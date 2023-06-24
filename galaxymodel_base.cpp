#include "galaxymodel_base.h"
#include "math_core.h"
#include "math_random.h"
#include "math_sample.h"
#include "math_specfunc.h"
#include "math_spline.h"
#include "math_linalg.h"
#include "potential_utils.h"
#include "actions_torus.h"
#include "smart.h"
#include "utils.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <cassert>
#include "alloca.h"

namespace galaxymodel{

const TrivialSelectionFunction trivialSF;

namespace{   // internal definitions

//------- HELPER ROUTINES -------//

/** map (0,1) to (0,infty) for distance and include s^2 factor to
 ** sample cone
*/
/*inline double unscaleDistance(const double x, double *jac){
	double s = x<1? x/(1-x) : 1e10;	//watch out for when var[0]=1 ,maybe account for this?
	if(jac != NULL) *jac = x<1? 1/pow_2(1-x) : 1e10;
	return s;
}*/

/** convert from scaled velocity variables to the actual velocity.
    \param[in]  vars are the scaled variables: chi, psi, phi/(2*Pi),
    where the magnitude of the velocity is v = v_esc * g(chi, zeta), g is some scaling function,
    and the two angles {theta(psi), phi} specify the orientation of velocity vector 
    in spherical coordinates centered at a given point.
    \param[in]  vesc is the maximum magnutude of velocity (equal to the escape velocity).
    \param[in]  zeta is the scaling factor for the velocity magnitude -
    the ratio of circular to escape velocity at the given radius.
    The non-trivial transformation is needed to accurately handle distribution functions of
    cold disks at large radii, which are very strongly peaked near {v_R,v_z,v_phi} = {0,0,v_circ}.
    To improve the robustness of integration, we make sure that a large proportion of unit cube
    in scaled variables maps onto a relatively small region around this circular velocity:
    the scaling function g(chi) is nearly horizontal for a large range of chi when its value
    is close to zeta, and the angle theta = pi * psi^2, i.e. again a large range of psi maps
    onto a small region of theta near zero, where the velocity is directed nearly azimuthally.
    \param[out] jac (optional) if not NULL, output the jacobian of transformation.
    \return  three components of velocity in cylindrical coordinates.
*/
inline coord::VelCyl unscaleVelocity(
    const double vars[], const double vesc, const double zeta, double* jac=0)
{
    double sintheta, costheta, sinphi, cosphi,
    eta = sqrt(1/zeta-1) + 1,
    chi = vars[0] * eta - 1,
    vel = vesc * zeta * (1 + math::sign(chi) * pow_2(chi));
    math::sincos(M_PI * pow_2(vars[1]), sintheta, costheta);
    math::sincos(2*M_PI * vars[2], sinphi, cosphi);
    if(jac)
        *jac = 8*M_PI*M_PI * vars[1] * zeta * fabs(chi) * eta * vesc * pow_2(vel) * sintheta;
    return coord::VelCyl(vel * sintheta * cosphi, vel * sintheta * sinphi, vel * costheta);
}

/** compute the escape velocity and the ratio of circular to escape velocity
    at a given position in the given ponential */
inline void getVesc(
    const coord::PosCyl& pos, const potential::BasePotential& poten, double& vesc, double& zeta)
{
    if(pow_2(pos.R) + pow_2(pos.z) == INFINITY) {
        vesc = 0.;
        zeta = 0.5;
        return;
    }
    double Phi;
    coord::GradCyl grad;
    poten.eval(pos, &Phi, &grad);
    vesc = sqrt(-2. * Phi);
    zeta = math::clip(sqrt(grad.dR * pos.R) / vesc, 0.1, 0.9);
    if(!isFinite(vesc)) {
        throw std::invalid_argument("Error in computing moments: escape velocity is undetermined at "
            "R="+utils::toString(pos.R)+", z="+utils::toString(pos.z)+", phi="+utils::toString(pos.phi)+
            " (Phi="+utils::toString(Phi)+")");
    }
}

/** convert from scaled position/velocity coordinates to the real ones.
    The coordinates in cylindrical system are scaled in the same way as for 
    the density integration; the velocity magnitude is scaled with local escape velocity.
    If needed, also provide the jacobian of transformation.
*/
inline coord::PosVelCyl unscalePosVel(const double vars[], 
    const potential::BasePotential& pot, double* jac=0)
{
    // 1. determine the position from the first three scaled variables
    double jacPos=0;
    const coord::PosCyl pos = potential::unscaleCoords(vars, jac==NULL ? NULL : &jacPos);
    // 2. determine the velocity from the second three scaled vars
    double vesc, zeta;
    getVesc(pos, pot, vesc, zeta);
    const coord::VelCyl vel = unscaleVelocity(vars+3, vesc, zeta, jac);
    if(jac!=NULL)
        *jac *= jacPos;
    return coord::PosVelCyl(pos, vel);
}

//------- HELPER CLASSES FOR MULTIDIMENSIONAL INTEGRATION OF DF -------//

/** Base helper class for integrating the distribution function over the position/velocity space.
    Various tasks in this module boil down to computing the integrals or sampling the values of DF
    over the (x,v) space, where the DF is expressed in terms of actions.
    This involves the following steps:
    1) scaled variables in N-dimensional unit cube are transformed to the actual (x,v);
    2) x,v are transformed to actions (J);
    3) the value(s) of DF f(J) is computed (for a multicomponent DF, one may either use the
    sum of all components, or treat them separately, depending on the flag 'separate');
    4) one or more quantities that are products of f(J) times something
    (e.g., velocity components) are returned to the integration or sampling routines.
    These tasks differ in the first and the last steps, and also in the number of dimensions
    that the integration/sampling is carried over. This diversity is handled by the class
    hierarchy descending from DFIntegrandNdim, where the base class performs the steps 2 and 3,
    and the derived classes implement virtual methods `unscaleVars()` and `outputValues()`,
    which are responsible for the steps 1 and 4, respectively.
    The derived classes also specify the dimensions of integration space (numVars)
    and the number of simultaneously computed quantities (numValues).
    The action finder for performing the step 2 and the DF used in the step 3 are provided
    as members of the GalaxyModel structure.
*/
class DFIntegrandNdim: public math::IFunctionNdim {
	public:
		double bright, faint;
		const obs::BaseLos* los;
		explicit DFIntegrandNdim(const GalaxyModel& _model, bool separate,
					 const obs::BaseLos* _los=NULL,
					 double _bright=NULL, double _faint=NULL)
				:
		    model(_model),
		    dflen(separate ? model.distrFunc.numValues() : 1),
		los(_los),
		bright(_bright),
		faint(_faint)
		{}

    /** compute one or more moments of distribution function. */
		    virtual void eval(const double vars[], double values[]) const
		    {
	// value(s) of distribution function - more than one for a multicomponent DF treated separately
	// (array allocated on the stack, no need to delete it manually)
			    double* dfval = static_cast<double*>(alloca(dflen * sizeof(double)));
			    std::fill(dfval, dfval+dflen, 0);  // initialize with zeroes (remain so in case of errors)
	// 1. get the position/velocity components in cylindrical coordinates
			    double jac;   // jacobian of variable transformation
			    coord::PosVelCyl posvel = unscaleVars(vars, &jac);
			    if(jac == 0) {  // we can't compute actions, but pretend that DF*jac is zero
				    outputValues(posvel, dfval, values);
				    return;
			    }

			    try{
	    // 2. determine the actions
				    actions::Actions act = model.actFinder.actions(posvel);

	    // 3. compute the value of distribution function times the jacobian
	    // FIXME: in some cases the Fudge action finder may fail and produce
	    // zero values of Jr,Jz instead of very large ones, which may lead to
	    // unrealistically high DF values. We therefore ignore these points
	    // entirely, but the real problem is with the action finder, not here.
				    bool valid = true;
				    if(isFinite(act.Jr + act.Jz + act.Jphi) && (act.Jr!=0 || act.Jz!=0)) {
					    if(los){
						    double sKpc = los->s0 + los->sKpc(posvel);
						    if(bright){//Compute abs mags
							    double shift = 5*log10(100 * sKpc);
							    shift += los->A_H(sKpc);
							    double Bright = bright - shift, Faint = faint - shift;//Now abs mags
							    model.distrFunc.withSF(act, dfval, Bright, Faint);
						    } else {//Just reduce star count for extinction
							    double fac = pow(10,-0.4*los->A_H(sKpc));
							    if(dflen==1) dfval[0]=model.distrFunc.value(act);
							    else model.distrFunc.eval(act, dfval);
							    for(int cpt=0; cpt<dflen; cpt++)
								    dfval[cpt] *= fac;
						    }
					    }
					    else{
						    if(dflen==1) dfval[0]=model.distrFunc.value(act);
						    else model.distrFunc.eval(act, dfval);  
					    }
				// multiply by jacobian, check for possibly invalid values and replace them with zeroes
					    for(unsigned int i=0; i<dflen; i++) {
						    if(!isFinite(dfval[i])) {
							    valid = false;
							    dfval[i] = 0;
						    } 
						    else
							    dfval[i] *= jac;
					    }
				    }
				    if(!valid)
					    throw std::runtime_error("DF is not finite");
			    }
			    catch(std::exception& e) {
	    // dump out the error at the highest debug logging level
				    std::cout << "DFIntegrandNdim" + std::string(e.what()) +
						    " at R="+utils::toString(posvel.R)  +", z="   +utils::toString(posvel.z)+
						    ", phi="+utils::toString(posvel.phi)+", vR="  +utils::toString(posvel.vR)+
						    ", vz=" +utils::toString(posvel.vz) +", vphi="+utils::toString(posvel.vphi)
						    << "\n";
				    if(utils::verbosityLevel >= utils::VL_VERBOSE) {
/*					    utils::msg(utils::VL_VERBOSE, "DFIntegrandNdim", std::string(e.what()) +
						    " at R="+utils::toString(posvel.R)  +", z="   +utils::toString(posvel.z)+
						    ", phi="+utils::toString(posvel.phi)+", vR="  +utils::toString(posvel.vR)+
						    ", vz="
						    +utils::toString(posvel.vz)
						    +", vphi="+utils::toString(posvel.vphi));*/
				    }
	    // ignore the error and proceed as if the DF was zero (keep the initially assigned zero)
			    }

	// 4. output the value(s) to the integration routine
			    outputValues(posvel, dfval, values);
		    }

    /** convert from scaled variables used in the integration routine 
        to the actual position/velocity point.
        \param[in]  vars  is the array of scaled variables;
        \param[out] jac (optional)  is the jacobian of transformation, if NULL it is not computed;
        \return  the position and velocity in cylindrical coordinates.
    */
		    virtual coord::PosVelCyl unscaleVars(const double vars[], double* jac=0) const = 0;

    /** output the value(s) computed at a given point to the integration routine.
        \param[in]  point  is the position/velocity point;
        \param[in]  dfval  is the value or array of values of distribution function at this point;
        \param[out] values is the array of one or more values that are computed
    */
		    virtual void outputValues(
					      const coord::PosVelCyl& point, const double dfval[], double values[]) const = 0;

		    const GalaxyModel& model;  ///< reference to the galaxy model to work with
    /// number of values in the DF in the case when a composite DF has more than component, and
    /// they are requested to be considered separately, otherwise 1 (a sum of all components)
		    unsigned int dflen;
};


/** helper class for computing the projected distribution function at a given point in x,y,vz space  */
class DFIntegrandProjected: public DFIntegrandNdim, public math::IFunctionNoDeriv {
public:
    DFIntegrandProjected(const GalaxyModel& _model, bool separate,
        double _R, double _vz, double _vz_error)
    :
        DFIntegrandNdim(_model, separate), R(_R), vz(_vz), vz_error(_vz_error) {}

    /// IFunctionNoDeriv: return v^2-vz^2 (used in setting the integration limits by root-finding)
    virtual double value(double zscaled) const
    {
        return -vz*vz + (zscaled==1 ? 0 : -2*model.potential.value(
            coord::PosCyl(R, math::unscale(math::ScalingSemiInf(), zscaled), 0)));
    }

    /// input variables define the missing components of position and velocity
    /// to be integrated over, suitably scaled: z, vx, vy
    virtual coord::PosVelCyl unscaleVars(const double vars[], double* jac=0) const
    {
        // map the input vars[0] ranging from 0 to 1 into z ranging from 0 to infinity
        double z   = math::unscale(math::ScalingSemiInf(), vars[0], jac);
        double vz1 = vz;
        if(vz_error!=0)  // add velocity error sampled from Gaussian c.d.f.
            vz1 += M_SQRT2 * vz_error * math::erfinv(2*vars[3]-1);
        // vx^2+vy^2 = -2 Phi(r) - vz^2
        double v2 = vars[0]==1 ? 0 : -2*model.potential.value(coord::PosCyl(R, z, 0)) - vz1*vz1;
        if(v2<=0) {    // we're outside the allowed range of z
            if(jac!=NULL)
                *jac = 0;
            return coord::PosVelCyl(R, 0, 0, 0, vz, 0);
        }
        double v = sqrt(v2) * vars[1], sinphi, cosphi;
        math::sincos(2*M_PI*vars[2], sinphi, cosphi);
        if(jac!=NULL)
            *jac *= 2*M_PI * v2 * vars[1] * 2;    // jacobian of velocity transformation
        return coord::PosVelCyl(R, z, 0, v * cosphi, vz1, v * sinphi);
    }

private:
    double R, vz, vz_error;
    virtual unsigned int numVars()   const { return vz_error==0 ? 3 : 4; }
    virtual unsigned int numValues() const { return dflen; }

    /// output array contains the value of DF (or all values for a multicomponent DF)
    virtual void outputValues(const coord::PosVelCyl& , const double dfval[], double values[]) const
    {
        std::copy(dfval, dfval+dflen, values);
    }
};


/** helper class for computing the moments of distribution function
    (surface density, scale height, line-of-sight velocity dispersion) at a given point in x,y plane  */
class DFIntegrandProjectedMoments: public DFIntegrandNdim {
public:
    DFIntegrandProjectedMoments(const GalaxyModel& _model, bool separate, double _R) :
        DFIntegrandNdim(_model, separate), R(_R) {}

    /// input variables define the z-coordinate and all three velocity components, suitably scaled
    virtual coord::PosVelCyl unscaleVars(const double vars[], double* jac=0) const
    {
        // map the input vars[0] ranging from 0 to 1 into z ranging from 0 to infinity
        coord::PosCyl pos(R, math::unscale(math::ScalingSemiInf(), vars[0], jac), 0);
        double vesc, zeta, jacVel;
        getVesc(pos, model.potential, vesc, zeta);
        const coord::VelCyl vel = unscaleVelocity(vars+1, vesc, zeta, &jacVel);
        if(jac!=NULL)
            *jac = vesc==0 ? 0 : *jac * jacVel * 2 /*since we only consider z>0*/;
        return coord::PosVelCyl(pos, vel);
    }

    virtual unsigned int numVars()   const { return 4; }
    virtual unsigned int numValues() const { return 3*dflen; }

private:
    double R;

    /// output array contains three elements per each DF component -
    /// the value of DF multiplied by 1, z^2 or vz^2
    virtual void outputValues(const coord::PosVelCyl& pv, const double dfval[], double values[]) const
    {
        for(unsigned int ic=0; ic<dflen; ic++) {
            values[ic + dflen*0] = dfval[ic];
            values[ic + dflen*1] = dfval[ic]!=0 ? dfval[ic] * pow_2(pv.z)  : 0;
            values[ic + dflen*2] = dfval[ic]!=0 ? dfval[ic] * pow_2(pv.vz) : 0;
        }
    }
};


/** helper class for integrating the distribution function over the entire 6d phase space */
class DFIntegrand6dim: public DFIntegrandNdim {
	public:
		DFIntegrand6dim(const GalaxyModel& _model) :
		    DFIntegrandNdim(_model, false) {}

    /// input variables define 6 components of position and velocity, suitably scaled
		virtual coord::PosVelCyl unscaleVars(const double vars[], double* jac=0) const
		{
			return unscalePosVel(vars, model.potential, jac);
		}

	private:
		virtual unsigned int numVars()   const { return 6; }
		virtual unsigned int numValues() const { return 1; }

    /// output contains just one value of DF (even for a multicomponent DF)
		virtual void outputValues(const coord::PosVelCyl& , const double dfval[], double values[]) const
		{
			values[0] = dfval[0];
		}
};


/// this will be redesigned
class DFIntegrandProjection: public math::IFunctionNdim {
    const GalaxyModel& model;
    const math::IFunctionNdim& fnc;  ///< spatial selection function
    const double* mat;               ///< orthogonal matrix for coordinate transformation
public:
    DFIntegrandProjection(const GalaxyModel& _model,
        const math::IFunctionNdim& _fnc, const double* _transformMatrix) :
    model(_model), fnc(_fnc), mat(_transformMatrix) {}

    virtual unsigned int numVars()   const { return 6; }
    virtual unsigned int numValues() const { return fnc.numValues(); }

    virtual void eval(const double vars[], double values[]) const
    {
        try{
            double X = vars[0], Y = vars[1], W = 2*vars[2]-1, Z, jac;  // W is scaled Z:
            if(W<0) {
                Z   = -exp(1/(1+W) + 1/W);
                jac = -Z * (1/pow_2(1+W) + 1/pow_2(W)) * 2;  // dZ/dW
            } else if(W>0) {
                Z   =  exp(1/(1-W) - 1/W);
                jac =  Z * (1/pow_2(1-W) + 1/pow_2(W)) * 2;
            } else {
                Z = jac = 0;
            }

            // transform the position from observed to intrinsic frame
            const coord::PosCyl pos(coord::toPosCyl(coord::PosCar(
                mat[0] * X + mat[3] * Y + mat[6] * Z,
                mat[1] * X + mat[4] * Y + mat[7] * Z,
                mat[2] * X + mat[5] * Y + mat[8] * Z) ) );

            // construct the full position/velocity in intrinsic frame
            double vesc, zeta, jacVel;
            getVesc(pos, model.potential, vesc, zeta);
            const coord::PosVelCyl posvel(pos, unscaleVelocity(vars+3, vesc, zeta, &jacVel));
            jac = isFinite(jacVel) && jac!=0 ?  jac * jacVel  : 0;

            // transform the velocity back to observed frame
            coord::VelCar velrot(toPosVelCar(posvel));
            const double posvelrot[6] = { X, Y, Z,
                mat[0] * velrot.vx + mat[1] * velrot.vy + mat[2] * velrot.vz,
                mat[3] * velrot.vx + mat[4] * velrot.vy + mat[5] * velrot.vz,
                mat[6] * velrot.vx + mat[7] * velrot.vy + mat[8] * velrot.vz };

            // query the spatial selection function
            fnc.eval(posvelrot, values);

            // check if there are any nonzero values reported by fnc
            // [skipped for the moment]

            // 2. determine the actions
            actions::Actions act = model.actFinder.actions(posvel);

            // 3. compute the value of distribution function times the jacobian
            // FIXME: in some cases the Fudge action finder may fail and produce
            // zero values of Jr,Jz instead of very large ones, which may lead to
            // unrealistically high DF values. We therefore ignore these points
            // entirely, but the real problem is with the action finder, not here.
            double dfval = isFinite(act.Jr + act.Jz + act.Jphi) && (act.Jr!=0 || act.Jz!=0) ?
                model.distrFunc.value(act) * jac : 0.;

            if(!isFinite(dfval))
                dfval = 0;

            // 4. output the value(s) to the integration routine
            for(unsigned int i=0, count=fnc.numValues(); i<count; i++)
                values[i] *= dfval;
        }
        catch(std::exception& e) {
            for(unsigned int i=0, count=fnc.numValues(); i<count; i++)
                values[i] = 0;
        }
    }
};


/** specification of the velocity moments of DF to be computed at a single point in space
    (a combination of them is given by bitwise OR) */
enum OperationMode {
    OP_VEL1MOM=1,  ///< first moments of velocity
    OP_VEL2MOM=2,   ///< second moments of velocity
    OP_VEL2FRQ=3   ///< mean frequencies
};

/** helper class for integrating the distribution function over velocity at a fixed position */
class DFIntegrandAtPoint: public DFIntegrandNdim {
	private:
		const coord::PosCyl point;    ///< fixed position
		double vesc;                  ///< escape velocity at this position
		double zeta;                  ///< the ratio of circular to escape velocity
		const OperationMode mode;     ///< determines which moments of DF to compute
		const unsigned int numOutVal; ///< number of output values for each component of DF
	public:
		DFIntegrandAtPoint(const GalaxyModel& _model, bool separate,
				   const coord::PosCyl& _point, OperationMode _mode,
				   obs::BaseLos* _los = NULL, double bright = NULL, double faint = NULL)
				:
		    DFIntegrandNdim(_model, separate, _los, bright, faint),
		    point(_point), mode(_mode),
		    numOutVal(1 + (mode&OP_VEL1MOM ? 1 : 0) + (mode&OP_VEL2MOM ? 6 : 0)
			      + (mode&OP_VEL2FRQ ? 3 : 0))
		{
			getVesc(point, model.potential, vesc, zeta);
		}

		    virtual coord::PosVelCyl unscaleVars(const double vars[], double* jac=0) const
		    {
			    return coord::PosVelCyl(point, unscaleVelocity(vars, vesc, zeta, jac));
		    }

		    virtual void outputValues(const coord::PosVelCyl& pv, const double dfval[], double values[]) const
		    {
	// 4. output the value(s) of DF, multiplied by various combinations of velocity components:
	// {f, f*vR, f*vz, f*vphi, f*vR^2, f*vz^2, f*vphi^2, f*vR*vz, f*vR*vphi, f*vz*vphi },
	// depending on the mode of operation.
			    for(unsigned int ic=0; ic<dflen; ic++) {  // loop over components of DF
				    values[ic]   = dfval[ic];
				    unsigned int im=1;      // index of the output moment, increases with each stored value
				    if(mode & OP_VEL1MOM) {
					    values[ic + dflen * (im++)] = dfval[ic] * pv.vphi;  // only <v_phi> may be nonzero
				    }
				    if(mode & OP_VEL2MOM) {
					    values[ic + dflen * (im++)] = dfval[ic] * pv.vR   * pv.vR;
					    values[ic + dflen * (im++)] = dfval[ic] * pv.vz   * pv.vz;
					    values[ic + dflen * (im++)] = dfval[ic] * pv.vphi * pv.vphi;
					    values[ic + dflen * (im++)] = dfval[ic] * pv.vR   * pv.vz;
					    values[ic + dflen * (im++)] = dfval[ic] * pv.vR   * pv.vphi;
					    values[ic + dflen * (im++)] = dfval[ic] * pv.vz   * pv.vphi;
				    }
				    if(mode & OP_VEL2FRQ) {
					    actions::Frequencies Freqs;
					    model.actFinder.actionAngles(pv,&Freqs);
					    values[ic + dflen * (im++)] = dfval[ic] * Freqs.Omegar;
					    values[ic + dflen * (im++)] = dfval[ic] * Freqs.Omegaz;
					    values[ic + dflen * (im++)] = dfval[ic] * fabs(Freqs.Omegaphi);
				    }
			    }
		    }

    /// dimension of the input array (3 scaled velocity components)
		    virtual unsigned int numVars()   const { return 3; }

    /// dimension of the output array
		    virtual unsigned int numValues() const { return dflen * numOutVal; }
};


/** Helper class for constructing histograms of velocity distribution */
template <int N>
class DFIntegrandVelDist: public DFIntegrandNdim {
public:
    DFIntegrandVelDist(const GalaxyModel& _model, bool separate,
        const coord::PosCyl& _point, bool _projected, 
        const math::BsplineInterpolator1d<N>& _bsplVR,
        const math::BsplineInterpolator1d<N>& _bsplVz,
        const math::BsplineInterpolator1d<N>& _bsplVphi)
    :
        DFIntegrandNdim(_model, separate),
        point(_point),
        projected(_projected),
        bsplVR(_bsplVR), bsplVz(_bsplVz), bsplVphi(_bsplVphi),
        NR(bsplVR.numValues()), Nz(bsplVz.numValues()), Ntotal(1 + NR + Nz + bsplVphi.numValues())
    {
        if(!projected)
            getVesc(point, model.potential, vesc, zeta);
    }

    virtual unsigned int numVars()   const { return projected ? 4 : 3; }
    virtual unsigned int numValues() const { return dflen * Ntotal; }
private:
    const coord::PosCyl point;  ///< position
    const bool projected;       ///< if true, only use R and phi and integrate over z
    double vesc, zeta;          ///< escape velocity at this position (if not projected)
    const math::BsplineInterpolator1d<N>& bsplVR, bsplVz, bsplVphi;
    const unsigned int NR, Nz, Ntotal;

    /// input variables define the z-coordinate and all three velocity components, suitably scaled
    virtual coord::PosVelCyl unscaleVars(const double vars[], double* jac=0) const
    {
        if(projected) {
            coord::PosCyl pos(point.R, math::unscale(math::ScalingInf(), vars[0], jac), point.phi);
            double vesc, zeta, jacVel;
            getVesc(pos, model.potential, vesc, zeta);
            const coord::VelCyl vel = unscaleVelocity(vars+1, vesc, zeta, &jacVel);
            if(jac!=NULL)
                *jac = vesc==0 ? 0 : *jac * jacVel;
            return coord::PosVelCyl(pos, vel);
        } else
            return coord::PosVelCyl(point, unscaleVelocity(vars, vesc, zeta, jac));
    }

    /// output the weighted integrals over basis functions;
    /// we scan only half of the (v_R, v_z) plane, and add the same contributions to (-v_R, -v_z),
    /// since the actions and hence the value of f(J) do not change with this inversion
    virtual void outputValues(const coord::PosVelCyl& pv, const double dfval[], double values[]) const
    {
        std::fill(values, values + dflen * Ntotal, 0);
        double valRp[N+1], valRm[N+1], valzp[N+1], valzm[N+1], valphi[N+1];
        unsigned int 
            iRp  = bsplVR.  nonzeroComponents( pv.vR,  0, valRp),
            iRm  = bsplVR.  nonzeroComponents(-pv.vR,  0, valRm),
            izp  = bsplVz.  nonzeroComponents( pv.vz,  0, valzp),
            izm  = bsplVz.  nonzeroComponents(-pv.vz,  0, valzm),
            iphi = bsplVphi.nonzeroComponents(pv.vphi, 0, valphi);
        for(unsigned int ic=0; ic<dflen; ic++) {  // loop over DF components
            // factor of two here because we only integrate over half-space in v_z
            values[ic * Ntotal] = dfval[ic] * 2;
            for(int ib=0; ib<=N; ib++) {  // loop over nonzero basis functions
                values[ic * Ntotal + 1 + ib + iRp]           += dfval[ic] * valRp[ib];
                values[ic * Ntotal + 1 + ib + iRm]           += dfval[ic] * valRm[ib];
                values[ic * Ntotal + 1 + ib + izp  + NR]     += dfval[ic] * valzp[ib];
                values[ic * Ntotal + 1 + ib + izm  + NR]     += dfval[ic] * valzm[ib];
                values[ic * Ntotal + 1 + ib + iphi + NR + Nz] = dfval[ic] * valphi[ib] * 2;
            }
        }
    }
};

/** Solve the equation for amplitudes of B-spline expansion of the velocity distribution function.
    The VDF is represented as  \f$  f(v) = \sum_i A_i B_i(v)  \f$,  where B_i(v) are the B-spline basis
    functions and A_i are the amplitudes to be found by solving the following linear system:
    \f$  \int f(v) B_j(v) dv = \sum_i A_i  [ \int B_i(v) B_j(v) dv ] = C_j  \f$,
    where C_j is the RHS vector computed through the integration of f(v) weighted with each
    basis function, and the overlap matrix in square brackets is provided by the B-spline object.
    Even though the RHS by definition is non-negative, the solution vector is not guaranteed to be so
    (unless the matrix is diagonal, which is the case only for N=0, i.e., a histogram representation);
    that is, the interpolated f(v) may attain unphysical negative values.
    We employ an additional measure that helps to reduce this effect:
    if the order of B-spline interpolator is larger than zero (i.e., it's not a simple histogram), and
    if the endpoints of the velocity interval are at the escape velocity (meaning that f(v) must be 0),
    we enforce the amplitudes of the first and the last basis functions to be zero.
    In this case, of course, the number of variables in the system is less than the number of equations,
    so it is solved in the least-square sense using the singular-value decomposition, instead of
    the standard Cholesky decomposition for a full-rank symmetric matrix.
    \param[in]  bspl  is the B-spline basis in the 1d velocity space;
    \param[in]  rhs   is the RHS of the linear system;
    \param[in]  vesc  is the escape velocity: if the endpoints of the B-spline interval are at or beyond
    the escape velocity, we force the corresponding amplitudes to be zero.
    \return  the vector of amplitudes.
*/
template<int N>
std::vector<double> solveForAmplitudes(const math::BsplineInterpolator1d<N>& bspl,
    const std::vector<double>& rhs, double vesc)
{
    math::BandMatrix<double> bandMat = math::FiniteElement1d<N>(bspl.xvalues()).computeProjMatrix();
    int size = bandMat.rows();
#if 0
    // another possibility is to use a linear or quadratic optimization solver that enforce
    // non-negativity constraints on the solution vector, at the same time minimizing the deviation
    // from the exact solution. This is perhaps the cleanest approach, but it requires
    // the optimization libraries to be included at the compile time, which we can't guarantee.
    return math::quadraticOptimizationSolveApprox(bandMat, rhs,
        std::vector<double>()          /*linear penalty for variables is absent*/,
        math::BandMatrix<double>()     /*quadratic penalty for variables is absent*/,
        std::vector<double>()          /*linear penalty for constraint violation is absent*/,
        std::vector<double>(size, 1e9) /*quadratic penalty for constraint violation*/,
        std::vector<double>(size, 0.)  /*lower bound on the solution vector*/ );
#endif
    int skipFirst = (N >= 1 && size>2 && math::fcmp(bspl.xmin(), -vesc, 1e-8) <= 0);
    int skipLast  = (N >= 1 && size>2 && math::fcmp(bspl.xmax(),  vesc, 1e-8) >= 0);
    if(skipFirst+skipLast==0)
        return math::solveBand(bandMat, rhs);
    // otherwise create another matrix with fewer columns (copy row-by-row from the original matrix)
    math::Matrix<double> fullMat(bandMat);
    math::Matrix<double> reducedMat(size, size-skipFirst-skipLast, 0.);
    for(int i=0; i<size; i++)
        std::copy(&fullMat(i, skipFirst), &fullMat(i, size-skipLast), &reducedMat(i, 0));
    // use the SVD to solve the rank-deficient system
    std::vector<double> sol = math::SVDecomp(reducedMat).solve(rhs);
    // append the skipped amplitudes
    if(skipFirst)
        sol.insert(sol.begin(), 0);
    if(skipLast)
        sol.insert(sol.end(),   0);
    return sol;
}

/** helper class for integrating the distribution function along a LOS */
class DFIntegrandLOS: public DFIntegrandNdim {
	private:
		double vesc;                  ///< escape velocity at this position
		double zeta;                  ///< the ratio of circular to escape velocity
		const OperationMode mode;     ///< determines which moments of DF to compute
		const unsigned int numOutVal; ///< number of output values for each component of DF
	public:
		DFIntegrandLOS(const GalaxyModel& _model, bool separate,	
			       const obs::BaseLos* _los, OperationMode _mode,
			       const double bright = NULL, const double faint= NULL)
				:
		    DFIntegrandNdim(_model, separate, _los, bright, faint),		//seperate means it does the model components seperately 
		    mode(_mode),			
		    numOutVal(1 + (mode&OP_VEL1MOM ? 1 : 0) + (mode&OP_VEL2MOM ? 6 : 0))
		{
			getVesc(los->deepest(), model.potential, vesc, zeta);
		}

		    virtual coord::PosVelCyl unscaleVars(const double vars[], double* jac=0) const
		    {
			    double jacVel, s, cone;
			    if(los->semi_inf){
				    s = math::unscale(math::ScalingSemiInf(), vars[0], jac);
				    cone=1;
			    } else {
				    s = math::unscale(math::ScalingInf(), vars[0], jac);
				    cone = s*s;
			    }
			    coord::PosCyl pos = los->Rzphi(s);
			    coord::VelCyl vel = unscaleVelocity(vars+1, vesc, zeta, &jacVel);
			    if(jac) (*jac) *= jacVel * cone;
			    return coord::PosVelCyl(pos,vel);
		    }

		    virtual void outputValues(const coord::PosVelCyl& pv, const double dfval[], double values[]) const	//pv is structure of position/velocity coords
		    {
	// 4. output the value(s) of DF, multiplied by various combinations of velocity components:
	// {f, f*vR, f*vz, f*vphi, f*vR^2, f*vz^2, f*vphi^2, f*vR*vz, f*vR*vphi, f*vz*vphi },
	// depending on the mode of operation.
			    for(unsigned int ic=0; ic<dflen; ic++) {  // loop over components of DF
				    values[ic]   = dfval[ic];
				    unsigned int im=1;      // index of the output moment, increases with each stored value
				    if(mode & OP_VEL1MOM) {
					    values[ic + dflen * (im++)] = dfval[ic] * pv.vphi;  // only <v_phi> may be nonzero
				    }
				    if(mode & OP_VEL2MOM) {
					    values[ic + dflen * (im++)] = dfval[ic] * pv.vR   * pv.vR;
					    values[ic + dflen * (im++)] = dfval[ic] * pv.vz   * pv.vz;
					    values[ic + dflen * (im++)] = dfval[ic] * pv.vphi * pv.vphi;
					    values[ic + dflen * (im++)] = dfval[ic] * pv.vR   * pv.vz;
					    values[ic + dflen * (im++)] = dfval[ic] * pv.vR   * pv.vphi;
					    values[ic + dflen * (im++)] = dfval[ic] * pv.vz   * pv.vphi;
				    }
			    }
		    }

    /// dimension of the input array (3 scaled velocity components)				//needs to be 4 instead of 3
		    virtual unsigned int numVars()   const { return 4; }

    /// dimension of the output array
		    virtual unsigned int numValues() const { return dflen * numOutVal; }
};

/* The next section contains classes enabling one to sample models with
 * composite DFs such that each star is assigned a magnitude as well as
 * phase-space coordinates. All classes descend from
 * DF_LFIntegrandNdim, which is slightly exteded from DFIntegrandNdim.
 * The inheriting classes are generally sipler than the classes that
 * inherit DFIntegrandNdim because they are only intended for sampling
 * rather than evaluating moments.
 */

inline double unscaleMag(const double x, const double bright, const double faint, double* jac){
	double d = (faint - bright);
	if(jac) *jac = d;
	return bright + x*d;
}

class DF_LFIntegrandNdim: public math::IFunctionNdim {
	public:
		double bright, faint;
		const obs::BaseLos* los;
		explicit DF_LFIntegrandNdim(const GalaxyModel& _model, bool separate,
					 const obs::BaseLos* _los=NULL, double _bright=NULL, double _faint=NULL)
				:
		    model(_model),
		    dflen(separate ? model.distrFunc.numValues() : 1),
		los(_los),
		bright(_bright),
		faint(_faint)
		{}

    /** compute one or more moments of distribution function. */
		virtual void eval(const double vars[], double values[]) const
		    {
	// value(s) of distribution function - more than one for a multicomponent DF treated separately
	// (array allocated on the stack, no need to delete it manually)
			    double* dfval = static_cast<double*>(alloca(dflen * sizeof(double)));
			    std::fill(dfval, dfval+dflen, 0);  // initialize with zeroes (remain so in case of errors)
	// 1. get the position/velocity components in cylindrical coordinates
			    double jac;   // jacobian of variable transformation
			    coord::PosVelCylMag posvel = unscaleVars(vars, &jac);
			    if(jac == 0) {  // we can't compute actions, but pretend that DF*jac is zero
				    outputValues(posvel, dfval, values);
				    return;
			    }

			    try{
	    // 2. determine the actions
				    actions::Actions act = model.actFinder.actions(posvel.first);

	    // 3. compute the value of distribution function times the jacobian
	    // FIXME: in some cases the Fudge action finder may fail and produce
	    // zero values of Jr,Jz instead of very large ones, which may lead to
	    // unrealistically high DF values. We therefore ignore these points
	    // entirely, but the real problem is with the action finder, not here.
				    bool valid = true;
				    if(isFinite(act.Jr + act.Jz + act.Jphi) && (act.Jr!=0 || act.Jz!=0)) {
					    if(bright){//proposed mag is apparent
						    double s = los->s(posvel.first);
						    double shift = 5*log10(100 * s);
						    shift += los->A_H(s);
						    double Mag=posvel.second - shift;
						    model.distrFunc.withLF(act, dfval, Mag);
					    } 
					    else{//proposed mag is absolute
						    model.distrFunc.withLF(act, dfval, posvel.second);  
					    }
				// multiply by jacobian, check for possibly invalid values and replace them with zeroes
					    for(unsigned int i=0; i<dflen; i++) {
						    if(!isFinite(dfval[i])) {
							    valid = false;
							    dfval[i] = 0;
						    } 
						    else
							    dfval[i] *= jac;
					    }
				    }
				    if(!valid)
					    throw std::runtime_error("DF is not finite");
			    }
			    catch(std::exception& e) {
	    // dump out the error at the highest debug logging level
				    if(utils::verbosityLevel >= utils::VL_VERBOSE) {
					    utils::msg(utils::VL_VERBOSE, "DFIntegrandNdim", std::string(e.what()) +
						    " at R="+utils::toString(posvel.first.R)  +", z="   +utils::toString(posvel.first.z)+
						    ", phi="+utils::toString(posvel.first.phi)+", vR="  +utils::toString(posvel.first.vR)+
						    ", vz=" +utils::toString(posvel.first.vz) +", vphi="+utils::toString(posvel.first.vphi));
				    }
	    // ignore the error and proceed as if the DF was zero (keep the initially assigned zero)
			    }

	// 4. output the value(s) to the integration routine
			    outputValues(posvel, dfval, values);
		    }

    /** convert from scaled variables used in the integration routine 
        to the actual position/velocity point.
        \param[in]  vars  is the array of scaled variables;
        \param[out] jac (optional)  is the jacobian of transformation, if NULL it is not computed;
        \return  the position and velocity in cylindrical coordinates.
    */
		    virtual coord::PosVelCylMag unscaleVars(const double vars[], double* jac=0) const = 0;

    /** output the value(s) computed at a given point to the integration routine.
        \param[in]  point  is the position/velocity point;
        \param[in]  dfval  is the value or array of values of distribution function at this point;
        \param[out] values is the array of one or more values that are computed
    */
		    virtual void outputValues(const coord::PosVelCylMag& point,
					      const double dfval[], double values[]) const = 0;

		    const GalaxyModel& model;  ///< reference to the galaxy model to work with
    /// number of values in the DF in the case when a composite DF has more than component, and
    /// they are requested to be considered separately, otherwise 1 (a sum of all components)
		    unsigned int dflen;
};

/* The next class is what's needed to obtain initial conditions for an
 * N-body simulation  with each star assigned an absolute magnitude.
 */
class DF_LFIntegrand6dim : public DF_LFIntegrandNdim {
	public:
		DF_LFIntegrand6dim(const GalaxyModel& _model) :
		    DF_LFIntegrandNdim(_model, false) {}

    /// input variables define 6 scaled components of position and velocity plus abs magnitude
		virtual coord::PosVelCylMag unscaleVars(const double vars[], double* jac=0) const
		{
			coord::PosVelCyl point(unscalePosVel(vars, model.potential, jac));
			const double Bright=-7, Faint=12;//brightest & faintest abs mags
			double jacMag, Mag = unscaleMag(vars[6], Bright, Faint, &jacMag);
			if(jac) (*jac) *= jacMag;
			return coord::PosVelCylMag(point, Mag);
		}

	private:
		virtual unsigned int numVars()   const { return 7; }
		virtual unsigned int numValues() const { return 1; }

    /// output contains just one value of DF (even for a multicomponent DF)
		virtual void outputValues(const coord::PosVelCylMag& ,
					  const double dfval[], double values[]) const
		{
			values[0] = dfval[0];
		}
};

/* The next class samples a LOS assigning to each star an apparent mag
 * that lies between the given limiting mags of the target survey. The
 * returned mags reflect extinction as well as distance.
 */
class DF_LFIntegrandLOS: public DF_LFIntegrandNdim {
	private:
		double vesc;                  ///< escape velocity at this position
		double zeta;                  ///< the ratio of circular to escape velocity
	public:
		DF_LFIntegrandLOS(const GalaxyModel& _model,
			       const obs::BaseLos* _los,
			       const double bright, const double faint)
				:
		    DF_LFIntegrandNdim(_model, false, _los, bright, faint)
		{
			getVesc(los->deepest(), model.potential, vesc, zeta);
		}

		    virtual coord::PosVelCylMag unscaleVars(const double vars[], double* jac=0) const
		    {
			    double jacVel, s;
			    if(los->semi_inf) s = math::unscale(math::ScalingSemiInf(), vars[0],jac);
			    else s = math::unscale(math::ScalingInf(), vars[0],jac);
			    double jacMag, mag = unscaleMag(vars[4],bright,faint,&jacMag);
			    coord::PosCyl pos = los->Rzphi(s);
			    coord::VelCyl vel = unscaleVelocity(vars+1, vesc, zeta, &jacVel);
			    if(jac) (*jac) *= jacMag * jacVel;
			    return coord::PosVelCylMag(coord::PosVelCyl(pos,vel), mag);
		    }

	// 4. output the value(s) of DF, multiplied by various combinations of velocity components:
	// {f, f*vR, f*vz, f*vphi, f*vR^2, f*vz^2, f*vphi^2, f*vR*vz, f*vR*vphi, f*vz*vphi },
	// depending on the mode of operation.
		    virtual void outputValues(const coord::PosVelCylMag& pv,
					      const double dfval[], double values[]) const
		    {
			    values[0] = dfval[0];//separate = false so only 1 value
		    }

    /// dimension of the input array (distance, 3 scaled velocity
    /// components, magnitude)
		    virtual unsigned int numVars()   const { return 5; }

    /// dimension of the output array
		    virtual unsigned int numValues() const { return 1; }
};

/* The next class enables sampling of the velocity distribution at a
 * given location and assigning to each an apparent magnitude star that
 * lies between given limiting mags (those of the target survey).
 * these magnitudes reflect extinction between the Sun and the given
 * location.
 */ 
class DF_LFIntegrandAtPoint: public DF_LFIntegrandNdim {
	private:
		const coord::PosCyl point;    ///< fixed position
		double vesc;                  ///< escape velocity at this position
		double zeta;                  ///< the ratio of circular to escape velocity
		const unsigned int numOutVal; ///< number of output values for each component of DF
	public:
		DF_LFIntegrandAtPoint(const GalaxyModel& _model,
				   const coord::PosCyl& _point,
				   obs::BaseLos* _los, double bright, double faint)
				:
		    DF_LFIntegrandNdim(_model, false, _los, bright, faint),
		    point(_point),
		    numOutVal(1)
		{
			getVesc(point, model.potential, vesc, zeta);
		}
		    
		    virtual coord::PosVelCylMag unscaleVars(const double vars[], double* jac=0) const
		    {
			    double jacMag, mag = unscaleMag(vars[3],bright,faint,&jacMag);
			    coord::PosVelCyl pv(point, unscaleVelocity(vars, vesc, zeta, jac));
			    if(jac) (*jac) *= jacMag;
			    return coord::PosVelCylMag(pv, mag);
		    }

		    virtual void outputValues(const coord::PosVelCylMag& pv,
					      const double dfval[], double values[]) const
		    {
	// 4. output the value(s) of DF, multiplied by various combinations of velocity components:
	// {f, f*vR, f*vz, f*vphi, f*vR^2, f*vz^2, f*vphi^2, f*vR*vz, f*vR*vphi, f*vz*vphi },
	// depending on the mode of operation.
			    values[0] = dfval[0];// separate = false so just one val
		    }

    /// dimension of the input array (3 scaled velocity components)
		    virtual unsigned int numVars()   const { return 4; }

    /// dimension of the output array
		    virtual unsigned int numValues() const { return 1; }
};

}  // unnamed namespace

//------- DRIVER ROUTINES -------//

EXP void computeMomentsLOS(const GalaxyModel& model, const obs::BaseLos* los,
    double* density,    double* velocityFirstMoment,    coord::Vel2Cyl* velocitySecondMoment,
    double* densityErr, double* velocityFirstMomentErr, coord::Vel2Cyl* velocitySecondMomentErr,
    bool separate, const double reqRelError, const int maxNumEval,const double bright, const double faint)
{
    OperationMode mode = static_cast<OperationMode>(
        (velocityFirstMoment !=NULL ? OP_VEL1MOM : 0) |
	    (velocitySecondMoment!=NULL ? OP_VEL2MOM : 0) );
    DFIntegrandLOS fnc(model, separate, los, mode, bright, faint);
    // the integration region in scaled velocities
    double xlower[4] = {0, 0, 0, 0};
    double xupper[4] = {1, 1, 1, 1};
    // the values of integrals and their error estimates: allocate a temp array on the stack
    double* result = static_cast<double*>(alloca(fnc.numValues() * 2 * sizeof(double)));
    double* error  = result + fnc.numValues();   // second half of the temporary array
    try{
	    math::integrateNdim(fnc, xlower, xupper, reqRelError, maxNumEval, result, error);
    }
    catch(std::exception& e)
    {
	    printf("%s",e.what());
    }
    // store the results
    for(unsigned int ic=0; ic<fnc.dflen; ic++) {
        if(density!=NULL) {
            density[ic] = result[ic];
            if(densityErr!=NULL)
                densityErr[ic] = error[ic];
        }
        double densVal = result[ic], densRelErr2 = pow_2(error[ic]/result[ic]);
        unsigned int im=1;  // index of the computed moment in the results array
        if(velocityFirstMoment!=NULL) {
            velocityFirstMoment[ic] = densVal==0 ? 0 : result[ic + im*fnc.dflen] / densVal;
            if(velocityFirstMomentErr!=NULL) {
                // relative errors in moments are summed in quadrature from errors in rho and rho*v
                velocityFirstMomentErr[ic] = fabs(velocityFirstMoment[ic]) *
                    sqrt(densRelErr2 + pow_2(error[ic + im*fnc.dflen] / result[ic + im*fnc.dflen]) );
            }
            im++;
        }
        if(velocitySecondMoment!=NULL) {
            velocitySecondMoment[ic].vR2    = densVal ? result[ic + (im+0)*fnc.dflen] / densVal : 0;
            velocitySecondMoment[ic].vz2    = densVal ? result[ic + (im+1)*fnc.dflen] / densVal : 0;
            velocitySecondMoment[ic].vphi2  = densVal ? result[ic + (im+2)*fnc.dflen] / densVal : 0;
            velocitySecondMoment[ic].vRvz   = densVal ? result[ic + (im+3)*fnc.dflen] / densVal : 0;
            velocitySecondMoment[ic].vRvphi = densVal ? result[ic + (im+4)*fnc.dflen] / densVal : 0;
            velocitySecondMoment[ic].vzvphi = densVal ? result[ic + (im+5)*fnc.dflen] / densVal : 0;
            if(velocitySecondMomentErr!=NULL) {
                velocitySecondMomentErr[ic].vR2 =
                    fabs(velocitySecondMoment[ic].vR2) * sqrt(densRelErr2 +
                    pow_2(error[ic + (im+0)*fnc.dflen] / result[ic + (im+0)*fnc.dflen]) );
                velocitySecondMomentErr[ic].vz2 =
                    fabs(velocitySecondMoment[ic].vz2) * sqrt(densRelErr2 +
                    pow_2(error[ic + (im+1)*fnc.dflen] / result[ic + (im+1)*fnc.dflen]) );
                velocitySecondMomentErr[ic].vphi2 =
                    fabs(velocitySecondMoment[ic].vphi2) * sqrt(densRelErr2 +
                    pow_2(error[ic + (im+2)*fnc.dflen] / result[ic + (im+2)*fnc.dflen]) );
                velocitySecondMomentErr[ic].vRvz =
                    fabs(velocitySecondMoment[ic].vRvz) * sqrt(densRelErr2 +
                    pow_2(error[ic + (im+3)*fnc.dflen] / result[ic + (im+3)*fnc.dflen]) );
                velocitySecondMomentErr[ic].vRvphi =
                    fabs(velocitySecondMoment[ic].vRvphi) * sqrt(densRelErr2 +
                    pow_2(error[ic + (im+4)*fnc.dflen] / result[ic + (im+4)*fnc.dflen]) );
                velocitySecondMomentErr[ic].vzvphi =
                    fabs(velocitySecondMoment[ic].vzvphi) * sqrt(densRelErr2 +
                    pow_2(error[ic + (im+5)*fnc.dflen] / result[ic + (im+5)*fnc.dflen]) );
            }
        }
    }
}

EXP void computeMoments(const GalaxyModel& model, const coord::PosCyl& point,
    double* density,    double* velocityFirstMoment,    coord::Vel2Cyl* velocitySecondMoment,
    actions::Frequencies* frequencyMoment,
    double* densityErr, double* velocityFirstMomentErr, coord::Vel2Cyl* velocitySecondMomentErr,
			bool separate, const double reqRelError, const int maxNumEval,
			const double bright, const double faint, obs::solarShifter* sun)
{
    OperationMode mode = static_cast<OperationMode>(
        (velocityFirstMoment !=NULL ? OP_VEL1MOM : 0) |
        (velocitySecondMoment!=NULL ? OP_VEL2MOM : 0) |
        (frequencyMoment     !=NULL ? OP_VEL2FRQ : 0));
    obs::BaseLos* ptrLos=NULL;
    if(sun != NULL){
	    obs::SunLos Los(coord::toPosCar(point), *sun);
	    ptrLos=&Los;
    }
    DFIntegrandAtPoint fnc(model, separate, point, mode, ptrLos, bright, faint);

    // the integration region in scaled velocities
    double xlower[3] = {0, 0, 0};
    double xupper[3] = {1, 1, 1};
    // the values of integrals and their error estimates: allocate a temp array on the stack
    double* result = static_cast<double*>(alloca(fnc.numValues() * 2 * sizeof(double)));
    double* error  = result + fnc.numValues();   // second half of the temporary array

    math::integrateNdim(fnc, xlower, xupper, reqRelError, maxNumEval, result, error);

    // store the results
    for(unsigned int ic=0; ic<fnc.dflen; ic++) {
        if(density!=NULL) {
            density[ic] = result[ic];
            if(densityErr!=NULL)
                densityErr[ic] = error[ic];
        }
        double densVal = result[ic], densRelErr2 = pow_2(error[ic]/result[ic]);
        unsigned int im=1;  // index of the computed moment in the results array
        if(velocityFirstMoment!=NULL) {
            velocityFirstMoment[ic] = densVal==0 ? 0 : result[ic + im*fnc.dflen] / densVal;
            if(velocityFirstMomentErr!=NULL) {
                // relative errors in moments are summed in quadrature from errors in rho and rho*v
                velocityFirstMomentErr[ic] = fabs(velocityFirstMoment[ic]) *
                    sqrt(densRelErr2 + pow_2(error[ic + im*fnc.dflen] / result[ic + im*fnc.dflen]) );
            }
            im++;
        }
        if(velocitySecondMoment!=NULL) {
            velocitySecondMoment[ic].vR2    = densVal ? result[ic + (im+0)*fnc.dflen] / densVal : 0;
            velocitySecondMoment[ic].vz2    = densVal ? result[ic + (im+1)*fnc.dflen] / densVal : 0;
            velocitySecondMoment[ic].vphi2  = densVal ? result[ic + (im+2)*fnc.dflen] / densVal : 0;
            velocitySecondMoment[ic].vRvz   = densVal ? result[ic + (im+3)*fnc.dflen] / densVal : 0;
            velocitySecondMoment[ic].vRvphi = densVal ? result[ic + (im+4)*fnc.dflen] / densVal : 0;
            velocitySecondMoment[ic].vzvphi = densVal ? result[ic + (im+5)*fnc.dflen] / densVal : 0;
            if(velocitySecondMomentErr!=NULL) {
                velocitySecondMomentErr[ic].vR2 =
                    fabs(velocitySecondMoment[ic].vR2) * sqrt(densRelErr2 +
                    pow_2(error[ic + (im+0)*fnc.dflen] / result[ic + (im+0)*fnc.dflen]) );
                velocitySecondMomentErr[ic].vz2 =
                    fabs(velocitySecondMoment[ic].vz2) * sqrt(densRelErr2 +
                    pow_2(error[ic + (im+1)*fnc.dflen] / result[ic + (im+1)*fnc.dflen]) );
                velocitySecondMomentErr[ic].vphi2 =
                    fabs(velocitySecondMoment[ic].vphi2) * sqrt(densRelErr2 +
                    pow_2(error[ic + (im+2)*fnc.dflen] / result[ic + (im+2)*fnc.dflen]) );
                velocitySecondMomentErr[ic].vRvz =
                    fabs(velocitySecondMoment[ic].vRvz) * sqrt(densRelErr2 +
                    pow_2(error[ic + (im+3)*fnc.dflen] / result[ic + (im+3)*fnc.dflen]) );
                velocitySecondMomentErr[ic].vRvphi =
                    fabs(velocitySecondMoment[ic].vRvphi) * sqrt(densRelErr2 +
                    pow_2(error[ic + (im+4)*fnc.dflen] / result[ic + (im+4)*fnc.dflen]) );
                velocitySecondMomentErr[ic].vzvphi =
                    fabs(velocitySecondMoment[ic].vzvphi) * sqrt(densRelErr2 +
                    pow_2(error[ic + (im+5)*fnc.dflen] / result[ic + (im+5)*fnc.dflen]) );
            }
	    im+=6;
        }
	if(frequencyMoment!=NULL){
		frequencyMoment[ic].Omegar   = densVal ? result[ic + (im+0)*fnc.dflen] / densVal : 0;
		frequencyMoment[ic].Omegaz   = densVal ? result[ic + (im+1)*fnc.dflen] / densVal : 0;
		frequencyMoment[ic].Omegaphi = densVal ? result[ic + (im+2)*fnc.dflen] / densVal : 0;
	}
    }
}


EXP void computeProjectedMoments(const GalaxyModel& model, const double R,
    double* surfaceDensity, double* rmsHeight, double* rmsVel,
    double* surfaceDensityErr, double* rmsHeightErr, double* rmsVelErr,
    const bool separate, const double reqRelError, const int maxNumEval)
{
    double xlower[4] = {0, 0, 0, 0};  // integration region in scaled variables
    double xupper[4] = {1, 1, 1, 1};
    DFIntegrandProjectedMoments fnc(model, separate, R);
    // the values of integrals and their error estimates: allocate a temp array on the stack
    double* result = static_cast<double*>(alloca(fnc.numValues() * 2 * sizeof(double)));
    double* error  = result + fnc.numValues();   // second half of the temporary array
    math::integrateNdim(fnc, xlower, xupper, reqRelError, maxNumEval, result, error);
    for(unsigned int ic=0; ic<fnc.dflen; ic++) {
        if(surfaceDensity)
            surfaceDensity[ic] = result[ic];
        if(rmsHeight)
            rmsHeight[ic] = result[ic]>0 ? sqrt(result[ic + fnc.dflen  ] / result[ic]) : 0;
        if(rmsVel)
            rmsVel[ic]    = result[ic]>0 ? sqrt(result[ic + fnc.dflen*2] / result[ic]) : 0;
        if(surfaceDensityErr)
            surfaceDensityErr[ic] = error[ic];
        if(rmsHeightErr)
            rmsHeightErr[ic] = result[ic]>0 ?
                sqrt(pow_2(error[ic] / result[ic] * result[ic + fnc.dflen]) +
                pow_2(error[ic + fnc.dflen])) : 0;
        if(rmsVelErr)
            rmsVelErr[ic] = result[ic]>0 ?
                sqrt(pow_2(error[ic] / result[ic] * result[ic + fnc.dflen*2]) +
                pow_2(error[ic + fnc.dflen*2])) : 0;
    }
}


EXP void computeProjectedDF(const GalaxyModel& model,
    const double R, const double vz, double result[], const double vz_error,
    bool separate, const double reqRelError, const int maxNumEval)
{
    double xlower[4] = {0, 0, 0, 0};  // integration region in scaled variables
    double xupper[4] = {1, 1, 1, 1};
    DFIntegrandProjected fnc(model, separate, R, vz, vz_error);
    if(vz_error==0) {
        // in this case we may put a tighter upper limit on the integration interval in z,
        // confining it to the region where v^2-vz^2>0 (lower limit is always at z=0)
        xupper[0] = math::findRoot(fnc, 0, 1, 1e-8);
    }
    math::integrateNdim(fnc, xlower, xupper, reqRelError, maxNumEval, result);
}


template <int N>
EXP void computeVelocityDistribution(const GalaxyModel& model,
				     const coord::PosCyl& point, bool projected,
				     const std::vector<double>& gridVR,
				     const std::vector<double>& gridVz,
				     const std::vector<double>& gridVphi,
    /*output arrays (one per each DF component if separate==true, otherwise just one element*/
				     double density[],
				     std::vector<double> amplVR[],
				     std::vector<double> amplVz[],
				     std::vector<double> amplVphi[],
				     bool separate, double reqRelError, int maxNumEval)
{
	double vesc = sqrt(-2 * model.potential.value(
		projected ? coord::PosCyl(point.R, 0, point.phi) : point) );   // escape velocity
	math::BsplineInterpolator1d<N> bsplVR(gridVR), bsplVz(gridVz), bsplVphi(gridVphi);
	const unsigned int NR = bsplVR.numValues(), Nz = bsplVz.numValues(), Nphi = bsplVphi.numValues();
	DFIntegrandVelDist<N> fnc(model, separate, point, projected, bsplVR, bsplVz, bsplVphi);

    // the integration region [scaled z, 3 components of scaled velocity]
	double xlower[4] = {0, 0, 0, 0};
	double xupper[4] = {1, 1, 1, 0.5};  // scan only half of {v_R,v_z} plane, since the VDF is symmetric
    // the values of integrals: allocate a temp array on the stack
	double* result = static_cast<double*>(alloca(fnc.numValues() * sizeof(double)));

	math::integrateNdim(fnc,
			    projected? xlower : xlower+1,   // the 0th dimension (z) only used in the case of projected VDF,
			    projected? xupper : xupper+1,   // otherwise only three components of scaled velocity
			    reqRelError, maxNumEval, result);

    // loop over elements of a multicomponent DF
	for(unsigned int ic=0; ic<fnc.dflen; ic++) {
	// beginning of storage for the current element
		const double* begin = result + fnc.numValues() / fnc.dflen * ic;
	// compute the amplitudes of un-normalized VDF
		amplVR[ic]   = solveForAmplitudes<N>(bsplVR,
			std::vector<double>(begin+1, begin+1+NR), vesc);
		amplVz[ic]   = solveForAmplitudes<N>(bsplVz,
			std::vector<double>(begin+1+NR, begin+1+NR+Nz), vesc);
		amplVphi[ic] = solveForAmplitudes<N>(bsplVphi,
			std::vector<double>(begin+1+NR+Nz, begin+1+NR+Nz+Nphi), vesc);

	// normalize by the value of density
		density[ic] = *begin;
		if(density[ic]!=0) {
			math::blas_dmul(1/density[ic], amplVR[ic]);
			math::blas_dmul(1/density[ic], amplVz[ic]);
			math::blas_dmul(1/density[ic], amplVphi[ic]);
		}
	}
}

// force compilation of several template instantiations
template<>
EXP void computeVelocityDistribution<0>(const GalaxyModel&, const coord::PosCyl&, bool,
    const std::vector<double>&, const std::vector<double>&, const std::vector<double>&,
    double[], std::vector<double>[], std::vector<double>[], std::vector<double>[], bool, double, int);
template<>
EXP void computeVelocityDistribution<1>(const GalaxyModel&, const coord::PosCyl&, bool,
    const std::vector<double>&, const std::vector<double>&, const std::vector<double>&,
    double[], std::vector<double>[], std::vector<double>[], std::vector<double>[], bool, double, int);
template<>
EXP void computeVelocityDistribution<3>(const GalaxyModel&, const coord::PosCyl&, bool,
    const std::vector<double>&, const std::vector<double>&, const std::vector<double>&,
    double[], std::vector<double>[], std::vector<double>[], std::vector<double>[], bool, double, int);

EXP void computeVelocityDistributionO3(const GalaxyModel& model,
				     const coord::PosCyl& point, bool projected,
				     const std::vector<double>& gridVR,
				     const std::vector<double>& gridVz,
				     const std::vector<double>& gridVphi,
    /*output arrays (one per each DF component if separate==true, otherwise just one element*/
				     double density[],
				     std::vector<double> amplVR[],
				     std::vector<double> amplVz[],
				     std::vector<double> amplVphi[],
				     bool separate, double reqRelError, int maxNumEval)
{
	double vesc = sqrt(-2 * model.potential.value(
		projected ? coord::PosCyl(point.R, 0, point.phi) : point) );   // escape velocity
	math::BsplineInterpolator1d<3> bsplVR(gridVR), bsplVz(gridVz), bsplVphi(gridVphi);
	const unsigned int NR = bsplVR.numValues(), Nz = bsplVz.numValues(), Nphi = bsplVphi.numValues();
	DFIntegrandVelDist<3> fnc(model, separate, point, projected, bsplVR, bsplVz, bsplVphi);

    // the integration region [scaled z, 3 components of scaled velocity]
	double xlower[4] = {0, 0, 0, 0};
	double xupper[4] = {1, 1, 1, 0.5};  // scan only half of {v_R,v_z} plane, since the VDF is symmetric
    // the values of integrals: allocate a temp array on the stack
	double* result = static_cast<double*>(alloca(fnc.numValues() * sizeof(double)));

	math::integrateNdim(fnc,
			    projected? xlower : xlower+1,   // the 0th dimension (z) only used in the case of projected VDF,
			    projected? xupper : xupper+1,   // otherwise only three components of scaled velocity
			    reqRelError, maxNumEval, result);

    // loop over elements of a multicomponent DF
	for(unsigned int ic=0; ic<fnc.dflen; ic++) {
	// beginning of storage for the current element
		const double* begin = result + fnc.numValues() / fnc.dflen * ic;
	// compute the amplitudes of un-normalized VDF
		amplVR[ic]   = solveForAmplitudes<3>(bsplVR,
			std::vector<double>(begin+1, begin+1+NR), vesc);
		amplVz[ic]   = solveForAmplitudes<3>(bsplVz,
			std::vector<double>(begin+1+NR, begin+1+NR+Nz), vesc);
		amplVphi[ic] = solveForAmplitudes<3>(bsplVphi,
			std::vector<double>(begin+1+NR+Nz, begin+1+NR+Nz+Nphi), vesc);

	// normalize by the value of density
		density[ic] = *begin;
		if(density[ic]!=0) {
			math::blas_dmul(1/density[ic], amplVR[ic]);
			math::blas_dmul(1/density[ic], amplVz[ic]);
			math::blas_dmul(1/density[ic], amplVphi[ic]);
		}
	}
}

EXP void computeProjection(const GalaxyModel& model,
    const math::IFunctionNdim& spatialSelection,
    const double Xlim[2], const double Ylim[2],
    const double transformMatrix[9],
    double* result, double* error,
    const double reqRelError, const int maxNumEval)
{
    const double xlower[6] = {Xlim[0], Ylim[0], 0, 0, 0, 0};
    const double xupper[6] = {Xlim[1], Ylim[1], 1, 1, 1, 1};
    DFIntegrandProjection fnc(model, spatialSelection, transformMatrix);
    math::integrateNdim(fnc, xlower, xupper, reqRelError, maxNumEval, result, error);
}


EXP particles::ParticleArrayCyl sampleActions(
    const GalaxyModel& model, const size_t nSamp, std::vector<actions::Actions>* actsOutput)
{
    // first sample points from the action space:
    // we use nAct << nSamp  distinct values for actions, and construct tori for these actions;
    // then each torus is sampled with nAng = nSamp/nAct  distinct values of angles,
    // and the action/angles are converted to position/velocity points
    size_t nAng = std::min<size_t>(nSamp/100+1, 16);   // number of sample angles per torus
    size_t nAct = nSamp / nAng + 1;

    // do the sampling in actions space
    double totalMass, totalMassErr;
    std::vector<actions::Actions> actions = df::sampleActions(
        model.distrFunc, nAct, &totalMass, &totalMassErr);
    assert(nAct == actions.size());
    double pointMass = totalMass / (nAct*nAng);

    // next sample angles from each torus
    particles::ParticleArrayCyl points;
    if(actsOutput!=NULL)
        actsOutput->clear();
    for(size_t t=0; t<nAct && points.size()<nSamp; t++) {
        actions::ActionMapperTorus torus(model.potential, actions[t]);
        for(size_t a=0; a<nAng; a++) {
            actions::Angles ang;
            ang.thetar   = 2*M_PI*math::random();
            ang.thetaz   = 2*M_PI*math::random();
            ang.thetaphi = 2*M_PI*math::random();
            points.add(torus.map(actions::ActionAngles(actions[t], ang)), pointMass);
            if(actsOutput!=NULL)
                actsOutput->push_back(actions[t]);
        }
    }
    return points;
}


EXP particles::ParticleArrayCyl samplePosVel(
					     const GalaxyModel& model, const size_t numSamples)
{
	DFIntegrand6dim fnc(model);
	math::Matrix<double> result;      // sampled scaled coordinates/velocities
	double totalMass, errorMass;      // total normalization of the DF and its estimated error
	double xlower[6] = {0,0,0,0,0,0}; // boundaries of sampling region in scaled coordinates
	double xupper[6] = {1,1,1,1,1,1};
	math::sampleNdim(fnc, xlower, xupper, numSamples, result, NULL, &totalMass, &errorMass);
	const double pointMass = totalMass / result.rows();
	particles::ParticleArrayCyl points;
	points.data.reserve(result.rows());
	for(size_t i=0; i<result.rows(); i++) {
		double scaledvars[6] = {result(i,0), result(i,1), result(i,2),
		result(i,3), result(i,4), result(i,5)};
	// transform from scaled vars (array of 6 numbers) to real pos/vel
		points.add(fnc.unscaleVars(scaledvars), pointMass);
	}
	return points;
}

EXP particles::ParticleArrayCyl samplePosVelLF(
					     const GalaxyModel& model, const size_t numSamples)
{
	DF_LFIntegrand6dim fnc(model);
	math::Matrix<double> result;      // sampled scaled coordinates/velocities
	double totalMass, errorMass;      // total normalization of the DF and its estimated error
	double xlower[7] = {0,0,0,0,0,0,0}; // boundaries of sampling region in scaled coordinates
	double xupper[7] = {1,1,1,1,1,1,1};
	math::sampleNdim(fnc, xlower, xupper, numSamples, result, NULL, &totalMass, &errorMass);
	const double pointMass = totalMass / result.rows();
	particles::ParticleArrayCyl points;
	points.data.reserve(result.rows());
	for(size_t i=0; i<result.rows(); i++) {
		double scaledvars[7] = {result(i,0), result(i,1), result(i,2),
		result(i,3), result(i,4), result(i,5), result(i,6)};
		coord::PosVelCylMag star(fnc.unscaleVars(scaledvars));
		points.add(star.first, star.second);
	}
	return points;
}

EXP particles::ParticleArrayCyl samplePosVel_LF(
		const GalaxyModel& model, const size_t numSamples, const double Bright, const double Faint)
{
	particles::ParticleArrayCyl points = samplePosVel(model, numSamples);
	for(unsigned int i=0; i<points.size(); i++)
		points[i].second = model.distrFunc.selectMag(Bright,Faint);
	return points;
}
	
EXP particles::ParticleArray<coord::PosCyl> sampleDensity(
    const potential::BaseDensity& dens, const size_t numPoints)
{
    if(!isFinite(dens.totalMass()))   // safety precautions
        throw std::runtime_error("sampleDensity: model has infinite mass");
    potential::DensityIntegrandNdim fnc(dens, /*require the values of density to be non-negative*/true);
    math::Matrix<double> result;      // sampled scaled coordinates
    double totalMass, errorMass;      // total mass and its estimated error
    double xlower[3] = {0,0,0};       // boundaries of sampling region in scaled coordinates
    double xupper[3] = {1,1,1};
    math::sampleNdim(fnc, xlower, xupper, numPoints, result, NULL, &totalMass, &errorMass);
    const double pointMass = totalMass / result.rows();
    particles::ParticleArray<coord::PosCyl> points;
    points.data.reserve(result.rows());
    for(size_t i=0; i<result.rows(); i++) {
        // if the system is axisymmetric, phi is not provided by the sampling routine
        double scaledvars[3] = {result(i,0), result(i,1), 
            fnc.axisym ? math::random() : result(i,2)};
        // transform from scaled coordinates to the real ones, and store the point into the array
        points.add(fnc.unscaleVars(scaledvars), pointMass);
    }
    return points;
}

EXP std::vector<coord::VelCyl> sampleVelocity(const GalaxyModel& model,
					      const coord::PosCyl loc, const size_t numSamples,
					      obs::BaseLos* Los, const double bright, const double faint
					      )
{
	double *velocityMoment=NULL;
	OperationMode mode = static_cast<OperationMode>(0);
	DFIntegrandAtPoint fnc(model, false, loc, mode, Los, bright, faint);
	math::Matrix<double> result;      // sampled scaled velocities
	double totalDens, errorDens;      // total normalization of the DF and its estimated error
	double xlower[3] = {0,0,0}; // boundaries of sampling region in scaled coordinates
	double xupper[3] = {1,1,1};
	math::sampleNdim(fnc, xlower, xupper, numSamples, result, NULL, &totalDens, &errorDens);
	const double pointDens = totalDens / result.rows();
	std::vector<coord::VelCyl> vels;
	vels.reserve(result.rows());
	for(size_t i=0; i<result.rows(); i++) {
		double scaledvars[3] = {result(i,0), result(i,1), result(i,2)};
	// transform from scaled vars (array of 3 numbers) to real pos/vel
		vels.push_back(fnc.unscaleVars(scaledvars));
	}
	return vels;
}

EXP std::vector<coord::VelCyl> sampleHalfVelocity(const GalaxyModel& model,
					      const coord::PosCyl loc, const size_t numSamples,
					      obs::BaseLos* Los, const double bright, const double faint
					      )
{
	double *velocityMoment=NULL;
	OperationMode mode = static_cast<OperationMode>(0);
	DFIntegrandAtPoint fnc(model, false, loc, mode, Los, bright, faint);
	math::Matrix<double> result;      // sampled scaled velocities
	double totalDens, errorDens;      // total normalization of the DF and its estimated error
	double xlower[3] = {0,0,0}; // boundaries of sampling region in scaled coordinates
	double xupper[3] = {1,1/sqrt(2),1};
	math::sampleNdim(fnc, xlower, xupper, numSamples, result, NULL, &totalDens, &errorDens);
	const double pointDens = totalDens / result.rows();
	std::vector<coord::VelCyl> vels;
	vels.reserve(result.rows());
	for(size_t i=0; i<result.rows(); i++) {
		double scaledvars[3] = {result(i,0), result(i,1), result(i,2)};
	// transform from scaled vars (array of 3 numbers) to real pos/vel
		vels.push_back(fnc.unscaleVars(scaledvars));
	}
	return vels;
}

EXP std::vector<coord::VelCylMag> sampleVelocityLF(const GalaxyModel& model,
					      const coord::PosCyl loc, const size_t numSamples,
					      obs::BaseLos* Los, const double bright, const double faint
					      )
{
	DF_LFIntegrandAtPoint fnc(model, loc, Los, bright, faint);
	math::Matrix<double> result;      // sampled scaled velocities
	double totalDens, errorDens;      // total normalization of the DF and its estimated error
	double xlower[4] = {0,0,0,0}; // boundaries of sampling region in scaled coordinates
	double xupper[4] = {1,1,1,1};
	math::sampleNdim(fnc, xlower, xupper, numSamples, result, NULL, &totalDens, &errorDens);
	const double pointDens = totalDens / result.rows();
	std::vector<coord::VelCylMag> vels;
	vels.reserve(result.rows());
	for(size_t i=0; i<result.rows(); i++) {
		double scaledvars[4] = {result(i,0), result(i,1), result(i,2), result(i,3)};
		coord::PosVelCylMag pvM(fnc.unscaleVars(scaledvars));
		coord::VelCylMag vM(coord::VelCyl(pvM.first),pvM.second);
		vels.push_back(vM);
	}
	return vels;
}

EXP std::vector<coord::VelCylMag> sampleVelocity_LF(const GalaxyModel& model,
					      const coord::PosCyl loc, const size_t numSamples,
					      obs::BaseLos* Los, const double bright, const double faint
					      )
{
	std::vector<coord::VelCyl> vels=sampleVelocity(model, loc, numSamples, Los, bright, faint);
	double s = Los->s(loc);
	double shift = 5*log10(100 * s) + Los->A_H(s);
	double Bright = bright - shift, Faint = faint - shift;//Now abs mags
	std::vector<coord::VelCylMag> stars;
	stars.reserve(vels.size());
	for(size_t i=0; i<vels.size(); i++) {
		//mag chosen from LF of component[0]
		double m = model.distrFunc.selectMag(Bright, Faint) + shift;//apparent mag
		coord::VelCylMag VM(vels[i],m);
		stars.push_back(VM);
	}
	return stars;
}

// Sample stars seen along a line of sight 
EXP std::vector<coord::PosVelCyl> sampleLOS(const GalaxyModel& model,
					    const obs::BaseLos* Los, const size_t numSamples,
					    const double bright, const double faint)
{
	double *velocityMoment=NULL;
	OperationMode mode = static_cast<OperationMode>(0);
	DFIntegrandLOS fnc(model, false, Los, mode, bright, faint);
	math::Matrix<double> result;  // sampled scaled distances & velocities
	double totalDens, errorDens;  // total normalization of the DF and its estimated error
	double xlower[4] = {0,0,0,0}; // boundaries of sampling region in scaled coordinates		
	double xupper[4] = {1,1,1,1};
	math::sampleNdim(fnc, xlower, xupper, numSamples, result, NULL, &totalDens, &errorDens);
	const double pointDens = totalDens / result.rows();
	std::vector<coord::PosVelCyl> posVels;						//change to PosVelCyl
	posVels.reserve(result.rows());
	for(size_t i=0; i<result.rows(); i++) {
		double scaledvars[4] = {result(i,0), result(i,1), result(i,2), result(i,3)};
	// transform from scaled vars (array of 4 numbers) to real pos/vel
		posVels.push_back(fnc.unscaleVars(scaledvars));
	}														
	return posVels;
}

// Sample stars seen along a line of sight returning only (s,Vlos)
EXP std::vector<std::pair<double,double> > sampleLOSsVlos(const GalaxyModel& model,
					    const obs::BaseLos* Los, const size_t numSamples,
					    const double bright, const double faint)
{
	double *velocityMoment=NULL;
	OperationMode mode = static_cast<OperationMode>(0);
	DFIntegrandLOS fnc(model, false, Los, mode, bright, faint);
	math::Matrix<double> result;  // sampled scaled distances & velocities
	double totalDens, errorDens;  // total normalization of the DF and its estimated error
	double xlower[4] = {0,0,0,0}; // boundaries of sampling region in scaled coordinates		
	double xupper[4] = {1,1,1,1};
	math::sampleNdim(fnc, xlower, xupper, numSamples, result, NULL, &totalDens, &errorDens);
	const double pointDens = totalDens / result.rows();
	std::vector<std::pair<double,double> >  sVs;
	sVs.reserve(result.rows());
	for(size_t i=0; i<result.rows(); i++) {
		double scaledvars[4] = {result(i,0), result(i,1), result(i,2), result(i,3)};
	// transform from scaled vars (array of 4 numbers) to real pos/vel
		sVs.push_back(Los->sVlos(fnc.unscaleVars(scaledvars)));
	}														
	return sVs;
}

EXP std::vector<coord::PosVelCylMag> sampleLOSLF(const GalaxyModel& model,
					    const obs::BaseLos* Los, const size_t numSamples,
					    const double bright, const double faint)
{
	DF_LFIntegrandLOS fnc(model, Los, bright, faint);
	math::Matrix<double> result;  // sampled scaled distances & velocities
	double totalDens, errorDens;  // total normalization of the DF and its estimated error
	double xlower[5] = {0,0,0,0,0}; // boundaries of sampling region in scaled coordinates		
	double xupper[5] = {1,1,1,1,1};
	math::sampleNdim(fnc, xlower, xupper, numSamples, result, NULL, &totalDens, &errorDens);
	const double pointDens = totalDens / result.rows();
	std::vector<coord::PosVelCylMag> posVels;						//change to PosVelCyl
	posVels.reserve(result.rows());
	for(size_t i=0; i<result.rows(); i++) {
		double scaledvars[5] = {result(i,0), result(i,1), result(i,2), result(i,3), result(i,4)};
	// transform from scaled vars (array of 4 numbers) to real pos/vel
		posVels.push_back(fnc.unscaleVars(scaledvars));
	}														
	return posVels;
}

EXP std::vector<coord::PosVelCylMag> sampleLOS_LF(const GalaxyModel& model,
	const obs::BaseLos* Los, const size_t numSamples,
	const double bright, const double faint)
{
	std::vector<coord::PosVelCyl> points = sampleLOS(model, Los, numSamples, bright, faint);
	std::vector<coord::PosVelCylMag> stars;
	stars.reserve(points.size());
	for(size_t i=0; i<points.size(); i++) {
		double s = Los->s(points[i]);
		double shift = 5*log10(100 * s) + Los->A_H(s);
		double Bright = bright - shift, Faint = faint - shift;//Now abs mags
		//apparent mag chosen from LF of component[0]
		double m = model.distrFunc.selectMag(Bright, Faint) + shift;
		coord::PosVelCylMag pvM(points[i],m);
		stars.push_back(pvM);
	}
	return stars;
}

}  // namespace
