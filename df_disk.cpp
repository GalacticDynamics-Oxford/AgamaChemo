#include "df_disk.h"
#include <cmath>
#include <stdexcept>

namespace df{

namespace{  // internal

/// compute the average of DF over stellar age:
/// ( \int_0^1 dt B^2(t) \exp[ t/t_0 - A*B(t) ] ) / ( \int_0^1 dt \exp[ t/t_0 ] ),
/// where B(t) = ( (t + t_1) / (1 + t_1) )^{-2\beta}
double averageOverAge(double A, const AgeVelocityDispersionParam& par)
{
    if(par.beta == 0 || par.sigmabirth == 1 || !isFinite(par.Tsfr))
        return exp(-A);
    // if we have a non-trivial age-velocity dispersion relation,
    // then we need to integrate over sub-populations convolved with star formation history
    static const int NT = 5;      // number of points in quadrature rule for integration over age
    static const double qx[NT] =  // nodes of quadrature rule
    { 0.04691007703066802, 0.23076534494715845, 0.5, 0.76923465505284155, 0.95308992296933198 };
    static const double qw[NT] =  // weights of quadrature rule
    { 0.11846344252809454, 0.23931433524968324, 64./225, 0.23931433524968324, 0.11846344252809454 };
    double s = std::pow(par.sigmabirth, 1./par.beta); 
    double integ=0, norm = 0;
    for(int i=0; i<NT; i++) {
        // t is the lookback time (stellar age) measured in units of galaxy time (ranges from 0 to 1)
        double t = qx[i];
        // star formation rate exponentially increases with look-back time
        double weight = exp(t / par.Tsfr) * qw[i];
        // velocity dispersions {sigma_r, sigma_z} scale as  [ t + s * (1-t) ]^beta
        double multsq = std::pow(t + (1-t) * s, -2*par.beta);  // multiplied by sigma^-2
        integ += weight * exp(-A * multsq) * pow_2(multsq);
        norm  += weight;
    }
    return integ / norm;
}
}

EXP QuasiIsothermal::QuasiIsothermal(const QuasiIsothermalParam &params, const potential::Interpolator& freqs) :
    par(params), freq(freqs)
{
    // sanity checks on parameters
    if(!(par.Sigma0>0))
        throw std::invalid_argument("QuasiIsothermal: surface density Sigma0 must be positive");
    if(!(par.Rdisk>0))
        throw std::invalid_argument("QuasiIsothermal: disk scale radius Rdisk must be positive");
    if(!(par.sigmar0>0))
        throw std::invalid_argument("QuasiIsothermal: velocity dispersion sigmar0 must be positive");
    if(!(par.Rsigmar>0))
        throw std::invalid_argument("QuasiIsothermal: velocity scale radius Rsigmar must be positive");
    if(!( (par.Hdisk>0) ^ (par.sigmaz0>0 && par.Rsigmaz>0) ))
        throw std::invalid_argument("QuasiIsothermal: should have either "
            "Hdisk>0 to assign the vertical velocity dispersion from disk scaleheight, or "
            "Rsigmaz>0, sigmaz0>0 to make it exponential in radius");
    if(par.Hdisk<0 || par.sigmaz0<0 || par.Rsigmaz<0)  // these are optional but non-negative
        throw std::invalid_argument("QuasiIsothermal: parameters cannot be negative");
    if(par.sigmabirth<=0 || par.sigmabirth>1)
        throw std::invalid_argument("QuasiIsothermal: invalid value for velocity dispersion at birth");
}

EXP double QuasiIsothermal::value(const actions::Actions &J) const
{
    // obtain the radius of in-plane motion with the given "characteristic" angular momentum
    double Rcirc = freq.R_from_Lz(sqrt(pow_2(par.Jmin) +
        pow_2(fabs(J.Jphi) + par.coefJr * J.Jr + par.coefJz * J.Jz)) );
    if(Rcirc > 20 * par.Rdisk)
        return 0;   // we're too far out, DF is negligibly small
    double kappa, nu, Omega;   // characteristic epicyclic freqs
    freq.epicycleFreqs(Rcirc, kappa, nu, Omega);
    // surface density follows an exponential profile in radius
    double Sigma = par.Sigma0 * exp( -Rcirc / par.Rdisk );
    // squared radial velocity dispersion is exponential in radius
    double sigmarsq = pow_2(par.sigmar0 * exp ( -Rcirc / par.Rsigmar ) ) + pow_2(par.sigmamin);
    // squared vertical velocity dispersion computed by either of the two methods: 
    double sigmazsq = pow_2(par.sigmamin) + (par.Hdisk>0 ?
        2 * pow_2(nu * par.Hdisk) :     // keep the disk thickness approximately equal to Hdisk, or
        pow_2(par.sigmaz0 * exp ( -Rcirc / par.Rsigmaz ) ) );  // make sigmaz exponential in radius
    // suppression factor for counterrotating orbits
    double negJphi = J.Jphi>0 ? 0. : 2*Omega * J.Jphi;
    double result = 1./(2*M_PI*M_PI) * Sigma * nu * Omega / (kappa * sigmarsq * sigmazsq) *
        averageOverAge( (kappa * J.Jr - negJphi) / sigmarsq + nu * J.Jz / sigmazsq, par);
    return isFinite(result) ? result : 0;
}


EXP Exponential::Exponential(const ExponentialParam& params) :
    par(params)
{
	if(!(par.norm>0))
		throw std::invalid_argument("Exponential: norm must be positive");
	if(!(par.Jr0>0) || !(par.Jz0>0) || !(par.Jphi0>0))
		throw std::invalid_argument("Exponential: scale actions must be positive");
	if(par.addJden<0 || par.addJden >= par.Jphi0)
		throw std::invalid_argument("Exponential: addJden must be in (0, Jphi0)");
	if(par.Fname!="")
		readBrighterThan(par.Fname);
}

EXP double Exponential::value(const actions::Actions &J) const
{
	if(J.Jphi<=0) return 0;
	double Jt=1.5*J.Jr+J.Jz+J.Jphi;
	double Jvel = Jt + par.addJvel, Jden = Jt + par.addJden;
	double xr = pow(Jvel/par.Jphi0,par.pr)/par.Jr0;
	double xz = pow(Jvel/par.Jphi0,par.pz)/par.Jz0;
	double fr = xr * exp(-xr*J.Jr), fz = xz * exp(-xz*J.Jz);
	double xp = Jden / par.Jphi0;
	double fp = par.norm/par.Jphi0 * J.Jphi / par.Jphi0 * exp(-xp);
	return fr * fz * fp;
}

/*EXP double Exponential::value(const actions::Actions &J) const
{
	double Jp = J.Jphi<=0 ? 0 : J.Jphi;
	if(Jp==0) return 0;
	double Jvel = fabs(Jp) + par.addJvel;
	double xr = pow(Jvel/par.Jphi0,par.pr)/par.Jr0;
	double xz = pow(Jvel/par.Jphi0,par.pz)/par.Jz0;
	double fr = xr * exp(-xr*J.Jr), fz = xz * exp(-xz*J.Jz);
	double Jden = Jp + par.addJden;
	double xp = Jden / par.Jphi0;
	double fp = par.norm/par.Jphi0 * fabs(J.Jphi) / par.Jphi0 * exp(-xp);
	double F = fr * fz * fp;
	if(J.Jphi >= 0)
		return F;
	else{
		double x=J.Jphi/par.addJden;
		return exp(x*(1-x)) * F;
	}
}*/

EXP void Exponential::write_params(std::ofstream &strm,const units::InternalUnits &intUnits) const{
	strm << "type\t = Exponential\n";
	strm << "mass\t = " << par.mass * intUnits.to_Msun << '\n';
	strm << "Jr0\t = "<< par.Jr0 * intUnits.to_Kpc_kms << '\n';
	strm << "Jz0\t = " << par.Jz0 * intUnits.to_Kpc_kms << '\n';
	strm << "Jphi0\t = " << par.Jphi0 * intUnits.to_Kpc_kms << '\n';
	strm << "pr\t = " << par.pr << '\n';
	strm << "pz\t = " << par.pz << '\n';
	strm << "addJden\t = " << par.addJden * intUnits.to_Kpc_kms << '\n';
	strm << "addJvel\t = " << par.addJvel * intUnits.to_Kpc_kms << '\n';
}

taperExp::taperExp(const taperExpParam& params) :
    par(params)
{
	if(!(par.norm>0))
		throw std::invalid_argument("taperExp: overall normalization must be positive");
	if(!(par.Jr0>0) || !(par.Jz0>0) || !(par.Jphi0>0))
		throw std::invalid_argument("taperExp: scale actions must be positive");
	if(par.Fname!="")
		readBrighterThan(par.Fname);
}

double taperExp::value(const actions::Actions &J) const
{
	if(J.Jphi<=0) return 0;
	//double Jt = J.Jr+J.Jz+J.Jphi;
	double Jt = J.Jr + J.Jz + J.Jphi;
	double Jvel = Jt + par.addJvel, Jden = Jt + par.addJden;
	double xr = pow(Jvel/par.Jphi0,par.pr)/par.Jr0;
	double xz = pow(Jvel/par.Jphi0,par.pz)/par.Jz0;
	double fr = xr * exp(-xr*J.Jr), fz = xz * exp(-xz*J.Jz);
	double xp = Jden / par.Jphi0;
	double fp = par.norm/par.Jphi0  * Jden / par.Jphi0 * exp(-xp);
	double taper,cut;
	if(par.Jtrans>0){
		double ep = exp((J.Jphi-par.Jtaper)/par.Jtrans);
		taper = ep/(ep+1/ep);//.5*(1+tanh(ep))
	} else
		taper = 1;
	if(par.Delta>0) {
		double ep = exp(-(J.Jphi-par.Jcut)/par.Delta);
		cut = ep/(ep+1/ep);
	} else
		cut = 1;
//	double taper = .5*fmax(0,1+tanh);
//	double cut   = .5*fmax(0,1-tanh((J.Jphi-par.Jcut)/par.Delta));
	return fr * fz * fp * taper * cut;
}

void taperExp::write_params(std::ofstream &strm,const units::InternalUnits &intUnits) const{
	strm << "type\t = taperExp\n";
	strm << "norm\t = " << par.norm * intUnits.to_Msun << '\n';
	strm << "mass\t = " << par.mass * intUnits.to_Msun << '\n';
	strm << "Jr0\t = " << par.Jr0 * intUnits.to_Kpc_kms << '\n';
	strm << "Jz0\t = " << par.Jz0 * intUnits.to_Kpc_kms << '\n';
	strm << "Jphi0\t = " << par.Jphi0 * intUnits.to_Kpc_kms << '\n';
	if(par.Jtrans>0){
		strm << "Jtaper\t = " << par.Jtaper * intUnits.to_Kpc_kms << '\n';
		strm << "Jtrans\t = " << par.Jtrans * intUnits.to_Kpc_kms << '\n';
	}
	if(par.Delta>0){
		strm << "Jcut\t = " << par.Jcut * intUnits.to_Kpc_kms << '\n';
		strm << "Delta\t = " << par.Delta * intUnits.to_Kpc_kms << '\n';
	}
	strm << "pr\t = " << par.pr << '\n';
	strm << "pz\t = " << par.pz << '\n';
	strm << "addJden\t = " << par.addJden * intUnits.to_Kpc_kms << '\n';
	strm << "addJvel\t = " << par.addJvel * intUnits.to_Kpc_kms << '\n';
}
//void taperExp::tabulate_params(std::ofstream& strm,const::InternalUnits& intUnits) const{
	
/*Exponential::Exponential(const ExponentialParam& params) :
    par(params)
{
    if(!(par.norm>0))
        throw std::invalid_argument("Exponential: overall normalization must be positive");
    if(!(par.Jr0>0) || !(par.Jz0>0) || !(par.Jphi0>0))
        throw std::invalid_argument("Exponential: scale actions must be positive");
    if(par.sigmabirth<=0 || par.sigmabirth>1)
        throw std::invalid_argument("Exponential: invalid value for velocity dispersion at birth");
}

double Exponential::value(const actions::Actions &J) const
{
    // weighted sum of actions
    double Jsum = fabs(J.Jphi) + par.coefJr * J.Jr + par.coefJz * J.Jz;
    double Jden = sqrt(pow_2(Jsum) + pow_2(par.addJden));
    double Jvel = sqrt(pow_2(Jsum) + pow_2(par.addJvel));
    // suppression factor for counterrotating orbits
    double negJphi = J.Jphi>0 ? 0. : J.Jphi;
    return 1. / TWO_PI_CUBE * par.norm / pow_2(par.Jr0 * par.Jz0 * par.Jphi0) *
        Jvel * Jvel * Jden * exp(-Jden / par.Jphi0) *
        averageOverAge(Jvel * ((J.Jr - negJphi) / pow_2(par.Jr0) + J.Jz / pow_2(par.Jz0)), par);
}
*/
}  // namespace df
