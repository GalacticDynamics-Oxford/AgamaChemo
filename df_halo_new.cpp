#include "df_halo.h"
#include "math_core.h"
#include "math_ode.h"
#include <cmath>
#include <stdexcept>

#define DE

namespace df{

namespace {  // internal ns

/// helper class used in the root-finder to determine the auxiliary
/// coefficient beta for a cored halo
/// We ned one for DoublePowerLaw and one for ModifiedDoublePowerLaw
class BetaFinder: public math::IFunctionNoDeriv{
	const DoublePowerLawParam& par;

    // return the difference between the non-modified and modified DF as a function of beta and
    // the appropriately scaled action variable (t -> hJ), weighted by d(hJ)/dt for the integration in t
	double deltaf(const double t, const double beta) const
	{
	// integration is performed in a scaled variable t, ranging from 0 to 1,
	// which is remapped to hJ ranging from 0 to infinity as follows:
		double hJ    = par.Jcore * t*t*(3-2*t) / pow_2(1-t) / (1+2*t); 
		double dhJdt = par.Jcore * 6*t / pow_3(1-t) / pow_2(1+2*t);
		return hJ * hJ * dhJdt *
				math::pow(1 + math::pow(par.J0 / hJ, par.steepness),  par.slopeIn  / par.steepness) *
				math::pow(1 + math::pow(hJ / par.J0, par.steepness), -par.slopeOut / par.steepness) *
				(math::pow(1 + par.Jcore/hJ * (par.Jcore/hJ - beta), -0.5*par.slopeIn) - 1);
	}

	public:
		BetaFinder(const DoublePowerLawParam& _par) : par(_par) {}

		virtual double value(const double beta) const
		{
			double result = 0;
	// use a fixed-order GL quadrature to compute the integrated difference in normalization between
	// unmodified and core-modified DF, which is sought to be zero by choosing an appropriate beta
			static const int GLORDER = 20;  // should be even, to avoid singularity in the integrand at t=0.5
			for(int i=0; i<GLORDER; i++)
				result += math::GLWEIGHTS[GLORDER][i] * deltaf(math::GLPOINTS[GLORDER][i], beta);
			return result;
		}
};

// helper function to compute the auxiliary coefficient beta in the case of a central core
double computeBeta(const DoublePowerLawParam &par)
{
	if(par.Jcore<=0)
		return 0;
	return math::findRoot(BetaFinder(par), 0.0, 2.0, /*root-finder tolerance*/ SQRT_DBL_EPSILON);
}

/// 2nd helper class used in the root-finder to determine the auxiliary coefficient beta for a cored halo
class ModifiedBetaFinder: public math::IFunctionNoDeriv{
	const ModifiedDoublePowerLawParam& par;

    // return the difference between the non-modified and modified DF as a function of beta and
    // the appropriately scaled action variable (t -> hJ), weighted by d(hJ)/dt for the integration in t
	double deltaf(const double t, const double beta) const
	{
	// integration is performed in a scaled variable t, ranging from 0 to 1,
	// which is remapped to hJ ranging from 0 to infinity as follows:
		double hJ    = par.Jcore * t*t*(3-2*t) / pow_2(1-t) / (1+2*t); 
		double dhJdt = par.Jcore * 6*t / pow_3(1-t) / pow_2(1+2*t);
		return hJ * hJ * dhJdt *
				math::pow(1 + par.J0 / hJ,  par.slopeIn) *
				math::pow(1 + hJ / par.J0, -par.slopeOut) *
				(math::pow(1 + par.Jcore/hJ * (par.Jcore/hJ - beta), -0.5*par.slopeIn) - 1);
	}

	public:
		ModifiedBetaFinder(const ModifiedDoublePowerLawParam& _par) : par(_par) {}

		virtual double value(const double beta) const
		{
			double result = 0;
	// use a fixed-order GL quadrature to compute the integrated difference in normalization between
	// unmodified and core-modified DF, which is sought to be zero by choosing an appropriate beta
			static const int GLORDER = 20;  // should be even, to avoid singularity in the integrand at t=0.5
			for(int i=0; i<GLORDER; i++)
				result += math::GLWEIGHTS[GLORDER][i] * deltaf(math::GLPOINTS[GLORDER][i], beta);
			return result;
		}
};

// helper function to compute the auxiliary coefficient beta in the case of a central core
double computeBeta(const ModifiedDoublePowerLawParam &par)
{
	if(par.Jcore<=0)
		return 0;
	return math::findRoot(ModifiedBetaFinder(par), 0.0, 2.0, /*root-finder tolerance*/ SQRT_DBL_EPSILON);
}

/// 3rd helper class used in the root-finder to determine the auxiliary coefficient beta for a cored halo
class SinBetaFinder: public math::IFunctionNoDeriv{
	const SinDoublePowerLawParam& par;

    // return the difference between the non-modified and modified DF as a function of beta and
    // the appropriately scaled action variable (t -> hJ), weighted by d(hJ)/dt for the integration in t
	double deltaf(const double t, const double beta) const
	{
	// integration is performed in a scaled variable t, ranging from 0 to 1,
	// which is remapped to hJ ranging from 0 to infinity as follows:
		double hJ    = par.Jcore * t*t*(3-2*t) / pow_2(1-t) / (1+2*t); 
		double dhJdt = par.Jcore * 6*t / pow_3(1-t) / pow_2(1+2*t);
		return hJ * hJ * dhJdt *
				math::pow(1 + par.J0 / hJ,  par.slopeIn) *
				math::pow(1 + hJ / par.J0, -par.slopeOut) *
				(math::pow(1 + par.Jcore/hJ * (par.Jcore/hJ - beta), -0.5*par.slopeIn) - 1);
	}

	public:
		SinBetaFinder(const SinDoublePowerLawParam& _par) : par(_par) {}

		virtual double value(const double beta) const
		{
			double result = 0;
	// use a fixed-order GL quadrature to compute the integrated difference in normalization between
	// unmodified and core-modified DF, which is sought to be zero by choosing an appropriate beta
			static const int GLORDER = 20;  // should be even, to avoid singularity in the integrand at t=0.5
			for(int i=0; i<GLORDER; i++)
				result += math::GLWEIGHTS[GLORDER][i] * deltaf(math::GLPOINTS[GLORDER][i], beta);
			return result;
		}

/*		double g(double x,double y){
			return (1+y*x*x)/(1+x*x);
		}*/
};

// helper function to compute the auxiliary coefficient beta in the case of a central core
double computeBeta(const SinDoublePowerLawParam &par)
{
	if(par.Jcore<=0)
		return 0;
	return math::findRoot(SinBetaFinder(par), 0.0, 2.0, /*root-finder tolerance*/ SQRT_DBL_EPSILON);
}

}  // internal ns

EXP DoublePowerLaw::DoublePowerLaw(const DoublePowerLawParam &inparams) :
    par(inparams), beta(computeBeta(par))
{
    // sanity checks on parameters
	if(!(par.norm>0))
		throw std::invalid_argument("DoublePowerLaw: normalization must be positive");
	if(!(par.J0>0))
		throw std::invalid_argument("DoublePowerLaw: break action J0 must be positive");
	if(!(par.Jcore>=0 && beta>=0))
		throw std::invalid_argument("DoublePowerLaw: core action Jcore is invalid");
	if(!(par.Jcutoff>=0))
		throw std::invalid_argument("DoublePowerLaw: truncation action Jcutoff must be non-negative");
	if(!(par.slopeOut>3) && par.Jcutoff==0)
		throw std::invalid_argument(
					    "DoublePowerLaw: mass diverges at large J (outer slope must be > 3)");
	if(!(par.slopeIn<3))
		throw std::invalid_argument(
					    "DoublePowerLaw: mass diverges at J->0 (inner slope must be < 3)");
	if(!(par.steepness>0))
		throw std::invalid_argument("DoublePowerLaw: transition steepness parameter must be positive");
	if(!(par.cutoffStrength>0))
		throw std::invalid_argument("DoublePowerLaw: cutoff strength parameter must be positive");
	if(!(par.coefJrIn>0 && par.coefJzIn >0 && par.coefJrIn +par.coefJzIn <3 &&
	     par.coefJrOut>0 && par.coefJzOut>0 && par.coefJrOut+par.coefJzOut<3) )
		throw std::invalid_argument(
					    "DoublePowerLaw: invalid weights in the linear combination of actions");
	if(!(fabs(par.rotFrac)<=1))
		throw std::invalid_argument(
					    "DoublePowerLaw: amplitude of odd-Jphi component must be between -1 and 1");
	if(par.Fname!=""){
		readBrighterThan(par.Fname);
	}
}

EXP double DoublePowerLaw::value(const actions::Actions &J) const
{
    // linear combination of actions in the inner part of the model (for J<J0)
	double hJ  = par.coefJrIn * J.Jr + par.coefJzIn * J.Jz +
		     (3-par.coefJrIn -par.coefJzIn) * fabs(J.Jphi);
    // linear combination of actions in the outer part of the model (for J>J0)
	double gJ  = par.coefJrOut* J.Jr + par.coefJzOut* J.Jz +
		     (3-par.coefJrOut-par.coefJzOut)* fabs(J.Jphi);
	double val = par.norm / pow_3(2*M_PI * par.J0) *
		     math::pow(1 + math::pow(par.J0 / hJ, par.steepness),  par.slopeIn  / par.steepness) *
		     math::pow(1 + math::pow(gJ / par.J0, par.steepness), -par.slopeOut / par.steepness);
	if(par.rotFrac!=0)  // add the odd part
		val *= 1 + par.rotFrac * tanh(J.Jphi / par.Jphi0);
	if(par.Jcutoff>0)   // exponential cutoff at large J
		val *= exp(-math::pow(gJ / par.Jcutoff, par.cutoffStrength));
	if(par.Jcore>0) {   // central core of nearly-constant f(J) at small J
		if(hJ==0) return par.norm / pow_3(2*M_PI * par.J0);
		val *= math::pow(1 + par.Jcore/hJ * (par.Jcore/hJ - beta), -0.5*par.slopeIn);
	}
	return val;
}

EXP ModifiedDoublePowerLaw::ModifiedDoublePowerLaw(const ModifiedDoublePowerLawParam &inparams) :
    par(inparams), beta(computeBeta(par))
{
    // sanity checks on parameters
	if(!(par.norm>0))
		throw std::invalid_argument("ModifiedDoublePowerLaw: normalization must be positive");
	if(!(par.J0>0))
		throw std::invalid_argument("ModifiedDoublePowerLaw: break action J0 must be positive");
	if(!(par.Jcore>=0 && beta>=0))
		throw std::invalid_argument("ModifiedDoublePowerLaw: core action Jcore is invalid");
	if(!(par.Jcutoff>=0))
		throw std::invalid_argument("ModifiedDoublePowerLaw: truncation action Jcutoff must be non-negative");
	if(!(par.slopeOut>3) && par.Jcutoff==0)
		throw std::invalid_argument(
					    "ModifiedDoublePowerLaw: mass diverges at large J (outer slope must be > 3)");
	if(!(par.slopeIn<3))
		throw std::invalid_argument(
					    "ModifiedDoublePowerLaw: mass diverges at J->0 (inner slope must be < 3)");
	if(!(par.cutoffStrength>0))
		throw std::invalid_argument("ModifiedDoublePowerLaw: cutoff strength parameter must be positive");
	if(!(fabs(par.rotFrac)<=1))
		throw std::invalid_argument(
					    "ModifiedDoublePowerLaw: amplitude of odd-Jphi component must be between -1 and 1");
	if(par.Fname!=""){
		readBrighterThan(par.Fname);
	}
}

EXP double ModifiedDoublePowerLaw::value(const actions::Actions &J) const{
	double Jp=fabs(J.Jphi), L=Jp+J.Jz, L1=0;
	double b=.5*(par.coefJzIn+1), c=.5*(par.coefJzIn-1);
	double cL=Jp+b*J.Jz+c*J.Jz*Jp/(L+L1);
	double h=par.coefJrIn*J.Jr+.5*cL;
	double f = par.norm / pow_3(2*M_PI * par.J0) *
		   math::pow(1 + par.J0 / h,  par.slopeIn) *
		   math::pow(1 + h / par.J0, -par.slopeOut);
	if(par.Jcutoff>0)   // exponential cutoff at large J
		f *= exp(-math::pow(h / par.Jcutoff, par.cutoffStrength));
	if(par.Jcore>0) {   // central core of nearly-constant f(J) at small J
		if(h==0) return par.norm / pow_3(2*M_PI * par.J0);
		f *= math::pow(1 + par.Jcore/h * (par.Jcore/h - beta), -0.5*par.slopeIn);
	}
	if(par.rotFrac!=0)  // add the odd part
		f *= 1 + par.rotFrac * tanh(J.Jphi / par.Jphi0);
	return f;
}

EXP SinDoublePowerLaw::SinDoublePowerLaw(const SinDoublePowerLawParam &inparams) :
    par(inparams), beta(computeBeta(par))
{
    // sanity checks on parameters
	if(!(par.norm>0))
		printf("SinDoublePowerLaw: normalization must be positive");
	if(!(par.J0>0))
		printf("SinDoublePowerLaw: break action J0 must be positive");
	if(!(par.Jcore>=0 && beta>=0))
		printf("SinDoublePowerLaw: core action Jcore is invalid");
	if(!(par.Jcutoff>=0))
		printf("SinDoublePowerLaw: truncation action Jcutoff must be non-negative");
	if(!(par.slopeOut>3) && par.Jcutoff==0)
		printf("SinDoublePowerLaw: mass diverges at large J (outer slope must be > 3)");
	if(!(par.slopeIn<3))
		printf("SinDoublePowerLaw: mass diverges at J->0 (inner slope must be < 3)");
	if(!(par.cutoffStrength>0))
		printf("SinDoublePowerLaw: cutoff strength parameter must be positive");
	if(!(fabs(par.rotFrac)<=1))
		printf("SinDoublePowerLaw: amplitude of odd-Jphi component must be between -1 and 1");
	if(par.Fname!=""){
		readBrighterThan(par.Fname);
	}
}
EXP double SinDoublePowerLaw::value(const actions::Actions &J) const
{
	double modJphi=fabs(J.Jphi);
	double L=J.Jz+modJphi;
	double c=L/(L+J.Jr);
	double jt=(1.5*J.Jr+L)/par.L0;
	double jta=pow(jt,par.alpha), xi=jta/(1+jta);
	double rat=(1-xi)*par.Fin+xi*par.Fout;
	double a=.5*(rat+1), b=.5*(rat-1);
	double cL= L>0? J.Jz*(a+b*modJphi/L)+modJphi: 0;
	double fac=exp(par.beta*sin(0.5*M_PI*c));
	double hJ=J.Jr/fac + .5*(1+c*xi)*fac*cL;
	double gJ=hJ;
	double val = par.norm / pow_3(2*M_PI * par.J0) *
		     math::pow(1 + par.J0 / hJ,  par.slopeIn) *
		     math::pow(1 + gJ / par.J0, -par.slopeOut);
//	printf("%g %g %g |",J.Jz,L,cL);
	if(par.Jcutoff>0){   // exponential cutoff at large J
		double fac=math::pow(gJ / par.Jcutoff, par.cutoffStrength);
		if(fac>25) return 0;
		else val *= exp(-fac);
	}
	if(par.Jcore>0) {   // central core of nearly-constant f(J) at small J
		if(hJ==0) return par.norm / pow_3(2*M_PI * par.J0);
		val *= math::pow(1 + par.Jcore/hJ * (par.Jcore/hJ - beta), -0.5*par.slopeIn);
	}
	if(par.rotFrac!=0)  // add the odd part
		val *= 1 + par.rotFrac * tanh(J.Jphi / par.Jphi0);
	if(std::isnan(val)){
		printf("SinDoublePowerLaw: val=NAN, J=(%g %g %g), hJ=%g\n",J.Jr,J.Jz,J.Jphi,hJ);
		printf("%g %g %g %g %g\n",cL,fac,hJ,par.norm,par.Jcore);
		printf("%g %g %g\n",par.J0,par.slopeIn,par.slopeOut); exit(1);
	}
	return val;
}

EXP void SinDoublePowerLaw::write_params(std::ofstream &strm,const units::InternalUnits &intUnits) const{
	strm << "type\t = SinDoublePowerLaw\n";
	strm << "mass\t = " << par.mass * intUnits.to_Msun << '\n';
	strm << "J0\t = " << par.J0 * intUnits.to_Kpc_kms << '\n';
	strm << "L0\t = " << par.L0 * intUnits.to_Kpc_kms << '\n';
	strm << "Jcore\t = " << par.Jcore * intUnits.to_Kpc_kms << '\n';
	strm << "Jcutoff\t = " << par.Jcutoff * intUnits.to_Kpc_kms << '\n';
	strm << "cutoffStrength\t = " << par.cutoffStrength << '\n';
	strm << "slopeIn\t = " << par.slopeIn << '\n';
	strm << "slopeOut\t = " << par.slopeOut << '\n';
	strm << "Fin\t = " << par.Fin <<'\n';
	strm << "Fout\t = " << par.Fout <<'\n';
	strm << "alpha\t = " << par.alpha << '\n';
	strm << "beta\t = " << par.beta << '\n';
}

/*
EXP double ModifiedDoublePowerLaw::value(const actions::Actions &J) const{
	const double pisq=pow_2(M_PI);
	double modJphi = fabs(J.Jphi);
	double L = J.Jz + modJphi;
	double c = L/(L+J.Jr), p = modJphi/L;//define circularity & planarity
	double cJ = J.Jr + .5*L;
	double jt = pow_2(cJ)/(pow_2(cJ)+pow_2(par.L0));
//	double jt = pow_2(L)/(pow_2(L)+pow_2(par.L0));
	double sinr = pow_2(sin(c*M_PI*.5));
	double facp = (1 + par.Fout*sinr*pow_2(cos(p*M_PI*.5)))/(1+.5*par.Fout*(.5+2/pisq));//Fout biases to planar orbits
//	double facr=par.Fin*pow_2(cos(.5*M_PI*modJphi/(L+J.Jr)));
//	double h = cJ/facr*facp;
//	double facr=pow_2(cos(.5*M_PI*modJphi/(L+J.Jr)));
//	double h = cJ*(1+par.Fin/facr)*facp;
		      double facr=(cosh(par.Fin*modJphi/(L+J.Jr)));
	double h = cJ*.5*(1+facr)*facp;
	double f = par.norm / pow_3(2*M_PI * par.J0) *
		   math::pow(1 + par.J0 / h,  par.slopeIn) *
		   math::pow(1 + h / par.J0, -par.slopeOut);
	if(par.Jcutoff>0)   // exponential cutoff at large J
		f *= exp(-math::pow(h / par.Jcutoff, par.cutoffStrength));
	if(par.Jcore>0) {   // central core of nearly-constant f(J) at small J
		if(h==0) return par.norm / pow_3(2*M_PI * par.J0);
		f *= math::pow(1 + par.Jcore/h * (par.Jcore/h - beta), -0.5*par.slopeIn);
	}
	if(par.rotFrac!=0)  // add the odd part
		f *= 1 + par.rotFrac * tanh(J.Jphi / par.Jphi0);
	return f;
}*/
}  // namespace df
