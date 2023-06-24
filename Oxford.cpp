#include "df_halo.h"
#include "math_core.h"
#include "math_ode.h"
#include "Oxford.h"
#include <cmath>
#include <stdexcept>

#define DE

namespace df{

/// A class and a function required by Oxford & Bologna
void Oxford_z_syst::eval(const double t,const double s[],double dsdt[]) const{
	//returns Om_z/Om_phi varies with inclination
	double i;//inclination
	if(s[1]>0) i=s[1]/(s[0]+s[1]);	else i=0;
	dsdt[0]=i+(1-i)*z_rat;//z_rat=Om_z/Om_phi planar
	dsdt[1]=-1;//t=const-Jz
}

void Oxford_r_syst::eval(const double t,const double s[],double dsdt[]) const{
	//returns Om_r/Om_phi as unction of circularity
	double c;//circularity
	if(s[1]>0) c=s[0]/(s[0]+s[1]);	else c=1;
	dsdt[0]=1/(.51+c*(r_rat-.51));//r_rat=Om/kappa
	dsdt[1]=-1;//t=const-Jr
}
/*double Oxford::z_rat(const actions::Actions &J) const{
	const double sixth=1./6., L=J.Jz+fabs(J.Jphi);
	const double jt=(1.5*J.Jr+L)/par.L0;
	return par.kIn+jt/(sixth+jt)*(par.kOut-par.kIn);// Omz/Omphi at this E
}*/
EXP double Oxford::h(const actions::Actions &J) const{
	//Move over E surface to planar orbit then to circular orbit &
	//return Lc
/*	const double sixth=1./6., L=J.Jz+fabs(J.Jphi), c0=L/(L+J.Jr);
	const double jt=(1.5*J.Jr+L)/par.L0,powj=pow(jt,par.gamma1);
	const double _z_rat
	const double _z_rat=par.kIn+jt/(sixth+jt)*(par.kOut-par.kIn);// Omz/Omphi at this E
	const double _r_rat=.5*(1+powj/(1+powj)+par.gamma0);*/
	const double L=J.Jz+fabs(J.Jphi), c0=L/(L+J.Jr);
	const double Lc=(1.4*J.Jr+L);//estimated Lc(E)
	double rats[2];	freq.epicycle_ratios(Lc,rats);//rats{Om/kappa,nu/Om} at this E
	Oxford_z_syst z_system(rats[1]);
	math::OdeSolverDOP853 z_solver(z_system);
	double ic[2]={fabs(J.Jphi)/par.L0,J.Jz/par.L0};
	z_solver.init(ic);
	double t=0, jz=ic[1]; int iz=0;
	while(iz<20 && jz>0){//until Jz<=0
		z_solver.doStep(); t=z_solver.getTime(); jz=z_solver.getSol(t,1);
		iz++;
	}
	double jp=ic[0];
	if(iz) jp=z_solver.getSol(ic[1],0);/*Jphi/L0 when Jz=0*/
	double icr[2]={jp, J.Jr/par.L0};
	if(iz==20 || icr[0]<0){
		printf("z_solver error: %g %g\n(%g %g)\n",
		       t,icr[0],ic[0],ic[1]); exit(0);
	}
//	if(jp>=0) return L0*jp;
	Oxford_r_syst r_system(rats[0]);
	math::OdeSolverDOP853 r_solver(r_system);
	r_solver.init(icr);
	t=0; int j=0; double jr=icr[1];
	while(j<20 && jr>0){//until Jr=0
		r_solver.doStep(0); t=r_solver.getTime(); jr=r_solver.getSol(t,1);
		j++;
	}
	jp=icr[0];
	if(j) jp=r_solver.getSol(icr[1],0);/*Jphi/L0 when Jr=0*/
	if(j==20){
		double s[2]={r_solver.getSol(t,0),r_solver.getSol(t,1)};
		printf("r_solver error: %g (%g %g) (%g %g)\n",t,icr[0],icr[1],s[0],s[1]);
		exit(0);
	}
	return par.L0*jp;//return Jphi when Jr=0
}

EXP double Oxford::value(const actions::Actions &J) const
{
#ifdef DE
	double hJ=h(J);
	double gJ=hJ;
//	if(par.coefLout==par.coefLin) gJ=hJ;
//	else gJ=h(J,par.coefLout);
#else
	double fJphi=fabs(J.Jphi), L=J.Jz+fJphi;
	double Jt=6*(1.5*J.Jr+L)/par.L0, zeta=Jt/(1+Jt);
	double k=par.kIn + zeta*(par.kOut-par.kIn);
	double a=.5*(k+1), b=.5*(k-1);
	double cL=(a*J.Jz+b*J.Jz*fJphi/L+fJphi);
	//double cL=Oxford_h(J,par.coefLout);
	double powJt=pow((1.5*J.Jr+L)/par.L0,par.gamma1), xi=powJt/(1+powJt);
	double c=L/(L+J.Jr), stheta=sin(.5*M_PI*c), fac=exp(-.5*par.coefLin*stheta);
	double hJ=J.Jr * fac + .5*(1+xi*c)/fac * cL;
	fac=exp(-.5*par.coefLout*stheta);
	double gJ=J.Jr * fac + .5*(1+xi*c)/fac * cL;
#endif
	double val =  norm *
		      math::pow(1 + math::pow(par.J0 / hJ, par.steepness),  par.slopeIn  / par.steepness) *
		      math::pow(1 + math::pow(gJ / par.J0, par.steepness), -par.slopeOut / par.steepness);
	if(par.rotFrac!=0)  // add the odd part
		val *= 1 + par.rotFrac * tanh(J.Jphi / par.Jphi0);
	if(par.Jcutoff>0)   // exponential cutoff at large J
		val *= exp(-math::pow((J.Jr+J.Jz+fabs(J.Jphi)) / par.Jcutoff, par.cutoffStrength));
	return g(hJ) * val;
}

EXP Oxford::Oxford(const OxfordParam &inparams, const potential::Interpolator& infreq) :  par(inparams), freq(infreq)
{
	norm=1;
    // sanity checks on parameters
	if(par.mass == 0)
		throw std::invalid_argument("Oxford: mass must be non-zero");
	if(!(par.J0>0))
		throw std::invalid_argument("Oxford: break action J0 must be positive");
	if(par.Jcutoff<0)
		throw std::invalid_argument("Oxford: truncation action Jcutoff must be non-negative");
	if(!(par.slopeOut>3) && par.Jcutoff==0)
		throw std::invalid_argument(
					    "Oxford: mass diverges at large J (outer slope must be > 3)");
	if(!(par.slopeIn<3))
		throw std::invalid_argument(
					    "Oxford: mass diverges at J->0 (inner slope must be < 3)");
	if(par.steepness<=0)
		throw std::invalid_argument(
					    "Oxford: transition steepness parameter must be positive");
	if(par.cutoffStrength<=0)
		throw std::invalid_argument("Oxford: cutoff strength parameter must be positive");
	if( par.L0 <=0 )
		throw std::invalid_argument("Oxford: L0 must be >0");
	if( par.kIn < 1 || par.kOut < 1  )
		throw std::invalid_argument("Oxford: Must have k >= 1");
	if( par.coefLin > 1.5 || par.coefLout > 1.5  )
		throw std::invalid_argument("Oxford: Must have coefL < 1.5");
	if(fabs(par.rotFrac)>1)
		throw std::invalid_argument(
					    "Oxford: amplitude of odd-Jphi component must be between -1 and 1");
	set_beta();
	printf("beta: %f\n",beta);
	/*printf("beta,mass: %f %f %g\n",beta,norm,totalMass());*/
	if(par.mass>0)
		norm*=par.mass/totalMass();
	else
		norm = -par.mass/pow_3(par.J0);

}
EXP void Oxford::set_norm(double fac){
	norm *= fac;
}
EXP double Oxford::g(const double hJ) const{
	double rat=par.Jcore/hJ;
	return pow(rat*rat-beta*rat+1,-.5*par.slopeIn);
}
EXP void Oxford::write_params(std::ofstream &strm,const units::InternalUnits &intUnits) const{
	strm << "mass\t\t" << par.mass*intUnits.to_Msun << '\n';
	strm << "J0\t\t" << par.J0*intUnits.to_Kpc_kms << '\n';
	strm << "Jcutoff\t\t" << par.Jcutoff*intUnits.to_Kpc_kms << '\n';
	strm << "Jphi0\t\t" << par.Jphi0*intUnits.to_Kpc_kms << '\n';
	strm << "Jcore\t\t" << par.Jcore*intUnits.to_Kpc_kms << '\n';
	strm << "slopeIn\t\t" << par.slopeIn << '\n';
	strm << "slopeOut\t" << par.slopeOut << '\n';
	strm << "steepness\t" << par.steepness << '\n';
	strm << "coefLin\t" << par.coefLin << '\n';
	strm << "coefLout\t" << par.coefLout << '\n';
	strm << "kIn\t" << par.kIn << '\n';
	strm << "kOut\t" << par.kOut << '\n';
	strm << "L0\t" << par.L0 << '\n';
	strm << "gamma0\t" << par.gamma0 << '\n';
	strm << "rotFrac\t\t" << par.rotFrac << '\n';
	strm << "cutoffStrength\t" << par.cutoffStrength << '\n';
}
EXP void Oxford::tab_params(std::ofstream &strm,const units::InternalUnits &intUnits) const{
	strm << "mass\t J0\t Jcutoff\t Jphi0\t Jcore\t L0\t slopeIn\t slopeOut\t steepness\t coefLin\n";
	strm << "kIn\t coefLout\t kOut\t rotFrac\t cutoffStrength\t gamma\n";
	strm << par.mass*intUnits.to_Msun << " & "
			<< par.J0*intUnits.to_Kpc_kms << " & "
			<< par.Jcutoff*intUnits.to_Kpc_kms << " & "
			<< par.Jphi0*intUnits.to_Kpc_kms << " & "
			<< par.Jcore*intUnits.to_Kpc_kms << " & "
			<< par.L0*intUnits.to_Kpc_kms << " & "
			<< par.slopeIn << " & "
			<< par.slopeOut << " & "
			<< par.coefLin << " & "
			<< par.kIn << " & "
			<< par.coefLout << " & "
			<< par.kOut << " & "
			<< par.rotFrac << " & "
			<< par.cutoffStrength << " & "
			<< par.gamma0 << " \\cr \n";
}
double Oxford::diff(const double hJ) const{
	double val = 
		    math::pow(1 + math::pow(par.J0 / hJ, par.steepness),  par.slopeIn  / par.steepness) *
		    math::pow(1 + math::pow(hJ / par.J0, par.steepness), -par.slopeOut / par.steepness);
	return hJ*hJ*val*(1-g(hJ));
}
double Oxford::intDiff(const double beta_) const{
	Oxford CB(par,freq,beta_);
	OxfordDiff Diff(CB);
	return math::integrateGL(Diff,1e-5,10*par.Jcore,10);
}
EXP void Oxford::set_beta(void){
	beta=0;
	OxfordInt Int(*this);
	beta=math::findRoot(Int,0,1,1e-5);
}

double Bologna_rat(const double s[],double c,double gamma){
	return .5*(1+c*tanh(s[0]+1.5*s[1]+gamma));
}

EXP Bologna::Bologna(const BolognaParam &inparams) :
    par(inparams)
{
	norm=1;
    // sanity checks on parameters
	if(par.mass==0)
		throw std::invalid_argument("Bologna: normalization must be non-zero");
	if(!(par.J0>0))
		throw std::invalid_argument("Bologna: break action J0 must be positive");
	if(par.alpha < 0 || par.kIn < 0 || par.kOut < 0) 
		throw std::invalid_argument("Bologna: steepness and k parameters must be positive");
	if(par.mass>0)
		norm *= par.mass/totalMass();
	else
		norm = -par.mass/pow_3(par.J0);
}

EXP double Bologna::value(const actions::Actions &J) const {
	double fJphi=fabs(J.Jphi), L=J.Jz+fJphi;
	double Jt=2*(1.5*J.Jr+L)/par.L0, zeta=Jt/(1+Jt);
	double k=par.kIn + zeta*(par.kOut-par.kIn);
	//Now that we have good estimate of k compute cL and xi
	double a=.5*(k+1), b=.5*(k-1);
	double cL=a*J.Jz+b*J.Jz*fJphi/(L+par.L1)+fJphi;
	double ic[2]={cL/par.L0,J.Jr/par.L0};
	double gJ=0;  //par.L0*get_h(1,ic,par.coefL,k,par.gamma,NAN);
/*	double powJt=pow((1.5*J.Jr+L)/par.L0,par.delta), xi=tanh(powJt);
	double c=L/(L+J.Jr);
	if(par.coefL!=0){
		double fac=exp(-par.coefL*sin(.5*M_PI*c));
		gJ=J.Jr * fac + .5*(1+xi*c)/fac * cL;
	} else {
		gJ=J.Jr + .5*(1+xi*c) * cL;
	}*/
	double val = norm * exp(-pow( gJ / par.J0, par.alpha));
	if(par.rotFrac!=0)  // add the odd part
		val *= 1 + par.rotFrac * tanh(J.Jphi / par.Jphi);
	return val; 
}
}