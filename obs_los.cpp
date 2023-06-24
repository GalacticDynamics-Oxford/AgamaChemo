#pragma once
#include "obs.h"
#include "math_ODE.h"
//#include <cmath>


namespace obs{

EXP SunLos::SunLos(const PosSky _pos, const solarShifter _sun, dust::dustModel* _dm)//angles in radians// added for sKpc
: pos(_pos), cosl(cos(pos.l)), cosb(cos(pos.b)), sinl(sin(pos.l)), sinb(sin(pos.b)),
sun(_sun), BaseLos(_sun.fromKpc(), true, 0, _dm)
{
	if(_dm) tab_extinct();
}

EXP SunLos::SunLos(const coord::PosCar xyz, const solarShifter _sun, dust::dustModel* _dm)
: sun(_sun), BaseLos(_sun.fromKpc(), true, 0, _dm)// added for sKpc
{
	double xhel = xyz.x-sun.xyz().x, yhel=xyz.y-sun.xyz().y, zhel=xyz.z-sun.xyz().z;//position wrt Sun
	double s = sqrt(xhel*xhel+yhel*yhel+zhel*zhel);
	pos.b = asin(zhel/s); pos.l = atan2(yhel,xhel);
	cosl = cos(pos.l); sinl = sin(pos.l);
	cosb = cos(pos.b); sinb = sin(pos.b);
	if(_dm) tab_extinct();
}

EXP coord::PosCyl SunLos::deepest() const{//locate clostest distance r_p of los to GC
	double r_p = sun.xyz().x*cosb*cosl + sun.xyz().y*cosb*sinl
		     + sun.xyz().z*sinb;
	return Rzphi(r_p);
}

//helper class for tabulating extinctions: we integrate dN/ds = rho
//along LOS. Distances in Kpc
class A_syst: public math::IOdeSystem{
	private:
		const obs::BaseLos* Los;
		const dust::BaseDustModel* dm;
		const double start;
		const double from_Kpc;
	public:
		A_syst(const obs::BaseLos* _los, double _start, double _from_Kpc, dust::BaseDustModel* _dm)
				: Los(_los), start(_start), from_Kpc(_from_Kpc), dm(_dm){}
		virtual void eval(const double sKpc,const double A[], double dAds[]) const{
			dAds[0] = dm->dens(Los->Rzphi((start+sKpc)*from_Kpc));
		}
		virtual unsigned int size() const {return 1;}
};

EXP void BaseLos::tab_extinct(void){
	std::vector<double> dists, Ns;
	const int ns=3;
	A_syst syst(this, start(), from_Kpc, dm);
	math::OdeSolverDOP853 solver(syst);
	double sKpc=start(), rho=fmax(0,dm->dens(Rzphi(sKpc*from_Kpc))), N=0, rhoPeak=0;
	int iz=0;
	solver.init(&N);
	dists.push_back(sKpc); Ns.push_back(N);
	while(iz<200 && rho >= 0.01*rhoPeak){
		solver.doStep();
		double s_last=sKpc, ds=(start()+solver.getTime()-s_last)/(double)ns;
		for(int i=0; i<ns; i++){// store ns points over step
			sKpc+=ds;
			dists.push_back(sKpc); Ns.push_back(solver.getSol(sKpc-start(),0));
		}
		rho=dm->dens(Rzphi(sKpc*from_Kpc)); rhoPeak= rho>rhoPeak? rho : rhoPeak;
		iz++;
	}
	extinct = math::CubicSpline(dists,Ns);
}

EXP extLos::extLos(double _xs, double _ys, double _incl, double _s0,
		   double _fromKpc, dust::BaseDustModel* _dm)
: xs(_xs*_fromKpc), ys(_ys*_fromKpc), incl(_incl*M_PI/180.), 
cosi(cos(_incl)), sini(sin(_incl)), BaseLos(_fromKpc, false, _s0*_fromKpc, _dm)
{
	if(_dm) tab_extinct();
}

EXP coord::PosCyl extLos::deepest() const{
	return Rzphi(-ys*cosi/sini);
}

} //namespace obs