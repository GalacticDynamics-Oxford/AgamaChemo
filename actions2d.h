#pragma once
#include "coord.h"
#include "actions_focal_distance_finder.h"
#include "math_fit.h"
#include "potential_factory.h"
#include "actions_isochrone.h"

#define EXP __declspec(dllexport)

namespace actions2d{

struct Actions2d {
	double Jr;
	double Jphi;
	Actions2d() {};
	Actions2d(double _Jr,double _Jphi) : Jr(_Jr), Jphi(_Jphi) {};
};
struct Angles2d {
	double thetar;
	double thetaphi;
	Angles2d() {};
	Angles2d(double _thetar,double _thetaphi) : thetar(_thetar), thetaphi(_thetaphi) {};
};
struct ActionAngles2d : Actions2d, Angles2d {
	ActionAngles2d() {};
	ActionAngles2d(const Actions2d _J,const Angles2d _theta) : Actions2d(_J), Angles2d(_theta) {};
};
struct Frequencies2d {
	double Omegar;
	double Omegaphi;
	Frequencies2d() {};
	Frequencies2d(double _Omegar,double _Omegaphi) : Omegar(_Omegar), Omegaphi(_Omegaphi) {};
};
/* Class to find the variance of Hisochrone along a trajectory for
 * varying b with M(b) so the circular orbit of given action has a
 * given radius
*/
class EXP varFinder : public math::IFunction{
	private:
		std::vector<double> rs,KEs;
		double Jz, Rshell;
	public:
		varFinder(std::vector<std::pair<coord::PosVelCyl, double> >& traj,
			  double _Jz, double _Rshell);
		virtual void evalDeriv(const double var, double* value, double *derivs=NULL, double *der2=NULL) const;
		double sigH(const double var) const;
		virtual unsigned int numDerivs() const { return 0; }
		double Mass(const double b) const {
			double ksh=sqrt(1+pow_2(Rshell/b));
			return pow_2(Jz*b*(1+ksh))*b/pow_4(Rshell)*ksh;
		}
};
EXP Actions2d actionsIsochrone2d(const double M,const double b,const double q, const coord::PosVelCar& point);
EXP ActionAngles2d actionAnglesIsochrone2d(const double M,const double b,const double q, const coord::PosVelCar& point);
EXP Actions2d actionsHO2d(const double omegax,const double omegaz,const coord::PosVelCar& point);
EXP ActionAngles2d actionAnglesHO2d(const double omegax,const double omegaz,const coord::PosVelCar& point);
EXP actions2d::Actions2d getJs(const potential::BasePotential& pot,
			       std::vector<std::pair<coord::PosVelCyl, double> >& traj,
			       double& Mass, double& b_scale, double& Rsh, double& qR, double& qz, double& Jcl,
			       std::vector<double>& Sn, std::vector<std::pair<int,int> >& n);
EXP double checkH(potential::PtrPotential pot, double M, double b, const double q, int N, actions2d::Actions2d& J,
	      std::vector<double>& Sn, std::vector<std::pair<int,int> >& n);


}