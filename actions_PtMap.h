#pragma once
#include "coord.h"
#include "math_base.h"
#include "math_core.h"
#include "math_spline.h"
#include "math_ode.h"
#include "potential_base.h"
#include "actions_base.h"
#include "actions_focal_distance_finder.h"

#define EXP __declspec(dllexport)

namespace actions2d{

class EXP PtClosedSyst: public math::IOdeSystem{
	private:
		math::CubicSpline u,pu,pv;
		double v_min;
		double L2, Lz;
		coord::UVSph CS;
	public:
		PtClosedSyst(std::vector<coord::PosVelCyl>& closedOrb,
			     const coord::UVSph& _CS, const double Ltot, const double _Lz)
				: CS(_CS), L2(Ltot*Ltot), Lz(_Lz) {
			std::vector<double> ui,vi,pui,pvi;
			for(size_t i=0;i<closedOrb.size();i++){
				coord::PosMomUVSph uvt(coord::toPosMom(closedOrb[i],CS));
				ui.push_back(uvt.u);
				vi.push_back(uvt.v);
				pui.push_back(uvt.pu);
				pvi.push_back(uvt.pv);
			}
			v_min=vi.back();
			u = math::CubicSpline(vi,ui);
			pu= math::CubicSpline(vi,pui);
			pv= math::CubicSpline(vi,pvi);
		}
		void eval(const double psi, const double x[], double dxdpsi[]) const{
			double cpsi=cos(psi), dvdpsi=-(.5*M_PI-v_min)*cpsi;
			double vc=.5*M_PI-(.5*M_PI-v_min)*sin(psi);
			double dudv,uc,puc,pvc;
			u.evalDeriv(vc,&uc,&dudv); pu.evalDeriv(vc,&puc); pv.evalDeriv(vc,&pvc);
			double pt=sqrt(L2-pow_2(Lz/sin(x[0]*x[1])));
			dxdpsi[0]=dvdpsi*puc*dudv/(x[1]*pt);
			dxdpsi[1]=dvdpsi*pvc/(x[0]*pt);
		}
		unsigned int size()  const {return 2;}
		double get_u(double v){
			double uc; u.evalDeriv(v,&uc);
			return uc;
		}
};
std::vector<coord::PosVelCyl> unpair(std::vector<std::pair<coord::PosVelCyl,double> >& traj2){
	std::vector<coord::PosVelCyl> traj;
	for(size_t i=0;i<traj2.size();i++)
		traj.push_back(traj2[i].first);
	return traj;
}
class EXP pointMapper{
	private:
		const potential::BasePotential& pot;
		double Rsh, Lz, Ltot, v_min;
		math::CubicSpline xs,ys,zs;
		coord::UVSph CS;
	public:
		pointMapper(const potential::BasePotential& _pot,const double _Rsh,const double _Lz)
				: pot(_pot), Rsh(_Rsh), Lz(_Lz){
			double Phi; pot.eval(coord::PosCyl(Rsh,0,0),&Phi);
			double Vesc=sqrt(-2*Phi);
			double timeCross;
			std::vector<std::pair<coord::PosVelCyl,double> > traj;
			actions::FindRzClosedOrbitV FcO(pot,Rsh,Lz,timeCross,traj);
			double Vsh = math::findRoot(FcO,0.01*Vesc,.95*Vesc,1e-3);
			Ltot = Lz + Rsh*Vsh;
			std::vector<coord::PosVelCyl> closedOrb=unpair(traj);
			CS.set(actions::fitFocalDistanceShellOrbit(closedOrb));
			PtClosedSyst yz_system(closedOrb,CS,Ltot,Lz);
			//The independent variable is psi st v=Pi/2-(Pi/2-v_min)*sin(psi)
			math::OdeSolverDOP853 yz_solver(yz_system);
			double ic[2]={1, .5*M_PI}; 
			yz_solver.init(ic);
			std::vector<double> u,v,x,y,z;
			//start with dummy point
			//r.push_back(Rsh/2); theta.push_back(-0.5);
			//xi.push_back(1); eta.push_back(ic[0]); zeta.push_back(ic[1]);
			//first real data point
			u.push_back(asinh(Rsh/CS.Delta)); v.push_back(.5*M_PI);
			x.push_back(1); y.push_back(ic[0]); z.push_back(ic[1]);
			double psi=0;
			do {
				yz_solver.doStep();
				psi=yz_solver.getTime();
				double vc=.5*M_PI-(.5*M_PI-v_min)*sin(psi);
				v.push_back(vc);
				double uc=yz_system.get_u(vc);
				u.push_back(uc);
				x.push_back(Rsh/(CS.Delta*sinh(uc)));
				y.push_back(yz_solver.getSol(psi,0));
				z.push_back(yz_solver.getSol(psi,1));
			} while (psi<.5*M_PI);
			//now add dunmy pointa at end
		/*	theta.push_back(t_max+.5);
			r.push_back(r[r.size()-1]);
			xi.push_back(xi[xi.size()-1]);
			eta.push_back(eta[eta.size()-1]);
			zeta.push_back(zeta[zeta.size()-1]);*/
			xs = math::CubicSpline(v,x);
			ys = math::CubicSpline(u,y);
			zs = math::CubicSpline(v,z);
		}
		coord::PosVelCyl FromToy(const coord::PosVelCyl Rzt){
			coord::PosMomUVSph uvt(coord::toPosMom(Rzt,CS));
			double y, dy; ys.evalDeriv(uvt.u,&y,&dy);
			double z, dz, colat=.5*M_PI-uvt.v;
			zs.evalDeriv(colat,&z,&dz);
			double x, dx; xs.evalDeriv(colat,&x,&dx);
			double uc=x*uvt.u, vc=.5*M_PI-y*z;
			double det=x*y*dzeta-dx*dy*uvt.u;
			double puc=(y*dz*uvt.pu-dy*z*uvt.u*uvt.pphi)/det;
			double pvp=-(dx*uvt.u*uvt.pu-x*uvt.u*uvt.pv)/det;
			coord::PosVelSph im(rp,thetap,uvt.phi,prp,ptp/rp,uvt.phi);
			return coord::PosVelCyl(coord::toPosVelCyl(im));
}
		coord::PosVelCyl ToToy(const coord::PosVelCyl& Rzt){
			coord::PosVelSph uvt(coord::toPosVelSph(Rzt));
			double eta, deta; etas.evalDeriv(uvt.r,&eta,&deta);
			double zeta, dzeta, colat=.5*M_PI-uvt.theta;
			zetas.evalDeriv(colat,&zeta,&dzeta);
			double xi, dxi; xis.evalDeriv(colat,&xi,&dxi);
			double rp=xi*uvt.r, thetap=.5*M_PI-eta*zeta;
			double det=xi*eta*dzeta-dxi*deta*uvt.r;
			double prp=(eta*dzeta*uvt.vr-deta*zeta*uvt.r*uvt.vphi)/det;
			double ptp=-(dxi*uvt.r*uvt.vr-xi*uvt.r*uvt.vtheta)/det;
			coord::PosVelSph im(rp,thetap,uvt.phi,prp,ptp/rp,uvt.phi);
			return coord::PosVelCyl(coord::toPosVelCyl(im));
		}
};
	
}//namespace