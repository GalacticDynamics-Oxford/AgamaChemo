#include "actions2d.h"

namespace actions2d {

EXP Actions2d actionsIsochrone2d(const double M,const double b, const double q, const coord::PosVelCar& point) {
	double R = sqrt(pow_2(point.x)+pow_2(point.z/q)), phi = atan2(point.z/q,point.x), z = 0;
	double vR = (point.x*point.vx + point.z*point.vz)/R;
	double vphi = (point.x*point.vz*q - point.z/q*point.vx)/R, vz = 0;
	coord::PosVelCyl pv(R,z,phi,vR,vz,vphi);
	actions::Actions J(actions::actionsIsochrone(M,b,pv));
	return Actions2d(J.Jr,J.Jphi);
}

EXP ActionAngles2d actionAnglesIsochrone2d(const double M,const double b, const double qR, const double qz, const coord::PosVelCar& point) {
	double R = sqrt(pow_2(point.x/qR)+pow_2(point.z/qz)), phi = atan2(point.z/qz,point.x/qR), z = 0;
	double vR = (point.x*point.vx + point.z*point.vz)/R;
	double vphi = (point.x/qR*point.vz*qz - point.z/qz*point.vx*qR)/R, vz = 0;
	coord::PosVelCyl pv(R,z,phi,vR,vz,vphi);
	actions::ActionAngles JA = actions::actionAnglesIsochrone(M,b,pv);
	Actions2d J(JA.Jr,JA.Jphi);
	Angles2d A(JA.thetar,JA.thetaphi);
	return ActionAngles2d(J,A);
}

EXP Actions2d actionsHO2d(const double omegax,const double omegaz,const coord::PosVelCar& point){
	return Actions2d(.5*(pow_2(point.vx)+pow_2(omegax*point.x))/omegax,
			 .5*(pow_2(point.vz)+pow_2(omegaz*point.z))/omegaz);
}

EXP ActionAngles2d actionAnglesHO2d(const double omegax,const double omegaz,const coord::PosVelCar& point){
	Actions2d J(.5*(pow_2(point.vx)+pow_2(omegax*point.x))/omegax,
		  .5*(pow_2(point.vz)+pow_2(omegaz*point.z))/omegaz);
	Angles2d A(atan2(omegax*point.x,point.vx),atan2(omegaz*point.z,point.vz));
	return ActionAngles2d(J,A);
}

/* Determine rms variation of H along torus defined by isochrone Phi
   * true actions J and the GF (Sn,n)
*/
EXP double checkH(potential::PtrPotential pot, double M, double b, const double q, int N, actions2d::Actions2d& J,
	      std::vector<double>& Sn, std::vector<std::pair<int,int> >& n){
	double* dH=new double[N*N];
	double varH=0, Hbar=0;
	int nlines=n.size();
	actions::Angles T(0,0,0);
	for(int i=0; i<N; i++){
		T.thetar=M_PI*((double)i+.5)/(double)N;
		for(int j=0; j<N; j++){
			T.thetaz=M_PI*((double)j+.5)/(double)N;
			actions::Actions Jtoy(J.Jr,J.Jphi,0);
			for(int k=0; k<nlines; k++){//compute isochrone actions
				double c = cos(n[k].first*T.thetar+n[k].second*T.thetaz);
				Jtoy.Jr+=2*Sn[k]*n[k].first * c;
				Jtoy.Jz+=2*Sn[k]*n[k].second* c;
			}
			Jtoy.Jr= Jtoy.Jr<0? 0:Jtoy.Jr;
			actions::ActionAngles aa(Jtoy,T);
			coord::PosVelCyl Rzphi(actions::mapIsochrone(M,b,aa));
			coord::PosVelCar xv(coord::toPosVelCar(Rzphi));
			xv.z*=q; xv.vz/=q;
			double H=potential::totalEnergy(*pot,xv);
			if(std::isnan(H)){
				//printf("ISNAN %g %g %g %g\n", T.thetar/M_PI,T.thetaz/M_PI,Jtoy.Jr,Jtoy.Jz);
				dH[i*N+j]=0;
			} else {
				varH+=pow_2(H); Hbar+=H;
				dH[i*N+j]=H;
			}
		}
	}
	varH /=N*N; Hbar/=N*N;
	double sigH = sqrt((varH-pow_2(Hbar)))/fabs(Hbar);
	printf("sigmaH/<H>: %g\n",sigH);
	for(int i=0;i<N*N;i++){
		float sz=40.0+fmin(.95,.5*fabs(dH[i]-Hbar));
		//if(i%N==0) printf("\n"); printf("%5.2f ",sz); 
		//pl.relocate(i/N,i%N); pl.point(sz);
	}
	return sigH;
}
double getJz(std::vector<std::pair<coord::PosVelCyl,double> >& traj){
	double Jz=0;
	for(size_t i=0;i<traj.size()-1;i++)
		Jz+=(traj[i].first.vR+traj[i+1].first.vR)*(traj[i+1].first.R-traj[i].first.R)
		    +(traj[i].first.vz+traj[i+1].first.vz)*(traj[i+1].first.z-traj[i].first.z);
	return Jz/M_PI;//we've only integrated 1/4 cycle but have twice that integral
}

class closedFitter: public math::IFunctionNdim{
	private:
		std::vector<std::pair<coord::PosVelCyl,double> >& closedOrb;
		double Jclosed, Rsh, q, tmax;
		size_t nt;
	public:
		closedFitter(std::vector<std::pair<coord::PosVelCyl,double> >& _closedOrb):
		    closedOrb(_closedOrb) {
			nt = closedOrb.size();
			tmax =closedOrb[nt-1].second;
			Jclosed = getJz(closedOrb);//Action of shell orbit
			Rsh=closedOrb[0].first.R;
			q = closedOrb[nt-1].first.z/Rsh;
			printf("Rsh, tmax, q, Jclosed: %f %f %f %f\n",Rsh,tmax,q,Jclosed);
		}
		virtual void eval(const double vars[], double values[]) const{
			double M=exp(vars[0]), b=exp(vars[1]), qR=exp(vars[2]);//vary logs to ensure >0
			double qz=q*qR;
		//	printf("qs: %f %f ",qR,qz);
			double var=0;
			for(size_t i=0;i<nt;i+=nt-1){
				double psi=.5*M_PI*closedOrb[i].second/tmax;
				double R=closedOrb[i].first.R, z=closedOrb[i].first.z;
				double vR=closedOrb[i].first.vR, vz=closedOrb[i].first.vz;
				//To avoid singularity on z axis, we
				//rotate so orbit is in xy plane of
				//spherical isochrone
				actions::Actions J(0,0,Jclosed); actions::Angles T(0,0,psi);
				actions::ActionAngles a(J,T);
				coord::PosVelCyl Rz(actions::mapIsochrone(M,b,a));
				coord::PosVelCar xv(coord::toPosVelCar(Rz));
		//		printf("psi: %f %f %f %f %f\n",psi,xv.x,xv.vx,xv.y,xv.vy);
				double add=pow_2(xv.x*qR-R)+pow_2(xv.y*qz-z)
					   +pow_2(xv.vx/qR-vR)+pow_2(xv.vy/qz-vz);
				var+=add;
				printf("%g ",add);

			}
			if(values) values[0]=var;
		}
		double give_q(){ return q; }
		double give_J(){ return Jclosed; }
		virtual unsigned int numVars() const { return 3; }
		virtual unsigned int numValues() const { return 1; }
};

varFinder::varFinder(std::vector<std::pair<coord::PosVelCyl, double> >& traj,
		     double _Jz, double _Rshell): Jz(_Jz), Rshell(_Rshell) {
	for(size_t k=0; k<traj.size(); k++){//save radii & KEs
		rs.push_back(sqrt(pow_2(traj[k].first.R)+pow_2(traj[k].first.z)));
		KEs.push_back(.5*(pow_2(traj[k].first.vR)+pow_2(traj[k].first.vz)+pow_2(traj[k].first.vphi)));
	}
}
//The value returned below is the derivative of varH wrt log b
void varFinder::evalDeriv(const double var, double* value, double* d1,double* d2) const{
	size_t nt=KEs.size();
	double b=exp(var), ksh=sqrt(1+pow_2(Rshell/b));//vary logb to ensure >0
	double M=pow_2(Jz*b*(1+ksh))*b/pow_4(Rshell)*ksh;
	double H2=0, dH2=0, Hbar=0, dHbar=0;
	for(size_t n=0;n<nt;n++){
		double k=sqrt(1+pow_2(rs[n]/b));
		double Bottom = b*(1 + k);
		double Phi=-M/Bottom;
		double H = KEs[n] + Phi;
		double dH=Phi*((1+2*ksh)/(ksh*ksh)-1+pow_2(rs[n]/b)/((1+k)*k));
		H2+=H*H; dH2+=2*H*dH; Hbar+=H; dHbar+=dH;
	}
	H2/=nt; dH2/=nt; Hbar/=nt; dHbar/=nt; 
	if(value) *value = dH2-2*Hbar*dHbar;
}
double varFinder::sigH(const double v) const{
	size_t nt=KEs.size();
	double b=exp(v);
	double ksh=sqrt(1+pow_2(Rshell/b));
	double M=pow_2(Jz*b*(1+ksh))*b/pow_4(Rshell)*ksh;
	double H2=0, Hbar=0;
	for(size_t n=0;n<nt;n++){
		double k=sqrt(1+pow_2(rs[n]/b));
		double Bottom = b*(1 + k);
		double Phi=-M/Bottom;
		double H = KEs[n] + Phi;
		H2+=H*H; Hbar+=H;
	}
	H2/=nt; Hbar/=nt;
	return (H2-Hbar*Hbar);
}
/* Helper class used to fit an Isochrone Phi to an orbits. We minimize
 * the variation of Hisochrone along the numerically integrated orbit
 */
class MLfunc: public math::IFunctionNdimDeriv{
	private:
		std::vector<double> rs,KEs;
	public:
		MLfunc(std::vector<std::pair<coord::PosVelCyl, double> >& traj){
			for(size_t k=0; k<traj.size(); k++){//save radii & KEs
				rs.push_back(sqrt(pow_2(traj[k].first.R)+pow_2(traj[k].first.z)));
				KEs.push_back(.5*(pow_2(traj[k].first.vR)+pow_2(traj[k].first.vz)+pow_2(traj[k].first.vphi)));
			}
		}
		virtual void evalDeriv(const double vars[], double values[], double *derivs) const{
			size_t nt=KEs.size();
			double M=exp(vars[0]), b=exp(vars[1]);//vary logM & logb to ensure >0
			double varH=0, Hbar=0, dHbar0=0, dHbar1=0;
			for(size_t k=0;k<nt;k++){
				double bottom = sqrt(pow_2(b)+pow_2(rs[k]));
				double Bottom = b + bottom;
				double H = KEs[k] - M/Bottom;
				if(values) values[k] = H; Hbar += H; varH+=pow_2(H);
				if(derivs){
					derivs[k*2] = -M/Bottom;
					derivs[k*2+1] = b*M/pow_2(Bottom)*(1+vars[1]/bottom);
					dHbar0 += derivs[k*2]; dHbar1 += derivs[k*2+1];
				}
			}
			Hbar /= nt; dHbar0 /= nt; dHbar1 /= nt; varH = varH/nt-pow_2(Hbar);
			for(size_t k=0;k<nt;k++){
				if(values) values[k] -= Hbar;
				if(derivs){
					derivs[2*k] -= dHbar0; derivs[2*k+1] -= dHbar1;
				}
			}
		}
		virtual unsigned int numVars() const { return 2; }
		virtual unsigned int numValues() const { return (unsigned int)KEs.size(); }

};

void getRminRmax(std::vector<std::pair<coord::PosVelCyl, double> >& traj,
		 double& Rmin, double& Rmax, double& Vmin, double& Vmax){
	Rmin=1e6; Rmax=-1e6; Vmin=1e6; Vmax=-1e6;
	for(size_t i=1; i<traj.size(); i++){
		if(traj[i].first.z * traj[i-1].first.z<0){//plane crossing
			double f = traj[i-1].first.z/(traj[i-1].first.z-traj[i].first.z);
			double Rc = fabs((1-f)*traj[i-1].first.R + f*traj[i].first.R);
			double Vc = fabs((1-f)*traj[i-1].first.vz + f*traj[i].first.vz);
			if(Rc>Rmax) Rmax=Rc; if(Rc<Rmin) Rmin=Rc;
			if(Vc>Vmax) Vmax=Vc; if(Vc<Vmin) Vmin=Vc;
		}
	}
}
	
void print_orb(std::vector<std::pair<coord::PosVelCyl,double> >& traj){
	FILE* ifile;
	if(fopen_s(&ifile,"orb.dat","w")) return;
	for(size_t i=0;i<traj.size();i++)
		fprintf(ifile,"%f %f %f %f\n",traj[i].first.R,traj[i].first.z,traj[i].first.vR,traj[i].first.vz);
	fclose(ifile);
}
	
			
/* Determine actions and GF cpts Sn from a numerically computed orbit
 * (traj) with Lz=0 by the Sanders&Binney algoritthm. On exit the Mass and b_scale relate to
 * the isochrone that fits traj best, Sn contains the GF cpts with the identification of the
 * lines given by (n.first,n.second)
 */
EXP actions2d::Actions2d getJs(const potential::BasePotential& pot,
			       std::vector<std::pair<coord::PosVelCyl, double> >& traj,
			       double& Mass, double& b_scale, double& Rsh, double& qR, double& qz, double& Jclosed,
			       std::vector<double>& Sn, std::vector<std::pair<int,int> >& n){
	size_t nt=traj.size();
	//Get R &Vz of plane crossings
	double Rmin,Rmax,Vmin,Vmax;
	getRminRmax(traj,Rmin,Rmax,Vmin,Vmax);
	//Get shell through Rsh
	Rsh=.5*(Rmin+Rmax);
	double Phi; pot.eval(coord::PosCyl(Rsh,0,0),&Phi);
	double Vesc=sqrt(-2*Phi);
	printf("Rmin/max: %f %f Vesc %f Vmin/max: %f %f\n",Rmin,Rmax,Vesc,Vmin,Vmax);
	double timeCross;
	std::vector<std::pair<coord::PosVelCyl,double> > closedOrb;
	actions::FindRzClosedOrbitV FcO(pot,Rsh,0,timeCross,closedOrb);
	double Vsh = math::findRoot(FcO,0.01*Vesc,.95*Vesc,1e-3);
	closedFitter cF(closedOrb);
	double xinit[3]={0,0,0}, xstep[3]={.2,.2,.2}, results[3];
	int nops=math::findMinNdim(cF,xinit,xstep,1e-3,100,results);
	Mass=exp(results[0]); b_scale=exp(results[1]);
	qR=exp(results[2]);
	qz=qR*cF.give_q(); Jclosed=cF.give_J();
	printf("\nnops etc: %d %f %f %f %f\n",nops,Mass,b_scale,qR,qz);
//	Jclosed = getJz(closedOrb);//Action of shell orbit
//	q = closedOrb[closedOrb.size()-1].z/Rsh;
//	printf("Rsh, q, Jclosed: %f %f %f\n",Rsh,q,Jclosed);
	print_orb(closedOrb);

/*
	//get b of isochrone that minimises varH(traj) given that circ
	//orbit radius Rsh has action Jclosed
	varFinder vF(traj,Jclosed,Rsh);
	double dera; vF.evalDeriv(0,&dera);
	double eps=1e-3, val=vF.sigH(eps);
	val-=vF.sigH(-eps);
	printf("vF %g %g\n",.5*val/eps,dera);
	double x;
	for(x=-3; x<5; x+=.2){
		vF.evalDeriv(x,&dera); printf("(%f %f) ",exp(x),dera);
		if(dera>0) break;
	}
	b_scale = exp(math::findRoot(vF,x-0.2,x,1e-4));
	Mass = vF.Mass(b_scale);
	printf("best M,b = %f %f\n",Mass,b_scale);
	if(std::isnan(b_scale)) exit(0);  */
/*	
	MLfunc dH(traj);
	double yinit[2]={0,0}, best[2];
	int niter = nonlinearMultiFit(dH, yinit, 1e-3, 100, best);
	Mass = exp(best[0]), b_scale = exp(best[1]);
	printf("After %d iterations, best M,b = %f %f\n",niter,Mass,b_scale);
	*/
	std::vector<actions2d::ActionAngles2d> AAs(nt);
	double addr=0,addphi=0;
	for(size_t k=0; k<nt; k++){//compute toy AAs
		AAs[k]=actionAnglesIsochrone2d(Mass,b_scale,qR,qz,coord::toPosVelCar(traj[k].first));
		AAs[k].thetar+=addr; AAs[k].thetaphi+=addphi;
		if(k>0){//eliminate discontinuities in angles
			if(fabs(AAs[k].thetar-AAs[k-1].thetar)>M_PI){
				if(AAs[k].thetar>AAs[k-1].thetar){
					addr-=2*M_PI; AAs[k].thetar-=2*M_PI;
				}else{
					addr+=2*M_PI; AAs[k].thetar+=2*M_PI;
				}
			}
			if(fabs(AAs[k].thetaphi-AAs[k-1].thetaphi)>M_PI){
				if(AAs[k].thetaphi>AAs[k-1].thetaphi){
					addphi-=2*M_PI; AAs[k].thetaphi-=2*M_PI;
				}else{
					addphi+=2*M_PI; AAs[k].thetaphi+=2*M_PI;
				}
			}
			//if(k<200)printf("%5.1f %5.1f\n",AAs[k].thetar,AAs[k].thetaphi);
		}
	}
	double Dthetar=AAs[nt-1].thetar-AAs[0].thetar;
	double Dthetaphi=AAs[nt-1].thetaphi-AAs[0].thetaphi;
	int nmax=8;// 8 implies going to nr=7, nphi=6
	for(int i=0; i<=nmax; i++){//pick well samples lines
		if(i>0){
			if(fabs(i*Dthetar)<(nt-1)*M_PI && fabs(i*Dthetar)>2*M_PI)
				n.push_back(std::make_pair(i,0));
		}
		for(int j=1;j<=nmax/2;j++){
			double ndtheta=fabs(i*Dthetar+2*j*Dthetaphi);
			if(i>0 && ndtheta<(nt-1)*M_PI && ndtheta>2*M_PI)
				n.push_back(std::make_pair(i, 2*j));
			ndtheta=fabs(i*Dthetar-2*j*Dthetaphi);
			if(ndtheta<(nt-1)*M_PI && ndtheta>2*M_PI)
				n.push_back(std::make_pair(i,-2*j));
		}
	}
	int nlines=n.size();
	printf("Using %d lines ",nlines);
	math::Matrix<double> A(2+nlines,2+nlines);
	std::vector<double> b(2+nlines);
	for(int i=0; i<2+nlines; i++) b[i]=0;
	for(size_t k=0; k<nt; k++){
		b[0]+=AAs[k].Jr; b[1]+=AAs[k].Jphi;
	}
	A(0,0)=nt; A(0,1)=0; A(1,0)=0; A(1,1)=nt;
	for(int i=0; i<nlines; i++){
		double sum1=0;
		for(size_t k=0; k<nt; k++){
			double ci=2*cos(n[i].first*AAs[k].thetar+n[i].second*AAs[k].thetaphi);
			sum1+=ci;
			b[2+i]+=ci*(n[i].first*AAs[k].Jr+n[i].second*AAs[k].Jphi);
		}
		A(0,2+i)=n[i].first*sum1;
		A(1,2+i)=n[i].second*sum1;
		A(2+i,0)=A(0,2+i);//A is symmetric
		A(2+i,1)=A(1,2+i);
		for(int j=i; j<nlines; j++){
			double sum2=0;
			for(size_t k=0; k<nt; k++){
				double ci=2*cos(n[i].first*AAs[k].thetar+n[i].second*AAs[k].thetaphi);
				double cj=2*cos(n[j].first*AAs[k].thetar+n[j].second*AAs[k].thetaphi);
				sum2+=ci*cj;
			}
			A(2+i,2+j)=sum2*(n[i].first*n[j].first + n[i].second*n[j].second);
			if(j>i) A(2+j,2+i)=A(2+i,2+j);
		}
	}
	math::LUDecomp K(A);
	Sn = K.solve(b);
	actions2d::Actions2d Js(Sn[0],Sn[1]);
	for(size_t i=0;i<nlines;i++) Sn[i]=Sn[i+2];//move GF cpts down
	return Js;
}

}//namespace