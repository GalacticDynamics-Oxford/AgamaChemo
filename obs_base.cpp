#pragma once
#include "obs.h"
#include <cmath>


namespace obs{


//Three constants defining Galactic coordinates
//#define	ap (192.85948*torad)
//#define dp (27.12825*torad)
//#define l0 (122.931918*torad)

#define	RA_GP (192.85948*torad)
#define decGP (27.12825*torad)
#define lCP (122.931918*torad)

const double torad=M_PI/180., sndecGP=sin(decGP), csdecGP=cos(decGP);

EXP PosSky from_RAdec(double RA,double dec){//in and out in degrees
	const double alpha=RA*torad, delta=dec*torad;
	double b=asin(sndecGP*sin(delta)+csdecGP*cos(delta)*cos(alpha-RA_GP));
	double l=lCP-atan2(cos(delta)*sin(alpha-RA_GP),
			  csdecGP*sin(delta)-sndecGP*cos(delta)*cos(alpha-RA_GP));
	PosSky p(l/torad, b/torad, false);
	return p;
}
EXP PosSky from_RAdec(PosSky p){
	assert(p.is_ra);
	return from_RAdec(p.l,p.b);
}
EXP PosSky to_RAdec(double l,double b){//in and out in degrees
	l*=torad; b*=torad;
	double delta=asin(sndecGP*sin(b)+csdecGP*cos(b)*cos(lCP-l));
	double alpha=RA_GP + atan2(cos(b)*sin(lCP-l),
			      csdecGP*sin(b)-sndecGP*cos(b)*cos(lCP-l));
	PosSky p(alpha/torad, delta/torad, true);
	return p;
}
EXP PosSky to_RAdec(PosSky p){
	assert(!p.is_ra);
	return to_RAdec(p.l,p.b);
}
//muRA =RA dot cos(delta) & similarly for mul
EXP PosVelSky from_muRAdec(double RA,double dec,double muRA,double mudec){
	PosSky poslb(from_RAdec(RA,dec));
	RA*=torad; dec*=torad;
	double sb,cb,l=poslb.l*torad,sd,cd;
	math::sincos(poslb.b*torad,sb,cb);
	math::sincos(dec,sd,cd);
	double mub=((sndecGP*cd-csdecGP*sd*cos(RA-RA_GP))*mudec
		    -csdecGP*sin(RA-RA_GP)*muRA)/cb;
	double mul;
	if(fabs(cos(lCP-l))>fabs(sin(lCP-l)))
		mul=(sd*sin(RA-RA_GP)*mudec-cos(RA-RA_GP)*muRA
		     -sb*sin(lCP-l)*mub)/cos(lCP-l);
	else
		mul=((csdecGP*cd+sndecGP*sd*cos(RA-RA_GP))*mudec
		     +sndecGP*sin(RA-RA_GP)*muRA+sb*cos(lCP-l)*mub)/sin(lCP-l);
	VelSky pm(mul,mub,false);
	PosVelSky p(poslb,pm);
	return p;
}
EXP PosVelSky from_muRAdec(PosSky radec,VelSky pmradec){
	assert(radec.is_ra && pmradec.is_ra);
	return from_muRAdec(radec.l,radec.b,pmradec.mul,pmradec.mub);
}
EXP PosVelSky from_muRAdec(PosVelSky p){//p in ra dec of course
	assert(p.is_ra);
	return from_muRAdec(p.pos,p.pm);
}
EXP PosVelSky to_muRAdec(double l,double b,double mul,double mub){
	PosSky posra(to_RAdec(l,b));
	l*=torad; b*=torad;
	double sb, cb, RA=posra.l*torad, dec=posra.b*torad, sd, cd;
	math::sincos(b,sb,cb);
	math::sincos(dec,sd,cd);
	double A=(sndecGP*cd-csdecGP*sd*cos(RA-RA_GP))/cb;
	double B=csdecGP*sin(RA-RA_GP)/cb;
	double C,D,E,clC=cos(lCP-l),slC=sin(lCP-l);
	if(fabs(clC)>fabs(slC)){
		C=sd*sin(RA-RA_GP)/clC;
		D=cos(RA-RA_GP)/clC;
		E=sb*sin(lCP-l)/clC;
	} else {
		C=(csdecGP*cd+sndecGP*sd*cos(RA-RA_GP))/slC;
		D=-sndecGP*sin(RA-RA_GP)/slC;
		E=-sb*cos(lCP-l)/slC;
	}
	double det=A*D-C*B;
	double mudec=((D-E*B)*mub-B*mul)/det;
	double mura =((C-A*E)*mub-A*mul)/det;
	VelSky pm(mura,mudec,true);
	PosVelSky p(posra,pm);
	return p;
}
EXP PosVelSky to_muRAdec(PosSky lb,VelSky pm){
	return to_muRAdec(lb.l,lb.b,pm.mul,pm.mub);
}
EXP PosVelSky to_muRAdec(PosVelSky pv){
	return to_muRAdec(pv.pos,pv.pm);
}

EXP solarShifter::solarShifter(const units::InternalUnits& intUnits,
			       coord::PosVelCar* _Sun):
    from_Kpc(intUnits.from_Kpc), from_kms(intUnits.from_kms), from_mas_per_yr(intUnits.from_mas_per_yr)
{
	torad = M_PI/180.;
	if(_Sun) Sun=*_Sun;
	//Gravity collab AA 657 L12 (2022) Read Brunthaler ApJ 892 39
	//(2020) Schoenrich MN 427 274 (2012)
	else Sun=coord::PosVelCar(-8.27 * from_Kpc,0,0.025 * from_Kpc,
				  14 * from_kms,251.3 * from_kms,7 * from_kms);
}
EXP coord::PosCar solarShifter::toCar(const PosSky pos, double sKpc) const{
	assert(!pos.is_ra);
	double cb,sb,cl,sl;
	math::sincos(torad*pos.b,sb,cb);
	math::sincos(torad*pos.l,sl,cl);
	double s = sKpc*from_Kpc;
	return coord::PosCar(Sun.x+s*cb*cl, Sun.y+s*cb*sl, Sun.z+s*sb);
}	
EXP coord::PosVelCar solarShifter::toCar(const PosSky pos, double sKpc,
					 const VelSky pm, double Vlos_kms) const{
	assert(!pos.is_ra && !pm.is_ra);
	double cb,sb,cl,sl;
	math::sincos(torad*pos.b,sb,cb);
	math::sincos(torad*pos.l,sl,cl);
	double s = sKpc*from_Kpc, Vlos = Vlos_kms*from_kms;
	double X = Sun.x+s*cb*cl;
	double Y = Sun.y+s*cb*sl;
	double Z = Sun.z+s*sb;
	double U = Sun.vx+Vlos*cb*cl-s*(sb*cl*pm.mub+sl*pm.mul)*from_mas_per_yr;
	double V = Sun.vy+Vlos*cb*sl-s*(sb*sl*pm.mub-cl*pm.mul)*from_mas_per_yr;
	double W = Sun.vz+Vlos*sb+s*cb*pm.mub*from_mas_per_yr;
	coord::PosVelCar O(X,Y,Z,U,V,W);//xv in internal units
	return O;
}
EXP coord::PosVelCyl solarShifter::toCyl(const PosSky pos, const double sKpc, const VelSky pm, double Vlos_kms) const{
	return toPosVelCyl(toCar(pos,sKpc,pm,Vlos_kms));
}
EXP PosSky solarShifter::toSky(const coord::PosCar p, double& sKpc) const{
	double xhel = p.x-Sun.x, yhel=p.y-Sun.y, zhel=p.z-Sun.z;//position wrt Sun
	double s = sqrt(xhel*xhel+yhel*yhel+zhel*zhel);
	sKpc = s/from_Kpc; //return in kpc
	double b=asin(zhel/s), l=atan2(yhel,xhel);
	return PosSky(l/torad, b/torad, false);
}
EXP PosVelSky solarShifter::toSky(const coord::PosVelCar pv, double& sKpc, double& Vlos_kms) const{
	double xhel = pv.x-Sun.x, yhel=pv.y-Sun.y, zhel=pv.z-Sun.z;//position wrt Sun
	double s = sqrt(xhel*xhel+yhel*yhel+zhel*zhel);
	sKpc = s/from_Kpc; //return in kpc
	double Vhelx = pv.vx-Sun.vx, Vhely = pv.vy-Sun.vy, Vhelz = pv.vz-Sun.vz;//V wrt Sun
	double Vlos = (xhel*Vhelx + yhel*Vhely + zhel*Vhelz)/s;
	Vhelx -= Vlos*xhel/s; Vhely -= Vlos*yhel/s; Vhelz -= Vlos*zhel/s; //sky velocity
	Vlos_kms = Vlos/from_kms; //return Vlos in kms
	double b=asin(zhel/s), l=atan2(yhel,xhel);
	double sb,cb,sl,cl;
	math::sincos(b,sb,cb); math::sincos(l,sl,cl);
	double mub = Vhelz/(s*cb);
	double mul;
	if(fabs(sl)>fabs(cl))
		mul=-(Vhelx+Vhelz*sb*cl/cb)/(s*sl);
	else
		mul=(Vhely+Vhelz*sb*sl/cb)/(s*cl);
	PosVelSky p(l/torad, b/torad, mul/from_mas_per_yr, mub/from_mas_per_yr, false);
	return p;
}

EXP VelSky solarShifter::toPM(const coord::PosVelCar pv, double& Vlos_kms) const{
	double sKpc;
	PosVelSky p = toSky(pv, sKpc, Vlos_kms);
	return p.pm;
}


} // namespace obs