/** \file    obs_base.h 
    \brief   Sky coordinates & functions relating to observations
    \author  Fred Thompson & James Binney
    \date    2021

This module provides code to facilitate interaction with observations.
   Polar sky coords (right acsension, declination) or Galactic (longitude, latitude)
   are defined alongside associated proper motions (with mul = dot\ell\cos(b), etc).
   A solarhifter moves data between Galactocentric and sky coordinates. Lines of
   sight are defined with functions to compute Galctocentric coords to/from distance down los.
   Extinction down the los is available if a dustModel is specfied.

*/
#pragma once

namespace obs {

/* Structs for sky coordinates and proper motions. Units are degrees
 * and mas/yr rather than internal units
 */
struct PosSky{
	double l, b;
	bool is_ra;
	PosSky(){};
	PosSky(double _l, double _b, bool _is_ra = false) : l(_l), b(_b), is_ra(_is_ra) {};
};

struct VelSky{
	double mul, mub;
	bool is_ra;
	VelSky(double _mul, double _mub, bool _is_ra = false) : mul(_mul), mub(_mub), is_ra(_is_ra){}
};

struct PosVelSky{
	bool is_ra;
	PosSky pos;
	VelSky pm;
	PosVelSky(PosSky _pos, VelSky _pm) : pos(_pos), pm(_pm){}
	PosVelSky(double l, double b, double mul, double mub, bool is_ra = false) : pos(l,b,is_ra), pm(mul,mub,is_ra){}
};

/* Transformations between equatorial & Galactic coords (units degreea)
 */
EXP PosSky from_RAdec(double ra,double dec);
EXP PosSky from_RAdec(PosSky p);
EXP PosSky to_RAdec(double l,double b);
EXP PosSky to_RAdec(PosSky p);
//PMs are mura=ra dot cos(delta) and mul=l dot cos(b)
EXP PosVelSky from_muRAdec(double ra,double dec,double mura,double mudec);
EXP PosVelSky from_muRAdec(PosSky radec,VelSky pmradec);
EXP PosVelSky from_muRAdec(PosVelSky p);//p in ra dec of course
EXP PosVelSky to_muRAdec(double l,double b,double mul,double mub);
EXP PosVelSky to_muRAdec(PosSky lb,VelSky pm);
EXP PosVelSky to_muRAdec(PosVelSky pv);

/* Class to convert obs coords (l,b in deg, mu in mas/yr, s in kpc,
 * Vlos in km/s) into 6d Galactocentric phase space coords in int units
 * and back. Sun's Galactocentric phase-space coords (if specified) should be in Kpc &
 * kms. Similarly s in kpc in and out; Vlos in kms in and out
 */
class EXP solarShifter{
	private:
		coord::PosVelCar Sun;
	public:
		double from_mas_per_yr, torad, from_Kpc, from_kms;
		solarShifter(const units::InternalUnits& intUnits, coord::PosVelCar* _Sun=NULL);
		// l, b in degrees, s in kpc, VelSky in mas/yr, Vlos in kms
		coord::PosCar toCar(PosSky pos, double sKpc) const;
		coord::PosVelCar toCar(const PosSky p,double sKpc,const VelSky pm,double Vlos_kms) const;
		coord::PosVelCyl toCyl(const PosSky p,double sKpc,const VelSky pm,double Vlos_kms) const;
		PosSky toSky(const coord::PosCar p, double& sKpc) const;
		PosSky toSky(const coord::PosCyl p, double& sKpc) const{
			return toSky(coord::toPosCar(p), sKpc);
		}
		PosVelSky toSky(const coord::PosVelCar pv, double& sKpc, double& Vlos_kms) const;
		PosVelSky toSky(const coord::PosVelCyl pv, double& sKpc, double& Vlos_kms) const{
			return toSky(coord::toPosVelCar(pv), sKpc, Vlos_kms);
		}
		VelSky toPM(const coord::PosVelCar p, double& Vlos_kms) const;
		VelSky toPM(const coord::PosVelCyl p, double& Vlos_kms) const{
			return toPM(coord::PosVelCar(coord::toPosVelCar(p)), Vlos_kms);
		}
		coord::PosCar xyz(void) const{// sun's location in internal units
			return coord::PosCar(Sun.x,Sun.y,Sun.z);
		}
		coord::VelCar Vxyz(void) const{// sun's velocity in internal units
			return coord::VelCar(Sun.vx,Sun.vy,Sun.vz);
		}
		double sKpc(const coord::PosCar p) const{// returns s in kpc from p in int units
			double S=pow_2(Sun.x-p.x)+pow_2(Sun.y-p.y)+pow_2(Sun.z-p.z);
			return (S>0? sqrt(S)/from_Kpc : 0);
		}
		double sKpc(const coord::PosCyl p) const{// returns s in kpc from p in int units
			return sKpc(coord::toPosCar(p));
		}
		double fromKpc() const {return from_Kpc;}
};

//You'll find obs::los in dust.h

}//namespace obs