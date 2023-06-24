/** \file    obs_los.h
    \brief   Defines lines of sight
    \author  James Binney
    \date    2022
*/
#pragma once

namespace obs {

class EXP BaseLos{
	private:
		math::CubicSpline extinct;
		dust::BaseDustModel* dm;
	public:
		double from_Kpc, s0;
		bool semi_inf;
		BaseLos(){};
		virtual ~BaseLos(){};
		BaseLos(const double _fromKpc, const bool _semi_inf, const double _s0,
			  dust::BaseDustModel* _dm=NULL):
		    from_Kpc(_fromKpc), semi_inf(_semi_inf), s0(_s0), dm(_dm) {
		}
		virtual coord::PosCyl deepest() const{
			printf("Erroneous call to deepest()");
			return coord::PosCyl(0,0,0);
		};
		virtual coord::PosCar xyz(const double s) const{
			printf("Error: xyz call in BaseLos\n");
			return coord::PosCar(0,0,0);
		}
		coord::PosCyl Rzphi(const double s) const{// s in internal units
			coord::PosCyl Cyl(coord::toPosCyl(xyz(s)));			
			return Cyl;
		}
		virtual double s(const coord::PosCar& xyz) const{// s in internal units
			printf("Error: s call in BaseLos\n");
			return 0;
		}
		double s(const coord::PosCyl& p) const{
			return s(coord::toPosCar(p));
		}
		double sKpc(const coord::PosCar& xyz) const{// s in Kpc
			return s(xyz)/from_Kpc;
		}
		double sKpc(const coord::PosCyl& p) const{
			return sKpc(coord::toPosCar(p));
		}
		virtual std::pair<double,double> sVlos(const coord::PosVelCyl& xv) const{// s Vlos int units
			return std::make_pair(0,0);
		}
		double A_V(const double sKpc) const{
			if(!dm) return 0;
			else{
				double A,dA;
				extinct.evalDeriv(sKpc,&A,&dA);
				return A;
			}
		}
		double A_B(const double sKpc) const{
			return 1.324 * A_V(sKpc);
		}
		double A_R(const double sKpc) const{
			return 0.748 * A_V(sKpc);
		}
		double A_H(const double sKpc) const{
			return 0.175 * A_V(sKpc);
		}
		double A_K(const double sKpc) const{
			return 0.112 * A_V(sKpc);
		}
		virtual double start() const{
			printf("Erroneous call to start\n");
			return 0;
		}
		virtual void tab_extinct(void);
};

/*
 *class for lines of sight through external galaxies
 *(xs,ys) coords (in kpc) wrt apparent maj/min axies
 *incl inclination (deg), s0 distance to galaxy (kpc)
 */
class EXP extLos: public obs::BaseLos {
	public: 
		const double xs, ys, incl, cosi, sini;
		extLos(double _xs, double _ys, double _incl, double s0,
		       double _fromKpc, dust::BaseDustModel* _dm=NULL);
		virtual coord::PosCyl deepest() const;// Point of closest approach to GC
		virtual coord::PosCar xyz(const double s) const{// s in internal units
			double x = xs;//s is measured from crossing of xz plane
			double y = s*sini;
			double z = ys/sini+s*cosi;
			return coord::PosCar(x, y, z);
		}
		virtual double s(const coord::PosCar& _xyz) const{// s in internal units
			double s = pow_2(_xyz.y) + pow_2(_xyz.z-ys/sini);
			s = s>0?  sqrt(s) : 0;
			return _xyz.y>0? s : -s;
		}
		virtual std::pair<double,double> sVlos(const coord::PosVelCyl& p) const{
			coord::PosVelCar xv(coord::toPosVelCar(p));
			return std::make_pair(s(xv), xv.vy*sini + xv.vz*cosi);
		}
		virtual double start() const{//starting s (Kpc) value for Av integration 
			return fabs(sini)>0? (1-ys/from_Kpc/sini)/cosi : -30;
		}
};

//class line of sight from Sun, specified by (l,b) and location of Sun.
//Extinction in BVRHK down it will be computed if a dustModel isspecified;
//otherwise zero extinction
class EXP SunLos: public obs::BaseLos {
	public: 
		const obs::solarShifter sun;
		PosSky pos;// l, b in radians
		double cosl, cosb, sinl, sinb;
		SunLos(const PosSky _pos, const obs::solarShifter _sun,
		    dust::dustModel* _dm=NULL);//angles in radians
		SunLos(const coord::PosCar xyz, const obs::solarShifter _sun,
		    dust::dustModel* _dm=NULL);
		virtual coord::PosCyl deepest() const;// Point of closest approach to GC
		virtual coord::PosCar xyz(const double s) const{// s in internal units
			double x = sun.xyz().x + s*cosb*cosl;
			double y = sun.xyz().y + s*cosb*sinl;
			double z = sun.xyz().z + s*sinb;
			coord::PosCar Car(x, y, z);
			return Car;
		}
		virtual double s(const coord::PosCar& xyz) const{// s in internal units
			double s = pow_2(xyz.x-sun.xyz().x) + pow_2(xyz.y-sun.xyz().y) + pow_2(xyz.z-sun.xyz().z);
			return (s>0?  sqrt(s) : 0);
		}
		virtual double start() const{
			return 0;
		}
};

}// namespace obs