/** \file    obs_dust.h
    \brief   Defines BaseDustModel and DustModel
    The BaseDustModel includes a warp but is otherwise v simple
    A DustModel is intended for the MW and can include bubbles, spirals
    and Julian-Toomre clouds
    \author  James Binney
    \date    2021-2022
*/

#pragma once
#include "potential_factory.h"

namespace dust {

class EXP BaseDustModel{
	private:
		const double Rd, zd, rho0, Rw, Hw;
		potential::PtrDensity ptr;
		double Zw(const coord::PosCyl& p) const{
			if(p.R<Rw) return 0;
			else return Hw*(p.R-Rw)/Rw*sin(p.phi);
		}
	public:
		const double from_Kpc;
		BaseDustModel(const double _Rd, const double _zd, const double _rho0,
			       const double _Rw, const double _Hw, const double _fromKpc):
		    Rd(_Rd*_fromKpc), zd(_zd*_fromKpc), rho0(_rho0),
		Rw(_Rw*_fromKpc), Hw(_Hw*_fromKpc), from_Kpc(_fromKpc) {}
		BaseDustModel(potential::PtrDensity _ptr, const double _rho0,
			       const double _fromKpc): Rd(0), zd(0), Rw(0), Hw(0),
		    ptr(_ptr), rho0(_rho0), from_Kpc(_fromKpc) {}
		virtual double dens(const coord::PosCyl& p) const;
		double BaseDens(const coord::PosCyl& p) const{
			return dens(p);
		}
};

/* Class introducing a global (logarithimc) spiral pattarn. Ampl peak
   * density, phase peak phi, alpha radial wavenumber, 1/kz (exponential) scale height.
   * Distances in internal units
 */
class EXP Spiral{
	double Ampl,phase,alpha,Narms,kz;
	public:
		Spiral(){}
		Spiral(double _Ampl,double _phase,double _alpha,int _Narms,double _kz) :
		    Ampl(_Ampl), phase(_phase), alpha(_alpha), Narms(_Narms), kz(_kz){}			
		double dens(coord::PosCyl Rzphi){
			return Ampl*cos(alpha*log(Rzphi.R+1.e-10)-Narms*(Rzphi.phi-phase))*exp(-kz*fabs(Rzphi.z));
		}
};

/* Class creating spiral-shaped dust clouds. (Rc,phic) Galactocentric
 * coords of cloud centre, z0 (exponential) scale height - all in
 * intUnits - norm central density, fname file from which to read surface
 * density computed per Julian & Toomre (1966)
*/
class EXP JTcloud{
	math::LinearInterpolator2d vals;
	double Rc,phic,z0,norm;
	public:
		JTcloud(){}
		JTcloud(double _Rc, double _phic, double _z0, double _norm,const std::string& fname);
		double dens(const coord::PosCyl& Rzphi) const;
};

/* Class creating blobs at p
*/
class EXP Blob{
	coord::PosCyl p;
	double Ampl,rad;
	public:
		Blob(){}
		Blob(const coord::PosCyl _p, const double _Ampl, const double _rad) :
		    p(_p), Ampl(_Ampl), rad(_rad){//ensure continuity around 180
			if(fabs(p.phi)>0.5*M_PI && p.phi<0) p.phi+=2*M_PI;
		}
		double dens(const coord::PosCyl& Rzphi) const{
			double phi=Rzphi.phi;
			if(fabs(p.phi)>0.5*M_PI && phi<0) phi+=2*M_PI;
			double x=(Rzphi.R-p.R), y=p.R*(phi-p.phi);//assume clockwise rotation
			//if(fabs(x)>3*rad || fabs(y)>3*rad) return 0;
			return Ampl * exp(-.5*(pow_2(x/rad)+pow_2(y/rad)+pow_2((Rzphi.z-p.z)/rad)));
		}
};

class EXP dustModel: public BaseDustModel {
	private:
		obs::solarShifter* shifter;
		potential::PtrDensity ptr;
		std::vector<JTcloud>* cl;
		std::vector<Spiral>* sp;
		std::vector<Blob>* bl;
	public:
		dustModel(const double _Rd, const double _zd, const double dAvds,
			  const double _Rw, const double _Hw,
			  obs::solarShifter* _shifter):
		    BaseDustModel(_Rd,_zd,dAvds/_shifter->from_Kpc/dens(coord::toPosCyl(_shifter->xyz())),
		    _Rw, _Hw, _shifter->from_Kpc),
		    shifter(_shifter){
			ptr = NULL; sp = NULL; cl = NULL;
		}
		dustModel(potential::PtrDensity _ptr, const double dAvds,
			  obs::solarShifter* _shifter):
		    ptr(_ptr), shifter(_shifter),
		    BaseDustModel(_ptr, dAvds/_shifter->from_Kpc/ptr->density(coord::toPosCyl(_shifter->xyz())),
		    _shifter->from_Kpc)
 {
			bl = NULL; sp = NULL; cl = NULL;
		}
		~dustModel(void){
			if(bl!=NULL) bl->~vector();
			if(sp!=NULL) sp->~vector();
			if(cl!=NULL) cl->~vector();
		}
		double dens(const coord::PosCyl& p) const;
		void addBlob(const coord::PosCyl p,double Ampl,double rad){
			Blob nbl(p,Ampl,rad);
			if(bl==NULL){
				std::vector<Blob>* blobs = new std::vector<Blob>;
				bl=blobs;
			}
			bl->push_back(nbl);
		}
		void deleteBlob(int n=-1){
			if(bl!=NULL){
				if(n==-1) bl->erase(bl->end());
				else bl->erase(bl->begin()+n);
			}
		}
		void addSpiral(double Ampl,double phase,double alpha,int Narms,double kz){
			Spiral nsp(Ampl,phase,alpha,Narms,kz);
			if(sp==NULL){
				std::vector<Spiral>* spirals = new std::vector<Spiral>;
				sp=spirals;
			}
			sp->push_back(nsp);
		}
		void deleteSpiral(int n=-1){
			if(sp!=NULL){
				if(n==-1) sp->erase(sp->end());
				else sp->erase(sp->begin()+n);
			}
		}
		void addCloud(double Rc, double phic, double z0, double norm,const std::string& fname){
			JTcloud ncl(Rc,phic,z0,norm,fname);
			if(cl==NULL){
				std::vector<JTcloud>* clouds = new std::vector<JTcloud>;
				cl=clouds;
			}
			cl->push_back(ncl);
		}
		void deleteCloud(int n=-1){
			if(cl!=NULL){
				if(n==-1) cl->erase(cl->end());
				else cl->erase(cl->begin()+n);
			}
		}
};

} // namespace dust
