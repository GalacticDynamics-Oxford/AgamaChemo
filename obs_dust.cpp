#pragma once
#include "math_ODE.h"
#include "obs.h"

namespace dust{

JTcloud::JTcloud(double _Rc, double _phic, double _z0, double _norm,const std::string& fname) :
    Rc(_Rc), phic(_phic), z0(_z0), norm(_norm) {
	if(phic>M_PI) phic-=2*M_PI;//move dscontinuity to M_PI
	FILE *ifile;
	if(fopen_s(&ifile,fname.c_str(),"r")){
		printf(" I can't open %s\n",fname.c_str());
		exit(0);
	}
	size_t ny,nx;
	float Q,amax,Xmax,Ymax;
	fscanf_s(ifile,"%zd %zd %f %f %f %f",&ny,&nx,&Q,&Xmax,&Ymax,&amax);
	std::vector<double> xs,ys;
	math::Matrix<double> array(nx,ny);
	for(int i=0;i<nx;i++) xs.push_back(Xmax*(2*i/(float)(nx-1)-1));
	for(int i=0;i<ny;i++) ys.push_back(Ymax*(2*i/(float)(nx-1)-1));
	float x,fac=norm/amax;
	for(int j=0; j<ny; j++){
		for(int i=0; i<nx; i++){
			fscanf_s(ifile,"%f",&x);
			array(i,j)=fac*x;
			if(std::isnan(array(i,j)))
				printf("error NAN in %s\n",fname.c_str());
		}
	}
	vals=math::LinearInterpolator2d(xs,ys,array);
}

EXP double JTcloud::dens(const coord::PosCyl& Rzphi) const{
	double phi=Rzphi.phi;
	if(fabs(phic)>0.5*M_PI && phi<0) phi+=2*M_PI;
	double x=(Rzphi.R/Rc-1), y=(phi-phic);//assume clockwise rotation
	if(fabs(x)>vals.xmax() || fabs(y)>vals.ymax()) return 0;
	return vals.value(x,y)*.5/z0*exp(-fabs(Rzphi.z)/z0);
}

EXP double BaseDustModel::dens(const coord::PosCyl& p) const{
	if(ptr)
		return rho0 * ptr->density(p);
	else
		return rho0 * exp(-p.R/Rd-fabs(p.z-Zw(p))/zd);
}

EXP double dustModel::dens(const coord::PosCyl& p) const{
	double log_rho=log(BaseDens(p));
	if(bl!=NULL){//add bubbles
		for(size_t i = 0; i < bl->size(); i++)
			log_rho += (*bl)[i].dens(p);
	}
	if(sp!=NULL){//ass spirals
		for(size_t i = 0; i < sp->size(); i++)
			log_rho += (*sp)[i].dens(p);
	}
	if(cl!=NULL){//add clouds
		for(size_t i = 0; i < cl->size(); i++)
			log_rho += (*cl)[i].dens(p);
	}
	return exp(log_rho);
}

}//namespace dust