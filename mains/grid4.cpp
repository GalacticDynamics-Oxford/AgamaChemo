//#include "/u/sm/mgo.h"
#include "/u/c/spress/hacked/press.h"
#include "grid2.h"
#include "grid4.h"
/* Grid centres are at (x0,y0), (x0+i*Deltax,y0+k*Deltay) etc so half
 * of edge cells lie outside plot. This ensures that contour plots are
 * correct because we are computing values on a regular grid that
 * touches frame*/

grid4::grid4(int _nw,int _nx,int _ny,int _nz,int _na,double wmin,double wmax,double xmin,double xmax,
	     double ymin,double ymax,double zmin,double zmax):
    nw(_nw), nx(_nx), ny(_ny), nz(_nz), na(_na), w0(wmin), x0(xmin), y0(ymin), z0(zmin) {
	greyit=false;
	Deltaw=(wmax-wmin)/(double)(nw-1);
	Deltax=(xmax-xmin)/(double)(nx-1);
	Deltay=(ymax-ymin)/(double)(ny-1);
	Deltaz=(zmax-zmin)/(double)(nz-1);
	m=fmatrix(nw,nx,ny,nz); av = na>0? fmatrix(nw,nx,ny,nz,na) : NULL;
	for(int h=0;h<nw;h++)
		for(int i=0;i<nx;i++)
			for(int j=0;j<ny;j++)
				for(int k=0;k<nz;k++){
					m[h][i][j][k]=0;
					for(int l=0;l<na;l++) av[h][i][j][k][l]=0;
				}
}
grid4::grid4(grid4& ogrid,bool copy) {//clone a grid4 & possible copy its data
	greyit=false;
	ogrid.disgorge(nw,nx,ny,nz,na,w0,x0,y0,z0,Deltaw,Deltax,Deltay,Deltaz);
	m=fmatrix(nw,nx,ny,nz); av = na>0? fmatrix(nw,nx,ny,nz,na) : NULL;
	if(copy){
		double* avs = na>0? new double[na] : NULL;
		for(int h=0;h<nw;h++){
			double w=w0+h*Deltaw;
			for(int i=0;i<nx;i++){
				double x=x0+i*Deltax;
				for(int j=0;j<ny;j++){
					double y=y0+j*Deltay;
					for(int k=0;k<nz;k++){
						double z=z0+k*Deltaz;
						m[h][i][j][k]=ogrid.dens(w,x,y,z,avs);
						for(int l=0;l<na;l++) av[h][i][j][k][l]=m[h][i][j][k]*avs[l];
					}
				}
			}
		}
		if(na>0) delete[] avs;
	}
}
/*
grid4::grid4(FILE *ifile,int na1) {
	greyit=false;
	if(10!=fscanf(ifile,"%d %d %d %d %lg %lg %lg %lg %lg %lg",&nx,&ny,&nz,&na,&x0,&y0,&z0,&Deltax,&Deltay,&Deltaz)){
		printf("grid4::grid4: something wrong with data file\n"); exit(0);
	}
	int na2=MAX(na1,na);
	m=fmatrix(nw,nx,ny,nz); av = na2>0? fmatrix(nw,nx,ny,nz,na2) : NULL;
	get(ifile,m,nw,nx,ny,nz);
	if(na>0) get(ifile,av,nw,nx,ny,nz,na);
	fclose(ifile);
	na=MAX(na1,na);
}*/
grid4::~grid4(){
	delmatrix(m,nw,nx,ny); if(na>0) delmatrix(av,nw,nx,ny,nz);
}
void grid4::get_centres(float *ws,float* xs,float* ys,float* zs){
	ws[0]=w0;
	for(int i=1;i<nw;i++) ws[i]=ws[i-1]+Deltaw;
	xs[0]=x0;
	for(int i=1;i<nx;i++) xs[i]=xs[i-1]+Deltax;
	ys[0]=y0;
	for(int i=1;i<ny;i++) ys[i]=ys[i-1]+Deltay;
	zs[0]=z0;
	for(int i=1;i<nz;i++) zs[i]=zs[i-1]+Deltaz;
}
void grid4::get_minmax(float& lmin,float& lmax) const{
	lmin=1e10; lmax=-1e10;
	for(int h=0;h<nw;h++)
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++)
			for(int k=0;k<nz;k++){
			if(m[h][i][j][k]<lmin) lmin=m[h][i][j][k];
			if(m[h][i][j][k]>lmax) lmax=m[h][i][j][k];
			}
}
void grid4::get_minmax(int l,float& lmin,float& lmax) const{
	lmin=1e10; lmax=-1e10;
	for(int h=0;h<nw;h++)
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++)
			for(int k=0;k<nz;k++){
				if(m[h][i][j][k]<=2) continue;
				float average=av[h][i][j][k][l]/m[h][i][j][k];
				if(average<lmin) lmin=average;
				if(average>lmax) lmax=average;
			}
}
void grid4::setup(int nw0,int nx0,int ny0,int nz0,int na0,
		  double wmin,double wmax,double xmin,double xmax,
		 double ymin,double ymax,double zmin,double zmax){
	nw=nw0; nx=nx0; ny=ny0; nz=nz0; na=na0;
	w0=wmin; x0=xmin; y0=ymin; z0=zmin; greyit=false;
	Deltaw=(wmax-wmin)/(double)(nw-1);
	Deltax=(xmax-xmin)/(double)(nx-1); Deltay=(ymax-ymin)/(double)(ny-1);
	Deltaz=(zmax-zmin)/(double)(nz-1);
	m=fmatrix(nw,nx,ny,nz); av = na>0? fmatrix(nw,nx,ny,nz,na) : NULL;
	for(int h=0;h<nw;h++)
		for(int i=0;i<nx;i++)
			for(int j=0;j<ny;j++)
				for(int k=0;k<nz;k++){
					m[h][i][j][k]=0;
					for(int l=0;l<na;l++) av[h][i][j][k][l]=0;
				}
}
void grid4::setup(grid4& grd){
	greyit=false;
	grd.disgorge(nw,nx,ny,nz,na,w0,x0,y0,z0,Deltaw,Deltax,Deltay,Deltaz);
//	printf("%d %d %d %d %d\n%f %f %f %f\n%f %f %f %f\n",nw,nx,ny,nz,na,w0,x0,y0,z0,Deltaw,Deltax,Deltay,Deltaz);
	m=fmatrix(nw,nx,ny,nz); av = na>0? fmatrix(nw,nx,ny,nz,na) : NULL;
	for(int h=0;h<nw;h++)
		for(int i=0;i<nx;i++)
			for(int j=0;j<ny;j++)
				for(int k=0;k<nz;k++){
					m[h][i][j][k]=0;
					for(int l=0;l<na;l++) av[h][i][j][k][l]=0;
				}
}
void grid4::disgorge(int& nw0,int& nx0,int& ny0,int& nz0,int& na0,
		     double& w00,double& x00,double& y00,double& z00,double& Deltaw0,
		     double& Deltax0,double& Deltay0,double& Deltaz0){
	nw0=nw, nx0=nx; ny0=ny; nz0=nz; na0=na; w00=w0; x00=x0; y00=y0; z00=z0;
	Deltaw0=Deltaw; Deltax0=Deltax; Deltay0=Deltay; Deltaz0=Deltaz;
}
void grid4::regrid(grid4& old,int in,int out){//put data from slot in into slot out
	int nap=old.navs();
	double* avs = new double[nap];
	for(int h=0;h<nw;h++){
		double w=w0+h*Deltaw;
	for(int i=0;i<nx;i++){
		double x=x0+i*Deltax;
		for(int j=0;j<ny;j++){
			double y=y0+j*Deltay;
			for(int k=0;k<nz;k++){
				double z=z0+k*Deltaz;
				m[h][i][j][k]=old.dens(w,x,y,z,avs);
				if(0!=m[h][i][j][k]) av[h][i][j][k][out]=m[h][i][j][k]*avs[in];
			}
		}
	}
	}
	delete[] avs;
}
/*double grid4::w(double x,double Dx){//w is fraction to home cell
	double f=1-fabs(x)/Dx;
	return MAX(0,f);
}*/
void grid4::wts(double w,double x,double y,double z,double& ww,double& wwp,
		double& wx,double& wxp,double& wy,double& wyp,double& wz,double& wzp){
	ww= MAX(0,1-fabs(w)/Deltaw); wwp=1-ww;
	wx= MAX(0,1-fabs(x)/Deltax); wxp=1-wx;
	wy= MAX(0,1-fabs(y)/Deltay); wyp=1-wy;
	wz= MAX(0,1-fabs(z)/Deltaz); wzp=1-wz;
}
void grid4::CIC(double weight,double w,double x,double y,double z){
	int h=(int)((w-w0)/Deltaw+.5);//index cell we're in
	int i=(int)((x-x0)/Deltax+.5);
	int j=(int)((y-y0)/Deltay+.5);
	int k=(int)((z-z0)/Deltaz+.5);
	double dw=w-(w0+((double)h)*Deltaw), dx=x-(x0+((double)i)*Deltax),
	dy=y-(y0+((double)j)*Deltay), dz=z-(z0+((double)k)*Deltaz);//ofsets from cell centre
	double ww[2],wx[2],wy[2],wz[2];
	wts(dw,dx,dy,dz,ww[0],ww[1],wx[0],wx[1],wy[0],wy[1],wz[0],wz[1]);
	int H= dw>0? h+1:h-1;
	int I= dx>0? i+1:i-1;
	int J= dy>0? j+1:j-1;
	int K= dz>0? k+1:k-1;
	int hw[2]={h,H},ix[2]={i,I},jy[2]={j,J},kz[2]={k,K};
	for(int d=0; d<2; d++)
		for(int a=0; a<2; a++)
			for(int b=0; b<2; b++)
				for(int c=0; c<2; c++){
					if(hw[d]>=0 && hw[d]<nw && ix[a]>=0 && ix[a]<nx && jy[b]>=0 && jy[b]<ny && kz[c]>=0 && kz[c]<nz)
						m[hw[d]][ix[a]][jy[b]][kz[c]]+=weight*ww[d]*wx[a]*wy[b]*wz[c];
				}
}
void grid4::CIC(double weight,double *value,double w,double x,double y,double z){
	int h=(int)((w-w0)/Deltaw+.5);//index cell we're in
	int i=(int)((x-x0)/Deltax+.5);
	int j=(int)((y-y0)/Deltay+.5);
	int k=(int)((z-z0)/Deltaz+.5);
	double dw=w-(w0+((double)h)*Deltaw), dx=x-(x0+((double)i)*Deltax),
	dy=y-(y0+((double)j)*Deltay), dz=z-(z0+((double)k)*Deltaz);//ofsets from cell centre
	double ww[2],wx[2],wy[2],wz[2];
	wts(dw,dx,dy,dz,ww[0],ww[1],wx[0],wx[1],wy[0],wy[1],wz[0],wz[1]);
	int H= dw>0? h+1:h-1;
	int I= dx>0? i+1:i-1;
	int J= dy>0? j+1:j-1;
	int K= dz>0? k+1:k-1;
	int hw[2]={h,H},ix[2]={i,I},jy[2]={j,J},kz[2]={k,K};
	for(int d=0; d<2; d++)
		for(int a=0; a<2; a++)
			for(int b=0; b<2; b++)
				for(int c=0; c<2; c++){
					if(hw[d]>=0 && hw[d]<nw && ix[a]>=0 && ix[a]<nx && jy[b]>=0 && jy[b]<ny && kz[c]>=0 && kz[c]<nz){
						double fac=ww[d]*wx[a]*wy[b]*wz[c];
						m[hw[d]][ix[a]][jy[b]][kz[c]]+=weight*fac;
						for(int l=0; l<na; l++) av[hw[d]][ix[a]][jy[b]][kz[c]][l]+=weight*value[l]*fac;
					}

				}
}
double grid4::dens(double w,double x,double y,double z){//Unnormalised: need to divide by cell volume
	//The virtual cell centres should be at (i)Delta
	if(w<w0 || x<x0 || y<y0 || z<z0) return 0;// off grid4
	int h=(int)((w-w0)/Deltaw+.5);//index cell we're in
	int i=(int)((x-x0)/Deltax+.5);//index cell we're in
	int j=(int)((y-y0)/Deltay+.5);
	int k=(int)((z-z0)/Deltaz+.5);
	if(h<0 || i<0 || j<0 || k<0 || h>= nw || i>=nx || j>=ny || k>=nz) return 0; //off grid
	double dw=w-(w0+((double)h)*Deltaw), dx=x-(x0+((double)i)*Deltax),
	dy=y-(y0+((double)j)*Deltay), dz=z-(z0+((double)k)*Deltaz);//ofsets from cell centre
	double ww[2],wx[2],wy[2],wz[2];
	wts(dw,dx,dy,dz,ww[0],ww[1],wx[0],wx[1],wy[0],wy[1],wz[0],wz[1]);
	int H= dw>0? h+1:h-1;
	int I= dx>0? i+1:i-1;
	int J= dy>0? j+1:j-1;
	int K= dz>0? k+1:k-1;
	int hw[2]={h,H},ix[2]={i,I},jy[2]={j,J},kz[2]={k,K};
	double mass=0;
	for(int d=0; d<2; d++)
		for(int a=0; a<2; a++)
			for(int b=0; b<2; b++)
				for(int c=0; c<2; c++){
					if(hw[d]>=0 && hw[d]<nw && ix[a]>=0 && ix[a]<nx && jy[b]>=0 && jy[b]<ny && kz[c]>=0 && kz[c]<nz){
						double fac=ww[d]*wx[a]*wy[b]*wz[c];
						mass+=m[hw[d]][ix[a]][jy[b]][kz[c]]*fac;
					}

				}
	return mass;
}
double grid4::dens(double w,double x,double y,double z,double *avs){//Unnormalised: need to divide by cell volume
	//The virtual cell centres should be at (i)Delta
	if(w<x0 || x<x0 || y<y0 || z<z0) return 0;// off grid4
	int h=(int)((w-w0)/Deltaw+.5);//index cell we're in
	int i=(int)((x-x0)/Deltax+.5);//index cell we're in
	int j=(int)((y-y0)/Deltay+.5);
	int k=(int)((z-z0)/Deltaz+.5);
	if(h<0 || i<0 || j<0 || k<0 || h>=nw || i>=nx || j>=ny || k>=nz) return 0; //off grid
	double dw=w-(w0+((double)h)*Deltaw), dx=x-(x0+((double)i)*Deltax),
	dy=y-(y0+((double)j)*Deltay), dz=z-(z0+((double)k)*Deltaz);//ofsets from cell centre
	double ww[2],wx[2],wy[2],wz[2];
	wts(dw,dx,dy,dz,wx[0],ww[0],ww[1],wx[1],wy[0],wy[1],wz[0],wz[1]);
	int H= dw>0? h+1:h-1;
	int I= dx>0? i+1:i-1;
	int J= dy>0? j+1:j-1;
	int K= dz>0? k+1:k-1;
	int hw[2]={h,H},ix[2]={i,I},jy[2]={j,J},kz[2]={k,K};
	double mass=0; for(int l=0;l<na;l++) avs[l]=0;
	for(int d=0; d<2; d++)
		for(int a=0; a<2; a++)
		for(int b=0; b<2; b++)
			for(int c=0; c<2; c++){
				if(hw[d]>=0 && hw[d]<nw && ix[a]>=0 && ix[a]<nx && jy[b]>=0 && jy[b]<ny && kz[c]>=0 && kz[c]<nz){
					double fac=ww[d]*wx[a]*wy[b]*wz[c];
					mass+=m[hw[d]][ix[a]][jy[b]][kz[c]]*fac;
					for(int l=0; l<na; l++) avs[l]+=av[hw[d]][ix[a]][jy[b]][kz[c]][l]*fac;
				}

			}
	return mass;
}
void grid4::get_values(int l,float blank,float ****mp) const{//returns av values with blank value of empty cells
	if(l<0){
		for(int h=0; h<nw; h++)
			for(int i=0;i<nx;i++)
			for(int j=0;j<ny;j++)
				for(int k=0; k<nz; k++) mp[h][i][j][k]=m[h][i][j][k];
	} else
		for(int h=0; h<nw; h++)
		for(int i=0;i<nx;i++){
			for(int j=0;j<ny;j++)
				for(int k=0;k<nz;k++){
				if(m[h][i][j][k]>2){
					mp[h][i][j][k]=av[h][i][j][k][l]/m[h][i][j][k];
				}else{
					mp[h][i][j][k]=blank;
				}
			}
		}
}

float grid4::get_total(void){
	float sum=0;
	for(int h=0; h<nw; h++)
	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
			for(int k=0; k<nz; k++)
				sum+=m[h][i][j][k];
	return sum;
}
/*
void grid4::dump(FILE *ofile){
	fprintf(ofile,"%d %d %d %d %d %g %g %g %g %g %g %g %g\n",nw,nx,ny,nz,na,w0,x0,y0,z0,Deltaw,Deltax,Deltay,Deltaz);
	compress(ofile,m,nw,nx,ny,nz);
	if(na>0) compress(ofile,av,nw,nx,ny,nz,na);
}*/
void grid4::zero(void){
	for(int h=0; h<nw; h++)
		for(int i=0; i<nx; i++)
			for(int j=0; j<ny; j++){
				for(int k=0; k<nz; k++){
					m[h][i][j][k]=0;
				}
			}
}
void grid4::divide(grid4 &grd){
	for(int h=0; h<nw; h++)
		for(int i=0; i<nx; i++)
			for(int j=0;j<ny;j++)
				for(int k=0; k<nz; k++){
					double M=grd.dens(w0+h*Deltaw,x0+i*Deltax,y0+j*Deltay,z0+k*Deltaz);
					if(fabs(M)>2) m[h][i][j][k]/=M;
				}
}
void grid4::divide(float fac){
	for(int h=0; h<nw; h++)
		for(int i=0;i<nx;i++)
			for(int j=0;j<ny;j++)
				for(int k=0; k<nz; k++){
					m[h][i][j][k]/=fac;
				}
}
void grid4::subtract(grid4 &grd){
	for(int h=0; h<nw; h++)
		for(int i=0; i<nx; i++)
			for(int j=0;j<ny;j++)
				for(int k=0; k<nz; k++){
					m[h][i][j][k]-=grd.dens(w0+h*Deltaw,x0+i*Deltax,y0+j*Deltay,z0+k*Deltaz);
				}
}
float grid4::reveal(int h,int i,int j,int k,float* avij){
	for(int l=0; l<na; l++) avij[l]=av[h][i][j][k][l];
	return m[h][i][j][k];
}
void grid4::add(grid4& grd){
	int iw,ix,iy,iz,ia;
	grd.get_ints(iw,ix,iy,iz,ia);
	if(iw!=nw || ix!=nx || iy!=ny || iz!=nz || ia!=na){
		printf("Error: illegal attempt to add incommensurable grid4s\n"); return;
	}
	float *avij=NULL;
	if(na>0) avij=new float[na];
	for(int h=0; h<nw; h++)
		for(int i=0; i<nx; i++){
			for(int j=0; j<ny; j++)
				for(int k=0; k<nz; k++){
					m[h][i][j][k]+=grd.reveal(h,i,j,k,avij);
					for(int l=0; l<na; l++) av[h][i][j][k][l]+=avij[l];
				}
		}
	if(na>0) delete avij;
}
void grid4::toggle_grey(void){
	greyit=(greyit+1)%2;
}

//Now routines involving mgo::plt
/*
void grid4::plot_slice(char c1,double h1,char c2,double h2,mgo::plt* pl,bool show,bool logit){
	if(!((c1=='w' || c1=='x' || c1=='y' || c1=='z') && (c2=='w' || c2=='x' || c2=='y' || c2=='z'))
	   || c1==c2){
		printf("ERROR: First 2 args %c %c of plot_slice must be one of w,x,y,z\n",c1,c2);
		return;
	}
	double x1,y1,Dx1,Dy1;
	int nx1,ny1;
	if((c1=='w' && c2=='x') || (c2=='w' && c1=='x')){
		nx1=ny; ny1=nz; x1=y0; y1=z0; Dx1=Deltay; Dy1=Deltaz;
	} else if((c1=='w' && c2=='y') || (c1=='y' && c2=='w')){
		nx1=nx; ny1=nz; x1=x0; y1=z0; Dx1=Deltax; Dy1=Deltaz;
	} else if((c1=='w' && c2=='z') || (c1=='z' && c2=='w')){
		nx1=nx; ny1=ny; x1=x0; y1=y0; Dx1=Deltax; Dy1=Deltay;
	} else if((c1=='x' && c2=='y') || (c1=='y' && c2=='x')){
		nx1=nw; ny1=nz; x1=w0; y1=z0; Dx1=Deltaw; Dy1=Deltaz;
	} else if((c1=='x' && c2=='z') || (c1=='z' && c2=='x')){
		nx1=nw; ny1=ny; x1=w0; y1=y0; Dx1=Deltaw; Dy1=Deltay;
	} else if((c1=='y' && c2=='z') || (c1=='z' && c2=='y')){
		nx1=nw; ny1=nx; x1=w0; y1=x0; Dx1=Deltaw; Dy1=Deltax;
	}
	grid2 grd2(2*nx1-1,2*ny1-1,0,x1,x1+(nx1-1)*Dx1,y1,y1+(ny1-1)*Dy1);
	for(int i=0; i<2*nx1-1; i++){
		double x=x1 + (double)i*Dx1*.5;
		for(int j=0; j<2*ny1-1; j++){
			double y=y1 + (double)j*Dy1*.5;
			if((c1=='w' && c2=='x'))
				grd2.CIC(dens(h1,h2,x,y),x,y);
			else if((c1=='x' && c2=='w'))
				grd2.CIC(dens(h2,h1,x,y),x,y);
			else if((c1=='w' && c2=='y'))
				grd2.CIC(dens(h1,x,h2,y),x,y);
			else if((c1=='y' && c2=='w'))
				grd2.CIC(dens(h2,x,h1,y),x,y);
			else if((c1=='w' && c2=='z'))
				grd2.CIC(dens(h1,x,y,h2),x,y);
			else if((c1=='z' && c2=='w'))
				grd2.CIC(dens(h2,x,y,h1),x,y);
			else if((c1=='x' && c2=='y'))
				grd2.CIC(dens(x,h1,h2,y),x,y);
			else if((c1=='y' && c2=='x'))
				grd2.CIC(dens(x,h2,h1,y),x,y);
			else if((c1=='x' && c2=='z'))
				grd2.CIC(dens(x,h1,y,h2),x,y);
			else if((c1=='z' && c2=='x'))
				grd2.CIC(dens(x,h2,y,h1),x,y);
			else if((c1=='y' && c2=='z'))
				grd2.CIC(dens(x,y,h1,h2),x,y);
			else if((c1=='z' && c2=='y'))
				grd2.CIC(dens(x,y,h2,h1),x,y);
			else printf("Error in grid4::plot_slice: %c %c\n",c1,c2);
		}
	}
	grd2.plot_dens(pl,show,logit);
}
void grid4::plot_projection(char c1,char c2,mgo::plt* pl,bool show,bool logit){
	if(!((c1=='w' || c1=='x' || c1=='y' || c1=='z') && (c2=='w' || c2=='x' || c2=='y' || c2=='z'))
	   || c1==c2){
		printf("ERROR: First 2 args %c %c of plot_projection must be one of w,x,y,z\n",c1,c2);
		return;
	}
	double x1,y1,Dx1,Dy1;
	int nx1,ny1;
	if((c1=='w' && c2=='x') || (c2=='w' && c1=='x')){
		nx1=ny; ny1=nz; x1=y0; y1=z0; Dx1=Deltay; Dy1=Deltaz;
	} else if((c1=='w' && c2=='y') || (c1=='y' && c2=='w')){
		nx1=nx; ny1=nz; x1=x0; y1=z0; Dx1=Deltax; Dy1=Deltaz;
	} else if((c1=='w' && c2=='z') || (c1=='z' && c2=='w')){
		nx1=nx; ny1=ny; x1=x0; y1=y0; Dx1=Deltax; Dy1=Deltay;
	} else if((c1=='x' && c2=='y') || (c1=='y' && c2=='x')){
		nx1=nw; ny1=nz; x1=w0; y1=z0; Dx1=Deltaw; Dy1=Deltaz;
	} else if((c1=='x' && c2=='z') || (c1=='z' && c2=='x')){
		nx1=nw; ny1=ny; x1=w0; y1=y0; Dx1=Deltaw; Dy1=Deltay;
	} else if((c1=='y' && c2=='z') || (c1=='z' && c2=='y')){
		nx1=nw; ny1=nx; x1=w0; y1=x0; Dx1=Deltaw; Dy1=Deltax;
	}
	grid2 grd2(2*nx1-1,2*ny1-1,0,x1,x1+(nx1-1)*Dx1,y1,y1+(ny1-1)*Dy1);
	for(int i=0; i<2*nx1-1; i++){
		float x=x1+i*Dx1*.5;
		for(int j=0; j<2*ny1-1; j++){
			float y=y1+j*Dy1*.5;
			if((c1=='w' && c2=='x') || (c2=='w' && c1=='x')){
				for(int s=0; s<nw; s++){
					float s1=w0+s*Deltaw;
					for(int t=0; t<nx; t++){
						float t1=x0+t*Deltax;
						grd2.CIC(dens(s1,t1,x,y),x,y);	
					}
				}
			}else if((c1=='w' && c2=='y') || (c1=='y' && c2=='w')){
				for(int s=0; s<nw; s++){
					float s1=w0+s*Deltaw;
					for(int t=0; t<ny; t++){
						float t1=y0+t*Deltay;
						grd2.CIC(dens(s1,x,t1,y),x,y);	
					}
				}
			}else if((c1=='w' && c2=='z') || (c1=='z' && c2=='w')){
				for(int s=0; s<nw; s++){
					float s1=w0+s*Deltaw;
					for(int t=0; t<nz; t++){
						float t1=z0+t*Deltaz;
						grd2.CIC(dens(s1,x,y,t1),x,y);	
					}
				}
			}else if((c1=='x' && c2=='y') || (c1=='y' && c2=='x')){
				for(int s=0; s<nx; s++){
					float s1=x0+s*Deltax;
					for(int t=0; t<ny; t++){
						float t1=y0+t*Deltay;
						grd2.CIC(dens(x,s1,t1,y),x,y);	
					}
				}
			}else if((c1=='x' && c2=='z') || (c1=='z' && c2=='x')){
				for(int s=0; s<nx; s++){
					float s1=x0+s*Deltax;
					for(int t=0; t<nz; t++){
						float t1=z0+t*Deltaz;
						grd2.CIC(dens(x,s1,y,t1),x,y);	
					}
				}
			}else if((c1=='y' && c2=='z') || (c1=='z' && c2=='y')){
				for(int s=0; s<ny; s++){
					float s1=y0+s*Deltay;
					for(int t=0; t<nz; t++){
						float t1=z0+t*Deltaz;
						grd2.CIC(dens(x,y,s1,t1),x,y);	
					}
				}
			}else printf("Error in grid4::plot_projection: %c %c\n",c1,c2);
		}
	}
	grd2.plot_dens(pl,show,logit);		
}
void grid4::plot_mean(char c1,char c2,double lmin,double lmax,mgo::plt* pl,bool show,bool logit){
	if(!((c1=='w' || c1=='x' || c1=='y' || c1=='z') && (c2=='w' || c2=='x' || c2=='y' || c2=='z'))
	   || c1==c2){
		printf("ERROR: First 2 args %c %c of plot_mean must be one of w,x,y,z\n",c1,c2);
		return;
	}
	double x1,y1,Dx1,Dy1;
	int nx1,ny1;
	if((c1=='w' && c2=='x') || (c2=='w' && c1=='x')){
		nx1=ny; ny1=nz; x1=y0; y1=z0; Dx1=Deltay; Dy1=Deltaz;
	} else if((c1=='w' && c2=='y') || (c1=='y' && c2=='w')){
		nx1=nx; ny1=nz; x1=x0; y1=z0; Dx1=Deltax; Dy1=Deltaz;
	} else if((c1=='w' && c2=='z') || (c1=='z' && c2=='w')){
		nx1=nx; ny1=ny; x1=x0; y1=y0; Dx1=Deltax; Dy1=Deltay;
	} else if((c1=='x' && c2=='y') || (c1=='y' && c2=='x')){
		nx1=nw; ny1=nz; x1=w0; y1=z0; Dx1=Deltaw; Dy1=Deltaz;
	} else if((c1=='x' && c2=='z') || (c1=='z' && c2=='x')){
		nx1=nw; ny1=ny; x1=w0; y1=y0; Dx1=Deltaw; Dy1=Deltay;
	} else if((c1=='y' && c2=='z') || (c1=='z' && c2=='y')){
		nx1=nw; ny1=nx; x1=w0; y1=x0; Dx1=Deltaw; Dy1=Deltax;
	}
	grid2 grd2(2*nx1-1,2*ny1-1,2,x1,x1+(nx1-1)*Dx1,y1,y1+(ny1-1)*Dy1);
	for(int i=0; i<2*nx1-1; i++){
		double x=x1+i*Dx1*.5;
		for(int j=0; j<2*ny1-1; j++){
			double y=y1+j*Dy1*.5;
			if(c1=='w' && c2=='x'){//average of w summed on x
				for(int s=0; s<nw; s++){
					double s1=w0+s*Deltaw;
					for(int t=0; t<nx; t++){
						double t1=x0+t*Deltax;
						double rho=dens(s1,t1,x,y);
						double v[2]={s1,t1};
						grd2.CIC(rho,v,x,y);
					}
				}
			} else if(c1=='x' && c2=='w'){//average of x summed on w
				for(int s=0; s<nx; s++){
					double s1=x0+s*Deltax;
					for(int t=0; t<nw; t++){
						double t1=w0+t*Deltaw;
						double rho=dens(t1,s1,x,y);
						double v[2]={s1,t1};
						grd2.CIC(rho,v,x,y);
					}
				}
			} else if(c1=='w' && c2=='y'){//average of w summed on y
				for(int s=0; s<nw; s++){
					double s1=w0+s*Deltaw;
					for(int t=0; t<ny; t++){
						double t1=y0+t*Deltay;
						double rho=dens(s1,x,t1,y);
						double v[2]={s1,t1};
						grd2.CIC(rho,v,x,y);
					}
				}
			} else if(c1=='y' && c2=='w'){//average  of y summed on w
				for(int s=0; s<ny; s++){
					double s1=y0+s*Deltay;
					for(int t=0; t<nw; t++){
						double t1=w0+t*Deltaw;
						double rho=dens(t1,x,s1,y);
						double v[2]={s1,t1};
						grd2.CIC(rho,v,x,y);
					}
				}
			} else if(c2=='w' && c1=='z'){//average of w sumed on z
				for(int s=0; s<nw; s++){
					double s1=w0+s*Deltaw;
					for(int t=0; t<nz; t++){
						double t1=z0+t*Deltaz;
						double rho=dens(s1,x,y,t1);
						double v[2]={s1,t1};
						grd2.CIC(rho,v,x,y);
					}
				}
			} else if(c1=='z' && c2=='w'){//average of z summed on w
				for(int s=0; s<nz; s++){
					double s1=z0+s*Deltaz;
					for(int t=0; t<nw; t++){
						double t1=w0+t*Deltaw;
						double rho=dens(t1,x,y,s1);
						double v[2]={s1,t1};
						grd2.CIC(rho,v,x,y);
					}
				}
			} else if(c1=='x' && c2=='y'){//average of x summed on y
				for(int s=0; s<nx; s++){
					double s1=x0+s*Deltax;
					for(int t=0; t<ny; t++){
						double t1=y0+t*Deltay;
						double rho=dens(x,s1,t1,y);
						double v[2]={s1,t1};
						grd2.CIC(rho,v,x,y);
					}
				}
			} else if(c1=='y' && c2=='x'){//average of y summed on x
				for(int s=0; s<ny; s++){
					double s1=y0+s*Deltay;
					for(int t=0; t<nx; t++){
						double t1=x0+t*Deltax;
						double rho=dens(x,t1,s1,y);
						double v[2]={s1,t1};
						grd2.CIC(rho,v,x,y);
					}
				}
			} else if(c1=='x' && c2=='z'){//average of x summed on z
				for(int s=0; s<nx; s++){
					double s1=x0+s*Deltax;
					for(int t=0; t<nz; t++){
						double t1=z0+t*Deltaz;
						double rho=dens(x,s1,y,t1);
						double v[2]={s1,t1};
						grd2.CIC(rho,v,x,y);
					}
				}
			} else if(c1=='z' && c2=='x'){//average of z summed on x
				for(int s=0; s<nz; s++){
					double s1=z0+s*Deltaz;
					for(int t=0; t<nx; t++){
						double t1=x0+t*Deltax;
						double rho=dens(x,t1,y,s1);
						double v[2]={s1,t1};
						grd2.CIC(rho,v,x,y);
					}
				}
			} else if(c1=='y' && c2=='z'){// average of y summed on z
				for(int s=0; s<ny; s++){
					double s1=y0+s*Deltay;
					for(int t=0; t<nz; t++){
						double t1=z0+t*Deltaz;
						double rho=dens(x,y,s1,t1);
						double v[2]={s1,t1};
						grd2.CIC(rho,v,x,y);
					}
				}
			} else if(c1=='z' && c2=='y'){//average of z summed on y
				for(int s=0; s<nz; s++){
					double s1=z0+s*Deltaz;
					for(int t=0; t<ny; t++){
						double t1=y0+t*Deltay;
						double rho=dens(x,y,t1,s1);
						double v[2]={s1,t1};
						grd2.CIC(rho,v,x,y);
					}
				}
			} else	printf("Error in grid4::plot_mean: %c %c\n",c1,c2);
		}
	}
	float blank= logit? -10 : -10;
	if(lmin==NAN)
		grd2.plot_values(0,blank,pl,show,logit);
	else
		grd2.plot_values(lmin,lmax,0,blank,pl,show,logit);
}
void grid4::plot_meanR(char c1,char c2,double lmin,double lmax,mgo::plt* pl,bool show,bool logit){
	if(!((c1=='w' || c1=='x' || c1=='y' || c1=='z') && (c2=='w' || c2=='x' || c2=='y' || c2=='z'))
	   || c1==c2){
		printf("ERROR: First 2 args %c %c of plot_mean must be one of w,x,y,z\n",c1,c2);
		return;
	}
	double x1,y1,Dx1,Dy1;
	int nx1,ny1;
	if((c1=='w' && c2=='x') || (c2=='w' && c1=='x')){
		nx1=ny; ny1=nz; x1=y0; y1=z0; Dx1=Deltay; Dy1=Deltaz;
	} else if((c1=='w' && c2=='y') || (c1=='y' && c2=='w')){
		nx1=nx; ny1=nz; x1=x0; y1=z0; Dx1=Deltax; Dy1=Deltaz;
	} else if((c1=='w' && c2=='z') || (c1=='z' && c2=='w')){
		nx1=nx; ny1=ny; x1=x0; y1=y0; Dx1=Deltax; Dy1=Deltay;
	} else if((c1=='x' && c2=='y') || (c1=='y' && c2=='x')){
		nx1=nw; ny1=nz; x1=w0; y1=z0; Dx1=Deltaw; Dy1=Deltaz;
	} else if((c1=='x' && c2=='z') || (c1=='z' && c2=='x')){
		nx1=nw; ny1=ny; x1=w0; y1=y0; Dx1=Deltaw; Dy1=Deltay;
	} else if((c1=='y' && c2=='z') || (c1=='z' && c2=='y')){
		nx1=nw; ny1=nx; x1=w0; y1=x0; Dx1=Deltaw; Dy1=Deltax;
	}
	grid2 grd2(2*nx1-1,2*ny1-1,2,x1,x1+(nx1-1)*Dx1,y1,y1+(ny1-1)*Dy1);
	for(int i=0; i<2*nx1-1; i++){
		double x=x1+i*Dx1*.5;
		for(int j=0; j<2*ny1-1; j++){
			double y=y1+j*Dy1*.5;
			if(c1=='w' && c2=='x'){//average of w summed on x
				for(int s=0; s<nw; s++){
					double s1=w0+s*Deltaw;
					for(int t=0; t<nx; t++){
						double t1=x0+t*Deltax;
						double rho=dens(s1,t1,x,y);
						double v[2]={s1,t1};
						grd2.CIC(rho,v,x,y);
					}
				}
			} else if(c1=='x' && c2=='w'){//average of x summed on w
				for(int s=0; s<nx; s++){
					double s1=x0+s*Deltax;
					for(int t=0; t<nw; t++){
						double t1=w0+t*Deltaw;
						double rho=dens(t1,s1,x,y);
						double v[2]={s1,t1};
						grd2.CIC(rho,v,x,y);
					}
				}
			} else if(c1=='w' && c2=='y'){//average of w summed on y
				for(int s=0; s<nw; s++){
					double s1=w0+s*Deltaw;
					for(int t=0; t<ny; t++){
						double t1=y0+t*Deltay;
						double rho=dens(s1,x,t1,y);
						double v[2]={s1,t1};
						grd2.CIC(rho,v,x,y);
					}
				}
			} else if(c1=='y' && c2=='w'){//average  of y summed on w
				for(int s=0; s<ny; s++){
					double s1=y0+s*Deltay;
					for(int t=0; t<nw; t++){
						double t1=w0+t*Deltaw;
						double rho=dens(t1,x,s1,y);
						double v[2]={s1,t1};
						grd2.CIC(rho,v,x,y);
					}
				}
			} else if(c2=='w' && c1=='z'){//average of w sumed on z
				for(int s=0; s<nw; s++){
					double s1=w0+s*Deltaw;
					for(int t=0; t<nz; t++){
						double t1=z0+t*Deltaz;
						double rho=dens(s1,x,y,t1);
						double v[2]={s1,t1};
						grd2.CIC(rho,v,x,y);
					}
				}
			} else if(c1=='z' && c2=='w'){//average of z summed on w
				for(int s=0; s<nz; s++){
					double s1=z0+s*Deltaz;
					for(int t=0; t<nw; t++){
						double t1=w0+t*Deltaw;
						double rho=dens(t1,x,y,s1);
						double v[2]={s1,t1};
						grd2.CIC(rho,v,x,y);
					}
				}
			} else if(c1=='x' && c2=='y'){//average of x summed on y
				for(int s=0; s<nx; s++){
					double s1=x0+s*Deltax;
					for(int t=0; t<ny; t++){
						double t1=y0+t*Deltay;
						double rho=dens(x,s1,t1,y);
						double v[2]={s1,t1};
						grd2.CIC(rho,v,x,y);
					}
				}
			} else if(c1=='y' && c2=='x'){//average of y summed on x
				for(int s=0; s<ny; s++){
					double s1=y0+s*Deltay;
					for(int t=0; t<nx; t++){
						double t1=x0+t*Deltax;
						double rho=dens(x,t1,s1,y);
						double v[2]={s1,t1};
						grd2.CIC(rho,v,x,y);
					}
				}
			} else if(c1=='x' && c2=='z'){//average of x summed on z
				for(int s=0; s<nx; s++){
					double s1=x0+s*Deltax;
					for(int t=0; t<nz; t++){
						double t1=z0+t*Deltaz;
						double rho=dens(x,s1,y,t1);
						double v[2]={s1,t1};
						grd2.CIC(rho,v,x,y);
					}
				}
			} else if(c1=='z' && c2=='x'){//average of z summed on x
				for(int s=0; s<nz; s++){
					double s1=z0+s*Deltaz;
					for(int t=0; t<nx; t++){
						double t1=x0+t*Deltax;
						double rho=dens(x,t1,y,s1);
						double v[2]={s1,t1};
						grd2.CIC(rho,v,x,y);
					}
				}
			} else if(c1=='y' && c2=='z'){// average of y summed on z
				for(int s=0; s<ny; s++){
					double s1=y0+s*Deltay;
					for(int t=0; t<nz; t++){
						double t1=z0+t*Deltaz;
						double rho=dens(x,y,s1,t1);
						double v[2]={s1,t1};
						grd2.CIC(rho,v,x,y);
					}
				}
			} else if(c1=='z' && c2=='y'){//average of z summed on y
				for(int s=0; s<nz; s++){
					double s1=z0+s*Deltaz;
					for(int t=0; t<ny; t++){
						double t1=y0+t*Deltay;
						double rho=dens(x,y,t1,s1);
						double v[2]={s1,t1};
						grd2.CIC(rho,v,x,y);
					}
				}
			} else	printf("Error in grid4::plot_mean: %c %c\n",c1,c2);
		}
	}
	float blank= logit? -10 : -10;
	if(lmin==NAN)
		grd2.plot_valuesR(0,blank,pl,show,logit);
	else
		grd2.plot_valuesR(lmin,lmax,0,blank,pl,show,logit);
}
void grid4::plot_plaque(double lmin,double lmax,double ymin,double ymax,double zmin,double zmax,mgo::plt* pl,bool show,bool logit){
	grid2 grd2(2*nw-1,2*nx-1,0,w0,w0+(nw-1)*Deltaw,x0,x0+(nx-1)*Deltax);
	for(int i=0;i<2*ny-1;i++){
		double y=y0+i*Deltay*.5;
		if(y<ymin || y>ymax) continue;
		for(int j=0;j<2*nz-1;j++){
			double z=z0+j*Deltaz*.5;
			if(z<zmin || z>zmax) continue;
			for(int s=0; s<2*nw-1; s++){
				double s1=w0+s*Deltaw*.5;
				for(int t=0; t<2*nx-1; t++){
					double t1=x0+t*Deltax*.5;
					grd2.CIC(dens(s1,t1,y,z),s1,t1);
				}
			}
		}
	}
//	printf("grid2.total: %f\n",grd2.get_total());
	if(lmin==NAN)
		grd2.plot_dens(pl,show,logit);
	else
		grd2.plot_dens(lmin,lmax,pl,show,logit);
}
*/