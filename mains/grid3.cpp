//#include "/u/sm/mgo.h"
#include "/u/c/spress/hacked/press.h"
#include "grid2.h"
#include "grid3.h"
/* Grid centres are at (x0,y0), (x0+i*Deltax,y0+k*Deltay) etc so half
 * of edge cells lie outside plot. This ensures that contour plots are
 * correct because we are computing values on a regular grid that
 * touches frame*/

grid3::grid3(int _nx,int _ny,int _nz,int _na,double xmin,double xmax,
	     double ymin,double ymax,double zmin,double zmax):
    nx(_nx), ny(_ny), nz(_nz), na(_na), x0(xmin), y0(ymin), z0(zmin) {
	greyit=false;
	Deltax=(xmax-xmin)/(double)(nx-1); Deltay=(ymax-ymin)/(double)(ny-1);
	Deltaz=(zmax-zmin)/(double)(nz-1);
	m=fmatrix(nx,ny,nz); av = na>0? fmatrix(nx,ny,nz,na) : NULL;
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++)
			for(int k=0;k<nz;k++){
			m[i][j][k]=0;
			for(int l=0;l<na;l++) av[i][j][k][l]=0;
			}
	for(int i=0; i<3; i++){
		means[i]=0; vars[i]=0;
	}
}
grid3::grid3(grid3& ogrid,bool copy) {//copy a grid3 & its data
	greyit=false;
	ogrid.disgorge(nx,ny,nz,na,x0,y0,z0,Deltax,Deltay,Deltaz);
	m=fmatrix(nx,ny,nz); av = na>0? fmatrix(nx,ny,nz,na) : NULL;
	if(copy){
		double* avs = na>0? new double[na] : NULL;
		for(int i=0;i<nx;i++){
			double x=x0+i*Deltax;
			for(int j=0;j<ny;j++)
			{
				double y=y0+j*Deltay;
				for(int k=0;k<nz;k++){
					double z=z0+k*Deltaz;
					m[i][j][k]=ogrid.dens(x,y,z,avs);
					for(int l=0;l<na;l++) av[i][j][k][l]=m[i][j][k]*avs[l];
				}
			}
		}
		if(na>0) delete[] avs;
	}
}
/*
grid3::grid3(FILE *ifile,int na1){
	greyit=false;
	if(10!=fscanf(ifile,"%d %d %d %d %lg %lg %lg %lg %lg %lg",&nx,&ny,&nz,&na,&x0,&y0,&z0,&Deltax,&Deltay,&Deltaz)){
		printf("grid3::grid3: something wrong with data file\n"); exit(0);
	}
	int na2=MAX(na1,na);
	m=fmatrix(nx,ny,nz); av = na2>0? fmatrix(nx,ny,nz,na2) : NULL;
	get(ifile,m,nx,ny,nz);
	if(na>0) get(ifile,av,nx,ny,nz,na);
	fclose(ifile);
	na=MAX(na1,na);
}*/
grid3::~grid3(){
	delmatrix(m,nx,ny); if(na>0) delmatrix(av,nx,ny,nz);
}
void grid3::get_centres(float* xs,float* ys,float* zs){
	xs[0]=x0;
	for(int i=1;i<nx;i++) xs[i]=xs[i-1]+Deltax;
	ys[0]=y0;
	for(int i=1;i<ny;i++) ys[i]=ys[i-1]+Deltay;
	zs[0]=z0;
	for(int i=1;i<nz;i++) zs[i]=zs[i-1]+Deltaz;
}
void grid3::get_minmax(float& lmin,float& lmax) const{
	lmin=1e10; lmax=-1e10;
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++)
			for(int k=0;k<nz;k++){
			if(m[i][j][k]<lmin) lmin=m[i][j][k];
			if(m[i][j][k]>lmax) lmax=m[i][j][k];
			}
}
void grid3::get_minmax(int l,float& lmin,float& lmax) const{
	lmin=1e10; lmax=-1e10;
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++)
			for(int k=0;k<nz;k++){
				if(m[i][j][k]<=2) continue;
				float average=av[i][j][k][l]/m[i][j][k];
				if(average<lmin) lmin=average;
				if(average>lmax) lmax=average;
			}
}
void grid3::setup(int nx0,int ny0,int nz0,int na0,double xmin,double xmax,
		 double ymin,double ymax,double zmin,double zmax){
	nx=nx0; ny=ny0; nz=nz0; na=na0;
	x0=xmin; y0=ymin; z0=zmin; greyit=false;
	Deltax=(xmax-xmin)/(double)(nx-1); Deltay=(ymax-ymin)/(double)(ny-1);
	Deltaz=(zmax-zmin)/(double)(nz-1);
	m=fmatrix(nx,ny,nz); av = na>0? fmatrix(nx,ny,nz,na) : NULL;
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++)
			for(int k=0;k<nz;k++){
				m[i][j][k]=0;
				for(int l=0;l<na;l++) av[i][j][k][l]=0;
			}
	for(int i=0; i<3; i++){
		means[i]=0; vars[i]=0;
	}
}
void grid3::disgorge(int &nx0,int &ny0,int &nz0,int &na0,
	       double &x00,double &y00,double &z00,double &Deltax0,double &Deltay0,double &Deltaz0){
	nx0=nx; ny0=ny; nz0=nz; na0=na; x00=x0; y00=y0; z00=z0;
	Deltax0=Deltax; Deltay0=Deltay; Deltaz0=Deltaz;
}
void grid3::regrid(grid3 &old,int in,int out){//put data from slot in into slot out
	int nap=old.navs();
	double* avs = new double[nap];
	for(int i=0;i<nx;i++){
		double x=x0+i*Deltax;
		for(int j=0;j<ny;j++){
			double y=y0+j*Deltay;
			for(int k=0;k<nz;k++){
				double z=z0+k*Deltaz;
				m[i][j][k]=old.dens(x,y,z,avs);
				if(0!=m[i][j][k]) av[i][j][k][out]=m[i][j][k]*avs[in];
			}
		}
	}
	delete[] avs;
}
double grid3::w(double x,double Dx){//w is fraction to home cell
	double f=1-fabs(x)/Dx;
	return MAX(0,f);
}
void grid3::wts(double x,double y,double z,double &wx,double &wxp,
	       double &wy,double &wyp,double &wz,double &wzp){
	wx= MAX(0,1-fabs(x)/Deltax); wxp=1-wx;
	wy= MAX(0,1-fabs(y)/Deltay); wyp=1-wy;
	wz= MAX(0,1-fabs(z)/Deltaz); wzp=1-wz;
}
void grid3::CIC(double weight,double x,double y,double z){
	int i=(int)((x-x0)/Deltax+.5);//index cell we're in
	int j=(int)((y-y0)/Deltay+.5);//
	int k=(int)((z-z0)/Deltaz+.5);//
	double dx=x-(x0+((double)i)*Deltax), dy=y-(y0+((double)j)*Deltay),
	dz=z-(z0+((double)k)*Deltaz);//ofsets from cell centre
	double wx[2],wy[2],wz[2];
	wts(dx,dy,dz,wx[0],wx[1],wy[0],wy[1],wz[0],wz[1]);
	int I= dx>0? i+1:i-1;
	int J= dy>0? j+1:j-1;
	int K= dz>0? k+1:k-1;
	int ix[2]={i,I},jy[2]={j,J},kz[2]={k,K};
	for(int a=0; a<2; a++)
		for(int b=0; b<2; b++)
			for(int c=0; c<2; c++){
				if(ix[a]>=0 && ix[a]<nx && jy[b]>=0 && jy[b]<ny && kz[c]>=0 && kz[c]<nz)
					m[ix[a]][jy[b]][kz[c]]+=weight*wx[a]*wy[b]*wz[c];
			}
	means[0]+=weight*x; means[1]+=weight*y; means[2]+=weight*z;
	vars[0]+=weight*x*x; vars[1]+=weight*y*y; vars[2]+=weight*z*z;
}
void grid3::CIC(double weight,double *value,double x,double y,double z){
	int i=(int)((x-x0)/Deltax+.5);//index cell we're in
	int j=(int)((y-y0)/Deltay+.5);
	int k=(int)((z-z0)/Deltaz+.5);
	double dx=x-(x0+((double)i)*Deltax), dy=y-(y0+((double)j)*Deltay),
	dz=z-(z0+((double)k)*Deltaz);//ofsets from cell centre
	double wx[2],wy[2],wz[2];
	wts(dx,dy,dz,wx[0],wx[1],wy[0],wy[1],wz[0],wz[1]);
	int I= dx>0? i+1:i-1;
	int J= dy>0? j+1:j-1;
	int K= dz>0? k+1:k-1;
	int ix[2]={i,I},jy[2]={j,J},kz[2]={k,K};
	for(int a=0; a<2; a++)
		for(int b=0; b<2; b++)
			for(int c=0; c<2; c++){
				if(ix[a]>=0 && ix[a]<nx && jy[b]>=0 && jy[b]<ny && kz[c]>=0 && kz[c]<nz){
					double fac=wx[a]*wy[b]*wz[c];
					m[ix[a]][jy[b]][kz[c]]+=weight*fac;
					for(int l=0; l<na; l++) av[ix[a]][jy[b]][kz[c]][l]+=weight*value[l]*fac;
				}

			}
	means[0]+=weight*x; means[1]+=weight*y; means[2]+=weight*z;
	vars[0]+=weight*x*x; vars[1]+=weight*y*y; vars[2]+=weight*z*z;
}
double grid3::dens(double x,double y,double z){//Unnormalised: need to divide by cell volume
	//The virtual cell centres should be at (i)Delta
	if(x<x0 || y<y0 || z<z0) return 0;// off grid3
	int i=(int)((x-x0)/Deltax+.5);//index cell we're in
	int j=(int)((y-y0)/Deltay+.5);
	int k=(int)((z-z0)/Deltaz+.5);
	if(i<0 || j<0 || k<0 || i>=nx || j>=ny || k>=nz) return 0; //off grid
	double dx=x-(x0+((double)i)*Deltax),dy=y-(y0+((double)j)*Deltay),
	dz=z-(z0+((double)k)*Deltaz);//ofsets from cell centre
	double wx[2],wy[2],wz[2];
	wts(dx,dy,dz,wx[0],wx[1],wy[0],wy[1],wz[0],wz[1]);
	int I= dx>0? i+1:i-1;
	int J= dy>0? j+1:j-1;
	int K= dz>0? k+1:k-1;
	int ix[2]={i,I},jy[2]={j,J},kz[2]={k,K};
	double mass=0;
	for(int a=0; a<2; a++)
		for(int b=0; b<2; b++)
			for(int c=0; c<2; c++){
				if(ix[a]>=0 && ix[a]<nx && jy[b]>=0 && jy[b]<ny && kz[c]>=0 && kz[c]<nz){
					double fac=wx[a]*wy[b]*wz[c];
					mass+=m[ix[a]][jy[b]][kz[c]]*fac;
				}

			}
	return mass;
}
double grid3::dens(double x,double y,double z,double *avs){//Unnormalised: need to divide by cell volume
	//The virtual cell centres should be at (i)Delta
	if(x<x0 || y<y0 || z<z0) return 0;// off grid3
	int i=(int)((x-x0)/Deltax+.5);//index cell we're in
	int j=(int)((y-y0)/Deltay+.5);
	int k=(int)((z-z0)/Deltaz+.5);
	if(i<0 || j<0 || k<0 || i>=nx || j>=ny || k>=nz) return 0; //off grid
	double dx=x-(x0+((double)i)*Deltax),dy=y-(y0+((double)j)*Deltay),
	dz=z-(z0+((double)k)*Deltaz);//ofsets from cell centre
	double wx[2],wy[2],wz[2];
	wts(dx,dy,dz,wx[0],wx[1],wy[0],wy[1],wz[0],wz[1]);
	int I= dx>0? i+1:i-1;
	int J= dy>0? j+1:j-1;
	int K= dz>0? k+1:k-1;
	int ix[2]={i,I},jy[2]={j,J},kz[2]={k,K};
	double mass=0; for(int l=0;l<na;l++) avs[l]=0;
	for(int a=0; a<2; a++)
		for(int b=0; b<2; b++)
			for(int c=0; c<2; c++){
				if(ix[a]>=0 && ix[a]<nx && jy[b]>=0 && jy[b]<ny && kz[c]>=0 && kz[c]<nz){
					double fac=wx[a]*wy[b]*wz[c];
					mass+=m[ix[a]][jy[b]][kz[c]]*fac;
					for(int l=0; l<na; l++) avs[l]+=av[ix[a]][jy[b]][kz[c]][l]*fac;
				}

			}
	return mass;
}
void grid3::get_values(int l,float blank,float ***mp) const{//returns av values with blank value of empty cells
	if(l<0){
		for(int i=0;i<nx;i++)
			for(int j=0;j<ny;j++)
				for(int k=0; k<nz; k++) mp[i][j][k]=m[i][j][k];
	} else
		for(int i=0;i<nx;i++){
			for(int j=0;j<ny;j++)
				for(int k=0;k<nz;k++){
					if(m[i][j][k]>2){
						mp[i][j][k]=av[i][j][k][l]/m[i][j][k];
					}else{
						mp[i][j][k]=blank;
					}
				}
		}
}

float grid3::get_total(void){
	float sum=0;
	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
			for(int k=0; k<nz; k++)
				sum+=m[i][j][k];
	return sum;
}
void grid3::get_stats(float* _means,float* disps){
	float tot=get_total();
	for(int i=0; i<3; i++) _means[i]=means[i]/tot;
	for(int i=0; i<3; i++){
		float sq=vars[i]/tot-pow(_means[i],2);
		disps[i]= sq>0? sqrt(sq) : 0;
	}
}
/*
void grid3::dump(FILE *ofile){
	fprintf(ofile,"%d %d %d %d %g %g %g %g %g %g\n",nx,ny,nz,na,x0,y0,z0,Deltax,Deltay,Deltaz);
	compress(ofile,m,nx,ny,nz);
	if(na>0) compress(ofile,av,nx,ny,nz,na);
}*/
void grid3::zero(void){
	for(int i=0; i<nx; i++){
		for(int j=0; j<ny; j++){
			for(int k=0; k<nz; k++){
				m[i][j][k]=0;
			}
		}
	}
}
void grid3::divide(grid3 &grd){
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++)
			for(int k=0; k<nz; k++){
				double M=grd.dens(x0+i*Deltax,y0+j*Deltay,z0+k*Deltaz);
				if(fabs(M)>2) m[i][j][k]/=M;
			}
}
void grid3::divide(float fac){
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++)
			for(int k=0; k<nz; k++){
				m[i][j][k]/=fac;
			}
	for(int i=0; i<3; i++){
		means[i]/=fac; vars[i]/=fac;
	}
}
void grid3::subtract(grid3 &grd){
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++)
			for(int k=0; k<nz; k++){
				m[i][j][k]-=grd.dens(x0+i*Deltax,y0+j*Deltay,z0+k*Deltaz);
			}
}
float grid3::reveal(int i,int j,int k,float* avij){
	for(int l=0; l<na; l++) avij[l]=av[i][j][k][l];
	return m[i][j][k];
}
void grid3::add(grid3& grd){
	int ix,iy,iz,ia;
	grd.get_ints(ix,iy,iz,ia);
	if(ix!=nx || iy!=ny || iz!=nz || ia!=na){
		printf("Error: illegal attempt to add incommensurate grid3s\n"); return;
	}
	float* avij=NULL;
	if(na>0) avij=new float[na];
	for(int i=0; i<nx; i++){
		for(int j=0; j<ny; j++)
			for(int k=0; k<nz; k++){
				m[i][j][k]+=grd.reveal(i,j,k,avij);
				for(int l=0; l<na; l++) av[i][j][k][l]+=avij[l];
			}
	}
	if(na>0) delete[] avij;
}
void grid3::toggle_grey(void){
	greyit=(greyit+1)%2;
}

//Now routines involving mgo::plt
/*
void grid3::plot_slice(char c,double h,mgo::plt* pl,bool show,bool logit){
	if(!(c=='x' || c=='y' || c=='z')){
		printf("First arg plot_slice must be x,y,z\n");
		return;
	}
	double x1,y1,z1,Dx1,Dy1,Dz1;
	int nx1,ny1,nz1;
	if(c=='x'){
		nx1=ny; ny1=nz; nz1=nx; x1=y0; y1=z0; z1=x0; Dx1=Deltay; Dy1=Deltaz; Dz1=Deltax;
	}else if(c=='y'){
		nx1=nx; ny1=nz; nz1=ny; x1=x0; y1=z0; z1=y0; Dx1=Deltax; Dy1=Deltaz; Dz1=Deltay;
	}else if(c=='z'){
		nx1=nx; ny1=ny; nz1=nz; x1=x0; y1=y0; z1=z0; Dx1=Deltax; Dy1=Deltay; Dz1=Deltaz;
	}
	grid2 grd2(2*nx1-1,2*ny1-1,0,x1,x1+(nx1-1)*Dx1,y1,y1+(ny1-1)*Dy1);
	float lmin=1e10,lmax=-1e10;
	float z=z1+h*Dz1;
	for(int i=0; i<2*nx1-1; i++){
		double x=x1 + i*Dx1*.5;
		for(int j=0; j<2*ny1-1; j++){
			double y=y1 + j*Dy1*.5;
			if(c=='z')
				grd2.CIC(dens(x,y,z),x,y);
			else if(c=='y')
				grd2.CIC(dens(x,z,y),x,y);
			else
				grd2.CIC(dens(z,x,y),x,y);
		}
	}
	grd2.plot_dens(pl,show,logit);
}
std::pair<float,float> grid3::plot_projection(char c,mgo::plt* pl,bool show,bool logit){
	if(!(c=='x' || c=='y' || c=='z')){
		printf("First arg of plot_slice must be one of x,y,z\n");
		return std::make_pair((float)0.,(float)0.);
	}
	double x1,y1,z1,Dx1,Dy1,Dz1;
	int nx1,ny1,nz1;
	if(c=='x'){      //project along x
		nx1=ny; ny1=nz; nz1=nx; x1=y0; y1=z0; z1=x0; Dx1=Deltay; Dy1=Deltaz; Dz1=Deltax;
	}else if(c=='y'){//project along y
		nx1=nx; ny1=nz; nz1=ny; x1=x0; y1=z0; z1=y0; Dx1=Deltax; Dy1=Deltaz; Dz1=Deltay;
	}else if(c=='z'){//project along z
		nx1=nx; ny1=ny; nz1=nz; x1=x0; y1=y0; z1=z0; Dx1=Deltax; Dy1=Deltay; Dz1=Deltaz;
	}
	grid2 grd2(2*nx1-1,2*ny1-1,0,x1,x1+(nx1-1)*Dx1,y1,y1+(ny1-1)*Dy1);
	for(int i=0; i<2*nx1-1; i++){
		float x=x1+i*Dx1*.5;
		for(int j=0; j<2*ny1-1; j++){
			float y=y1+j*Dy1*.5;
			for(int k=0; k<nz1; k++){
				float z=z1+k*Dz1;
				if(c=='z')
					grd2.CIC(dens(x,y,z),x,y);
//					mp[i][j] += m[i][j][k];
				else if(c=='y')
					grd2.CIC(dens(x,z,y),x,y);
//					mp[i][j] += m[i][k][j];
				else
					grd2.CIC(dens(z,x,y),x,y);
//					mp[i][j] += m[k][i][j];
			}
		}
	}
	return grd2.plot_dens(pl,show,logit);
}
void grid3::plot_projection(float lmin,float lmax,char c,mgo::plt* pl,bool show,bool logit){
	if(!(c=='x' || c=='y' || c=='z')){
		printf("First arg of plot_slice must be one of x,y,z\n");
		return;
	}
	double x1,y1,z1,Dx1,Dy1,Dz1;
	int nx1,ny1,nz1;
	if(c=='x'){      //project along x
		nx1=ny; ny1=nz; nz1=nx; x1=y0; y1=z0; z1=x0; Dx1=Deltay; Dy1=Deltaz; Dz1=Deltax;
	}else if(c=='y'){//project along y
		nx1=nx; ny1=nz; nz1=ny; x1=x0; y1=z0; z1=y0; Dx1=Deltax; Dy1=Deltaz; Dz1=Deltay;
	}else if(c=='z'){//project along z
		nx1=nx; ny1=ny; nz1=nz; x1=x0; y1=y0; z1=z0; Dx1=Deltax; Dy1=Deltay; Dz1=Deltaz;
	}
	grid2 grd2(2*nx1-1,2*ny1-1,0,x1,x1+(nx1-1)*Dx1,y1,y1+(ny1-1)*Dy1);
	for(int i=0; i<2*nx1-1; i++){
		float x=x1+i*Dx1*.5;
		for(int j=0; j<2*ny1-1; j++){
			float y=y1+j*Dy1*.5;
			for(int k=0; k<nz1; k++){
				float z=z1+k*Dz1;
				if(c=='z')
					grd2.CIC(dens(x,y,z),x,y);
//					mp[i][j] += m[i][j][k];
				else if(c=='y')
					grd2.CIC(dens(x,z,y),x,y);
//					mp[i][j] += m[i][k][j];
				else
					grd2.CIC(dens(z,x,y),x,y);
//					mp[i][j] += m[k][i][j];
			}
		}
	}
	grd2.plot_dens(lmin,lmax,pl,show,logit);
}
void grid3::x_hist(mgo::plt* pl,bool logit){
	float nmax=-30,nmin=30;
	float* n=new float[nx];
	float* x=new float[nx];
	for(int i=0;i<nx;i++){//for each column
		n[i]=0; x[i]=x0+((double)i)*Deltax;
		for(int j=0; j<ny; j++)
			for(int k=0;k<nz;k++) n[i]+=m[i][j][k];//add cell height j
		nmin=MIN(nmin,n[i]); nmax=MAX(nmax,n[i]);
	}
	if(logit){
		for(int i=0;i<nx;i++) n[i]=log10(MAX(1e-20,n[i]));
		nmin=log10(MAX(1e-20,nmin)); nmax=log10(MAX(1e-20,nmax));
	}
	pl->setlimits(x0,0,x0+(nx)*Deltax,1.05*(nmax)); pl->histogram(x,n,nx);
	delete[] n; delete[] x;
}
void grid3::y_hist(mgo::plt* pl,bool logit){
	float nmax=-30,nmin=30;
	float* n=new float[ny];
	float* y=new float[ny];
	for(int i=0;i<ny;i++){//for each row
		n[i]=0; y[i]=y0+((double)i)*Deltay;
		for(int j=0; j<nx; j++)
			for(int k=0; k<nz; k++) n[i]+=m[j][i][k];//add cell position j
		nmin=MIN(nmin,n[i]); nmax=MAX(nmax,n[i]);
	}
	if(logit){
		for(int i=0;i<ny;i++) n[i]=log10(MAX(1e-20,n[i]));
		nmin=log10(MAX(1e-20,nmin)); nmax=log10(MAX(1e-20,nmax));
	}
	//pl->setlimits(y0,.5,y0+(ny)*Deltay,1.05*log10(nmax)); //
	pl->histogram(y,n,ny);
	delete[] n; delete[] y;
}
*/