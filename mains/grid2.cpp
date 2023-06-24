#include "grid2.h"
/* Grid centres are at (x0,y0), (x0+i*Deltax,y0+k*Deltay) etc so half
 * of edge cells lie outside plot. This ensures that contour plots are
 * correct because we are computing values on a regular grid that
 * touches frame*/

grid2::grid2(int _nx,int _ny,int _na,double xmin,double xmax,double ymin,double ymax):
    nx(_nx), ny(_ny), na(_na), x0(xmin), y0(ymin) {
	greyit=false;
	Deltax=(xmax-xmin)/(double)(nx-1); Deltay=(ymax-ymin)/(double)(ny-1);
	m=fmatrix(nx,ny); av = na>0? fmatrix(nx,ny,na) : NULL;
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++){
			m[i][j]=0;
			for(int k=0;k<na;k++) av[i][j][k]=0;
		}
}
grid2::grid2(grid2& ogrid,bool copy) {//copy a grid & its data
	greyit=false;
	ogrid.disgorge(nx,ny,na,x0,y0,Deltax,Deltay);
	m=fmatrix(nx,ny); av = na>0? fmatrix(nx,ny,na) : NULL;
	if(copy){
		double* avs = na>0? new double[na] : NULL;
		for(int i=0;i<nx;i++){
			double x=x0+i*Deltax;
			for(int j=0;j<ny;j++){
				double y=y0+j*Deltay;		
				m[i][j]=ogrid.dens(x,y,avs);
				for(int k=0;k<na;k++) av[i][j][k]=m[i][j]*avs[k];
			}
		}
		if(na>0) delete[] avs;
	}
}
/*
grid2::grid2(FILE *ifile,int na1) {
	greyit=false;
	if(7!=fscanf(ifile,"%d %d %d %lg %lg %lg %lg",&nx,&ny,&na,&x0,&y0,&Deltax,&Deltay)){
		printf("grid2::grid2: something wrong with data file\n"); exit(0);
	}
	int na2=MAX(na1,na);
	m=fmatrix(nx,ny); av = na2>0? fmatrix(nx,ny,na2) : NULL;
	get(ifile,m,nx,ny);
	if(na>0) get(ifile,av,nx,ny,na);
	fclose(ifile);
	na=MAX(na1,na);
}*/
grid2::~grid2(){
	delmatrix(m,nx); if(na>0) delmatrix(av,nx,ny);
}
void grid2::get_centres(float* xs,float* ys){
	xs[0]=x0;
	for(int i=1;i<nx;i++) xs[i]=xs[i-1]+Deltax;
	ys[0]=y0;
	for(int i=1;i<ny;i++) ys[i]=ys[i-1]+Deltay;
}
void grid2::get_minmax(float& lmin,float& lmax) const{
	lmin=1e10; lmax=-1e10;
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++){
			if(m[i][j]<lmin) lmin=m[i][j];
			if(m[i][j]>lmax) lmax=m[i][j];
		}
}
void grid2::get_minmax(int k,float& lmin,float& lmax) const{
	lmin=1e10; lmax=-1e10;
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++){
			if(m[i][j]<.0001) continue;
			float average=av[i][j][k]/m[i][j];
			if(average<lmin) lmin=average;
			if(average>lmax) lmax=average;
		}
}
void grid2::setup(int _nx,int _ny,int _na,double xmin,double xmax,
		 double ymin,double ymax) {
	nx=_nx; ny=_ny; na=_na; x0=xmin; y0=ymin;
	greyit=false;
	Deltax=(xmax-xmin)/(double)(nx-1); Deltay=(ymax-ymin)/(double)(ny-1);
	m=fmatrix(nx,ny); av = na>0? fmatrix(nx,ny,na) : NULL;
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++){
			m[i][j]=0;
			for(int k=0;k<na;k++) av[i][j][k]=0;
		}
}
void grid2::disgorge(int& _nx,int& _ny,int& _na,
	       double& _x0,double& _y0,double& _Deltax,double& _Deltay){
	_nx=nx; _ny=ny; _na=na; _x0=x0; _y0=y0; _Deltax=Deltax; _Deltay=Deltay;
}
void grid2::regrid(grid2 &old,int in,int out){//put data from slot in into slot out
	int ona=old.navs();
	double* avs = new double[ona];
	for(int i=0;i<nx;i++){
		double x=x0+i*Deltax;
		for(int j=0;j<ny;j++){
			double y=y0+j*Deltay;
			m[i][j]=old.dens(x,y,avs);
			if(0!=m[i][j]) av[i][j][out]=m[i][j]*avs[in];
		}
	}
	delete[] avs;
}
double grid2::w(double x,double Dx){//w is fraction to home cell
	double f=1-fabs(x)/Dx;
	return fmax(0,f);
}
void grid2::wts(double x,double y,double &wx,double &wxp,double &wy,double &wyp){
	wx= fmax(0,1-fabs(x)/Deltax); wxp=1-wx;
	wy= fmax(0,1-fabs(y)/Deltay); wyp=1-wy;
}
void grid2::CIC(double weight,double x,double y){
	int i=(int)((x-x0)/Deltax+.5);//index cell we're in
	int j=(int)((y-y0)/Deltay+.5);//
	double dx=x-(x0+((double)i)*Deltax),dy=y-(y0+((double)j)*Deltay);//ofsets from cell centre
	double wx[2],wy[2];
	wts(dx,dy,wx[0],wx[1],wy[0],wy[1]);
	int I= dx>0? i+1:i-1;
	int J= dy>0? j+1:j-1;
	int ix[2]={i,I},jy[2]={j,J};
	for(int a=0; a<2; a++)
		for(int b=0; b<2; b++){
			if(ix[a]>=0 && ix[a]<nx && jy[b]>=0 && jy[b]<ny)
				m[ix[a]][jy[b]]+=weight*wx[a]*wy[b];
		}
}
void grid2::CIC(double weight,double* value,double x,double y){
	int i=(int)((x-x0)/Deltax+.5);//index cell we're in
	int j=(int)((y-y0)/Deltay+.5);//
	double dx=x-(x0+((double)i)*Deltax),dy=y-(y0+((double)j)*Deltay);//ofsets from cell centre
	double wx[2],wy[2];
	wts(dx,dy,wx[0],wx[1],wy[0],wy[1]);
	int I= dx>0? i+1:i-1;
	int J= dy>0? j+1:j-1;
	int ix[2]={i,I},jy[2]={j,J};
	for(int a=0; a<2; a++)
		for(int b=0; b<2; b++){
			if(ix[a]>=0 && ix[a]<nx && jy[b]>=0 && jy[b]<ny){
				double fac=wx[a]*wy[b];
				m[ix[a]][jy[b]]+=weight*fac;
				for(int k=0;k<na;k++) av[ix[a]][jy[b]][k]+=weight*value[k]*fac;
			}
		}
}
double grid2::test_CIC(double weight,double *value,double x,double y){
	double total=0; int nc=0;
	int i=(int)((x-x0)/Deltax+.5);//index cell we're in
	int j=(int)((y-y0)/Deltay+.5);//
	double dx=x-(x0+((double)i)*Deltax),dy=y-(y0+((double)j)*Deltay);//ofsets from cell centre
	double wx,wy,wxp,wyp; wts(dx,dy,wx,wxp,wy,wyp);
	printf("i, j: %d %d dx, dy: %f %f %f %f ",i,j,dx/Deltax,dy/Deltay,wx,wy);
	bool iok = (i>=0 && i<nx), jok = (j>=0 && j<ny);
	bool ipok = (i+1>=0 && i+1<nx), jpok = (j+1>=0 && j+1<ny);
	bool imok = (i-1>=0 && i-1<nx), jmok = (j-1>=0 && j-1<ny);
	if(iok && jok){
		m[i][j]+=weight*wx*wy; total+=wx*wy; nc++;
		for(int k=0;k<na;k++) av[i][j][k]+=weight*value[k]*wx*wy;
	}
	if(dy>0 && iok && jpok){
		m[i][j+1]+=weight*wx*wyp; total+=wx*wyp; nc++;
		for(int k=0;k<na;k++) av[i][j+1][k]+=weight*value[k]*wx*wyp;
	} else if(dy<0 && iok && jmok){
		m[i][j-1]+=weight*wx*wyp; total+=wx*wyp; nc++;
		for(int k=0;k<na;k++) av[i][j-1][k]+=weight*value[k]*wx*wyp;
	}		
	if(dx>0 && ipok){
		if(jok){
			m[i+1][j]+=weight*wxp*wy; total+=wxp*wy; nc++;
			for(int k=0;k<na;k++) av[i+1][j][k]+=weight*value[k]*wxp*wy;
		}
		if(dy>0 && jpok){
			m[i+1][j+1]+=weight*wxp*wyp; total+=wxp*wyp; nc++;
			for(int k=0;k<na;k++) av[i+1][j+1][k]+=weight*value[k]*wxp*wyp;
		} else if(dy<0 && jmok){
			m[i+1][j-1]+=weight*wxp*wyp; total+=wxp*wyp; nc++;
			for(int k=0;k<na;k++) av[i+1][j-1][k]+=weight*value[k]*wxp*wyp;
		}
	} else if(dx<0 && imok){
		if(jok){
			m[i-1][j]+=weight*wxp*wy; total+=wxp*wy; nc++;
			for(int k=0;k<na;k++) av[i-1][j][k]+=weight*value[k]*wxp*wy;
		}
		if(dy>0 && jpok){
			m[i-1][j+1]+=weight*wxp*wyp; total+=wxp*wyp; nc++;
			for(int k=0;k<na;k++) av[i-1][j+1][k]+=weight*value[k]*wxp*wyp;
		} else if(dy<0 && jmok){
			m[i-1][j-1]+=weight*wxp*wyp; total+=wxp*wyp; nc++;
			for(int k=0;k<na;k++) av[i-1][j-1][k]+=weight*value[k]*wxp*wyp;
		}
	}
	printf("nc: %d ",nc);
	return total;
}
double grid2::dens(double x,double y){//Unnormalised: divide by cell volume
	if(x<x0 || y<y0) return 0;// off grid
	int i=(int)((x-x0)/Deltax+.5);//index cell we're in
	int j=(int)((y-y0)/Deltay+.5);//
	if(i<0 || j<0 || i>=nx || j>=ny) return 0;//off grid
	double dx=x-(x0+((double)i)*Deltax),dy=y-(y0+((double)j)*Deltay);//ofsets from cell centre
	double wx[2],wy[2];
	wts(dx,dy,wx[0],wx[1],wy[0],wy[1]);
	int I= dx>0? i+1:i-1;
	int J= dy>0? j+1:j-1;
	int ix[2]={i,I},jy[2]={j,J};
	double mass=0;
	for(int a=0; a<2; a++)
		for(int b=0; b<2; b++){
			if(ix[a]>=0 && ix[a]<nx && jy[b]>=0 && jy[b]<ny)
				mass+=m[ix[a]][jy[b]]*wx[a]*wy[b];
		}
	return mass;
}
double grid2::dens(double x,double y,double *avs){//Unnormalised: need to divide by cell volume
	if(x<x0 || y<y0) return 0;// off grid
	int i=(int)((x-x0)/Deltax+.5);//index cell we're in
	int j=(int)((y-y0)/Deltay+.5);//
	if(i<0 || j<0 || i>=nx || j>=ny) return 0;//off grid
	double dx=x-(x0+((double)i)*Deltax),dy=y-(y0+((double)j)*Deltay);//ofsets from cell centre
	double wx[2],wy[2];
	wts(dx,dy,wx[0],wx[1],wy[0],wy[1]);
	int I= dx>0? i+1:i-1;
	int J= dy>0? j+1:j-1;
	int ix[2]={i,I},jy[2]={j,J};
	double mass=0; for(int k=0;k<na;k++) avs[k]=0;
	for(int a=0; a<2; a++)
		for(int b=0; b<2; b++){
			if(ix[a]>=0 && ix[a]<nx && jy[b]>=0 && jy[b]<ny){
				double fac=wx[a]*wy[b];
				mass+=m[ix[a]][jy[b]]*fac;
				for(int k=0;k<na;k++) avs[k]+=av[i][j][k]*fac;
			}
		}
	return mass;
}
void grid2::find_peak(float *x,float *n){
	//find peak & surrounding 68% confidence interval. n should
	//contain probabilities not logs
	float peak=n[0],s; s=n[0];
	float* N=new float[nx];
	int ip=0;
	for(int i=1;i<nx;i++){
		s+=n[i];
		if(n[i]>peak){
			ip=i; peak=n[i];
		}
	}
	for(int i=0;i<nx;i++) N[i]=n[i]/s;
	int i=ip,plus=ip+1,minus=ip-1; float prob=N[ip];
	while(plus<nx && minus>=0 && prob<0.68){
		if(N[plus]>N[minus]){
			prob+=N[plus]; plus++;
		}else{
			prob+=N[minus]; minus--;
		}
	}
	if(prob<0.68 && plus<nx){
		while(prob<0.68 && plus<nx){
			prob+=N[plus]; plus++;
		}
	}
	if(prob<0.68 && minus>=0){
		while(prob<0.68){
			prob+=N[minus]; minus--;
		}
	}
	plus--;	minus++;
	printf("peak/confidence: %f + %f - %f\n",x[ip],x[plus]-x[ip],x[ip]-x[minus]);
	delete[] N;
}
void grid2::print_dens(void){
	for(int i=0;i<ny;i++){
		for(int j=0;j<nx;j++) printf("%6.2f ",m[j][i]);
		printf("\n");
	}
}
void grid2::fillup(void){//fill in gaps in values
	float ***avp; avp=fmatrix(nx,ny,na);
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			if(m[i][j]>1)
				for(int k=0;k<na;k++) avp[i][j][k]=av[i][j][k]/m[i][j];
			else
				for(int k=0;k<na;k++) avp[i][j][k]=0;
		}
	}
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			if(m[i][j]==0){
				int ns=0;
				for(int ix=i-1;ix<=i+1;ix+=2){
					if(ix>=0 && ix<nx && m[ix][j]>1){
						for(int k=0;k<na;k++) avp[i][j][k]+=avp[ix][j][k];
						ns++;
					}
				}
				for(int iy=j-1;iy<=j+1;iy+=2){
					if(iy>=0 && iy<ny && m[i][iy]>1){
						for(int k=0;k<na;k++) avp[i][j][k]+=avp[i][iy][k];
						ns++;
					}
				}
				if(ns>0)  for(int k=0;k<na;k++) avp[i][j][k]/=ns;
			}
		}
	}
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			if(m[i][j]==0 && avp[i][j][0]!=0){
				for(int k=0;k<na;k++) av[i][j][k]=avp[i][j][k];
				m[i][j]=1;
			}
		}
	}
	delmatrix(avp,nx,ny);
}	
void grid2::get_values(int k,float blank,float **mp) const{//returns av values with blank value of empty cells
	if(k<0){
		for(int i=0;i<nx;i++)
			for(int j=0;j<ny;j++) mp[i][j]=m[i][j];
	} else
		for(int i=0;i<nx;i++){
			for(int j=0;j<ny;j++){
				if(m[i][j]>2){
					mp[i][j]=av[i][j][k]/m[i][j];
				}else{
					mp[i][j]=blank;
				}
			}
		}
}
float grid2::get_total(void){
	float sum=0;
	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
			sum+=m[i][j];
	return sum;
}
void grid2::zero(void){
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			m[i][j]=0;
//			for(int k=0;k<na;k++) av[i][j][k]=0;
		}
	}
}
void grid2::divide(grid2 &grd){
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++){
			double M=grd.dens(x0+i*Deltax,y0+j*Deltay);
			if(fabs(M)>0) m[i][j]/=M;
		}
}
void grid2::subtract(grid2 &grd){
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++)
			m[i][j]-=grd.dens(x0+i*Deltax,y0+j*Deltay);
}
float grid2::reveal(int i,int j,float* avij){
	for(int k=0; k<na; k++) avij[k]=av[i][j][k];
	return m[i][j];
}
void grid2::add(grid2& grd){
	int ix,iy,iz;
	grd.get_ints(ix,iy,iz);
	if(ix!=nx || iy!=ny || iz!=na){
		printf("Error: illegal attempt to add incommensurate grids\n"); return;
	}
	float* avij=NULL;
	if(na>0) avij=new float[na];
	for(int i=0; i<nx; i++){
		for(int j=0; j<ny; j++){
			m[i][j]+=grd.reveal(i,j,avij);
			for(int k=0; k<na; k++) av[i][j][k]+=avij[k];
		}
	}
	if(na>0) delete[] avij;
}
void grid2::filter(void){
	float **mp; mp=fmatrix(nx,ny);
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++) mp[i][j]=m[i][j];
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			if(m[i][j]>1e-30){
				double sum=0,val=m[i][j]; int ns=0,nh=0;//compute mean of neighbours
				bool change=true;
				for(int k=-1;k<2;k++){
					int ik=i+k;
					if(ik>=0 && ik<nx){//keep within grid
						for(int l=-1;l<2;l++){
							if(k!=0 || l!=0){
								int jl=j+l;
								if(jl>=0 && jl<ny){
									if(m[ik][jl]>1e-30){
										if(fabs(m[ik][jl]-val)>fabs(m[ik][jl])){
											nh++;; sum+=m[ik][jl]; ns++;
										}
									}
								}
							}
						}
					}
				}
				if(sum!=0){
					sum/=(double)ns;
					if(nh>4){
						//printf("(%d %d) %g %g\n",i,j,m[i][j],sum);
						mp[i][j]=sum;//replace with mean
					}
				}
			}
		}
	}
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++) m[i][j]=mp[i][j];
	delmatrix(mp,nx);
}
void grid2::toggle_grey(void){
	greyit=(greyit+1)%2;
}

//routines involving mgo::plt
/*
void grid2::plot_contour(float frac,mgo::plt* pl){
	//contour inside which lies fraction frac of probability
	float h[1],**mp; h[0]=log10(get_contour(frac));
	mp=fmatrix(nx,ny);
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++) mp[i][j]=log10(m[i][j]);
	pl->setcolour("white"); pl->setlweight(1);
	pl->contour(mp,nx,ny,h,1);
	pl->setcolour("black"); pl->setlweight(0);
	delmatrix(mp,nx);
}
void grid2::contours(float* h,int nc,mgo::plt* pl) const{
	pl->contour(m,nx,ny,h,nc);
}
void grid2::x_hist(mgo::plt* pl,bool logit){
	float nmax=-30,nmin=30;
	float* n=new float[nx];
	float* x=new float[nx];
	for(int i=0;i<nx;i++){//for each column
		n[i]=0; x[i]=x0+((double)i)*Deltax;
		for(int j=0; j<ny; j++) n[i]+=m[i][j];//add cell height j
		nmin=MIN(nmin,n[i]); nmax=MAX(nmax,n[i]);
	}
	find_peak(x,n);
	if(logit){
		for(int i=0;i<nx;i++) n[i]=log10(MAX(1e-20,n[i]));
		nmin=log10(MAX(1e-20,nmin)); nmax=log10(MAX(1e-20,nmax));
	}
	pl->setlimits(x0,0,x0+(nx)*Deltax,1.05*(nmax)); pl->histogram(x,n,nx);
	delete[] n; delete[] x;
}
void grid2::x_hist(FILE *ofile,mgo::plt* pl,bool logit){
	float nmax=-30,nmin=30;
	float* n=new float[nx];
	float* x=new float[nx];
	for(int i=0;i<nx;i++){//for each column
		n[i]=0; x[i]=x0+((double)i)*Deltax;
		for(int j=0; j<ny; j++) n[i]+=m[i][j];//add cell height j
		nmin=MIN(nmin,n[i]); nmax=MAX(nmax,n[i]);
	}
	find_peak(x,n);
	if(logit){
		for(int i=0;i<nx;i++) n[i]=log10(MAX(1e-20,n[i]));
		nmin=log10(MAX(1e-20,nmin)); nmax=log10(MAX(1e-20,nmax));
	}
	pl->setlimits(x0,0,x0+(nx)*Deltax,1.05*(nmax)); pl->histogram(x,n,nx);
	compress(ofile,x,nx); compress(ofile,n,nx);
	delete[] n; delete[] x;
}
void grid2::y_hist(mgo::plt* pl,bool logit){
	float nmax=-30,nmin=30;
	float* n=new float[ny];
	float* y=new float[ny];
	for(int i=0;i<ny;i++){//for each row
		n[i]=0; y[i]=y0+((double)i)*Deltay;
		for(int j=0; j<nx; j++) n[i]+=m[j][i];//add cell position j
		nmin=MIN(nmin,n[i]); nmax=MAX(nmax,n[i]);
	}
	find_peak(y,n);
	if(logit){
		for(int i=0;i<ny;i++) n[i]=log10(MAX(1e-20,n[i]));
		nmin=log10(MAX(1e-20,nmin)); nmax=log10(MAX(1e-20,nmax));
	}
	//pl->setlimits(y0,.5,y0+(ny)*Deltay,1.05*log10(nmax)); //
	pl->histogram(y,n,ny);
	delete[] n; delete[] y;
}
void grid2::y_hist(FILE *ofile,mgo::plt* pl,bool logit){
	float nmax=-30,nmin=30;
	float* n=new float[ny];
	float* y=new float[ny];
	for(int i=0;i<ny;i++){//for each column
		n[i]=0; y[i]=y0+((double)i)*Deltay;
		for(int j=0; j<nx; j++) n[i]+=m[j][i];//add cell position j
		nmin=MIN(nmin,n[i]); nmax=MAX(nmax,n[i]);
	}
	find_peak(y,n);
	if(logit){
		for(int i=0;i<ny;i++) n[i]=log10(MAX(1e-20,n[i]));
		nmin=log10(MAX(1e-20,nmin)); nmax=log10(MAX(1e-20,nmax));
	}
	pl->setlimits(y0,.5,y0+(ny)*Deltay,1.05*log10(nmax)); pl->histogram(y,n,ny);
	compress(ofile,y,ny); compress(ofile,n,ny);
	delete[] n; delete[] y;
}
std::pair<float,float> grid2::plot_dens(mgo::plt* pl,bool show,bool logit){
	float lmin=1e10,lmax=-1e10,**mp;
	mp=fmatrix(nx,ny);
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++){
			mp[i][j]= logit? log10(MAX(1e-30,m[i][j])) : m[i][j];
			if((logit && mp[i][j]>-2) || (!logit && mp[i][j]> 0)){
				lmin=MIN(lmin,mp[i][j]); lmax=MAX(lmax,mp[i][j]);
			}
		}
	//lmin=MAX(lmin,lmax-2.5);
//	printf("lmin/max: %g %g\n",lmin,lmax);
	if(greyit)
//		pl->greywhite(mp,nx,ny,lmax,lmin+.1*(lmax-lmin),lmin+.1*(lmax-lmin),1);
		pl->greywhite(mp,nx,ny,lmin,lmax,lmin,show);
	else
		pl->colourwh(mp,nx,ny,lmin,lmax,lmin,show);
	//pl->colourwh(mp,nx,ny,lmin+.1*(lmax-lmin),lmax,lmin+.2*(lmax-lmin));
	pl->box(0,0,0,0); delmatrix(mp,nx);
	return std::make_pair(lmin,lmax);
}
void grid2::plot_dens(grid2 &grd,mgo::plt* pl){//plot log of density divided by that of grd
	float lmin=1e10,lmax=-1e10,**mp;
	mp=fmatrix(nx,ny);
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++){
			float bot=grd.dens(x0+i*Deltax,y0+j*Deltay);
			mp[i][j]=log10(MAX(1e-30,m[i][j]/bot));
			if(m[i][j]>5){
				lmin=MIN(lmin,mp[i][j]); lmax=MAX(lmax,mp[i][j]);
			}
		}
	lmin=MAX(lmin,lmax-2.5);
	//printf("lmin/max: %g %g\n",lmin,lmax);
	if(greyit)
		pl->greywhite(mp,nx,ny,lmax,lmin+.1*(lmax-lmin),lmin+.1*(lmax-lmin),false);
	else
		pl->colourwh(mp,nx,ny,lmin+.1*(lmax-lmin),lmax,lmin+.1*(lmax-lmin));
	pl->box(0,0,0,0); delmatrix(mp,nx);
}
void grid2::plot_dens(float lmin,float lmax,mgo::plt* pl,bool show,bool logit){
	float **mp;
	mp=fmatrix(nx,ny);
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			mp[i][j]= logit? log10(MAX(1e-30,m[i][j])) : m[i][j];
		}
	}
	if(greyit)
		pl->greywhite(mp,nx,ny,lmin,lmax,lmin,show);
	else
		pl->colourwh(mp,nx,ny,lmin,lmax,lmin,show);
	pl->box(0,0,0,0); delmatrix(mp,nx);
}
void grid2::plot_densR(mgo::plt* pl,bool show,bool logit){
	float lmin=1e10,lmax=-1e10,**mp;
	mp=fmatrix(nx,ny);
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++){
			mp[i][j]= logit? log10(MAX(1e-30,m[i][j])) : m[i][j];
			if((logit|| mp[i][j]>-29) || (!logit && mp[i][j]>1)){
				lmin=MIN(lmin,mp[i][j]); lmax=MAX(lmax,mp[i][j]);
			}
		}
	//printf("lmin, max: %f %f\n",lmin,lmax);
	pl->colourwhR(mp,nx,ny,lmin+.1*(lmax-lmin),lmax,lmin+.1*(lmax-lmin),show);
	pl->box(0,0,0,0); delmatrix(mp,nx);
}
void grid2::plot_values(int k,float blank,mgo::plt* pl,bool show,bool logit) const{//plots av values with blank value of empty cells
	float lmin=1e6,lmax=-1e6,**avp;
	avp=fmatrix(nx,ny);
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			if(m[i][j]>0){
				avp[i][j]=av[i][j][k]/m[i][j];
				if(logit) avp[i][j]=log10(avp[i][j]);
				if(m[i][j]>.0001){
					lmin=MIN(lmin,avp[i][j]);
					lmax=MAX(lmax,avp[i][j]);
				}
			}else{
				avp[i][j]=blank;
			}
		}
	}
	//printf("lmin/max: %g %g\n",lmin,lmax);
	if(greyit)
		pl->greywhite(avp,nx,ny,lmin,lmax,lmin,show);
	else
		pl->colourwh(avp,nx,ny,lmin,lmax,lmin,show);
	pl->box(0,0,0,0); delmatrix(avp,nx);
}
void grid2::plot_valuesR(int k,float blank,mgo::plt* pl,bool show,bool logit) const{//plots av values with blank value of empty cells
	float lmin=1e6,lmax=-1e6,**avp;
	avp=fmatrix(nx,ny);
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			if(m[i][j]>0){
				avp[i][j]=av[i][j][k]/m[i][j];
				if(logit) avp[i][j]=log10(avp[i][j]);
				if(m[i][j]>.0001){
					lmin=MIN(lmin,avp[i][j]);
					lmax=MAX(lmax,avp[i][j]);
				}
			}else{
				avp[i][j]=blank;
			}
		}
	}
	//printf("lmin/max: %g %g\n",lmin,lmax);
	if(greyit)
		pl->greywhite(avp,nx,ny,lmin,lmax,lmin,show);
	else
		pl->colourwhR(avp,nx,ny,lmin,lmax,lmin,show);
	pl->box(0,0,0,0); delmatrix(avp,nx);
}
void grid2::plot_values(float lmin,float lmax,int k,float blank,mgo::plt* pl,bool show,bool logit) const{//plots av values with blank value of empty cells
	float **avp;
	avp=fmatrix(nx,ny);
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			if(m[i][j]>1){
				avp[i][j]= logit? log10(avp[i][j]/m[i][j]) : av[i][j][k]/m[i][j];
			}else{
				avp[i][j]=blank;
			}
		}
	}
	if(greyit)
		pl->greywhite(avp,nx,ny,lmin,lmax,lmin,show);
	else
		pl->colourwh(avp,nx,ny,lmin,lmax,lmin,show);
	pl->box(0,0,0,0); delmatrix(avp,nx);
}
void grid2::plot_valuesR(float lmin,float lmax,int k,float blank,mgo::plt* pl,bool show,bool logit) const{//plots av values with blank value of empty cells
	float **avp;
	avp=fmatrix(nx,ny);
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			if(m[i][j]>0){
				avp[i][j]= logit? log10(avp[i][j]/m[i][j]) : av[i][j][k]/m[i][j];
			}else{
				avp[i][j]=blank;
			}
		}
	}
	pl->colourwhR(avp,nx,ny,lmin,lmax,lmin,show);
	pl->box(0,0,0,0); delmatrix(avp,nx);
}
void grid2::plot_fn(float (*fn)(float,float),mgo::plt* pl,bool show){//plot provided function
	float lmin=1e6,lmax=-1e6,**avp;
	avp=fmatrix(nx,ny);
	for(int i=0;i<nx;i++){
		float x=x0+i*Deltax;
		for(int j=0;j<ny;j++){
			float y=y0+j*Deltay;
			avp[i][j]=fn(x,y);
			lmin=MIN(lmin,avp[i][j]);
			lmax=MAX(lmax,avp[i][j]);
		}
	}
	//printf("lmin/max: %g %g\n",lmin,lmax);
	if(greyit)
		pl->greywhite(avp,nx,ny,lmin,lmax,lmin,show);
	else
		pl->colourwh(avp,nx,ny,lmin,lmax,lmin,show);
	pl->box(0,0,0,0); delmatrix(avp,nx);
}	
void grid2::contour_fn(float (*fn)(float,float),float* h,int nc,mgo::plt* pl,bool show){//plot provided function
	float **avp; avp=fmatrix(nx,ny);
	for(int i=0;i<nx;i++){
		float x=x0+i*Deltax;
		for(int j=0;j<ny;j++){
			float y=y0+j*Deltay;
			//if(i==0) printf("%f ",y);
			avp[i][j]=fn(x,y);
		}
	}
	pl->contour(avp,nx,ny,h,nc);
	delmatrix(avp,nx);
}	
void grid2::plot_fn_values(int k,float blank,float (*fn)(float),mgo::plt* pl,bool show){//plots av values with blank value of empty cells
	float lmin=1e6,lmax=-1e6,**avp;
	avp=fmatrix(nx,ny);
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			if(m[i][j]>2){
				avp[i][j]=(*fn)(av[i][j][k]/m[i][j]);
				if(m[i][j]>.0001){
					lmin=MIN(lmin,avp[i][j]);
					lmax=MAX(lmax,avp[i][j]);
				}
			}else{
				avp[i][j]=blank;
			}
		}
	}
	//printf("lmin/max: %g %g\n",lmin,lmax);
	if(greyit)
		pl->greywhite(avp,nx,ny,lmin,lmax,lmin,show);
	else
		pl->colourwh(avp,nx,ny,lmin,lmax,lmin,show);
	pl->box(0,0,0,0); delmatrix(avp,nx);
}
void grid2::percentiles(int n,float *p,mgo::plt* pl){
	unsigned long np=nx*ny;
	unsigned long* indx=new unsigned long[nx*ny];
	double* arr=new double[nx*ny];
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++) arr[ny*i+j]=m[i][j];
	indexx(np,arr,indx);
	float sum=0,**mp; mp=fmatrix(nx,ny);
	for(int k=0;k<np;k++){
		int i=indx[k]/ny,j=indx[k]%ny;
		sum+=m[i][j]; mp[i][j]=sum;
	}
	sum/=100;
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++) mp[i][j]/=sum;
	pl->colourwh(mp,nx,ny,2,100,2,false);
	if(n>0) pl->contour(mp,nx,ny,p,n);
	delmatrix(mp,nx); delete[] indx; delete[] arr;
}
void grid2::contour_values(int k,float* h,int nc,float blank,mgo::plt* pl){
	float **avp;
	avp=fmatrix(nx,ny);
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			if(m[i][j]>2){
				avp[i][j]=av[i][j][k]/m[i][j];
			}else{
				avp[i][j]=blank;
			}
		}
	}
	pl->contour(avp,nx,ny,h,nc);
}

void grid2::contour_values(int k,int nc,float blank,mgo::plt* pl){//plots av values with blank value of empty cells
	float lmin=1e6,lmax=-1e6,**avp;
	avp=fmatrix(nx,ny);
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			if(m[i][j]>2){
				avp[i][j]=av[i][j][k]/m[i][j];
				lmin=MIN(lmin,avp[i][j]);
				lmax=MAX(lmax,avp[i][j]);
			}else{
				//if(i==nx/2) printf("%d %g %g\n",j,m[i][j],av[i][j]);
				avp[i][j]=blank;
			}
		}
	}
	//printf("lmin/max: %g %g\n",lmin,lmax);
	float* h=new float[nc];
	for(int i=0;i<nc;i++) h[i]=lmin+i*(lmax-lmin)/(float)(nc-1);
	pl->contour(avp,nx,ny,h,nc);
	delete[] h;
}
void grid2::plot_fn_values(float lmin,float lmax,int k,float blank,float (*fn)(float),mgo::plt* pl,bool show){//plots av values with blank value of empty cells
	float **mp;
	mp=fmatrix(nx,ny);
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			if(m[i][j]>2){
				mp[i][j]=fn(av[i][j][k]/m[i][j]);
			}else{
				printf("(%d %d) %g %g\n",i,j,m[i][j],av[i][j][k]);
				mp[i][j]=blank;
			}
		}
	}
	if(greyit)
		pl->greywhite(mp,nx,ny,lmin,lmax,lmin,show);
	else
		pl->colourwh(mp,nx,ny,lmin,lmax,lmin,show);
	pl->box(0,0,0,0); delmatrix(mp,nx);
}
void grid2::plot_dispersion(int mn,int sd,float blank,mgo::plt* pl,bool show) const{
	float **mp,lmin=1e10,lmax=0;
	mp=fmatrix(nx,ny);
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			if(m[i][j]>1){
				float d=av[i][j][sd]/m[i][j]-pow(av[i][j][mn]/m[i][j],2);
				if(d>10)
					printf("%d %d %g %g %g\n",i,j,d,m[i][j],av[i][j][mn],av[i][j][sd]);
				mp[i][j]= d>0? sqrt(d) : 0;
				lmin=MIN(lmin,mp[i][j]);
				lmax=MAX(lmax,mp[i][j]);
			} else {
				mp[i][j]=blank;
			}
		}
	}
	if(greyit)
		pl->greywhite(mp,nx,ny,lmin,lmax,lmin,show);
	else
		pl->colourwh(mp,nx,ny,lmin,lmax,lmin,show);
	pl->box(0,0,0,0);
	delmatrix(mp,nx);
}
void grid2::plot(float frac,mgo::plt* pl){
	plot_dens(pl);	plot_contour(frac,pl);
}
void grid2::dump(FILE *ofile){
	fprintf(ofile,"%d %d %d %g %g %g %g\n",nx,ny,na,x0,y0,Deltax,Deltay);
	compress(ofile,m,nx,ny);
	if(na>0) compress(ofile,av,nx,ny,na);
}
void grid2::test_pl(float x,float y,char* msg,mgo::plt* pl){
	pl->relocate(x,y);
	pl->label(msg);
}
*/