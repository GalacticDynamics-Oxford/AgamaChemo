#include "grid.h"
#define SIGN(A,B) ((B)>0?(A):(-A))
int *imatrix(int n){
	int *m1 = new int[n];
	for(int i=0;i<n;i++) m1[i]=0;
	return m1;
}
int **imatrix(int n,int m){
	int **m1 = new int*[n];
	for(int i=0; i<n; i++) m1[i] = imatrix(m);
	return m1;
}
int ***imatrix(int n,int m,int l){
	int ***m1 = new int**[n];
	for(int i=0; i<n; i++) m1[i] = imatrix(m,l);
	return m1;
}
void delmatrix(int **m1,int n){
	for(int i=0;i<n;i++) delete[] m1[i];
	delete [] m1;
}
void delmatrix(int ***m1,int n,int m){
	for(int i=0;i<n;i++) delmatrix(m1[i],m);
	delete [] m1;
}
double *dmatrix(int n){
	double *m1 = new double[n];
	for(int i=0;i<n;i++) m1[i]=0;
	return m1;
}
double **dmatrix(int n,int m){
	double **m1 = new double*[n];
	for(int i=0; i<n; i++) m1[i] = dmatrix(m);
	return m1;
}
double ***dmatrix(int n,int m,int l){
	double ***m1 = new double**[n];
	for(int i=0; i<n; i++) m1[i] = dmatrix(m,l);
	return m1;
}
double ****dmatrix(int n,int m,int l,int k){
	double ****m1 = new double***[n];
	for(int i=0; i<n; i++) m1[i] = dmatrix(m,l,k);
	return m1;
}
void delmatrix(double **m1,int n){
	for(int i=0;i<n;i++) delete[] m1[i];
	delete [] m1;
}
void delmatrix(double ***m1,int n,int m){
	for(int i=0;i<n;i++) delmatrix(m1[i],m);
	delete [] m1;
}
void delmatrix(double ****m1,int n,int m,int l){
	for(int i=0;i<n;i++) delmatrix(m1[i],m,l);
	delete [] m1;
}
float *fmatrix(int n){
	float *m1 = new float[n];
	for(int i=0;i<n;i++) m1[i]=0;
	return m1;
}
float **fmatrix(int n,int m){
	float **m1 = new float*[n];
	for(int i=0; i<n; i++) m1[i] = fmatrix(m);
	return m1;
}
float ***fmatrix(int n,int m,int l){
	float ***m1 = new float**[n];
	for(int i=0; i<n; i++) m1[i] = fmatrix(m,l);
	return m1;
}
float ****fmatrix(int n,int m,int l,int k){
	float ****m1 = new float***[n];
	for(int i=0; i<n; i++) m1[i] = fmatrix(m,l,k);
	return m1;
}
float *****fmatrix(int o,int n,int m,int l,int k){
	float *****m1 = new float****[o];
	for(int i=0; i<o; i++) m1[i] = fmatrix(n,m,l,k);
	return m1;
}
void delmatrix(float **m1,int n){
	for(int i=0;i<n;i++) delete[] m1[i];
	delete [] m1;
}
void delmatrix(float ***m1,int n,int m){
	for(int i=0;i<n;i++) delmatrix(m1[i],m);
	delete [] m1;
}
void delmatrix(float ****m1,int n,int m,int l){
	for(int i=0;i<n;i++) delmatrix(m1[i],m,l);
	delete [] m1;
}
void delmatrix(float *****m1,int o,int n,int m,int l){
	for(int i=0;i<o;i++) delmatrix(m1[i],n,m,l);
	delete [] m1;
}
void compress(FILE *tmpf,float *x,int npt){

	int i,j,k,lnx,rem,i1 = 90,i2 = 8100,i3 = 729000;
	float at2t23 = 8388608., big=1.e28,small=1.e-28;
	char m[6];
	m[5]='\0';

	const double aln2 = log(2.);
	for(k=0 ; k<npt ; k++){
		if (x[k]==0.){
			m[4] = 80;
			m[3] = 0;
			m[2] = 0;
			m[1] = 0;
			m[0] = 0;
		}else{
			if (fabs(x[k])>big) x[k]= SIGN(big, x[k]);
			if (fabs(x[k])<small) x[k] = SIGN(small, x[k]);
			if (fabs(x[k])>=1.)
				lnx = (int)(log(fabs(x[k]))/aln2) + 129;
			else
				lnx = (int)(log(fabs(x[k]))/aln2) + 128;
			j = lnx/i1;
			m[3] = lnx - j*i1;
			lnx = lnx - 129;
			if (lnx>0)
				rem = (int)((fabs(x[k])/pow(2.,(double)(lnx)) - 1.)*at2t23);
			else
				rem = (int)((fabs(x[k])*pow(2.,(double)(fabs(lnx))) - 1.)*at2t23);
			m[4] = rem/i3;
			rem = rem - m[4]*i3;
			m[2] = rem/i2;
			rem = rem - m[2]*i2;
			m[1] = rem/i1;
			m[0] = rem - m[1]*i1;
			m[4] = m[4] + j*12;
			if (x[k]<0.) m[4] = m[4] + 40;
		}
		for (i=0;i<5;i++){
			m[i] = m[i] + 33;
			if (m[i]>=94) m[i] = m[i] + 1;
			if (m[i]>=123) m[i] = m[i] + 1;
		}
		fprintf(tmpf,"%s",m);
		if(15*((k+1)/15)==k+1) fprintf(tmpf,"\n");
	}
	if(15*((npt)/15)!=npt) fprintf(tmpf,"\n");
}
void compress(FILE *tmpf,double *x,int npt){

	int i,j,k,lnx,rem,i1 = 90,i2 = 8100,i3 = 729000;
	float at2t23 = 8388608., big=1.e28,small=1.e-28;
	char m[6];
	m[5]='\0';

	const double aln2 = log(2.);
	for(k=0 ; k<npt ; k++){
		if (x[k]==0.){
			m[4] = 80;
			m[3] = 0;
			m[2] = 0;
			m[1] = 0;
			m[0] = 0;
		}else{
			if (fabs(x[k])>big) x[k]= SIGN(big, x[k]);
			if (fabs(x[k])<small) x[k] = SIGN(small, x[k]);
			if (fabs(x[k])>=1.)
				lnx = (int)(log(fabs(x[k]))/aln2) + 129;
			else
				lnx = (int)(log(fabs(x[k]))/aln2) + 128;
			j = lnx/i1;
			m[3] = lnx - j*i1;
			lnx = lnx - 129;
			if (lnx>0)
				rem = (int)((fabs(x[k])/pow(2.,(double)(lnx)) - 1.)*at2t23);
			else
				rem = (int)((fabs(x[k])*pow(2.,(double)(fabs(lnx))) - 1.)*at2t23);
			m[4] = rem/i3;
			rem = rem - m[4]*i3;
			m[2] = rem/i2;
			rem = rem - m[2]*i2;
			m[1] = rem/i1;
			m[0] = rem - m[1]*i1;
			m[4] = m[4] + j*12;
			if (x[k]<0.) m[4] = m[4] + 40;
		}
		for (i=0;i<5;i++){
			m[i] = m[i] + 33;
			if (m[i]>=94) m[i] = m[i] + 1;
			if (m[i]>=123) m[i] = m[i] + 1;
		}
		fprintf(tmpf,"%s",m);
		if(15*((k+1)/15)==k+1) fprintf(tmpf,"\n");
	}
	if(15*((npt)/15)!=npt) fprintf(tmpf,"\n");
}

bool get(FILE *tmpf,float *xg,int npt)//encode numbers t base 90
{
	float expon,at2t23= 8388608.;
	long int i1=90;
	int nread=0;
	char lin[100],m[5];
	int nline=npt/15;
	if(15*nline!=npt) nline=nline+1;//# of lines with data, including short lines
	for(int line=0; line<nline; line++){
		fgets(lin,100,tmpf);
		while(strlen(lin)<5 && !feof(tmpf)) fgets(lin,100,tmpf);
		int np=strlen(lin)/5;//we'll read np numbers from lin
		for(int kp=0; kp<np; kp++){
			int k=line*15+kp;
			for(int i = 0; i<5; i++){
				m[i] = (unsigned int)lin[5*kp+i];
				if (m[i]>=124) m[i] = m[i] - 1;
				if (m[i]>=95) m[i] = m[i] - 1;
				m[i] = m[i] - 33;}
			if (m[4]==80)
				xg[k] = 0.;
			else{
				if (m[4]>=40){
					m[4] = m[4] - 40;
					xg[k] = -1.;}
				else xg[k] = 1.;
				int j = m[4]/12;
				m[4] = m[4] - j*12;
				expon=j*i1+m[3]-129;
				xg[k]=xg[k]*((((m[4]*i1+m[2])*i1+m[1])*i1
					      +m[0])/at2t23+ 1.);
				if (expon>0.)
					xg[k]=xg[k]*pow(2.,expon);
				else
					xg[k]=xg[k]/pow(2.,-expon);
			}
			nread++;
		}
	}
	return nread==npt? true : false;
}
bool get(FILE *tmpf,double *xg,int npt)//encode numbers t base 90
{
	double expon,at2t23= 8388608.;
	long int i1=90;
	int nread=0;
	char lin[100],m[5];
	int nline=npt/15;
	if(15*nline!=npt) nline=nline+1;//# of lines with data, including short lines
	for(int line=0; line<nline; line++){
		fgets(lin,100,tmpf);
		while(strlen(lin)<5 && !feof(tmpf)) fgets(lin,100,tmpf);
//		printf("str: %zd\n",strlen(lin));
		int np=strlen(lin)/5;//we'll read np numbers from lin
		for(int kp=0; kp<np; kp++){
			int k=line*15+kp;
			for(int i = 0; i<5; i++){
				m[i] = (unsigned int)lin[5*kp+i];
				if (m[i]>=124) m[i] = m[i] - 1;
				if (m[i]>=95) m[i] = m[i] - 1;
				m[i] = m[i] - 33;}
			if (m[4]==80)
				xg[k] = 0.;
			else{
				if (m[4]>=40){
					m[4] = m[4] - 40;
					xg[k] = -1.;}
				else xg[k] = 1.;
				int j = m[4]/12;
				m[4] = m[4] - j*12;
				expon=j*i1+m[3]-129;
				xg[k]=xg[k]*((((m[4]*i1+m[2])*i1+m[1])*i1
					      +m[0])/at2t23+ 1.);
				if (expon>0.)
					xg[k]=xg[k]*pow(2.,expon);
				else
					xg[k]=xg[k]/pow(2.,-expon);
			}
			nread++;
		}
	}
	return nread==npt? true : false;
}
bool get(FILE *tmpf,float **array,int nx,int ny){
	bool res=true;
	for(int i=0; i<nx; i++) res=res && get(tmpf,array[i],ny);
	return res;
}
bool get(FILE *tmpf,float ***array,int nx,int ny,int nz){
	bool res=true;
	for(int i=0; i<nx; i++) res=res && get(tmpf,array[i],ny,nz);
	return res;
}
bool get(FILE *tmpf,float ****array,int nx,int ny,int nz,int na){
	bool res=true;
	for(int i=0; i<nx; i++) res=res && get(tmpf,array[i],ny,nz,na);
	return res;
}
bool get(FILE *tmpf,float *****array,int nw,int nx,int ny,int nz,int na){
	bool res=true;
	for(int i=0; i<nw; i++) res=res && get(tmpf,array[i],nx,ny,nz,na);
	return res;
}
bool get(FILE *tmpf,double **array,int nx,int ny){
	bool res=true;
	for(int i=0; i<nx; i++) res=res && get(tmpf,array[i],ny);
	return res;
}
bool get(FILE *tmpf,double ***array,int nx,int ny,int nz){
	bool res=true;
	for(int i=0; i<nx; i++) res=res && get(tmpf,array[i],ny,nz);
	return res;
}
bool get(FILE *tmpf,double ****array,int nx,int ny,int nz,int na){
	bool res=true;
	for(int i=0; i<nx; i++) res=res && get(tmpf,array[i],ny,nz,na);
	return res;
}

void compress(FILE *tmpf,float **array,int nx,int ny){
	for(int i=0; i<nx;i++) compress(tmpf,array[i],ny);
}
void compress(FILE *tmpf,float ***array,int nx,int ny,int nz){
	for(int i=0; i<nx; i++) compress(tmpf,array[i],ny,nz);
}
void compress(FILE *tmpf,float ****array,int nx,int ny,int nz,int na){
	for(int i=0; i<nx; i++) compress(tmpf,array[i],ny,nz,na);
}
void compress(FILE *tmpf,float *****array,int nw,int nx,int ny,int nz,int na){
	for(int i=0; i<nw; i++) compress(tmpf,array[i],nx,ny,nz,na);
}
void compress(FILE *tmpf,double **array,int nx,int ny){
	for(int i=0; i<nx;i++) compress(tmpf,array[i],ny);
}
void compress(FILE *tmpf,double ***array,int nx,int ny,int nz){
	for(int i=0; i<nx; i++) compress(tmpf,array[i],ny,nz);
}
void compress(FILE *tmpf,double ****array,int nx,int ny,int nz,int na){
	for(int i=0; i<nx; i++) compress(tmpf,array[i],ny,nz,na);
}
void compress(FILE* tmpf,std::vector<float>& x){
	int n=(int)x.size();
	float* xp = new float[n];
	for(int i=0; i<n; i++) xp[i]=x[i];
	compress(tmpf,xp,n);
	delete[] xp;
}
void compress(FILE* tmpf,std::vector<double>& x){
	int n=(int)x.size();
	double* xp = new double[n];
	for(int i=0; i<n; i++) xp[i]=x[i];
	compress(tmpf,xp,n);
	delete[] xp;
}
