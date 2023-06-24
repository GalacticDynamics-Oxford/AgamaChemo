#ifndef __grid__
#define __grid__
#include "grid.h"
#endif
#ifndef __grid4__
#define __grid4__
class grid4{
	private:
		double w(double,double);
		void wts(double,double,double,double,double&,double&,double&,double&,double&,double&,double&,double&);
		float ****m, *****av;
		double w0,x0,y0,z0,Deltaw,Deltax,Deltay,Deltaz;
		int nw,nx,ny,nz,na;
		bool greyit;
	public:
		grid4(void){};
		grid4(int,int,int,int,int,double,double,double,double,double,double,double,double);
		grid4(grid4&,bool copy=true);//copy grid & its data
//		grid4(FILE*,int=0);
		~grid4(void);
		int navs(){
			return na;
		}
		double CellVolume(void){
			return Deltaw*Deltax*Deltay*Deltaz;
		}
		void get_centres(float*,float*,float*,float*);//returns cell centres
		void setup(int,int,int,int,int,double,double,double,double,double,double,double,double);
		void setup(grid4&);
		void CIC(double,double,double,double,double);//add contrib to density
		void CIC(double,double*,double,double,double,double);//add contrib to density and averages
		double dens(double,double,double,double);//recover density at (x,y)
		double dens(double,double,double,double,double*);//recover density and averages at (x,y)
		float get_total(void);//returns total number on grid
		void get_minmax(float&,float&) const;//return lmin, lmax
		void get_minmax(int,float&,float&) const;//return min, max of kth value
		void get_ints(int& iw,int& ix,int& iy,int& iz,int &ia){
			iw=nw; ix=nx; iy=ny; iz=nz; ia=na;
		}
		void get_values(int,float,float****) const;// recovers totals (int<0) or average of kth value; float if m=0
//		void dump(FILE*);//write m and av to file
		void zero(void);//initialise density and averages
		void disgorge(int&,int&,int&,int&,int&,double&,double&,double&,double&,double&,double&,double&,double&);//returns nx,ny,x0,y0,Deltas
		void regrid(grid4&,int,int);//put data from slot in into slot out 
		void divide(grid4&);
		void divide(float);
		void subtract(grid4&);//
		float reveal(int h,int i,int j,int k,float* avij);
		void add(grid4& grd);
		void toggle_grey(void);
//Now routines involving mgo::plt
/*		void plot_slice(char,double,char,double,mgo::plt*,bool show=true,bool logit=true);//slice is at char axis = double
		void plot_projection(char,char,mgo::plt*,bool show=true,bool logit=true);//summed over  char axes
		void plot_mean(char,char,double,double,mgo::plt*,bool show=true,bool logit=false);//summed over  char axes
		void plot_meanR(char,char,double,double,mgo::plt*,bool show=true,bool logit=false);//summed over  char axes
		void plot_plaque(double lmin,double lmax,double ymin,double ymax,double zmin,double zmax,mgo::plt*,bool show=true,bool logit=false);
*/
};

#endif
