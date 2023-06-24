#ifndef __grid__
#define __grid__
#include "grid.h"
#endif
#ifndef __grid3__
#define __grid3__
class grid3{
	private:
		double w(double,double);
		void wts(double,double,double,double&,double&,double&,double&,double&,double&);
		float ***m,****av;
		double means[3],vars[3];
		double x0,y0,z0,Deltax,Deltay,Deltaz;
		int nx,ny,nz,na;
		bool greyit;
	public:
		grid3(void){};
		grid3(int,int,int,int,double,double,double,double,double,double);
		grid3(grid3&,bool copy=true);//copy grid & its data
//		grid3(FILE*,int=0);
		~grid3(void);
		int navs(){
			return na;
		}
		double CellVolume(void){
			return Deltax*Deltay*Deltaz;
		}
		void get_centres(float*,float*,float*);//returns cell centres
		void setup(int,int,int,int,double,double,double,double,double,double);
		void CIC(double,double,double,double);//add contrib to density
		void CIC(double,double*,double,double,double);//add contrib to density and averages
		double dens(double,double,double);//recover density at (x,y)
		double dens(double,double,double,double*);//recover density and averages at (x,y)
		double add_dens(grid3&,double (*)(double,double),int,int,bool);
		double add_backgnd(double (*)(double,double,double),double,bool);
		float get_total(void);//returns total number on grid
		void get_stats(float* means,float* disps);//returns means and dispersions
		void get_minmax(float&,float&) const;//return lmin, lmax
		void get_minmax(int,float&,float&) const;//return min, max of kth value
		void get_ints(int& ix,int& iy,int& iz,int &ia){
			ix=nx; iy=ny; iz=nz; ia=na;
		}
		void get_values(int,float,float***) const;// recovers totals (int<0) or average of kth value; float if m=0
//		void dump(FILE*);//write m and av to file
		void zero(void);//initialise density and averages
		void disgorge(int&,int&,int&,int&,double&,double&,double&,double&,double&,double&);//returns nx,ny,x0,y0,Deltas
		void regrid(grid3&,int,int);//put data from slot in into slot out 
		void divide(grid3&);
		void divide(float);
		void subtract(grid3&);//
		float reveal(int i,int j,int k,float* avij);
		void add(grid3& grd);
		void toggle_grey(void);
//Now routines involving mgo::plt
/*		std::pair<float,float> plot_projection(char c,mgo::plt* pl,bool show=true,bool logit=true);//plot density projected along c axis
		void plot_slice(char c,double h,mgo::plt* pl,bool show=true,bool logit=true);//plot density of slice at c=h
		void plot_projection(float lmin,float lmax,char c,mgo::plt* pl,bool show=true,bool logit=true);//plot density projected along c axis
		void plot_dens(grid3&,mgo::plt*);//plot density ratio
		void x_hist(mgo::plt*,bool logit=false);//project density onto x axis
		void y_hist(mgo::plt*,bool logit=false);//project density onto y axis
*/
};

#endif
