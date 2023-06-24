#ifndef __grid__
#define __grid__
#include "grid.h"
#endif
#ifndef __grid2__
#define __grid2__
class grid2{
	private:
		double w(double,double);
		void wts(double,double,double&,double&,double&,double&);
		void find_peak(float*,float*);
		float **m,***av;
		double x0,y0,Deltax,Deltay;
		int nx,ny,na;
		bool greyit;
	public:
		grid2(void){};
		grid2(int,int,int,double,double,double,double);
		grid2(grid2&,bool copy=true);//copy grid & its data
//		grid2(FILE*,int=0);
		~grid2(void);
		int navs(){
			return na;
		}
		double CellVolume(void){
			return Deltax*Deltay;
		}
		void get_centres(float*,float*);//returns cell centres
		void setup(int,int,int,double,double,double,double);
		void CIC(double,double,double);//add contrib to density
		void CIC(double,double*,double,double);//add contrib to density and averages
		double test_CIC(double,double*,double,double);
		double dens(double,double);//recover density at (x,y)
		double dens(double,double,double*);//recover density and averages at (x,y)
		void get_minmax(float&,float&) const;//return lmin, lmax
		void get_minmax(int,float&,float&) const;//return min, max of kth value
		void get_ints(int& _nx,int& _ny,int& _na){
			_nx=nx; _ny=ny; _na=na;
		}
		void get_values(int,float,float**) const;// recovers totals (int<0) or average of kth value; float if m=0
		void print_dens(void);
		float get_total(void);//returns total number on grid
		void fillup(void);//patches gaps
//		void dump(FILE*);//write m and av to file
		void zero(void);//initialise density and averages
		void filter(void);
		void disgorge(int&,int&,int&,double&,double&,double&,double&);//returns nx,ny,x0,y0,Deltas
		void regrid(grid2&,int,int);//put data from slot in into slot out 
		void divide(grid2&);
		void subtract(grid2&);//
		float reveal(int i,int j,float* avij);
		void add(grid2& grd);
		void toggle_grey(void);
// Now routines involving mgo::ply
/*		std::pair<float,float> plot_dens(mgo::plt*,bool show=true,bool logit=true);//plot density
		void plot_dens(grid2&,mgo::plt*);//plot density ratio
		void plot_dens(float,float,mgo::plt*,bool show=true,bool logit=true);//plot density with set limits
		void plot_densR(mgo::plt*,bool show=true,bool logit=true);//as above reversed colours
		void plot(float,mgo::plt*);//as above plus contour at float 
		void plot_fn(float (*fn)(float,float),mgo::plt*,bool);//plots fn(x,y) 
		void contour_fn(float (*fn)(float,float),float*,int,mgo::plt*,bool show=true);//contour fn(x,y) 
		void plot_fn_values(int,float,float (*)(float),mgo::plt*,bool show=true);//plot an average
		void contour_values(int,int,float,mgo::plt*);//contour an average, blank below float
		void contour_values(int,float*,int,float,mgo::plt*);//contour an average, blank below float
		void percentiles(int,float*,mgo::plt*);
		void plot_dispersion(int i,int j,float blank,mgo::plt*,bool show=true) const;//plot dispersion with av[i]=<x>, av[j]=<x^2>
		void plot_values(int,float,mgo::plt*,bool show=true,bool logit=true) const;//plot an average with float value of empty cells
		void plot_valuesR(int,float,mgo::plt*,bool show=true,bool logit=true) const;//plot an average with float value of empty cells
		void plot_values(float,float,int,float,mgo::plt*,bool show=true,bool logit=true) const;//plot an average with lmin,lmax first args
		void plot_valuesR(float,float,int,float,mgo::plt*,bool show=true,bool logit=true) const;//as above reversed colours
		void plot_fn_values(float,float,int,float,float (*)(float),mgo::plt*,bool show=true);//plot an average with lmin,lmax first args
		void contours(float*,int,mgo::plt*) const;//contours of values
		void plot_contour(float,mgo::plt*);
		void x_hist(mgo::plt* pl,bool logit=false);//project density onto x axis
		void y_hist(mgo::plt* pl,bool logit=false);//project density onto y axis
		void x_hist(FILE*,mgo::plt* pl,bool logit=false);
		void y_hist(FILE*,mgo::plt* pl,bool logit=false);
		void test_pl(float,float,char*,mgo::plt*);
*/
};

#endif
