/* \file    getPhi.cpp
    \author  James Binney & Eugene Vasiliev
    \date    2015-2023

    This is edited from the AGAMA distribution file example_self_consistent_model.cpp
    It determines the self-consistent potential of a model specified by
    a gas disc and DFs specified in m.ini
    It outputs the potential in 3 files that can be used to reproduce
    the potential by other programs.
    The function printoutInfo outputs various other diagnostics, including the
    circular-speed curve and density profiles
    The function writeHayden reads a file that contains the barycentre (R,z) of a
    spatial cell and the number of stars in that cell. It outputs for each cell
    samples of stars  and after computing the then for each cell drawn from the
    stellar components in the correct proportion.
    By calling writeSnapshot one can create N-body representations of all mass components: dark matter halo,
    stars (bulge, thin and thick disks and stellar halo combined), and gas disk.
*/
#include "galaxymodel_base.h"
#include "galaxymodel_selfconsistent.h"
#include "galaxymodel_velocitysampler.h"
#include "df_factory.h"
#include "potential_composite.h"
#include "potential_factory.h"
#include "potential_multipole.h"
#include "potential_utils.h"
#include "particles_io.h"
#include "math_core.h"
#include "math_spline.h"
#include "utils.h"
#include "utils_config.h"
#include "actions_staeckel.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <string.h>
#include <sys/stat.h>
#include "/u/c/spress/hacked/press.h"

using potential::PtrDensity;
using potential::PtrPotential;

// define internal unit system - arbitrary numbers here! the result should not depend on their choice
const units::InternalUnits intUnits(2.7183*units::Kpc, 3.1416*units::Myr);

// define external unit system describing the data (including the parameters in INI file)
const units::ExternalUnits extUnits(intUnits, 1.*units::Kpc, 1.*units::kms, 1.*units::Msun);

// used for outputting the velocity distribution (the value is read from the ini file)
double solarRadius = NAN;

// various auxiliary functions for printing out information are non-essential
// for the modelling itself; the essential workflow is contained in main()

/// print the rotation curve for a collection of potential components into a text file
void writeRotationCurve(const std::string& filename, const std::vector<PtrPotential>& potentials)
{
    std::ofstream strm(filename.c_str());
    strm << "# radius[Kpc]\tv_circ,total[km/s]\tdisk\tbulge\thalo\n";
    // print values at certain radii, expressed in units of Kpc
    std::vector<double> radii = math::createExpGrid(81, 0.01, 100);
    for(unsigned int i=0; i<radii.size(); i++) {
        strm << radii[i];  // output radius in kpc
        double v2sum = 0;  // accumulate squared velocity in internal units
        double r_int = radii[i] * intUnits.from_Kpc;  // radius in internal units
        std::string str;
        for(unsigned int j=0; j<potentials.size(); j++) {
            double vc = v_circ(*potentials[j], r_int);
            if(vc>0) v2sum += pow_2(vc); else v2sum -= pow_2(vc);
            str += "\t" + utils::toString(vc * intUnits.to_kms);  // output in km/s
	}
	double x = v2sum>=0? sqrt(v2sum) : -sqrt(fabs(v2sum));
        strm << '\t' << (x * intUnits.to_kms) << str << '\n';
    }
}

/// print surface density profiles to a file
void writeSurfaceDensityProfile(const std::string& filename, const galaxymodel::GalaxyModel& model)
{
    std::cout << "Writing surface density profile on "+filename+"\n";
    std::vector<double> radii;
    // convert radii to internal units
    for(double r=1./8; r<=30; r<1 ? r*=2 : r<16 ? r+=0.5 : r+=2)
        radii.push_back(r * intUnits.from_Kpc);
    int nr = radii.size();
    int nc = model.distrFunc.numValues();  // number of DF components
    std::vector<double> surfDens(nr*nc), rmsHeight(nr*nc), rmsVel(nr*nc);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for(int ir=0; ir<nr; ir++) {
	    computeProjectedMoments(model, radii[ir], &surfDens[ir*nc], &rmsHeight[ir*nc], &rmsVel[ir*nc],
				   NULL, NULL, NULL, /*separate*/ true);
    }

    std::ofstream strm(filename.c_str());
    strm << "# Radius[Kpc]\tsurfaceDensity[Msun/pc^2]\trmsHeight[kpc]\tsigma[km/s]\n";
    strm << nr << '\n';
    for(int ir=0; ir<nr; ir++){
        strm << radii[ir] * intUnits.to_Kpc << '\n';
	for(int ic=0; ic<nc; ic++)
		strm << surfDens[ir*nc+ic]  * intUnits.to_Msun_per_pc2 << '\t';
	strm << '\n';
	for(int ic=0; ic<nc; ic++)
		strm << rmsHeight[ir*nc+ic] * intUnits.to_Kpc << '\t';
	strm << '\n';
	for(int ic=0; ic<nc; ic++)
		strm << rmsVel[ir*nc+ic] * intUnits.to_kms << '\t';
	strm << '\n';
    }
}
/*
 Write into a file (x,v) of a sample of stars along the line of sight of an
 observer external to the galaxy
 */
void do_los(char* filename,const galaxymodel::GalaxyModel& model,obs::extLos* Los){
	size_t N=10000;
	std::vector<coord::PosVelCyl> sample=sampleLOS(model,Los,N);
	FILE* ofile;
	fopen_s(&ofile,filename,"w");
	for(size_t i=0;i<N;i++) fprintf(ofile,"%f %f %f %f %f %f\n",
					sample[i].R*intUnits.to_Kpc,
					sample[i].z*intUnits.to_Kpc,
					sample[i].phi,
					sample[i].vR*intUnits.to_kms,
					sample[i].vz*intUnits.to_kms,
					sample[i].vphi*intUnits.to_kms);
	fclose(ofile);
}
/*
 Write just (s,Vlos) for stars sampled on LOS
*/
void do_sVlos(char* filename,const galaxymodel::GalaxyModel& model,obs::extLos* Los){
	size_t N=10000;
	std::vector<std::pair<double,double> > sample=sampleLOSsVlos(model,Los,N);
	FILE* ofile;
	fopen_s(&ofile,filename,"w");
	for(size_t i=0;i<N;i++) fprintf(ofile,"%f %f\n",
					sample[i].first*intUnits.to_Kpc,
					sample[i].second*intUnits.to_kms);
	fclose(ofile);
}
/*
 Sample velocity space at locations specified in the input file /APOGEE/DR17_nos
 which specifies the number of real stars. The sample contains a multiple
 of this number distributed between components according to the model
 On return the filename contains the total masses of the components and
 their fractional contributions to each point. The sampled stars are
 left in the files filenamex.y with x
 the number of the location and y the component
*/
void writeHayden(const std::string& filename,const galaxymodel::GalaxyModel& model,
		 const std::vector<df::PtrDistributionFunction>& dfStellarArray)
{
	int nc = dfStellarArray.size();
	std::vector<double> Ms;
	for(int ic=0; ic<nc; ic++)
		Ms.push_back(dfStellarArray[ic]->totalMass());
	FILE *ofile;// whither to write masses of cpts
	fopen_s(&ofile,filename.c_str(),"w");
	for(size_t i=0; i<Ms.size(); i++) fprintf(ofile,"%g ",Ms[i]*intUnits.to_Msun);
	fclose(ofile);
	FILE *ifile;// whence to read star numbers in each cell
	int nh;
	if(fopen_s(&ifile,"/APOGEE/DR17_nos","r")){
		printf("I can't open /APOGEE/DR17_nos\n"); exit(0);
	}
	fscanf(ifile,"%d",&nh);// read # of cells
	std::vector<double> Rs(nh),zs(nh);// for their barycentres
	std::vector<int> N(nh);// for # of stars in cell
	for(int i=0;i<nh;i++){
		fscanf(ifile,"%lf %lf %d",&Rs[i],&zs[i],&N[i]);
		Rs[i]*=intUnits.from_Kpc; zs[i]*=intUnits.from_Kpc;
	}
	fclose(ifile);
	std::cout << "Writing Hayden data for " << nc <<" components at " << nh << " locations  on "+filename+"\n";
	std::vector<float> fracs(nh*nc);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for(int ih=0; ih<nh; ih++) {//run over locations
		coord::PosCyl pos(Rs[ih],zs[ih],0);
		double* dens = new double[nc];
		galaxymodel::computeMoments(model, coord::PosCyl(Rs[ih],zs[ih],0), dens, NULL, NULL,
					    NULL, NULL, NULL, NULL, true);
		double total=0;
		for(int ic=0; ic<nc; ic++){
			total+=dens[ic];
		}
		for(int ic=0; ic<nc; ic++) fracs[ih*nc+ic]=dens[ic]/total;
		for(int ic=0; ic<nc; ic++){//run over components
			size_t n=5*N[ih]*dens[ic]/total;
			char stuff[10]; sprintf(stuff,"%d_%d",ih,ic);
			std::string fname=filename+stuff;
			std::ofstream strm(fname.c_str());
			strm << Rs[ih]*intUnits.to_Kpc << " " << zs[ih]*intUnits.to_Kpc << " " << ic << " " << n << "\n";
			if(n>0){
				std::vector<coord::VelCyl> Vs = ic!=4?
					galaxymodel::sampleHalfVelocity(galaxymodel::GalaxyModel(
					model.potential, model.actFinder, *dfStellarArray[ic]),pos,n)
					: galaxymodel::sampleVelocity(galaxymodel::GalaxyModel(
					model.potential, model.actFinder, *dfStellarArray[ic]),pos,n);
				for(int i=0;i<n;i++){
					coord::PosVelCyl pv(pos,Vs[i]);
					actions::Actions J; double f;
					if(totalEnergy(model.potential,pv)>0){
						J.Jr=-1; J.Jz=-1; f=-1;
					} else {					
						J = model.actFinder.actions(pv);
						f = dfStellarArray[ic]->value(J);
					}
					strm << Vs[i].vR*intUnits.to_kms << " "
							<< Vs[i].vz*intUnits.to_kms << " "
							<< Vs[i].vphi*intUnits.to_kms << " "
							<< J.Jr*intUnits.to_Kpc*intUnits.to_kms<< " "
							<< J.Jz*intUnits.to_Kpc*intUnits.to_kms<< " "
							<< J.Jphi*intUnits.to_Kpc*intUnits.to_kms<< " "
							<< f <<"\n";
				}
			}
			strm.close();
		}
		delete[] dens;
	}
	fopen_s(&ofile,filename.c_str(),"a");
	fprintf(ofile,"\nFractions\n");
	for(int i=0; i<nh; i++){
		fprintf(ofile,"(%6.3f %6.3f) ",Rs[i]*intUnits.to_Kpc,zs[i]*intUnits.to_Kpc);
		for(int j=0; j<nc; j++) fprintf(ofile,"%6.3f ",fracs[i*nc+j]);
		fprintf(ofile,"\n");
	}
	fclose(ifile);
}
/// print radial density profile for several sub-components of the stellar DF
void writeRadialDensityProfile(const std::string& filename, const galaxymodel::GalaxyModel& model)
{
	std::cout << "Writing radial density profile on "+filename+"\n";
	std::vector<double> Rs;
    // convert height to internal units
	for(double R=0.5; R<=20; R<2.5 ? R+=0.25 : R+=1)
		Rs.push_back(R * intUnits.from_Kpc);
	double z = 0.0 * intUnits.from_Kpc;
	int nh = Rs.size();
	int nc = model.distrFunc.numValues();
	int n = nh*nc;
	std::vector<double> dens(n),vbar(n);
	std::vector<coord::Vel2Cyl> sigmas(n);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for(int ih=0; ih<nh; ih++) {
		int ic=ih*nc;
		//printf("Radius %d %f\n",ih,Rs[ih]*intUnits.to_Kpc); 
		computeMoments(model, coord::PosCyl(Rs[ih],z,0), &dens[ic], &vbar[ic],
			       &sigmas[ic], NULL, NULL, NULL, NULL, true);
	}
	std::ofstream strm(filename.c_str());
	strm << "# R[Kpc]\tdarkHalo\tThinDisk\tThickDisk\tStellarHalo\t[Msun/pc^3]\n";
	strm << nh << '\n';
	for(int ih=0; ih<nh; ih++) {
		strm << (Rs[ih] * intUnits.to_Kpc);
		for(int ic=0; ic<nc; ic++)
			strm << '\t' << (dens[ih*nc+ic] * intUnits.to_Msun_per_pc3);
		strm << '\n';
	}
	strm << "<Vphi>\n";
	for(int ih=0; ih<nh; ih++) {
		strm << (Rs[ih] * intUnits.to_Kpc);
		for(int ic=0; ic<nc; ic++)
			strm << '\t' << (vbar[ih*nc+ic] * intUnits.to_kms);
		strm << '\n';
	}
	strm << "sigma_R\n";
	for(int ih=0; ih<nh; ih++) {
		strm << (Rs[ih] * intUnits.to_Kpc);
		for(int ic=0; ic<nc; ic++)
			strm << '\t' << (sqrt(sigmas[ih*nc+ic].vR2) * intUnits.to_kms);
		strm << '\n';
	}
	strm << "sigma_z\n";
	for(int ih=0; ih<nh; ih++) {
		strm << (Rs[ih] * intUnits.to_Kpc);
		for(int ic=0; ic<nc; ic++)
			strm << '\t' << (sqrt(sigmas[ih*nc+ic].vz2) * intUnits.to_kms);
		strm << '\n';
	}
	strm << "sigma_phi\n";
	for(int ih=0; ih<nh; ih++) {
		strm << (Rs[ih] * intUnits.to_Kpc);
		for(int ic=0; ic<nc; ic++)
			strm << '\t' << (sqrt(sigmas[ih*nc+ic].vphi2-pow_2(vbar[ih*nc+ic])) * intUnits.to_kms);
		strm << '\n';
	}
}

/// print vertical density profile for several sub-components of the stellar DF
void writeVerticalDensityProfile(const std::string& filename, double frac,
				 const galaxymodel::GalaxyModel& model)
{
	std::cout << "Writing vertical density profile on "+filename+"\n";
	std::vector<double> heights;
    // convert height to internal units
	for(double h=0; h<=3; h<.5 ? h+=0.05 : h+=0.25)
		heights.push_back(h * intUnits.from_Kpc);
	double R = frac * solarRadius * intUnits.from_Kpc;
	int nh = heights.size();
	int nc = model.distrFunc.numValues();
	int n = nh*nc;
	std::vector<double> dens(n),vbar(n);
	std::vector<coord::Vel2Cyl> sigmas(n);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for(int ih=0; ih<nh; ih++) {
//		printf("z=%5.2f\n",heights[ih]*intUnits.to_Kpc);
		computeMoments(model, coord::PosCyl(R,heights[ih],0), &dens[ih*nc], &vbar[ih*nc], &sigmas[ih*nc],
			       NULL, NULL, NULL, NULL, /*separate*/ true);
	}

	std::ofstream strm(filename.c_str());
	strm << "# z[Kpc]\tDM\tbulge\tThinDisk\tThickDisk\tStellarHalo[Msun/pc^3] at\n" << R*intUnits.to_Kpc;
	strm << "  " << nh << '\n';
	for(int ih=0; ih<nh; ih++) {
		strm << (heights[ih] * intUnits.to_Kpc);
		for(int ic=0; ic<nc; ic++)
			strm << '\t' << (dens[ih*nc+ic] * intUnits.to_Msun_per_pc3);
		strm << '\n';
	}
	strm << "<Vphi>\n";
	for(int ih=0; ih<nh; ih++) {
		strm << (heights[ih] * intUnits.to_Kpc);
		for(int ic=0; ic<nc; ic++)
			strm << '\t' << (vbar[ih*nc+ic] * intUnits.to_kms);
		strm << '\n';
	}
	strm << "sigma_R	<Omega.r>\n";
	for(int ih=0; ih<nh; ih++) {
		strm << (heights[ih] * intUnits.to_Kpc);
		for(int ic=0; ic<nc; ic++)
			strm << '\t' << (sqrt(sigmas[ih*nc+ic].vR2) * intUnits.to_kms);
		strm << '\n';
	}
	strm << "sigma_z\n";
	for(int ih=0; ih<nh; ih++) {
		strm << (heights[ih] * intUnits.to_Kpc);
		for(int ic=0; ic<nc; ic++)
			strm << '\t' << (sqrt(sigmas[ih*nc+ic].vz2) * intUnits.to_kms);
		strm << '\n';
	}
	strm << "sigma_phi\n";
	for(int ih=0; ih<nh; ih++) {
		strm << (heights[ih] * intUnits.to_Kpc);
		for(int ic=0; ic<nc; ic++)
			strm << '\t' << (sqrt(sigmas[ih*nc+ic].vphi2-pow_2(vbar[ih*nc+ic])) * intUnits.to_kms);
		strm << '\n';
	}
	strm << "sigma2_Rz\n";
	for(int ih=0; ih<nh; ih++) {
		strm << (heights[ih] * intUnits.to_Kpc);
		for(int ic=0; ic<nc; ic++)
			strm << '\t' << (sigmas[ih*nc+ic].vRvz * intUnits.to_kms);
		strm << '\n';
	}
}

/// print velocity distributions at the given point to a file
void writeVelocityDistributions(const std::string& filename, const double RoverRsun,
				const galaxymodel::GalaxyModel& model)
{
	const coord::PosCyl point(RoverRsun*solarRadius * intUnits.from_Kpc, 0.015 * intUnits.from_Kpc, 0);
	std::cout << "Writing velocity distributions at "
			"(R=" << point.R * intUnits.to_Kpc << ", z=" << point.z * intUnits.to_Kpc << ")\n";
    // create grids in velocity space
	double vR_max = 100 * intUnits.from_kms, vR_step=   1 * intUnits.from_kms;
	double vz_max =  50 * intUnits.from_kms, vz_step=   0.5 * intUnits.from_kms;
	double vp_max = 350 * intUnits.from_kms;
	std::vector<double> gridvR   = math::createSymmetricGrid(75, vR_step, vR_max);
	std::vector<double> gridvz   = math::createSymmetricGrid(75, vz_step, vz_max);
	std::vector<double> gridvphi = math::createUniformGrid(75, -.25*vp_max, vp_max);
	std::vector<double> amplvR, amplvz, amplvphi;
	double density;
    // compute the distributions
	const int ORDER = 3;
	math::BsplineInterpolator1d<ORDER> intvR(gridvR), intvz(gridvz), intvphi(gridvphi);
	galaxymodel::computeVelocityDistributionO3(model, point, false /*not projected*/,
		gridvR, gridvz, gridvphi, /*output*/ &density, &amplvR, &amplvz, &amplvphi);

	std::ofstream strm(filename.c_str());
	strm << "# V\tf(V_R)\tf(V_z) [1/(km/s)]\n";
	for(int i=-100; i<=100; i++) {
		double v = i*vR_max/100;
	// unit conversion: the VDF has a dimension 1/V, so that \int f(V) dV = 1;
	// therefore we need to multiply it by 1/velocityUnit
		strm << utils::toString(v * intUnits.to_kms)+'\t'+
				utils::toString(intvR.  interpolate(v, amplvR)   / intUnits.to_kms)+'\t'+
				utils::toString(intvz.  interpolate(v*vz_max/vR_max, amplvz)   / intUnits.to_kms)+'\n';
	}
	for(int i=-25; i<=100; i++) {
		double v = i*vp_max/100;
		strm << utils::toString(v * intUnits.to_kms)+'\t'+
				utils::toString(intvphi.interpolate(v, amplvphi) / intUnits.to_kms)+'\n';
	}
}
// write to file velocity distributions at locations read from file
void writeVelocityDistributions(const std::string &dir, const std::string &dname, const galaxymodel::GalaxyModel& model)
{
	FILE *ifile; int err=fopen_s(&ifile,dname.c_str(),"r");
	if(err){
		printf("I cannot find %s\n",dname.c_str()); return;
	}
	std::cout << "Computing velocity distributions\n";
	int Nloc=0; float Rloc,zloc;
	std::vector<coord::PosCyl> points;
	while(Nloc>=0){
		if(2!=fscanf(ifile,"%f %f",&Rloc,&zloc)) break;
		const coord::PosCyl point(Rloc * intUnits.from_Kpc, zloc * intUnits.from_Kpc, 0);
		points.push_back(point);
		Nloc++;
	}
	fclose(ifile);
	printf("%d locations read from file\n",Nloc);
    // create grids in velocity space
	double vR_max = 120 * intUnits.from_kms;
	double vp_max = 350 * intUnits.from_kms;
	std::vector<double> gridvR   = math::createUniformGrid(75, -vR_max, vR_max);
	std::vector<double> gridvz   = gridvR;  // for simplicity use the same grid for two dimensions
	std::vector<double> gridvphi = math::createUniformGrid(75, -.25*vp_max, vp_max);
	const unsigned int nc_max=7;
	unsigned int nc=model.distrFunc.numValues();
	if(nc>nc_max){
		printf("Error: nc_max (%d) < nc (%d)",nc_max,nc); exit(0);
	}
	double ***Ramps,***zamps,***pamps;
	Ramps=dmatrix(Nloc,nc,78); zamps=dmatrix(Nloc,nc,78); pamps=dmatrix(Nloc,nc,78);
	const int ORDER = 3;
	math::BsplineInterpolator1d<ORDER> intvR(gridvR), intvz(gridvz), intvphi(gridvphi);
	double **density; density=dmatrix(Nloc,nc);
	int **Ni; Ni=imatrix(Nloc,nc);
    // compute the distributions
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for(int loc=0; loc<Nloc; loc++){
		std::vector<double> amplvR[nc_max], amplvz[nc_max], amplvphi[nc_max];
		double* dens=new double[nc];
		coord::PosCyl pt(points[loc]);
		galaxymodel::computeVelocityDistributionO3(model, pt, false /*not projected*/,
			gridvR, gridvz, gridvphi, /*output*/ dens, amplvR, amplvz, amplvphi, /*separate*/ true);
		if(amplvR[0].size()>78) printf("Error size: %d\n",(int)amplvR[0].size());
		for(unsigned int ic=0; ic<nc; ic++){
			density[loc][ic]=dens[ic];
			Ni[loc][ic]=amplvR[ic].size();
			for(unsigned i=0; i<Ni[loc][ic]; i++){
				Ramps[loc][ic][i]=amplvR[ic][i];
				zamps[loc][ic][i]=amplvz[ic][i];
				pamps[loc][ic][i]=amplvphi[ic][i];
			}
		}
		delete[] dens;
	}
	printf("profiles computed\n");
#ifdef _OPENMP
//#pragma omp parallel for schedule(dynamic)
#endif
	for(int loc=0; loc<Nloc; loc++){
		char fname[50],locname[30];
		strcpy(fname,dir.c_str());
		if(Nloc<12 )sprintf(locname,"loc_%d",loc); else sprintf(locname,"loc%d",loc);
		strcat(fname,locname);
		//printf("opening %s\n",fname);
		FILE *ofile;
		if(fopen_s(&ofile,fname,"w")) printf("cannot open %s\n",fname);
		for(unsigned int ic=0; ic<nc; ic++) fprintf(ofile,"%g ",density[loc][ic]);
		fprintf(ofile,"\n");
		for(int i=-100; i<=100; i++) {
			double v = i*vR_max/100;
			fprintf(ofile,"%f",v * intUnits.to_kms);
			for(unsigned int ic=0;ic<nc;ic++){
				std::vector<double> amplvR, amplvz;
				for(unsigned i=0; i<Ni[loc][ic]; i++){
					amplvR.push_back(Ramps[loc][ic][i]); amplvz.push_back(zamps[loc][ic][i]);
				}
	// unit conversion: the VDF has a dimension 1/V, so that \int f(V) dV = 1;
	// therefore we need to multiply it by 1/velocityUnit
				fprintf(ofile," %g %g",
					  intvR.interpolate(v, amplvR)   / intUnits.to_kms,
					  intvz.interpolate(v, amplvz)   / intUnits.to_kms);
			}
			fprintf(ofile,"\n");
		}
		for(int i=-25; i<=100; i++) {
			double v = i*vp_max/100;
			fprintf(ofile,"%f",v * intUnits.to_kms);
			for(unsigned int ic=0;ic<nc;ic++){
				std::vector<double> amplvphi;
				for(unsigned i=0; i<Ni[loc][ic]; i++){
					amplvphi.push_back(pamps[loc][ic][i]);
				}
				fprintf(ofile," %g",intvphi.interpolate(v, amplvphi) / intUnits.to_kms);
			}
			fprintf(ofile,"\n");
		}
		fclose(ofile);
	}
	delmatrix(Ramps,Nloc,nc); delmatrix(zamps,Nloc,nc); delmatrix(pamps,Nloc,nc);
	delmatrix(density,Nloc); delmatrix(Ni,Nloc);
	printf("Profiles left in %s\n",dir.c_str());
}

/// report progress after an iteration
void printoutInfo(const galaxymodel::SelfConsistentModel& model, const std::string& dir)
{
	const potential::BaseDensity& compDisk = *model.components[0]->getDensity();
	//const potential::BaseDensity& compBulge= *model.components[1]->getDensity();
	const potential::BaseDensity& compHalo = *model.components[1]->getDensity();
	const potential::BaseDensity& compGas= *model.components[2]->getDensity();
	coord::PosCyl pt0(solarRadius * intUnits.from_Kpc, 0, 0);
	coord::PosCyl pt1(solarRadius * intUnits.from_Kpc, 1 * intUnits.from_Kpc, 0);
	//printf("Entered printoutInfo; R0 = %f Mgas = %g\n",solarRadius,compGas.totalMass()*intUnits.to_Msun);
	std::cout <<
			"Disk total mass="      << (compDisk.totalMass()  * intUnits.to_Msun) << " Msun"
			", rho(Rsolar,z=0)="    << (compDisk.density(pt0) * intUnits.to_Msun_per_pc3) <<
			", rho(Rsolar,z=1kpc)=" << (compDisk.density(pt1) * intUnits.to_Msun_per_pc3) << " Msun/pc^3\n"
			"Halo total mass="      << (compHalo.totalMass()  * intUnits.to_Msun) << " Msun"
			", rho(Rsolar,z=0)="    << (compHalo.density(pt0) * intUnits.to_Msun_per_pc3) <<
			", rho(Rsolar,z=1kpc)=" << (compHalo.density(pt1) * intUnits.to_Msun_per_pc3) << " Msun/pc^3\n"
			"Potential at origin=-("<<
			(sqrt(-model.totalPotential->value(coord::PosCyl(0,0,0))) * intUnits.to_kms) << " km/s)^2"
			", total mass=" << (model.totalPotential->totalMass() * intUnits.to_Msun) << " Msun\n";
	writeDensity(dir + "dens_disk", compDisk, extUnits);
//	writeDensity(dir + "dens_bulge", compBulge, extUnits);
	writeDensity(dir + "dens_halo", compHalo, extUnits);
	writePotential(dir + "potential", *model.totalPotential, extUnits);
	std::vector<PtrPotential> potentials(2);
	potentials[0] = dynamic_cast<const potential::CompositeCyl&>(*model.totalPotential).component(1);
//	potentials[1] = potential::Multipole::create(compBulge, /*lmax*/6, /*mmax*/0, /*gridsize*/25);
	potentials[1] = potential::Multipole::create(compHalo,  /*lmax*/6, /*mmax*/0, /*gridsize*/25);
	writeRotationCurve(dir + "rotcurve", potentials);
}

/// perform one iteration of the model
void doIteration(galaxymodel::SelfConsistentModel& model,const int iterationIndex, const std::string &dir)
{
    std::cout << "Starting iteration #" << iterationIndex << "\n";
    bool error=false;
    try {
	    galaxymodel::doIteration(model);
    }
    catch(std::exception& ex) {
        error=true;  // report the error and allow to save the results of the last iteration
        std::cout << "==== Exception occurred: \n" << ex.what();
    }
    printoutInfo(model, dir);
    if(error)
        exit(1);  // abort in case of problems
}
int main(int narg,char **args)
{
	if(narg!=3){
		printf("You must specify #iteration & directory\n"); return 0;
	}
	printf("|%s\n",args[2]);
	std::time_t now, start_t = std::time(NULL);
	//Specify what outputs you require
	bool doHayden=true, vertProfile=false, Nbody=false, newModel=true;
	int noIterations;
	sscanf(args[1],"%d",&noIterations); if(noIterations==0) newModel=false;
	std::string Mno(args[2]);
	const std::string dir = Mno+"/";
    // read parameters from the INI file
	const std::string iniFileName = dir + "m.ini";
	utils::ConfigFile ini(iniFileName);
	utils::KeyValueMap
		    iniPotenThinDisk = ini.findSection("Potential thin disk"),
    iniPotenThickDisk= ini.findSection("Potential thick disk"),
    iniPotenGasDisk  = ini.findSection("Potential gas disk"),
    iniPotenBulge    = ini.findSection("Potential bulge"),
    iniPotenDarkHalo = ini.findSection("Potential dark halo"),
	iniDFyoungDisk   = ini.findSection("DF young disk"),
	iniDFmiddleDisk  = ini.findSection("DF middle disk"),
	iniDFoldDisk     = ini.findSection("DF old disk"),
	iniDFhighADisk   = ini.findSection("DF highA disk"),
        iniDFStellarHalo = ini.findSection("DF stellar halo"),
        iniDFBulge       = ini.findSection("DF bulge"),
	iniDFDarkHalo    = ini.findSection("DF dark halo"),
        iniSCMDisk       = ini.findSection("SelfConsistentModel disk"),
        iniSCMBulge      = ini.findSection("SelfConsistentModel bulge"),
        iniSCMHalo       = ini.findSection("SelfConsistentModel halo"),
        iniSCM           = ini.findSection("SelfConsistentModel");
    if(!iniSCM.contains("rminSph")) {  // most likely file doesn't exist
        std::cout << "Invalid INI file " << iniFileName << "\n";
        return -1;
    }
    solarRadius = ini.findSection("Data").getDouble("SolarRadius", solarRadius);

    // set up parameters of the entire Self-Consistent Model
    galaxymodel::SelfConsistentModel model;
    model.rminSph         = iniSCM.getDouble("rminSph") * extUnits.lengthUnit;
    model.rmaxSph         = iniSCM.getDouble("rmaxSph") * extUnits.lengthUnit;
    model.sizeRadialSph   = iniSCM.getInt("sizeRadialSph");
    model.lmaxAngularSph  = iniSCM.getInt("lmaxAngularSph");
    model.RminCyl         = iniSCM.getDouble("RminCyl") * extUnits.lengthUnit;
    model.RmaxCyl         = iniSCM.getDouble("RmaxCyl") * extUnits.lengthUnit;
    model.zminCyl         = iniSCM.getDouble("zminCyl") * extUnits.lengthUnit;
    model.zmaxCyl         = iniSCM.getDouble("zmaxCyl") * extUnits.lengthUnit;
    model.sizeRadialCyl   = iniSCM.getInt("sizeRadialCyl");
    model.sizeVerticalCyl = iniSCM.getInt("sizeVerticalCyl");
    model.useActionInterpolation = iniSCM.getBool("useActionInterpolation");

    // initialize density profiles of various components
    std::vector<PtrDensity> densityStellarDisk(3);
    PtrDensity densityDarkHalo = potential::createDensity(iniPotenDarkHalo, extUnits);
    densityStellarDisk[0]      = potential::createDensity(iniPotenThinDisk, extUnits);
    densityStellarDisk[1]      = potential::createDensity(iniPotenThickDisk,extUnits);
    densityStellarDisk[2]      = potential::createDensity(iniPotenBulge,    extUnits);
    PtrDensity densityGasDisk  = potential::createDensity(iniPotenGasDisk,  extUnits);

    // add components to SCM - at first, all of them are static density profiles
    model.components.push_back(galaxymodel::PtrComponent(
	    new galaxymodel::ComponentStatic(PtrDensity(
	    new potential::CompositeDensity(densityStellarDisk)), true)));
    model.components.push_back(galaxymodel::PtrComponent(
	    new galaxymodel::ComponentStatic(densityDarkHalo, false)));
    model.components.push_back(galaxymodel::PtrComponent(
	    new galaxymodel::ComponentStatic(densityGasDisk, true)));

    // initialize total potential of the model (first guess)
    updateTotalPotential(model);
    printoutInfo(model, "init");
    printf("back from printoutInfo\n");

    std::cout << "**** STARTING MODELLING ****\nInitial masses of density components: "
		    "Mdisk + MstellHalo + Mbulge="  << (model.components[0]->getDensity()->totalMass() * intUnits.to_Msun) << " Msun, "
		    "MdarHalo="  << (densityDarkHalo->totalMass() * intUnits.to_Msun) << " Msun, "
		    "Mgas="   << (densityGasDisk->totalMass() * intUnits.to_Msun) << " Msun\n";

    // create the dark halo DF
    df::PtrDistributionFunction dfHalo = df::createDistributionFunction(
	    iniDFDarkHalo, model.totalPotential.get(), NULL, extUnits);
    printf("dark halo DF created\n");
    // same for the bulge
    df::PtrDistributionFunction dfBulge = df::createDistributionFunction(
	    iniDFBulge, model.totalPotential.get(), NULL, extUnits);
    printf("bulge DF created..");
    // same for the stellar components (thin/thick disks and stellar halo)
    std::vector<df::PtrDistributionFunction> dfStellarArray;
    dfStellarArray.push_back(df::createDistributionFunction(
	    iniDFyoungDisk, model.totalPotential.get(), NULL, extUnits));
    printf("young disc DF created..");
    dfStellarArray.push_back(df::createDistributionFunction(
	    iniDFmiddleDisk, model.totalPotential.get(), NULL, extUnits));
    printf("middle disc DF created..");
    dfStellarArray.push_back(df::createDistributionFunction(
	    iniDFoldDisk, model.totalPotential.get(), NULL, extUnits));
    printf("old disc DF created..");
    dfStellarArray.push_back(df::createDistributionFunction(
	    iniDFhighADisk, model.totalPotential.get(), NULL, extUnits));
    printf("highA disc DF created\n");
    dfStellarArray.push_back(df::createDistributionFunction(
	    iniDFStellarHalo, model.totalPotential.get(), NULL, extUnits));
    printf("Stellar halo DF created\n");
    dfStellarArray.push_back(dfBulge);
// composite DF of all stellar components
    df::PtrDistributionFunction dfStellar(new df::CompositeDF(dfStellarArray));
    printf("Total stellar mass: %g Msun\n",dfStellar->totalMass()*intUnits.to_Msun);
//assemble total DF
    std::vector<df::PtrDistributionFunction> dfAllArray;
    dfAllArray.push_back(dfHalo);
    for(int i=0; i<(int)dfStellarArray.size(); i++) dfAllArray.push_back(dfStellarArray[i]);
    df::PtrDistributionFunction dfAll(new df::CompositeDF(dfAllArray));
    printf("Inserting DFs into galaxyModel. ");
    int ncmp=0;
    // replace the static disk density component of SCM with a DF-based disk component
    model.components[ncmp] = galaxymodel::PtrComponent(
	    new galaxymodel::ComponentWithDisklikeDF(dfStellar, PtrDensity(),
	    iniSCMDisk.getInt("mmaxAngularCyl"),
	    iniSCMDisk.getInt("sizeRadialCyl"),
	    iniSCMDisk.getDouble("RminCyl") * extUnits.lengthUnit,
	    iniSCMDisk.getDouble("RmaxCyl") * extUnits.lengthUnit,
	    iniSCMDisk.getInt("sizeVerticalCyl"),
	    iniSCMDisk.getDouble("zminCyl") * extUnits.lengthUnit,
	    iniSCMDisk.getDouble("zmaxCyl") * extUnits.lengthUnit));
    printf("disc's, halo & bulge in...");
// same for the dark halo
    ncmp++;
    model.components[ncmp] = galaxymodel::PtrComponent(
	    new galaxymodel::ComponentWithSpheroidalDF(dfHalo, potential::PtrDensity(),
	    iniSCMHalo.getInt("lmaxAngularSph"),
	    iniSCMHalo.getInt("mmaxAngularSph"),
	    iniSCMHalo.getInt("sizeRadialSph"),
	    iniSCMHalo.getDouble("rminSph") * extUnits.lengthUnit,
	    iniSCMHalo.getDouble("rmaxSph") * extUnits.lengthUnit));
    printf("DM's in\n");
    
    // we can compute the masses even though we don't know the density profile yet
    std::cout <<
        "Masses of DF components:"
		    "\n Mdisk + Mstel.halo + M.bulge=" << (dfStellar->totalMass() * intUnits.to_Msun) <<
		    " Msun\n (Myoung=" << (dfStellarArray[0]->totalMass() * intUnits.to_Msun) <<
		    ", Mmiddle="     << (dfStellarArray[1]->totalMass() * intUnits.to_Msun) <<
		    ", Mold="     << (dfStellarArray[2]->totalMass() * intUnits.to_Msun) <<
		    ", Mthick="     << (dfStellarArray[3]->totalMass() * intUnits.to_Msun);
    std::cout <<
		    ", Mstel.halo=" << (dfStellarArray[4]->totalMass() * intUnits.to_Msun) <<
		    ")\n Mbulge="    << (dfBulge->totalMass() * intUnits.to_Msun) << " Msun"
		    "; Mdark="      << (dfHalo ->totalMass() * intUnits.to_Msun) << " Msun\n";
    if(newModel){
	    if(utils::fileExists(dir + "potential"))
		    model.totalPotential = potential::readPotential(dir + "potential", extUnits);
	    std::cout << "Potential loaded: value at origin=-("<<
			    (sqrt(-model.totalPotential->value(coord::PosCyl(0,0,0))) * intUnits.to_kms) << " km/s)^2\n";
	    // update the action finder
	    std::cout << "Updating action finder..."<<std::flush;
	    model.actionFinder.reset(new actions::ActionFinderAxisymFudge(model.totalPotential,
		    model.useActionInterpolation));
	    std::cout << "done"<<std::endl;
	    std::ofstream strm0(dir + "masses");
	    strm0 << "\n Mdisk + Mstel.halo=" << (dfStellar->totalMass() * intUnits.to_Msun) <<
			    " Msun\n (Myoung=" << (dfStellarArray[0]->totalMass() * intUnits.to_Msun) <<
			    ", Mmiddle="     << (dfStellarArray[1]->totalMass() * intUnits.to_Msun) <<
			    ", Mold="     << (dfStellarArray[2]->totalMass() * intUnits.to_Msun) <<
			    ", MHighA="     << (dfStellarArray[3]->totalMass() * intUnits.to_Msun) <<
			    ", Mstel.halo=" << (dfStellarArray[4]->totalMass() * intUnits.to_Msun) <<
			    ")\n Mbulge="    << (dfBulge->totalMass() * intUnits.to_Msun) << " Msun"
			    "; Mdark="      << (dfHalo ->totalMass() * intUnits.to_Msun) << " Msun\n";
	    strm0.close();
    // do a few more iterations to obtain the self-consistent density profile for both disks
	    for(int iteration=1; iteration<=noIterations; iteration++)
		    doIteration(model, iteration, dir);
	    now = std::time(NULL);
	    std::cout << now-start_t << " secs to build\n";
	    start_t=now;
    } else {
	    model.totalPotential = potential::readPotential(dir + "potential", extUnits);
	    // update the action finder
	    std::cout << "Updating action finder..."<<std::flush;
	    model.actionFinder.reset(new actions::ActionFinderAxisymFudge(model.totalPotential,
		    model.useActionInterpolation));
	    std::cout << "done"<<std::endl;
    }
    double dR=.1, Phi; coord::GradCyl K;
    model.totalPotential->eval(coord::PosCyl(solarRadius*intUnits.from_Kpc,1.1*intUnits.from_Kpc,0),
			       &Phi, &K);
    double Vc1=v_circ(*(model.totalPotential),(solarRadius+dR)*intUnits.from_Kpc);
    double Vc0=v_circ(*(model.totalPotential),(solarRadius-dR)*intUnits.from_Kpc);
    double dVcsqdR=.5*(pow_2(Vc1)-pow_2(Vc0))/(dR*intUnits.from_Kpc);//internal units
    printf("Vc etc %f %f %f\n",Vc1,Vc0,dVcsqdR);
    printf("Kz(1.1)   : %f km/s/Myr\n",K.dz*intUnits.to_kms/intUnits.to_Myr);
    //convert to astronomical units
    K.dz*=intUnits.to_kms*1e5/(intUnits.to_Myr*units::Myr)*pow_2(units::pc)/(2*M_PI*units::Grav*units::Msun);
    dVcsqdR*=intUnits.to_kms*1e5/(intUnits.to_Myr*units::Myr)*pow_2(units::pc)/(2*M_PI*units::Grav*units::Msun);
    printf("Kz/2piG   : %f Msun/pc^2\n",K.dz);
    printf("Sigma     : %f Msun/pc^2\n",K.dz+(1.1/solarRadius)*dVcsqdR);
    // output various profiles
    galaxymodel::GalaxyModel modelStars(*model.totalPotential, *model.actionFinder, *dfStellar);
    galaxymodel::GalaxyModel modelAll(*model.totalPotential, *model.actionFinder, *dfAll);

    if(doHayden)
	    writeHayden(dir + "Hayden", modelStars, dfStellarArray);
    if(newModel)
	    writeSurfaceDensityProfile(dir + "stellar.surfdens", modelStars);
    if(newModel)
	    writeRadialDensityProfile(dir + "stellar_dens.radial", modelAll);
    if(vertProfile){
	    double RoverRsun=1;
	    writeVerticalDensityProfile(dir + "dens.vertical", RoverRsun, modelAll);
    }
	    if(newModel){
		    writePotential(dir + "potential", *model.totalPotential, extUnits);
		    printf("Potential written\n");
	    }
	    now = std::time(NULL);
	    std::cout << now-start_t << " secs for diagnostics\n";
    if(Nbody){
	    std::cout << "Writing N-body sample of solar neighbourhood\n";
    // export model to an N-body snapshot
	    std::cout << "Creating an N-body representation of the model\n";
	    std::string format = "CMP";   // could use "text", "nemo" or "gadget" here

    // first create a representation of density profiles without velocities
    // (just for demonstration), by drawing samples from the density distribution
	    std::cout << "Writing N-body sampled density profile for the dark matter halo\n";
	    particles::writeSnapshot("dens_dm_final", galaxymodel::sampleDensity(
		    *model.components[ncmp]->getDensity(), 800000), format, extUnits);
	    std::cout << "Writing N-body sampled density profile for the stellar bulge, disk and halo\n";
	    std::vector<PtrDensity> densityStars(ncmp);
	    densityStars[0] = model.components[0]->getDensity();  // stellar disks and halo
	    if(ncmp==2); densityStars[1] = model.components[1]->getDensity();  // bulge
	    particles::writeSnapshot("dens_stars_final", galaxymodel::sampleDensity(
		    potential::CompositeDensity(densityStars), 200000), format, extUnits);

    // now create genuinely self-consistent models of both components,
    // by drawing positions and velocities from the DF in the given (self-consistent) potential
	    std::cout << "Writing a complete DF-based N-body model for the dark matter halo\n";
	    particles::writeSnapshot("model_dm_final", galaxymodel::samplePosVel(
		    galaxymodel::GalaxyModel(*model.totalPotential, *model.actionFinder, *dfHalo),
		    800000), format, extUnits);
	    std::cout << "Writing a complete DF-based N-body model for the stellar bulge, disk and halo\n";
	    //dfStellarArray.push_back(dfBulge); now already included
	    //dfStellar.reset(new df::CompositeDF(dfStellarArray));  // all stellar components incl. bulge
	    particles::writeSnapshot("model_stars_final", galaxymodel::samplePosVel(
		    galaxymodel::GalaxyModel(*model.totalPotential, *model.actionFinder, *dfStellar), 200000),
				     format, extUnits);
    // we didn't use an action-based DF for the gas disk, leaving it as a static component;
    // to create an N-body representation, we sample the density profile and assign velocities
    // from the axisymmetric Jeans equation with equal velocity dispersions in R,z,phi
	    std::cout << "Writing an N-body model for the gas disk\n";
	    particles::writeSnapshot("model_gas_final", galaxymodel::assignVelocity(
		    galaxymodel::sampleDensity(*model.components[2]->getDensity(), 24000),
	/*parameters for the axisymmetric Jeans velocity sampler*/
		    *model.components[3]->getDensity(), *model.totalPotential, /*beta*/ 0., /*kappa*/ 1.),
				     format, extUnits);
    }
    std::cout << "Results left in " + dir << '\n'; 
}
