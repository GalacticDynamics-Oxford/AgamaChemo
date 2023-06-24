/* \file    relax.cpp
    \author  James Binney & Eugene Vasiliev
    \date    2022 - 2023

   Nelder-Mead is used to ajust DFs (at constant Phi) and a chemical model to fit
   observed kinematics in 72 spatial bins and the chemistry in 30 bins in J space.
   It reads a potential and star samples written by getPhi.cpp. The DFs are specified
   by m.ini and on exit a prposal for new DFs is written to dyn_params.dat in a form
   that's ready to edit in to m.ini prior to interating the potential to self-consistency
   with get,Phi.cpp.
*/


#include "relax.h"

const int nxJ=60,nyJ=60;
const float Jphi_min=-1000,Jphi_max=4000,Jz_min=0,Jz_max=300;

/*
   Computes lnL in the case that the densities in velocity space are oomputed by distributing
   realrather than mock stars - i.e. the largest lnL that's logically possible
*/
double perfect_dyn(std::vector<df::PtrDistributionFunction>& dfs, cell& Cell,
		   grid3& grd3){
	for(size_t i=0; i<Cell.reals.size(); i++){
		grd3.CIC(1, Cell.reals[i].VR, Cell.reals[i].Vz, -Cell.reals[i].Vphi);
	}
	float total=Cell.reals.size();
	grd3.divide(total);
	double lnL=0;
	for(size_t i=0;i<Cell.reals.size();i++)
		lnL+=log(fmax(1e-20,grd3.dens(Cell.reals[i].VR,Cell.reals[i].Vz,-Cell.reals[i].Vphi)));
	return lnL;
}

/*
   Computes lnL in the case that the densities in chemical space are oomputed by distributing
   realrather than mock stars - i.e. the largest lnL that's logically possible
*/
double perfect_chem(std::vector<cell>& cells, std::vector<double> weights, int nmock, int nreal,
		std::vector<df::PtrDistributionFunction>& dfs,
		chem_model* M, grid4& grd4){
	const int NPROC=8;
	std::vector<grid4> grds(NPROC);
	for(int I=0; I<NPROC; I++) grds[I].setup(grd4);
#pragma omp parallel for schedule(dynamic)
	for(int I=0; I<NPROC; I++){
		for(size_t i=I; i<cells.size(); i+=NPROC){ 
			for(size_t j=0; j<cells[i].reals.size(); j++){// run over reals
				std::pair<double,double> chm = std::make_pair(cells[i].reals[j].FeH,cells[i].reals[j].MgFe);
				double Jphi=-cells[i].reals[j].R*cells[i].reals[j].Vphi;
				grds[I].CIC(1,chm.first,chm.second,Jphi,cells[i].reals[j].Jz);
			}
		}
	}
	for(int I=0; I<NPROC; I++) grd4.add(grds[I]);
	grd4.divide((float)nreal);
	double lnL=0;
	for(size_t i=0; i<cells.size(); i++){
		for(size_t j=0;j<cells[i].reals.size(); j++){
			double Jphi=-cells[i].reals[j].R*cells[i].reals[j].Vphi;
			lnL+=log(fmax(1e-20,grd4.dens(cells[i].reals[j].FeH,cells[i].reals[j].MgFe,
				Jphi,cells[i].reals[j].Jz)));
		}
	}
	return lnL;
}

/*
   Computes the lnL from velocity-space densities densities
*/
double val_dyn(std::vector<df::PtrDistributionFunction>& dfs, cell& Cell,
	       grid3& grd3){
	float total = 0;
	for(size_t i=0; i<Cell.mocks.size(); i++){
		int cpt=Cell.mocks[i].cpt;
		// f0 was computed with true actions, so here we must
		// also use true actions
		actions::Actions J(Cell.mocks[i].Jt.Jr*from_kpckms,Cell.mocks[i].Jt.Jz*from_kpckms,Cell.mocks[i].Jt.Jphi*from_kpckms);
		double frat=dfs[cpt]->value(J)/Cell.mocks[i].f0;
		total += frat;
		if(frat>100) printf("%d (%3.0f,%3.0f,%3.0f) %4.0f |",cpt,Cell.mocks[i].V.vR,
				    Cell.mocks[i].V.vz, Cell.mocks[i].V.vphi,frat);
		grd3.CIC(frat, Cell.mocks[i].V.vR, Cell.mocks[i].V.vz, Cell.mocks[i].V.vphi);
	}
	double lnL=0;
	grd3.divide(total);
	for(size_t i=0;i<Cell.reals.size();i++)
		lnL+=log(fmax(1e-20,grd3.dens(Cell.reals[i].VR,Cell.reals[i].Vz,-Cell.reals[i].Vphi)));
	return lnL;
}
/*
   Computes lnL from densities in cheical space
*/
double val_chem(std::vector<cell>& cells, std::vector<double> weights, int nmock, int nreal,
		std::vector<df::PtrDistributionFunction>& dfs,
		chem_model* M, grid4& grd4){
	const int NPROC=8;
	float total[NPROC] = {NPROC*0};
	std::vector<grid4> grds(NPROC);
	for(int I=0; I<NPROC; I++) grds[I].setup(grd4);
#pragma omp parallel for schedule(dynamic)
	for(int I=0; I<NPROC; I++){
		total[I]=0;
		for(size_t i=I; i<cells.size(); i+=NPROC){ 
			double EFeH=cells[i].EFeH, EMgFe=cells[i].EMgFe;
			for(size_t j=0; j<cells[i].mocks.size(); j++){// run over mocks
				std::pair<double,double> chm = M->chem(cells[i].mocks[j],EFeH,EMgFe);
				int cpt=cells[i].mocks[j].cpt;
				// f0 was computed with true actions, so here we must
				// also use true actions
				actions::Actions J(cells[i].mocks[j].Jt.Jr*from_kpckms,
					cells[i].mocks[j].Jt.Jz*from_kpckms,
					cells[i].mocks[j].Jt.Jphi*from_kpckms);
				double frat=dfs[cpt]->value(J)/cells[i].mocks[j].f0;
				total[I] += frat;
				grds[I].CIC(frat,chm.first,chm.second,cells[i].mocks[j].J.Jphi,cells[i].mocks[j].J.Jz);
			}
		}
	}
	float Total=0;
	for(int I=0; I<NPROC; I++){
		grd4.add(grds[I]); Total += total[I];
	}
	grd4.divide(Total);
	double lnL=0;
	for(size_t i=0; i<cells.size(); i++){
		for(size_t j=0;j<cells[i].reals.size(); j++){
			double Jphi=-cells[i].reals[j].R*cells[i].reals[j].Vphi;
			lnL+=log(fmax(1e-20,grd4.dens(cells[i].reals[j].FeH,cells[i].reals[j].MgFe,
				Jphi,cells[i].reals[j].Jz)));
		}
	}
	grd4.divide(1/(float)nreal);//change norm if plotting
	return lnL + M->prior();
}

/*
   Computes total lnL
*/
void Lhd::eval(const double vars[], double values[]) const{
// Update chemistry
	M->change_params(vars);
	M->unload_params();
//Update DFs
	int n0 = nvars_chem;
	if(dyn){
		for(size_t i=0; i<Kvs->size(); i++){//update KeyValueMaps
			if(i!=4) update((*Kvs)[i],*changing[i],n0,vars,true);
			else update((*Kvs)[i],*changing[i],n0,vars,false);
		}
		//adjust old disc to give required total mass
		double mass=0;
		for(int i=0;i<6;i++){
			if(i==2) continue;
			mass+=(*Kvs)[i].getDouble("mass",NAN);
		}
		(*Kvs)[2].set("mass",stellarMass-mass);
	}
	std::vector<df::PtrDistributionFunction> dfs;
	for(size_t i=0; i<Kvs->size(); i++)//create DFs
		dfs.push_back(df::createDistributionFunction((*Kvs)[i], pot.get(), NULL, extUnits));
// Compute LnLhd
	std::vector<double> vals(cells.size());
	double nc_av=0, nc_max=0, nc_min=100;
#pragma omp parallel for schedule(dynamic)
	for(int icell=0; icell<(int)cells.size(); icell++){
		int nc=(int)(.8*pow((double)cells[icell].reals.size(),.33333));
		nc_max=fmax(nc,nc_max); nc_min=fmin(nc,nc_min);
		nc_av+=nc;
		double VRmax=3*cells[icell].sigR, Vzmax=3*cells[icell].sigz;
		double Vphimin=-50, Vphimax=fabs(cells[icell].vphibar) + 2.5*cells[icell].sigphi;
		grid3 grd3(nc,nc,nc,0,-VRmax,VRmax,-Vzmax,Vzmax,Vphimin,Vphimax);
		vals[icell] = val_dyn(dfs,cells[icell],grd3);
	}
	nc_av/=cells.size();
	double val_dyn=0;
	for(size_t icell=0; icell<cells.size(); icell++)
		val_dyn += vals[icell];//Minus LnL for minimisation
	val_dyn /= nreal; //cells.size();
	if(std::isnan(val_dyn)) printf("Lkh: val_dyn is NAN\n");
	grid4 grd4(15,15,20,20,0,Femin,Femax,Mgmin,Mgmax,JphiMin,JphiMax,0,JzMax);
	double val_c = val_chem(cells,weights,nmock,nreal,dfs,M,grd4)/nreal;
	fprintf(log_file,"%f %f\n",val_dyn,val_c);
	values[0] = -(val_dyn + val_c);
	printf(".");
	if(plot){
		printf("nc_bar %5.1f, nc_min %5.1f, nc_max %5.1f\n",nc_av,nc_min,nc_max);
#pragma omp parallel for schedule(dynamic)
		for(int icell=0; icell<(int)cells.size(); icell++){
			int nc=(int)(.8*pow((double)cells[icell].reals.size(),.33333));
			double VRmax=3*cells[icell].sigR, Vzmax=3*cells[icell].sigz;
			double Vphimin=-50, Vphimax=fabs(cells[icell].vphibar) + 2.5*cells[icell].sigphi;
			grid3 grd3(nc,nc,nc,0,-VRmax,VRmax,-Vzmax,Vzmax,Vphimin,Vphimax);
			vals[icell] = perfect_dyn(dfs,cells[icell],grd3);
		}
		double lnLperf=0;
		for(int icell=0; icell<(int)cells.size(); icell++) lnLperf+=vals[icell];
		lnLperf /= nreal;
		grid4 grd4p(15,15,20,20,0,Femin,Femax,Mgmin,Mgmax,JphiMin,JphiMax,0,JzMax);
		double perf_c = perfect_chem(cells,weights,nmock,nreal,dfs,M,grd4p)/nreal;
		printf("Actual (perfect) fits: %g (%g) %g (%g)\n",val_dyn,lnLperf,val_c,perf_c);
// plot projections of 4d space ([Fe/H],[Mg/Fe],Jphi,Jz)
//		plot_chem(grd4); 
	}
}
void Lhd::plot_chem(grid4& grd4m) const{
	//First plot model <[Mg/Fe]> in JJ plane
/*
	pl1->new_plot(JphiMin,JphiMax,0,JzMax,"J\\d\\gf[kpc km s\\u-\\u1]","J\\dz[kpc km s\\u-\\u1]",1.1);
	grd4m.plot_mean('x','w',.01,.35,pl1,true,false);//average of model Mg in J plane
	pl1->relocate(JphiMax+.08*(JphiMax-JphiMin),.5*JzMax);
	pl1->setangle(90); pl1->putlabel(5,"[Mg/Fe]");
	pl1->grend();
*/
//Compute density of real stars in 4d space
	grid4 grd4o(15,15,20,20,0,Femin,Femax,Mgmin,Mgmax,JphiMin,JphiMax,0,JzMax);
	for(size_t i=0; i<cells.size(); i++){
		for(size_t j=0;j<cells[i].reals.size(); j++){
			double Jphi=-cells[i].reals[j].R*cells[i].reals[j].Vphi;
			grd4o.CIC(1,cells[i].reals[j].FeH,cells[i].reals[j].MgFe,
				  Jphi,cells[i].reals[j].Jz);
		}
	}
// Plot observed <[Mg/Fe]> in JJ plane
/*
	pl1->new_plot(JphiMin,JphiMax,0,JzMax,"J\\d\\gf[kpc km s\\u-\\u1]","J\\dz[kpc km s\\u-\\u1]",1.1);
	grd4o.plot_mean('x','w',.01,.35,pl1,true,false);//average of Mg in J plane
	pl1->relocate(JphiMax+.08*(JphiMax-JphiMin),.5*JzMax);
	pl1->setangle(90); pl1->putlabel(5,"[Mg/Fe]");
	pl1->grend();
// Plot model <[Fe/H]> in JJ plane
	pl1->new_plot(JphiMin,JphiMax,0,JzMax,"J\\d\\gf[kpc km s\\u-\\u1]","J\\dz[kpc km s\\u-\\u1]",1.1);
	grd4m.plot_meanR('w','x',-1.4,.2,pl1,true,false);//average of model Fe in J plane
	pl1->relocate(JphiMax+.08*(JphiMax-JphiMin),.5*JzMax);
	pl1->setangle(90); pl1->putlabel(5,"[Fe/H]");
	pl1->grend();
// Plot observed <[Fe/H]> in JJ plane
	pl1->new_plot(JphiMin,JphiMax,0,JzMax,"J\\d\\gf[kpc km s\\u-\\u1]","J\\dz[kpc km s\\u-\\u1]",1.1);
	grd4o.plot_meanR('w','x',-1.4,.2,pl1,true,false);//average of Fe in J plane
	pl1->relocate(JphiMax+.08*(JphiMax-JphiMin),.5*JzMax);
	pl1->setangle(90); pl1->putlabel(5,"[Fe/H]");
	pl1->grend();
//Plot c space for 30 bins in action space
	double lmin=-.2, lmax[6]={3,3,3.5,4.3,4.3,4};
	float mag=.4;
	for(int j=0;j<6;j++){
		pl1->multipane(-4,-5); 
		for(int k=0;k<5;k++){
			double _Jr=20,_Jphi=500*j,_Jz=5+35*k;
			pl1->pane(0,k);
			if(j%2==0) pl1->new_plot(Femin,Femax,Mgmin,Mgmax,"[Fe/H]","[Mg/Fe]",mag);
			else pl1->new_plot(Femin,Femax,Mgmin,Mgmax,"[Fe/H]"," ",mag);
			//plot model distribution
			grd4m.plot_plaque(lmin,lmax[j],_Jphi-200,_Jphi+200,_Jz-15,_Jz+15,pl1,false,true);
			char lab[50];
			if(k==0) sprintf(lab," (J\\d\\gf,J\\dz)=(%4.0f,%3.0f)",_Jphi,.5*(_Jz+15));
			else sprintf(lab," (J\\d\\gf,J\\dz)=(%4.0f,%3.0f)",_Jphi,_Jz);
			pl1->relocate(Femin,Mgmin+.08*(Mgmax-Mgmin)); pl1->label(lab);
			sprintf(lab,"model");
			pl1->relocate(Femin+.96*(Femax-Femin),Mgmin+.9*(Mgmax-Mgmin)); pl1->putlabel(4,lab);
			pl1->pane(1,k);
			if(j%2==0) pl1->new_plot(Femin,Femax,Mgmin,Mgmax,"[Fe/H]","[Mg/Fe]",mag);
			else pl1->new_plot(Femin,Femax,Mgmin,Mgmax,"[Fe/H]"," ",mag);
			//plot observed distribution
			grd4o.plot_plaque(lmin,lmax[j],_Jphi-200,_Jphi+200,_Jz-15,_Jz+15,pl1,true,true);
			sprintf(lab,"data");
			pl1->relocate(Femin+.96*(Femax-Femin),Mgmin+.9*(Mgmax-Mgmin));
			pl1->putlabel(4,lab);
		}
		pl1->grend();
	}
*/
}	
void Lhd::plot_grid(const double vars[], int icell) const{
//Plot projected V-distributions in a spatial cell
	int n0 = nvars_chem;
	if(dyn){
		for(int i=0;i<6;i++){//update KeyValueMaps
			if(i!=4) update((*Kvs)[i],*changing[i],n0,vars,true);
			else update((*Kvs)[i],*changing[i],n0,vars,false);
		}
	}
	std::vector<df::PtrDistributionFunction> dfs;
	for(size_t i=0; i<Kvs->size(); i++)//create DFs
		dfs.push_back(df::createDistributionFunction((*Kvs)[i], pot.get(), NULL, extUnits));
	std::vector<double> Ms;
	for(int i=0;i<6;i++) Ms.push_back(dfs[i]->totalMass());
	float total=0;
	int nc=(int)(.8*pow((double)cells[icell].reals.size(),.33333));
	double VRmax=3*cells[icell].sigR, Vzmax=3*cells[icell].sigz;
	double Vphimin=-50, Vphimax=fabs(cells[icell].vphibar)+2.5*cells[icell].sigphi;
	grid3 grd(nc,nc,nc,0,Vphimin,Vphimax,-VRmax,VRmax,-Vzmax,Vzmax);
	grid3 grdE0(nc,nc,nc,0,Vphimin,Vphimax,-VRmax,VRmax,-Vzmax,Vzmax);
	for(size_t i=0; i<cells[icell].mocks.size(); i++){
		int cpt=cells[icell].mocks[i].cpt;
		// f0 was computed with true actions, so here we must
		// also use true actions
		actions::Actions J(cells[icell].mocks[i].Jt.Jr*from_kpckms,
				   cells[icell].mocks[i].Jt.Jz*from_kpckms,
				   cells[icell].mocks[i].Jt.Jphi*from_kpckms);
		double frat=dfs[cpt]->value(J)/cells[icell].mocks[i].f0;
		total+=frat;
		grd.CIC(frat, cells[icell].mocks[i].V.vphi, cells[icell].mocks[i].V.vR,
			  cells[icell].mocks[i].V.vz);
		grdE0.CIC(frat, cells[icell].mocks[i].Vt.vphi, cells[icell].mocks[i].Vt.vR,
			  cells[icell].mocks[i].Vt.vz);
	}
	grd.divide(total/(float)cells[icell].reals.size());
	grdE0.divide(total/(float)cells[icell].reals.size());
	grid3 grdo(nc,nc,nc,0,Vphimin,Vphimax,-VRmax,VRmax,-Vzmax,Vzmax);
	for(int i=0; i<cells[icell].reals.size(); i++)
		grdo.CIC(1,-cells[icell].reals[i].Vphi,cells[icell].reals[i].VR,cells[icell].reals[i].Vz);
	char lab[100],lab_t[100],lab_s[100];
	sprintf(lab,"(%4.1f %4.2f)                            data  %zd stars",
		  cells[icell].R,cells[icell].z,cells[icell].reals.size());
	float total_m=grd.get_total(), total_o=grdo.get_total();
	sprintf(lab_t,"%5.0f  %5.0f",total_m,total_o);
	
// Now plot data and model grids
/*	pl->multipane(-2,-2);
	pl->pane(0,0);
	pl->new_plot(Vphimin,Vphimax,-Vzmax,Vzmax,"v\\d\\gf","v\\dz",.8);
	std::pair<float,float> lims = grd.plot_projection('y',pl,false);
	pl->pane(1,0);
	pl->new_plot(Vphimin,Vphimax,-Vzmax,Vzmax,"v\\d\\gf","v\\dz",.8);
	grdo.plot_projection(lims.first,lims.second,'y',pl,true);
	pl->pane(0,0);
	pl->relocate(Vphimin,-.92*Vzmax); pl->label(lab);
	pl->relocate(Vphimin+.78*(Vphimax-Vphimin),.92*Vzmax); pl->label(lab_t);
	pl->pane(0,1);
	pl->new_plot(Vphimin,Vphimax,-Vzmax,Vzmax,"v\\d\\gf","v\\dR",.8);
	lims = grd.plot_projection('z',pl,false);
	float means[3],disps[3];
	grd.get_stats(means,disps);
	sprintf(lab_s,"%5.0f %5.1f %5.1f %5.1f",
		  means[0],disps[0],disps[1],disps[2]);
	pl->relocate(Vphimin,.9*Vzmax); pl->label(lab_s);
	pl->pane(1,1);
	pl->new_plot(Vphimin,Vphimax,-Vzmax,Vzmax,"v\\d\\gf","v\\dR",.8);
	grdo.plot_projection(lims.first,lims.second,'z',pl,true);
	grdo.get_stats(means,disps);
	sprintf(lab_s,"%5.0f %5.1f %5.1f %5.1f",
		  means[0],disps[0],disps[1],disps[2]);
	pl->relocate(Vphimin,.9*Vzmax); pl->label(lab_s);
	pl->grend();

	pl2->multipane(-2,-2);
	pl2->pane(0,0);
	pl2->new_plot(Vphimin,Vphimax,-Vzmax,Vzmax,"v\\d\\gf","v\\dz",.8);
	lims = grdE0.plot_projection('y',pl2,false);
	pl2->pane(1,0);
	pl2->new_plot(Vphimin,Vphimax,-Vzmax,Vzmax,"v\\d\\gf","v\\dz",.8);
	grdo.plot_projection(lims.first,lims.second,'y',pl2,true);
	pl2->pane(0,0);
	pl2->relocate(Vphimin,-.92*Vzmax); pl2->label(lab);
	pl2->relocate(Vphimin+.78*(Vphimax-Vphimin),.92*Vzmax); pl2->label(lab_t);
	pl2->pane(0,1);
	pl2->new_plot(Vphimin,Vphimax,-Vzmax,Vzmax,"v\\d\\gf","v\\dR",.8);
	lims = grdE0.plot_projection('z',pl2,false);
	pl2->pane(1,1);
	pl2->new_plot(Vphimin,Vphimax,-Vzmax,Vzmax,"v\\d\\gf","v\\dR",.8);
	grdo.plot_projection(lims.first,lims.second,'z',pl2,true);
	pl2->grend();
*/
}

/*
 Read in star sample, either from raw APOGEE data or from pre-digested
 data for cells
*/
int read_data(const std::string dir,int ncpt,float Rs[],float zs[],int nR,int nz,
	      std::vector<double>& Ms0, std::vector<cell>& cells, potential::PtrPotential& pot,errors& errs){
	int ncell=cells.size();
	std::string fname(dir + "cells");
	std::cout << "reading "+fname+"\n";
	FILE* ifile;
	if(!fopen_s(&ifile,fname.c_str(),"r")){//file cells can beread
		errs.read(ifile,Ms0);
		errs.report("Gaia_errors");
		for(int i=0; i<ncell; i++) cells[i].read(ifile);
		fclose(ifile);
		printf("%zd cells read in\n",cells.size());
	} else {
/*
 file cells not there so read raw data from /APOGEE/APOGEE_DR17_extra.txt
 Creat files cells and /APOGEE/DR17_nos that's read by getPhi.cpp
*/
		printf("I can't open file %s so reading raw data\n",fname.c_str());
		obs::solarShifter sun(intUnits);
		actions::PtrActionFinder AF;
		std::cout << "Updating action finder..."<<std::flush;
		AF.reset(new actions::ActionFinderAxisymFudge(pot,false));
		std::cout << "done" << std::endl;
		std::vector<long_real> reals;
		int nreal = read_17(true,1,sun,intUnits,reals,errs,pot,AF);
		if(nreal==0){
			printf("nreal==0!\n"); return 0;
		}
		printf("%8.0f real stars on grids\n", pack_cells(nR,nz,Rs,zs,reals,cells));
		FILE* ofile; //write locations & pops of cells for getPhi.cpp
		fopen_s(&ofile,"/APOGEE/DR17_nos","w");
		fprintf(ofile,"%zd\n",cells.size());
		for(size_t i=0; i<cells.size(); i++)
			fprintf(ofile,"%f %f %zd\n",cells[i].R,cells[i].z,cells[i].reals.size());
		fclose(ofile);
		int nmock = read_mocks(sun,ncpt,dir,Ms0,cells,errs,pot,AF);
		if(nmock==0){
			printf("You need to run getPhi\n"); return 0;
		}
		fopen_s(&ofile,fname.c_str(),"w");
		errs.write(ofile,Ms0);
		for(int i=0; i<ncell; i++)
			cells[i].write(ofile);
		fclose(ofile);
	}
	return cells.size();
}
/*
 "relax 200 M7" will run Nelder-Mead for 200 iterations on the model
 specified in the ../M7 directory 
*/
int main(int narg,char **args){
	if(narg!=3){
		printf("you must enter a number of iterations and Mx\n"); return 0;
	}
// if cont-true we are continuing to iterate & can read a chemical model from chem_params.dat"	
	bool cont = true;
// if dyn=true we adjust the DFs at fixed chemistry & vice versa if
// dyn=false 
	bool dyn = false, chem = !dyn;
	//cpts ordered: Yd,Md,Od,Thk,SH,B
	int mod_n=2, nstep; sscanf(args[1],"%d",&nstep);
	std::string Mno(args[2]);
	std::string dir(Mno+"/");
	potential::PtrPotential pot;
	if(utils::fileExists(dir + "potential"))
		pot = potential::readPotential(dir + "potential", extUnits);
	std::cout << "Potential loaded: value at origin=-("<<
			(sqrt(-pot->value(coord::PosCyl(0,0,0))) * intUnits.to_kms) << " km/s)^2\n";
	const std::string iniFileName = dir + "m.ini";
	utils::ConfigFile ini(iniFileName);
	utils::KeyValueMap
			iniDM = ini.findSection("DF dark halo"),
	iniDFyoungDisk   = ini.findSection("DF young disk"),
	iniDFmiddleDisk  = ini.findSection("DF middle disk"),
	iniDFoldDisk     = ini.findSection("DF old disk"),
	iniDFhighADisk   = ini.findSection("DF highA disk"),
	iniDFStellarHalo = ini.findSection("DF stellar halo"),
	iniDFBulge       = ini.findSection("DF bulge");
	std::vector<utils::KeyValueMap> Kvs;
	Kvs.push_back(iniDFyoungDisk);
	Kvs.push_back(iniDFmiddleDisk);
	Kvs.push_back(iniDFoldDisk);
	Kvs.push_back(iniDFhighADisk);
	Kvs.push_back(iniDFStellarHalo);
	Kvs.push_back(iniDFBulge);

//Declare variables to adjust. d0i specifies which parameters of the
//young disc are adjusted, d1i ditto the middle disc, etc. xstep0
//specifies the largest permissable steps to be taken in the parameters
	const int nd0=9, nd1=9, nd2=8, nd3=8, nd4=5, nd5=8, nvars_dyn = nd0+nd1+nd2+nd3+nd4+nd5; 
	int d0i[nd0]={0,1,2,3,4,5,8,9,10};
	double xstep0[nd0]={.1,.1,.1,.1,.1,.1,.1,.1,.1};//extUnits
	int d1i[nd1]={0,1,2,3,4,5,8,9,10};
	double xstep1[nd1]={.1,.1,.1,.1,.1,.1,.1,.1,.1};
	int d2i[nd2]={1,2,3,4,5,8,9,10};
	double xstep2[nd2]={.1,.1,.1,.1,.1,.1,.1,.1};
	int d3i[nd3]={0,1,2,3,6,7,8,9};//High A disc
	double xstep3[nd3]={.1,.1,.1,.1,.1,.1,.1,.1};
	int d4i[nd4]={0,1,2,3,8};//S Halo
	double xstep4[nd4]={.1,.1,.1,.1,.1};//extUnits
	int d5i[nd5]={0,1,2,3,6,7,8,9};//Bulge
	double xstep5[nd5]={.1,.1,.1,.1,.1,.1,.1,.1};//extUnits
	std::vector<int> d0 (d0i,d0i+nd0);
	std::vector<int> d1 (d1i,d1i+nd1);
	std::vector<int> d2 (d2i,d2i+nd2);
	std::vector<int> d3 (d3i,d3i+nd3);
	std::vector<int> d4 (d4i,d4i+nd4);
	std::vector<int> d5 (d5i,d5i+nd5);
	const int nchange=6;
	std::vector<int>* changing[nchange] = {&d0,&d1,&d2,&d3,&d4,&d5};

	std::vector<double> xstep_dyn;
	for(size_t i=0; i<d0.size(); i++) xstep_dyn.push_back(xstep0[i]);
	for(size_t i=0; i<d1.size(); i++) xstep_dyn.push_back(xstep1[i]);
	for(size_t i=0; i<d2.size(); i++) xstep_dyn.push_back(xstep2[i]);
	for(size_t i=0; i<d3.size(); i++) xstep_dyn.push_back(xstep3[i]);
	for(size_t i=0; i<d4.size(); i++) xstep_dyn.push_back(xstep4[i]);
	for(size_t i=0; i<d5.size(); i++) xstep_dyn.push_back(xstep5[i]);
//Pick out starting values of variable parameters
	double* xi_dyn = new double[nvars_dyn];
	double stellarMass=0;
	int n0=0;
	FILE* ofile;
// We always start from the DFs specified in m.ini.
// If we are varying the parameters, we copy them from the KeyValueMaps into xi_dyn
	if(dyn){
		for(int i=0; i<6; i++){
			if(i!=4) copyin(Kvs[i],*changing[i],n0,xi_dyn,true);
			else copyin(Kvs[i],*changing[i],n0,xi_dyn,false);
		}
		for(int i=0;i<6;i++) stellarMass+=Kvs[i].getDouble("mass",NAN);
		print_dyn(changing,nchange,xi_dyn);
	}
	const float Femin=-2.5, Femax=0.7, Mgmin=-0.15, Mgmax=0.5,
	JphiMax=3000, JphiMin=-300, JzMax=200;
	const int nx=40, ny=30;
// Specify grid of spatial cells by nR + nz boundary lines 
	const int nR=13,nz=7,ncell=(nR-1)*(nz-1);
	float Rs[nR]={.5,1.5,3,4,5,6,7,8,9,10,11.5,13,14.5}, zs[nz]={0,0.3,.7,1,1.5,2,3};
	const int ncpt=6;// we are using 6 stellar DFs
	std::vector<cell> cells(ncell);
	errors errs(26,25);// Errors are assumed to depend on distance: 26 dist bins out to 25 kpc
	std::vector<double> Ms0;
	if(0==read_data(dir,ncpt,Rs,zs,nR,nz,Ms0,cells,pot,errs)) return 1;
	FILE* hfile;
	std::string fname = dir + "histograms";
	fopen_s(&hfile,fname.c_str(),"w");
	for(int i=0; i<cells.size(); i++){
		float VRmax=3*cells[i].sigR, Vzmax=3*cells[i].sigz, Vphimax=fabs(cells[i].vphibar)+3*cells[i].sigphi;
		V_hists(hfile,cells[i],VRmax,Vzmax,-50,Vphimax,50);
	}
	fclose(hfile);
	fopen_s(&hfile,fname.c_str(),"r");
	
/* plot 1d histograms of velocities for human consumption
	pl.multipane(-9,-12); float mag=.3;
	for(int i=0; i<cells.size(); i++){
		plot_V_hists(hfile,(i/12)%2,i%12,mag,&pl);
		if(i%12==11) pl.grend();
	}
	fclose(hfile);
*/
	std::vector<grid2> grdso(cells.size());
	for(size_t icell=0; icell<cells.size(); icell++)
		grdso[icell].setup(nx,ny,2,Femin,Femax,Mgmin,Mgmax);
	std::vector<std::string> cpts;
	cpts.push_back("young disk"); cpts.push_back("middle disk"); cpts.push_back("old disk");
	cpts.push_back("high-$\\alpha$ disk");
	cpts.push_back("stellar halo"); cpts.push_back("bulge");
	cpts.push_back("Thick d., low $\\alpha$");
	//cpts ordered: Yd,Md,Od,Thk,SH,B
	const int nvars_chem=70;

// map specifies which parameters of the chemical model are varied
	int map[nvars_chem]={
		0, 1, 2, 3, 4, 5, //FeH
		6, 7, 8, 9, 10,11,//MgFe
		12,13,14,15,16,17,//sigx
		18,19,20,21,22,23,//sigy
		24,25,26,27,      //th
		30,31,32,33,34,35, //grads thin discs
		36,37,38,39,40,41,
		42,43,44,45,46,47,
		48,49,50,51,52,53, //grads high-A
		54,55,56,57,58,59,//stellar halo
		60,61,62,63,64,65  //grads bulge
	};
// xstep_chem sets upper limits on steps of varied paramtes
	double xstep_chem[nvars_chem]={
		.05,.05,.02,.02,.02,.02,//Fbar
		.02,.02,.02,.02,.02,.02,//Mgbar
		.02,.02,.03,.03,.05,.05,//sigx
		.02,.02,.02,.02,.02,.02,//sigy
		.5,.5,.5,.5,//th
		.005,.005,.005,.002,.002,.002,//grads young disc
		.005,.005,.005,.002,.002,.002,//grads neium disc
		.005,.005,.005,.002,.002,.002,//grads old disc
		.005,.005,.005,.002,.002,.002,//grads high-A disc
		.005,.005,.005,.002,.002,.002,//grads stellar halo
		.005,.005,.005,.002,.002,.002 //grads bulge
	};
	int nvars_chem1 = chem? nvars_chem : 0;
	chem_model M(map,nvars_chem1);
// Load default values for the chemical parameters. If cont-true these
// will be overwritten by values in chem_params.dat
	M.load_Fbar(-.06,-.1,-.4,-0.5,-1.1,.4);//0-5
	M.load_Mgbar(.036,.04,.1,.33,.3,.03);//6-11
	M.load_sigx(.1,.12,.16,.2,.15,.3);//12-17
	M.load_sigy(.035,.03,.04,.05,.05,.04);//18-23
	M.load_th(-7,-7,-8,-8,-3,-6);//24-29
	float C0[6]={0,0,-.07/240,  0,0,0},//young 30-35
	C1[6]={0,0,-.07/240.,  0,0,0},//middle 36-41
	C2[6]={0,0,-.07/240.,  .0001,.0001,0},//old 42-47
	C3[6]={0,-.01,0,  0,0,0},//high-A  48-53
	C4[6]={0,0,-.07/240,  0,0,0},//SH 54-59
	C5[6]={0,0,0,  0,0,0};//bulge  60-65
	M.load_grad(0,C0);
	M.load_grad(1,C1);
	M.load_grad(2,C2);
	M.load_grad(3,C3);
	M.load_grad(4,C4);
	M.load_grad(5,C5);
	M.load_params();
// Write the initial chemical parameters on a TeX file
	fname = dir + "chem_table0.tex";
	tab_chemParams(fname,cpts,M);
	if(cont){//read in latest values of chemical parameters
		fname = dir + "chem_params.dat";
		if(!M.read(fname))	return 0;
		M.load_params();
	}
	int nvars_dyn1 = dyn? nvars_dyn : 0;
	int nvars= nvars_dyn1 + nvars_chem1;
// Initialise the likelihood object
	Lhd F(&Kvs, stellarMass, changing, dyn, &M, nvars_chem1, nvars, pot, grdso, nx, ny, ncpt,
	      Femin, Femax, Mgmin, Mgmax, JphiMin, JphiMax, JzMax, Ms0, cells);
	double Lstart,Lstop,*xi, *xstep, *result;
	xi = new double[nvars]; xstep = new double[nvars]; result = new double[nvars];
	if(chem){
		for(int i=0; i<nvars_chem1; i++) xstep[i]=xstep_chem[i];
		M.get_xi(xi);
	}
	if(dyn){
		for(int i=0; i<nvars_dyn;  i++) xstep[nvars_chem1+i]=xstep_dyn[i];
		for(int i=0; i<nvars_dyn1; i++)    xi[nvars_chem1+i]=xi_dyn[i];
	}
	printf("Calling F with %d chem + %d dyn = %d parameters\n",
	       nvars_chem1,nvars_dyn1,nvars);
	F.eval(xi,&Lstart); Lstart*=-1;//Starting likelihood
	printf("Lstart: %f\n",Lstart);
	if(isnan(Lstart)) return 0;
// Cal NelderMead
	int niter = nstep>0? math::findMinNdim(F,xi,xstep,.0000001,nstep,result) : 0;
	printf("%d iterations done\n",niter);
	if(niter == 0)	for(int i=0; i<nvars; i++) result[i]=xi[i];
// Compute final likelihood and produce diagnostic plots
//	F.toggle_plot();
	F.eval(result,&Lstop); Lstop*=-1;
// Plot a selection of velocity spces
	const int n1=1,n2=3,n3=5,n4=7,n5=10,ngr=5*5;
	int ngrd[ngr]={n1,n1+12,n1+24,n1+36,n1+48,n2,n2+12,n2+24,n2+36,n2+48,
	n3,n3+12,n3+24,n3+36,n3+48,n4,n4+12,n4+24,n4+36,n4+48,n5,n5+12,n5+24,n5+36,n5+48};
	
	for(int i=0;i<ngr;i++)
		F.plot_grid(result,ngrd[i]);
	if(nstep>0){
		printf("L & increase in log L: %f %f\n",Lstop,Lstop-Lstart);
		fname = dir +  "chem_params.dat";
		M.write(fname.c_str(),result);
		fname = dir +  "dyn_params.dat";
		std::ofstream ostrm(fname);
		n0=nvars_chem1;
		for(int i=0; i<Kvs.size(); i++){
		    if(dyn){
			    if(i!=4) update((Kvs)[i],*changing[i],n0,result,true);
			    else update((Kvs)[i],*changing[i],n0,result,false);
		    }
		    ostrm << "[DF " << cpts[i] << "]\n"; 
		    df::PtrDistributionFunction DF = df::createDistributionFunction(Kvs[i], pot.get(), NULL, extUnits);
		    DF->write_params(ostrm,intUnits);
		}
		ostrm.close();
		fname = dir + "DFparams.tex";
		std::ofstream tabstrm(fname);
		tabstrm << "Component \t&$M$\t&$J_{\\phi0}$\t&$J_{r0}$\t&$J_{z0}$\t&$J_{\\rm int}$\t&$D_{\\rm int}$\t&" <<
				"$J_{\\rm ext}$\t&$D_{\\rm ext}$\t&$p_r$\t&$p_z$\t&$J_{\\rm v0}$\t&$J_{\\rm d0}$\\cr\n \\hline \n";
		for(int i=0; i<Kvs.size(); i++){
			if(i==0) tab_params(tabstrm,cpts[i],Kvs[i],true,false);
			else if(i<3) tab_params(tabstrm,cpts[i],Kvs[i],false,false);
			else if(i==3 || i==5) tab_params(tabstrm,cpts[i],Kvs[i],false,true);
		}
		tabstrm.close();
		fname = dir + "DFhalo.tex";
		std::ofstream tabstrm2(fname);
		tabstrm2 << "Component \t&$M$\t& $J_0$\t& $J_{\\rm core}$\t& $J_{\\rm cutoff}$\t& $\\alpha_{\\rm in}$\t& $\\alpha_{\\rm out}$\t& $F_{\\rm in}$\t& $F_{\\rm out}$\\cr\n\\hline\n";
		tab_params(tabstrm2,std::string("dark halo"),iniDM);
		tab_params(tabstrm2,cpts[4],Kvs[4]);
	}
	fname = dir + "chem_table.tex";
	tab_chemParams(fname,cpts,M);
	return 1;
}
