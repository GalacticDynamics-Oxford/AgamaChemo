#include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "units.h"
#include "utils.h"
#include "coord.h"
#include "obs.h"
#include "math_base.h"
#include "math_fit.h"
#include "math_random.h"
#include "potential_composite.h"
#include "potential_factory.h"
#include "potential_multipole.h"
#include "potential_utils.h"
#include "actions_staeckel.h"
#include "df_base.h"
#include "df_factory.h"
#include "grid.h"
#include "utils_config.h"


// define internal unit system - arbitrary numbers here! the result should not depend on their choice
const units::InternalUnits intUnits(2.7183*units::Kpc, 3.1416*units::Myr);

// define external unit system describing the data (including the parameters in INI file)
const units::ExternalUnits extUnits(intUnits, 1.*units::Kpc, 1.*units::kms, 1.*units::Msun);

double from_kpckms=intUnits.from_Kpc*intUnits.from_kms;

//Note:: KeyValueMaps in extUnits, pars in intUnits. So val[] in
//extUnits

void tab_params(std::ofstream& strm,const std::string& cpt,const utils::KeyValueMap& Kvm,const bool taper, const bool cut){
	double mass=Kvm.getDouble("mass")/1e10;
	double Jphi0=Kvm.getDouble("Jphi0");
	double Jr0=Kvm.getDouble("Jr0");
	double Jz0=Kvm.getDouble("Jz0");
	double Jtaper=Kvm.getDouble("Jtaper");
	double Jtrans=Kvm.getDouble("Jtrans");
	double Jcut=Kvm.getDouble("Jcut");
	double Delta=Kvm.getDouble("Delta");
	double pr=Kvm.getDouble("pr");
	double pz=Kvm.getDouble("pz");
	double addJvel=Kvm.getDouble("addJvel");
	double addJden=Kvm.getDouble("addJden");
	strm << cpt << "\t&";
	strm << std::setprecision(3) << "$" << mass << "$\t&";
	strm << std::setprecision(4) << "$" << Jphi0 << "$\t&";
	strm << std::setprecision(4) << "$" << Jr0 << "$\t&";
	strm << std::setprecision(4) << "$" << Jz0 << "$\t&";
	if(taper){
		strm << std::setprecision(4) << "$" << Jtaper << "$\t&";
		strm << std::setprecision(4) << "$" << Jtrans << "$\t&";
	} else strm << "--\t&--\t&";
	if(cut){
		strm << std::setprecision(4) << "$" << Jcut << "$\t&";
		strm << std::setprecision(4) << "$" << Delta << "$\t&";
	} else strm << "--\t&--\t&";		
	strm << std::setprecision(2) << "$" << pr << "$\t&";
	strm << std::setprecision(2) << "$" << pz << "$\t&";
	strm << std::setprecision(4) << "$" << addJden << "$\t&";
	strm << std::setprecision(4) << "$" << addJvel << "$\\cr\n";
}
void tab_params(std::ofstream& strm,std::string& cpt,const utils::KeyValueMap& Kvm){
	double mass=Kvm.getDouble("mass")/1e10;
	double J0=Kvm.getDouble("J0");
	double Jcore=Kvm.getDouble("Jcore");
	double Jcutoff=Kvm.getDouble("Jcutoff");
	double slopeIn=Kvm.getDouble("slopeIn");
	double slopeOut=Kvm.getDouble("slopeOut");
	double Fin=Kvm.getDouble("Fin");
	double Fout=Kvm.getDouble("Fout");
	strm << cpt << "\t&";
	strm << std::setprecision(4) << "$" << mass << "$\t&";
	strm << std::setprecision(4) << "$" << J0 << "$\t&";
	strm << std::setprecision(4) << "$" << Jcore << "$\t&";
	strm << std::setprecision(4) << "$" << Jcutoff << "$\t&";
	strm << std::setprecision(2) << "$" << slopeIn << "$\t&";
	strm << std::setprecision(2) << "$" << slopeOut << "$\t&";
	strm << std::setprecision(2) << "$" << Fin << "$\t&";
	strm << std::setprecision(2) << "$" << Fout << "$\\cr\n";
}
void change(utils::KeyValueMap& Kv,int n,const double val){
	switch(n){//12 parameters of taperExp disc (all intUnits)
		case 0:
			Kv.set("mass",exp(val+18)); break;
		case 1:
			Kv.set("Jr0",exp(val)); break;
		case 2:
			Kv.set("Jz0",exp(val)); break;
		case 3:
			Kv.set("Jphi0",exp(val)*10); break;
		case 4:
			Kv.set("Jtaper",exp(val)*10); break;
		case 5:
			Kv.set("Jtrans",exp(val)*5); break;
		case 6:
			Kv.set("Jcut",exp(val)*10); break;
		case 7:
			Kv.set("Delta",exp(val)*5); break;
		case 8:
			Kv.set("pr",val); break;
		case 9:
			Kv.set("pz",val); break;
		case 10:
			Kv.set("addJden",exp(val)*5); break;
		case 11:
			Kv.set("addJvel",exp(val)); break;
	}
}
void extract(utils::KeyValueMap& Kv,int n,double& val){
	switch(n){//12 parameters of taperExp disc (all intUnits)
		case 0:
			val = log(Kv.getDouble("mass",NAN))-18; break;
		case 1:
			val = log(Kv.getDouble("Jr0",NAN)); break;
		case 2:
			val = log(Kv.getDouble("Jz0",NAN)); break;
		case 3:
			val = log(Kv.getDouble("Jphi0",NAN)/10); break;
		case 4:
			val = log(Kv.getDouble("Jtaper",NAN)/10); break;
		case 5:
			val = log(Kv.getDouble("Jtrans",NAN)/5); break;
		case 6:
			val = log(Kv.getDouble("Jcut",NAN)/10); break;
		case 7:
			val = log(Kv.getDouble("Delta",NAN)/5); break;
		case 8:
			val = Kv.getDouble("pr",NAN); break;
		case 9:
			val = Kv.getDouble("pz",NAN); break;
		case 10:
			val = log(Kv.getDouble("addJden",NAN)/5); break;
		case 11:
			val = log(Kv.getDouble("addJvel",NAN)); break;
	}
}
void changeH(utils::KeyValueMap& Kv,int n,const double val){
	switch(n){//9 parameters of spheroid DF
		case 0:
			Kv.set("mass",exp(val+18)); break;
		case 1:
			Kv.set("J0",exp(val)*5); break;
		case 2:
			Kv.set("slopeIn",exp(val)); break;
		case 3:
			Kv.set("slopeOut",exp(val)); break;
		case 4:
			Kv.set("Fin",val); break;
		case 5:
			Kv.set("Fout",val); break;
		case 6:
			Kv.set("alpha",val); break;
		case 7:
			Kv.set("beta",val); break;
		case 8:
			Kv.set("Jcore",exp(val)); break;
	}
}
void extractH(utils::KeyValueMap& Kv,int n,double& val){
	switch(n){//9 parameters of spheroid DF
		case 0:{
			val = log(Kv.getDouble("mass",NAN))-18; ;break;}
		case 1:
			val = log(Kv.getDouble("J0",NAN)/5); break;
		case 2:
			val = log(Kv.getDouble("slopeIn",NAN)); break;
		case 3:
			val = log(Kv.getDouble("slopeOut",NAN)); break;
		case 4:
			val = Kv.getDouble("Fin",val); break;
		case 5:
			val = Kv.getDouble("Fout",val); break;
		case 6:
			val = Kv.getDouble("alpha",val); break;
		case 7:
			val = Kv.getDouble("beta",val); break;
		case 8:
			val = log(Kv.getDouble("Jcore",NAN)); break;
	}
}
void reportD(utils::KeyValueMap& Kv){
	printf("M: %6.2g, Jr0:%5.1f, Jz0:%5.1f, Jphi0: %5.1f, ",Kv.getDouble("mass",NAN),
	       Kv.getDouble("Jr0",NAN), Kv.getDouble("Jz0",NAN), Kv.getDouble("Jphi0",NAN));
	printf("pr: %5.2f, pz: %5.2f, addJden:%5.1f, addJvel:%5.1f\n", Kv.getDouble("pr",NAN)
	       , Kv.getDouble("pz",NAN), Kv.getDouble("addJden",NAN), Kv.getDouble("addJvel",NAN));
}
void reportH(utils::KeyValueMap& Kv){
	printf("M: %6.2g, J0: %5.1f, slopeIn: %5.3f, slopeOut: %5.3f, ",Kv.getDouble("mass",NAN),
	       Kv.getDouble("J0",NAN),Kv.getDouble("slopeIn",NAN),Kv.getDouble("slopeOut",NAN));
	printf("Fin: %5.3f, Fout: %5.3f, alpha: %5.3f, beta: %5.3f, Jc:%5.1f\n", Kv.getDouble("Fin",NAN),
	       Kv.getDouble("Fout",NAN), Kv.getDouble("alpha",NAN), Kv.getDouble("beta",NAN), Kv.getDouble("Jcore",NAN));
}
void update(utils::KeyValueMap& Kv,std::vector<int>& d,int& n0,const double val[],bool disc){
	int n = d.size();
	for(int i=0; i<n; i++){
		if(disc) change(Kv,d[i],val[n0+i]);
		else changeH(Kv,d[i],val[n0+i]);
	}
	n0+=n;
}
void copyin(utils::KeyValueMap& Kv,std::vector<int>& d,int& n0,double val[],bool disc){
	int n = d.size();
	for(int i=0; i<n; i++){
		if(disc) extract(Kv,d[i],val[n0+i]);
		else extractH(Kv,d[i],val[n0+i]);
	}
	n0+=n;
}
int print_dyn(std::vector<int>* changing[],const int nchange,double vars[]){
	printf("DF parameters that we're changing by cpt:\n");
	int n=0;
	for(int i=0;i<nchange;i++){//run over discs
		for(size_t j=0;j<changing[i]->size();j++)  printf("(%d %g)",(*changing[i])[j],vars[n+j]);
		printf("\n");
		n+=changing[i]->size();
	}
	return n;
}
class errors{
	public:
		int n;//# of distance bins
		int* Ns;//# of stars in bin
		float* s;//distance of bin
		float *EFeH, *EMgFe, *Epmra, *Epmdec, *EVlos, *Es;
		errors(int _n,float smax) : n(_n){
			s=new float[n]; Ns=new int[n];
			Epmra=new float[n]; Epmdec=new float[n];
			EVlos=new float[n]; Es=new float[n];
			EFeH=new float[n]; EMgFe=new float[n];
			float ds=smax/(float)(n-1);
			for(int i=0;i<n;i++){
				s[i]=i*ds; Ns[i]=0;
				Epmra[i]=0; Epmdec[i]=0;
				Es[i]=0; EVlos[i]=0;
				EFeH[i]=0; EMgFe[i]=0;
			}
		}
		~errors(void){
			delete[] s; delete[] Ns;
			delete[] Epmra; delete[] Epmdec;
			delete[] EVlos; delete[] Es;
			delete[] EFeH; delete[] EMgFe;
		}
		void add(float _s,float _Epmra,float _Epmdec,float _Es,float _EVlos,float _EFeH,float _EMgFe){
			int k=0;
			while(_s>s[k] && k<=n) k++;
			k--;
			if(k==n) return;
			Ns[k]++;
			double f1=1./(double)Ns[k],f0=1-f1;
			Epmra[k]=f0*Epmra[k]+f1*_Epmra;
			Epmdec[k]=f0*Epmdec[k]+f1*_Epmdec;
			Es[k]=f0*Es[k]+f1*_Es;
			EVlos[k]=f0*EVlos[k]+f1*_EVlos;
			EFeH[k]=f0*EFeH[k]+f1*_EFeH;
			EMgFe[k]=f0*EMgFe[k]+f1*_EMgFe;
		}
		void get(float _s,float& _Epmra,float& _Epmdec,float& _Es,
			 float& _EVlos,float& _EFeH,float& _EMgFe){
			int k=0;
			while(_s>s[k] && k<n) k++;
			if(k==n) return;
			k--;
			_Epmra=Epmra[k]; _Epmdec=Epmdec[k]; _Es=Es[k];
			_EVlos=EVlos[k]; _EFeH=EFeH[k]; _EMgFe=EMgFe[k];
		}			
		bool read(FILE* ifile,std::vector<double>& Ms0){
			double M; int nc;
			fscanf(ifile,"%d",&nc);
			for(int i=0; i<nc; i++){
				fscanf(ifile,"%lg",&M); Ms0.push_back(M);
			}
			bool ok;
			int m;
			fscanf(ifile,"%d\n",&m);
			ok = m==n;
			for(int i=0; i<n; i++){
				ok=ok && 8==fscanf(ifile,"%f %d %f %f %f %f %f %f\n",
					  s+i,Ns+i,Epmra+i,Epmdec+i,Es+i,EVlos+i,EFeH+i,EMgFe+i);
			}
			if(!ok) printf("Error in errors.read\n");
			return ok;
		}
		void write(FILE* ofile,std::vector<double> Ms0){
			fprintf(ofile,"%zd\n",Ms0.size());
			for(size_t i=0; i<Ms0.size(); i++) fprintf(ofile,"%g ",Ms0[i]);
			fprintf(ofile,"\n");
			fprintf(ofile,"%d\n",n);
			for(int i=0; i<n; i++){
				fprintf(ofile,"%f %d %f %f %f %f %f %f\n",
					  s[i],Ns[i],Epmra[i],Epmdec[i],Es[i],EVlos[i],EFeH[i],EMgFe[i]);
			}
		}
		void report(std::string fname){
			FILE* ofile; fopen_s(&ofile,fname.c_str(),"w");
			fprintf(ofile,"   s   N    Es     Epmra  Epmdec EVlos  EFeH  EMgFe\n");
			for(int i=0;i<n;i++) fprintf(ofile,"%4.0f %5d %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",
				s[i],Ns[i],Es[i],Epmra[i],Epmdec[i],EVlos[i],EFeH[i],EMgFe[i]);
			fprintf(ofile,"\n");
			fclose(ofile);
		}
};
struct  star{
	char id[30];
	double ra,dec;
	float Teff, logg;//StarHorse
	float Vlos, Vlos_err, DVlos;
	float FeH, FeH_err, MgFe, MgFe_err;//from APOCAST
	float FeH_nn, FeH_nn_err, MgFe_nn, MgFe_nn_err;//from Bovy
	float pi, pi_err;
	float pm_ra, pm_ra_err, pm_dec, pm_dec_err;
	float s_nn, s_nn_err, s_sh, s_sh_err;//StarHorse
};
class long_real{
	public:
		float R, z, phi;
		float FeH, MgFe;
		float VR, Vz, Vphi;
		float Jr, Jz;
		float Epmra,Epmdec,Es,EVlos;
		float EFeH,EMgFe;
		long_real(float _FeH, float _MgFe, float _R, float _z, float _phi,
			  float _VR, float _Vz, float _Vphi, float _Jr, float _Jz,
			  float _Epmra, float _Epmdec, float _Es, float _EVlos, float _EFeH, float _EMgFe)  :
		    FeH(_FeH), MgFe(_MgFe), R(_R), z(_z), phi(_phi), VR(_VR), Vz(_Vz), Vphi(_Vphi), Jr(_Jr), Jz(_Jz),
		    Epmra(_Epmra), Epmdec(_Epmdec), Es(_Es), EVlos(_EVlos),
		    EFeH(_EFeH), EMgFe(_EMgFe) {};
};
class real{
	public:
		float R, z, phi;
		float FeH, MgFe;
		float VR, Vz, Vphi;
		float Jr, Jz;
		real(float _FeH, float _MgFe, float _R, float _z, float _phi, float _VR, float _Vz, float _Vphi, float _Jr, float _Jz) :
		    FeH(_FeH), MgFe(_MgFe), R(_R), z(_z), phi(_phi), VR(_VR), Vz(_Vz), Vphi(_Vphi), Jr(_Jr), Jz(_Jz) {};
		void write(FILE* ofile){//we always write Js
			float w[10]={FeH,MgFe,R,z,phi,VR,Vz,Vphi,Jr,Jz};
			compress(ofile,w,10);
		};
};
real read_real(FILE* ifile){
	float w[10]; if(!get(ifile,w,10))
		printf("Error reading real\n%f %f %f %f %f\n",w[0],w[1],w[2],w[3],w[4]);
	return real(w[0],w[1],w[2],w[3],w[4],w[5],w[6],w[7],w[8],w[9]);
};
class mock{
	public:
		int cpt;
		actions::Actions Jt,J;
		coord::VelCyl V, Vt;
//		double Jtr,Jtz,Jtphi;//stored in extUnits
//		double VR,Vz,Vphi;//stored in extUnits
		double f0;
		mock(int _i,actions::Actions _Jt,actions::Actions _J,coord::VelCyl _V,
		     coord::VelCyl _Vt, double _f0):
		    cpt(_i), Jt(_Jt), J(_J), V(_V), Vt(_Vt), f0(_f0) {};
		void write(FILE* ofile){
			fprintf(ofile,"%d\n",cpt);
			double w[13]={Jt.Jr,Jt.Jz,Jt.Jphi,J.Jr,J.Jz,J.Jphi,V.vR,V.vz,V.vphi,Vt.vR,Vt.vz,Vt.vphi,f0};
			compress(ofile,w,13);
		};
};
mock read_mock(FILE* ifile){
	int _cpt;
	if(1!=fscanf(ifile,"%d",&_cpt)){
		printf("Err reading cpt: %d\n",_cpt); _cpt=-1;
		char line[200];
		fgets(line,200,ifile); printf("%s\n",line);
		fgets(line,200,ifile); printf("%s\n",line);
	}
	double w[13];
	if(!get(ifile,w,13))
		printf("Err reading mock: %f %f %f %f %f\n",w[0],w[1],w[2],w[3],w[4]);
	return mock(_cpt,actions::Actions(w[0],w[1],w[2]),actions::Actions(w[3],w[4],w[5]),
		    coord::VelCyl(w[6],w[7],w[8]),coord::VelCyl(w[9],w[10],w[11]),w[12]);
};
class cell {
	public:
		double R,z;
		double Epmra,Epmdec,Es,EVlos,EFeH,EMgFe;
		double sigR,sigz,sigphi,vphibar;
		std::vector<real> reals;//stores (Fe,a)
		std::vector<double> phis;
		std::vector<mock> mocks;//store (Jphi,Jz)
		cell(void){
			R=0; z=0; sigR=0; sigz=0; sigphi=0; vphibar=0;
			Epmra=0; Epmdec=0; Es=0; EVlos=0; EFeH=0; EMgFe=0;
		}
		void load_real(long_real& st){
			reals.push_back(real(st.FeH,st.MgFe,st.R,st.z,st.phi,st.VR,st.Vz,st.Vphi,st.Jr,st.Jz));
		}
		void load_mock(actions::Actions Jt, actions::Actions J, coord::VelCyl V, coord::VelCyl Vt,
			       double f0, int ic){
			mocks.push_back(mock(ic,Jt,J,V,Vt,f0));
		}
		void read(FILE* ifile){
			int nr,nm;
			if(2!=fscanf(ifile,"%d %d",&nr,&nm))
				printf("Err reading nr %d %d\n",nr,nm);
			//else printf("%d %d\n",nr,nm);
			if(6!=fscanf(ifile,"%lf %lf %lf %lf %lf %lf",&R,&z,&sigR,&sigz,&sigphi,&vphibar))
				printf("Err reading Rz..: %f %f\n",R,z);
			if(6!=fscanf(ifile,"%lf %lf %lf %lf %lf %lf",&Epmra,&Epmdec,&Es,&EVlos,&EFeH,&EMgFe))
				printf("Err reading Epmra..: %f",Epmra);
			for(int i=0;i<nr;i++) reals.push_back(read_real(ifile));
			for(int i=0;i<nm;i++){
				mock m=read_mock(ifile);
				if(m.cpt==-1){
					printf("Err nm,i: %d %d\n",nm,i); exit(0);
				}
				mocks.push_back(m);
			}
		}
		void write(FILE* ofile){
			fprintf(ofile,"%zd %zd\n",reals.size(),mocks.size());
			fprintf(ofile,"%lf %lf %lf %lf %lf %lf\n",R,z,sigR,sigz,sigphi,vphibar);
			fprintf(ofile,"%lf %lf %lf %lf %lf %lf\n",Epmra,Epmdec,Es,EVlos,EFeH,EMgFe);
			for(int i=0;i<reals.size();i++) reals[i].write(ofile);
			for(int i=0;i<mocks.size();i++) mocks[i].write(ofile);
		}
};

#define  NC 6
#define NP 66

class chem_model{
	public:
		float Fbar[NC];
		float Mgbar[NC];
		float sigx[NC];
		float sigy[NC];
		float th[NC],c[NC],s[NC];
		float*** grads;
		double params[NP];
		std::vector<std::string> names;
		int *map;
		const int nvars;
		chem_model(int* _map, const int _nvars):
		    map(_map), nvars(_nvars) {
			grads = fmatrix(6,2,3);
			names.push_back("Fbar"); names.push_back("Mgbar"); names.push_back("sigx");
			names.push_back("sigy"); names.push_back("th"); names.push_back("FeGrad");
		}
		~chem_model(void){
			names.clear();
			delmatrix(grads,6,2);
		}
		void load_Fbar(float f0,float f1,float f2,float f3,float f4,float f5){
		    Fbar[0]=f0; Fbar[1]=f1; Fbar[2]=f2; Fbar[3]=f3; Fbar[4]=f4; Fbar[5]=f5;}
		void load_Mgbar(float f0,float f1,float f2,float f3,float f4,float f5){
		    Mgbar[0]=f0,Mgbar[1]=f1,Mgbar[2]=f2; Mgbar[3]=f3; Mgbar[4]=f4; Mgbar[5]=f5;}
		void load_sigx(float f0,float f1,float f2,float f3,float f4,float f5){
		    sigx[0]=f0; sigx[1]=f1; sigx[2]=f2; sigx[3]=f3; sigx[4]=f4; sigx[5]=f5;}
		void load_sigy(float f0,float f1,float f2,float f3,float f4,float f5){
		    sigy[0]=f0; sigy[1]=f1; sigy[2]=f2; sigy[3]=f3; sigy[4]=f4; sigy[5]=f5;}
		void load_th(float f0,float f1,float f2,float f3,float f4,float f5){
			th[0]=f0; th[1]=f1; th[2]=f2; th[3]=f3; th[4]=f4; th[5]=f5;}
		void load_grad(int n,float* C){
			for(int i=0; i<6; i++)
					grads[n][i/3][i%3]=C[i];
		}
		void write(std::string fname,const double vars[]){
			change_params(vars);
			unload_params();
			FILE* ofile;
			fopen_s(&ofile,fname.c_str(),"w");
			compress(ofile,Fbar,NC);
			compress(ofile,Mgbar,NC);
			compress(ofile,sigx,NC);
			compress(ofile,sigy,NC);
			compress(ofile,th,NC);
			compress(ofile,grads,NC,2,3);
			fclose(ofile);
		}
		bool read(std::string fname){
			FILE* ifile;
			if(fopen_s(&ifile,fname.c_str(),"r")){
				printf("I can't open %s\n",fname.c_str()); return false;
			}
			get(ifile,Fbar,NC);
			get(ifile,Mgbar,NC);
			get(ifile,sigx,NC);
			get(ifile,sigy,NC);
			get(ifile,th,NC);
			get(ifile,grads,NC,2,3);
			//printf("C(2,phi) %f\n",grads[0][1][2]);
			//grads[0][1][2]*=-1;
			fclose(ifile);
			for(int i=0; i<NC; i++){
				c[i]=cos(th[i]*M_PI/180); s[i]=sin(th[i]*M_PI/180);
			}
			return true;
		}
		void repeats(void){//equate some gradient parameters
			for(int i=1; i<3; i++)
				for(int j=0; j<2; j++)
					for(int k=0; k<3; k++)
						grads[i][j][k]=grads[0][j][k];
		}
		void load_params(void){
			int j=0;
			for(int i=0;i<NC;i++){//0-5
				params[j]=Fbar[i]; j++;
			}
			for(int i=0;i<NC;i++){//6-11
				params[j]=Mgbar[i]; j++;
			}
			for(int i=0;i<NC;i++){//12-17
				params[j]=sigx[i]; j++;
			}
			for(int i=0;i<NC;i++){//18-23
				params[j]=sigy[i]; j++;
			}
			for(int i=0;i<NC;i++){//24-29
				params[j]=th[i]; j++;
			}
			//repeats();//impose specified equality of gradients
			for(int i=0;i<NC;i++){//30-65
				for(int m=0;m<2;m++)
					for(int n=0;n<3;n++){
						params[j]=grads[i][m][n]; j++;
					}
			}
			if(j!=NP) printf("Error %d != %d parameters..",j,NP);
		}
		void unload_params(void){
			int j=0;
			for(int i=0;i<NC;i++){//0-5
				Fbar[i]=params[j]; j++;
			}
			for(int i=0;i<NC;i++){//6-11
				Mgbar[i]=params[j]; j++;
			}
			for(int i=0;i<NC;i++){//12-17
				sigx[i]=params[j]; j++;
			}
			for(int i=0;i<NC;i++){//18-23
				sigy[i]=params[j]; j++;
			}
			for(int i=0;i<NC;i++){//24-29
				th[i]=params[j];
				c[i]=cos(th[i]*M_PI/180); s[i]=sin(th[i]*M_PI/180);
				j++;
			}
			for(int i=0;i<NC;i++){//30-65
				for(int m=0;m<2;m++)
					for(int n=0;n<3;n++){
						grads[i][m][n]=params[j]; j++;
					}
			}
			//repeats();//required because some params out of date 
		}
		void print(void){
			for(int i=0;i<NC;i++){
			//	if(i<NC) continue;
				printf("\n%d\t",i);
				for(int j=0;j<2;j++){
					for(int k=0;k<3;k++)
						printf("%8.4f ",1000*grads[i][j][k]);
								}
			}
			printf("\n");
			for(int i=0; i<NC; i++) printf("%6.3f ",Fbar[i]); printf("\n");
			for(int i=0; i<NC; i++) printf("%6.3f ",Mgbar[i]); printf("\n");
			for(int i=0; i<NC; i++) printf("%6.3f ",sigx[i]); printf("\n");
			for(int i=0; i<NC; i++) printf("%6.3f ",sigy[i]); printf("\n");
			for(int i=0; i<NC; i++) printf("%6.3f ",th[i]); printf("\n");
		}
		void print_params(void){
			for(int i=0; i<NP; i++){
				printf("%8.4f ",params[i]);
				if(i%6==5) printf("\n");
			}
		}
		void change_params(const double vars[]){
			for(int i=0; i<nvars; i++) params[map[i]]=vars[i];
		}
		void get_xi(double vars[]){
			for(int i=0; i<nvars; i++) vars[i]=params[map[i]];
		}
		bool variable(int n){
			for(int i=0; i<nvars; i++) if(map[i]==n) return true;
			return false;
		}
		std::pair<double,double> chem(mock& star,double EFeH,double EMgFe){
			//Randomly chosen chemistry
			const double R0=8.3, Vc=240, Jsun=R0*Vc;
			double Jtr=star.Jt.Jr;
			double Jtz=star.Jt.Jz;
			double Jtphi=star.Jt.Jphi;
			int cpt=star.cpt;
			double r1,r2,r3,r4;
			math::getNormalRandomNumbers(r1,r2);
			math::getNormalRandomNumbers(r3,r4);
			double FebR =  Fbar[cpt] + r3*EFeH, x=r1*sigx[cpt];
			double MgbR = Mgbar[cpt] + r4*EMgFe, y=r2*sigy[cpt];
			FebR+=grads[cpt][0][0]*Jtr + grads[cpt][0][1]*Jtz + grads[cpt][0][2]*(Jtphi-Jsun);
			MgbR+=grads[cpt][1][0]*Jtr + grads[cpt][1][1]*Jtz + grads[cpt][1][2]*(Jtphi-Jsun);
			double Fe=FebR + x*c[cpt] - y*s[cpt];
			double Mg=MgbR + x*s[cpt] + y*c[cpt];
			return std::make_pair(Fe,Mg);
		}
		double P(double Fe,double Mg,mock& star){//PDF in (Fe,Mg) plane
			const double R0=8.3, Vc=240, Jsun=R0*Vc;
			double Jtr=star.Jt.Jr;
			double Jtz=star.Jt.Jz;
			double Jtphi=star.Jt.Jphi;
			int cpt=star.cpt;
			float Feb = Fbar[cpt];
			float Mgb = Mgbar[cpt];
			Feb+=grads[cpt][0][0]*Jtr + grads[cpt][0][1]*Jtz + grads[cpt][0][2]*(Jtphi-Jsun);
			Mgb+=grads[cpt][1][0]*Jtr + grads[cpt][1][1]*Jtz + grads[cpt][1][2]*(Jtphi-Jsun);
			float x = (Fe-Feb)*c[cpt]+(Mg-Mgb)*s[cpt];
			float y =-(Fe-Feb)*s[cpt]+(Mg-Mgb)*c[cpt];
			return exp(-.5*(pow_2(x/sigx[cpt])+pow_2(y/sigy[cpt])));
		}
		double prior(void){
			double p=0;
			for(int i=0; i<NC; i++){
				p+=pow(sigx[i],2);
				p+=pow(10*sigy[i],2);
				for(int k=0; k<2; k++)
					for(int j=0; j<2; j++)
						p+=pow(1000*grads[i][k][j],2);
			}
			return -p;
		}
};
void tab_chemParams(std::string& fname,	std::vector<std::string>& cpts, chem_model& M){
	std::vector<std::string> nm;
	nm.push_back("$x_0$"); nm.push_back("$y_0$"); nm.push_back("$\\sigma_x$");
	nm.push_back("$\\sigma_y$"); nm.push_back("$\\theta$");
	int map2[5]={4,0,1,2,3};//order to tabulate parameters
	std::ofstream strm(fname);
	strm << "Component ";
	for(int i=0; i<5; i++) strm <<"\t & " << nm[map2[i]];
	strm <<"\t& $C_{1,J_r}$\t& $C_{1,J_z}$\t& $C_{1,J_\\phi}$\t& $C_{2,J_r}$ \t& $C_{2,J_z}$ \t& $C_{2,J_\\phi}$\\cr \n \\hline \n";
	for(int ic=0; ic<NC; ic++){//cpts
		for(int ip=0; ip<5; ip++){//parameters
			int n=map2[ip]*NC+ic;
/*			if(M.variable(n)){
				if(ip==0) strm << cpts[ic] << std::setprecision(2) << " \t&${\\bf " << M.params[n] <<"}$";
				else strm << "\t&${\\bf " << M.params[n] << "}$";
			}else{*/
				if(ip==0) strm << cpts[ic] << std::setprecision(3) << "\t&$" << M.params[n] <<"$";
				else strm << "\t&$" << M.params[n] << "$";
//			}
		}
		for(int row=0; row<2; row++){
			for (int col=0; col<3; col++) strm << std::setprecision(3) << "\t&$" <<  1000*M.grads[ic][row][col] << "$";
			//if(row==1) strm << "\\cr\n";
		}
		strm << "\\cr\n";
	}
}
class Lhd: public math::IFunctionNdim {
	private:
		std::vector<utils::KeyValueMap>* Kvs;
		std::vector<int>** changing;
		chem_model* M;
		FILE* log_file;
		const int nvars, nvars_chem; 
		potential::PtrPotential pot;
		const int nx,ny,nc;
		int nreal, nmock;//total number of real or mock stars
		const float Femin, Femax, Mgmin, Mgmax, JphiMax, JphiMin, JzMax;
		std::vector<double>& Ms0;//masses of the cpts from DFs
		const double stellarMass;
		std::vector<cell>& cells;
		std::vector<double> weights;//ratio # reals/# mocks in each cell
		std::vector<grid2>& grdso;
		bool plot,dyn;
//		mgo::plt *pl, *pl1, *pl2;
		void do_plots(std::vector<grid2>&) const;
	public:
		Lhd(std::vector<utils::KeyValueMap>* _Kvs,const double _stellarMass,
		    std::vector<int>** _changing,
		    bool _dyn, chem_model* _M,
		    const int _nvars_chem, const int _nvars, potential::PtrPotential _pot,
		    std::vector<grid2>& _grdso,const int _nx, const int _ny,const int _nc,
		    float _Femin,float _Femax, float _Mgmin,float _Mgmax,
		    float _JphiMin, float _JphiMax, float _JzMax,
		    std::vector<double>& _Ms0, std::vector<cell>& _cells) :
		    Kvs(_Kvs), stellarMass(_stellarMass), changing(_changing),
		    dyn(_dyn), M(_M), nvars_chem(_nvars_chem), nvars(_nvars),
		    pot(_pot), grdso(_grdso), nx(_nx), ny(_ny), nc(_nc), 
		    Femin(_Femin), Femax(_Femax), Mgmin(_Mgmin), Mgmax(_Mgmax),
		    JphiMin(_JphiMin), JphiMax(_JphiMax), JzMax(_JzMax),
		    Ms0(_Ms0), cells(_cells)
		{
			plot = false;
			fopen_s(&log_file,"Lhd.log","w");
			nreal=0; nmock=0;
			for(size_t i=0; i<cells.size(); i++){
				weights.push_back((double)cells[i].reals.size()/(double)cells[i].mocks.size());
				nreal+=cells[i].reals.size();
				nmock+=cells[i].mocks.size();
			}
		}
		virtual void eval(const double vars[], double values[]) const;
		virtual unsigned int numVars() const{
			return nvars;
		}
		virtual unsigned int numValues() const{
			return 1;
		}
		void plot_grid(const double*,int icell) const;
		double do_plot(double vars[]){
			plot=true;
			double values[1];
			eval(vars,values);
			plot=false;
			return values[0];
		}
		void plot_chem(grid4&) const;
		void toggle_plot(void){
			plot=!plot;
		}
		void contourM(int nx, mock& star) const{
			float** m; m=fmatrix(nx,nx);
			for(int i=0;i<nx;i++){
				float FeH=Femin+i/(float)(nx-1)*(Femax-Femin);
				for(int j=0;j<nx;j++){
					float MgFe=Mgmin+j/(float)(nx-1)*(Mgmax-Mgmin);
					m[i][j]=M->P(FeH,MgFe,star);
				}
			}
			const int nc=5; float h[nc]={.05,.1,.2,.4,.8};
			//pl1->contour(m,nx,nx,h,nc);
			delmatrix(m,nx);
		}
		void testM(void) const{
			double _Jr=20,_Jphi=1000,_Jz=35;
			for(int cpt=0; cpt<6; cpt++){
				grid2 grd(50,50,0,Femin,Femax,Mgmin,Mgmax);
				mock star(cpt,actions::Actions(_Jr,_Jz,_Jphi),
					  actions::Actions(_Jr,_Jz,_Jphi),coord::VelCyl(0,0,0),coord::VelCyl(0,0,0),0);
				//pl1->new_plot(Femin,Femax,Mgmin,Mgmax,"[Fe/H}","[Mg/Fe]");
				for(int i=0;i<10000;i++){
					std::pair<double,double> chm = M->chem(star,0,0);
					grd.CIC(1,chm.first,chm.second);
				}
				//grd.plot_dens(pl1,true,false);
				contourM(50,star);
				//pl1->grend();
			}
		}		

};
coord::PosVelCyl scatter(obs::solarShifter& sun,errors& errs,coord::PosVelCyl xv,double& sKpctrue){
	double Vlos_kms;
	obs::PosVelSky pospm(sun.toSky(xv,sKpctrue,Vlos_kms));//get l,b and mul etc
	obs::PosVelSky radec=to_muRAdec(pospm);//convert to ra,dec,mura,mudec
	float Epmra,Epmdec,Es,EVlos,EFeH,EMgFe;
	errs.get(sKpctrue,Epmra,Epmdec,Es,EVlos,EFeH,EMgFe);
	Es*=.5;
	double r1,r2;
	math::getNormalRandomNumbers(r1,r2);
	//now add errors to radec cpts
	radec.pm.mul+=r1*Epmra; radec.pm.mub+=r2*Epmdec;
	obs::PosVelSky xvSky(obs::from_muRAdec(radec));//convert scattered cords back
	//add errors to skpc and Vlos_kms
	math::getNormalRandomNumbers(r1,r2);
	double sKpc=sKpctrue+r1*Es; Vlos_kms+=r2*EVlos;
	return coord::PosVelCyl(sun.toCyl(xvSky.pos,sKpc,xvSky.pm,Vlos_kms));
}

void read_header(char *line){
	std::size_t Pos=0, pos=0;
	int col=0;
	while(pos!=std::string::npos){
		std::string ln(line+Pos);
		pos=ln.find(" ");
		Pos+=pos+1;
		printf("%d |%s|",col-1,(ln.substr(0,pos)).c_str());
		col++;
	}
}
int read_extra(char* line,star& Star,long int& flag,float& fidelity,float &SNR,char* GC){
	std::size_t Pos=0, pos=0;
	int col=0;
	while(pos!=std::string::npos){
		std::string ln(line+Pos);
		pos=ln.find(" ");
		Pos+=pos+1;
		if     (col==0) sscanf((ln.substr(0,pos)).c_str(),"%s",&(Star.id));
		else if(col==1) sscanf((ln.substr(0,pos)).c_str(),"%lf",&(Star.ra));
		else if(col==2) sscanf((ln.substr(0,pos)).c_str(),"%lf",&Star.dec);
		else if(col==3) sscanf((ln.substr(0,pos)).c_str(),"%f",&(Star.Vlos));
		else if(col==4) sscanf((ln.substr(0,pos)).c_str(),"%f",&(Star.DVlos));
		else if(col==5) sscanf((ln.substr(0,pos)).c_str(),"%f",&(Star.Vlos_err));
		else if(col==6) sscanf((ln.substr(0,pos)).c_str(),"%f",&(Star.Teff));
		else if(col==7) sscanf((ln.substr(0,pos)).c_str(),"%f",&(Star.logg));
		else if(col==8) sscanf((ln.substr(0,pos)).c_str(),"%f",&(Star.FeH));
		else if(col==9) sscanf((ln.substr(0,pos)).c_str(),"%f",&(Star.FeH_err));
		else if(col==10) sscanf((ln.substr(0,pos)).c_str(),"%f",&(Star.MgFe));
		else if(col==11) sscanf((ln.substr(0,pos)).c_str(),"%f",&(Star.MgFe_err));
		else if(col==12) sscanf((ln.substr(0,pos)).c_str(),"%f",&(Star.FeH_nn));
		else if(col==13) sscanf((ln.substr(0,pos)).c_str(),"%f",&(Star.FeH_nn_err));
		else if(col==14) sscanf((ln.substr(0,pos)).c_str(),"%f",&(Star.MgFe_nn));
		else if(col==15) sscanf((ln.substr(0,pos)).c_str(),"%f",&(Star.MgFe_nn_err));
		else if(col==16) sscanf((ln.substr(0,pos)).c_str(),"%f",&(Star.s_sh));
		else if(col==17) sscanf((ln.substr(0,pos)).c_str(),"%f",&(Star.s_sh_err));
		else if(col==18) sscanf((ln.substr(0,pos)).c_str(),"%f",&(Star.s_nn));
		else if(col==19) sscanf((ln.substr(0,pos)).c_str(),"%f",&(Star.s_nn_err));
		else if(col==20) sscanf((ln.substr(0,pos)).c_str(),"%d",&flag);
		else if(col==21) sscanf((ln.substr(0,pos)).c_str(),"%f",&SNR);
		else if(col==22) sscanf((ln.substr(0,pos)).c_str(),"%f",&fidelity);
		else if(col==23) sscanf((ln.substr(0,pos)).c_str(),"%s",GC);
		else if(col==24) sscanf((ln.substr(0,pos)).c_str(),"%f",&(Star.pi));
		else if(col==25) sscanf((ln.substr(0,pos)).c_str(),"%f",&(Star.pi_err));
		else if(col==26) sscanf((ln.substr(0,pos)).c_str(),"%f",&(Star.pm_ra));
		else if(col==27) sscanf((ln.substr(0,pos)).c_str(),"%f",&(Star.pm_ra_err));
		else if(col==28) sscanf((ln.substr(0,pos)).c_str(),"%f",&(Star.pm_dec));
		else if(col==29) sscanf((ln.substr(0,pos)).c_str(),"%f",&(Star.pm_dec_err));
		col++;
	}
	return col;
}
bool set7(long int n){
	long int N = n & 1<<7;
	return N!=0? true : false;
}
// read_17 reads DR17 data
int read_17(bool bothways,int sense,obs::solarShifter& sun,const units::InternalUnits& intUnits,
	    std::vector<long_real>& reals,errors& errs,potential::PtrPotential& pot,actions::PtrActionFinder& AF){
	FILE* ifile;
	if(fopen_s(&ifile,"/APOGEE/APOGEE_DR17_extra.txt","r")){
		printf("I can't open data file\n"); return 0;
	}
	const int N=1000;
	char line[N];
	int il=0, n_sh=0, N_GC=0, Nfree=0;
	fgets(line,N,ifile);
	read_header(line);
	while(!feof(ifile)){
		if(!fgets(line,N,ifile)) break;
		star Star;
		float fidelity,SNR;
		long int flag;
		char GC[30];
		read_extra(line,Star,flag,fidelity,SNR,GC);
		int nGC=strcmp("-",GC);
		if(set7(flag) || isnan(Star.s_sh) || SNR<50 || isnan(Star.logg)
		   || isnan(Star.Teff) || isnan(Star.FeH_err) || isnan(Star.MgFe_err) || nGC){//grounds for rejection
			if(nGC) N_GC++;
		} else {
			Star.MgFe_nn-=Star.FeH_nn;
			n_sh++;
			if(Star.logg<3.5 && Star.Teff<5500 && fidelity>0.5 && (Star.s_sh_err < .75)){
				obs::PosVelSky skyPosVel(obs::from_muRAdec(Star.ra,Star.dec,Star.pm_ra,Star.pm_dec));
				coord::PosVelCyl xv(sun.toCyl(skyPosVel.pos,Star.s_sh,skyPosVel.pm,Star.Vlos));
				if(!isnan(xv.vphi) && (bothways || sense*xv.vphi<0)){
					double E=totalEnergy(*pot,xv);
					if(E>=0){
						Nfree++; continue;
					}
					actions::Actions J(AF->actions(xv));
					if(!isnan(J.Jz)){
						reals.push_back(long_real(Star.FeH,Star.MgFe,
							xv.R*intUnits.to_Kpc,xv.z*intUnits.to_Kpc,
							xv.phi,
						xv.vR*intUnits.to_kms,xv.vz*intUnits.to_kms,
						xv.vphi*intUnits.to_kms,
						J.Jr*(intUnits.to_kms*intUnits.to_Kpc),
						J.Jz*(intUnits.to_kms*intUnits.to_Kpc),
						Star.pm_ra_err,Star.pm_dec_err,Star.s_sh_err/Star.s_sh,
							Star.Vlos_err,Star.FeH_err,Star.MgFe_err));
						errs.add(Star.s_sh,Star.pm_ra_err,Star.pm_dec_err,
							Star.s_sh_err,Star.Vlos_err,Star.FeH_err,Star.MgFe_err);
					}
				}
			}
		}
		il++;
	}
	printf("%d stars read, %d in StarHorse, of which %d free and %zd in bound giants\n",il,n_sh,Nfree,reals.size());
	return reals.size();
}
float pack_cells(int nR, int nz, float Rs[], float zs[],
		   std::vector<long_real>& reals, std::vector<cell>& cells){
	for(size_t i=0; i<reals.size(); i++){
		int iR=0,iz=0;
		while(reals[i].R>Rs[iR] && iR<nR) iR++;
		if(iR<1 || iR==nR) continue; iR--;
		while(fabs(reals[i].z)>zs[iz] && iz<nz) iz++;
		if(iz<1 || iz==nz) continue; iz--;
		int indx=(nR-1)*iz+iR;
		cells[indx].R += reals[i].R; cells[indx].z += fabs(reals[i].z);
		cells[indx].sigR += pow_2(reals[i].VR); cells[indx].sigz += pow_2(reals[i].Vz);
		cells[indx].vphibar += reals[i].Vphi; cells[indx].sigphi += pow_2(reals[i].Vphi);
		cells[indx].phis.push_back(reals[i].phi);
		cells[indx].Epmra += reals[i].Epmra; cells[indx].Epmdec += reals[i].Epmdec;
		cells[indx].Es += reals[i].Es; cells[indx].EVlos += reals[i].EVlos;
		cells[indx].EFeH += reals[i].EFeH; cells[indx].EMgFe += reals[i].EMgFe;
		cells[indx].load_real(reals[i]);
	}
	float Nfiled=0;
	for(size_t i=0; i<cells.size(); i++){
		cells[i].R /= cells[i].reals.size();
		cells[i].z /= cells[i].reals.size();
		cells[i].sigR /= cells[i].reals.size(); cells[i].sigR = sqrt(cells[i].sigR);
		cells[i].sigz /= cells[i].reals.size(); cells[i].sigz = sqrt(cells[i].sigz);
		cells[i].vphibar /= cells[i].reals.size();
		cells[i].sigphi /= cells[i].reals.size();
		cells[i].sigphi = sqrt(fmax(0,cells[i].sigphi-pow_2(cells[i].vphibar)));
		cells[i].Epmra /= cells[i].reals.size();
		cells[i].Epmdec /= cells[i].reals.size();
		cells[i].Es /= cells[i].reals.size();
		cells[i].EVlos /= cells[i].reals.size();
		cells[i].EFeH /= cells[i].reals.size();
		cells[i].EMgFe /= cells[i].reals.size();
		Nfiled+=cells[i].reals.size();
	}
	return Nfiled;
}
int read_mocks(obs::solarShifter& sun, int nc, const std::string dir, std::vector<double>& Ms0,
	       std::vector<cell>& cells,
	       errors& errs,potential::PtrPotential& pot,actions::PtrActionFinder& AF){
//	char fname[30];
//	sprintf(fname,"../%s/Hayden",dir);
	std::string fname = dir + "Hayden";
	FILE *ifile;
	if(fopen_s(&ifile,fname.c_str(),"r")){
		printf("I can't open %s to read cpt masses\n",fname.c_str()); return 0;
	}
	for(int i=0; i<nc; i++){
		double M;
		if(1!=fscanf(ifile,"%lg",&M)){
			printf("Not enough mases in %s\n",fname.c_str()); return 0;
		}
		Ms0.push_back(M);
	}
	fclose(ifile);
	const int nh=cells.size();
	int nstar=0, nfree=0;
	std::vector<double> sbar(nh);
	errs.report("Gaia_errors");
//#pragma omp parallel for schedule(dynamic)
	for(int ih=0; ih<nh; ih++){//run over locations
		int nphis=cells[ih].phis.size();
		sbar[ih]=0;		
		for(int ic=0; ic<nc; ic++){//run over cpts
			char fn[80];
			sprintf(fn,"%s%d_%d",fname.c_str(),ih,ic);
			FILE *ifile;
			if(fopen_s(&ifile,fn,"r")){
				printf("I can't open %s\n",fn); return 0;
			}
			double R,z;
			int cp,np;
			fscanf(ifile,"%lf %lf %d %d",&R,&z,&cp,&np); nstar+=np;
			double Rint=R*intUnits.from_Kpc, zint=z*intUnits.from_Kpc;
			if(ic!=cp){
				printf("Error: ic,cp = %d %d\n",ic,cp); return 0;
			}
			std::vector<coord::PosVelCyl> xvs(np),xvt(np);
			std::vector<actions::Actions> Jt(np),J(np);
			std::vector<double> f0(np);
			std::vector<double> sKpc(np);
			for(int ip=0;ip<np;ip++){
				int nphi=nphis*math::random();
				double phi=cells[ih].phis[nphi];//assign phi of randomly chosen real
				double VR,Vz,Vphi,Jtr,Jtz,Jtphi;
				if(7!=fscanf(ifile,"%lf %lf %lf %lf %lf %lf %lf",&VR,&Vz,&Vphi,&Jtr,&Jtz,&Jtphi,&f0[ip])){
					printf("failed to read V,Jt in %s\n",fn);//read V anf f for this star
					return 0;
				}
				xvt[ip]=coord::PosVelCyl(Rint,zint,phi,
					VR*intUnits.from_kms,Vz*intUnits.from_kms,Vphi*intUnits.from_kms);
				Jt[ip] = actions::Actions(Jtr,Jtz,Jtphi);//true actions
				xvs[ip]= scatter(sun,errs,xvt[ip],sKpc[ip]);
			}
			//printf("computing apparent Js\n");
			for(int ip=0; ip<np; ip++){
				J[ip] =  AF->actions(xvs[ip]);//apparent actions
				if(Jt[ip].Jz<0 || J[ip].Jz<0){
					nfree++; continue;
				}
				J[ip].Jr/=from_kpckms; J[ip].Jz/=from_kpckms; J[ip].Jphi/=from_kpckms;
				coord::VelCyl V(xvs[ip].vR*intUnits.to_kms, xvs[ip].vz*intUnits.to_kms,
						xvs[ip].vphi*intUnits.to_kms);
				coord::VelCyl Vt(xvt[ip].vR*intUnits.to_kms, xvt[ip].vz*intUnits.to_kms,
						xvt[ip].vphi*intUnits.to_kms);
				cells[ih].load_mock(Jt[ip], J[ip], V, Vt, f0[ip], ic);
				if(isnan(sKpc[ip])) printf("NaN: %g ",xvs[ip].R);
				sbar[ih]+=sKpc[ip];
			}
			fclose(ifile);
		}
		sbar[ih]/=(double)cells[ih].mocks.size();
	}
	printf("%d mock stars loaded of which %d free. Mean distances:\n",nstar,nfree);
	return nstar;
}

void V_hists(FILE* ofile,cell& Cell,
	     float VRmax,float Vzmax,float Vphimin,float Vphimax,int nb){
	std::vector<float> NR(nb),Nz(nb),Nphi(nb);
	float dVR=2*VRmax/(float)nb,dVz=2*Vzmax/(float)nb;
	float dVphi=(Vphimax-Vphimin)/(float)nb;
	for(int k=0; k<nb; k++){
		NR[k]=0; Nz[k]=0; Nphi[k]=0;
	}
	for(int i=0; i<Cell.reals.size(); i++){
		int iR=(Cell.reals[i].VR+VRmax)/dVR, iz=(Cell.reals[i].Vz+Vzmax)/dVz;
		int ip=(-Cell.reals[i].Vphi-Vphimin)/dVphi;
		if(iR>-1 && iR<nb) NR[iR]+=1;
		if(iz>-1 && iz<nb) Nz[iz]+=1;
		if(ip>-1 && ip<nb) Nphi[ip]+=1;
	}
	fprintf(ofile,"%d %f %f %f %f %f %f\n",nb,Cell.R,Cell.z,VRmax,Vzmax,Vphimin,Vphimax);
	compress(ofile,NR); compress(ofile,Nz);	compress(ofile,Nphi);//data histograms
	for(int k=0; k<nb; k++){
		NR[k]=0; Nz[k]=0; Nphi[k]=0;
	}
	float fac = (float)Cell.reals.size()/(float)Cell.mocks.size();
	for(int i=0; i<Cell.mocks.size(); i++){
//		int iR=(Cell.mocks[i].Vt.vR+VRmax)/dVR, iz=(Cell.mocks[i].Vt.vz+Vzmax)/dVz;
//		int ip=(Cell.mocks[i].Vt.vphi-Vphimin)/dVphi;
		int iR=(Cell.mocks[i].V.vR+VRmax)/dVR, iz=(Cell.mocks[i].V.vz+Vzmax)/dVz;
		int ip=(Cell.mocks[i].V.vphi-Vphimin)/dVphi;
		int cpt=Cell.mocks[i].cpt;
		//if(cpt==5) continue;
		if(iR>-1 && iR<nb) NR[iR]+=fac;
		if(iz>-1 && iz<nb) Nz[iz]+=fac;
		if(ip>-1 && ip<nb) Nphi[ip]+=fac;
	}
	compress(ofile,NR); compress(ofile,Nz);	compress(ofile,Nphi);//model histograms
}
