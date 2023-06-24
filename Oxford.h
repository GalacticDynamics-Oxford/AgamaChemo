#pragma once
#include "df_base.h"
#include "units.h"
#include "math_core.h"
#include "math_ode.h"
#include "potential_utils.h"
#include <fstream>
#include <cmath>
#define EXP __declspec(dllexport)

/// Now the Pascale, Binney, Nipoti (2020 DF
EXP double Oxford_r_rat(const double[],double);
EXP double Bologna_rat(const double[],double);

class EXP Oxford_z_syst: public math::IOdeSystem{
	private:
		double z_rat;
	public:
		Oxford_z_syst(double _z_rat) : z_rat(_z_rat) {};
		void eval(const double,const double[],double[]) const;
		unsigned int size()  const {return 2;}
};

class EXP Oxford_r_syst: public math::IOdeSystem{
	private:
		double r_rat;
	public:
		Oxford_r_syst(double _r_rat) : r_rat(_r_rat) {};
		void eval(const double,const double[],double[]) const;
		unsigned int size()  const {return 2;}
};

struct OxfordParam{
	double
			mass,      ///< total mass in Msun
			J0,        ///< break action (defines the transition between inner and outer regions)
			Jcore,        ///< controls Cole&Binney core size
			Jcutoff,   ///< cutoff action (sets exponential suppression at J>Jcutoff, 0 to disable)
			slopeIn,   ///< power-law index for actions below the break action (Gamma)
			slopeOut,  ///< power-law index for actions above the break action (Beta)
			steepness, ///< steepness of the transition between two asymptotic regimes (eta)
			cutoffStrength, ///< steepness of exponential suppression at J>Jcutoff (zeta)
			coefLin,  ///< radial anisotropy
			kIn,      ///< flattening
			coefLout,  ///< radial anisotropy
			kOut,      ///< flattening
			L0,        ///< epicycle freqs vs r
			gamma0,     ///< epicycle freqs as r->0
			gamma1,     ///< epicycle freqs as r->0
			bz,	   ///< flattening
			rotFrac,   ///< relative amplitude of the odd-Jphi component (-1 to 1, 0 means no rotation)
			Jphi0;     ///< controls the steepness of rotation and the size of non-rotating core

	OxfordParam() :  ///< set default values for all fields (NAN means that it must be set manually)
	    mass(NAN), J0(NAN), Jcore(0), Jcutoff(0), slopeIn(NAN), slopeOut(NAN), steepness(1), cutoffStrength(2),
	    coefLin(1), kIn(1), coefLout(1), kOut(1), L0(NAN), gamma0(0), gamma1(1), bz(0), rotFrac(0), Jphi0(0) {}
};
class EXP Oxford: public BaseDistributionFunction{
	private:
		const OxfordParam par;  ///< parameters of DF
		const potential::Interpolator freq;  ///< interface providing the epicyclic frequencies and Rcirc
		double norm,beta;
		double g(const double hJ) const;
		void set_beta(void);
		Oxford(const OxfordParam &params_, const potential::Interpolator& freq_, const double beta_) :
		    par(params_), freq(freq_), beta(beta_){}
	public:
    /** Create an instance of ColeBinney double-power-law distribution function with given parameters
        \param[in] params  are the parameters of DF
        \throws std::invalid_argument exception if parameters are nonsense
    */
		Oxford(const OxfordParam&, const potential::Interpolator&);

    /** return value of DF for the given set of actions.
        \param[in] J are the actions  */
		virtual double value(const actions::Actions &J) const;
		virtual void set_norm(double);
		void write_params(std::ofstream&,const units::InternalUnits&) const;
		void tab_params(std::ofstream&,const units::InternalUnits&) const;
		double diff(const double) const;
		double intDiff(const double) const;
		double h(const actions::Actions&) const;
		double z_rat(const actions::Actions&) const;
		double r_rat(const actions::Actions&) const;
};

//We start with two helper classes used to determine Cole&Binney
//parameter beta

class EXP OxfordDiff: public math::IFunctionNoDeriv{
	public:
		const Oxford& Sft;
		OxfordDiff(Oxford &Sft_): Sft(Sft_) {}
		virtual double value(const double hJ) const{
			return Sft.diff(hJ);
		}
};
class EXP OxfordInt: public math::IFunctionNoDeriv{
	public:
		const Oxford& Sft;
		OxfordInt(Oxford &Sft_): Sft(Sft_) {}
		virtual double value(const double x) const{
			return Sft.intDiff(x);
		}
};

struct BolognaParam {
	double
			mass,      ///< total mass of a model
			J0,        ///< controls core size
			Jphi,      ///< controls the steepness of rotation and the size of non-rotating core
			L0,        ///< controls approx to Omega_phi/Omega_r
			gamma,     ///< epicycle ratio r->0
			L1,        ///< introduces harmonic core
			coefL,      ///< control radial anisotropy
			kIn, kOut,  ///< control flattening
			alpha,      ///< controls out density profile 
			rotFrac;    ///< relative amplitude of the odd-Jphi component (-1 to 1, 0 means no rotation)
	BolognaParam():     ///< set default values for all fields (NAN means that it must be set manually)
	    mass(NAN), J0(NAN), Jphi(0), L0(NAN), gamma(0), L1(0), coefL(0), kIn(1), kOut(1),
	    alpha(0.5), rotFrac(0) {}
};

class EXP Bologna: public BaseDistributionFunction{
	const BolognaParam par;  ///< parameters of DF
	private:
		double norm;
	public:
    /** Create an instance of Pascale Exponential distribution function with given parameters
        \param[in] params  are the parameters of DF
        \throws std::invalid_argument exception if parameters are nonsense
    */
		Bologna(const BolognaParam &params);

    /** return value of DF for the given set of actions.
        \param[in] J are the actions  */
		virtual double value(const actions::Actions &J) const;
};

