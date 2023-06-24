#include "df_factory.h"
#include "df_disk.h"
#include "df_halo.h"
#include "df_spherical.h"
#include "utils.h"
#include <cassert>
#include <stdexcept>

namespace df {

EXP DoublePowerLawParam parseDoublePowerLawParam(
					     const utils::KeyValueMap& kvmap,
					     const units::ExternalUnits& conv)
{
	DoublePowerLawParam par;
	par.norm      = kvmap.getDouble("norm",      par.norm)    * conv.massUnit;
	par.J0        = kvmap.getDouble("J0",        par.J0)      * conv.lengthUnit * conv.velocityUnit;
	par.Jcutoff   = kvmap.getDouble("Jcutoff",   par.Jcutoff) * conv.lengthUnit * conv.velocityUnit;
	par.Jphi0     = kvmap.getDouble("Jphi0",     par.Jphi0)   * conv.lengthUnit * conv.velocityUnit;
	par.Jcore     = kvmap.getDouble("Jcore",     par.Jcore)   * conv.lengthUnit * conv.velocityUnit;
	par.slopeIn   = kvmap.getDouble("slopeIn",   par.slopeIn);
	par.slopeOut  = kvmap.getDouble("slopeOut",  par.slopeOut);
	par.steepness = kvmap.getDouble("steepness", par.steepness);
	par.coefJrIn  = kvmap.getDouble("coefJrIn",  par.coefJrIn);
	par.coefJzIn  = kvmap.getDouble("coefJzIn",  par.coefJzIn);
	par.coefJrOut = kvmap.getDouble("coefJrOut", par.coefJrOut);
	par.coefJzOut = kvmap.getDouble("coefJzOut", par.coefJzOut);
	par.rotFrac   = kvmap.getDouble("rotFrac",   par.rotFrac);
	par.cutoffStrength = kvmap.getDouble("cutoffStrength", par.cutoffStrength);
	par.Fname = kvmap.getString("PopFile", par.Fname);
	return par;
}

EXP ModifiedDoublePowerLawParam parseModifiedDoublePowerLawParam(
	const utils::KeyValueMap& kvmap,
	const units::ExternalUnits& conv)
{
	ModifiedDoublePowerLawParam par;
	par.norm      = kvmap.getDouble("norm",      par.norm)    * conv.massUnit;
	par.J0        = kvmap.getDouble("J0",        par.J0)      * conv.lengthUnit * conv.velocityUnit;
	par.Jcutoff   = kvmap.getDouble("Jcutoff",   par.Jcutoff) * conv.lengthUnit * conv.velocityUnit;
	par.Jphi0     = kvmap.getDouble("Jphi0",     par.Jphi0)   * conv.lengthUnit * conv.velocityUnit;
	par.Jcore     = kvmap.getDouble("Jcore",     par.Jcore)   * conv.lengthUnit * conv.velocityUnit;
	par.L0	      = kvmap.getDouble("L0",        par.L0)      * conv.lengthUnit * conv.velocityUnit;
	par.slopeIn   = kvmap.getDouble("slopeIn",   par.slopeIn);
	par.slopeOut  = kvmap.getDouble("slopeOut",  par.slopeOut);
	par.cutoffStrength = kvmap.getDouble("cutoffStrength", par.cutoffStrength);
	par.coefJrIn  = kvmap.getDouble("coefJrIn",  par.coefJrIn);
	par.coefJzIn  = kvmap.getDouble("coefJzIn",  par.coefJzIn);
	par.rotFrac   = kvmap.getDouble("rotFrac",   par.rotFrac);
	par.Fname = kvmap.getString("PopFile", par.Fname);
	return par;
}

EXP SinDoublePowerLawParam parseSinDoublePowerLawParam(
	const utils::KeyValueMap& kvmap,
	const units::ExternalUnits& conv)
{
	SinDoublePowerLawParam par;
	par.norm      = kvmap.getDouble("norm",      par.norm)    * conv.massUnit;
	par.mass      = kvmap.getDouble("mass",      par.mass)    * conv.massUnit;
	par.J0        = kvmap.getDouble("J0",        par.J0)      * conv.lengthUnit * conv.velocityUnit;
	par.Jcutoff   = kvmap.getDouble("Jcutoff",   par.Jcutoff) * conv.lengthUnit * conv.velocityUnit;
	par.Jphi0     = kvmap.getDouble("Jphi0",     par.Jphi0)   * conv.lengthUnit * conv.velocityUnit;
	par.Jcore     = kvmap.getDouble("Jcore",     par.Jcore)   * conv.lengthUnit * conv.velocityUnit;
	par.L0	      = kvmap.getDouble("L0",        par.L0)      * conv.lengthUnit * conv.velocityUnit;
	par.slopeIn   = kvmap.getDouble("slopeIn",   par.slopeIn);
	par.slopeOut  = kvmap.getDouble("slopeOut",  par.slopeOut);
	par.cutoffStrength = kvmap.getDouble("cutoffStrength", par.cutoffStrength);
	par.alpha     = kvmap.getDouble("alpha",     par.alpha);
	par.beta      = kvmap.getDouble("beta",      par.beta);
	par.Fin       = kvmap.getDouble("Fin",       par.Fin);
	par.Fout      = kvmap.getDouble("Fout",      par.Fout);
	par.rotFrac   = kvmap.getDouble("rotFrac",   par.rotFrac);
	par.Fname     = kvmap.getString("PopFile",   par.Fname);
	return par;
}

EXP QuasiIsothermalParam parseQuasiIsothermalParam(
    const utils::KeyValueMap& kvmap,
    const units::ExternalUnits& conv)
{
    QuasiIsothermalParam par;
    par.Sigma0  = kvmap.getDouble("Sigma0",  par.Sigma0)  * conv.massUnit / pow_2(conv.lengthUnit);
    par.Rdisk   = kvmap.getDouble("Rdisk",   par.Rdisk)   * conv.lengthUnit;
    par.Hdisk   = kvmap.getDouble("Hdisk",   par.Hdisk)   * conv.lengthUnit;
    par.sigmar0 = kvmap.getDouble("sigmar0", par.sigmar0) * conv.velocityUnit;
    par.sigmaz0 = kvmap.getDouble("sigmaz0", par.sigmaz0) * conv.velocityUnit;
    par.sigmamin= kvmap.getDouble("sigmamin",par.sigmamin)* conv.velocityUnit;
    par.Rsigmar = kvmap.getDouble("Rsigmar", par.Rsigmar) * conv.lengthUnit;
    par.Rsigmaz = kvmap.getDouble("Rsigmaz", par.Rsigmaz) * conv.lengthUnit;
    par.coefJr  = kvmap.getDouble("coefJr",  par.coefJr);
    par.coefJz  = kvmap.getDouble("coefJz",  par.coefJz);
    par.Jmin    = kvmap.getDouble("Jmin",    par.Jmin)    * conv.lengthUnit * conv.velocityUnit;
    par.beta    = kvmap.getDouble("beta",    par.beta);
    par.Tsfr    = kvmap.getDouble("Tsfr",    par.Tsfr);  // dimensionless! in units of galaxy age
    par.sigmabirth = kvmap.getDouble("sigmabirth", par.sigmabirth);  // dimensionless ratio
    return par;
}
/*
ExponentialParam parseExponentialParam(
				       const utils::KeyValueMap& kvmap,
				       const units::ExternalUnits& conv)
{
	ExponentialParam par;
	par.norm   = kvmap.getDouble("norm",   par.norm)   * conv.massUnit;
	par.Jr0    = kvmap.getDouble("Jr0",    par.Jr0)    * conv.lengthUnit * conv.velocityUnit;
	par.Jz0    = kvmap.getDouble("Jz0",    par.Jz0)    * conv.lengthUnit * conv.velocityUnit;
	par.Jphi0  = kvmap.getDouble("Jphi0",  par.Jphi0)  * conv.lengthUnit * conv.velocityUnit;
	par.addJden= kvmap.getDouble("addJden")* conv.lengthUnit * conv.velocityUnit;
	par.addJvel= kvmap.getDouble("addJvel")* conv.lengthUnit * conv.velocityUnit;
	par.coefJr = kvmap.getDouble("coefJr", par.coefJr);
	par.coefJz = kvmap.getDouble("coefJz", par.coefJz);
	par.beta   = kvmap.getDouble("beta",   par.beta);
	par.Tsfr   = kvmap.getDouble("Tsfr",   par.Tsfr);  // dimensionless! in units of Hubble time
	par.sigmabirth = kvmap.getDouble("sigmabirth", par.sigmabirth);  // dimensionless ratio
	return par;
}*/

EXP ExponentialParam parseExponentialParam(
				       const utils::KeyValueMap& kvmap,
				       const units::ExternalUnits& conv)
{
	ExponentialParam par;
	par.norm   = kvmap.getDouble("norm",	par.norm)   * conv.massUnit;
	par.mass   = kvmap.getDouble("mass",	par.mass)   * conv.massUnit;
	par.Jr0    = kvmap.getDouble("Jr0",	par.Jr0)    * conv.lengthUnit * conv.velocityUnit;
	par.Jz0    = kvmap.getDouble("Jz0",	par.Jz0)    * conv.lengthUnit * conv.velocityUnit;
	par.Jphi0  = kvmap.getDouble("Jphi0",	par.Jphi0)  * conv.lengthUnit * conv.velocityUnit;
	par.pr     = kvmap.getDouble("pr",	par.pr);
	par.pz     = kvmap.getDouble("pz",	par.pz);
	par.addJden= kvmap.getDouble("addJden",	par.addJden)* conv.lengthUnit * conv.velocityUnit;
	par.addJvel= kvmap.getDouble("addJvel",	par.addJvel)* conv.lengthUnit * conv.velocityUnit;
	par.Fname = kvmap.getString("PopFile",	par.Fname);
	return par;
}

EXP taperExpParam parsetaperExpParam(
				  const utils::KeyValueMap& kvmap,
				  const units::ExternalUnits& conv)
{
	taperExpParam par;
	par.norm   = kvmap.getDouble("norm",	par.norm)   * conv.massUnit;
	par.mass   = kvmap.getDouble("mass",	par.mass)   * conv.massUnit;
	par.Jr0    = kvmap.getDouble("Jr0",	par.Jr0)    * conv.lengthUnit * conv.velocityUnit;
	par.Jz0    = kvmap.getDouble("Jz0",	par.Jz0)    * conv.lengthUnit * conv.velocityUnit;
	par.Jtaper = kvmap.getDouble("Jtaper",	par.Jtaper) * conv.lengthUnit * conv.velocityUnit;
	par.Jtrans = kvmap.getDouble("Jtrans",	par.Jtrans) * conv.lengthUnit * conv.velocityUnit;
	par.Jcut   = kvmap.getDouble("Jcut",	par.Jcut) * conv.lengthUnit * conv.velocityUnit;
	par.Delta  = kvmap.getDouble("Delta",	par.Delta) * conv.lengthUnit * conv.velocityUnit;
	par.Jphi0  = kvmap.getDouble("Jphi0",	par.Jphi0)  * conv.lengthUnit * conv.velocityUnit;
	par.pr	   = kvmap.getDouble("pr",	par.pr);
	par.pz	   = kvmap.getDouble("pz",	par.pz);
	par.addJden= kvmap.getDouble("addJden",	par.addJden)* conv.lengthUnit * conv.velocityUnit;
	par.addJvel= kvmap.getDouble("addJvel", par.addJvel)* conv.lengthUnit * conv.velocityUnit;
	par.Fname  = kvmap.getString("PopFile", par.Fname);
	return par;
}

EXP OxfordParam parseOxfordParams(
			      const utils::KeyValueMap& kvmap,
			      const units::ExternalUnits& conv)
{
	OxfordParam par;
	par.mass      = kvmap.getDouble("mass")    * conv.massUnit;
	par.J0        = kvmap.getDouble("J0")      * conv.lengthUnit * conv.velocityUnit;
	par.Jcutoff   = kvmap.getDouble("Jcutoff") * conv.lengthUnit * conv.velocityUnit;
	par.Jphi0     = kvmap.getDouble("Jphi0")   * conv.lengthUnit * conv.velocityUnit;
	par.Jcore     = kvmap.getDouble("Jcore")   * conv.lengthUnit * conv.velocityUnit;
	par.L0        = kvmap.getDouble("L0")      * conv.lengthUnit * conv.velocityUnit;
	par.gamma0    = kvmap.getDouble("gamma0",    par.gamma0);
	par.gamma1    = kvmap.getDouble("gamma1",    par.gamma1);
	par.bz        = kvmap.getDouble("bz",        par.bz);
	par.slopeIn   = kvmap.getDouble("slopeIn",   par.slopeIn);
	par.slopeOut  = kvmap.getDouble("slopeOut",  par.slopeOut);
	par.steepness = kvmap.getDouble("steepness", par.steepness);
	par.coefLin   = kvmap.getDouble("coefLin",   par.coefLin);
	par.kIn       = kvmap.getDouble("kIn",       par.kIn);
	par.coefLout  = kvmap.getDouble("coefLout",  par.coefLout);
	par.kOut      = kvmap.getDouble("kOut",      par.kOut);
	par.rotFrac   = kvmap.getDouble("rotFrac",   par.rotFrac);
	par.cutoffStrength = kvmap.getDouble("cutoffStrength", par.cutoffStrength);
	return par;
}

EXP BolognaParam parseBolognaParams(
			      const utils::KeyValueMap& kvmap,
			      const units::ExternalUnits& conv)
{
	BolognaParam par;
	par.mass      = kvmap.getDouble("mass")    * conv.massUnit;
	par.J0        = kvmap.getDouble("J0")      * conv.lengthUnit * conv.velocityUnit;
	par.Jphi      = kvmap.getDouble("Jphi")   * conv.lengthUnit * conv.velocityUnit;
	par.L0        = kvmap.getDouble("L0")   * conv.lengthUnit * conv.velocityUnit;
	par.gamma     = kvmap.getDouble("gamma", par.gamma);
	par.L1        = kvmap.getDouble("L1")   * conv.lengthUnit * conv.velocityUnit;
	par.alpha     = kvmap.getDouble("alpha", par.alpha);
	par.coefL     = kvmap.getDouble("coefL", par.coefL);
	par.kIn       = kvmap.getDouble("kIn", par.kIn);
	par.kOut      = kvmap.getDouble("kOut", par.kOut);
	par.rotFrac   = kvmap.getDouble("rotFrac", par.rotFrac);
	return par;
}

inline void checkNonzero(const potential::BasePotential* potential, const std::string& type)
{
    if(potential == NULL)
        throw std::invalid_argument("Need an instance of potential to initialize "+type+" DF");
}

EXP PtrDistributionFunction createDistributionFunction(
    const utils::KeyValueMap& kvmap,
    const potential::BasePotential* potential,
    const potential::BaseDensity* density,
    const units::ExternalUnits& converter)
{
    std::string type = kvmap.getString("type");
    // for some DF types, there are two alternative ways of specifying the normalization:
    // either directly as norm, Sigma0, etc., or as the total mass, from which the norm is computed
    // by creating a temporary instance of a corresponding DF class, and computing its mass
    double mass = kvmap.getDouble("mass", NAN)* converter.massUnit;
    if(utils::stringsEqual(type, "DoublePowerLaw")) {
	    DoublePowerLawParam par = parseDoublePowerLawParam(kvmap, converter);
	    if(mass>0) {
		    par.norm = 1.0;
		    par.norm = mass / DoublePowerLaw(par).totalMass();
	    }
	    return PtrDistributionFunction(new DoublePowerLaw(par));
    }
    if(utils::stringsEqual(type, "ModifiedDoublePowerLaw")) {
	    ModifiedDoublePowerLawParam par = parseModifiedDoublePowerLawParam(kvmap, converter);
	    if(mass>0) {
		    par.norm = 1.0;
		    par.norm = mass / ModifiedDoublePowerLaw(par).totalMass();
	    }
	    return PtrDistributionFunction(new ModifiedDoublePowerLaw(par));
    }
    if(utils::stringsEqual(type, "SinDoublePowerLaw")) {
	    SinDoublePowerLawParam par = parseSinDoublePowerLawParam(kvmap, converter);
	    if(mass>0) {
		    par.norm = 1.0;
		    par.norm = mass / SinDoublePowerLaw(par).totalMass();
	    }
	    return PtrDistributionFunction(new SinDoublePowerLaw(par));
    }
    else if(utils::stringsEqual(type, "Exponential")) {
	    ExponentialParam par = parseExponentialParam(kvmap, converter);
	    if(mass>0) {
		    par.norm = 1.0;
		    par.norm = mass / Exponential(par).totalMass();
	    }
	    return PtrDistributionFunction(new Exponential(par));
    }
    else if(utils::stringsEqual(type, "taperExp")) {
	    taperExpParam par = parsetaperExpParam(kvmap, converter);
	    if(mass>0) {
		    par.norm = 1.0;
		    par.norm = mass / taperExp(par).totalMass();
	    }
	    return PtrDistributionFunction(new taperExp(par));
    }
    else if(utils::stringsEqual(type, "QuasiIsothermal")) {
	    checkNonzero(potential, type);
	    potential::Interpolator pot_interp(*potential);
        QuasiIsothermalParam par = parseQuasiIsothermalParam(kvmap, converter);
        if(mass>0) {
            par.Sigma0 = 1.0;
            par.Sigma0 = mass / QuasiIsothermal(par, pot_interp).totalMass();
        }
        return PtrDistributionFunction(new QuasiIsothermal(par, pot_interp));
    }
    else if(utils::stringsEqual(type, "QuasiSpherical")) {
        checkNonzero(potential, type);
        if(density == NULL)
            density = potential;
        double beta0 = kvmap.getDoubleAlt("beta", "beta0", 0);
        double r_a   = kvmap.getDoubleAlt("anisotropyRadius", "r_a", INFINITY) * converter.lengthUnit;
        return PtrDistributionFunction(new QuasiSphericalCOM(
            potential::DensityWrapper(*density), potential::PotentialWrapper(*potential), beta0, r_a));
    }
    else if(utils::stringsEqual(type, "Oxford")) {
	    checkNonzero(potential, type);
	    potential::Interpolator pot_interp(*potential);
	    return PtrDistributionFunction(new Oxford(parseOxfordParams(kvmap, converter),pot_interp));
    }
    else if(utils::stringsEqual(type, "Bologna")) {
	    return PtrDistributionFunction(new Bologna(parseBolognaParams(kvmap, converter)));
    }

    else{
	    printf("Unknown type of distribution function: %s xx\n",type.c_str()); exit(1);
    }
}

}  // namespace df
