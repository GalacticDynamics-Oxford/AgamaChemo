/** \file    potential_galpot.h
    \brief   a reimplementation of Walter Dehnen's GalaxyPotential code
    \author  Eugene Vasiliev, based on the earlier work of Walter Dehnen, Paul McMillan
    \date    2015-2017

The original GalPot code is written by W.Dehnen:

Copyright Walter Dehnen, 1996-2005
e-mail:   walter.dehnen@astro.le.ac.uk
address:  Department of Physics and Astronomy, University of Leicester
          University Road, Leicester LE1 7RH, United Kingdom

Put into the Torus code (with a minimum of fuss) by Paul McMillan, Oxford 2010
email: p.mcmillan1@physics.ox.ac.uk

Modifications by Eugene Vasiliev, 2015-2017:
The present code is a completely redesigned implementation,
with only a few class names retaining the traces of the original code.


The method, explained in Dehnen & Binney (1998, MNRAS, 294, 429) and based
on the approach of Kuijken & Dubinski (1994, MNRAS, 269, 13), is applicable
to any disk density profile which is separable in cylindrical coordinates.

Let the density profile of the disk be

\f$  \rho_d(R,z) = f(R) h(z)  \f$,

and let H(z) be the second antiderivative of h(z), so that d^2 H / dz^2 = h.
Then the potential of the disk can be written as a sum of 'main' and 'residual' parts:

\f$  \Phi(R,z) = 4\pi f(r) H(z) + \Phi_{res}  \f$,

where the argument of f is spherical rather than cylindrical radius,
and the residual potential is generated by the following density profile:

\f$  \rho_{res} = [f(R)-f(r)] h(z) - f''(r) H(z) - 2 f'(r) [H(z) + z H'(z)]/r  \f$.

This residual potential is not strongly confined to the disk plane, and can be
efficiently approximated by a multipole expanion, which, in turn, is represented
by a two-dimensional quintic spline in (R,z) plane.

The original GalaxyPotential uses this method for any combination of disk components
and additional, possibly flattened spheroidal components: the residual density of all
disks and the entire density of spheroids serves as the source to the Multipole potential
approximation.

In the present modification, the GalaxyPotential class is replaced by a more generic Composite
potential, which contains one Multipole potential and possibly several DiskAnsatz components.
The latter come in pairs with DiskDensity density components, so that the difference between
the full input density and the one provided by DiskAnsatz is used in the multipole expansion.
A composite density model with all SpheroidDensity components and all pairs of DiskDensity minus
DiskAnsatz components is used to initialize the Multipole potential.
Of course this input may be generalized to contain other density components, and the Composite
potential may also contain some other potential models apart from DiskAnsatz and Multipole.

The Multipole potential solves the Poisson equation using the spherical-harmonic expansion
of its input density profile, and then stores the values and derivatives of potential on
a 2d grid in (r,theta) plane, so that the potential evaluation uses 2d spline-interpolated
values; however, if the radius lies outside the grid definition region, the potential is computed
by summing up appropriately extrapolated multipole components (unlike the original GalPot).

For compatibility with the original implementation, an utility function `readGalaxyPotential`
is provided in potential_factory.h, taking the name of parameter file and the Units object as parameters.
*/

#pragma once
#include "potential_base.h"
#include "smart.h"
#include <vector>
#define EXP __declspec(dllexport)

namespace potential{

/// \name  Separable disk density profile
///@{

/** Parameters that describe a disk component.

    Specification of a disk density profile separable in R and z requires two auxiliary function,
    f(R) and H(z)  (the former essentially describes the surface density of the disk,
    and the latter is the second antiderivative of vertical density profile h(z) ).
    They are used by both DiskAnsatz potential and DiskDensity density classes.
    In the present implementation they are almost the same as in GalPot,
    generalized for an arbitrary Sersic index 'n' of the surface density profile:

    \f$  \rho = f(R) h(z)  \f$,

    \f$  f(R) = \Sigma_0  \exp [ -(R/R_d)^{1/n} - R_0/R + \epsilon \cos(R/R_d) ]  \f$,

    \f$  h(z) = \delta(z)                 \f$  for  h=0, or 
    \f$  h(z) = 1/(2 h)  * exp(-|z/h|)    \f$  for  h>0, or
    \f$  h(z) = 1/(4|h|) * sech^2(|z/2h|) \f$  for  h<0.

    The corresponding second antiderivatives of h(z) are given in Table 2 of Dehnen&Binney 1998.
    Alternatively, one may provide two arbitrary 1d functions to be used in the separable profile.
*/
struct EXP DiskParam{
    double surfaceDensity;      ///< surface density normalisation Sigma_0
    double scaleRadius;         ///< scale length R_d (NOT the effective radius!)
    double scaleHeight;         ///< scale height h: 
    ///< For h<0 an isothermal (sech^2) profile is used, for h>0 an exponential one, 
    ///< and for h=0 the disk is infinitesimally thin
    double innerCutoffRadius;   ///< if nonzero, specifies the radius of a hole at the center R_0
    double modulationAmplitude; ///< a term eps*cos(R/R_d) is added to the radial exponent
    double sersicIndex;         ///< Sersic index n (by default n=1, corresponding to an exponential disk)
    /// set up default values for all parameters
    DiskParam(
        double _surfaceDensity=0,
        double _scaleRadius=1,
        double _scaleHeight=1,
        double _innerCutoffRadius=0,
        double _modulationAmplitude=0,
        double _sersicIndex=1)
    :
        surfaceDensity(_surfaceDensity), scaleRadius(_scaleRadius), scaleHeight(_scaleHeight),
        innerCutoffRadius(_innerCutoffRadius), modulationAmplitude(_modulationAmplitude),
        sersicIndex(_sersicIndex)
    {}
    double mass() const;        ///< return the total mass of a density profile with these parameters
};

/** helper routine to create an instance of radial density function */
EXP math::PtrFunction createRadialDiskFnc(const DiskParam& params);

/** helper routine to create an instance of vertical density function */
EXP math::PtrFunction createVerticalDiskFnc(const DiskParam& params);

/** Density profile of a separable disk model */
class EXP DiskDensity: public BaseDensity {
public:
    /// construct the density profile with provided parameters
    DiskDensity(const DiskParam& _params) : 
        radialFnc  (createRadialDiskFnc(_params)),
        verticalFnc(createVerticalDiskFnc(_params)) {};

    /// construct a generic profile with user-specified radial and vertical functions
    DiskDensity(const math::PtrFunction& _radialFnc, const math::PtrFunction& _verticalFnc) :
        radialFnc(_radialFnc), verticalFnc(_verticalFnc) {}

    virtual coord::SymmetryType symmetry() const { return coord::ST_AXISYMMETRIC; }
    virtual const char* name() const { return myName(); }
    static const char* myName() { static const char* text = "Disk"; return text; }
private:
    math::PtrFunction radialFnc;     ///< function describing radial dependence of surface density
    math::PtrFunction verticalFnc;   ///< function describing vertical density profile
    virtual double densityCyl(const coord::PosCyl &pos) const;
    virtual double densityCar(const coord::PosCar &pos) const
    {  return densityCyl(toPosCyl(pos)); }
    virtual double densitySph(const coord::PosSph &pos) const
    {  return densityCyl(toPosCyl(pos)); }
};

/** Part of the disk potential provided analytically as  4 pi f(r) H(z).
    This potential model is intended for the internal usage in the GalPot approach,
    not as an independent potential; its total mass is zero. */
class EXP DiskAnsatz: public BasePotentialCyl {
public:
    DiskAnsatz(const DiskParam& _params) : 
        radialFnc  (createRadialDiskFnc(_params)),
        verticalFnc(createVerticalDiskFnc(_params)) {};
    DiskAnsatz(const math::PtrFunction& _radialFnc, const math::PtrFunction& _verticalFnc) :
        radialFnc(_radialFnc), verticalFnc(_verticalFnc) {};
    virtual coord::SymmetryType symmetry() const { return coord::ST_AXISYMMETRIC; }
    virtual const char* name() const { return myName(); }
    static const char* myName() { static const char* text = "DiskAnsatz"; return text; }
    virtual double totalMass() const { return 0; }  // all the mass is contained in the residual density
private:
    math::PtrFunction radialFnc;     ///< function describing radial dependence of surface density
    math::PtrFunction verticalFnc;   ///< function describing vertical density profile
    /** Compute _part_ of the disk potential: f(r)*H(z) */
    virtual void evalCyl(const coord::PosCyl &pos,
        double* potential, coord::GradCyl* deriv, coord::HessCyl* deriv2) const;
    virtual double densityCyl(const coord::PosCyl &pos) const;
};

///@}
/// \name  Spheroical density profile
///@{

/** Parameters describing a double-power-law spheroidal density profile.
    The density is given by the Zhao(1996) alpha-beta-gamma model multiplied by an optional
    exponential cutoff, and includes many popular profiles (e.g., Dehnen, Prugniel-Simien, Gaussian):
    \f$  \rho = \rho_0  (r/r_0)^{-\gamma} ( 1 + (r/r_0)^\alpha )^{(\gamma-\beta) / \alpha}
    \exp[ -(r/r_{cut})^\xi], \f$.
*/
struct EXP SpheroidParam{
    double densityNorm;         ///< density normalization rho_0
    double axisRatioY;          ///< axis ratio p (y/x)
    double axisRatioZ;          ///< axis ratio q (z/x)
    double alpha;               ///< steepness of transition alpha
    double beta;                ///< outer power slope beta
    double gamma;               ///< inner power slope gamma
    double scaleRadius;         ///< transition radius r_0
    double outerCutoffRadius;   ///< outer cut-off radius r_{cut}
    double cutoffStrength;      ///< steepness of the exponential cutoff xi
    /// set up default values for all parameters
    SpheroidParam(
        double _densityNorm=0,
        double _axisRatioY=1,
        double _axisRatioZ=1,
        double _alpha=1,
        double _beta=4,
        double _gamma=1,
        double _scaleRadius=1,
        double _outerCutoffRadius=0,
        double _cutoffStrength=2)
    :
        densityNorm(_densityNorm), axisRatioY(_axisRatioY), axisRatioZ(_axisRatioZ),
        alpha(_alpha), beta(_beta), gamma(_gamma),
        scaleRadius(_scaleRadius), outerCutoffRadius(_outerCutoffRadius), cutoffStrength(_cutoffStrength)
    {}
    double mass() const;        ///< return the total mass of a density profile with these parameters
    static const char* myName() ///< the name of the density model in potential_factory routines
    { static const char* text = "Spheroid"; return text; }

};

/** Parameters describing a Sersic density profile.
    In the spherical case, the projected density is given by
    \f$  \Sigma(R) = \Sigma_0  \exp[ -b (R/R_e)^{1/n} ]  \f$,
    where Sigma_0 is the central surface density, R_e is the effective radius,
    n is the shape parameter (Sersic index), and b is the internally computed numerical constant,
    approximately equal to 2n - 1/3.
    The 3d density profile is obtained by deprojecting this expression and interpolating the result.
    In non-spherical cases, the flattening with constant axis ratios is applied to the 3d density,
    so that equidensity surfaces are concentric ellipsoids.
    For consistency with other density models, we use the central value of surface density,
    not the one at the effective radius (they are related by \f$  \Sigma_e = \Sigma_0 \exp(-b)  \f$);
    moreover, it is possible to provide the total mass rather than the surface density
    when constructing this model via the `createDensity()` routine.
*/
struct EXP SersicParam{
    double surfaceDensity; ///< central surface density Sigma_0
    double scaleRadius;    ///< effective radius containing 1/2 of the total mass in projection
    double sersicIndex;    ///< shape parameter `n` (Sersic index), should be positive
    double axisRatioY;     ///< axis ratio p (Y/X)
    double axisRatioZ;     ///< axis ratio q (Z/X)
    /// set up default values for all parameters
    SersicParam(
        double _surfaceDensity=0,
        double _scaleRadius=1,
        double _sersicIndex=4,
        double _axisRatioY=1,
        double _axisRatioZ=1)
    :
        surfaceDensity(_surfaceDensity), scaleRadius(_scaleRadius), sersicIndex(_sersicIndex),
        axisRatioY(_axisRatioY), axisRatioZ(_axisRatioZ)
    {}
    double b() const;      ///< compute the numerical coefficient b as a function of n
    double mass() const;   ///< return the total mass of a density profile with these parameters
    static const char* myName() ///< the name of the density model in potential_factory routines
    { static const char* text = "Sersic"; return text; }
};

/// helper routine to construct a one-dimensional function describing a double-power-law profile
EXP math::PtrFunction createSpheroidDensity(const SpheroidParam& params);

/// helper routine to construct a one-dimensional function describing a Sersic profile
EXP math::PtrFunction createSersicDensity(const SersicParam& params);

/** Density profile described by an arbitrary function f of ellipsoidal radius:
    \f$  \rho(x,y,z) = f(\tilde r),  \tilde r = \sqrt{x^2 + (y/p)^2 + (z/q)^2} \f$,
    where p = y/x and q = z/x  are two axis ratios.
*/
class EXP SpheroidDensity: public BaseDensity{
public:
    /// construct a generic profile with a user-specified one-dimensional function and two axis ratios
    SpheroidDensity(const math::PtrFunction& fnc, const double axisRatioY=1, const double axisRatioZ=1) :
        p2(pow_2(axisRatioY)), q2(pow_2(axisRatioZ)), rho(fnc) {}

    /// construct a model with the provided Spheroid density parameters (convenience overload)
    SpheroidDensity(const SpheroidParam& params) :
        p2(pow_2(params.axisRatioY)), q2(pow_2(params.axisRatioZ)),
        rho(createSpheroidDensity(params)) {}

    /// construct a model with the provided Sersic density parameters (convenience overload)
    SpheroidDensity(const SersicParam& params) :
        p2(pow_2(params.axisRatioY)), q2(pow_2(params.axisRatioZ)),
        rho(createSersicDensity(params)) {}

    virtual coord::SymmetryType symmetry() const {
        return p2==1 ? (q2==1 ? coord::ST_SPHERICAL : coord::ST_AXISYMMETRIC) : coord::ST_TRIAXIAL; }
    virtual const char* name() const { return myName(); }
    /// the name reported by the density model is always 'Spheroid',
    /// regardless of whether it was initialized from SersicParam or SpheroidParam
    static const char* myName() { return SpheroidParam::myName(); }
private:
    const double p2, q2;    ///< squared axis ratios p=y/x, q=z/x
    math::PtrFunction rho;  ///< one-dimensional density as a function of elliptical radius
    virtual double densityCar(const coord::PosCar &pos) const;
    virtual double densityCyl(const coord::PosCyl &pos) const;
    virtual double densitySph(const coord::PosSph &pos) const;
};

///@}

} // namespace potential
