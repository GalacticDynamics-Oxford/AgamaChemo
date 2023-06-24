/** \file   math_specfunc.h
    \brief  various special functions
    \date   2015-2016
    \author Eugene Vasiliev
*/
#pragma once
#define EXP __declspec(dllexport)

// most of the functions here return NAN in case of invalid arguments or out-of-range results;
// the error message is stored in the string variable math::exceptionText
namespace math {

/** Gegenbauer (ultraspherical) polynomial:  \f$ C_n^{(\lambda)}(x) \f$ */
EXP double gegenbauer(const int n, double lambda, double x);

/** Array of Gegenbauer (ultraspherical) polynomials for n=0,1,...,nmax */
EXP void gegenbauerArray(const int nmax, double lambda, double x, double* result);

/** Inverse error function (defined for -1<x<1) */
EXP double erfinv(const double x);

/** Gauss's hypergeometric function 2F1(a, b; c; x) */
EXP double hypergeom2F1(const double a, const double b, const double c, const double x);

/** Associate Legendre function of the second kind, together with its derivative if necessary */
EXP double legendreQ(const double m, const double x, double* deriv=0/*NULL*/);

/** Factorial of an integer number */
EXP double factorial(const unsigned int n);

/** Logarithm of factorial of an integer number (doesn't overflow quite that easy) */
EXP double lnfactorial(const unsigned int n);

/** Double-factorial n!! of an integer number */
EXP double dfactorial(const unsigned int n);

/** Gamma function */
EXP double gamma(const double x);

/** Logarithm of gamma function (doesn't overflow quite that easy) */
EXP double lngamma(const double n);

/** Upper incomplete Gamma function: the ordinary Gamma(x) = gammainc(x, y=0) */
EXP double gammainc(const double x, const double y);

/** Psi (digamma) function */
EXP double digamma(const double x);

/** Psi (digamma) function for integer argument */
EXP double digamma(const int x);

/** Complete elliptic integrals of the first kind K(k) = F(pi/2, k);
    in this and other elliptic integrals, default mode is single-precision accuracy,
    to obtain higher accuracy for a larger cost, pass true as the last argument */
EXP double ellintK(const double k, bool accurate=false);

/** Complete elliptic integrals of the second kind K(k) = E(pi/2, k) */
EXP double ellintE(const double k, bool accurate=false);

/** Incomplete elliptic integrals of the first kind:
    \f$  F(\phi,k) = \int_0^\phi d t \, \frac{1}{\sqrt{1-k^2\sin^2 t}}  \f$  */
EXP double ellintF(const double phi, const double k, bool accurate=false);

/** Incomplete elliptic integrals of the second kind:
    \f$  E(\phi,k) = \int_0^\phi d t \, \sqrt{1-k^2\sin^2 t}  \f$  */
EXP double ellintE(const double phi, const double k, bool accurate=false);

/** Incomplete elliptic integrals of the third kind:
    \f$  \Pi(\phi,k,n) = \int_0^\phi d t \, \frac{1}{(1+n\sin^2 t)\sqrt{1-k^2\sin^2 t}}  \f$  */
EXP double ellintP(const double phi, const double k, const double n, bool accurate=false);

/** Bessel J_n(x) function (regular at x=0) */
EXP double besselJ(const int n, const double x);

/** Bessel Y_n(x) function (singular at x=0) */
EXP double besselY(const int n, const double x);

/** Modified Bessel function I_n(x) (regular) */
EXP double besselI(const int n, const double x);

/** Modified Bessel function K_n(x) (singular) */
EXP double besselK(const int n, const double x);

/** Lambert W function W_0,W_{-1}(x), satisfying the equation W exp(W) = x;
    for -exp(-1)<=x<0 there are two branches, W_0>=-1, and W_{-1}<=-1;
    the second one is selected by setting the second argument to true */
EXP double lambertW(const double x, bool Wminus1branch=false);

/** solve the Kepler equation:  phase = eta - ecc * sin(eta)  for eta (eccentric anomaly) */
EXP double solveKepler(double ecc, double phase);

}  // namespace
