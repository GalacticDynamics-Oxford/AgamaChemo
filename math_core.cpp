#include "math_core.h"
#include "math_glquadrature.h"
#include "utils.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_version.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_fft_real.h>
#include <stdexcept>
#include <cassert>
#include <cmath>
#include <iostream>
#include "alloca.h"

#if not defined(GSL_MAJOR_VERSION) || (GSL_MAJOR_VERSION == 1) && (GSL_MINOR_VERSION < 15)
#error "GSL version is too old (need at least 1.15)"
#endif

#ifdef HAVE_CUBA
#include <cuba.h>
#else
#include "cubature.h"
#endif

namespace math{

/// upper limit on the number of iterations in root-finders, minimizers, etc.
static const int MAXITER = 64;

/// size of workspace for adaptive integration
static const int MAX_INTEGR_POINTS = 1000;

// ------ error handling ------ //

/// stores the exception text to be propagated through external C code that does not support exceptions
std::string exceptionText;

#define CALL_FUNCTION_OR_THROW(x) \
    exceptionText.clear(); \
    x; \
    if(!exceptionText.empty()) throw std::runtime_error(exceptionText);

/// callback function invoked by GSL in case of error; stores the error text in a global variable
/// (not thread-safe! assumed that these events don't occur often)
EXP void GSLerrorHandler(const char *reason, const char* file, int line, int gsl_errno)
{
    if( // list error codes that are non-critical and don't need to be reported
        gsl_errno == GSL_ETOL ||
        gsl_errno == GSL_EROUND ||
        gsl_errno == GSL_ESING ||
        gsl_errno == GSL_EDIVERGE )
        return;  // do nothing
    if(exceptionText.empty())
        exceptionText = (
        gsl_errno == GSL_ERANGE || gsl_errno == GSL_EOVRFLW ? "GSL range error" :
        gsl_errno == GSL_EDOM ? "GSL domain error" :
        gsl_errno == GSL_EINVAL ? "GSL invalid argument error" :
        "GSL error " + utils::toString(gsl_errno) ) +
        " in " + file + ", line " + utils::toString(line) + ": " + reason + "\n" + utils::stacktrace();
}

// a hacky way to initialize our error handler on module startup
EXP bool gsl_error_handler_set = gsl_set_error_handler(&GSLerrorHandler);

// ------ math primitives -------- //

/// wrapper for a user-defined function to provide to GSL
/// (any exception arising in the C++ code should not propagate back to the C code in GSL,
/// and are instead stored in a global variable, which is unfortunately not thread-safe)
EXP double functionWrapper(double x, void* param)
{
    try{
        return static_cast<IFunction*>(param)->value(x);
    }
    catch(std::exception& e){
	    std::cout << e.what() << "\n";
        return NAN;
    }
}

EXP int fcmp(double x, double y, double eps)
{
    if(x==0)
        return y<-eps ? -1 : y>eps ? +1 : 0;
    if(y==0)
        return x<-eps ? -1 : x>eps ? +1 : 0;
    if(x!=x)
        return -2;
    if(y!=y)
        return +2;
    return gsl_fcmp(x, y, eps);
}

EXP double pow(double x, int n)
{
    if(n<0) {
        n = -n;
        x = 1/x;
    }
    double result = 1;
    do {
        if(n%2) result *= x;
        n >>= 1;
        x *= x;
    } while(n);
    return result;
}

EXP double pow(double x, double n)
{
    if(n == 0.0) return 1;
    if(n == 1.0) return x;
    if(n ==-1.0) return 1/x;
    if(n == 2.0) return x*x;
    if(n ==-2.0) return 1/(x*x);
    if(n == 0.5) return sqrt(x);
    if(n ==-0.5) return 1/sqrt(x);
    if(n == 3.0) return x*x*x;
    if(n ==-3.0) return 1/(x*x*x);
    return std::pow(x, n);
}

EXP double wrapAngle(double x)
{
    // make no attempt to "accurately" subtract the exact value of 2pi
    // (not representable as a finite-mantissa floating-point number).
    // on the contrary, take the modulus w.r.t the floating-point approximation of 2pi,
    // so that 2*M_PI, 4*M_PI, -2*M_PI, etc., produce exactly 0
    x -= 2*M_PI * floor(1/(2*M_PI) * x);
    return x>=0 && x<2*M_PI ? x : 0;
    // x==2*M_PI may occur only if x<0 and |x| is very small, so it is rounded to zero
}

EXP double unwrapAngle(double x, double xprev)
{
    double diff=(x-xprev)/(2*M_PI);
    double nwraps=0;
    if(diff>0.5) 
        modf(diff+0.5, &nwraps);
    else if(diff<-0.5) 
        modf(diff-0.5, &nwraps);
    return x - 2*M_PI * nwraps;
}

EXP void sincos(double x, double& s, double& c)
{
    double y = x>=0 ? x : -x;     // fabs(x)
    double z = x>=0 ? 1 : -1;     // sign(x)
    long quad = long(4/M_PI * y); // floor(...), non-negative (!!no overflow check!!)
    quad = (quad+1) >> 1;         // 0 => 0, 1 => 1, 2 => 1, 3 => 2, 4 => 2, 5 => 3, etc.
    y -= M_PI/2 * quad;           // bring the range to [-pi/4 .. pi/4];
    // note that multiples of M_PI/2 are exactly mapped to zero (this is a deliberate tweak).
    double y2 = y * y;
    // use a Chebyshev approximation for sin and cos on this interval
    double sy = y + y * (((((
        +1.5896230157654657e-10 * y2
        -2.5050747762857807e-8) * y2
        +2.7557313621385725e-6) * y2
        -1.9841269829589539e-4) * y2
        +8.3333333333221186e-3) * y2
        -0.1666666666666663073) * y2;
    double cy = 1.0 + (((((
        +2.064357075039214e-09  * y2
        -2.755549453909573e-07) * y2
        +2.480158051577225e-05) * y2
        -0.0013888888877062660) * y2
        +0.0416666666665909000) * y2
        -0.5000000000000000000) * y2;
    // assign the output values, depending on the quadrant (also change the sign of sin(x) if x<0)
    switch(quad & 3) {
        case 0: s = sy * z; c = cy; return;
        case 1: s = cy * z; c =-sy; return;
        case 2: s =-sy * z; c =-cy; return;
        case 3: s =-cy * z; c = sy; return;
        default: s=NAN; c=NAN; /*shouldn't occur*/
    }
}


template<typename NumT>
ptrdiff_t binSearch(const NumT x, const NumT arr[], size_t size)
{
    if(size<1 || !(x>=arr[0]))
        return -1;
    if(x>arr[size-1] || size<2)
        return size-1;
    // first guess the likely location in the case that the input grid is equally-spaced
    ptrdiff_t index = static_cast<ptrdiff_t>( (x-arr[0]) / (arr[size-1]-arr[0]) * (size-1) );
    ptrdiff_t indhi = size-1;
    if(index==static_cast<ptrdiff_t>(size)-1)
        return size-2;     // special case -- we are exactly at the end of array, return the previous node
    if(x>=arr[index]) {
        if(x<arr[index+1])
            return index;  // guess correct, exiting
        // otherwise the search is restricted to [ index .. indhi ]
    } else {
        indhi = index;     // search restricted to [ 0 .. index ]
        index = 0;
    }
    // this will always end up with one grid node in O(log(N)) steps,
    // even if the grid nodes were not monotonic (we don't check this assertion to avoid wasting time)
    while(indhi > index + 1) {
        ptrdiff_t i = (indhi + index)/2;
        if(arr[i] > x)
            indhi = i;
        else
            index = i;
    }
    return index;
}

// template instantiations
template ptrdiff_t binSearch(const double,    const double[],    size_t);
template ptrdiff_t binSearch(const float,     const float[],     size_t);
template ptrdiff_t binSearch(const int,       const int[],       size_t);
template ptrdiff_t binSearch(const long,      const long[],      size_t);
template ptrdiff_t binSearch(const long long, const long long[], size_t);
template ptrdiff_t binSearch(const unsigned int,       const unsigned int[],       size_t);
template ptrdiff_t binSearch(const unsigned long,      const unsigned long[],      size_t);
template ptrdiff_t binSearch(const unsigned long long, const unsigned long long[], size_t);


/* ------ algebraic transformations of functions ------- */

// u in (-inf,inf) but not too large (e.g. a logarithm of something else)
template<> EXP double scale(const ScalingInf& /*scaling*/, double u) {
    return  fabs(u) < 1 ? // two cases depending on whether |u| is small or large
        1 / (1 + sqrt(1 + 0.25 * u*u) - 0.5 * u) :      // u is close to zero
        0.5 + sqrt(0.25 + pow_2(1/u)) * sign(u) - 1/u;  // u is large or even infinite
}

template<> EXP double unscale(const ScalingInf& /*scaling*/, double s, double* duds) {
    if(duds)
        *duds = 1 / pow_2(1-s) + 1 / pow_2(s);
    return 1 / (1-s) - 1 / s;
}

// u in (-inf, u0] or [u0, +inf), when u0 is +-zero or the sign of u0 is the same as the sign of infinity
template<> EXP double scale(const ScalingSemiInf& scaling, double u) {
    if(scaling.u0 == 0)  // transform u to log(u) and then use the scaling on a doubly-infinite interval
        return scale(ScalingInf(), std::signbit(scaling.u0) ? -log(-u) : log(u));
    double l = log(u / scaling.u0);  // expected to be >=0, but this must be ensured by the calling code
    return l / (1+l);
}

template<> EXP double unscale(const ScalingSemiInf& scaling, double s, double* duds) {
    if(scaling.u0 == 0) {
        double u = exp( 1 / (1-s) - 1 / s );
        if(std::signbit(scaling.u0))
            u = -1/u;
        if(duds)
            *duds = fabs(u) * (1 / pow_2(1-s) + 1 / pow_2(s));
        return u;
    } else {
        double u = scaling.u0 * exp( s / (1-s) );
        if(duds)
            *duds = u / pow_2(1-s);
        return u;
    }
}

// u in [uleft, uright], trivial linear transformation
template<> EXP double scale(const ScalingLin& scaling, double u) {
    return (u - scaling.uleft) / (scaling.uright - scaling.uleft);
}

template<> EXP double unscale(const ScalingLin& scaling, double s, double* duds) {
    if(duds)
        *duds = scaling.uright - scaling.uleft;
    return scaling.uleft * (1-s) + scaling.uright * s;
}
// u in [uleft, uright], trig transformation 
template<> EXP double scale(const ScalingCos2& scaling, double u) {
	if(u == scaling.uleft)  return 0;
	if(u == scaling.uright) return 1;
	return 2/M_PI*acos(sqrt((scaling.uright-u)/(scaling.uright-scaling.uleft)));
}

template<> EXP double unscale(const ScalingCos2& scaling, double s, double* duds) {
	if(duds)
		*duds = .5*M_PI*sin(M_PI*s)*(scaling.uright - scaling.uleft);
	return scaling.uright-(scaling.uright-scaling.uleft)*pow_2(cos(.5*M_PI*s));
}

template<> EXP double scale(const ScalingCos2Inv& scaling, double u){
	double s = 1/(u-scaling.Delta2);
	return 2/M_PI*acos(sqrt((scaling.sright-s)/(2*scaling.Ds)));
}

template<> EXP double unscale(const ScalingCos2Inv& scaling, double s, double* duds){
	double t=scaling.sright-2*scaling.Ds*pow_2(cos(.5*M_PI*s));
	if(duds)
		*duds = -M_PI*sin(M_PI*s)*scaling.Ds/(t*t);
	return scaling.Delta2 + 1/t;
}

template<> EXP double scale(const ScalingSin2Inv& scaling, double u){
	double s = 1/(u-scaling.Delta2);
	double t = (scaling.sbar-s)/scaling.Ds;
	return acos((scaling.sbar - s)/scaling.Ds)/M_PI;
}

template<> EXP double unscale(const ScalingSin2Inv& scaling, double s, double* duds){
	double t = scaling.sbar - scaling.Ds*cos(s*M_PI);
	if(duds)
		*duds = -M_PI*scaling.Ds*sin(s*M_PI)/(t*t);
	return scaling.Delta2 + 1/t;
}

template<> EXP double scale(const ScalingSin2& scaling, double u) {
	if(u <= scaling.uleft)  return 0;
	if(u >= scaling.uright) return 1;
	return acos((scaling.ubar - u)/scaling.Du)/M_PI;
}

template<> EXP double unscale(const ScalingSin2& scaling, double s, double* duds) {
	if(duds)
		*duds = M_PI*scaling.Du*sin(s*M_PI);
	return scaling.ubar-scaling.Du*cos(s*M_PI);
}

// u in [uleft, uright], cubic transformation 
template<> EXP double scale(const ScalingCub& scaling, double u) {
    if(u == scaling.uleft)  return 0;
    if(u == scaling.uright) return 1;
    double  half = 0.5 * (scaling.uleft + scaling.uright);
    if(u == half) return 0.5;
    // rescale the input value into the range [0..0.5],
    // considering only the lower half of the symmetric transformation
    double w = u < half ?
        (u - scaling.uleft)  / (scaling.uright - scaling.uleft) :  // choose the more accurate expression
        (scaling.uright - u) / (scaling.uright - scaling.uleft);   // depending on the lower/upper half
    // initial guess, accurate to <0.1%
    double s = 1/M_SQRT3 * sqrt(w) + w * (0.1184 + 0.1291*w);
    // two iterations of the Newton method to reach almost machine accuracy (faster than analytic solution)
    for(int i=0; i<2; i++)
        s -= (1./6) * (s*s * (3-2*s) - w) / (s * (1-s));
    // if the rescaled input value was in the upper half of the unit interval, do the same for the output
    return u < half ? s : 1-s;
}

template<> EXP double unscale(const ScalingCub& scaling, double s, double* duds) {
    if(duds)
        *duds = (scaling.uright - scaling.uleft) * 6 * s * (1-s);
    return scaling.uleft * pow_2(1-s) * (1 + 2*s) + scaling.uright * pow_2(s) * (3 - 2*s);
}

namespace{
/// fast approximation for cubic root of a number 0<x<1, accurate to better than 1%
inline double fastcbrt(double x)
{
    int ex;
    double m = frexp(x, &ex);  // get mantissa and the exponent (assumed to be <= 0)
    switch(ex % 3) {           // linearly rescale the mantissa which lies in the range [0.5:1]
        case -2: m = 0.374 + 0.260 * m; break;  // multiplied by 0.5^(-2/3)
        case -1: m = 0.471 + 0.327 * m; break;  // multiplied by 0.5^(-1/3)
        default: m = 0.593 + 0.413 * m; break;  // no further multiplication
    }
    return ldexp(m, ex / 3);   // m * 2 ^ (ex [integer_divide_by] 3) 
}
}

// u in [uleft, uright], quintic transformation
template<> EXP double scale(const ScalingQui& scaling, double u) {
    // rescale the input value into the range [0..1]
    double v = (u - scaling.uleft) / (scaling.uright - scaling.uleft);
    // consider only the lower half of the symmetric transformation
    double w = v<=0.5 ? v : 1-v;
    // initial guess, accurate to ~1%
    double s = fastcbrt( w * (sqrt(w * (1.54 * w + 0.067)) - w + 0.1) );
    // two iterations of the Halley method to reach machine accuracy
    for(int i=0; i<2; i++) {
        double s2 = s*s, s3 = s2*s;
        s *= (s3 * (12*s3 - 33*s2 + 30*s - 10) + (3*s-2) * w) /
             (s3 * (18*s3 - 54*s2 + 55*s - 20) + (2*s-1) * w);
    }
    // if the rescaled input value was in the upper half of the unit interval, do the same for the output
    return v>0.5 ? 1-s : s;
}

template<> EXP double unscale(const ScalingQui& scaling, double s, double* duds) {
    if(duds)
        *duds = (scaling.uright - scaling.uleft) * 30 * pow_2(s * (1-s));
    return scaling.uleft * pow_3(1-s) * (1  + (6*s + 3 ) * s)
        + scaling.uright * pow_3( s ) * (10 + (6*s - 15) * s);
}


EXP void FncProduct::evalDeriv(const double x, double *val, double *der, double *der2) const
{
    double v1, v2, d1, d2, dd1, dd2;
    bool needDer = der!=NULL || der2!=NULL, needDer2 = der2!=NULL;
    f1.evalDeriv(x, &v1, needDer ? &d1 : 0, needDer2 ? &dd1 : 0);
    f2.evalDeriv(x, &v2, needDer ? &d2 : 0, needDer2 ? &dd2 : 0);
    if(val)
        *val = v1 * v2;
    if(der)
        *der = v1 * d2 + v2 * d1;
    if(der2)
        *der2 = v1 * dd2 + 2 * d1 * d2 + v2 * dd1;
}

EXP void LogLogScaledFnc::evalDeriv(const double logx,
    /*output*/ double* logf, double* der, double* der2) const
{
    double x = exp(logx), fval, fder=0;
    fnc.evalDeriv(x, &fval, der || der2 ? &fder : NULL, der2);
    double logder = fder * x / fval;  // logarithmic derivative d[ln(f)] / d[ln(x)]
    if(logf)
        *logf = log(fval);
    if(der)
        *der  = logder;
    if(der2) {
        fder *= logder-1;
        // if fder is very close to zero, der2 = (logder-1) * fder / x  may be rounded to zero,
        // in which case we won't be able to compute it correctly, so should return zero
        if(fabs(fder) < DBL_MIN || fabs(*der2) < DBL_MIN)
            *der2 = 0;
        else
            *der2 = (*der2 * x - fder) * x / fval;
    }
}


EXP void Monomial::evalDeriv(const double x, double *val, double *der, double *der2) const
{
    if(val)
        *val = pow(x, m);
    if(der)
        *der = m==0 ? 0 : m * pow(x, m-1);
    if(der2)
        *der2 = m==0 || m==1 ? 0 : m * (m-1) * pow(x, m-2);
}

EXP double Monomial::integrate(double x1, double x2, int n) const
{
    return m+n+1==0 ? log(x2/x1) : (pow(x2, m+n+1) - pow(x1, m+n+1)) / (m+n+1);
}

EXP void Gaussian::evalDeriv(const double x, double *val, double *der, double *der2) const
{
    double v = 1/M_SQRT2/M_SQRTPI / sigma * exp(-0.5 * pow_2(x / sigma));  // easy, yeah?
    if(val)
        *val = v;
    if(der)
        *der = -x / pow_2(sigma) * v;
    if(der2)
        *der2 = (pow_2(x / sigma) - 1) / pow_2(sigma) * v;
}

EXP double Gaussian::integrate(double x1, double x2, int n) const
{
    if(n<0)
        return NAN;  // undefined
    if(sigma==0)
        return (n==0 && x1<=0 && x2>=0) ? 1 : (n==0 && x1>=0 && x2<=0) ? -1 : 0;
    if(x1+x2==0 && n%2==1)
        return 0;  // avoid roundoff errors if the integral is symmetrically zero
    // not easy at all...
    double y1 = x1 / sigma, y2 = x2 / sigma;  // first, normalize the input
    double e1 = exp(-0.5*y1*y1), e2 = exp(-0.5*y2*y2);
    double f1 = gsl_sf_erf(y1/M_SQRT2) * (M_SQRTPI/M_SQRT2);  // note: erf from cmath is not accurate!
    double f2 = gsl_sf_erf(y2/M_SQRT2) * (M_SQRTPI/M_SQRT2);
    double val= NAN;
    switch(n) {
        case 0: val = f2 - f1;  break;
        case 1: val = e1 - e2;  break;
        case 2: val = y1 * e1 - y2 * e2  +  f2 - f1;  break;
        case 3: val = (2 + y1 * y1) * e1  -  (2 + y2 * y2) * e2;  break;
        case 4: val = y1 * (3 + y1 * y1) * e1  -  y2 * (3 + y2 * y2) * e2  +  3 * (f2 - f1);  break;
        case 5: val = (8 + y1 * y1 * (4 + y1 * y1)) * e1  -  (8 + y2 * y2 * (4 + y2 * y2)) * e2;  break;
        case 6: val = y1 * (15 + y1 * y1 * (5 + y1 * y1)) * e1
                    - y2 * (15 + y2 * y2 * (5 + y2 * y2)) * e2  +  15 * (f2 - f1);  break;
        case 7: val = (48 + y1 * y1 * (24 + y1 * y1 * (6 + y1 * y1))) * e1
                    - (48 + y2 * y2 * (24 + y2 * y2 * (6 + y2 * y2))) * e2; break;
        default: {  // general case
            exceptionText.clear();
            double v0 = gsl_sf_gamma(0.5*(n+1));
            double v1 = (v0 - gsl_sf_gamma_inc(0.5*(n+1), 0.5*y1*y1)) * (y1>0 || n%2==1 ? 1 : -1);
            double v2 = (v0 - gsl_sf_gamma_inc(0.5*(n+1), 0.5*y2*y2)) * (y2>0 || n%2==1 ? 1 : -1);
            if(exceptionText.empty())
                val = pow(M_SQRT2, n-1) * (v2 - v1);
        }
    }
    return val * pow(sigma, n) / (M_SQRTPI * M_SQRT2);
}

// ------- tools for analyzing the behaviour of a function around a particular point ------- //
// this comes handy in root-finding and related applications, when one needs to ensure that 
// the endpoints of an interval strictly bracked the root: 
// if f(x) is exactly zero at one of the endpoints, and we want to locate the root inside the interval,
// then we need to shift slightly the endpoint to ensure that f(x) is strictly positive (or negative).

EXP PointNeighborhood::PointNeighborhood(const IFunction& fnc, double x0) : absx0(fabs(x0))
{
    // small offset used in computing numerical derivatives, if the analytic ones are not available
    double delta = fmax(fabs(x0) * ROOT3_DBL_EPSILON, 16*DBL_EPSILON);
    // we assume that the function can be computed at all points, but not necessarily with derivatives
    double fplusd = NAN, fderplusd = NAN, fminusd = NAN;
    f0 = fder = fder2=NAN;
    if(fnc.numDerivs()>=2) {
        fnc.evalDeriv(x0, &f0, &fder, &fder2);
        if(isFinite(fder+fder2))
            return;  // no further action necessary
    }
    if(!isFinite(f0))  // haven't called it yet
        fnc.evalDeriv(x0, &f0, fnc.numDerivs()>=1 ? &fder : NULL);
    fnc.evalDeriv(x0+delta, &fplusd, fnc.numDerivs()>=1 ? &fderplusd : NULL);
    if(isFinite(fder)) {
        if(isFinite(fderplusd)) {  // have 1st derivative at both points
            fder2 = (6*(fplusd-f0)/delta - (4*fder+2*fderplusd))/delta;
            return;
        }
    } else if(isFinite(fderplusd)) {  // have 1st derivative at one point only
        fder = 2*(fplusd-f0)/delta - fderplusd;
        fder2= 2*( fderplusd - (fplusd-f0)/delta )/delta;
        return;
    }
    // otherwise we don't have any derivatives computed
    fminusd= fnc(x0-delta);
    fder = (fplusd-fminusd)/(2*delta);
    fder2= (fplusd+fminusd-2*f0)/(delta*delta);
}

EXP double PointNeighborhood::dxToPosneg(double sgn) const
{
    // safety factor to make sure we overshoot in finding the value of opposite sign
    double s0 = sgn*f0 * 1.1;
    double sder = sgn*fder, sder2 = sgn*fder2;
    // offset should be no larger than the scale of variation of the function,
    // but no smaller than the minimum resolvable distance between floating point numbers
    const double delta = fmin(fabs(fder/fder2)*0.5,
        fmax(1000*DBL_EPSILON * absx0, fabs(f0)) / fabs(sder));  //TODO!! this is not satisfactory
    if(s0>0)
        return 0;  // already there
    if(sder==0) {
        if(sder2<=0)
            return NAN;  // we are at a maximum already
        else
            return fmax(sqrt(-s0/sder2), delta);
    }
    // now we know that s0<=0 and sder!=0
    if(sder2>=0)  // may only curve towards zero, so a tangent is a safe estimate
        return -s0/sder + delta*sign(sder);
    double discr = sder*sder - 2*s0*sder2;
    if(discr<=0)
        return NAN;  // never cross zero
    return sign(sder) * (delta - 2*s0/(sqrt(discr)+fabs(sder)) );
}

EXP double PointNeighborhood::dxToNearestRoot() const
{
    if(f0==0) return 0;  // already there
    double dx_nearest_root = -f0/fder;  // nearest root by linear extrapolation, if fder!=0
    if(f0*fder2<0) {  // use quadratic equation to find two nearest roots
        double discr = sqrt(fder*fder - 2*f0*fder2);
        if(fder<0) {
            dx_nearest_root = 2*f0/(discr-fder);
            //dx_farthest_root = (discr-fder)/fder2;
        } else {
            dx_nearest_root = 2*f0/(-discr-fder);
            //dx_farthest_root = (-discr-fder)/fder2;
        }
    }
    return dx_nearest_root;
}

EXP double PointNeighborhood::dxBetweenRoots() const
{
    if(f0==0 && fder==0) return 0;  // degenerate case
    return sqrt(fder*fder - 2*f0*fder2) / fabs(fder2);  // NaN if discriminant<0 - no roots
}

EXP void hermiteDerivs(double x0, double x1, double x2, double f0, double f1, double f2,
    double df0, double df1, double df2, double& der2, double& der3, double& der4, double& der5)
{
    // construct a divided difference table to evaluate 2nd to 5th derivatives via Hermite interpolation:
    const double
    // differences between grid points in x
    dx10   = x1-x0,
    dx21   = x2-x1,
    dx20   = x2-x0,
    sx20   = dx21-dx10,
    // 2nd column
    f00    = df0,
    f01    = (f1    - f0)    / dx10,
    f11    = df1,
    f12    = (f2    - f1)    / dx21,
    f22    = df2,
    // 3rd column
    f001   = (f01   - f00)   / dx10,
    f011   = (f11   - f01)   / dx10,
    f112   = (f12   - f11)   / dx21,
    f122   = (f22   - f12)   / dx21,
    // 4th column
    f0011  = (f011  - f001)  / dx10,
    f0112  = (f112  - f011)  / dx20,
    f1122  = (f122  - f112)  / dx21,
    // 5th column
    f00112 = (f0112 - f0011) / dx20,
    f01122 = (f1122 - f0112) / dx20;
    // 6th column - the tip of triangle, equal to (1/5!) * 5th derivative
    der5   = (f01122- f00112)/ dx20 * 120;
    // start unwinding back
    der4   = (f01122+ f00112) * 12 - sx20 * der5 * 0.3;
    der3   =  f0112 * 6 - sx20 * (der4 + der5 * (sx20 + dx10*dx21/sx20) * 0.2) * 0.25;
    der2   =  f011  + f112 - f0112 * sx20 - dx10 * dx21 * (der4 + der5 * sx20 * 0.2) / 12;

    // alternative expression:
    //der2 = ( -2 * (pow_2(dx21)*(f001-2*f011) + pow_2(dx10)*(f122-2*f112)) +
    //    4*dx10*dx21 * (dx10*f112 + dx21*f011) / dx20 ) / pow_2(dx20);
}


// ------ root finder and minimization routines ------//

namespace {  // internal
/// used in hybrid root-finder to predict the root location by Hermite interpolation:
/// compute the value of f(x) given its values and derivatives at two points x1,x2
/// (x1<=x<=x2 or x1>=x>=x2 is implied but not checked), if the function is expected to be
/// monotonic on this interval (i.e. its derivative does not have roots on x1..x2),
/// otherwise return NAN
inline double hermiteInterpMonotone(double x, double x1, double x2,
    double f1, double f2, double dfdx1, double dfdx2)
{
    // derivatives must exist and have the same sign
    if(!isFinite(dfdx1+dfdx2) || (dfdx1>=0 && dfdx2<0) || (dfdx2>=0 && dfdx1<0))
        return NAN;
    const double dx = x2-x1,
    // check if the interpolant is monotonic on t=[0:1] by solving a quadratic equation:
    // derivative is  a * (t^2 + b/a * t + c/a),  with 0<=t<=1 on the given interval.
    ai = 1 / (-6 * (f2-f1) / dx + 3 * (dfdx1 + dfdx2)),  ba = -1 + (dfdx2 - dfdx1) * ai,  ca = dfdx1 * ai,
    D  = pow_2(ba) - 4 * ca;
    if(!isFinite(ba+D))
        return NAN;  // prudently resign in case of troubles
    if(D >= 0) {     // need to check roots
        double sqD= sqrt(D);
        double t1 = 0.5 * (-ba-sqD);
        double t2 = 0.5 * (-ba+sqD);
        if( (t1>=0 && t1<=1) || (t2>=0 && t2<=1) )
            return NAN;    // there is a root ( y'=0 ) somewhere on the given interval
    }  // otherwise there are no roots
    return hermiteInterp(x, x1, x2, f1, f2, dfdx1, dfdx2);
}

/// the interpolation is accepted if the result differs from one of the endpoints by more than DELTA
static const double DELTA_HERMITE = 1e-12;
}  // internal ns

/// a hybrid between Brent's method and interpolation of root using function derivatives;
/// it is based on the implementation from GSL, original authors: Reid Priedhorsky, Brian Gough
double findRoot(const IFunction& fnc, 
    const double xlower, const double xupper, const double reltoler)
{
	if(reltoler<=0){
		printf("findRoot: relative tolerance must be positive");
		throw std::invalid_argument("findRoot: relative tolerance must be positive");
	}
	if(!isFinite(xlower+xupper)){
		printf("findRoot: endpoints must be finite,\n otherwise need to apply an appropriate scaling transformation manually");
		throw std::invalid_argument("findRoot: endpoints must be finite, "
					    "otherwise need to apply an appropriate scaling transformation manually");
	}
    double a = xlower;
    double b = xupper;
    double fa, fb;
    double fdera = NAN, fderb = NAN;
    bool have_derivs = fnc.numDerivs()>=1;
    fnc.evalDeriv(a, &fa, have_derivs? &fdera : NULL);
    fnc.evalDeriv(b, &fb, have_derivs? &fderb : NULL);
    if((fa < 0.0 && fb < 0.0) || (fa > 0.0 && fb > 0.0) || !isFinite(fa+fb))
        return NAN;   // endpoints do not bracket root
    /*  b  is the current estimate of the root,
        c  is the counter-point (i.e. f(b) * f(c) < 0, and |f(b)| < |f(c)| ),
        a  is the previous estimate of the root:  either
           (1) a==c, or 
           (2) f(a) has the same sign as f(b), |f(a)|>|f(b)|, and a, b, c form a monotonic sequence.
    */
    double c = a;
    double fc = fa;
    double fderc = fdera;
    double d = b - c;   // always holds the (signed) length of current interval
    double e = b - c;   // this is used to estimate roundoff (?)
    if (fabs(fc) < fabs(fb)) {  // swap b and c so that |f(b)| < |f(c)|
        a = b;
        b = c;
        c = a;
        fa = fb;
        fb = fc;
        fc = fa;
        fdera = fderb;
        fderb = fderc;
        fderc = fdera;
    }
    int numIter = 0;
    bool converged = false;
    double abstoler = fabs(xlower-xupper) * reltoler;
    do {
        double tol = 0.5 * DBL_EPSILON * fabs(b);
        double cminusb = c-b, m = 0.5 * cminusb;
        if(fb == 0 || fabs(m) <= tol) 
            return b;  // the ROOT
        if(fabs(e) < tol || fabs(fa) <= fabs(fb)) {  // use bisection
            d = m;
            e = m;
        } else {
            double dd = NAN;
            if(have_derivs && fderb * fderc > 0)  // derivs exist and have the same sign
            {   // attempt to obtain the approximation by Hermite interpolation
                dd = hermiteInterpMonotone(0, fb, fc, 0, cminusb, 1/fderb, 1/fderc);
            }
            if(isFinite(dd) && std::min(fabs(dd), fabs(cminusb-dd)) > fabs(cminusb) * DELTA_HERMITE) {
                d = dd;           // Hermite interpolation is successful
            } else {              // otherwise proceed as usual in the Brent method
                double p, q, r, s = fb / fa;
                if (a == c) {     // secant method (linear interpolation)
                    p = 2 * m * s;
                    q = 1 - s;
                } else {          // inverse quadratic interpolation
                    q = fa / fc;
                    r = fb / fc;
                    p = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
                    q = (q - 1) * (r - 1) * (s - 1);
                }
                if(p > 0)
                    q = -q;
                else
                    p = -p;
                if(2 * p < std::min(3 * m * q - fabs(tol * q), fabs(e * q))) { 
                    e = d;
                    d = p / q;
                } else {
                    /* interpolation failed, fall back to bisection */
                    d = m;
                    e = m;
                }
            }
        }

        a = b;
        fa = fb;
        fdera = fderb;
        if (fabs(d) > tol)
            b += d;
        else
            b += (m > 0 ? +tol : -tol);

        fnc.evalDeriv(b, &fb, have_derivs? &fderb : NULL);
        if(!isFinite(fb))
            return NAN;

        /* Update the best estimate of the root and bounds on each iteration */
        if((fb < 0 && fc < 0) || (fb > 0 && fc > 0)) {   // the root is between 'a' and the new 'b'
            c = a;       // so the new counterpoint is moved to the old 'a'
            fc = fa;
            fderc = fdera;
            d = b - c;
            e = b - c;
        }
        if(fabs(fc) < fabs(fb)) {   // ensure that 'b' is close to zero than 'c'
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
            fdera = fderb;
            fderb = fderc;
            fderc = fdera;
        }

        numIter++;
        if(fabs(b-c) <= abstoler) { // convergence criterion from bracketing algorithm
            converged = true;
            double offset = fb*(b-c)/(fc-fb);  // offset from b to the root via secant
            if((offset>0 && offset<c-b) || (offset<0 && offset>c-b))
                b += offset;        // final secant step
        } else if(have_derivs) {    // convergence from derivatives
            double offset = -fb / fderb;       // offset from b to the root via Newton's method
            bool bracketed = (offset>0 && offset<c-b) || (offset<0 && offset>c-b);
            if(bracketed && fabs(offset) < abstoler) {
                converged = true;
                b += offset;        // final Newton step
            }
        }
        if(numIter >= MAXITER) {
		converged = true;  // not quite ready, but can't loop forever
		printf("MAXITER exceeded in findRoot\n");
            utils::msg(utils::VL_WARNING, "findRoot", "max # of iterations exceeded: "
                "x="+utils::toString(b,15)+" +- "+utils::toString(fabs(b-c))+
                " on interval ["+utils::toString(xlower,15)+":"+utils::toString(xupper,15)+
                "], req.toler.="+utils::toString(abstoler));
        }
    } while(!converged);
    return b;  // best approximation
}

namespace{  // internal
/// choose a point inside an interval based on the golden section rule
inline double minGuess(double x1, double x2, double y1, double y2)
{
    const double golden = 0.618034;
    if(y1<y2)
        return x1 * golden + x2 * (1-golden);
    else
        return x2 * golden + x1 * (1-golden);
}
}  // internal ns

/// invoke the minimization routine with a valid initial guess (if it was not provided,
/// try to come up with a plausible one)
double findMin(const IFunction& fnc, double xlower, double xupper, double xinit, double reltoler)
{
    if(reltoler<=0)
        throw std::invalid_argument("findMin: relative tolerance must be positive");
    if(xlower>xupper)
        std::swap(xlower, xupper);
    double ylower = fnc(xlower);
    double yupper = fnc(xupper);
    if(xinit == xinit) {
        if(! (xinit >= xlower && xinit <= xupper) )
            throw std::invalid_argument("findMin: initial guess is outside the search interval");
    } else {    // initial guess not provided - choose a plausible point inside the interval
        xinit = minGuess(xlower, xupper, ylower, yupper);
    }
    double yinit  = fnc(xinit);
    if(!isFinite(ylower+yupper+yinit))
        return NAN;
    int iter = 0;
    double abstoler = reltoler * fabs(xupper-xlower);

    // if the initial guess does not enclose minimum, provide a new guess inside a smaller range
    while( (yinit>=ylower || yinit>=yupper) && iter < MAXITER && fabs(xlower-xupper) > abstoler) {
        if(ylower<yupper) {
            xupper = xinit;
            yupper = yinit;
        } else {
            xlower = xinit;
            ylower = yinit;
        }
        xinit = minGuess(xlower, xupper, ylower, yupper);
        yinit = fnc(xinit);
        if(!isFinite(yinit))
            return NAN;
        iter++;
    }
    if(yinit>=ylower || yinit>=yupper)
        // couldn't locate a minimum strictly inside the interval, so return one of the endpoints
        return ylower<yupper ? xlower : xupper;

    // use the GSL minimizer with a known initial point
    gsl_function F;
    F.function = &functionWrapper;
    F.params = const_cast<IFunction*>(&fnc);
    gsl_min_fminimizer *minser = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
    double xroot = NAN;
    gsl_min_fminimizer_set_with_values(minser, &F, xinit, yinit, xlower, ylower, xupper, yupper);
    iter=0;
    exceptionText.clear();
    while(iter < MAXITER && fabs(xlower-xupper) > abstoler) {
        iter++;
        gsl_min_fminimizer_iterate(minser);
        if(!exceptionText.empty()) {
            xroot = NAN;
            break;
        }
        xroot  = gsl_min_fminimizer_x_minimum(minser);
        xlower = gsl_min_fminimizer_x_lower(minser);
        xupper = gsl_min_fminimizer_x_upper(minser);
    }
    gsl_min_fminimizer_free(minser);
    return xroot;
}


// ------- integration routines ------- //

double integrate(const IFunction& fnc, double x1, double x2, double reltoler, 
    double* error, int* numEval)
{
    if(x1==x2)
        return 0;
    gsl_function F;
    F.function = functionWrapper;
    F.params = const_cast<IFunction*>(&fnc);
    double result, dummy;
    size_t neval;
    CALL_FUNCTION_OR_THROW(
    gsl_integration_qng(&F, x1, x2, 0, reltoler, &result, error!=NULL ? error : &dummy, &neval))
    if(numEval!=NULL)
        *numEval = neval;
    return result;
}

double integrateAdaptive(const IFunction& fnc, double x1, double x2, double reltoler, 
    double* error, int* numEval)
{
    if(x1==x2)
        return 0;
    gsl_function F;
    F.function = functionWrapper;
    F.params = const_cast<IFunction*>(&fnc);
    double result, dummy;
    size_t neval;
    gsl_integration_cquad_workspace* ws=gsl_integration_cquad_workspace_alloc(MAX_INTEGR_POINTS);
    CALL_FUNCTION_OR_THROW(
    gsl_integration_cquad(&F, x1, x2, 0, reltoler, ws, &result, error!=NULL ? error : &dummy, &neval))
    gsl_integration_cquad_workspace_free(ws);
    if(numEval!=NULL)
        *numEval = neval;
    return result;
}

// this routine is intended to be fast, so only works with pre-computed integration tables
double integrateGL(const IFunction& fnc, double x1, double x2, int N)
{
    if(x1==x2)
        return 0;
    if(N < 1 || N > MAX_GL_ORDER)
        throw std::invalid_argument("integrateGL: order is too high (not implemented)");
    // use pre-computed tables of points and weights (they are not available for every N,
    // so take the closest implemented one with at least the requested number of points)
    while(GLPOINTS[N] == NULL) N++;
    const double *points = GLPOINTS[N], *weights = GLWEIGHTS[N];
    double result = 0;
    for(int i=0; i<N; i++)
        result += weights[i] * fnc(x2 * points[i] + x1 * (1-points[i]));
    return result * (x2-x1);
}

// compute the node and weight of k-th node in the Gauss-Legendre quadrature rule with n points (0<=k<n)
void getGLnodeAndWeight(int n, int k, /*output: node*/ double& x, /*output: weight*/ double& w)
{
    // initial guess for the k-th root of Legendre polynomial P_n(x)
    x = (n==k*2+1) ? 0 : (1 - 0.125 * (n-1) / (n*n*n)) * cos(M_PI * (k+0.75) / (n+0.5));
    while(1) {
        // compute Legendre polynomial P_n(x) and its two derivatives
        double Pm2, Pm1 = x, P = 1.5 * x*x - 0.5;
        for(int l=3; l<=n; l++) {
            // recurrence relation assuming n>=2
            Pm2 = Pm1;
            Pm1 = P;
            P   = x * Pm1 + (1-1./l) * (x * Pm1 - Pm2);
        }
        double
        der = (Pm1 - x*P) * n / (1-x*x),           // first derivative of P_n(x)
        der2= (2*x*der - n*(n+1)*P) / (1-x*x),     // second derivative
        dx  = -P / (der - 0.5 * P * der2 / der);   // Halley root-finding method for x
        x  += dx;
        if(fabs(dx) <= 1e-10) {
            // if dx is so small, the next iteration will make it less than DBL_EPSILON, so is unnecessary.
            // apply a final correction for the first derivative, and compute the weight
            w = 2 / ( (1-x*x) * pow_2(der + dx*der2) );
            return;
        }
    }
}

void prepareIntegrationTableGL(double x1, double x2, int N, double* nodes, double* weights)
{
    if(N <= MAX_GL_ORDER && GLPOINTS[N] != NULL) {
        // retrieve a pre-computed table
        for(int i=0; i<=N/2; i++) {
            nodes  [i]     = x1 + GLPOINTS [N][i] * (x2-x1);
            nodes  [N-1-i] = x2 - GLPOINTS [N][i] * (x2-x1);
            weights[i] = weights[N-1-i] = GLWEIGHTS[N][i] * (x2-x1);
        }
    } else {
        // compute a table on-the-fly
        double x, w;
        for(int i=0; i<=N/2; i++) {
            getGLnodeAndWeight(N, i, /*output*/ x, w);
            nodes  [i]     = 0.5 * (x1 + x2 - x * (x2-x1) );
            nodes  [N-1-i] = 0.5 * (x1 + x2 + x * (x2-x1) );
            weights[i] = weights[N-1-i] = 0.5 * w * (x2-x1);
        }
    }
}


// ------- multidimensional integration ------- //
namespace {
#ifdef HAVE_CUBA
// wrapper for the Cuba library
struct CubaParams {
    const IFunctionNdim& F; ///< the original function
    const double* xlower;   ///< lower limits of integration
    const double* xupper;   ///< upper limits of integration
    std::string error;      ///< store error message in case of exception
    CubaParams(const IFunctionNdim& _F, const double* _xlower, const double* _xupper) :
        F(_F), xlower(_xlower), xupper(_xupper) {}
};

int integrandNdimWrapperCuba(const int *ndim, const double xscaled[],
    const int *fdim, double fval[], void *v_param, const int *npoints)
{
    CubaParams* param = static_cast<CubaParams*>(v_param);
    assert(*ndim == (int)param->F.numVars() && *fdim == (int)param->F.numValues());
    try {
        // un-scale the input point(s) from [0:1]^N to the original range:
        // allocate a temporary array for un-scaled values on the stack (no need to free it)
        double* xval = (double*) alloca( (*ndim) * (*npoints) * sizeof(double));
        for(int i=0; i<*npoints; i++) {
            for(int d=0, s=(*ndim)*i; d< *ndim; d++, s++)
                xval[s] = param->xlower[d] * (1-xscaled[s]) + param->xupper[d] * xscaled[s];
        }
        param->F.evalmany(*npoints, xval, fval);
        // check if the result is not finite (not performed unless in debug mode)
        if(utils::verbosityLevel >= utils::VL_WARNING) {
            for(int i=0; i< *npoints; i++)
                for(int f=0; f< *fdim; f++)
                    if(!isFinite(fval[f + (*fdim)*i])) {
                        param->error = "integrateNdim: invalid function value encountered at";
                        for(int d=0; d< *ndim; d++)
                            param->error += ' ' + utils::toString(xval[d + (*ndim)*i], 15);
                        param->error += '\n' + utils::stacktrace();
                        return -1;
                    }
        }
        return 0;   // success
    }
    catch(std::exception& e) {
	    param->error = std::string("integrateNdim: ") + e.what() + '\n' + utils::stacktrace();
	    std::cout << param->error << "\n";
	    printf("%s",(param->error).c_str());
        return -999;  // signal of error
    }
}

#else
// wrapper for the Cubature library
struct CubatureParams {
    const IFunctionNdim& F; ///< the original function
    int numEval;            ///< count the number of function evaluations
    std::string error;      ///< store error message in case of exception
    explicit CubatureParams(const IFunctionNdim& _F) :
        F(_F), numEval(0) {}
};

int integrandNdimWrapperCubature(unsigned int ndim, unsigned int npoints, const double *xval,
    void *v_param, unsigned int fdim, double *fval)
{
    CubatureParams* param = static_cast<CubatureParams*>(v_param);
    assert(ndim == param->F.numVars() && fdim == param->F.numValues());
    try {
        param->F.evalmany(npoints, xval, fval);
        param->numEval += npoints;
        // check if the result is not finite (only performed in debug mode)
        if(utils::verbosityLevel >= utils::VL_WARNING) {
            for(unsigned int i=0; i<npoints; i++)
                for(unsigned int f=0; f<fdim; f++)
			if(!isFinite(fval[f + i*fdim])) {
				printf("integrateNdim: invalid function value encountered at");
                        param->error = "integrateNdim: invalid function value encountered at";
                        for(unsigned int d=0; d<ndim; d++)
                            param->error += ' ' + utils::toString(xval[d + i*ndim], 15);
                        param->error += '\n' + utils::stacktrace();
                        return -1;
                    }
        }
        return 0;   // success
    }
    catch(std::exception& e) {
	    param->error = std::string("integrateNdim: ") + e.what() + '\n' + utils::stacktrace();
	    std::cout << param->error << "\n";
        return -1;  // signal of error
    }
}
#endif
}  // namespace

void integrateNdim(const IFunctionNdim& F, const double xlower[], const double xupper[], 
    const double relToler, const unsigned int maxNumEval, 
    double result[], double outError[], int* numEval)
{
    const unsigned int numVars = F.numVars();
    const unsigned int numValues = F.numValues();
    const double absToler = 0;  // the only possible way to stay invariant under scaling transformations
    // storage for errors: in the case that user doesn't need them, allocate a temp.array on the stack,
    // which will be automatically freed on exit (NOTE: this assumes that numValues is not too large!)
    double* error = outError!=NULL ? outError : static_cast<double*>(alloca(numValues * sizeof(double)));
#ifdef HAVE_CUBA
    CubaParams param(F, xlower, xupper);
    // allocate another temp.array, unused
    double* tempProb = static_cast<double*>(alloca(numValues * sizeof(double)));
    int nregions, neval, fail;
    const int NVEC = 1000, FLAGS = 0, KEY = 7, minNumEval = 0;
    cubacores(0, 0);  // disable parallelization at the CUBA level
    Cuhre(numVars, numValues, (integrand_t)&integrandNdimWrapperCuba, &param, NVEC,
        relToler, absToler, FLAGS, minNumEval, maxNumEval, 
        KEY, NULL/*STATEFILE*/, NULL/*spin*/,
        &nregions, numEval!=NULL ? numEval : &neval, &fail, 
        result, error, tempProb);
    if(fail==-1){
	    printf("integrateNdim: number of dimensions is too large\n");
	    throw std::runtime_error("integrateNdim: number of dimensions is too large");
    }if(fail==-2){
	    printf("integrateNdim: number of components is too large");
	    throw std::runtime_error("integrateNdim: number of components is too large");
    }
    // need to scale the result to account for coordinate transformation [xlower:xupper] => [0:1]
    double scaleFactor = 1.;
    for(unsigned int n=0; n<numVars; n++)
        scaleFactor *= (xupper[n]-xlower[n]);
    for(unsigned int m=0; m<numValues; m++) {
        result[m] *= scaleFactor;
        error [m] *= scaleFactor;
    }
    if(!param.error.empty()){
	    printf("integrateNdim: param.error\n");
	    throw std::runtime_error(param.error);
    }
#else
    CubatureParams param(F);
    int e=hcubature_v(numValues, &integrandNdimWrapperCubature, &param,
        numVars, xlower, xupper, maxNumEval, absToler, relToler,
		      ERROR_INDIVIDUAL, result, error);
    if(numEval!=NULL)
        *numEval = param.numEval;
    if(!param.error.empty()){
	    printf("integrateNdim: param.error\n");
	    throw std::runtime_error(param.error);
    }
#endif
}

EXP void getDFT(double data[],const int N){
	gsl_fft_real_radix2_transform(data,1,N);
}

EXP line& line::operator= (const line& aline){
	nu=aline.nu; A=aline.A; phi=aline.phi; resid=aline.resid; diag=aline.diag;
	return *this;
}
EXP FrequencyFinder::FrequencyFinder(double* ts,double* tc,const int _NF,double _deltat,
		       const double _frqmin,const double _frqmax) :
    NF(_NF), deltat(_deltat), frqmin(_frqmin), frqmax(_frqmax){
	NF5=NF/2;
	z = new double[NF5+2];
	Z=new cmplx[NF5+1];
	for(int i=0;i<NF5;i++)
		Z[i]=cmplx(tc[i],ts[i]);
	Z[NF5]=cmplx(tc[NF5],0);
	discrm=0.01;
	sin0=sin(M_PI/(double)NF); cos0=cos(M_PI/(double)NF);
	wmin=2.*M_PI/((double)NF*deltat);
}
EXP FrequencyFinder::FrequencyFinder(double* data,const int _NF,double _deltat,
		       const double _frqmin,const double _frqmax) :
    NF(_NF), deltat(_deltat), frqmin(_frqmin), frqmax(_frqmax){
	NF5=NF/2;
	z = new double[NF5+2];
	Z = new cmplx[NF5+1];
	math::getDFT(data,NF);
	for(int i=1; i<NF5; i++){
		Z[i]=cmplx(data[i],data[NF-i]);
	}
	Z[0]=cmplx(data[0],0); Z[NF5]=cmplx(data[NF5],0);
	discrm=0.01;
	sin0=sin(M_PI/(double)NF); cos0=cos(M_PI/(double)NF);
	wmin=2.*M_PI/((double)NF*deltat);
}
EXP std::vector<line> FrequencyFinder::analyse(double& resid){
	std::vector<line> lines;
	pwr0=get_pwr(); pwr_best=1e6;
	double pwr=pwr0, zmx, zstop;
	int nlines=0;
	while(pwr>.001*pwr0 && pwr<pwr_best){
		pwr_best=pwr;
		int jmx=difference(zmx);//compute z, the second difference of the spectrum
		if(lines.size() == 0) zstop=1e-8*zmx;
		if(jmx <= 2)  jmx=0;
		if(z[jmx] < zstop) break;
		for(int i=0;i<4;i++){//evaluate quantities near point of peak curvature
			int k=jmx-1+i;
			zl[i]=(Zt(k-1)-Zt(k))+(Zt(k+1)-Zt(k));
			if(k>=0) z[k]=sqrt(z[k]);
		}
		bool isol=true;
		for(int i=0; i<3; i+=2){
			isol=(isol && fabs(std::imag(zl[i]/zl[1]))<discrm);
		}
		if(isol || jmx<=1){
			if(isol && jmx>0)
				pwr=isolated(jmx,lines);
			else
				pwr=lowFreq(jmx,lines);
		} else
			pwr=pair(jmx,lines);
	}
	order(lines);
	resid=pwr_best/pwr0;
	return lines;
}
double FrequencyFinder::get_pwr(void){//pwr is a measure of the unclaimed portion of the spectrum
	double pwr=0,Zpwr=0;
	for(int i=0;i<NF5;i++)
		Zpwr+=std::norm(Z[i]);
	return sqrt(Zpwr/(double)(NF5));
}
int FrequencyFinder::difference(double& zmx){
	zmx=0;
	int jmx=0;
	for(int j=0;j<NF5-1;j++){
		z[j]=std::norm((Zt(j)-Zt(j-1))+(Zt(j)-Zt(j+1)));
		if(z[j]>zmx){
			zmx=z[j]; jmx=j;
		}
	}
	return jmx;
}
void FrequencyFinder::Zsubtract(int jmx,double x,cmplx zox){
	double sp = sin(M_PI*(x-(double)jmx)/(double)NF);//sin(pi*x_k/N) at origin
	double cp = cos(M_PI*(x-(double)jmx)/(double)NF);//cos same
	double sm = -sp, cm = cp;//sin & cos of pi*xk/N of line at -ve frequency
	cmplx zp = zox*x, zm = -std::conj(zp);
	for(int j=0; j<NF5; j++){
		if(j == jmx && fabs(x) < .01){
			if(jmx>0 || fabs(x)>.001)
				Z[j]-=cmplx(zox*cmplx(cp*(double)NF/M_PI,x)+zm*cmplx(cm/sm,1));
			else
				Z[j]-=cmplx(zox*cmplx(cp*(double)NF/M_PI,x));
		}else
			Z[j]-=zp*cmplx(cp/sp,1) + zm*cmplx(cm/sm,1);
//  Now add pi/N to args of sines & cosines
		double save = sp*cos0 + cp*sin0;//sp advanced
		cp = cp*cos0 - sp*sin0;
		sp = save;
		save = sm*cos0 + cm*sin0;//sm advanced
		cm = cm*cos0 - sm*sin0;
		sm = save;
	}
}
double FrequencyFinder::isolated(int jmx,std::vector<line>& lines){
	line aline;
	int k= z[jmx+1] > z[jmx-1]? k=jmx+1 : jmx-1;
	double wcor=(double)(jmx-k)*(2*z[k]-z[jmx])/(z[jmx]+z[k]);
	double freq=((double)jmx-wcor)*wmin;
	cmplx zox=.5*M_PI/(double)NF*(wcor*wcor-1)*zl[1];
	if(freq < frqmin || freq >= frqmax){
		Zsubtract(jmx,wcor,zox);
		return get_pwr();
	}
	aline.diag=0;
	aline.nu=freq;
	aline.A=fabs(wcor*wcor-1.)*z[jmx]/(double)(NF);
	if(fabs(wcor) > 1e-9) aline.A*=M_PI*wcor/sin(M_PI*wcor);
	aline.phi=M_PI*wcor+std::arg(zox);
	if(wcor*sin(M_PI*wcor) < 0) aline.phi+=M_PI;
	Zsubtract(jmx,wcor,zox);
	double pwr=get_pwr(); aline.resid=pwr/pwr0;
	if(pwr<pwr_best) lines.push_back(aline);
	return pwr;
}
double FrequencyFinder::lowFreq(int jmx,std::vector<line>& lines){//algorithm for low-frequency component
	line aline;
	double wcor,freq;
	cmplx zox;
	if(fabs(std::imag(zl[2])) < 1e-2*fabs(std::real(zl[2]))){
		wcor=0; zox=-.5*std::real(zl[1]);
	}else{
		wcor=std::imag(2.*zl[2]+9.*zl[3])/std::imag(zl[3]-2.*zl[2]);
		wcor=-sqrt(fmax(wcor,0));
		if(wcor == 0){//there's only the zero-freq line
			zox=-std::real(zl[1]);
		}else{
			zox=std::real(zl[2]*(wcor*wcor-4)-zl[1]*(wcor*wcor+2))/(3*wcor*wcor);
		}
	}
	if(frqmin <= 0){
		aline.diag=1; aline.nu=0;
		aline.A=std::abs(zox)/(double)(NF);
		aline.phi=std::arg(zox);
	}
	zl[1]+=2.*zox; zl[2]-=2.*zox;
	zox*=M_PI/(double)NF;
	Zsubtract(jmx,0,zox);
	double pwr=get_pwr();
	if(frqmin <= 0){
		aline.resid=pwr/pwr0;
		if(pwr<pwr_best) lines.push_back(aline);
	}			
	if(fabs(wcor) > 4 || -wcor < 1e-5){
		return pwr;
	}
	double ampr,ampi,tmod1,tmod2,fac;
	if(wcor < -1){
		double wcor2=wcor*wcor;
		ampr=-std::real(zl[3])*(9-wcor2)*(4-wcor2)*(1-wcor2)/(2*wcor2*(11+wcor2));
		ampi=-(4-wcor2)*(1-wcor2)/(6.*wcor)*
		     ((wcor+2)*std::imag(zl[2])+.5*(wcor+1)*
		      (9-wcor2)*std::imag(zl[3])/(1+wcor2));
		tmod2=.5*M_PI/(double)NF*ampr;
		fac=wcor*M_PI/sin(wcor*M_PI);
		ampr*=fac; ampi*=fac;
	}else{
		ampr=.5*std::real(zl[1])*(wcor*wcor-1);
		tmod2=.5*M_PI/(double)NF*ampr;
		if(-wcor > 1e-5){
			ampr*=M_PI*wcor/sin(M_PI*wcor);
			ampi=ampr*std::imag(zl[2])*(4-wcor*wcor)/(3*std::real(zl[1])*wcor);
		} else
			ampi=0;
	}
	tmod1=tmod2*ampi/ampr;
	freq=-wcor*wmin;
	if(freq < frqmin){
		Zsubtract(jmx,wcor,cmplx(tmod2,tmod1));
		return get_pwr(); 
	}
	aline.diag=2; aline.nu=freq;
	aline.A=sqrt(pow_2(ampr)+pow_2(ampi))/(double)(NF);
	aline.phi=M_PI*wcor+atan2(tmod1,tmod2);
	if(wcor*sin(M_PI*wcor) < 0) aline.phi+=M_PI;
	Zsubtract(jmx,wcor,cmplx(tmod2,tmod1));
	pwr=get_pwr(); aline.resid=pwr/pwr0;
	if(pwr<pwr_best) lines.push_back(aline);
	return pwr;
}
double FrequencyFinder::pair(int jmx,std::vector<line>& lines){// algorithm for examining close pairs of lines
	line aline;
	double a=std::imag((zl[1]-zl[2])*std::conj(zl[1]-zl[0]));
	double b=std::imag((zl[1]-zl[2])*std::conj(zl[1]+2.*zl[0])
			   -(zl[1]+2.*zl[2])*std::conj(zl[1]-zl[0]));
	double c=-std::imag((zl[1]+2.*zl[2])*std::conj(zl[1]+2.*zl[0]));
	double disc=(.25*b*b-a*c);
	if(disc < 0 || a == 0)  return isolated(jmx, lines);
	disc=sqrt(disc);
	double pwr, x[2]={(-.5*b+disc)/a, (-.5*b-disc)/a};
	if((double)jmx < fmin(x[0],x[1]))  return get_pwr();
	if(fabs(x[1]) < fabs(x[0])) std::swap(x[0],x[1]);
	for(int i=0;i<2;i++){
		if(fabs(x[i])>0.9) continue;
		int j=(1+i)%2;//now get 2 estimates of the amplitude
		double fac1=(x[i]+1)/((x[j]-1)/(x[i]-1)-(x[j]+2)/(x[i]+2));
		cmplx amp1=fac1*(zl[1]*(x[j]-1)-zl[2]*(x[j]+2));
		double fac2=(x[i]-1)/((x[j]+1)/(x[i]+1)-(x[j]-2)/(x[i]-2));
		cmplx amp2=fac2*(zl[1]*(x[j]+1)-zl[0]*(x[j]-2));
		cmplx amp=.5*(amp2*(1+2*x[i])+amp1*(1-2*x[i]));//merged estimate
		double freq=((double)jmx-x[i])*wmin;
		amp*=.5*M_PI/(double)NF;
		Zsubtract(jmx,x[i],amp);
		pwr=get_pwr();
		if(freq >= frqmin && freq < frqmax){
			aline.diag=3; aline.nu=freq;
			aline.A=sqrt(std::norm(amp))*2/M_PI;
			if(fabs(x[i]) > 1e-5)  aline.A *=fabs(M_PI*x[i]/sin(M_PI*x[i]));
			aline.phi=M_PI*x[i]+std::arg(amp);
			aline.resid=pwr/pwr0;
			if(x[i]*sin(M_PI*x[i]) < 0) aline.phi+=M_PI;
			if(pwr<pwr_best) lines.push_back(aline);
		}
	}
	return pwr;
}
void FrequencyFinder::order(std::vector<line>& lines){
	bool ordered;
	do{
		int i=lines.size()-1; ordered=true;
		while(i>0){
			if(lines[i].A>lines[i-1].A){
				std::swap(lines[i-1],lines[i]);
				ordered=false;
			}
			i--;
		}
	} while(!ordered);
}

}  // namespace
