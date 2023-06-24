/** \file   math_simple_cubature.h
    \brief  fixed-order rules for 2d and 3d integration
    \date   May 2017
    \author Eugene Vasiliev

    Computation of N-dimensional integrals using tailored quadrature rules.
    This file contains only the data for the case of 2d and 3d quadrature rules
    with 12 and 34 points, correspondingly. They are taken from Witherden&Vincent(2015),
    are exact for polynomials up to a total degree of 7 (or up to degree 3 in each dimension),
    evaluated using only interior points and have positive weight coefficients.
    They should be suitable for evaluating integrals over basis functions in finite-element methods.
*/
#pragma once

namespace math{

/// number of points for the 2d rule
static const int NUM_INT_POINTS_2D = 12;
/// points (x,y) and weights of the 2d rule for unit rectangle
static const double  INT_POINTS_2D[NUM_INT_POINTS_2D][3] = {
{9.629100498862757096e-01, 5.000000000000000000e-01, 6.049382716049382658e-02},
{5.000000000000000000e-01, 9.629100498862757096e-01, 6.049382716049382658e-02},
{3.708995011372429040e-02, 5.000000000000000000e-01, 6.049382716049382658e-02},
{5.000000000000000000e-01, 3.708995011372429040e-02, 6.049382716049382658e-02},
{9.029898914592993098e-01, 9.029898914592993098e-01, 5.935794367265755794e-02},
{9.029898914592993098e-01, 9.701010854070063472e-02, 5.935794367265755794e-02},
{9.701010854070063472e-02, 9.029898914592993098e-01, 5.935794367265755794e-02},
{9.701010854070063472e-02, 9.701010854070063472e-02, 5.935794367265755794e-02},
{6.902772166041578616e-01, 6.902772166041578616e-01, 1.301482291668486224e-01},
{6.902772166041578616e-01, 3.097227833958421384e-01, 1.301482291668486224e-01},
{3.097227833958421384e-01, 6.902772166041578616e-01, 1.301482291668486224e-01},
{3.097227833958421384e-01, 3.097227833958421384e-01, 1.301482291668486224e-01}};

/// number of points for the 3d rule
static const int NUM_INT_POINTS_3D = 34;
/// points (x,y,z) and weights of the 3d rule for unit cube
static const double  INT_POINTS_3D[NUM_INT_POINTS_3D][4] = {
{5.956919405661731393e-03, 5.000000000000000000e-01, 5.000000000000000000e-01, 2.501623637294587629e-02},
{5.000000000000000000e-01, 5.000000000000000000e-01, 9.940430805943383241e-01, 2.501623637294587629e-02},
{5.000000000000000000e-01, 9.940430805943383241e-01, 5.000000000000000000e-01, 2.501623637294587629e-02},
{5.000000000000000000e-01, 5.000000000000000000e-01, 5.956919405661731393e-03, 2.501623637294587629e-02},
{9.940430805943383241e-01, 5.000000000000000000e-01, 5.000000000000000000e-01, 2.501623637294587629e-02},
{5.000000000000000000e-01, 5.956919405661731393e-03, 5.000000000000000000e-01, 2.501623637294587629e-02},
{7.039775836791597463e-01, 2.960224163208403092e-01, 2.960224163208403092e-01, 5.714423200873159497e-02},
{2.960224163208403092e-01, 7.039775836791597463e-01, 7.039775836791597463e-01, 5.714423200873159497e-02},
{2.960224163208403092e-01, 7.039775836791597463e-01, 2.960224163208403092e-01, 5.714423200873159497e-02},
{2.960224163208403092e-01, 2.960224163208403092e-01, 2.960224163208403092e-01, 5.714423200873159497e-02},
{2.960224163208403092e-01, 2.960224163208403092e-01, 7.039775836791597463e-01, 5.714423200873159497e-02},
{7.039775836791597463e-01, 7.039775836791597463e-01, 2.960224163208403092e-01, 5.714423200873159497e-02},
{7.039775836791597463e-01, 7.039775836791597463e-01, 7.039775836791597463e-01, 5.714423200873159497e-02},
{7.039775836791597463e-01, 2.960224163208403092e-01, 7.039775836791597463e-01, 5.714423200873159497e-02},
{8.905514105020593263e-01, 1.094485894979407292e-01, 1.094485894979407292e-01, 1.922451750824483510e-02},
{1.094485894979407292e-01, 8.905514105020593263e-01, 8.905514105020593263e-01, 1.922451750824483510e-02},
{1.094485894979407292e-01, 8.905514105020593263e-01, 1.094485894979407292e-01, 1.922451750824483510e-02},
{1.094485894979407292e-01, 1.094485894979407292e-01, 1.094485894979407292e-01, 1.922451750824483510e-02},
{1.094485894979407292e-01, 1.094485894979407292e-01, 8.905514105020593263e-01, 1.922451750824483510e-02},
{8.905514105020593263e-01, 8.905514105020593263e-01, 1.094485894979407292e-01, 1.922451750824483510e-02},
{8.905514105020593263e-01, 8.905514105020593263e-01, 8.905514105020593263e-01, 1.922451750824483510e-02},
{8.905514105020593263e-01, 1.094485894979407292e-01, 8.905514105020593263e-01, 1.922451750824483510e-02},
{7.597386215798063169e-02, 7.597386215798063169e-02, 5.000000000000000000e-01, 1.991271546887610905e-02},
{9.240261378420193683e-01, 5.000000000000000000e-01, 7.597386215798063169e-02, 1.991271546887610905e-02},
{5.000000000000000000e-01, 9.240261378420193683e-01, 7.597386215798063169e-02, 1.991271546887610905e-02},
{9.240261378420193683e-01, 9.240261378420193683e-01, 5.000000000000000000e-01, 1.991271546887610905e-02},
{9.240261378420193683e-01, 5.000000000000000000e-01, 9.240261378420193683e-01, 1.991271546887610905e-02},
{5.000000000000000000e-01, 7.597386215798063169e-02, 9.240261378420193683e-01, 1.991271546887610905e-02},
{5.000000000000000000e-01, 7.597386215798063169e-02, 7.597386215798063169e-02, 1.991271546887610905e-02},
{7.597386215798063169e-02, 5.000000000000000000e-01, 9.240261378420193683e-01, 1.991271546887610905e-02},
{7.597386215798063169e-02, 9.240261378420193683e-01, 5.000000000000000000e-01, 1.991271546887610905e-02},
{9.240261378420193683e-01, 7.597386215798063169e-02, 5.000000000000000000e-01, 1.991271546887610905e-02},
{5.000000000000000000e-01, 9.240261378420193683e-01, 9.240261378420193683e-01, 1.991271546887610905e-02},
{7.597386215798063169e-02, 5.000000000000000000e-01, 7.597386215798063169e-02, 1.991271546887610905e-02}};

} // namespace
