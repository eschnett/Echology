#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include "stdio.h"
#include "math.h"
using namespace std;

// solve for r 
/*range of r from 2m to the domain*/

double fun(double rho, double r) {
  DECLARE_CCTK_PARAMETERS;
return pow(sqrt(1-2*param_m/r)*r + 2 * param_m * atanh(sqrt(1-2* param_m /r))+param_c,2)-pow(rho,2);
}

double solve(double rho)
{
DECLARE_CCTK_PARAMETERS;
double low=2 * param_m ,up=20,mid=(low+up)/2;
while(up-low>1e-6)
{
double y = fun(rho,mid);
if(y>0)up=mid;
else if(y<0)low=mid;
else break;
mid=(low+up)/2;
}
return mid;
}

extern "C" void FakeMatter_DefineMetricexp(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

// Loop over all grid points
for (int k = 0; k < cctk_lsh[2]; ++k) {
  for (int j = 0; j < cctk_lsh[1]; ++j) {
    for (int i = 0; i < cctk_lsh[0]; ++i) {
        int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
        double rho1 = r[ijk];
        double rr = solve(rho1);
// cerr << "r=" << rr << "\n";
// cerr << "rho=" << rho1 << "\n";
// define metric, gxx, alp, beta

        gxx[ijk] = (pow(x[ijk], 2)+pow(h * rr, 2) * (pow(y[ijk], 2)+pow(z[ijk], 2)))/pow(h * (pow(x[ijk], 2)+pow(y[ijk], 2)+pow(z[ijk], 2)), 2);
        gyy[ijk] = (pow(y[ijk], 2)+pow(h * rr, 2) * (pow(z[ijk], 2)+pow(x[ijk], 2)))/pow(h * (pow(x[ijk], 2)+pow(y[ijk], 2)+pow(z[ijk], 2)), 2);
        gzz[ijk] = (pow(z[ijk], 2)+pow(h * rr, 2) * (pow(y[ijk], 2)+pow(x[ijk], 2)))/pow(h * (pow(x[ijk], 2)+pow(y[ijk], 2)+pow(z[ijk], 2)), 2);
        gxy[ijk] = (( 1- pow(h*rr, 2) ) * x[ijk] * y[ijk])/pow(h * (pow(x[ijk], 2)+pow(y[ijk], 2)+pow(z[ijk], 2)), 2);
        gxz[ijk] = (( 1- pow(h*rr, 2) ) * x[ijk] * z[ijk])/pow(h * (pow(x[ijk], 2)+pow(y[ijk], 2)+pow(z[ijk], 2)), 2);
        gyz[ijk] = (( 1- pow(h*rr, 2) ) * z[ijk] * y[ijk])/pow(h * (pow(x[ijk], 2)+pow(y[ijk], 2)+pow(z[ijk], 2)), 2);

        alp[ijk] = sqrt(1-2 * param_m / rr + param_a);
cerr << "alp=" << alp[ijk] << "\n";
/*
do I need set beta=0?
 */
// define extrinsic kxx of 3+1 decomposition
            }}}
}

// remember to add it now make code routine, then run make.
