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
/*return pow(sqrt(1-2*param_m/r)*r + 2 * param_m * atanh(sqrt(1-2* param_m /r))+param_c,2)-pow(rho,2);*/
double res=1;
if ( rho > param_c ) { 
res=sqrt(1-2*param_m/r)*r + 2 * param_m * atanh(sqrt(1-2* param_m /r))+param_c-rho;}
else {
//cerr << "r=" << rr << "\n";
 //cerr << "rho=" << rho1 << "\n";
//cerr << "you enter here";
res=sqrt(1-2*param_m/r)*r + 2 * param_m * atanh(sqrt(1-2* param_m /r))-param_c+rho;
//cerr << "res" << res << "\n"; 
}
return res;
}

double solve(double rho)
{
DECLARE_CCTK_PARAMETERS;
double low=2 * param_m ,up=10,mid=(low+up)/2;
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

extern "C" void FakeMatter_DefineMetric(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

// Loop over all grid points
for (int k = 0; k < cctk_lsh[2]; ++k) {
  for (int j = 0; j < cctk_lsh[1]; ++j) {
    for (int i = 0; i < cctk_lsh[0]; ++i) {
        int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
        double rho1 = r[ijk];
        if(rho1<1e-2) rho1=1e-2;
        double rr = solve(rho1);
// define metric, gxx, alp, beta

        gxx[ijk] = (pow(x[ijk], 4)+(pow(x[ijk], 2)+pow( rr, 2)) * (pow(y[ijk], 2)+pow(z[ijk], 2)))/pow( rho1, 2);
        gyy[ijk] = (pow(y[ijk], 4)+(pow(y[ijk], 2)+pow( rr, 2)) * (pow(z[ijk], 2)+pow(x[ijk], 2)))/pow( rho1, 2);
        gzz[ijk] = (pow(z[ijk], 4)+(pow(z[ijk], 2)+pow( rr, 2) )* (pow(y[ijk], 2)+pow(x[ijk], 2)))/pow( rho1, 2);
        gxy[ijk] = (( (pow(x[ijk], 2)+pow(y[ijk], 2)+pow(z[ijk], 2)) - pow(rr, 2) ) * x[ijk] * y[ijk])/pow( rho1, 2);
        gxz[ijk] = (( (pow(x[ijk], 2)+pow(y[ijk], 2)+pow(z[ijk], 2)) - pow(rr, 2) ) * x[ijk] * z[ijk])/pow( rho1, 2);
        gyz[ijk] = (( (pow(x[ijk], 2)+pow(y[ijk], 2)+pow(z[ijk], 2)) - pow(rr, 2) ) * z[ijk] * y[ijk])/pow( rho1, 2);

        alp[ijk] = sqrt(1-2 * param_m / rr + param_a);
  if(isnan(gzz[ijk])){printf("nan in gzz at %d,%d,%d,%g,%g,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk],rho1,rr);}
  if(isnan(gxx[ijk])){printf("nan in gxx at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  if(isnan(gyy[ijk])){printf("nan in gyy at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  if(isnan(gxz[ijk])){printf("nan in gxz at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  if(isnan(gxy[ijk])){printf("nan in gxy at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  if(isnan(gyz[ijk])){printf("nan in gyz at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  if(isnan(alp[ijk])){printf("nan in alp at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  
  if(isinf(gzz[ijk])){printf("inf in gzz at %d,%d,%d,%g,%g,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk],rho1,rr);}
  if(isinf(gxx[ijk])){printf("inf in gxx at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  if(isinf(gyy[ijk])){printf("inf in gyy at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  if(isinf(gxz[ijk])){printf("inf in gxz at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  if(isinf(gxy[ijk])){printf("inf in gxy at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  if(isinf(gyz[ijk])){printf("inf in gyz at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  if(isinf(alp[ijk])){printf("inf in alp at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}


  if(x[ijk]==0&&y[ijk]==0&&z[ijk]==0){printf("metric gxx, gyy, gzz, gxy, gyz, gxz, alp at origin %d,%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk],gxx[ijk],gyy[ijk],gzz[ijk],gxy[ijk],gyz[ijk],gxz[ijk],alp[ijk]);};

 if(x[ijk]<1e-2&&y[ijk]<1e-2&&z[ijk]<1e-2)alp[ijk]=1e-2;


  if(gzz[ijk]<1e-2)gzz[ijk]=1e-2;
  if(gxx[ijk]<1e-2)gxx[ijk]=1e-2;
  if(gyy[ijk]<1e-2)gyy[ijk]=1e-2;
  if(gzz[ijk]<1e-2){printf("0 in gzz at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  if(gxx[ijk]<1e-2){printf("0 in gxx at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  if(gyy[ijk]<1e-2){printf("0 in gyy at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}

  if(alp[ijk]==0){printf("0 in alp at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}

  if(x[ijk]==0&&y[ijk]==0&&z[ijk]==0){printf("metric gxx, gyy, gzz, gxy, gyz, gxz, alp at origin after smooth %d,%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk],gxx[ijk],gyy[ijk],gzz[ijk],gxy[ijk],gyz[ijk],gxz[ijk],alp[ijk]);};

// define extrinsic kxx of 3+1 decomposition
if (rho1 < param_c) 
{
double res1=-sqrt(1-2*param_m/rr)*rr - 2 * param_m * atanh(sqrt(1-2* param_m /rr))+param_c;
//cerr << "rho from r" << res1 << "\n";
//cerr << "rho" << rho1 << "\n";
//cerr << "r" << rr << "\n";
//cerr << "alp" << alp[ijk] << "\n";
}
            }}}
}

// remember to add it now make code routine, then run make.
