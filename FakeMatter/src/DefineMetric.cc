#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include "stdio.h"
#include "math.h"
using namespace std;

extern "C" void FakeMatter_DefineMetric(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;


// Loop over all grid points
for (int k = 0; k < cctk_lsh[2]; ++k) {
  for (int j = 0; j < cctk_lsh[1]; ++j) {
    for (int i = 0; i < cctk_lsh[0]; ++i) {

        int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
        double xx=x[ijk];
        double zz=z[ijk];
        double yy=y[ijk];
        double rho1=sqrt((pow(yy, 2)+pow(zz, 2)+pow(xx, 2)));
        double dd=1e-2;
        double psi4=pow(1+param_a/(2*rho1+dd),4);

        gxx[ijk]=psi4;gyy[ijk]=psi4;gzz[ijk]=psi4;
        alp[ijk] = (1-param_m / (2*rho1) + param_a) ;

  if(isnan(gzz[ijk])){printf("nan in gzz at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  if(isnan(gxx[ijk])){printf("nan in gxx at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  if(isnan(gyy[ijk])){printf("nan in gyy at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  if(isnan(gxz[ijk])){printf("nan in gxz at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  if(isnan(gxy[ijk])){printf("nan in gxy at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  if(isnan(gyz[ijk])){printf("nan in gyz at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  if(isnan(alp[ijk])){printf("nan in alp at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  
  if(!isfinite(gzz[ijk])){printf("inf/ind in gzz at %d,%d,%d,%g,%g,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk],rho1,psi4);}
  if(!isfinite(gxx[ijk])){printf("inf/ind in gxx at %d,%d,%d,%g,%g,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk],rho1,psi4);}
  if(!isfinite(gyy[ijk])){printf("inf/ind in gyy at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  if(!isfinite(gxz[ijk])){printf("inf/ind in gxz at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  if(!isfinite(gxy[ijk])){printf("inf/ind in gxy at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  if(!isfinite(gyz[ijk])){printf("inf/ind in gyz at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  if(!isfinite(alp[ijk])){printf("inf/ind in alp at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}

  if(x[ijk]==0&&y[ijk]==0&&z[ijk]==0){printf("metric gxx, gyy, gzz, gxy, gyz, gxz, alp, rho, psi4 at origin %d,%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk],gxx[ijk],gyy[ijk],gzz[ijk],gxy[ijk],gyz[ijk],gxz[ijk],alp[ijk],rho1,psi4);}

}
}
}
}

// remember to add it now make code routine, then run make.
