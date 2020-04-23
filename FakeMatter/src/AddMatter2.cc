#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cmath>
#include <iostream>
using namespace std;



// Matrix determinant
static CCTK_REAL determinant(const CCTK_REAL (&g)[3][3]) {
  return g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[2][1]) -
         g[1][0] * (g[0][1] * g[2][2] - g[0][2] * g[2][0]) +
         g[2][0] * (g[0][1] * g[1][2] - g[0][2] * g[1][0]);
}

// Matrix inverse
static void inverse(const CCTK_REAL (&g)[3][3], const CCTK_REAL detg,
                    CCTK_REAL (&gu)[3][3]) {
  gu[0][0] = (g[1][1] * g[2][2] - g[1][2] * g[2][1]) / detg;
  gu[1][0] = (g[0][1] * g[2][2] - g[0][2] * g[2][1]) / detg;
  gu[2][0] = (g[0][1] * g[1][2] - g[0][2] * g[1][1]) / detg;
  gu[1][1] = (g[0][0] * g[2][2] - g[0][2] * g[2][0]) / detg;
  gu[2][1] = (g[0][0] * g[1][2] - g[0][2] * g[1][0]) / detg;
  gu[2][2] = (g[0][0] * g[1][1] - g[0][1] * g[1][0]) / detg;
  gu[0][1] = gu[1][0];
  gu[0][2] = gu[2][0];
  gu[1][2] = gu[2][1];
}

// Length of a vector
static CCTK_REAL length(const CCTK_REAL (&g)[3][3], const CCTK_REAL (&x)[3]) {
  CCTK_REAL len = 0;
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
      len += g[a][b] * x[a] * x[b];
  len = sqrt(len);
  return len;
}

// Normalize a vector (to unit length)
static void normalize(CCTK_REAL (&x)[3], CCTK_REAL len) {
  for (int a = 0; a < 3; ++a)
    x[a] /= len;
}

extern "C" void FakeMatter_AddMatter2(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

 for (int k = 0; k < cctk_lsh[2]; ++k) {
   for (int j = 0; j < cctk_lsh[1]; ++j) {
     for (int i = 0; i < cctk_lsh[0]; ++i) {
        int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);


          const CCTK_REAL g[3][3] = {{gxx[ijk], gxy[ijk], gxz[ijk]},
                                     {gxy[ijk], gyy[ijk], gyz[ijk]},
                                     {gxz[ijk], gyz[ijk], gzz[ijk]}};
          // Metric determinant
          const CCTK_REAL detg = determinant(g);
          // Inverse metric
          CCTK_REAL gu[3][3];
          inverse(g, detg, gu);
	  
        double rho1 = r[ijk];
        if(rho1<1e-1)rho1=0.1;

        double xx=x[ijk];
        double zz=z[ijk];
        double yy=y[ijk];


          double g11= - 8* param_turnon * param_a * param_m /  rho1*(pow(param_m,2)+4*(-1+2*param_a)*param_m*rho1+4*pow(rho1,2));
          double g22= 4*param_turnon* param_a * param_m * rho1 * ( pow(param_m,2) + 4 * param_a * param_m * rho1 + 4* pow(rho1,2))/ pow( pow(param_m,2) + 4 * (-1+ 2 * param_a) * param_m *rho1+4 * pow(rho1,2),2) ;
          
          double s1=sqrt(pow(xx,2)+pow(yy,2))/sqrt(pow(xx,2)+pow(yy,2)+pow(zz,2)); 
          double c1=zz/sqrt(pow(xx,2)+pow(yy,2)+pow(zz,2));
          double s2=xx/sqrt(pow(xx,2)+pow(yy,2));
          double c2=yy/sqrt(pow(xx,2)+pow(yy,2)); 
 if(s1<1e-3)s1=1e-3;
          double trans[4][4] = {{1,0,0,0},{0,s1*s2,c2*s1,c1},{0,c1*s2/rho1,c1*c2/rho1,-s1/rho1},{0,c2/(rho1*s1),-s2/(rho1*s1),0}};
          double Trho[4][4]  = {{0,0,0,0},{0,g11,0,0},{0,0,g22,0},{0,0,0,g22*pow(s1,2)}};
          double Txyz[4][4];
    
          for (int a = 0; a<4; ++a) 
              for (int b = 0; b<4; ++b)
          {
                for (int c = 0; c<4; ++c) 
                    for (int d = 0; d<4; ++d)
                      Txyz[a][b] += trans[c][a]* Trho[c][d]* trans[d][b];
          }
          
          
          eTxx[ijk] += Txyz[1][1];
          eTyy[ijk] += Txyz[2][2];
          eTzz[ijk] += Txyz[3][3];
          eTxy[ijk] += Txyz[1][2];
          eTxz[ijk] += Txyz[1][3];
          eTyz[ijk] += Txyz[2][3]; 
       
   if(isnan(eTxx[ijk])){printf("nan in eTxx at %g,%g,%g,%g,%g,%g,%g\n",x[ijk],y[ijk],z[ijk], s1,s2,c1,c2);}

       }
    }
  }
}
