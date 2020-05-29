#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cmath>
#include <iostream>
using namespace std;




extern "C" void FakeMatter_AddMatter2(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

 for (int k = 0; k < cctk_lsh[2]; ++k) {
   for (int j = 0; j < cctk_lsh[1]; ++j) {
     for (int i = 0; i < cctk_lsh[0]; ++i) {
        int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

        double xx=x[ijk];
        double zz=z[ijk];
        double yy=y[ijk];
     
        double rho1 = sqrt(pow(xx,2)+pow(yy,2))+pow(zz,2);
     //   double rho2= sqrt(pow(xx,2)+pow(yy,2));
     //   if(rho2<1e-1)rho2=0.1;

  double g112= (pow(param_m,2)+4*(-1+2*param_a)*param_m*rho1+4*pow(rho1,2));
  double g111= - 8* param_turnon * param_a * param_m ;
  double g11=g111/g112;
  double g22= 4*param_turnon* param_a * param_m * rho1 * ( pow(param_m,2) + 4 * param_a * param_m * rho1 + 4* pow(rho1,2))/ pow( pow(param_m,2) + 4 * (-1+ 2 * param_a) * param_m *rho1+4 * pow(rho1,2),2) ;
          

          eTxx[ijk] += (g22*(pow(zz,2) + pow(yy,2))+g11*pow(xx,2)*rho1)/(pow(rho1,4)+0.01);
          eTyy[ijk] += (g22*(pow(xx,2) + pow(zz,2))+g11*pow(yy,2)*rho1)/(pow(rho1,4)+0.01);
          eTzz[ijk] += (g22*(pow(xx,2) + pow(yy,2))+g11*pow(zz,2)*rho1)/(pow(rho1,4)+0.01);
          eTxy[ijk] += xx *yy *(-g22 + g11*rho1)/(pow(rho1,4)+0.01);
          eTxz[ijk] += xx *zz *(-g22 + g11*rho1)/(pow(rho1,4)+0.01);
          eTyz[ijk] += zz *yy *(-g22 + g11*rho1)/(pow(rho1,4)+0.01) ;          
  /*
     
  double g112= rho1*(pow(param_m,2)+4*(-1+2*param_a)*param_m*rho1+4*pow(rho1,2));
  double g111= - 8* param_turnon * param_a * param_m ;
  double g11=g111/g112;
          double g22= 4*param_turnon* param_a * param_m * rho1 * ( pow(param_m,2) + 4 * param_a * param_m * rho1 + 4* pow(rho1,2))/ pow( pow(param_m,2) + 4 * (-1+ 2 * param_a) * param_m *rho1+4 * pow(rho1,2),2) ;
    if(rho1<1e-1)rho1=0.1;
          eTxx[ijk] += (g22*(pow(zz,2) + pow(yy,2)))/pow(rho1,4) + g11*pow(xx,2)/pow(rho1,2) ;
          eTyy[ijk] += (g22*(pow(xx,2) + pow(zz,2)))/pow(rho1,4) + g11*pow(yy,2)/pow(rho1,2) ;
          eTzz[ijk] += (g22*(pow(xx,2) + pow(yy,2)))/pow(rho1,4) + g11*pow(zz,2)/pow(rho1,2) ;
          eTxy[ijk] += xx *yy *(-g22 + g11*pow(rho1,2))/pow(rho1,4) ;
          eTxz[ijk] += xx *zz *(-g22 + g11*pow(rho1,2))/pow(rho1,4) ;
          eTyz[ijk] += zz *yy *(-g22 + g11*pow(rho1,2))/pow(rho1,4) ;


*/
/*
          double s1=sqrt(pow(xx,2)+pow(yy,2))/rho1; 
          double c1=zz/rho1;
          double s2=xx/rho2;
          double c2=yy/rho2; 
          double trans[4][4] = {{1,0,0,0},{0,s1*s2,c2*s1,c1},{0,c1*s2/rho1,c1*c2/rho1,-s1/rho1},{0,c2/(rho1),-s2/(rho1),0}};
          double Trho[4][4]  = {{0,0,0,0},{0,g11,0,0},{0,0,g22,0},{0,0,0,g22}};
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
  
*/

        
  //        eTxx[ijk] += 0;
  //        eTyy[ijk] += 0;
  //        eTzz[ijk] += 0;
  //        eTxy[ijk] += 0;
  //        eTxz[ijk] += 0;
  //        eTyz[iji] += 0; 
       



//if(x[ijk]==0&&y[ijk]==0&&z[ijk]==0){printf("at origin i%d,j%d,k%d,x%g,y%g,z%g,eTxx%g,eTyy%g,eTzz%g,eTxy%g,eTyz%g,eTxz%g,g111%g,g1112%g,g11%g,g22%g,rho%g\n",i,j,k,x[ijk],y[ijk],z[ijk],eTxx[ijk],eTyy[ijk],eTzz[ijk],eTxy[ijk],eTyz[ijk],eTxz[ijk],g111, g112,g11,g22,rho1);}


       }
    }
  }
}
