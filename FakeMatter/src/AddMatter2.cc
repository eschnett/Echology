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
 // return;

/*
 // Horizon surface
  const int hi = horizon_index;
  if (!sf_active[hi]) {
    CCTK_WARN(CCTK_WARN_ALERT, "No horizon found at this time");
    return;
  }
  // Origin (centre) of horizon
  const CCTK_REAL hc[3] = {sf_origin_x[hi], sf_origin_y[hi], sf_origin_z[hi]};
  // Constant expansion surface
  const int si = surface_index;
  if (!sf_active[si]) {
    CCTK_WARN(CCTK_WARN_ALERT,
              "No constant expansion surface found at this time");
    return;
  }
  // Origin (centre) of constant expansion surface
  const CCTK_REAL sc[3] = {sf_origin_x[si], sf_origin_y[si], sf_origin_z[si]};
*/

 for (int k = 0; k < cctk_lsh[2]; ++k) {
   for (int j = 0; j < cctk_lsh[1]; ++j) {
     for (int i = 0; i < cctk_lsh[0]; ++i) {
        int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
// cerr << smask[ijk] << endl;
// Are we inside the constant expansion surface?

// since the horizon finder is scheduled later, we dont send matter until it 
        if (smask[ijk] == 1 && cctk_iteration > 0 ) {
//here?
        double xx=x[ijk];
        double zz=z[ijk];
        double yy=y[ijk];
     
        double rho1 = sqrt(pow(xx,2)+pow(yy,2)+pow(zz,2));
 // add expansion
  double rho2 = surface_expansion * pow(param_m,2) / sqrt(2) + param_m / 2; 
       
  double g112= rho2 * (pow(param_m,2)+4*(-1+2*param_a)*param_m*rho2+4*pow(rho2,2));
  double g111= - 8* param_turnon * param_a * param_m ;
  double g11=g111/g112;
  double g22= 4*param_turnon* param_a * param_m * rho2 * ( pow(param_m,2) + 4 * param_a * param_m * rho2 + 4* pow(rho2,2))/ pow( pow(param_m,2) + 4 * (-1 + 2 * param_a) * param_m * rho2 + 4 * pow(rho2,2),2) ;

  double co = rho2 / 0.1 * rho1;
//expansion I forgot.... co should be the linear interopolation.. anyway I will check it later
         
//i want to check the grammar for assert statement 
if(rho1 == 0){printf("rho shouldn't be zero in the loop!"); assert(0);}
          eTxx[ijk] += co * (g22*(pow(zz,2) + pow(yy,2))+g11*pow(xx,2)* pow(rho1,2))/(pow(rho1,4));
          eTyy[ijk] += co * (g22*(pow(xx,2) + pow(zz,2))+g11*pow(yy,2)* pow(rho1,2))/(pow(rho1,4));
          eTzz[ijk] += co * (g22*(pow(xx,2) + pow(yy,2))+g11*pow(zz,2)* pow(rho1,2))/(pow(rho1,4));
          eTxy[ijk] += co * xx * yy * (-g22 + g11*pow(rho1,2))/(pow(rho1,4));
          eTxz[ijk] += co * xx * zz * (-g22 + g11*pow(rho1,2))/(pow(rho1,4));
          eTyz[ijk] += co * zz * yy * (-g22 + g11*pow(rho1,2))/(pow(rho1,4)) ;
 if(eTxx[ijk] > 10e5){printf("Txx > 10^5 at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
 if(eTyy[ijk] > 10e5){printf("Tyy > 10^5 at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  if(eTzz[ijk] > 10e5){printf("Tzz > 10^5 at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  if(eTxz[ijk] > 10e5){printf("Txz > 10^5 at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
 if(eTxy[ijk] > 10e5){printf("Txy > 10^5 at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
 if(eTyz[ijk] > 10e5){printf("Tyz > 10^5 at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
 if(co > 10e5){printf("co > 10^5 at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
 if(g11 > 10e5){printf("g11 > 10^5 at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
 if(g22 > 10e5){printf("g22 > 10^5 at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
 if(rho1 > 10e5){printf("rho1 > 10^5 at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
 if(rho2 > 10e5){printf("rho2 > 10^5 at %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}

 if(!isfinite(eTxx[ijk])){printf("Txx nanat %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
 if(!isfinite(eTyy[ijk])){printf("Tyy nanat %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  if(!isfinite(eTzz[ijk])){printf("Tzz nanat %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
  if(!isfinite(eTxz[ijk])){printf("Txz nanat %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
 if(!isfinite(eTxy[ijk])){printf("Txy nanat %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
 if(!isfinite(eTyz[ijk])){printf("Tyz nanat %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
 if(!isfinite(co)){printf("co nanat %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
 if(!isfinite(g11)){printf("g11 nanat %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
 if(!isfinite(g22)){printf("g22 nanat %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
 if(!isfinite(rho1)){printf("rho1 nanat %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
 if(!isfinite(rho2)){printf("rho2 nanat %d,%d,%d,%g,%g,%g\n",i,j,k,x[ijk],y[ijk],z[ijk]);}
//cerr << "Txx, co, rho2, rho1 " << eTxx << co << rho2 << rho1  << endl;
/*  
  double g112= (pow(param_m,2)+4*(-1+2*param_a)*param_m*rho1+4*pow(rho1,2));
  double g111= - 8* param_turnon * param_a * param_m ;
  double g11=g111/g112;
  double g22= 4*param_turnon* param_a * param_m * rho1 * ( pow(param_m,2) + 4 * param_a * param_m * rho1 + 4* pow(rho1,2))/ pow( pow(param_m,2) + 4 * (-1+ 2 * param_a) * param_m *rho1+4 * pow(rho1,2),2) ;
*/

//smooth 3 set a cut off
/*
          eTxx[ijk] += (g22*(pow(zz,2) + pow(yy,2))+g11*pow(xx,2)*rho1)/(pow(rho1,4));
          eTyy[ijk] += (g22*(pow(xx,2) + pow(zz,2))+g11*pow(yy,2)*rho1)/(pow(rho1,4));
          eTzz[ijk] += (g22*(pow(xx,2) + pow(yy,2))+g11*pow(zz,2)*rho1)/(pow(rho1,4));
          eTxy[ijk] += xx *yy *(-g22 + g11*rho1)/(pow(rho1,4));
          eTxz[ijk] += xx *zz *(-g22 + g11*rho1)/(pow(rho1,4));
          eTyz[ijk] += zz *yy *(-g22 + g11*rho1)/(pow(rho1,4)) ;
  */        
//smooth 2 give a dd in the denominantor

/*
          eTxx[ijk] += (g22*(pow(zz,2) + pow(yy,2))+g11*pow(xx,2)*rho1)/(pow(rho1,4)+0.01);
          eTyy[ijk] += (g22*(pow(xx,2) + pow(zz,2))+g11*pow(yy,2)*rho1)/(pow(rho1,4)+0.01);
          eTzz[ijk] += (g22*(pow(xx,2) + pow(yy,2))+g11*pow(zz,2)*rho1)/(pow(rho1,4)+0.01);
          eTxy[ijk] += xx *yy *(-g22 + g11*rho1)/(pow(rho1,4)+0.01);
          eTxz[ijk] += xx *zz *(-g22 + g11*rho1)/(pow(rho1,4)+0.01);
          eTyz[ijk] += zz *yy *(-g22 + g11*rho1)/(pow(rho1,4)+0.01) ;          

*/

// smooth 1 set a mim for rho
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
}
