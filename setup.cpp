/*!  \file setup.cpp
 *   \brief Sets up variables specifying the initial conditions.
 *
 *   Fixed parameters, in contrast, are specified in constants.h.
 *   Users should mostly only need to alter setup.cpp and constants.h
 *
 *   \author Julian Pittard (Original version 02.06.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 16.06.11 (JMP)
 */

#undef MAIN

//  Header files:
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "constants.h"
#include "gas.h"
#include "global.h"
#include "grid.h"
#include "physconst.h"
#include "star.h"

using namespace std;

// Function declarations
void cellCentreOfMass();
void checkSetup();
double findForwardShock();
void grid(int nzones, double xmin, double xmax, double xa[], double xc[], double dx[]);
void gridSetup();
void mapWind(bool init);
void setupCWB();
void setupSOD();
void setupSNR();
void setupWBB();
extern void computeTimestep();
extern void quit();

/*!  \brief Sets up variables specifying the initial conditions.
 *
 *   Fixed parameters, in contrast, are specified in constants.h.
 *   Users should mostly only need to alter setup.cpp and constants.h
 *
 *   \author Julian Pittard (Original version 02.06.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 02.06.11 (JMP)
 */
void setup(){
  //////////////////////////////////////////////////////////////////////
  //                                                                  //
  //                   SET UP CODE VARIABLES                          //
  // ...students can (carefully) alter the values specified below...  //
  //                                                                  //
  //       NOTE: Remember to recompile after making changes!          // 
  //                                                                  //
  //////////////////////////////////////////////////////////////////////

  /*
  problem = "SOD";

  lg.geom = "XYZ";
  lg.ngeomx = 0;    
  lg.ngeomy = 0;
  lg.ngeomz = 0;

  lg.nleftx = 0; // reflecting
  lg.nlefty = 0;
  lg.nleftz = 0;

  lg.nrightx = 0;
  lg.nrighty = 0;
  lg.nrightz = 0;

  lg.xmin = 0.0;
  lg.xmax = 1.0;
  lg.ymin = 0.0;
  lg.ymax = 1.0;
  lg.zmin = 0.0;
  lg.zmax = 1.0;
  */
  
  
  //problem = "SNR";
  //problem = "WBB";
  problem = "CWB";
  lg.geom = "ZRP";
  lg.ngeomx = 0;    
  lg.ngeomy = 1;
  lg.ngeomz = 3;

  lg.nleftx = 1; 
  lg.nlefty = 0;
  lg.nleftz = 0;

  lg.nrightx = 1;
  lg.nrighty = 1;
  lg.nrightz = 0;

  lg.xmin = -1.0e12;
  lg.xmax = 1.0e12;
  lg.ymin = 0.0;
  lg.ymax = 2.0e12;
  lg.zmin = 0.0;
  lg.zmax = 1.0;
  
  
  ncycle = 0;
  ncycend = 500000;
  nprin1d = 1000000;
  nprin2d = 100;
  endtime = 1.0e4;
  ncycp1d = 0;
  ncycp2d = 0;
  ncycp3d = 0;
  timep1d = 0.0;
  timep2d = 0.0;
  timep3d = 0.0;
  tprin1d = 1.0e99;
  tprin2d = 1.0e99;
  tprin3d = 1.0e99;
  
  t = 0.0;
  svel = 1.0;

  // Check inconsistencies and setup grid...
  checkSetup();
  gridSetup();

  if (reStart == false){
    if (problem == "SOD"){
      setupSOD();
    }
    else if (problem == "SNR"){
      setupSNR();
    }
    else if (problem == "WBB"){
      setupWBB();
    }
    else if (problem == "CWB"){
      setupCWB();
    }
    else{
      cout << "Unknown problem in setup!\n";
      quit();
    }
  }
  else{
    //setupRestart();
  }
    
  computeTimestep();
  cout << "Initial time = " << t << "; dt = " << dt << "; endtime = " << endtime << "\n";
  //for (int i = 0; i < maxsweep; ++i){
  //  cout << i << " " << zparax[1][i] << " " << zparax[2][i] << " " << zparax[3][i] << " " << zparax[4][i] << " " << zparax[5][i] << "\n";
  //}
  //quit();
  return;
}


/*!  \brief Check that the setup is consistent.
 *
 *   \author Julian Pittard (Original version 05.06.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 10.08.11 (JMP)
 */
void checkSetup(){

  if ((nd == 1) && ((jmax != 1)||(kmax != 1))){
    cout << "nd = " << nd << ", jmax = " << jmax << ", kmax = " << kmax << "\n";
    cout << "Aborting...\n";
    exit(EXIT_FAILURE);
  }
  if ((nd != 3) && (kmax != 1)){
    cout << "nd = " << nd << ", kmax = " << kmax << "\n";
    cout << "Aborting...\n";
    exit(EXIT_FAILURE);
  }

  if ((lg.ngeomx == 0) && (lg.ngeomy == 0) && (lg.ngeomz == 0)) lg.geom = "XYZ";
  else if ((lg.ngeomx == 0) && (lg.ngeomy == 1) && (lg.ngeomz == 3)) lg.geom = "ZRP";
  else if ((lg.ngeomx == 2) && (lg.ngeomy == 4) && (lg.ngeomz == 5)) lg.geom = "RTP";
  else{
    if (procRank == 0){
      cout << "Error: geometry unrecognized!\n";
      cout << " lg.ngeomx = " << lg.ngeomx << "\n";
      cout << " lg.ngeomy = " << lg.ngeomy << "\n";
      cout << " lg.ngeomz = " << lg.ngeomz << "\n";
    }
    quit();
  }
  
  if ((nd == 3) && (lg.geom != "XYZ")){
    if (procRank == 0) cout << "Only XYZ geometry is currently supported in 3D. Aborting!\n";
    quit();
  }

  if ((nd == 1 && cn > 0.4)||(nd == 2 && cn > 0.25)||(nd == 3 && cn > 0.2)){
    if (procRank == 0){
      cout << "courant number is too high (cn = " << cn << ")!\n";
      cout << " Maximum permissible values are:\n";
      cout << "   Eulerian computation - 1D: cn = 0.4; 2D: cn = 0.25; 3D: cn = 0.2\n";
    }
    quit();
  }

  return;
}


/*!  \brief Setup grid coordinates and calculate parabolic coefficients.
 *
 *   \author Julian Pittard (Original version 02.06.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 23.08.11 (JMP)
 */
void gridSetup(){
  ncells = imax*jmax*kmax;

  // Set range for real cells on the grid (irs->ire inclusive)
  lg.irs = nghost;
  lg.ire = imax+nghost-1;
  if (nd > 1){
    lg.jrs = nghost;
    lg.jre = jmax+nghost-1;
  }
  else{
    lg.jrs = 0;
    lg.jre = 0;
  }
  if (nd == 3){
    lg.krs = nghost;
    lg.kre = kmax+nghost-1;
  }
  else{
    lg.krs = 0;
    lg.kre = 0;
  }
    
  grid(imax,lg.xmin,lg.xmax,lg.zxa,lg.zxc,lg.zdx);
  grid(jmax,lg.ymin,lg.ymax,lg.zya,lg.zyc,lg.zdy);
  grid(kmax,lg.zmin,lg.zmax,lg.zza,lg.zzc,lg.zdz);
  cellCentreOfMass();
  return;
}


/*!  \brief Function to create coords for the grid in a given direction.
 *
 *   \param xmin (double) is the left edge of the grid.
 *   \param xmax (double) is the right edge of the grid.
 *   \param xa (double) returns the left edge coord of each grid cell.
 *   \param xa (double) returns the coord of the centre of each grid cell.
 *   \param dx (double) returns the width of each grid cell.
 *
 *   xa(girs) is left boundary location (at xmin). xa(gire+1) is right boundary 
 *   location (at xmax).
 *
 *   \author Julian Pittard (Original version 31.05.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 04.08.11 (JMP)
 */
void grid(int nzones, double xmin, double xmax, double xa[], double xc[], double dx[])
{
  double dxfac = (xmax - xmin) / double(nzones);
  int ntot = nzones + 2*nghost;

  for (int n = 0; n <= ntot; n++){
    xa[n] = xmin + double(n-nghost)*dxfac;
    dx[n] = dxfac;
    xc[n] = xa[n] + 0.5*dx[n];
  }
  return;
}



/*!  \brief Calculate the cell centre of mass arrays. Based on centis in mg.
 *
 *   \author Julian Pittard (Original version 02.06.11)
 *   \version 0.1-development (Undomiel):
 *   \date Last modified: 21.09.11 (JMP)
 */
void cellCentreOfMass(){

  double x1,x2,y1,y2,y1abs,y2abs;
  
  if (lg.geom == "XYZ"){
    for (int i = 0; i <= imaxd; ++i){ 
      lg.zxg[i] = 0.5*(lg.zxa[i] + lg.zxa[i+1]);
    }
    if (nd > 1){
      for (int j = 0; j <= jmaxd; ++j){ 
	lg.zyg[j] = 0.5*(lg.zya[j] + lg.zya[j+1]);
      }
    }
    if (nd == 3){
      for (int k = 0; k <= kmaxd; ++k){ 
	lg.zzg[k] = 0.5*(lg.zza[k] + lg.zza[k+1]);
      }
    }
  }
  else if (lg.geom == "ZRP"){
    for (int i = 0; i <= imaxd; ++i){ 
      lg.zxg[i] = 0.5*(lg.zxa[i] + lg.zxa[i+1]);
    }
    if (nd > 1){
      for (int j = 0; j <= jmaxd; ++j){ 
	y1 = lg.zya[j];
	y2 = lg.zya[j+1];
	if (y1 >= 0.0){
	  lg.zyg[j] = (2.0/3.0)*(y2*y2 + y1*y2 + y1*y1)/(y1+y2);
	}
	else{
	  y1abs = abs(y2);
	  y2abs = abs(y1);
	  lg.zyg[j] = -(2.0/3.0)*(y2abs*y2abs + y1abs*y2abs + y1abs*y1abs)/(y1abs+y2abs);
	}
      }
    }
    if (nd == 3){
      for (int k = 0; k <= kmaxd; ++k){ 
	lg.zzg[k] = 0.5*(lg.zza[k] + lg.zza[k+1]);
      }
    }
  }
  else if (lg.geom == "RTP"){
    for (int i = 0; i <= imaxd; ++i){ 
      x1 = lg.zxa[i];
      x2 = lg.zxa[i+1];
      lg.zxg[i] = 0.75*(x2*x2*x2 + x1*x2*(x2 + x1) + x1*x1*x1)/(x2*x2 + x1*x2 + x1*x1);
    }
    if (nd > 1){
      for (int j = 0; j <= jmaxd; ++j){ 
	y1 = lg.zya[j];
	y2 = lg.zya[j+1];
	lg.zyg[j] = (y1*cos(y1) - y2*cos(y2) + sin(y2) - sin(y1))/(cos(y1) - cos(y2));
      }
    }
    if (nd == 3){
      for (int k = 0; k <= kmaxd; ++k){ 
	lg.zzg[k] = 0.5*(lg.zza[k] + lg.zza[k+1]);
      }
    }
  }
  else{
    if (procRank == 0) cout << "cellCentreOfMass: error - not coded for geometry " << lg.geom << ". Aborting!\n";
    quit();
  }

  return;
}



/*!  \brief Sets up a SOD shock tube problem for the initial conditions.
 *
 *   \author Julian Pittard (Original version 02.06.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 21.09.11 (JMP)
 */
void setupSOD(){

  prefx = "sod";
  
  //tmin = 1.0e4;
  //tmax = 1.0e9;
  
  // Set up ambient medium and ejecta
  for (int k = lg.krs; k <= lg.kre; ++k){
    for (int j = lg.jrs; j <= lg.jre; ++j){
      for (int i = lg.irs; i <= lg.ire; ++i){
	lg.P0[iqd][k][j][i] = 1.0;
        lg.P0[iqe][k][j][i] = 1.0;
	for (int n = 0; n < nd; ++n){
	  lg.P0[iqu0+n][k][j][i] = 0.0;
	}
	for (int n = 0; n < smax; ++n){
	  lg.P0[iqal0+n][k][j][i] = 0.0;
	}      
	if (nd == 1){
	  if (i < imax/2 + nghost){
  	    lg.P0[iqe][k][j][i] = 10.0;
	    lg.P0[iqal0][k][j][i] = 1.0;
	  }
	}
	else if (nd == 2){
	  int img = i - nghost;
	  int jmg = j - nghost;
	  if (img*img + jmg*jmg <= 0.5*imax*imax){
  	    lg.P0[iqe][k][j][i] = 10.0;
	    lg.P0[iqal0][k][j][i] = 1.0;
	  }
	}
      }
    }
  }
  
  return;
}


/*!  \brief Sets up a supernova remnant problem for the initial conditions.
 *
 *   \author Julian Pittard (Original version 02.06.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 21.09.11 (JMP)
 */
void setupSNR(){

  // 1D spherically symmetric
  //lg.geom = "RTP";
  //lg.ngeomx = 2;    
  //lg.ngeomy = 4;
  //lg.ngeomz = 5;

  // 2D ZRP
  lg.geom = "ZRP";
  lg.ngeomx = 0;    
  lg.ngeomy = 1;
  lg.ngeomz = 3;

  lg.nleftx = 0;
  lg.nlefty = 0;
  lg.nleftz = 0;

  lg.nrightx = 1;
  lg.nrighty = 1;
  lg.nrightz = 1;

  // Grid extent
  lg.xmin = 0.0;
  lg.xmax = 15.0*pc;
  lg.ymin = 0.0;
  //lg.ymax = 1.0;
  lg.ymax = 15.0*pc;
  lg.zmin = 0.0;
  lg.zmax = 1.0;

  ncycle = 0;
  ncycend = 100;
  nprin1d = 500;
  nprin2d = 10;
  endtime = 1.0e20;
  timep1d = 0.0;
  timep2d = 0.0;
  timep3d = 0.0;
  t = 0.0;
  svel = 1.0;
    
  tmin = 1.0e4;
  tmax = 1.0e9;

  esn = 1.0e51;      // ergs
  msn = 10.0;        // Msol
  ramb = 1.0e-24;    // g/cm^-3
  tamb = 1.0e4;      // K
  uamb = 0.0;
  avgmassamb = 1.0e-24; // g
  pamb = (ramb/avgmassamb)*boltzman*tamb;
    
  // Check inconsistencies and setup grid...
  checkSetup();
  gridSetup();

  remapRadius = 6.0*lg.zdx[0];
  prefx = "snr";

  double vol = (4.0/3.0)*pi*pow(remapRadius,3);
  double drho = msn*msol/vol;
  double dp = gamm*esn/vol;
  
  // Set up ambient medium and ejecta. For ZRP geometry, j is R and i is Z.
  for (int k = lg.krs; k <= lg.kre; ++k){
    double zc = lg.zzc[k];
    double zc2 = zc*zc;
    for (int j = lg.jrs; j <= lg.jre; ++j){
      double yc = lg.zyc[j];
      double yc2 = yc*yc;
      for (int i = lg.irs; i <= lg.ire; ++i){
        double xc = lg.zxc[i];
        double xc2 = xc*xc;
	double r = sqrt(xc2 + yc2 + zc2);
	lg.P0[iqd][k][j][i] = ramb;
        lg.P0[iqe][k][j][i] = pamb;
	for (int n = 0; n < nd; ++n){
	  lg.P0[iqu0+n][k][j][i] = 0.0;
	}
	for (int n = 0; n < smax; ++n){
	  lg.P0[iqal0+n][k][j][i] = 0.0;
	}
	if (r < remapRadius){
	  lg.P0[iqd][k][j][i] += drho;
	  lg.P0[iqe][k][j][i] += dp;
  	  lg.P0[iqal0][k][j][i] = 1.0;
	}
      }
    }
  }
  
  return;
}


/*!  \brief Sets up a wind-blown-bubble problem for the initial conditions.
 *
 *   \author Julian Pittard (Original version 02.06.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 21.09.11 (JMP)
 */
void setupWBB(){

  lg.geom = "RTP";
  lg.ngeomx = 2;    
  lg.ngeomy = 4;
  lg.ngeomz = 5;

  lg.nleftx = 0;
  lg.nlefty = 0;
  lg.nleftz = 0;

  lg.nrightx = 1;
  lg.nrighty = 1;
  lg.nrightz = 1;

  // Grid extent
  lg.xmin = 0.0;
  lg.xmax = 15.0*pc;
  lg.ymin = 0.0;
  lg.ymax = 1.0;
  lg.zmin = 0.0;
  lg.zmax = 1.0;

  ncycle = 0;
  ncycend = 5000;
  nprin1d = 500;
  endtime = 1.0e20;
  timep1d = 0.0;
  timep2d = 0.0;
  timep3d = 0.0;
  t = 0.0;
  svel = 1.0;
    
  tmin = 1.0e4;
  tmax = 1.0e9;

  stars[0].mdot = 1.0e-6*msol/yr;
  stars[0].vinf = 2.0e8;
  stars[0].xpos = 0.0;
  stars[0].ypos = 0.0;
  stars[0].zpos = 0.0;

  ramb = 1.0e-26;    // g/cm^-3
  tamb = 1.0e4;      // K
  uamb = 0.0;
  avgmassamb = 1.0e-24; // g
  pamb = (ramb/avgmassamb)*boltzman*tamb;
    
  // Check inconsistencies and setup grid...
  checkSetup();
  gridSetup();

  remapRadius = 6.0*lg.zdx[0];
  prefx = "wbb";
       
  // Set up ambient medium. For ZRP geometry, j is R and i is Z.
  for (int k = lg.krs; k <= lg.kre; ++k){
    for (int j = lg.jrs; j <= lg.jre; ++j){
      for (int i = lg.irs; i <= lg.ire; ++i){
	lg.P0[iqd][k][j][i] = ramb;
        lg.P0[iqe][k][j][i] = pamb;
	for (int n = 0; n < nd; ++n){
	  lg.P0[iqu0+n][k][j][i] = 0.0;
	}
	for (int n = 0; n < smax; ++n){
	  lg.P0[iqal0+n][k][j][i] = 0.0;
	}
      }
    }
  }

  // Set up wind
  mapWind(false);
  
  return;
}


/*!  \brief Set up the CWB problem.
 *
 *   \author Julian Pittard (Original version 02.06.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 21.09.11 (JMP)
 */
void setupCWB(){

  t = 0.0;
  svel = 1.0;
    
  tmin = 1.0e4;
  tmax = 1.0e9;

  stars[0].mdot = 1.0e-6*msol/yr;
  stars[0].vinf = 2.0e8;
  stars[0].xpos = -5.0e11;
  stars[0].ypos = 0.0;
  stars[0].zpos = 0.0;

  stars[1].mdot = 1.0e-6*msol/yr;
  stars[1].vinf = 2.0e8;
  stars[1].xpos = 5.0e11;
  stars[1].ypos = 0.0;
  stars[1].zpos = 0.0;

  tamb = 1.0e4;         // K
  avgmassamb = 1.0e-24; // g
  remapi = 10;
    
  remapRadius = remapi*lg.zdx[0];

  dsep = sqrt(pow(stars[1].xpos - stars[0].xpos,2) + pow(stars[1].ypos - stars[0].ypos,2) + pow(stars[1].zpos - stars[0].zpos,2));
  eta = stars[1].mdot*stars[1].vinf/(stars[0].mdot*stars[0].vinf);        // wind mtm ratio
  rob = (sqrt(eta)/(1.0 + sqrt(eta)))*dsep;//distance of stagnation point from star 1 (distance from star 0 is rwr)

  // Check that the WCR has enough room not to interfere with the remap regions
  if (int(rob/lg.zdx[0]) < remapi + 10){
    cout << "setupCWB: WCR too close to star 1. Aborting!";
    quit();
  }
  if (int((dsep - rob)/lg.zdx[0]) < remapi + 10){
    cout << "setupCWB: WCR too close to star 0. Aborting!";
    quit();
  }
  
  // Check inconsistencies and setup grid...
  checkSetup();
  gridSetup();

  prefx = "cwb";

  // Initialize winds
  mapWind(true);
    
  return;
}


/*!  \brief Map a wind solution onto the grid.
 *
 *   \author Julian Pittard (Original version 05.06.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 10.08.11 (JMP)
 */
void mapWind(bool init){

  double xc, yc, zc, xc2, yc2, zc2, r2, r, xy;
  double sinphi, cosphi, costhta, sinthta;
  double xmaxst1;
  double mdot1 = stars[0].mdot;
  double vinf1 = stars[0].vinf;
  double xpos1 = stars[0].xpos;
  double ypos1 = stars[0].ypos;
  double zpos1 = stars[0].zpos;
  double mdot2,vinf2,xpos2,ypos2,zpos2;
  double mdot,vinf,xpos,ypos,zpos;
  int istar,jstar,kstar,istl,istu,jstl,jstu,kstl,kstu;

  if (nstars > 1){
    mdot2 = stars[1].mdot;
    vinf2 = stars[1].vinf;
    xpos2 = stars[1].xpos;
    ypos2 = stars[1].ypos;
    zpos2 = stars[1].zpos;
  }

  if (!init){
    for (int n = 0; n < nstars; ++n){
      mdot = stars[n].mdot;
      vinf = stars[n].vinf;
      xpos = stars[n].xpos;
      ypos = stars[n].ypos;
      zpos = stars[n].zpos;
      // Determine extent of remap region.
      // lg.irs is the index of the first real cell, lg.ire is the last real cell
      istar = int((xpos - lg.zxa[0])/lg.zdx[0]);
      jstar = 0;
      kstar = 0;
      istl = max(lg.irs,istar - remapi - 2);
      istu = min(lg.ire,istar + remapi + 2);
      jstl = 0;
      jstu = 0;
      kstl = 0;
      kstu = 0;
      if (nd > 1){
	jstar = int((ypos - lg.zya[0])/lg.zdy[0]);
	jstl = max(lg.jrs,jstar - remapi - 2);
	jstu = min(lg.jre,jstar + remapi + 2);
      }
      if (nd == 3){
	kstar = int((zpos - lg.zza[0])/lg.zdz[0]); 
	kstl = max(lg.krs,kstar - remapi - 2);
	kstu = min(lg.kre,kstar + remapi + 2);
      }
      
      // for (int k = lg.krs; k <= lg.kre; ++k){
      for (int k = kstl; k <= kstu; ++k){
        zc = lg.zzc[k] - zpos;
        zc2 = zc*zc;
        //for (int j = lg.jrs; j <= lg.jre; ++j){
        for (int j = jstl; j <= jstu; ++j){
          yc = lg.zyc[j] - ypos;
          yc2 = yc*yc;
          //for (int i = lg.irs; i <= lg.ire; ++i){
          for (int i = istl; i <= istu; ++i){
	    xc = lg.zxc[i] - xpos;
	    xc2 = xc*xc;
	    r2 = xc2 + yc2 + zc2;
	    r = sqrt(r2);
	    xy = sqrt(xc2 + yc2);
	    sinphi = xy/r;
	    cosphi = zc/r;
	    costhta = xc/xy;
	    sinthta = yc/xy;
	    if (r <= remapRadius){
	      lg.P0[iqd][k][j][i] = mdot/(4.0*pi*r2*vinf);
	      lg.P0[iqe][k][j][i] = (lg.P0[iqd][k][j][i]/avgmassamb)*boltzman*tamb;
	      lg.P0[iqu0][k][j][i] = vinf*sinphi*costhta;
	      if (nd > 1)  lg.P0[iqu0+1][k][j][i] = vinf*sinphi*sinthta;
	      if (nd == 3) lg.P0[iqu0+2][k][j][i] = vinf*cosphi;
	      for (int m = 0; m < smax; ++m) lg.P0[iqal0+m][k][j][i] = 1.0 - float(n);
	    }
          }
        }
      }
    } // nstars
  }   // init T/F?
  else{
    if (problem == "CWB") xmaxst1 = xpos1 + dsep - rob; // maximum x-extent of wind
    else xmaxst1 = lg.zxa[imax+1];
    //cout << "xmaxst1 = " << xmaxst1 << "\n";
    // First map the wind of star 0
    for (int k = lg.krs; k <= lg.kre; ++k){
      zc = lg.zzc[k] - zpos1;
      zc2 = zc*zc;
      for (int j = lg.jrs; j <= lg.jre; ++j){
        yc = lg.zyc[j] - ypos1;
        yc2 = yc*yc;
        for (int i = lg.irs; i <= lg.ire; ++i){
	  if (lg.zxc[i] < xmaxst1){
	    xc = lg.zxc[i] - xpos1;
	    xc2 = xc*xc;
	    r2 = xc2 + yc2 + zc2;
	    r = sqrt(r2);
	    xy = sqrt(xc2 + yc2);
	    sinphi = xy/r;
	    cosphi = zc/r;
	    costhta = xc/xy;
	    sinthta = yc/xy;
	    lg.P0[iqd][k][j][i] = mdot1/(4.0*pi*r2*vinf1);
	    //cout << mdot1 << " " << r << " " << vinf1 << " " << lg.P0[iqd][k][j][i] << "\n";
	    //quit();
	    lg.P0[iqe][k][j][i] = (lg.P0[iqd][k][j][i]/avgmassamb)*boltzman*tamb;
	    lg.P0[iqu0][k][j][i] = vinf1*sinphi*costhta;
            if (nd > 1)  lg.P0[iqu0+1][k][j][i] = vinf1*sinphi*sinthta;
	    if (nd == 3) lg.P0[iqu0+2][k][j][i] = vinf1*cosphi;
	    for (int m = 0; m < smax; ++m) lg.P0[iqal0+m][k][j][i] = 1.0;
          }
        }
      }
    }
    // Now map the wind of star 1
    for (int k = lg.krs; k <= lg.kre; ++k){
      zc = lg.zzc[k] - zpos2;
      zc2 = zc*zc;
      for (int j = lg.jrs; j <= lg.jre; ++j){
        yc = lg.zyc[j] - ypos2;
        yc2 = yc*yc;
        for (int i = lg.irs; i <= lg.ire; ++i){
	  if (lg.zxc[i] >= xmaxst1){
	    xc = lg.zxc[i] - xpos2;
	    xc2 = xc*xc;
	    r2 = xc2 + yc2 + zc2;
	    r = sqrt(r2);
	    xy = sqrt(xc2 + yc2);
	    sinphi = xy/r;
	    cosphi = zc/r;
	    costhta = xc/xy;
	    sinthta = yc/xy;
	    lg.P0[iqd][k][j][i] = mdot2/(4.0*pi*r2*vinf2);
	    lg.P0[iqe][k][j][i] = (lg.P0[iqd][k][j][i]/avgmassamb)*boltzman*tamb;
	    lg.P0[iqu0][k][j][i] = vinf2*sinphi*costhta;
	    if (nd > 1)  lg.P0[iqu0+1][k][j][i] = vinf2*sinphi*sinthta;
	    if (nd == 3) lg.P0[iqu0+2][k][j][i] = vinf2*cosphi;
	    for (int m = 0; m < smax; ++m) lg.P0[iqal0+m][k][j][i] = 0.0;
          }
        }
      }
    }
    
  }
      
  return;
}




/*!  \brief Find the forward shock position.
 *
 *   \author Julian Pittard (Original version 05.06.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 10.08.11 (JMP)
 */
double findForwardShock(){
  int j = 0;
  int k = 0;
  double rshk = 0.0;
  for (int i = lg.ire; i >= lg.irs; --i){
    if (lg.P0[iqd][k][j][i] > 1.1*ramb){
      rshk = lg.zxc[i];
      break;
    }
  }
  return rshk;
}
