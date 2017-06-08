/*!  \file dtcon.cpp
 *   \brief Function for determining the safe timestep.
 * 
 *   The code can either call quit() or simply exit the compute_loop
 *   if the timestep becomes an INF, a NAN, or is zero or less.
 *
 *   \todo Double check the correct timestep for a Lagrangian code.
 *
 *   \author Julian Pittard (Original version 31.05.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 31.05.11 (JMP)
 */

#undef MAIN

//  Header files:
#include <cmath>              //needed to check for infinities and NaNs
#include <iomanip> 
#include <iostream>           //  The C++ iostream.
#include <math.h>
#include "constants.h"
#include "gas.h" 
#include "global.h"           // Global parameters
#include "grid.h"
#include "physconst.h"

using namespace std;

// Function declarations 
double ridt();
extern void quit();

void computeTimestep(){

  double dtx,dt3;

  if (problem == "SNR"){
    if      (ncycle < 10) cn = 0.000025;
    else if (ncycle < 20) cn = 0.00025;
    else if (ncycle < 30) cn = 0.0025;
    else if (ncycle < 40) cn = 0.025;
    else                  cn = 0.25;
  }
  else if (problem == "WBB"){
    if      (ncycle < 50) cn = 0.000025;
    else if (ncycle < 200) cn = 0.00025;
    else if (ncycle < 500) cn = 0.0025;
    else if (ncycle < 800) cn = 0.025;
    else                   cn = 0.25;
  }
  
  // global time constraint for given courant number

  dtx = cn / ridt();
                                        
  // limiting constraint on rate of increase of dt
  // dt3 = 1.1 * dt; // this causes problems if dt is adjusted for FITS file dumps etc.
  dt3 = dtx; 

  // use smallest required timestep                                        
  dt = min( dt3, dtx );

  // if timestep becomes too small, stop the program!

  if (t/dt > 1.0e10){
    if (procRank == 0){
      cout << "Timestep has become too small: dt = " << dt << "\n";
      cout << "                             time = " << t << "\n";
    }
    quit();
  }

  return;

}



double ridt(){

  int i,j,k;
  double xvel,yvel,zvel,widthy,widthz;
  double ri_dt = 0.0;
  j = 0;
  k = 0;
      
  cout << setprecision(4);
  //cout << "ridt start: procRank = " << procRank << "\n";

  if (nd == 1){
    for (int i = lg.irs; i <= lg.ire; ++i){
      svel = sqrt(gam*lg.P0[iqe][k][j][i]/lg.P0[iqd][k][j][i]);
      xvel = (fabs(lg.P0[iqu0][k][j][i])+svel)/lg.zdx[i];
      ri_dt = max(xvel,ri_dt);
    }
  }
  else if (nd == 2){
    for (int j = lg.jrs; j <= lg.jre; ++j){
      for (int i = lg.irs; i <= lg.ire; ++i){
        widthy = lg.zdy[j];
        if (lg.ngeomy > 2) widthy = widthy*lg.zxc[i];
        svel = sqrt(gam*lg.P0[iqe][k][j][i]/lg.P0[iqd][k][j][i]);
	xvel = (fabs(lg.P0[iqu0][k][j][i])+svel)/lg.zdx[i];
	yvel = (fabs(lg.P0[iqu0+1][k][j][i])+svel)/widthy;
	ri_dt = max(xvel,max(yvel,ri_dt));
      }
    }
  }
  else{ //nd == 3
    for (int k = lg.krs; k <= lg.kre; ++k){
      for (int j = lg.jrs; j <= lg.jre; ++j){
        for (int i = lg.irs; i <= lg.ire; ++i){
	  widthy = lg.zdy[j];
	  widthz = lg.zdz[k];
  	  if(lg.ngeomy > 2) widthy = widthy*lg.zxc[i];
	  if(lg.ngeomz > 2) widthz = widthz*lg.zxc[i];
	  if(lg.ngeomz == 5) widthz = widthz*sin(lg.zyc[j]);
	  svel = sqrt(gam*lg.P0[iqe][k][j][i]/lg.P0[iqd][k][j][i]);
	  xvel = (fabs(lg.P0[iqu0]  [k][j][i])+svel) / lg.zdx[i];
	  yvel = (fabs(lg.P0[iqu0+1][k][j][i])+svel) / widthy;
	  zvel = (fabs(lg.P0[iqu0+2][k][j][i])+svel) / widthz;
	  ri_dt = max(xvel,max(yvel,max(zvel,ri_dt)));
	} 
      }
    }
  }

  return ri_dt;

}
