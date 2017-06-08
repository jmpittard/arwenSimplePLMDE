/*!  \file update.cpp
 *   \brief Update the solution. 
 *
 *   See ~/work/stellarcluster/SSCwind_withIndividualSNe/ARWENv1.0
 *   for a more sophisticated version of this code, which deals
 *   with the grid being split across multiple processors by MPI. 
 *
 *   \author Julian Pittard (Original version 09.08.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 09.08.11 (JMP)
 */

#undef MAIN

//  Header files:
#include <cmath>
#include "constants.h"
#include "funcs.h"
#include "gas.h"
#include "global.h"
#include "grid.h"
#include "gridlet.h"
#include "physconst.h"

// Function declarations
void calculateAreas(int i, int j, int k);
void calculateFluxes(int i, int j, int k);
void calculateP(int i, int j, int k);
void calculatePstar(int i, int j, int k);
void calculateSourceTerms(int i, int j, int k);
void calculateVolumes(int i, int j, int k);
void constructGridlet(int i, int j, int k, double P[][kmaxd][jmaxd][imaxd]);
void consToPrim(double c[], double p[]);
void fillGhostCells(double P[][kmaxd][jmaxd][imaxd]);
void operatorsplit();
double pminFixUp(double avgmass, double tamb, double r, double p);
void primToCons(double p[], double c[]);
void riemann(int direc, double pL[], double pR[], double flux[]);
void riemann_nl(double dl, double ul, double pl, double cl, double dr, double ur, double pr, double cr, double rs[3]);
void setupGridletIndices();
void unsplit();
void update();
double wave(double p, double p0);
extern void mapWind(bool init);
extern double pminFixUp(double avgmassamb, double tamb, double r, double p);
extern void quit();

using namespace std;

// Global arrays with file scope
double p[iqmax][kk][jj][ii]; // primitive variables on gridlet
double c[iqmax][kk][jj][ii]; // conserved
double s[iqmax][kk][jj][ii]; // source terms
double diff[nd][kk][jj][ii][iqmax]; // primitive differences in each direction
double flux[nd][kk][jj][ii][iqmax]; // conserved fluxes      ---------"-------
double pL[iqmax];
double pR[iqmax];
double vol[kk][jj][ii]; // cell volumes on gridlet
double area[3][kk][jj][ii]; // cell areas in each direction on gridlet


/*!  \brief Set the values of the gridlet indices, which are only
 *   visible to the functions in this file.
 *
 *   \author Julian Pittard (Original version 09.08.11)
 *   \version 0.1-development (Und�miel):
 *   \date Last modified: 09.08.11 (JMP)
 */
void setupGridletIndices(){
  irs = nghost;
  ire = nghost+ig;
  if (nd > 1){
    jrs = nghost;
    jre = nghost+jg;
  }
  else{
    jrs = 0;
    jre = 1;
  }
  if (nd == 3){
    krs = nghost;
    kre = nghost+kg;
  }
  else{
    krs = 0;
    kre = 1;
  }
  
  return;
}


/*!  \brief Perform operator split activities (e.g. cooling, non-hydro evolution of scalars etc.)
 *
 *   \author Julian Pittard (Original version 09.08.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 09.08.11 (JMP)
 */
void operatorsplit(){

  if (problem == "WBB"){
    mapWind(false); 
  }
  else if (problem == "CWB"){
    mapWind(false); 
  }
  
  return;
}

/*!  \brief Update the fluid variables by advancing the solution.
 *
 *   A directionally unsplit update is performed on specific regions 
 *   (gridlets) on the grid.
 * 
 *   \author Julian Pittard (Original version 02.06.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 09.08.11 (JMP)
 */
void update(){

  setupGridletIndices();
  
  // Alternate operator split activities to achieve 2nd order in time
  if (ncycle%2 == 0){
    unsplit();
    operatorsplit();
  }
  else{
    operatorsplit();
    unsplit();
  }
  
  t       += dt;
  timep1d += dt;
  timep2d += dt;
  timep3d += dt;

  return;
}


/*!  \brief Advance the hydro solution using a directionally unsplit scheme.
 *
 *   \author Julian Pittard (Original version 28.03.12)
 *   \version 1.0-stable (Evenstar):
 *   \date Last modified: 19.04.12 (JMP)
 */
void unsplit(){

  fillGhostCells(lg.P0);

  // Predictor step: calculates Pstar at t^(n+1)
  // Loop over gridlets
  for (int k = 0; k < nkg; ++k){
    for (int j = 0; j < njg; ++j){
      for (int i = 0; i < nig; ++i){
	constructGridlet(i,j,k,lg.P0); // the gridlet p is equal to the grid P0
	calculateVolumes(i,j,k);
	calculateAreas(i,j,k);
	calculateSourceTerms(i,j,k);
	calculateFluxes(i,j,k);
	calculatePstar(i,j,k);
      }
    }
  }
  
  fillGhostCells(lg.Pstar);

  // Corrector step: 
  // Loop over gridlets
  for (int k = 0; k < nkg; ++k){
    for (int j = 0; j < njg; ++j){
      for (int i = 0; i < nig; ++i){
	constructGridlet(i,j,k,lg.Pstar); // the gridlet p is equal to the grid Pstar
	calculateVolumes(i,j,k);
	calculateAreas(i,j,k);
	calculateSourceTerms(i,j,k);
	calculateFluxes(i,j,k);
	calculateP(i,j,k);
      }
    }
  }
  
  return;
}


/*!  \brief Fill the ghost cells on the grid with their necessary values.
 *
 *   This function acts on the Grid, not a gridlet.
 *   P is passed in. This is either P0 (during the predictor step) or Pstar
 *   (during the corrector step). The boundaries are filled in the order XYZ 
 *   in order to fill the corner cells properly.
 *
 *   Each edge of the Grid represents a physical boundary at the edge of the
 *   global grid, and the BCs are specified by an integer flag set by the user.
 *   
 *   The values of the integer flags (nleftx, etc.) are:
 *    - 1 = reflecting; 2 = outflow; 3 = fixed, 4 = periodic; 5 = outflow ONLY
 *   These are set in setup().
 *
 *   \author Julian Pittard (Original version 28.03.12)
 *   \version 1.0-stable (Evenstar):
 *   \date Last modified: 19.04.12 (JMP)
 */
void fillGhostCells(double P[][kmaxd][jmaxd][imaxd]){

  // Leftx ghost cells    
  if (lg.nleftx == 0){ // reflecting
    for (int k = lg.krs; k <= lg.kre; ++k){
      for (int j = lg.jrs; j <= lg.jre; ++j){
        for (int i = 0; i < nghost; ++i){
          for (int n = 0; n < iqmax; ++n){
	    P[n][k][j][i] = P[n][k][j][2*nghost-1-i];
	  }
	  P[iqu0][k][j][i] *= -1.0; // reverse x-velocity
	}
      }
    }
  }
  else if (lg.nleftx == 1){ // outflow
    for (int k = lg.krs; k <= lg.kre; ++k){
      for (int j = lg.jrs; j <= lg.jre; ++j){
        for (int i = 0; i < nghost; ++i){
          for (int n = 0; n < iqmax; ++n){
	    P[n][k][j][i] = P[n][k][j][lg.irs];
	  }
	}
      }
    }
  }
  else if (lg.nleftx == 2){ // fixed
    for (int k = lg.krs; k <= lg.kre; ++k){
      for (int j = lg.jrs; j <= lg.jre; ++j){
        for (int i = 0; i < nghost; ++i){
          for (int n = 0; n < iqmax; ++n){
	    P[n][k][j][i] = lg.Pbc[n][0][0];
	  }
	}
      }
    }
  }
  else if (lg.nleftx == 3){ // periodic
    for (int k = lg.krs; k <= lg.kre; ++k){
      for (int j = lg.jrs; j <= lg.jre; ++j){
        for (int i = 0; i < nghost; ++i){
          for (int n = 0; n < iqmax; ++n){
	    P[n][k][j][i] = P[n][k][j][lg.ire+1-nghost+i];
	  }
	}
      }
    }
  }
  else if (lg.nleftx == 4){ // outflow only
    for (int k = lg.krs; k <= lg.kre; ++k){
      for (int j = lg.jrs; j <= lg.jre; ++j){
        for (int i = 0; i < nghost; ++i){
          for (int n = 0; n < iqmax; ++n){
	    P[n][k][j][i] = P[n][k][j][lg.irs];
	  }
	  if (P[iqu0][k][j][i] > 0.0) P[iqu0][k][j][i] *= -1.0;
	}
      }
    }
  }

  // Rightx ghost cells
  if (lg.nrightx == 0){ // reflecting
    for (int k = lg.krs; k <= lg.kre; ++k){
      for (int j = lg.jrs; j <= lg.jre; ++j){
        for (int i = 0; i < nghost; ++i){
          for (int n = 0; n < iqmax; ++n){
	    P[n][k][j][lg.ire+1+i] = P[n][k][j][lg.ire-i];
	  }
	  P[iqu0][k][j][lg.ire+1+i] *= -1.0; // reverse x-velocity
	}
      }
    }
  }
  else if (lg.nrightx == 1){ // outflow
    for (int k = lg.krs; k <= lg.kre; ++k){
      for (int j = lg.jrs; j <= lg.jre; ++j){
        for (int i = 0; i < nghost; ++i){
          for (int n = 0; n < iqmax; ++n){
	    P[n][k][j][lg.ire+1+i] = P[n][k][j][lg.ire];
	  }
	}
      }
    }
  }
  else if (lg.nrightx == 2){ // fixed
    for (int k = lg.krs; k <= lg.kre; ++k){
      for (int j = lg.jrs; j <= lg.jre; ++j){
        for (int i = 0; i < nghost; ++i){
          for (int n = 0; n < iqmax; ++n){
	    P[n][k][j][lg.ire+1+i] = lg.Pbc[n][0][1];
	  }
	}
      }
    }
  }
  else if (lg.nrightx == 3){ // periodic
    for (int k = lg.krs; k <= lg.kre; ++k){
      for (int j = lg.jrs; j <= lg.jre; ++j){
        for (int i = 0; i < nghost; ++i){
          for (int n = 0; n < iqmax; ++n){
	    P[n][k][j][lg.ire+1+i] = P[n][k][j][lg.irs+i];
	  }
	}
      }
    }
  }
  else if (lg.nrightx == 4){ // outflow only
    for (int k = lg.krs; k <= lg.kre; ++k){
      for (int j = lg.jrs; j <= lg.jre; ++j){
        for (int i = 0; i < nghost; ++i){
          for (int n = 0; n < iqmax; ++n){
	    P[n][k][j][lg.ire+1+i] = P[n][k][j][lg.ire];
	  }
	  if (P[iqu0][k][j][lg.ire+1+i] > 0.0) P[iqu0][k][j][lg.ire+1+i] *= -1.0;
	}
      }
    }
  }
  
  if (nd > 1){
  // Lefty ghost cells
  if (lg.nlefty == 0){ // reflecting
    for (int k = lg.krs; k <= lg.kre; ++k){
      for (int i = 0; i < imaxd; ++i){
        for (int j = 0; j < nghost; ++j){
          for (int n = 0; n < iqmax; ++n){
	    P[n][k][j][i] = P[n][k][2*nghost-1-j][i];
	  }
	  P[iqu0+1][k][j][i] *= -1.0; // reverse y-velocity
	}
      }
    }
  }
  else if (lg.nlefty == 1){ // outflow
    for (int k = lg.krs; k <= lg.kre; ++k){
      for (int i = 0; i < imaxd; ++i){
        for (int j = 0; j < nghost; ++j){
          for (int n = 0; n < iqmax; ++n){
	    P[n][k][j][i] = P[n][k][lg.jrs][i];
	  }
	}
      }
    }
  }
  else if (lg.nlefty == 2){ // fixed
    for (int k = lg.krs; k <= lg.kre; ++k){
      for (int i = 0; i < imaxd; ++i){
        for (int j = 0; j < nghost; ++j){
          for (int n = 0; n < iqmax; ++n){
	    P[n][k][j][i] = lg.Pbc[n][1][0];
	  }
	}
      }
    }
  }
  else if (lg.nlefty == 3){ // periodic
    for (int k = lg.krs; k <= lg.kre; ++k){
      for (int i = 0; i < imaxd; ++i){
        for (int j = 0; j < nghost; ++j){
          for (int n = 0; n < iqmax; ++n){
	    P[n][k][j][i] = P[n][k][lg.jre+1-nghost+j][i];
	  }
	}
      }
    }
  }
  else if (lg.nlefty == 4){ // outflow only
    for (int k = lg.krs; k <= lg.kre; ++k){
      for (int i = 0; i < imaxd; ++i){
        for (int j = 0; j < nghost; ++j){
          for (int n = 0; n < iqmax; ++n){
	    P[n][k][j][i] = P[n][k][lg.jrs][i];
	  }
	  if (P[iqu0+1][k][j][i] > 0.0) P[iqu0+1][k][j][i] *= -1.0;
	}
      }
    }
  }
  // Righty ghost cells
  if (lg.nrighty == 0){ // reflecting
    for (int k = lg.krs; k <= lg.kre; ++k){
      for (int i = 0; i < imaxd; ++i){
        for (int j = 0; j < nghost; ++j){
          for (int n = 0; n < iqmax; ++n){
	    P[n][k][lg.jre+1+j][i] = P[n][k][lg.jre-j][i];
	  }
	  P[iqu0+1][k][lg.jre+1+j][i] *= -1.0; // reverse y-velocity
	}
      }
    }
  }
  else if (lg.nrighty == 1){ // outflow
    for (int k = lg.krs; k <= lg.kre; ++k){
      for (int i = 0; i < imaxd; ++i){
        for (int j = 0; j < nghost; ++j){
          for (int n = 0; n < iqmax; ++n){
	    P[n][k][lg.jre+1+j][i] = P[n][k][lg.jre][i];
	  }
	}
      }
    }
  }
  else if (lg.nrighty == 2){ // fixed
    for (int k = lg.krs; k <= lg.kre; ++k){
      for (int i = 0; i < imaxd; ++i){
        for (int j = 0; j < nghost; ++j){
          for (int n = 0; n < iqmax; ++n){
	    P[n][k][lg.jre+1+j][i] = lg.Pbc[n][1][1];
	  }
	}
      }
    }
  }
  else if (lg.nrighty == 3){ // periodic
    for (int k = lg.krs; k <= lg.kre; ++k){
      for (int i = 0; i < imaxd; ++i){
        for (int j = 0; j < nghost; ++j){
          for (int n = 0; n < iqmax; ++n){
	    P[n][k][lg.jre+1+j][i] = P[n][k][lg.jrs+j][i];
	  }
	}
      }
    }
  }
  else if (lg.nrighty == 4){ // outflow only
    for (int k = lg.krs; k <= lg.kre; ++k){
      for (int i = 0; i < imaxd; ++i){
        for (int j = 0; j < nghost; ++j){
          for (int n = 0; n < iqmax; ++n){
	    P[n][k][lg.jre+1+j][i] = P[n][k][lg.jre][i];
	  }
	  if (P[iqu0+1][k][lg.jre+1+j][i] > 0.0) P[iqu0+1][k][lg.jre+1+j][i] *= -1.0;
	}
      }
    }
  }  
  } // nd > 1

  if (nd == 3){
  // Leftz ghost cells
  if (lg.nleftz == 0){ // reflecting
    for (int j = 0; j < jmaxd; ++j){
      for (int i = 0; i < imaxd; ++i){
        for (int k = 0; k < nghost; ++k){
          for (int n = 0; n < iqmax; ++n){
	    P[n][k][j][i] = P[n][2*nghost-1-k][j][i];
	  }
	  P[iqu0+2][k][j][i] *= -1.0; // reverse z-velocity
	}
      }
    }
  }
  else if (lg.nleftz == 1){ // outflow
    for (int j = 0; j < jmaxd; ++j){
      for (int i = 0; i < imaxd; ++i){
        for (int k = 0; k < nghost; ++k){
          for (int n = 0; n < iqmax; ++n){
	    P[n][k][j][i] = P[n][lg.krs][j][i];
	  }
	}
      }
    }
  }
  else if (lg.nleftz == 2){ // fixed
    for (int j = 0; j < jmaxd; ++j){
      for (int i = 0; i < imaxd; ++i){
        for (int k = 0; k < nghost; ++k){
          for (int n = 0; n < iqmax; ++n){
	    P[n][k][j][i] = lg.Pbc[n][2][0];
	  }
	}
      }
    }
  }
  else if (lg.nleftz == 3){ // periodic
    for (int j = 0; j < jmaxd; ++j){
      for (int i = 0; i < imaxd; ++i){
        for (int k = 0; k < nghost; ++k){
          for (int n = 0; n < iqmax; ++n){
	    P[n][k][j][i] = P[n][lg.kre+1-nghost+k][j][i];
	  }
	}
      }
    }
  }
  else if (lg.nleftz == 4){ // outflow only
    for (int j = 0; j < jmaxd; ++j){
      for (int i = 0; i < imaxd; ++i){
        for (int k = 0; k < nghost; ++k){
          for (int n = 0; n < iqmax; ++n){
	    P[n][k][j][i] = P[n][lg.krs][j][i];
	  }
	  if (P[iqu0+2][k][j][i] > 0.0) P[iqu0+2][k][j][i] *= -1.0;
	}
      }
    }
  }
  // Rightz ghost cells
  if (lg.nrightz == 0){ // reflecting
    for (int j = 0; j < jmaxd; ++j){
      for (int i = 0; i < imaxd; ++i){
        for (int k = 0; k < nghost; ++k){
          for (int n = 0; n < iqmax; ++n){
	    P[n][lg.kre+1+k][j][i] = P[n][lg.kre-k][j][i];
	  }
	  P[iqu0+2][lg.kre+1+k][j][i] *= -1.0; // reverse z-velocity
	}
      }
    }
  }
  else if (lg.nrightz == 1){ // outflow
    for (int j = 0; j < jmaxd; ++j){
      for (int i = 0; i < imaxd; ++i){
        for (int k = 0; k < nghost; ++k){
          for (int n = 0; n < iqmax; ++n){
	    P[n][lg.kre+1+k][j][i] = P[n][lg.kre][j][i];
	  }
	}
      }
    }
  }
  else if (lg.nrightz == 2){ // fixed
    for (int j = 0; j < jmaxd; ++j){
      for (int i = 0; i < imaxd; ++i){
        for (int k = 0; k < nghost; ++k){
          for (int n = 0; n < iqmax; ++n){
	    P[n][lg.kre+1+k][j][i] = lg.Pbc[n][2][1];
	  }
	}
      }
    }
  }
  else if (lg.nrightz == 3){ // periodic
    for (int j = 0; j < jmaxd; ++j){
      for (int i = 0; i < imaxd; ++i){
        for (int k = 0; k < nghost; ++k){
          for (int n = 0; n < iqmax; ++n){
	    P[n][lg.kre+1+k][j][i] = P[n][lg.krs+k][j][i];
	  }
	}
      }
    }
  }
  else if (lg.nrightz == 4){ // outflow only
    for (int j = 0; j < jmaxd; ++j){
      for (int i = 0; i < imaxd; ++i){
        for (int k = 0; k < nghost; ++k){
          for (int n = 0; n < iqmax; ++n){
	    P[n][lg.kre+1+k][j][i] = P[n][lg.kre][j][i];
	  }
	  if (P[iqu0+2][lg.kre+1+k][j][i] > 0.0) P[iqu0+2][lg.kre+1+k][j][i] *= -1.0;
	}
      }
    }
  }  
  } // nd == 3
  return;
}


/*!  \brief Construct gridlet.
 *
 *   Pgrid is passed in. This is either P0 (during the predictor step) or Pstar
 *   (during the corrector step).
 *
 *   The gridlet indices igg, jgg and kgg are passed in.
 *
 *   \author Julian Pittard (Original version 08.03.13)
 *   \version 1.0-stable (Evenstar):
 *   \date Last modified: 05.09.13 (JMP)
 */
void constructGridlet(int igg, int jgg, int kgg, double Pgrid[][kmaxd][jmaxd][imaxd]){

  // Loop over all (real and boundary) gridlet cells, extracting
  // data from real and boundary cells on the local grid.
  // igg = gridlet index
  // ig = number of real cells in gridlet
  // ii = total number of cells in gridlet
  // il = index of cell on grid (which stores ghost and real cells)

  int kl,jl,il;
    
  for (int n = 0; n < iqmax; ++n){
    for (int k = 0; k < kk; ++k){
      kl = kgg*kg + k;
      for (int j = 0; j < jj; ++j){
	jl = jgg*jg + j;
	for (int i = 0; i < ii; ++i){
	  il = igg*ig + i;
	  p[n][k][j][i] = Pgrid[n][kl][jl][il];
	}
      }
    }
  }

  //for (int i = 0; i < ii; ++i){
  //  cout << "Gridlet: i = " << i << "; p[iqe] = " << p[iqe][0][2][i] << "\n";
  //}
  //quit();
  
  return;
}


/*!  \brief Calculate the volume of cells on the gridlet.
 *
 *   The gridlet indices igg, jgg and kgg are passed in.
 *
 *   \author Julian Pittard (Original version 08.03.13)
 *   \version 1.0-stable (Evenstar):
 *   \date Last modified: 05.09.13 (JMP)
 */
void calculateVolumes(int igg, int jgg, int kgg){

  int gi,gj,gk;
  double r1,r2;

  // Default values
  double dx = 1.0;
  double dy = 1.0;
  double dz = 1.0;
    
  if (lg.geom == "XYZ"){
    for (int k = 0; k < kk; ++k){
      gk = kgg*kg + k; // index of cell on grid
      dz = lg.zdz[gk];
      for (int j = 0; j < jj; ++j){
	gj = jgg*jg + j; // index of cell on grid
	dy = lg.zdy[gj];
	for (int i = 0; i < ii; ++i){
	  gi = igg*ig + i; // index of cell on grid
	  dx = lg.zdx[gi];
	  vol[k][j][i] = dx*dy*dz;
	}
      }
    }
  }
  else if (lg.geom == "ZRP"){
    for (int k = 0; k < kk; ++k){
      gk = kgg*kg + k; // index of cell on grid
      dz = lg.zdz[gk];
      if (nd == 2) dz = pi;
      for (int j = 0; j < jj; ++j){
	gj = jgg*jg + j;   // index of cell on grid
	r1 = lg.zya[gj];   // inner radius of cell
	r2 = lg.zya[gj+1]; // outer radius
	for (int i = 0; i < ii; ++i){
	  gi = igg*ig + i; // index of cell on grid
	  dx = lg.zdx[gi];
	  vol[k][j][i] = mabs(r2*r2 - r1*r1)*dx*dz;
          if (nd == 3){
	    cout << "calculateVolumes: check this!\n";
            quit();
	  }
	}
      }
    }
  }
  else if (lg.geom == "RTP"){
    for (int k = 0; k < kk; ++k){
      gk = kgg*kg + k; // index of cell on grid
      dz = lg.zdz[gk];
      for (int j = 0; j < jj; ++j){
	gj = jgg*jg + j; // index of cell on grid
	dy = lg.zdy[gj];
	for (int i = 0; i < ii; ++i){
	  gi = igg*ig + i; // index of cell on grid
	  r1 = lg.zxa[gi];    // inner radius of cell
	  r2 = lg.zxa[gi+1];  // outer radius
	  vol[k][j][i] = 4.0*pi*mabs(r2*r2*r2 - r1*r1*r1)*dy*dz/3.0; // dy and dz should be fractions of pi if 2D/3D, and =1.0 if 1D
          if (nd > 1){
	    cout << "calculateVolumes: check this. Only coded for 1D at the mo!\n";
            quit();
	  }
	}
      }
    }
  }

  //for (int k = 0; k < kk; ++k){
  //  for (int j = 0; j < jj; ++j){
  //    for (int i = 0; i < ii; ++i){
  //	cout << i << "\t" << j << "\t" << k << "\t" << vol[k][j][i] << "\n";
  //    }
  //  }
  //}
  //quit();
  return;
}


/*!  \brief Calculate the area of the faces of the cells on the gridlet.
 *
 *   The gridlet indices igg, jgg and kgg are passed in.
 *
 *   \author Julian Pittard (Original version 08.03.13)
 *   \version 1.0-stable (Evenstar):
 *   \date Last modified: 05.09.13 (JMP)
 */
void calculateAreas(int igg, int jgg, int kgg){

  int gi,gj,gk;
  double r1,r2;

  // Default values
  double dx = 1.0;
  double dy = 1.0;
  double dz = 1.0;
    
  if (lg.geom == "XYZ"){
    for (int k = 0; k < kk; ++k){
      gk = kgg*kg + k; // index of cell on grid
      dz = lg.zdz[gk];
      for (int j = 0; j < jj; ++j){
	gj = jgg*jg + j; // index of cell on grid
	dy = lg.zdy[gj];
	for (int i = 0; i < ii; ++i){
	  gi = igg*ig + i; // index of cell on grid
	  dx = lg.zdx[gi];
	  area[0][k][j][i] = dy*dz; // face normal along x
	  area[1][k][j][i] = dx*dz; // face normal along y
	  area[2][k][j][i] = dx*dy; // face normal along z
	}
      }
    }
  }
  else if (lg.geom == "ZRP"){
    for (int k = 0; k < kk; ++k){
      gk = kgg*kg + k; // index of cell on grid
      dz = lg.zdz[gk];
      for (int j = 0; j < jj; ++j){
	gj = jgg*jg + j;   // index of cell on grid
	r1 = lg.zya[gj];   // inner radius of cell
	r2 = lg.zya[gj+1]; // outer radius
	for (int i = 0; i < ii; ++i){
	  gi = igg*ig + i; // index of cell on grid
	  dx = lg.zdx[gi];
	  area[0][k][j][i] = mabs((r2*r2 - r1*r1)*pi);
	  area[1][k][j][i] = 2.0*pi*mabs(r1)*dx; // face normal along "r" (dx is actually dz)
           if (nd == 3){
	    cout << "calculateAreas: check this!\n";
            quit();
	  }
	}
      }
    }
  }
  else if (lg.geom == "RTP"){
    for (int k = 0; k < kk; ++k){
      gk = kgg*kg + k; // index of cell on grid
      dz = lg.zdz[gk];
      for (int j = 0; j < jj; ++j){
	gj = jgg*jg + j; // index of cell on grid
	dy = lg.zdy[gj];
	for (int i = 0; i < ii; ++i){
	  gi = igg*ig + i; // index of cell on grid
	  r1 = lg.zxa[gi];    // inner radius of cell
	  r2 = lg.zxa[gi+1];  // outer radius
	  area[0][k][j][i] = 4.0*pi*r1*r1; // dy and dz should be fractions of pi if 2D/3D, and = 1.0 if 1D
          if (nd > 1){
	    cout << "calculateAreas: check this too!\n";
            quit();
	  }
	}
      }
    }
  }
  
  //for (int k = 0; k < kk; ++k){
  //  for (int j = 0; j < jj; ++j){
  //    for (int i = 0; i < ii; ++i){
  //	cout << i << "\t" << j << "\t" << k << "\t" << area[0][k][j][i] << "\t" << area[1][k][j][i] << "\t" << area[2][k][j][i] << "\n";
  //    }
  //  }
  //}
  //quit();
  return;
}


/*!  \brief Calculate the source terms.
 *
 *   The gridlet indices igg, jgg and kgg are passed in.
 *
 *   Note that the spatial variation of geometric source term's
 *   across the cell has to be taken account of. See the Appendix in 
 *   Falle (1991, MN, 250, 581) and p85-86 of Hydro Book III for
 *   further details.
 *
 *   \author Julian Pittard (Original version 08.03.13)
 *   \version 1.0-stable (Evenstar):
 *   \date Last modified: 05.09.13 (JMP)
 */
void calculateSourceTerms(int igg, int jgg, int kgg){

  int gi,gj,gk,dicell;
  double rg,rc,dr,fac,fac2,pre,difflp,diffrp,dpdr;
    
  // Calculate the geometric source term for each real cell on the gridlet.
  // Initialize...
  for (int n = 0; n < iqmax; ++n){
    for (int k = krs; k <= kre; ++k){
      for (int j = jrs; j <= jre; ++j){
	for (int i = irs; i <= ire; ++i){
	  s[n][k][j][i] = 0.0;
	}
      }
    }
  }

  if (lg.geom == "ZRP"){
    for (int k = krs; k <= kre; ++k){
      for (int j = jrs; j <= jre; ++j){
	gj = jgg*jg + j; // j index on the grid
	rc = lg.zyg[gj]; // Cell COM r-coordinate
	dr = lg.zdy[gj];
	dicell = gj-nghost+1; // j-index such that first real cell has index of 1 (to match Eq. A14 in the mg method paper) - see p85 of Hydro Book III
	fac = 2.0/(2.0*dicell - 1.0);
	fac2 = dr*(float(dicell) - 0.5) - rc;
	for (int i = irs; i <= ire; ++i){
	  pre = p[iqe][k][j][i]; // initial pressure
	  difflp = (p[iqe][k][j][i] - p[iqe][k][j-1][i])/(rc - lg.zyg[gj-1]);
	  diffrp = (p[iqe][k][j+1][i] - p[iqe][k][j][i])/(lg.zyg[gj+1] - rc);
	  dpdr = av(diffrp, difflp);
	  s[iqu0+1][k][j][i] += (fac/dr)*(pre + fac2*dpdr); // r-mtm geometric source term in ZRP geometry
	}
      }
    }
  }
  else if (lg.geom == "RTP"){
    for (int k = krs; k <= kre; ++k){
      for (int j = jrs; j <= jre; ++j){
	for (int i = irs; i <= ire; ++i){
	  pre = p[iqe][k][j][i]; // initial pressure
	  gi = igg*ig + i; // i index on the grid
	  rg = lg.zxg[gi]; // Cell COM r-coordinate
	  rc = lg.zxc[gi]; // Cell centre r-coordinate
	  dr = lg.zdx[gi];
	  dicell = gi-nghost+1; // i-index such that first real cell has index of 1 (to match Eq. A14 in the mg method paper) - see p85/86 of Hydro Book III
	  fac = 3.0*(2.0*dicell - 1.0)/(3.0*dicell*dicell - 3.0*dicell + 1.0);
	  fac2 = dr*(dicell - 0.5) - rg;
	  difflp = (p[iqe][k][j][i] - p[iqe][k][j][i-1])/(rg - lg.zxg[gi-1]);
	  diffrp = (p[iqe][k][j][i+1] - p[iqe][k][j][i])/(lg.zxg[gi+1] - rg);
	  dpdr = av(diffrp, difflp);
	  s[iqu0][k][j][i] += (fac/dr)*(pre + fac2*dpdr); // r-mtm geometric source term in RTP geometry 
	}
      }
    }
  }
  
  //for (int k = krs; k <= kre; ++k){
  //  for (int j = jrs; j <= jre; ++j){
  //    for (int i = irs; i <= ire; ++i){
  //	cout << i << "\t" << j << "\t" << k << "\t" << s[iqd][k][j][i] << "\t" << s[iqu0+1][k][j][i] << "\t" << s[iqe][k][j][i] << "\n";
  //    }
  //  }
  //}
  //quit();

  return;
}


/*!  \brief Calculate the fluxes.
 *
 *   The gridlet indices igg, jgg and kgg are passed in.
 *
 *   The following procedure occurs:
 *        i)   Calculate the average spatial derivative (monotonicity preserved
 *             using Eq. A18 from Falle 1991, MN, 250, 581) in each direction
 *             (store in the file scope array diff)
 *        ii)  Using diff, calculate the L and R states at each face (pL and pR)
 *        iii) Calculate the fluxes from the Riemann solver (store in the
 *             file scope array flux)
 *
 *   The value stored in flux[0][kg][jg][ig] is the flux across the lower-face 
 *   of cell[kg][jg][ig] (i.e. between cells [ig-1] and [i]).
 *
 *   It makes sense to perform the flux calculations one direction at a time, since
 *   this avoids nasty conditional statements within nested loops. Because of this,
 *   it also made sense to calculate diffl and diffr directly in the loops (a function
 *   call would have reintroduced a conditional statement within the loops).
 *   [UPDATE: In fact I don't think this is the case - probably I could do something
 *   more like Sam's mg code where diffl/r are actually gradients rather than just
 *   differences].
 *
 *   To avoid mass building up on symmetry axes (e.g the r=0 axis in ZRP geometry),
 *   two things are required:
 *     1) The initial conditions/winds must be set up using the cell centre-of-mass coordinates.
 *     2) 2nd order spatial derivatives must be calculated using the cell COM coords.
 *   Note that problems still occur if only one of the above conditions is met. 
 * 
 *   NOTE: See p187 of Hydro Book II.
 *
 *   \author Julian Pittard (Original version 08.03.13)
 *   \version 1.0-stable (Evenstar):
 *   \date Last modified: 05.09.13 (JMP)
 */
void calculateFluxes(int igg, int jgg, int kgg){

  int id;
  int igmin,igmax,igmin2,jgmin,kgmin,jgmin2,kgmin2,jgmax,kgmax;
  int gi,gj,gk;
  double diffl,diffr;
  double xint,yint,zint,fracL,fracR;
    
  // Initialize diff and flux
  for (int id = 0; id < nd; ++id){
    for (int k = 0; k < kk; ++k){
      for (int j = 0; j < jj; ++j){
	for (int i = 0; i < ii; ++i){
          for (int n = 0; n < iqmax; ++n){
	    diff[id][k][j][i][n] = 0.0;
	    flux[id][k][j][i][n] = 0.0;
	  }
	}
      }
    }
  }

  // Use a piecewise linear interpolation to determine pL and pR, using the same
  // code as in mg. Do not make any distinction between the 1st and 2nd parts of the timestep.

  // Set default loop ranges 
  igmin = irs-1;
  igmax = ire+1;
  igmin2 = irs;
  jgmin = 0;
  kgmin = 0;
  jgmin2 = 0;
  kgmin2 = 0;
  jgmax = 1;
  kgmax = 1;
  if (nd > 1){
    jgmin = jrs-1;
    jgmin2 = jrs;
    jgmax = jre+1;
  }
  if (nd == 3){
    kgmin = krs-1;
    kgmin2 = krs;
    kgmax = kre+1;
  }

  // For each direction, calculate the interpolated L/R state, then
  // pass to the Riemann solver to determine the fluxes.
  // The value in flux[id][k][j][i] is the net flux through the lower faces of cell[k,j,i] 

  // 1) Direction "0"

  id = 0;
  
  // Calculate the average gradient between cells. diff needs to be calculated
  // in the ghost zones adjacent to the real zones.
  for (int k = kgmin; k < kgmax; ++k){
    for (int j = jgmin; j < jgmax; ++j){
      for (int i = igmin; i < igmax; ++i){
	// Calculate the left and right differences, and their monotonicity preserving average
        for (int n = 0; n < iqmax; ++n){
	  diffl = p[n][k][j][i] - p[n][k][j][i-1];
	  diffr = p[n][k][j][i+1] - p[n][k][j][i];
	  diff[id][k][j][i][n] = av(diffl,diffr);
	}
      }
    }
  }

  // Calculate the L and R states at each real face.
  // pL is the interpolated value of p at the left side
  //            of the lower face of cell k,j,i
  // pL is the interpolated value of p at the right side
  //            of the lower face of cell k,j,i
  // Fluxes need to be calculated for faces irs to ire+1

  for (int k = kgmin2; k < kgmax; ++k){
    for (int j = jgmin2; j < jgmax; ++j){
      for (int i = igmin2; i < igmax; ++i){
	gi = igg*ig + i;    // i index on the grid
	xint = lg.zxa[gi];  // x-coord of left face of cell i on gridlet
	fracL = (xint - lg.zxg[gi-1])/lg.zdx[gi-1];
	fracR = (lg.zxg[gi] - xint)/lg.zdx[gi];
        for (int n = 0; n < iqmax; ++n){
	  pL[n] = p[n][k][j][i-1] + fracL*diff[id][k][j][i-1][n];
	  pR[n] = p[n][k][j][i] - fracR*diff[id][k][j][i][n];
	}
	riemann(id,pL,pR,flux[id][k][j][i]);
      }
    }
  }

  // 2) Direction "1"

  if (nd > 1){
    id = 1;
    // Calculate the average gradient between cells
    for (int k = kgmin; k < kgmax; ++k){
      for (int j = jgmin; j < jgmax; ++j){
        for (int i = igmin; i < igmax; ++i){
	  // Calculate the left and right differences, and their monotonicity preserving average
          for (int n = 0; n < iqmax; ++n){
	    diffl = p[n][k][j][i] - p[n][k][j-1][i];
	    diffr = p[n][k][j+1][i] - p[n][k][j][i];
	    diff[id][k][j][i][n] = av(diffl,diffr);
	  }
        }
      }
    }

    // Calculate the L and R states at each real face.
    // pL is the interpolated value of p at the left side
    //            of the lower face of cell k,j,i
    // pL is the interpolated value of p at the right side
    //            of the lower face of cell k,j,i
    // Fluxes need to be calculated for faces irs to ire+1

    for (int k = kgmin2; k < kgmax; ++k){
      for (int j = jgmin2; j < jgmax; ++j){
	gj = jgg*jg + j;    // j index on the grid
	yint = lg.zya[gj];  // y-coord of left face of cell j on gridlet
	fracL = (yint - lg.zyg[gj-1])/lg.zdy[gj-1];
	fracR = (lg.zyg[gj] - yint)/lg.zdy[gj];
        for (int i = igmin2; i < igmax; ++i){
          for (int n = 0; n < iqmax; ++n){
	    pL[n] = p[n][k][j-1][i] + fracL*diff[id][k][j-1][i][n];
	    pR[n] = p[n][k][j][i] - fracR*diff[id][k][j][i][n];
	  }
  	  riemann(id,pL,pR,flux[id][k][j][i]);
	}
      }
    }
  }

  
  // 3) Direction "2"

  if (nd == 3){
    id = 2;
    // Calculate the average gradient between cells
    for (int k = kgmin; k < kgmax; ++k){
      for (int j = jgmin; j < jgmax; ++j){
        for (int i = igmin; i < igmax; ++i){
	  // Calculate the left and right differences, and their monotonicity preserving average
          for (int n = 0; n < iqmax; ++n){
	    diffl = p[n][k][j][i] - p[n][k-1][j][i];
	    diffr = p[n][k+1][j][i] - p[n][k][j][i];
	    diff[id][k][j][i][n] = av(diffl,diffr);
	  }
        }
      }
    }

    // Calculate the L and R states at each real face.
    // pL is the interpolated value of p at the left side
    //            of the lower face of cell k,j,i
    // pL is the interpolated value of p at the right side
    //            of the lower face of cell k,j,i
    // Fluxes need to be calculated for faces irs to ire+1

    for (int k = kgmin2; k < kgmax; ++k){
      gk = kgg*kg + k;    // k index on the grid
      zint = lg.zza[gk];  // z-coord of left face of cell k on gridlet
      fracL = (zint - lg.zzg[gk-1])/lg.zdz[gk-1];
      fracR = (lg.zzg[gk] - zint)/lg.zdz[gk];
      for (int j = jgmin2; j < jgmax; ++j){
        for (int i = igmin2; i < igmax; ++i){
          for (int n = 0; n < iqmax; ++n){
	    pL[n] = p[n][k-1][j][i] + fracL*diff[id][k-1][j][i][n];
	    pR[n] = p[n][k][j][i] - fracR*diff[id][k][j][i][n];
	  }
  	  riemann(id,pL,pR,flux[id][k][j][i]);
	}
      }
    }
  }
  
  //for (int k = kgmin2; k < kgmax; ++k){
    //for (int j = jgmin2; j < jgmax; ++j){
  //  for (int j = 25; j < 26; ++j){
  //    for (int i = igmin2; i < igmax; ++i){
  //	cout << k << "\t" << j << "\t" << i << "\t" << flux[0][k][j][i][iqd] << "\t" << flux[0][k][j][i][iqu0] << "\t" << flux[0][k][j][i][iqu0+1] << "\t" << flux[0][k][j][i][iqe] << "\t" << flux[1][k][j][i][iqd] << "\t" << flux[1][k][j][i][iqu0] << "\t" << flux[1][k][j][i][iqu0+1] << "\t" << flux[1][k][j][i][iqe] << "\n";
  //    }
  //  }
  //}
  //quit();

  return;
}


/*!  \brief Calculate the intermediate state Pstar.
 *
 *   The gridlet indices igg, jgg and kgg are passed in.
 *   Update the real cells. flux has correct values for i=irs to i=ire, inclusive.
 *   
 *   \author Julian Pittard (Original version 08.03.13)
 *   \version 1.0-stable (Evenstar):
 *   \date Last modified: 05.09.13 (JMP)
 */
void calculatePstar(int igg, int jgg, int kgg){

  int gi,gj,gk;
  int igmin,jgmin,kgmin,igmax,jgmax,kgmax;
  double dtvol,area1,area2;
  double corig[iqmax];
  double cstar[iqmax];
  double porig[iqmax];
  double pstar[iqmax];

  // Pstar is calculated at t^(n+1) = t^(n) + dt
  double delt = dt;
        
  // Set default loop ranges. irs -> ire-1 inclusive will update only the real cells on the gridlet
    
  igmin = irs;
  igmax = ire;
  jgmin = 0;
  kgmin = 0;
  jgmax = 1;
  kgmax = 1;
  if (nd > 1){
    jgmin = jrs;
    jgmax = jre;
  }
  if (nd == 3){
    kgmin = krs;
    kgmax = kre;
  }

  for (int k = kgmin; k < kgmax; ++k){
    gk = kgg*kg + k; // k index on the grid
    for (int j = jgmin; j < jgmax; ++j){
      gj = jgg*jg + j; // j index on the grid
      for (int i = igmin; i < igmax; ++i){
	gi = igg*ig + i; // i index on the grid
	dtvol = delt/vol[k][j][i];

	for (int n = 0; n < iqmax; ++n){
	  porig[n] = p[n][k][j][i];
	}
	primToCons(porig,corig);  // calculate the conservative variables
	for (int n = 0; n < iqmax; ++n){
	  cstar[n] = corig[n];
	}
                
	// Multiply fluxes by cell areas
	// For cell[k,j,i], the net fluxes through the left and right faces
	// in direction "0" are (+flux[id][k][j][i] and -flux[id][k][j][i+1]). 
	// The left and right faces have an area of area[id][k][j][i] and area[id][k][j][i+1].

	// Direction "0"
	area1 = area[0][k][j][i];   // area of left face 
	area2 = area[0][k][j][i+1]; // area of right face 

	for (int n = 0; n < iqmax; ++n){
	  cstar[n] += (flux[0][k][j][i][n]*area1 - flux[0][k][j][i+1][n]*area2)*dtvol;
	}

	if (nd > 1){
	  // Direction "1"
	  area1 = area[1][k][j][i];   // area of left face 
	  area2 = area[1][k][j+1][i]; // area of right face 

  	  for (int n = 0; n < iqmax; ++n){
	    cstar[n] += (flux[1][k][j][i][n]*area1 - flux[1][k][j+1][i][n]*area2)*dtvol;
	  }
	}

	if (nd == 3){
	  // Direction "2"
	  area1 = area[2][k][j][i];   // area of left face 
	  area2 = area[2][k+1][j][i]; // area of right face 

  	  for (int n = 0; n < iqmax; ++n){
	    cstar[n] += (flux[2][k][j][i][n]*area1 - flux[2][k+1][j][i][n]*area2)*dtvol;
	  }
	}
	
	// Add source terms
	for (int n = 0; n < iqmax; ++n){
	  cstar[n] += delt*s[n][k][j][i];
	}

	// Now calculate Pstar on the grid
	consToPrim(cstar,pstar);
  	for (int n = 0; n < iqmax; ++n){
	  lg.Pstar[n][gk][gj][gi] = pstar[n];
	}
      }
    }
  }

  //for (int k = kgmin; k < kgmax; ++k){
  //  gk = kgg*kg + k; // k index on the grid
    //for (int j = jgmin; j < jgmax; ++j){
  //  for (int j = 25; j < 26; ++j){
  //    gj = jgg*jg + j; // j index on the grid
  //    for (int i = igmin; i < igmax; ++i){
  //	gi = igg*ig + i; // i index on the grid
  //	cout << i << "\t" << j << "\t" << k << "\t" << lg.Pstar[iqe][gk][gj][gi] << "\n";
  //  }
  //}
  //}
  //quit();
  return;
}


/*!  \brief Calculate the updated state P at time t + dt.
 *
 *   The gridlet indices igg, jgg and kgg are passed in.
 *   Update the real cells. flux has correct values for i=irs to i=ire, inclusive.
 *   
 *   \author Julian Pittard (Original version 08.03.13)
 *   \version 1.0-stable (Evenstar):
 *   \date Last modified: 05.09.13 (JMP)
 */
void calculateP(int igg, int jgg, int kgg){

  int gi,gj,gk;
  int igmin,jgmin,kgmin,igmax,jgmax,kgmax;
  double dtvol,area1,area2;
  double c[iqmax];
  double cstar[iqmax];
  double cnew[iqmax];
  double gridP[iqmax];
  double pstar[iqmax];
  double pnew[iqmax];

  // P is calculated by using 0.5*dt
  double delt = 0.5*dt;
        
  // Set default loop ranges. irs -> ire-1 inclusive will update only the real cells on the gridlet
    
  igmin = irs;
  igmax = ire;
  jgmin = 0;
  kgmin = 0;
  jgmax = 1;
  kgmax = 1;
  if (nd > 1){
    jgmin = jrs;
    jgmax = jre;
  }
  if (nd == 3){
    kgmin = krs;
    kgmax = kre;
  }

  for (int k = kgmin; k < kgmax; ++k){
    gk = kgg*kg + k; // k index on the grid
    for (int j = jgmin; j < jgmax; ++j){
      gj = jgg*jg + j; // j index on the grid
      for (int i = igmin; i < igmax; ++i){
	gi = igg*ig + i; // i index on the grid
	dtvol = delt/vol[k][j][i];

	for (int n = 0; n < iqmax; ++n){
	  pstar[n] = p[n][k][j][i];
	}
	primToCons(pstar,cstar);  // conservative variables from Pstar at t^(n+1)
	for (int n = 0; n < iqmax; ++n){
          gridP[n] = lg.P0[n][gk][gj][gi];
	}
	primToCons(gridP,c);      // conservative variables from P0 at t^(n)

	for (int n = 0; n < iqmax; ++n){
	  cnew[n] = 0.5*(cstar[n] + c[n]);
	}
                
	// Multiply fluxes by cell areas
	// For cell[k,j,i], the net fluxes through the left and right faces
	// in direction "0" are (+flux[id][k][j][i] and -flux[id][k][j][i+1]). 
	// The left and right faces have an area of area[id][k][j][i] and area[id][k][j][i+1].

	// Direction "0"
	area1 = area[0][k][j][i];   // area of left face 
	area2 = area[0][k][j][i+1]; // area of right face 

	for (int n = 0; n < iqmax; ++n){
	  cnew[n] += (flux[0][k][j][i][n]*area1 - flux[0][k][j][i+1][n]*area2)*dtvol;
	}

	if (nd > 1){
	  // Direction "1"
	  area1 = area[1][k][j][i];   // area of left face 
	  area2 = area[1][k][j+1][i]; // area of right face 

  	  for (int n = 0; n < iqmax; ++n){
	    cnew[n] += (flux[1][k][j][i][n]*area1 - flux[1][k][j+1][i][n]*area2)*dtvol;
	  }
	}

	if (nd == 3){
	  // Direction "2"
	  area1 = area[2][k][j][i];   // area of left face 
	  area2 = area[2][k+1][j][i]; // area of right face 

  	  for (int n = 0; n < iqmax; ++n){
	    cnew[n] += (flux[2][k][j][i][n]*area1 - flux[2][k+1][j][i][n]*area2)*dtvol;
	  }
	}
	
	// Add source terms
	for (int n = 0; n < iqmax; ++n){
	  cnew[n] += delt*s[n][k][j][i];
	}

	// Now calculate the new P on the grid, at t^(n+1)
	consToPrim(cnew,pnew);
  	for (int n = 0; n < iqmax; ++n){
	  lg.P0[n][gk][gj][gi] = pnew[n];
	}
      }
    }
  }

  //for (int k = kgmin; k < kgmax; ++k){
  //  gk = kgg*kg + k; // k index on the grid
    //for (int j = jgmin; j < jgmax; ++j){
  //  for (int j = 25; j < 26; ++j){
  //    gj = jgg*jg + j; // j index on the grid
  //    for (int i = igmin; i < igmax; ++i){
  //	gi = igg*ig + i; // i index on the grid
  //	cout << i << "\t" << j << "\t" << k << "\t" << lg.P0[iqe][gk][gj][gi] << "\n";
  //    }
  //  }
  //}
  //quit();
  return;
}


/*!  \brief Convert an array of primitives to conserved values.
 *
 *   \author Julian Pittard (Original version 08.03.13)
 *   \version 1.0-stable (Evenstar):
 *   \date Last modified: 05.09.13 (JMP)
 */
void primToCons(double p[iqmax], double c[iqmax]){

  c[iqd] = p[iqd];
  c[iqe] = 0.0;
  for (int n = 0; n < nd; ++n){
    c[iqu0+n] = c[iqd]*p[iqu0+n];
    c[iqe] += c[iqu0+n]*p[iqu0+n];
  }
  c[iqe] *= 0.5;
  c[iqe] += p[iqe]/gamm;

  for (int n = 0; n < smax; ++n){
    c[iqal0+n] = p[iqd]*p[iqal0+n];
  }
  
  return;
}

/*!  \brief Convert an array of conservatives to primitive values.
 *
 *   \author Julian Pittard (Original version 08.03.13)
 *   \version 1.0-stable (Evenstar):
 *   \date Last modified: 05.09.13 (JMP)
 */
void consToPrim(double c[iqmax], double p[iqmax]){

  double u,pre;
  double ke = 0.0;
  double rho = c[iqd];
    
  p[iqd] = rho;
  for (int n = 0; n < nd; ++n){
    u = c[iqu0+n]/p[iqd];
    p[iqu0+n] = u;
    ke += rho*u*u;
  }
  ke *= 0.5;
  pre = gamm*(c[iqe] - ke);
  p[iqe] = pre;

  for (int n = 0; n < smax; ++n){
    p[iqal0+n] = c[iqal0+n]/rho;
  }

  // Use additional criteria to stop p going too small. This is problem specific.
  p[iqe] = pminFixUp(avgmassamb,tamb,p[iqd],p[iqe]);
  return;
}


/*!  \brief Given left and right states (pL,pR), calculate the time-averaged flux.
 *
 *   Uses Sam's method in mg_g (i.e. in the hflux routine).
 *
 *   Note that direc=0 for the 1st axis (i.e. "X" in "XYZ" geometry), direc=1 for the 2nd axis, etc.
 *   This means that direc is different from the value passed to hflux in mg.
 *
 *   \todo Write a fixup routine more like those in mg_g.
 *
 *   \author Julian Pittard (Original version 18.03.13)
 *   \version 1.0-stable (Evenstar):
 *   \date Last modified: 21.03.13 (JMP)
 */
void riemann(int direc, double pL[iqmax], double pR[iqmax], double flux[iqmax]){

  double linear = 0.1;
  bool nlrieman = false;
  bool cwrieman = false;
  double rs[3];
  double vel[3];
  double cl,cr,ul,dl,pl,ur,dr,pr,cel,cer,p,u,d;
  double mflux,ud,ev,cvis,vis,tv;
    
  cl = 1.0; // defaults
  cr = 1.0;
    
  // Zero fluxes (safety)
  for (int n = 0; n < iqmax; ++n)
    flux[n] = 0.0;

  // Check that the values in pL and pR are sensible using problem-specific checks.
  pL[iqd] = max(pL[iqd],dfloor);
  pR[iqd] = max(pR[iqd],dfloor);
  pL[iqe] = pminFixUp(avgmassamb,tamb,pL[iqd],pL[iqe]);
  pR[iqe] = pminFixUp(avgmassamb,tamb,pR[iqd],pR[iqe]);

  int iqu = iqu0 + direc;

  if (cwrieman){
    // Use the Collea & Woodward (1984) 2-shock riemann solver
    //riemann_cw84(direc,pL,pR,d,u,p,vel)
    cout << "Not converted yet!\n";
    quit();
  }
  else{    
    // Use the riemann solver in mg
    // Left state
    ul = pL[iqu]; // normal velocity
    dl = pL[iqd];
    pl = pL[iqe];

    // Right state
    ur = pR[iqu]; // normal velocity
    dr = pR[iqd];
    pr = pR[iqe];

    // Calculate sound speeds
    cel = sqrt(gam*pl/dl);
    cl = dl*cel;
    cer = sqrt(gam*pr/dr);
    cr = dr*cer;

    // Linear Riemann solver
    p = (cr*pl + cl*pr - cr*cl*(ur - ul))/(cr + cl);
    u = (cr*ur + cl*ul - (pr - pl))/(cr + cl);
    if ((abs(p - pl) > linear*pl) || (abs(p - pr) > linear*pr)){
      // Use non-linear Riemann solver
      nlrieman = true;
      rs[0] = 1.0; // dummy argument (d)
      rs[1] = u;
      rs[2] = p;
      riemann_nl(dl, ul, pl, cl, dr, ur, pr, cr, rs);
      d = rs[0];
      u = rs[1];
      p = rs[2];
    }
     
    if (u >= 0.0){
      if (ul - cel > 0.0){
	d = dl;
        u = ul;
	p = pl;
      }
      else{
	d = cl*dl/(cl + (u - ul)*dl);
      }
      for (int n = 0; n < nd; ++n){
        vel[n] = pL[iqu0+n];
      }
    }
    else{
      if (ur + cer <= 0.0){
	d = dr;
	u = ur;
	p = pr;
      }
      else{
	d = cr*dr/(cr - (u - ur)*dr);
      }
      for (int n = 0; n < nd; ++n){
        vel[n] = pR[iqu0+n];
      }
    }
    vel[direc] = u;
  }
  
  // Make sure that the time-averaged fluid variables at the cell interface
  // (d, u and p) have sensible values. The pressure in particular is problem-specific
  // and so we load the interface values into a primitive variable for checking
  // by the routine pminFixUp.
  d = fmax(d,dfloor);
  p = pminFixUp(avgmassamb,tamb,d,p);

  // Mass flux and energy flux
  mflux = d*u;
  flux[iqd] = mflux;
  flux[iqe] = p;
  // Add velocity bits into momentum fluxes and energy
  for (int n = 0; n < nd; ++n){
    ud = vel[n];
    flux[iqu0+n] = mflux*ud;
    flux[iqe] += 0.5*d*ud*ud;
  }
  // Calculate internal energy
  ev = p/gamm;
  // Add this into energy flux, add pressure into normal momentum
  flux[iqe] = u*(flux[iqe] + ev);
  flux[iqu] += p;
  // Scalars
  if (u > 0.0){
    for (int n = 0; n < smax; ++n){
      flux[iqal0+n] = mflux*pL[iqal0+n];
    }
  }
  else{
    for (int n = 0; n < smax; ++n){
      flux[iqal0+n] = mflux*pR[iqal0+n];
    }
  }

  // Artificial dissipation
  cvis = fmax(cl,cr);
  vis = avisp*cvis;
  for (int n = 0; n < nd; ++n){
    tv = vis*(pL[iqu0+n] - pR[iqu0+n]);
    flux[iqu0+n] += tv;
    flux[iqe] += tv*vel[n];
  }
  flux[iqe] += avise*cvis*(pL[iqe]/pL[iqd] - pR[iqe]/pR[iqd]);
 
  return;

}


/*!  \brief Non-linear Riemann solver.
 *
 *   Calculates the resolved state which has dl,ul,pl on the left and dr,ur,pr
 *   on the right.
 *
 *   \author Julian Pittard (Original version 18.03.13)
 *   \version 1.0-stable (Evenstar):
 *   \date Last modified: 21.03.13 (JMP)
 */
void riemann_nl(double dl, double ul, double pl, double cl, double dr, double ur,
  double pr, double cr, double rs[3]){

  int nit = 0;
  double wl, wr, c, ci, p0;
  const double tol=1.0e-4;
  bool doit = true;
  double d = rs[0];
  double u = rs[1];
  double p = rs[2];
  
  if ((p < pl) && (p < pr)){
    // Non-iterative solution for two rarefactions
    ci = cr/(dr*(double)pow((double)pr, (double)g4));
    p = (0.5*g1*(ul - ur) + cr/dr + cl/dl)/(ci + cl/
      (dl*(double)pow((double)pl, (double)g4)));
    if (p < pfloor)
      p = pfloor;
    else
      p = (double)pow((double)p, (double)g5);
  }
  else{
    // Iterative solution for shock/rarefaction
    p0 = p;
    wl = -cl*wave(p,pl);
    wr = cr*wave(p,pr);
    // Iteration
    while (doit){
      p = (wr*pl - wl*pr + wr*wl*(ur - ul))/(wr - wl);
      if ((mabs(p - p0) > tol*p) && (p > 0.0)){
	p0 = p;
	wl = -cl*wave(p,pl);
	wr = cr*wave(p,pr);
	nit++;
	if (nit > 50){
	  cout << "convergence failure in riemann solver" << '\n';
	  // Use linear guess
	  p = (cr*pl + cl*pr - cr*cl*(ur - ul))/(cr + cl);
	  p = max(p, pfloor);
	  wl = -cl;
	  wr = cr;
	  doit = false;
	}
	else if (nit > 20)
	  cout << "slow convergence in riemann solver nit " <<
	    nit << '\n';
      }
      else
	doit = false;
    }
    // Calculate resolved density and velocity
    u = (wr*ur - wl*ul - pr + pl)/(wr -wl);
    // Check for position of contact
    if (u > 0.0){
      // Contact is on right of interface
      d = wl*dl/(wl - (u -ul)*dl);
      // Check velocity of left waves
      if (p < pl){
	// Left wave is rarefaction
	c = sqrt(gam*p/d);
	if ((u - c) > 0.0){
	  // Tail is on right of interface
	  cl /= dl;
	  if ((ul - cl) < 0.0){
	    // Head is on left of interface -- spans interface 
	    u = g8*(ul +g7*cl);
	    p = pl*pow((double)(u/cl), (double)g5);
	    d = gam*p/(u*u);
	  }
	  else{
	    // Rarefaction does not span interface
	    d = dl;
	    u = ul;
	    p = pl;
	  }
	}
      }
      else{
	// Left wave is shock
	if ((wl/dl+ul) > 0.0){
	  // On right of interface
	  d = dl;
	  u = ul;
	  p = pl;
	}
      }
    }
    else{
      //     Contact is on left of interface 
      d = wr*dr/(wr - (u - ur)*dr);
      // Check velocity of right waves
      if (p < pr){
	// Right wave is rarefaction 
	c = sqrt(gam*p/d);
	if ((u + c) < 0.0){
	  // Tail is on left of interface
	  cr /= dr;
	  if ((ur + cr) > 0.0){
	    // Head is on right of interface -- spans interface
	    u = g8*(ur - g7*cr);
	    p = pr*pow((double)(-u/cr), (double)g5);
	    d = gam*p/(u*u);
	  }
	  else{
	    // Rarefaction does not span interface
	    d = dr;
	    u = ur;
	    p = pr;
	  }
	}
      }
      else{
	// Right wave is shock 
	if ((wr/dr+ur) < 0.0){
	  // On left of interface 
	  d = dr;
	  u = ur;
	  p = pr;
	}
      }
    }
  }

  rs[0] = d;
  rs[1] = u;
  rs[2] = p;
  
  return;

}


/*!  \brief Calculate the nonlinear wave speed. 
 *
 *   This function calculates the wave speed for a wave connecting
 *   states with pressures p0 and p ahead of and behind respectively.
 *   From Sam Falle's wave function.
 *
 *   \author Julian Pittard (Original version 21.03.13)
 *   \version 1.0-stable (Evenstar):
 *   \date Last modified: 22.03.13 (JMP)
 */
double wave(double p, double p0){

  double x, w;
  x = p/p0;
  if (mabs(x - 1.0) < 1.0e-03)
    // Use linear expression
    w = 1.0 + 0.5*g3*(x - 1.0);
  else{
    // Use non-linear expression
    if (x >= 1.0)
      // Shock
      w = sqrt(1.0 + g3*(x - 1.0));
    else
      // Rarefaction
      w = g4*(1.0 - x)/(1.0 - pow((double)x, (double)g4));
  }

  return w;

}


/*!  \brief Calculate the minimum pressure, given the density. 
 *
 *   This is problem specific.
 *
 *   \author Julian Pittard (Original version 21.03.13)
 *   \version 1.0-stable (Evenstar):
 *   \date Last modified: 22.03.13 (JMP)
 */
double pminFixUp(double avgmassamb, double tamb, double r, double p){
  double const1 = boltzman*tamb/avgmassamb;
  double pmin;
    
  if (isnan(p) || isinf(p)) p = pfloor;
  if (isnan(r) || isinf(r)) r = dfloor;
        
  if ((problem == "WBB") || (problem == "SNR") || (problem == "CWB")){
    pmin = max(r*const1,pfloor);
    p = max(p,pmin);
  }
  else{
    p = max(p,pfloor);
  }
  return p;
}
