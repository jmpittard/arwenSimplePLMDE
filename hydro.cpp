/*!  \file hydro.cpp
 *   \brief Perform the hydro update.
 *
 *   sweepbc.cpp, volume.cpp, ppm.cpp and remap.cpp moved out of sweepx/y/z.cpp
 *   routines. Many variables moved into hydro.h.
 *
 *   \author Julian Pittard (Original version 01.08.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 08.08.11 (JMP)
 */


#undef MAIN

//  Header files:
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include "constants.h"
#include "gas.h"
#include "global.h"     // Global parameters
#include "grid.h"        
#include "hydro.h"
#include "physconst.h"
#include "sweep.h"

// Function declarations
void evolve(int n1, int n2, int ngeom);
void fict_forces(double xf[], double uf[], double vf[], double wf[], double fict[]);
void flatten();
void pminFixUp(int nmin, int nmax);
void riemann_cw84(int lmin, int lmax);
void sweepbc( int nleft, int nright);
void volume ( int ngeom);
extern void ppm(int n1, int n2, int nleft, int nright, int ngeom, double para[][maxsweep] );
extern void quit(); 
extern void remap( int nleft, int nright, int ngeom);


using namespace std;


/*!  \brief Function to perform the hydro update
 *
 * Set boundary conditions, compute volume elements, flatten (smooth)
 * the flow near shocks, calculate the input states to the Riemann solver,
 * obtain the zone face fluxes, and evolve the flow. Remap to a 
 * fixed or evolving grid if necessary
 *
 */
void hydroUpdate(int n1, int n2, int nleft, int nright, int ngeom, int ntot, double para[][maxsweep] ){

  // Set up boundary conditions
  sweepbc( nleft, nright);

  // Compute cell volumes
  volume ( ngeom);

  // Flatten (smooth) the flow near shocks 
  flatten();

  // Calculate the input states to the Riemann solver
  ppm( n1, n2, nleft, nright, ngeom, para );

  // Call the Riemann solver to obtain the zone face averages. For the
  // lagrangian solver we only require umid and pmid, while for the Eulerian
  // solver we require in addition rmid, vmid, wmid, and asmid

  riemann_cw84( nmin-1, nmax+2);      // OK for Lagrangian and Eulerian

  // Do lagrangian update using umid and pmid, or Eulerian update using
  // umid, rmid, pmid, vmid, wmid and asmid. 
  evolve( n1, n2, ngeom);

  // Remap to a fixed or expanding grid if necessary
  //if (lagrangian)
  remap( nleft, nright, ngeom);

  return;
}


/*!  \brief Adds appropriate boundary cells to 1D sweep arrays.
 *
 *   Boundary condition flags : nleft, nright
 *    \li = 0 : reflecting
 *    \li = 1 : outflow (zero gradients)
 *    \li = 2 : fixed inflow (eg, uinflo,pinflo,...)
 *    \li = 3 : periodic (eg, u(nmin-1) = u(nmax))
 *    \li = 4 : outflow ONLY (no inflow allowed)
 *
 *   \author Julian Pittard (Original version 16.05.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 30.09.11 (JMP)
 */
void sweepbc( int nleft, int nright){

  int n,m;
  if ( nleft == 0 ){
    for (n = 1; n <= nghost; n++){
      dx[nmin-n] = dx[nmin+n-1];
      xa[nmin-n] = xa[nmin-n+1] - dx[nmin-n];
      dx0[nmin-n]= dx0[nmin+n-1];
      xa0[nmin-n]= xa0[nmin-n+1] - dx0[nmin-n];
      r [nmin-n] = r [nmin+n-1];
      u [nmin-n] = -u[nmin+n-1];
      v [nmin-n] = v [nmin+n-1];
      w [nmin-n] = w [nmin+n-1];
      p [nmin-n] = p [nmin+n-1];
      e [nmin-n] = e [nmin+n-1];
      for (m = 0; m < smax; m++){ //advected scalars
        as[m][nmin-n] = as[m][nmin+n-1];
      }
    }
  }
  else if ( nleft == 1 ){
    for (n = 1; n <= nghost; n++){
      dx[nmin-n] = dx[nmin];
      xa[nmin-n] = xa[nmin-n+1] - dx[nmin-n];
      dx0[nmin-n]= dx0[nmin+n-1];
      xa0[nmin-n]= xa0[nmin-n+1] - dx0[nmin-n];
      r [nmin-n] = r [nmin];
      u [nmin-n] = u [nmin];
      v [nmin-n] = v [nmin];
      w [nmin-n] = w [nmin];
      p [nmin-n] = p [nmin];
      e [nmin-n] = e [nmin];
      for (m = 0; m < smax; m++){ //advected scalars
        as[m][nmin-n] = as[m][nmin];
      }
    }
  }
  else if ( nleft == 2 ){
    for (n = 1; n <= nghost; n++){
      dx[nmin-n] = dx[nmin];
      xa[nmin-n] = xa[nmin-n+1] - dx[nmin-n];
      dx0[nmin-n]= dx0[nmin+n-1];
      xa0[nmin-n]= xa0[nmin-n+1] - dx0[nmin-n];
      r [nmin-n] = dinflo;
      u [nmin-n] = uinflo;
      v [nmin-n] = vinflo;
      w [nmin-n] = winflo;
      p [nmin-n] = pinflo;
      e [nmin-n] = pinflo/(dinflo*gamm) + 0.5*
                      ( uinflo*uinflo + vinflo*vinflo + winflo*winflo );
      for (m = 0; m < smax; m++){ //advected scalars
        as[m][nmin-n] = asinflo[m];
      }
    }
  }
  else if ( nleft == 3 ){
    for (n = 1; n <= nghost; n++){
      dx[nmin-n] = dx[nmax+1-n];
      xa[nmin-n] = xa[nmin-n+1] - dx[nmin-n];
      dx0[nmin-n]= dx0[nmax+1-1];
      xa0[nmin-n]= xa0[nmin-n+1] - dx0[nmin-n];
      r [nmin-n] = r [nmax+1-n];
      u [nmin-n] = u [nmax+1-n];
      v [nmin-n] = v [nmax+1-n];
      w [nmin-n] = w [nmax+1-n];
      p [nmin-n] = p [nmax+1-n];
      e [nmin-n] = e [nmax+1-n];
      for (m = 0; m < smax; m++){ //advected scalars
        as[m][nmin-n] = as[m][nmax+1-n];
      }
    }
  }
  else if ( nleft == 4 ){  // outflow only
    for (n = 1; n <= nghost; n++){
      dx[nmin-n] = dx[nmin];
      xa[nmin-n] = xa[nmin-n+1] - dx[nmin-n];
      dx0[nmin-n]= dx0[nmin+n-1];
      xa0[nmin-n]= xa0[nmin-n+1] - dx0[nmin-n];
      r [nmin-n] = r [nmin];
      // We don't need to worry about grid splitting in the parallel code, because
      // the transpose means that we always sweep along a complete side of the grid
      if (u[nmin] > 0.0){
	u [nmin-n] = -u[nmin];  // reflect velocity to give zero flux
      }
      else u [nmin-n] = u [nmin];  // outflow
      v [nmin-n] = v [nmin];
      w [nmin-n] = w [nmin];
      p [nmin-n] = p [nmin];
      e [nmin-n] = e [nmin];
      for (m = 0; m < smax; m++){ //advected scalars
	as[m][nmin-n] = as[m][nmin];
      }
    }
  }

  if (nright == 0){
    for (n = 1; n <= nghost; n++){
      dx[nmax+n] = dx[nmax+1-n];
      xa[nmax+n] = xa[nmax+n-1] + dx[nmax+n-1];
      dx0[nmax+n]= dx0[nmax+1-n];
      xa0[nmax+n]= xa0[nmax+n-1] + dx0[nmax+n-1];
      r [nmax+n] = r [nmax+1-n];
      u [nmax+n] = -u[nmax+1-n];
      v [nmax+n] = v [nmax+1-n];
      w [nmax+n] = w [nmax+1-n];
      p [nmax+n] = p [nmax+1-n];
      e [nmax+n] = e [nmax+1-n];
      for (m = 0; m < smax; m++){ //advected scalars
        as[m][nmax+n] = as[m][nmax+1-n];
      }
    }
  }
  else if (nright == 1){
    for (n = 1; n <= nghost; n++){
      dx[nmax+n] = dx[nmax];
      xa[nmax+n] = xa[nmax+n-1] + dx[nmax+n-1];
      dx0[nmax+n]= dx0[nmax+1-n];
      xa0[nmax+n]= xa0[nmax+n-1] + dx0[nmax+n-1];
      r [nmax+n] = r [nmax];
      u [nmax+n] = u [nmax];
      v [nmax+n] = v [nmax];
      w [nmax+n] = w [nmax];
      p [nmax+n] = p [nmax];
      e [nmax+n] = e [nmax];
      for (m = 0; m < smax; m++){ //advected scalars
        as[m][nmax+n] = as[m][nmax];
      }
    }
  }
  else if (nright == 2){
    for (n = 1; n <= nghost; n++){
      dx[nmax+n] = dx [nmax];
      xa[nmax+n] = xa [nmax+n-1] + dx[nmax+n-1];
      dx0[nmax+n]= dx0[nmax+1-n];
      xa0[nmax+n]= xa0[nmax+n-1] + dx0[nmax+n-1];
      r [nmax+n] = dotflo;
      u [nmax+n] = uotflo;
      v [nmax+n] = votflo;
      w [nmax+n] = wotflo;
      p [nmax+n] = potflo;
      e [nmax+n] = potflo/(dotflo*gamm) + 0.5*
	           ( uotflo*uotflo + votflo*votflo + wotflo*wotflo );
      for (m = 0; m < smax; m++){ //advected scalars
        as[m][nmax+n] = asotflo[m];
      }
    }
  }
  else if (nright == 3){
    for (n = 1; n <= nghost; n++){
      dx[nmax+n] = dx[nmin+n-1];
      xa[nmax+n] = xa[nmax+n-1] + dx[nmax+n-1];
      dx0[nmax+n]= dx0[nmin+n-1];
      xa0[nmax+n]= xa0[nmax+n-1] + dx0[nmax+n-1];
      r [nmax+n] = r [nmin+n-1];
      u [nmax+n] = u [nmin+n-1];
      v [nmax+n] = v [nmin+n-1];
      w [nmax+n] = w [nmin+n-1];
      p [nmax+n] = p [nmin+n-1];
      e [nmax+n] = e [nmin+n-1];
      for (m = 0; m < smax; m++){ //advected scalars
        as[m][nmax+n] = as[m][nmin+n-1];
      }
    }
  }
  else if ( nright == 4 ){     
    for (n = 1; n <= nghost; n++){
      dx[nmax+n] = dx[nmax];
      xa[nmax+n] = xa[nmax+n-1] + dx[nmax+n-1];
      dx0[nmax+n]= dx0[nmax+1-n];
      xa0[nmax+n]= xa0[nmax+n-1] + dx0[nmax+n-1];
      r [nmax+n] = r [nmax];
      if (u[nmax] < 0.0){
	u [nmax+n] = -u[nmax]; // set opposite to ensure zero flux
      }
      else u [nmax+n] = u [nmax]; // outlfow
      v [nmax+n] = v [nmax];
      w [nmax+n] = w [nmax];
      p [nmax+n] = p [nmax];
      e [nmax+n] = e [nmax];
      for (m = 0; m < smax; m++){ //advected scalars
	as[m][nmax+n] = as[m][nmax];
      }
    }
  }

  return;
}




/*!  \brief Performs 1D sweeps in the x-direction.
 *
 *   For an expanding grid simply calculate a new xa0 and dx0 for the expansion.
 *   The parabolic coefficients in zparax should be calculated  
 *   from the original grid before the expansion, and updated
 *   each time step. However, if we are scaling a uniform grid, 
 *   the coeffs. remain constant, in which case there is no need
 *   to recompute...
 *
 *   \author Julian Pittard (Original version 02.06.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 01.08.11 (JMP)
 */
void sweepx(){
 
  //Locals
  int i,j,k,n,m,j1,j2;
  int ntot;
  double dr,u2,v2,w2;

  sweep = "x";

  // Add ghost zones to 1D sweep arrays
  nmin = nghost;                               // first true cell
  nmax = imax - 1 + nghost;                    // last true cell
  ntot = maxsweep;                             // MUST equal the value in zparax in grid.h

  // Expand the grid for each i. dr0 and dsep0 are the initial values and are 
  // not updated!
  /*
  if (expand){
    dr = dr0*pow((dsep/dsep0),expande);
    for (i=0; i<imax; i++){
      n = i + nghost;
      dx0[n] = dr;
      xa0[n] = float(i-ipivot)*dx0[n];
      // xa0(n) = xa0(n) + 0.4*sin(float(ncycle)/1.0)*dr; ! wiggle the grid
    }
  }
  */
  
  // Loop over each row...
  for (k = 0; k < kmax; k++){
    for (j = 0; j < jmax; j++){
      // Put state variables into 1D arrays, padding with ghost zones
      // Also extract the required component of the effective force(s)
      // (gravity plus radiative continuum) and the radiative line force(s). 
      for (i = 0; i < imax; i++){
        n = i + nghost;
	r  [n] = zro[k][j][i];
	p  [n] = zpr[k][j][i];
	u  [n] = zux[k][j][i];
	v  [n] = zuy[k][j][i];
	w  [n] = zuz[k][j][i];
        for (m = 0; m < smax; m++){
  	  as[m][n] = zas[k][j][i][m];  
        }
        //if (!expand){
          xa0[n] = zxa[i];
          dx0[n] = zdx[i];
	//}
        xa [n] = zxa[i];
        dx [n] = zdx[i];
        p  [n] = max(pfloor,p[n]);
	u2 = u[n]*u[n];
	v2 = v[n]*v[n];
	w2 = w[n]*w[n]; 
	e  [n] = p[n]/(r[n]*gamm)+0.5*(u2 + v2 + w2);
	/*
        if (force){
          if (rank == 0) cout << "Not coded for force in sweepx yet\n";
          quit();
          //act[n] = fx[k][j][i];
        }
        else
	*/
        act[n] = 0.0;
        //printf("%3d %15.6e %15.6e %15.6e\n",n,p[n],r[n],u[n]);
      }

      // Perform the hydro update
      hydroUpdate(j, k, nleftx, nrightx, ngeomx, ntot, zparax);

      // Put updated values back into 3D arrays, dropping ghost zones

      for (i = 0; i < imax; i++){
	n = i + nghost;
	zro[k][j][i] = r[n];
	zpr[k][j][i] = p[n];
	zux[k][j][i] = u[n];
	zuy[k][j][i] = v[n];
	zuz[k][j][i] = w[n];
        for (m = 0; m < smax; m++){
  	  zas[k][j][i][m] = as[m][n];
        }
      }

    }  // j loop
  }    // k loop

  /*
  // Store the new coords
  if (expand){
    for (i=0;i < imax; i++){
      n = i+nghost;
      zxa[i] = xa0[n];
      zdx[i] = dx0[n];
    }
    xmin = zxa[0];
    xmax = zxa[imax-1] + zdx[imax-1];
    xming = xmin;
    xmaxg = xmax;
    para_setup(); // compute new parabolic coeffs
  }
  */
  
  return;
}


/*!  \brief Performs 1D sweeps in the y-direction.
 *
 *   For an expanding grid simply calculate a new xa0 and dx0 for the expansion.
 *   The parabolic coefficients in zparay should be calculated  
 *   from the original grid before the expansion, and updated
 *   each time step. However, if we are scaling a uniform grid, 
 *   the coeffs. remain constant, in which case there is no need
 *   to recompute...
 *
 *   \author Julian Pittard (Original version 02.06.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 02.06.11 (JMP)
 */
void sweepy(){

  //Locals
  int i,j,k,n,m,i1,i2;
  int ntot;
  double dr;

  sweep = "y";
  
  // Add ghost zones to 1D sweep arrays
  nmin = nghost;                               // first true cell
  nmax = jmax - 1 + nghost;                    // last true cell
  ntot = maxsweep;                             // MUST equal the value in zparay in grid.h

  // Expand the grid for each j. dr0 and dsep0 are the initial values and are 
  // not updated!
  /*
  if (expand){
    dr = dr0*pow((dsep/dsep0),expande);
    for (j=0; j<jmax; j++){
      n = j + nghost;
      dx0[n] = dr;
      xa0[n] = float(j-jpivot)*dx0[n];
      // xa0[n] = xa0[n] + 0.4*sin(float(ncycle)/1.0)*dr; ! wiggle the grid
    }
  }
  */
  
  // Loop over each row...
  for (k = 0; k < kmax; k++){
    for (i = 0; i < imax; i++){
      radius = zxc[i];

      // Put state variables into 1D arrays, padding with ghost zones
      // Also extract the required component of the effective force(s)
      // (gravity plus radiative continuum) and the radiative line force(s). 
      for (j = 0; j < jmax; j++){
        n = j + nghost;
	r  [n] = zro[k][j][i];
	p  [n] = zpr[k][j][i];
	u  [n] = zuy[k][j][i];
	v  [n] = zuz[k][j][i];
	w  [n] = zux[k][j][i];
        for (m = 0; m < smax; m++){
  	  as[m][n] = zas[k][j][i][m];  
        }
        //if (!expand){
          xa0[n] = zya[j];
          dx0[n] = zdy[j];
	//}
        xa [n] = zya[j];
        dx [n] = zdy[j];
        p  [n] = max(pfloor,p[n]);
	e  [n] = p[n]/(r[n]*gamm)+0.5*(u[n]*u[n]+v[n]*v[n]+w[n]*w[n]);

	/*
        if (force){
          if (rank == 0) cout << "Not coded for force in sweepy yet\n";
          quit();
          //act[n] = fy[k][j][i];
        }
        else 
	*/
	act[n] = 0.0;
      }

      // Perform the hydro update

      hydroUpdate(i, k, nlefty, nrighty, ngeomy, ntot, zparay);

      // Put updated values back into 3D arrays, dropping ghost zones

      for (j = 0; j < jmax; j++){
	n = j + nghost;
	zro[k][j][i] = r[n];
	zpr[k][j][i] = p[n];
	zuy[k][j][i] = u[n];
	zuz[k][j][i] = v[n];
	zux[k][j][i] = w[n];
        for (m = 0; m < smax; m++){
  	  zas[k][j][i][m] = as[m][n];
        }
      }

    }  // i loop
  }    // k loop

  // Store the new coords
  /*
  if (expand){
    for (j=0;j < jmax; j++){
      n = j+nghost;
      zya[j] = xa0[n];
      zdy[j] = dx0[n];
    }
    ymin = zya[0];
    ymax = zya[jmax-1] + zdy[jmax-1];
    yming = ymin;
    ymaxg = ymax;
    para_setup();  // compute new parabolic coeffs
  }
  */
  
  return;
}


/*!  \brief Performs 1D sweeps in the z-direction.
 *
 *   For an expanding grid simply calculate a new xa0 and dx0 for the expansion.
 *   The parabolic coefficients in zparaz should be calculated  
 *   from the original grid before the expansion, and updated
 *   each time step. However, if we are scaling a uniform grid, 
 *   the coeffs. remain constant, in which case there is no need
 *   to recompute...
 *
 *   \author Julian Pittard (Original version 02.06.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 04.08.11 (JMP)
 */
void sweepz(){

  //Locals
  int i,j,k,n,m,k1,k2;
  int ntot;
  double dr;

  sweep = "z";
  
  // Add ghost zones to 1D sweep arrays
  nmin = nghost;                               // first true cell
  nmax = kmax - 1 + nghost;                    // last true cell
  ntot = maxsweep;                             // MUST equal the value in zparaz in grid.h

  // Expand the grid for each k. dr0 and dsep0 are the initial values and are 
  // not updated!
  /*
  if (expand){
    dr = dr0*pow((dsep/dsep0),expande);
    for (k=0; k<kmax; k++){
      n = k + nghost;
      dx0[n] = dr;
      xa0[n] = float(k-kpivot)*dx0[n];
      // xa0[n] = xa0[n] + 0.4*sin(float(ncycle)/1.0)*dr; ! wiggle the grid
    }
  }
  */
  
  // Loop over each row...

  for (j = 0; j < jmax; j++){
    for (i = 0; i < imax; i++){
      radius = zxc[i];
      theta  = zyc[j];
      stheta = sin(theta);
      radius = radius * stheta;

      // Put state variables into 1D arrays, padding with ghost zones
      // Also extract the required component of the effective force(s)
      // (gravity plus radiative continuum) and the radiative line force(s). 

      for (k = 0; k < kmax; k++){
        n = k + nghost;
	r  [n] = zro[k][j][i];
	p  [n] = zpr[k][j][i];
	u  [n] = zuz[k][j][i];
	v  [n] = zux[k][j][i];
	w  [n] = zuy[k][j][i];
        for (m = 0; m < smax; m++){
  	  as[m][n] = zas[k][j][i][m];  
        }
        //if (!expand){
          xa0[n] = zza[k];
          dx0[n] = zdz[k];
	//}
        xa [n] = zza[k];
        dx [n] = zdz[k];
        p  [n] = max(pfloor,p[n]);
	e  [n] = p[n]/(r[n]*gamm)+0.5*(u[n]*u[n]+v[n]*v[n]+w[n]*w[n]);

	/*
        if (force){
          if (rank == 0) cout << "Not coded for force in sweepz yet\n";
          quit();
          //act[n] = fz[k][j][i];
        }
        else */
	act[n] = 0.0;
      }

      // Perform the hydro update

      hydroUpdate(i, j, nleftz, nrightz, ngeomz, ntot, zparaz);

      // Put updated values back into 3D arrays, dropping ghost zones

      for (k = 0; k < kmax; k++){
	n = k + nghost;
	zro[k][j][i] = r[n];
	zpr[k][j][i] = p[n];
	zuz[k][j][i] = u[n];
	zux[k][j][i] = v[n];
	zuy[k][j][i] = w[n];
        for (m = 0; m < smax; m++){
  	  zas[k][j][i][m] = as[m][n];
        }
      }

    }  // i loop
  }    // j loop

  /*
  // Store the new coords
  if (expand){
    for (k=0;k < kmax; k++){
      n = k+nghost;
      zza[k] = xa0[n];
      zdz[k] = dx0[n];
    }
    zmin = zza[0];
    zmax = zza[kmax-1] + zdz[kmax-1];
    zming = zmin;
    zmaxg = zmax;
    para_setup(); // compute new parabolic coeffs
  }
  */
  
  return;
}


/*!  \brief Updates the fluid variables using the time averaged fluxes from
 *          the Riemann solver.
 *
 *   For the lagrangian solver, umid and pmid are used to update the 
 *   velocity, density, and total energy. For the Eulerian solver,
 *   rmid, vmid, wmid, and asmid are also used.
 *   Physical zones are from nmin to nmax.  Zone boundary numbers run from
 *   nmin to nmax+1.
 *
 *   smallp is used to prevent negative temperatures - the value is
 *   problem-specific. ewnd*rho(n) also sets a temperature floor, for
 *   instance preventing adiabatic cooling - this may not be desirable
 *   for some problems. 
 *
 *   Various error checking occurs: e.g. that there are no NaNs or INFs in
 *   the updated fluid variables, that advected scalars have sensible values,
 *   etc.
 *
 *   04.08.11 (JMP) - Replaced 1D Matrices with built-in arrays
 *
 *   \author Julian Pittard (Original version 05.06.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 04.08.11 (JMP)
 */
void evolve(int n1, int n2, int ngeom){

  // Locals
  int n,m;

  double pmin, scrch1, artd, arte, csp, hdt, dtdvol, dtdx, aux1, dtheta, smallp;

  double as_ex[smax];
  double avis[maxsweep]; // artifical viscosity
  double fict0[maxsweep], fict1[maxsweep],
    amid[maxsweep], rold[maxsweep], uold[maxsweep], rinv[maxsweep], 
      rhoflx[maxsweep], cflx[maxsweep], uflx[maxsweep], vflx[maxsweep], 
    wflx[maxsweep], eflx[maxsweep], dr[maxsweep], dm[maxsweep], dtbdm[maxsweep],
    xa1[maxsweep], xa2[maxsweep], xa3[maxsweep], dvol1[maxsweep], upmid[maxsweep];
  double asflx[smax][maxsweep];
  double as_exn1[smax], as_exn2[smax], as_flx[smax];

  //bool errorcheck = true;
  bool errorcheck = false;
  const double third = 1./3.;
  const double rndoff = 1.0e-09;

  //if (lagrangian){
    //grid position evolution
    for (n = nmin-2; n <= nmax+2; n++){
      dm   [n] = r[n] * dvol[n];
      dtbdm[n] = dt / dm[n];
      xa1  [n] = xa[n];
      dvol1[n] = dvol[n];                          // store original zone vol
      xa   [n] = xa[n] + dt * umid[n] / radius;
      upmid[n] = umid[n] * pmid[n];
    }

    xa1[nmin-3] = xa[nmin-3];
    xa1[nmax+3] = xa[nmax+3];
    for (n = nmin-2; n <= nmax+2; n++){
      xa2[n]   = xa1[n] + 0.5*dx[n];
      dx [n]   = xa[n+1] - xa[n];    
      xa3[n]   = xa [n] + 0.5*dx[n];
    }

    // Calculate forces using coordinates at t. Note that fictitious
    //   forces at t+dt depend on updated velocity, but we ignore this
    fict_forces(xa2,u,v,w,fict0); // at t    (xa1/xa2 is the original grid)
    fict_forces(xa3,u,v,w,fict1); // at t+dt (xa/xa3  is the evolved  grid)
    
    // Calculate dvolume and average area based on geometry of sweep
    // Maybe loop here can be from nmin-1 to nmax+2?
    if (ngeom == 0){
      for (n = nmin-2; n <= nmax+2; n++){
        dvol[n] = dx[n];
        amid[n] = 1.;
      }
    }
    else if (ngeom == 1){
      for (n = nmin-2; n <= nmax+2; n++){
        dvol[n] = dx[n]*(xa[n]+.5*dx[n]);
        amid[n] = .5*(xa[n] + xa1[n]); // average area of interface 
                                       // between t and t+dt
      }
    }
    else if (ngeom == 2){
      for (n = nmin-2; n <= nmax+2; n++){
        dvol[n] = dx[n]*(xa[n]*(xa[n]+dx[n])+dx[n]*dx[n]*third);
        amid[n] = (xa[n]-xa1[n])*(third*(xa[n]-xa1[n])+xa1[n])+xa1[n]*xa1[n];
      }
    }
    else if (ngeom == 3){
      for (n = nmin-2; n <= nmax+2; n++){
        dvol[n] = dx[n]*radius;
        amid[n] = 1.;
      }
    }
    else if (ngeom == 4){
      for (n = nmin-2; n <= nmax+2; n++){
        dvol[n] = (cos(xa[n])-cos(xa[n+1]))*radius;
        dtheta  = xa[n] - xa1[n];
        if (dtheta == 0.0) amid[n] = sin(xa[n]);
        else amid[n] = (cos(xa1[n])-cos(xa[n]))/dtheta;
      }
    }
    else if (ngeom == 5){
      for (n = nmin-2; n <= nmax+2; n++){
        dvol[n] = dx[n]*radius;
        amid[n] = 1.;
      }
    }
    else{
      if (procRank == 0) cout << "Geometry" << ngeom << " not implemented\n";
      quit();
    }

    // Evolve for the lagrangian code. For the density evolution, all we  
    // have to do is watch the change in the geometry.
    // velocity evolution due to pressure acceleration and forces, both 
    // fictious and external.
    // Total energy evolution (do not include fictious forces)
    // The fictious forces do not change the total energy, but rather 
    // redistribute the momentum between the different coordinate directions.
    // Note we rely on smallp to prevent negative temperatures. ewnd*rho[n]
    // will reheat material which may not be wanted if we have cooling on.
    // Because this is a lagrangian code, advected scalars are unchanged, and 
    // are not updated here.

    for (n = nmin-1; n <= nmax+1; n++){
      r[n] *= ( dvol1[n] / dvol[n] );
      r[n] = max(r[n],dfloor);
      uold [n] = u[n];
      u[n] -= dtbdm[n] * (pmid[n+1]-pmid[n])
                         * 0.5 * (amid[n+1]+amid[n])
	   + 0.5*dt*(act[n]+fict0[n]+act[n]+fict1[n]);
      double eold = e[n];
      e[n] -= dtbdm[n]*
                 (amid[n+1]*upmid[n+1] - amid[n]*upmid[n])
	   + .5*dt*(uold[n]*act[n] + u[n]*act[n]);
      p[n]  = r[n]*gamm*(e[n]-.5*(u[n]*u[n] + v[n]*v[n] + w[n]*w[n]));

      //smallp = max((r[n]/avgmassamb)*boltzman*tamb,pfloor);
      //p[n] = max(p[n],smallp);
      //cout << "evolve:" << n << " " << r[n] << " " << u[n] << " " << p[n] << " " << xa[n] << " " << dx[n] << " " << dvol[n] << " " << dvol1[n] << " " << e[n] << " " << eold << "\n";
    }

    pminFixUp(nmin,nmax);
    //quit();
    /*
  }
  else{

// Evolve for the Eulerian code
  } //Eulerian
    */
    
  return;
}


/*!  \brief Check the minimum pressure value. This is problem dependent.
 *
 *   \author Julian Pittard (Original version 05.06.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 10.08.11 (JMP)
 */
void pminFixUp(int nmin, int nmax){

  double pmin;
  double const1 = boltzman*tamb/avgmassamb;
  
  if (problem == "WBB" || problem == "SNR"){
    for (int n = nmin-1; n <= nmax+1; ++n){
      pmin = max(r[n]*const1,pfloor);
      p[n] = max(p[n],pmin);
    }
  }
  else{
    for (int n = nmin-1; n <= nmax+1; ++n){
      p[n] = max(p[n],pfloor);
    }

  }  
  return;
}


/*!  \brief Solve the Riemann shock tube problem for the left and right 
 *    input states
 *
 *   Uses the Colella & Woodard (1984) iterative, approximate 2-shock solver.
 *   The Newton iteration procedure described in van Leer (1979).
 *
 *   Because of the way the code is constructed, some of the ghost cell values 
 *   will be NaNs, but this is OK as long as this information is far enough 
 *   removed to not propagate into the real grid cells.
 *
 *   Note that the Euleian solver needs the information in w_lft and w_rgh.
 *  
 *   \author Julian Pittard (Original version 05.06.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 10.08.11 (JMP)
 */
void riemann_cw84(int lmin, int lmax){

  //Locals

  int m,n,l;
  double gamfac1, gamfac2;
  double rs,ps,us,vs,ws,w_s,c_s,sig,c_mid,ls,lmid,zeta;
  double tol,small_p;
  double w_lft[maxsweep], w_rgh[maxsweep], c_lft[maxsweep], 
    c_rgh[maxsweep], z_lft[maxsweep], z_rgh[maxsweep];
  double umidl, umidr, pmold, p_mid;
  double r_r, p_r, u_r, r_l, p_l, u_l;
  double w_r, c_r, z_r, w_l, c_l, z_l;
  double v_r, v_l;
  double ass[smax];

  //double gamma = gam;
  gamfac1 = 0.5*(gam+1.0)/gam;
  gamfac2 = gam + 1.0;
  tol     = 1.0e-3;
  small_p = 1.0e-20;

//
// Input variables are:
//
//   lmin   = zone number of first physical zone
//   lmax   = zone number of first ghost zone on right (lmax=nmax+1)
//   gam    = equation of state gamma
//   r_r     = density state on the right side of the boundary
//   r_l     = density state on the left side of the boundary
//   p_r     = pressure state on the right side of the boundary
//   p_l     = pressure state on the left side of the boundary
//   u_r     = velocity state on the right side of the boundary
//   u_l     = velocity state on the left side of the boundary
//  and for Eulerian part only,
//   v_r     = transverse velocity state on the right side of the boundary
//   v_l     = transverse velocity state on the left side of the boundary
//   w_r     = transverse velocity state on the right side of the boundary
//   w_l     = transverse velocity state on the left side of the boundary
//   asr(m) = (m)th advected scalar on the right side of the boundary
//   asl(m) = (m)th advected scalar on the left side of the boundary
//
// Output variables are:
// 
//   umid   = time averaged velocity at the zone interface
//   pmid   = time averaged pressure at the zone interface
//   rmid   = time averaged density at the zone interface
//   vmid   = time averaged transverse velocity at the zone interface
//   wmid   = time averaged transverse velocity at the zone interface
//   asmid  = time averaged (m)th advected scalar at the zone interface
//
//


  for (l = lmin; l <= lmax; l++){
    p_l = plft[l];
    u_l = ulft[l];
    r_l = rlft[l];
    p_r = prgh[l];
    u_r = urgh[l];
    r_r = rrgh[l];

// Obtain first guess for Pmid by assuming W_lft, W_rgh = C_lft, C_rgh
//
    c_l = sqrt(gam*p_l*r_l);
    c_r = sqrt(gam*p_r*r_r);
    // See Eq. 72 in Fryxell et al. (2000, ApJSS, 131, 273) for the next two lines
    p_mid = p_r - p_l - c_r*(u_r-u_l);     
    p_mid = p_l + p_mid * c_l/(c_l+c_r);      
    p_mid = max(small_p,p_mid);

    v_l = 1.0/r_l; // specific volume (saves divisions in the iterator)
    v_r = 1.0/r_r;
//
// Iterate up to 8 times using Newton's method to converge on correct Pmid
//     -use previous guess for pmid to get wavespeeds: w_l, w_r
//     -find the slope in the u-P plane for each state: z_l, z_r
//     -use the wavespeeds and pmid to guess umid on each side: umidl, umidr
//     -project tangents from (pmid,umidl) and (pmid,umidr) to get new pmid
//     -make sure pmid does not fall below floor value for pressure
//

    for (n = 1; n <= 8; n++){
      pmold = p_mid;
      w_l  = 1.0 + gamfac1*(p_mid - p_l) / p_l;   //Eq 77
      w_r  = 1.0 + gamfac1*(p_mid - p_r) / p_r;
      w_l  = c_l * sqrt(w_l);      
      w_r  = c_r * sqrt(w_r);   
      //w_lft[l] = w_l; // Save - these are needed in the Eulerian part 
      //w_rgh[l] = w_r;
      z_l  = 4.0 * w_l * w_l * v_l;      
      z_r  = 4.0 * w_r * w_r * v_r;   
      z_l  = -z_l * w_l/(z_l - gamfac2*(p_mid - p_l));      
      z_r  =  z_r * w_r/(z_r - gamfac2*(p_mid - p_r));     
      umidl = u_l - (p_mid - p_l) / w_l;      
      umidr = u_r + (p_mid - p_r) / w_r;   
      p_mid  = p_mid + (umidr - umidl)
                * (z_l  * z_r ) / (z_r  - z_l );   
      p_mid  = max(small_p,p_mid);
      if (fabs(p_mid-pmold)/p_mid < tol) break;
    }
    // 
    // Calculate umid by averaging umidl, umidr based on new pmid
    //
    umidl = u_l - (p_mid - p_l) / w_l;      
    umidr = u_r + (p_mid - p_r) / w_r;   
    umid[l]  = 0.5*(umidl + umidr);
    pmid[l] = p_mid;
    //cout << "riemann: " << l << " " << r_l << " " << p_l << " " << u_l << " " << r_r << " " << p_r << " " << u_r << " " << umid[l] << " " << pmid[l] << "\n";
  }

// At this point the Lagrangian solver returns umid and pmid. For the
// Eulerian solver, some extra computation needs to be done (see Flash 
// paper p 293)


  return;
}


/*!  \brief Look for signs of strong shocks and contact discontinuities.
 *
 *   Computes values for flat, jump and shck. Global numerical diffusion can
 *   be added by increasing the minimum value of the flattening parameter from
 *   0.0 to say 0.2 (see Sutherland 2010, Ap&SS, 327, 173, section 6.1.7).
 *   However, tests have shown that this can seriously mess up flow aligned with
 *   the grid, while it does not do as good a job at stopping the carbuncle
 *   instability in 2D axisymmetric sims of CWBs as the code in the Fortran version
 *   of vh1xe.
 *
 *   \var flat (double) returns a value between 0.0 (smooth flow) 
 *        and 1.0 (strong shock) using a simplified method of C&W: eqn. 
 *        A.1 and A.2.
 *
 *   \var jump (double) returns a measure of the pressure jump,
 *    which is needed in order to determine if a CD exists and therefore
 *    whether to apply steepening in the Eulerian solver.
 *
 *   10.08.11 (JMP) - reduced loops from +/- 4 to +/- 2, and sweep[n] to +/- 3.
 *
 *   \author Julian Pittard (Original version 05.06.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 10.08.11 (JMP)
 */
void flatten(){

  //Locals
  int n;
  double steep[maxsweep];
  double omega1, omega2, epsilon, delp1, delp2, shock;
  double temp1, temp2;

  // Parameters as suggested in C&W

  omega1 = 0.75;
  omega2 = 5.0;
  //if (lagrangian) omega2  = 5.0;
  //else omega2 = 10.0;
  epsilon = 0.33;

  // Look for presence of a shock using pressure gradient and sign of
  // velocity jump:  shock = 1 if there is a shock in the zone, else shock = 0
  // Compute steepness parameter based on steepness of pressure jump IF 
  // there is a shock.

  for (n = nmin-2; n <= nmax+2; n++){
    //shck[n] = 0;                                   // initialize
    delp1 = p[n+1] - p[n-1];
    delp2 = p[n+2] - p[n-2];
    if(fabs(delp2) < smallvalue) delp2 = smallvalue;
    jump[n] = fabs(delp1)/min(p[n+1],p[n-1]);
    shock = jump[n]-epsilon;                       // Eq 46 in Fryxell etal
    shock = max(0.0,shock);
    if(shock > 0.0) shock = 1.0;
    if(u[n-1] < u[n+1]) shock = 0.0;               // Eq 47 in Fryxell etal
    //if(shock == 1.0) shck[n] = 1;
    temp1 = ( delp1 / delp2 - omega1 ) * omega2;
    steep[n] = shock * max( 0., min(1.0,temp1) );  // Eq 45 including Eqs 46
                                                   // and 47 
  }

  // Set phony boundary conditions for the steepness parameter

  steep[nmin-3] = steep[nmin-2];
  steep[nmax+3] = steep[nmax+2];

  // Set flatening coefficient based on the steepness in neighboring zones
  for (n = nmin-2; n <= nmax+2; n++){
    //if (lagrangian){
      temp2   = max( steep[n-1], max(steep[n], steep[n+1]) );
      flat[n] = max( 0.0, min( 0.5, temp2 ) );
      //flat[n] = max( 1.0, min( 0.5, temp2 ) ); // global diffusion can be added by increasing the
                                               // minimum value of the flattening parameter
      /*}
    else{
      if (p[n+1] < p[n-1]){             // Eq 48 in Fryxell
        flat[n] = max(steep[n],steep[n+1]);
      }
      else flat[n] = max(steep[n],steep[n-1]);
      //flat[n] = max(1.0,flat[n]);   // adjust the numerical value to add/increase numerical diffusion
      }*/
  }

  return;
}


/*!  \brief Calculate volume of cells.
 *
 *   \author Julian Pittard (Original version 05.06.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 05.06.11 (JMP)
 */
void volume(int ngeom){

  int l;
  const float third = 0.333333333;

  // Calculate dvolume and average area based on geometry of sweep

  if (ngeom == 0){
    radius = 1.0;
    for (l = nmin-2; l <= nmax+2; l++){
      dvol [l] = dx [l];
      dvol0[l] = dx0[l];
    }
  }
  else if (ngeom == 1){
    radius = 1.0;
    for (l = nmin-2; l <= nmax+2; l++){
      dvol [l] = dx [l]*(xa [l]+0.5*dx [l]);
      dvol0[l] = dx0[l]*(xa0[l]+0.5*dx0[l]);
    }
  }
  else if (ngeom == 2){
    radius = 1.0;
    for (l = nmin-2; l <= nmax+2; l++){
      dvol [l]=dx [l]*(xa [l]*(xa [l]+dx [l])+dx [l]*dx [l]*third);
      dvol0[l]=dx0[l]*(xa0[l]*(xa0[l]+dx0[l])+dx0[l]*dx0[l]*third);
    }
  }
  else if (ngeom == 3){
    for (l = nmin-2; l <= nmax+2; l++){
      dvol [l] = dx [l]*radius;
      dvol0[l] = dx0[l]*radius;
    }
  }
  else if (ngeom == 4){
    for (l = nmin-2; l <= nmax+2; l++){
      dvol [l] = (cos(xa [l])-cos(xa [l+1]))*radius;
      dvol0[l] = (cos(xa0[l])-cos(xa0[l+1]))*radius;
    }
  }
  else if (ngeom == 5){
    for (l = nmin-2; l <= nmax+2; l++){
      dvol [l] = dx [l] * radius;
      dvol0[l] = dx0[l] * radius;
    }
  }

  return;
}


/*!  \brief Function which calculates fictitious (ie. coriolis and 
 *          centrifugal) forces. 
 *
 *   \var fict[] (double) returns with the fictitious forces.
 *
 *   Actual forces (eg. gravitational and radiative) should be calculated 
 *   elsewhere and extracted into the array act(n) in sweepx.cpp, 
 *   sweepy.cpp and sweepz.cpp
 *
 *   \author Julian Pittard (Original version 05.06.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 04.08.11 (JMP)
 */
void fict_forces(double xf[], double uf[], double vf[], double wf[], double fict[]){

  int n;
  float xf0;

  // Initialize the fictious array

  for (n = 0; n < maxsweep; n++){
    fict[n] = 0.0;
  }

  // Perform the force calculations
 
  if (sweep == "x"){
    if ((ngeomx == 0) && ((ngeomy == 0)||(ngeomy == 1))){

      // Cartesian or cyclindrical radial
      // Do nothing, since already initialized to zero...
    }
    else if ((ngeomx == 1) && (ngeomy == 3)){

      // Cylindrical polar

      for (n = 0; n < maxsweep; n++){
        fict[n] = vf[n]*vf[n]/xf[n];
      }
    }
    else if (ngeomx == 2){

      // Spherical

      for (n = 0; n < maxsweep; n++){
        if (xf[n] == 0) fict[n] = 0.;
        else fict[n] = (wf[n]*wf[n]+vf[n]*vf[n])/xf[n];
      }
    }
  }


  if (sweep == "y"){
     if ((ngeomx == 0) && ((ngeomy == 0)||(ngeomy == 1))){

      // Cartesian or cyclindrical radial (do nothing)

    }
    else if ((ngeomx == 1) && (ngeomy == 3)){

      // Cylindrical polar (2D)

      for (n = 0; n < maxsweep; n++){
        fict[n] = -uf[n]*wf[n]/radius;
      }
    }
    else if ((ngeomy == 3) && (ngeomx == 2)){

      // Spherical equator (2D)

      for (n = 0; n < maxsweep; n++){
        fict[n] = -uf[n]*wf[n]/radius;
      }
    }
    else if ((ngeomy == 4) && (ngeomx == 2)){

      // Spherical

      for (n = 0; n < maxsweep; n++){
        if (xf[n] == 0.0) xf0 = smallvalue;
        else xf0 = xf[n];
        fict[n] = -uf[n]*wf[n]/radius
	  +vf[n]*vf[n]/radius * cos(xf[n])/sin(xf0);
        if (xf[n] == 0.0) fict[n] = 0.0;
      }
    }
  }

  if (sweep == "z"){
    if (ngeomz == 0){

      // Cartesian (do nothing)
      
    }
    else if (ngeomz == 3){

      // Cylindrical

      for (n = 0; n < maxsweep; n++){
        fict[n] = -uf[n]*wf[n]/radius;
      }
    }
    else if (ngeomz == 5){

      // Spherical

      for (n = 0; n < maxsweep; n++){
        fict[n] = -uf[n]*wf[n]/radius
	  -uf[n]*vf[n]/radius * cos(theta)/stheta;
      }
    }
  }

  return;
}
