/*!  \file arwen.cpp
 *   \brief Main driver for ARWEN hydrodynamics code.
 *   
 *   ARWEN stands for Astrophysical Research With Enhanced Numerics.
 *   ARWEN is intended to be a simple-to-use, yet powerful, astrophysical
 *   hydrodynamics code.
 *   
 *   This version of ARWEN includes the bare minimum of features. It is not
 *   parallelized, for instance. It is useful as a basic frame on which 
 *   new code/features can be tested/developed. It is also useful as a 
 *   "starter" code for students. 
 *
 *   At present ARWEN solves the inviscid equations of hydrodynamics on
 *   a uniform grid, using unsplit PLMDE. This is not as robust as PPMLR (and
 *   may require courant numbers substantially below 1.0 for stability - e.g.,
 *   cn = 2.5e-5 in the initial stages of a high pressure SN explosion),
 *   but has the advantage of being able to easily include MHD. 
 *
 *   \warning Any problem-specific ".h" files which use an 
 *      "#ifdef MAIN #else" struture MUST also be "#included" in arwen.cpp 
 *      where "MAIN" is defined (otherwise the code may compile successfully but
 *      will fail to link).
 *
 *   \notes Compilation: make
 *          Run:         ./arwen
 * 
 *   \author Julian Pittard (Original version 25.05.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 18.11.11 (JMP)
 */

#define MAIN

//  Header files:
#include <cstdlib>
#include <limits.h>
#include <iomanip>
#include <iostream>
#include <string>
#include "constants.h"
#include "gas.h"
#include "global.h"
#include "grid.h"
#include "gridlet.h"
#include "star.h"
//#include <gperftools/profiler.h>

using namespace std;

// Function declarations
void quit();
void timing(bool init);
extern void computeLoop();
extern void prin1D();
extern void prin2D();
extern void prin2Dfits();
extern void setup();

int main(int argc, char *argv[]){
  // Locals

  // Setup the initial conditions and write to disk
  setup();
  if (nd == 1) prin1D();
  else if (nd == 2){
    //prin2D();
    prin2Dfits();
  }
  
  // Evolve...
  //ProfilerStart("arwen.prof");
  computeLoop();
  //ProfilerStop();
  
  // Write to disk the final solution
  //prin1D();
  
  quit();
}


/*!  \brief Exits out of programme.
 *
 *   Serial code simply stops, while parallel code also exits MPI.
 *
 *   \author Julian Pittard (Original version 02.06.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 02.06.11 (JMP)
 */
void quit(){

  exit(EXIT_SUCCESS);

}


/*!  \brief Times the compute loops.
 *
 *   \author Julian Pittard (Original version 02.06.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 02.06.11 (JMP)
 */
void timing(bool init)
{

  double elapsed_time = 0.0;

  if (init) time0=clock();
  else{
    time1=clock()-time0;
    elapsed_time = (double)time1 / ((double)CLOCKS_PER_SEC);
    if (procRank == 0){
      ec = ncycle;
      double cycps = double(ec-sc)/elapsed_time;
      double cellps = cycps*ncells;
      cout << setprecision(4) << "That took " << elapsed_time << " seconds; total number of steps = " << ncycle << "; cycles/s = " << cycps << "; cells/s = " << cellps << "\n";
    }
  }

  return;
}
