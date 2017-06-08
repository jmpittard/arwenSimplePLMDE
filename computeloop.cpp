/*!  \file computeloop.cpp
 *   \brief Main computational loop.
 * 
 *   \author Julian Pittard (Original version 01.06.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 01.06.11 (JMP)
 */

#undef MAIN

//  Header files:

#include <iostream>
#include "constants.h"
#include "global.h"
#include "grid.h"


// Function declarations
extern bool computeTimestep();
extern void prin1D();
extern void prin2D();
extern void prin2Dfits();
extern void quit();
extern void timing(bool init);
extern void update();

using namespace std;


/*!  \brief Main computational loop.
 * 
 *   Marches the solution forward in time until the desired number of steps
 *   or the target time is achieved. During this period files are output
 *   as desired and the safe timestep is computed. The actual solution update
 *   is computed in update(). Problem-specific code after each update is
 *   handled in endofcomputeloop().
 *
 *   \author Julian Pittard (Original version 01.06.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 24.08.11 (JMP)
 */
void computeLoop(){

  // Locals
  bool loop = true;
  int step = 0;
  double smalltime = 1.0e-6;
  double dtinit;

  // Begin timing of compute loop
  timing(true);
  sc = ncycle;
  
  while (loop){

    //Increase counters
    ncycle++;
    ncycp1d++;
    ncycp2d++;
    ncycp3d++;

    // Adjust timestep if needed (2 timesteps are performed per cycle)
    if (t + 2.0*dt > endtime){       //set dt to land on endtime
      dt = 0.5*(endtime - t);
    }
    if ( (timep1d + 2.0*dt - tprin1d)/tprin1d > -smalltime ){ // set dt to land on tprin
      dt = 0.5*(tprin1d - timep1d);
      ncycp1d = nprin1d;  //changing the counters ensures file output
    }
    if ( (timep2d + 2.0*dt - tprin2d)/tprin2d > -smalltime ){ 
      dt = 0.5*(tprin2d - timep2d);
      ncycp2d = nprin2d;
    }
    if ( (timep3d + 2.0*dt - tprin3d)/tprin3d > -smalltime ){ 
      dt = 0.5*(tprin3d - timep3d);
      ncycp3d = nprin3d;
    }

    // Update the fluid variables
    update();

    // Output marching to terminal
    if (ncycle%ncycshow == 0){
        if (procRank == 0) cout << "ncycle = " << ncycle << "; time = " << t << "; dt = " << dt << "\n";
    }
    
    // Dump files
    if ((ncycp1d >= nprin1d) || (timep1d >= tprin1d )){
      //cout << "ncycp1d = " << ncycp1d << "; nprin1d = " << nprin1d << "; timep1d = " << timep1d << "; tprin1d = " << tprin1d << "\n";
      //quit();
      timep1d = 0.;
      ncycp1d = 0;
      prin1D();
    }
    if ((ncycp2d >= nprin2d) || (timep2d >= tprin2d )){
      timep2d = 0.;
      ncycp2d = 0;
      //prin2D();
      prin2Dfits();
    }
    if ((ncycp3d >= nprin3d) || (timep3d >= tprin3d )){
      timep3d = 0.;
      ncycp3d = 0;
      // prin3D();
    }
    
    // Barrier is necessary here to prevent some processes exiting before
    // all the files are written

    //barrier();

    // Calculate shock radius and write out
    //fshk_radius = findForwardShock();
    //prinShkRadius();
    
    // Update dt
    computeTimestep();
      
    // Check if the loop is finished

    if ((ncycle >= ncycend)||(t >= endtime)){
      if (procRank==0) cout << "Target reached, total number of steps = " << ncycle << "\n";
      loop = false;
    }
    
  } // do while

  // Output the elapsed time of the compute loop
  timing(false);

  return;
}





