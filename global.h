/*!  \file global.h
 *   \brief  This is a repository for global variables and constants.
 *
 *   \param problem (string) defines the problem being tackled
 *   \param ncycle (int) specifies the current cycle number
 *   \param t1 (double) specifes the start time of the current march
 *   \param t2 (double) specifes the end time of the current march
 *   \param t (double) specifes the current time (s) of the calculation
 *   \param dt (double) specifes the current timestep (s) of the calculation
 *   \param svel (double) specifes the sound speed
 *   \param ramb (double) is the ambient density (g/cm^-3)
 *   \param pamb (double) is the ambient pressure (dyn/cm^-2)
 *   \param tamb (double) is the ambient temperature (K)
 *   \param uamb (double) is the ambient velocity (cm/s)
 *   \param camb (double) is the ambient sound speed (cm/s)
 *   \param avgmassamb (double) is the ambient average mass (g)
 *   \param muamb (double) is the ambient average mass (in units of mH)
 *   \param tmin (double) is the temperature floor (K)
 *   \param tmax (double) is the maximum permitted (K)
 *   \param dotflo (double) is a fixed density outflow at the grid boundary
 *
 *   \author Julian Pittard (Original version 31.05.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 13.09.11 (JMP)
 */

#include <time.h>
#include <iostream>
#include <string>

using namespace std;

// Variables
#ifdef MAIN
  clock_t time0, time1;       // compute loop timing
  bool reStart = false;       // reStart the current problem?
  string problem;             // the problem being tackled
  string prefx = "cwb";       // prefix to dump files
  int procRank = 0;           // processor rank
  int ncycle;                 // the current number of hydro steps taken
  int ncells;                 // the number of cells on the grid
  int nprin1d = 1e9;          // cycles between output of 1D dump files 
  int nprin2d = 1e9;          // cycles between output of 2D dump files
  int nprin3d = 1e9;          // cycles between output of 3D dump files
  int ncycp1d,ncycp2d,ncycp3d;// elapsed steps since last 1D/2D/3D output
  int ncycshow = 100;          // number of steps before outputing cycle
  int ncycend = 1e9;          // number of steps before ending code
  int nfile1d;                // file index for 1d dumps
  int nfile2d;                // file index for 2d dumps
  int nfile3d;                // file index for 3d dumps
  int sc,ec;                  // start and end cycles when timing
  int remapi;                 // number of cells in remap radius
  double tprin1d = 1.0e99;    // time between 1D dumps
  double tprin2d = 1.0e99;    // time between 2D dumps
  double tprin3d = 1.0e99;    // time between 3D dumps
  double timep1d,timep2d,timep3d; // elapsed time since last 1D/2D/3D output
  double endtime = 1.0e99;    // endtime
  double t;                   // current time
  double dt;                  // current timestep
  double svel;                // sound speed
  double cn = 0.25;            //courant number
  double ramb,pamb,tamb,uamb,vamb,wamb,camb,avgmassamb,muamb,tmin,tmax;
  double dotflo,uotflo,votflo,wotflo,potflo,dinflo,uinflo,vinflo,winflo,pinflo,einflo;
  double asotflo[smax],asinflo[smax];
  double remapRadius;
  double mdot,vinf,esn,msn;
#else
extern clock_t time0, time1; 
extern bool reStart;
extern string problem,prefx;
extern int procRank,ncycle,ncells,nprin1d,nprin2d,nprin3d,ncycp1d,ncycp2d,ncycp3d,nfile1d,nfile2d,nfile3d,ncycend,ncycshow,sc,ec,remapi;
extern double endtime,t,dt,timep1d,timep2d,timep3d,tprin1d,tprin2d,tprin3d;
extern double svel,cn;
extern double ramb,pamb,tamb,uamb,vamb,wamb,camb,avgmassamb,muamb,tmin,tmax;
extern double dotflo,uotflo,votflo,wotflo,potflo,dinflo,uinflo,vinflo,winflo,pinflo,einflo;
extern double asotflo[],asinflo[];
extern double remapRadius;
extern double mdot,vinf,esn,msn;
#endif
