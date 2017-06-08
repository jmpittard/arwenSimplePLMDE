/*!  \file constants.h
 *   \brief This is a repository for global constants.
 *
 *   These set/define the sizes of arrays, etc. This file must be
 *   "#included" in any ".cpp" which in turn "#includes" any other 
 *   ".h" files which uses arrays making use of the values below.
 *
 *   \param nd (const int) defines the dimensionality of the grid.
 *   \param imax (const int) is the number of zones in 1st direction.
 *   \param jmax (const int) is the number of zones in 2nd direction.
 *            If nd = 1, jmax MUST equal 1.
 *   \param kmax (const int) is the number of zones in 3rd direction.
 *            If nd < 3, kmax MUST equal 1.
 *   \param smax (const int) is the total number of advected scalars (min=1)
 *   \param cmax (const int) is the number of advected colours (min=1)
 *   \param nghost (const int) is the number of ghost cells to use
 *   \param cn (const double) is the courant number
 *   \param smallvalue (const double) specifes what stands for a small value in the code
 *   \param debug_lev (const int) sets the debug level for non-MPI code
 *            ( 0 = no messages, 1 = some messages, 3 = more messages, 5 = ALL messages)
 *   \param debug_mpi (const bool) sets the debug option for MPI code 
 *            (false = no messages, true = debug messages)
 *
 *   \author Julian Pittard (Original version 31.05.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 31.05.11 (JMP)
 */

  //////////////////////////////////////////////////////////////////////
  //                                                                  //
  //                   SET UP CODE PARAMETERS                         //
  //     ...students can (carefully) alter these values...            //
  //                                                                  //
  //       NOTE: Remember to recompile after making changes!          // 
  //                                                                  //
  //////////////////////////////////////////////////////////////////////

#define NDIM 2                // First, define the number of dimensions

#if NDIM == 1
const int nd   = 1;           // dimensionality of grid
#elif NDIM == 2
const int nd   = 2;           // dimensionality of grid
#else
const int nd   = 3;           // dimensionality of grid
#endif
const int imax = 200;         // number of zones in 1st direction
const int jmax = 200;         // number of zones in 2nd direction
const int kmax = 1;           // number of zones in 3rd direction
const int smax = 1;           // total number of advected scalars (min=1)
const int cmax = 1;           // number of advected colours (min=1)

  //////////////////////////////////////////////////////////////////////
  //                                                                  //
  //    STUDENTS SHOULD NOT CHANGE THE CODE BELOW THIS POSITION       //
  //                                                                  //
  //////////////////////////////////////////////////////////////////////

const int iqmax = 2+nd+smax;  // total number of hydro variables per cell (hydro: density, pressure, nd velocity, smax)
const int nghost = 2;         // number of ghost zones at each boundary, surrounding the grid
const int imaxd = imax + 2*nghost; // total number of zones in i-direction
#if NDIM > 1
const int jmaxd = jmax + 2*nghost; // total number of zones in j-direction
#else
const int jmaxd = 1;
#endif
#if NDIM == 3
const int kmaxd = kmax + 2*nghost; // total number of zones in k-direction
#else
const int kmaxd = 1;
#endif

const double smallvalue = 1.0e-40;  // small value in code

const int debug_lev = 0;      // Set debug level for non-MPI code
                              // 0 = no messages
                              // 1 = some messages
                              // 3 = more messages
                              // 5 = ALL messages
const bool debug_mpi = false;  // false = no messages, true = debug messages

const int iqd = 0;  // density
const int iqu0 = 1; // velocity0/mtm0 (in 2D the velocity components are iqu0 and iqu0+1)
const int iqe = iqu0+nd; // pressure/energy
const int iqal0 = iqe+1; // advected scalar 0 (with smax=2, the scalar components are iqal0 and iqal0+1)

const double avisp = 0.2; // artificial viscosity    
const double avise = 0.2; // artificial conductivity    
