/*!  \file grid.h
 *   \brief This is a repository for the grid variables and constants.
 *
 *   \param nc (const int array) stores the dimensionality of the grid.
 *   \param geom (string) specifies the grid geometry (XYZ or ZRP)
 *   \param ngeomx (int) stores the geometry of the 1st grid dimension
 *   \param ngeomy (int) stores the geometry of the 2nd grid dimension
 *   \param ngeomz (int) stores the geometry of the 3rd grid dimension
 *   \param nleftx (int) stores the left grid boundary type of the 1st dimension
 *   \param nlefty (int) stores the left grid boundary type of the 2nd dimension
 *   \param nleftz (int) stores the left grid boundary type of the 3rd dimension
 *   \param nrightx (int) stores the right grid boundary type of the 1st dimension
 *   \param nrighty (int) stores the right grid boundary type of the 2nd dimension
 *   \param nrightz (int) stores the right grid boundary type of the 3rd dimension
 *   \param xmin (double) is the left edge of the 1st grid dimension
 *   \param ymin (double) is the left edge of the 2nd grid dimension
 *   \param zmin (double) is the left edge of the 3rd grid dimension
 *   \param xmax (double) is the right edge of the 1st grid dimension
 *   \param ymax (double) is the right edge of the 2nd grid dimension
 *   \param zmax (double) is the right edge of the 3rd grid dimension
 *   \param P0 (<double,4>) stores the cell averaged primitives
 *   \param Pstar (<double,4>) stores the cell averaged primitives at the end of the predictor step
 *   \param zxa (<double>) stores the left edge coord of cells in the 1st dimension
 *   \param zya (<double>) stores the left edge coord of cells in the 2nd dimension
 *   \param zza (<double>) stores the left edge coord of cells in the 3rd dimension
 *   \param zxc (<double>) stores the coord of cell centres in the 1st dimension
 *   \param zyc (<double>) stores the coord of cell centres in the 2nd dimension
 *   \param zzc (<double>) stores the coord of cell centres in the 3rd dimension
 *   \param zxg (<double>) stores the coord of cell centre of mass in the 1st dimension
 *   \param zyg (<double>) stores the coord of cell centre of mass in the 2nd dimension
 *   \param zzg (<double>) stores the coord of cell centre of mass in the 3rd dimension
 *   \param zdx (<double>) stores the width of cells in the 1st dimension
 *   \param zdy (<double>) stores the width of cells in the 2nd dimension
 *   \param zdz (<double>) stores the width of cells in the 3rd dimension
 * 
 *   \author Julian Pittard (Original version 31.05.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 18.08.11 (JMP)
 */

#include <string>      

using namespace std;

// Constants. Most are defined in constants.h
const int nc[3] = {imax,jmax,kmax}; //  size of grid

// Struct for grid local to each processor
struct lgrid{
  string geom;      // specify geometry (XYZ or ZRP)
  int ngeomx, ngeomy, ngeomz, nleftx, nlefty, nleftz, nrightx, nrighty, nrightz;

  double xmin;   ///< min coord of 1st axis on local grid (real cells ONLY). global grid extent is defined in structure ggrid
  double xmax;   ///< max coord of 1st axis on local grid (real cells ONLY)
  double ymin;   ///< min coord of 2nd axis on local grid (real cells ONLY)
  double ymax;   ///< max coord of 2nd axis on local grid (real cells ONLY)
  double zmin;   ///< min coord of 3rd axis on local grid (real cells ONLY)
  double zmax;   ///< max coord of 3rd axis on local grid (real cells ONLY)

  // Construct grid arrays
  double P0   [iqmax][kmaxd][jmaxd][imaxd]; // cell averaged primitive variables
  double Pstar[iqmax][kmaxd][jmaxd][imaxd]; 

  double zxa[imaxd+1]; // left edges of cells in x-direction
  double zxc[imaxd+1]; // centre 
  double zxg[imaxd+1]; // centre of mass
  double zdx[imaxd+1]; // width
  double zya[jmaxd+1]; // left edges of cells in y-direction
  double zyc[jmaxd+1]; // centre 
  double zyg[jmaxd+1]; // centre of mass
  double zdy[jmaxd+1]; // width
  double zza[kmaxd+1]; // left edges of cells in z-direction
  double zzc[kmaxd+1]; // centre 
  double zzg[kmaxd+1]; // centre of mass
  double zdz[kmaxd+1]; // width

  double Pbc[iqmax][3][2]; // for storing fixed boundary conditions: (iq,ijk,lr)
  
  int irs, ire, jrs, jre, krs, kre; // inclusive range of real cells on the grid
};

// Variables
#ifdef MAIN
lgrid lg; // grid local to each processor
//ggrid gg; // global grid
#else
extern lgrid lg;
//extern ggrid gg;
#endif

