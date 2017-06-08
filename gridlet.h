/*!  \file gridlet.h
 *   \brief This is a repository for the gridlet variables and constants.
 *
 */

using namespace std;

// Constants...
const int nig = 1; // number of gridlets in i-direction
const int njg = 1; // number of gridlets in j-direction
const int nkg = 1; // number of gridlets in k-direction
 
const int ig = imax/nig; // number of real cells on the gridlet in the i-direction
const int jg = jmax/njg; // number of real cells on the gridlet in the j-direction
const int kg = kmax/nkg; // number of real cells on the gridlet in the k-direction

// Total number of zones on the gridlet[kk][jj][ii]
const int ii = ig + 2*nghost;
#if NDIM > 1
const int jj = jg + 2*nghost;
#else
const int jj = 1;
#endif
#if NDIM == 3
const int kk = kg + 2*nghost;
#else
const int kk = 1;
#endif

// Variables
#ifdef MAIN
 int irs,ire,jrs,jre,krs,kre;  
#else
 extern int irs,ire,jrs,jre,krs,kre;
#endif

