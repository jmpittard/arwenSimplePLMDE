/*!  \file stars.h
 *   \brief This is a repository for the Star class.
 *
 */

using namespace std;

const int nstars = 2;

struct Star{
  double mdot; // mass-loss rate (g/s)
  double vinf; // wind speed (cm/s)
  double xpos; // x-position (cm)
  double ypos; // y-position (cm)
  double zpos; // z-position (cm)
};

// Variables
#ifdef MAIN
Star stars[nstars];
double dsep = 0.0;
double rob = 0.0;
double eta = 0.0;
#else
extern Star stars[];
extern double dsep,rob,eta;
#endif

