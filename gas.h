/*!  \file gas.h
 *   \brief  This is a repository for gas variables and constants.
 *
 *   \author Julian Pittard (Original version 31.05.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 31.05.11 (JMP)
 */

// Constants
const double gam = 1.666666666667;  // gamma
const double gamm = gam - 1.0;      
const double dfloor = 1.0e-40;      // density floor in code (g/cm^3)
const double pfloor = 1.0e-22;      // pressure floor in code (dyn/cm^2)
//const double tfloor = 1.0e4;        // temperature floor in code (K) - currently not used (use Tmin instead)
const double avgmass = 1.0e-24;     // average particle mass (g)
const double g1 = gam - 1.0;
const double g2 = gam/g1;
const double g3 = 0.5*(gam + 1.0)/gam;
const double g4 = 0.5*g1/gam;
const double g5 = 1.0/g4;
const double g6 = 1.0/gam;
const double g7 = 2.0/g1;
const double g8 = g1/(gam + 1.0);

// Variables
#ifdef MAIN
  double avgmass1;
  double avgmass2;
  double avgmass3;
  double gg_Alpha = ((gam + 1.0)/4.0);
  double gg_Beta1 = ((gam - 1.0)/2.0);
  double gg_Beta2 = (2.0/(gam - 1.0)); //inverse
  double gg_Delta2 = (gam - 1.0)/(2.0*gam); //inverse
#else
extern double avgmass1, avgmass2, avgmass3;
extern double gg_Alpha, gg_Beta1, gg_Beta2, gg_Delta2;
#endif
