/*!  \file funcs.h
 *   \brief This is a repository for inline functions.
 *
 *   \author Julian Pittard (Original version 04.08.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 04.08.11 (JMP)
 */

#define small 1.0e-30

inline float min(float a, float b)
{
 return (a < b ? a : b);
}

inline float max(float a, float b)
{
 return (a > b ? a : b);
}

inline int min(int a, int b)
{
 return (a < b ? a : b);
}

inline int max(int a, int b)
{
 return (a > b ? a : b);
}

inline float av(float a, float b)
{
  return(b*a < small ? 0 : a*b*(a + b)/ (a*a + b*b));
}

inline int ipow (int j)
{
  return(1 << j);
}

inline float mabs(float a)
{
  return(a > 0.0 ? a : -a);
}

inline double min(double a, double b)
{
 return (a < b ? a : b);
}

inline double max(double a, double b)
{
 return (a > b ? a : b);
}

inline double av(double a, double b)
{
  return(b*a < small ? 0 : a*b*(a + b)/ (a*a + b*b));
}

inline double mabs(double a)
{
  return(a > 0.0 ? a : -a);
}
