#include <math.h>

#define LNSQRT2PI 0.9189385332046727417803296
#define ZLOW -87.0

double npdf(double z)
/* Compute the standard normal density. */
{
   double arg;

   /* return exp(-0.5*z*z-LNSQRT2PI); */
   return ((arg=(-0.5*z*z-LNSQRT2PI)) > ZLOW ? exp(arg) : 0.0);
}
