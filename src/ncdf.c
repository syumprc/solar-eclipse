#include <math.h>

#define LTONE 7.0
#define UTZRO 18.66
#define CONST 1.28

double ncdf(double x)
/*
Returns the standard normal cdf as the area to the left of x, using the
algorithm ALNORM of Hill (1973) written to return a lower tail area only.
The function can be made completely analogous to ALNORM by including `up'
in the function header.  Defined constants are as in Hill (1973).
Reference:  I.D. Hill. 1973.  Algorithm AS 66.  The Normal Integral.
Applied Statistics 22(3):424--427.
*/
{
   int up=0;
   double a,y,z;

   z=x;
   if (z < 0.0) {
      z=(-z);
      up=1-up;
   }
   if (z <= LTONE || (up && z <= UTZRO)) {
      y=0.5*z*z;
      if (z < CONST) a=0.5-z*(0.398942280444 - 0.399903438504*y/
         (y + 5.75885480458 - 29.8213557808/
         (y + 2.62433121679 + 48.6959930692/(y + 5.92885724438))));
      else a=0.398942280385 * exp(-y)/
         (z - 3.8052e-8 + 1.00000615302/
         (z + 3.98064794e-4 + 1.98615381364/
         (z - 0.151679116635 + 5.29330324926/
         (z + 4.8385912808 - 15.1508972451/
         (z + 0.742380924027 + 30.789933034/(z + 3.99019417011))))));
      return (up ? a : 1.0-a);
   } else return (up ? 0.0 : 1.0);
}
