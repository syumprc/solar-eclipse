// Determines closeness in approximate "significant digits" 0-100

// UPDATE: this now uses isfinite instead of the original finite
// to be compatible with -std=c++0x

#include <math.h>

#ifdef __SUNPRO_CC
#include <ieeefp.h>
#endif

#ifdef GCC_ON_SUN
#include <ieeefp>
#endif

int howclose (const double d1, const double d2)
{
    if (d1 == d2) return 100;

    if (!isfinite(d1) || !isfinite(d2)) return 100;

    double ad1 = fabs(d1);
    double ad2 = fabs(d2);

    if (ad1 > 2*ad2 || ad2 > 2*ad1) return 0;

    double delta = fabs (d1 - d2);
    double smallest_magnitude = (ad1 > ad2) ? ad2 : ad1;
    if (smallest_magnitude == 0.0) return 0;


    double fraction = delta / smallest_magnitude;
    double logarithm = - log10 (fraction);
    if (logarithm < 1.0) return 0;
    return (int) floor (logarithm);

}

extern "C" int howclose_ (const double *d1, const double *d2)
{
    return howclose (*d1, *d2);
}

#ifdef HOWCLOSE_STANDALONE
#include <iostream>
int main () {
    for (;;) {
	double d1, d2;
	std::cout << "enter number 1: ";
	std::cin >> d1;
	std::cout << "enter number 2: ";
	std::cin >> d2;
	std::cout << "They are the same for " << howclose(d1,d2) << " digits.\n";
    }
    return 0;
}
#endif
