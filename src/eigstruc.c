#include "safelib.h"
#include <math.h>

extern void symeig_ (int*, double*, double*, double*, double*, int*);

void eigstruc_ (double *x, int *n, int *info)
{
    int i, j, nd = (*n);
    double* z = (double*) Calloc (nd*nd, sizeof (double));
    double* d =  (double*) Calloc (nd, sizeof (double));
    double* e =  (double*) Calloc (nd, sizeof (double));

    symeig_(n, x, d, e, z, info);
    if (*info) return;

    for (i = 0; i < nd; i++) {
        for (j = 0; j < nd; j++) {
            if (d[j] < 0) d[j] = 0;
            x[i*nd+j] = z[j*nd+i] * sqrt(d[j]);
        }
    }

    free (e);
    free (d);
    free (z);
}
