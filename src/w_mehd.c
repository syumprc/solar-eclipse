#include <math.h>
#include <stdio.h>
#include "nrutil.h"
#include "mvn.h"

double w_mehd_(int *ndim, int *nind, double *mean, double *cov,
	double *a, double *b)
/*
This is an F77 function call wrapper for the C procedure `mehd'.
*/
{
    int i,j,n,nd;
    double *mu,p,**r,*s,*t;

    n=(*nind);
    nd=(*ndim);

    /* malloc C arrays */
    mu=dvector(1,n);
    r=dmatrix(1,n,1,n);
    s=dvector(1,n);
    t=dvector(1,n);

    /* load mean vector and vector limits of integration */
    for (i=1;i<=n;i++) {
	mu[i]=mean[i-1];
	s[i]=a[i-1];
	t[i]=b[i-1];
    }

    /* load covariance matrix */
    for (j=1;j<=n;j++)
	for (i=1;i<=n;i++)
	    r[j][i]=cov[(i-1)+(j-1)*nd];

    /* get MVN distribution */
    p=mehd(n,mu,r,s,t);

    free_dvector(t,1,n);
    free_dvector(s,1,n);
    free_dmatrix(r,1,n,1,n);
    free_dvector(mu,1,n);

    /* return ln-likelihood */
    return log(p);
}
