#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
#include "mvn.h"

double mehd(int n, double mu[], double **r, double *s, double *t)
/* Mendell-Elston-Hasstedt evaluation of the MVN distribution. */
{
    int i,j,k;
    double scdf,spdf,tcdf,tpdf;
    double a=0.0,detr=1.0,fac,fac1,fac2,p,tmp,v2=0.0;
    double *v;

    v=dvector(1,n);

    /* standardize the problem */
    for (i=1;i<=n;i++) {
	if (mu[i]) {
	    s[i] -= mu[i];
	    t[i] -= mu[i];
	    mu[i]=0.0;
	}
	if (r[i][i] != 1.0) {
	    tmp=sqrt(r[i][i]);
	    s[i] /= tmp;
	    t[i] /= tmp;
	    r[i][i]=1.0;
	    for (j=1;j<=n;j++) {
		if (j != i) r[i][j] = (r[j][i] /= tmp);
	    }
	}
    }

    /* initialize v[i] */
    for (i=1;i<=n;i++) {
	scdf=ncdf(s[i]);
	tcdf=ncdf(t[i]);
	spdf=npdf(s[i]);
	tpdf=npdf(t[i]);
	a=(spdf-tpdf)/(fac=tcdf-scdf);
	v[i]=(s[i]*spdf-t[i]*tpdf)/fac-a*a+1.0;
    }

    /* evaluate the probability for variable i,
       then condition remainder of integral on i */
    for (p=1.0,i=1;i<=n;i++) {
	scdf=ncdf(s[i]);
	tcdf=ncdf(t[i]);
	spdf=npdf(s[i]);
	tpdf=npdf(t[i]);
	a=(spdf-tpdf)/(fac=tcdf-scdf);
	p*=fac;
	v2=1.0-v[i];
	for (j=i+1;j<=n;j++) {
	    fac1=a*r[i][j];
	    fac2=sqrt(1.0-v2*r[i][j]*r[i][j]);
	    detr*=fac2;
	    s[j]=(s[j]-fac1)/fac2;
	    t[j]=(t[j]-fac1)/fac2;
	    v[j]/=(fac2*fac2);
	    for (k=j+1;k<=n;k++) {
		r[j][k]=r[k][j]=(r[j][k]-v2*r[i][j]*r[i][k])/
		    sqrt((1.0-v2*r[i][j]*r[i][j])*(1.0-v2*r[i][k]*r[i][k]));
		if (fabs(r[j][k]) > 1.0) {
		    fprintf(stderr,"[mehd]: r[%d,%d]=%f\n",j,k,r[j][k]);
		    /* return p; */
		    exit(1);
		}
	    }
	}
    }
    free_dvector(v,1,n);
    /* (*dr)=detr*detr; */
    return p;
}
