#ifndef _MVN_H_
#define _MVN_H_

/* MVN procedures */
void genz1(int n, double mu[], double **v, double a[], double b[],
	unsigned long nmax, double eps, double *fw, double *dw,
	unsigned long *nw);
double mehc(int n, double x[], double mu[], double **r, double *d);
double mehd(int n, double mu[], double **r, double *s, double *t);

/* univariate normal procedures */
double ncdf(double x);
double nicdf(double p);
double npdf(double z);

/* wrappers */
double w_genz1(int *ndim, int *nind, double *mean, double *cov,
	double *a, double *b);
double w_mehc(int *ndim, int *nind, double *trait, double *mean, double *cov);
double w_mehd(int *ndim, int *nind, double *mean, double *cov,
	double *a, double *b);

#endif /* _MVN_H_ */
