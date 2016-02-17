#ifndef  _GEMDEFS_H

#define  _GEMDEFS_H  1

/* Output from p2c, the Pascal-to-C translator */
/* From input file "ilink.p" */
/* This file contains definitions for a modified version of */
/* the gemini package used in the LODSCORE and ILINK programs*/
/* The modifications are described in the papers: */
/* R. W. Cottingham, Jr., R. M. Idury, and A. A. Schaffer, */
/* Faster Sequential Genetic Linkage Computations */
/* American Journal of Human Genetics, 53(1993), pp. 252-263 */
/* and A. A. Schaffer, S. K. Gupta, K. Shriram, and R. W. Cottingham, Jr. */
/* Avoiding Recomputation in Linkage Analysis, Human Heredity */
/* 44(1994), pp. 225-237. */

#include "commondefs.h"

/*GEMINI*/

#if !defined(iterationMultiple)
/* Maximum number of iterations will be
   iterationMultiple * n; n is number of dimensions */
#define iterationMultiple  50   
#endif

#ifndef maxn
#define maxn            20   /* MAXIMUM NUMBER OF ITERATED PARAMETERS */
#endif

#define nbit            23   /*NUMBER OF BITS OF MACHINE PRECISION*/

#define tolconst        1.0e-3   /*TOLERANCE ON LIKELIHOOD FOR STOP*/
/*DO NOT CHANGE THE FOLLOWING CONSTANTS*/
#define tolg            1.0e-35
#define xpmcon          1.0e30
#define clbcon          1.0e10
#define tbndcon         1.0e5
#define curv            0.75
#define small           1.0e-10

/*GEMINI*/
typedef int itertype[maxn];
typedef double vector[maxn];
typedef vector matrix[maxn];

/*VARIABLES FROM COMMON MIN1*/
 matrix tl;
 vector d, g, gsave, y, p;
#if PARALLEL
 vector *gg;   /*Tmk_malloc later, Sandeep*/
#endif

/*VARIABLE FROM COMMON MIN2*/
 int nit, nfe, idg, idif, isw, iret, ibnd, ihess, ivar, ihx, maxit;
 double tol, tmin, h, trupb, ptg;
#if PARALLEL
 int *gnfe;  /*Sandeep*/
#endif

/*DIMENSIONED IN GEMINI*/
 vector xall, x, v, se;
#if PARALLEL
 vector *gx, *gxall;  /*Tmk_malloc later, Sandeep*/
#endif
 itertype itp, holditp;
#if PARALLEL
 itertype *gitp;      /*Tmk_malloc later, Sandeep*/
#endif
 double bnd[maxn][2];
 int nall, n, icall;
 double tbnd, f, fsmf, fsav2, t, hx, xsave, fxph, fxmh, xp, xpm, ytp;

/*ARRAYS IN UPDATE*/
 vector wtil, ztil, w, z, wtjp1, ztjp1, s, dp;

 double nu, muj, lambj, lambj2, sb, sc, sd, fbcd, alpha, sa, thet1,
	      thet2, aa, bb, cc, del2, alph1, alph2, rdiv, eta, aj, thj, bj,
	      gamlj, betlj, del;

 int iswup;

/*FROM STEP*/
 vector xt;

 double fsave, sumt, twot, ft, f2t, ft2, scal;

/*FROM CHKBND */
 double clb, check, eh, teq;

/*OTHERS*/
 int itsys;
/*USUALLY IN INIB*/
 matrix bmat;

 boolean active;

/*Added to avoid extra function evaluations*/
vector outsavex;
double outsavefvalue;
double savedf[2 * maxn];
double likebyped[maxped];
double outsavelike[maxped];

#if PARALLEL

/*Prototypes */

#if defined(KNR_PROTO)
void OneMaster();
void AssignThetas();
unsigned MakeMask();
#else
void OneMaster(void);
void AssignThetas(void);
unsigned MakeMask(int start, int end);
#endif  /* defined(KNR_PROTO) */

#endif
#endif  /* _GEMDEFS_H */





