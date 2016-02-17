#ifndef  _ILDEFS_H

#define  _ILDEFS_H  1

/* Output from p2c, the Pascal-to-C translator */
/* From input file "ilink.p" */
/* This file cocntains definitions for a modified version of the ILINK program*/
/* The modifications are described in the papers: */
/* R. W. Cottingham, Jr., R. M. Idury, and A. A. Schaffer */
/* Faster Sequential Genetic Linkage Computations */
/* American Journal of Human Genetics, 53(1993), pp. 252-263*/
/* and A. A. Schaffer, S. K. Gupta, K. Shriram, and R. W. Cottingham, Jr. */
/* Avoiding Recomputation in Linkage Analysis, Human Heredity */
/* 44(1994), pp. 225-237. */


#ifndef _COMMONDEFS_H
#include <commondefs.h>
#endif

#define scale           3.0   /*SCALE FACTOR*/

#ifndef byfamily
#define byfamily        false
#endif

#define epsilon         1.0e-6


boolean inconsistent;
extern FILE *outfile, *ipedfile, *datafile, *stream, *speedfile;
/*ILINK*/
extern FILE *final;
boolean penalty;
double penlike;

/*GEMINI*/
 enum {
  go, quit, restart
} continue_;

#endif  /* _ILDEFS_H */
