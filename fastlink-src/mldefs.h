#ifndef _MLDEFS_H

#define _MLDEFS_H 1

/* Output from p2c, the Pascal-to-C translator */
/* From input file "mlink.p" */
/* This file contains definitions for a modified version of MLINK */
/* The modifications are described in the papers: */
/* R. W. Cottingham, Jr., R. M. Idury, and A. A. Schaffer */
/* Faster Sequential Genetic Linkage Computations */
/* American Journal of Human Genetics, 53(1993), pp. 252-263 */
/* and A. A. Schaffer, S. K. Gupta, K. Shriram, and R. W. Cottingham, Jr. */
/* Avoiding Recomputation in Linkage Analysis, Human Heredity */
/* 44(1994), pp. 225-237 */

#ifndef _COMMONDEFS_H
#include <commondefs.h>
#endif

#define scale           2.0   /*SCALE FACTOR*/
#define score           true   /*CALCULATE LOD SCORES*/

#ifndef byfamily
#define byfamily        true   /*GIVE LIKELIHOODSS BY FAMILY*/
#endif

#ifndef lodbyfamily
#define lodbyfamily     true
#endif

#ifndef BIGNEGATIVE
#define BIGNEGATIVE     (-150)
#endif

#define epsilon         1.0e-6


boolean zeromale[maxlocus], zerofemale[maxlocus];
/*OTHERS*/
long j, nlocus, which, thissystem;
double  tlike, finish, inc, scorevalue, holdtheta;

int continue_;  /* cgh */

#if defined(DOS)
typedef enum
{
  normalRun ,
  checkpointedRun
}
checkpointType;

checkpointType checkpointStatus;
#endif
#endif

