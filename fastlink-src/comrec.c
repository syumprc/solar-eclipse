/* Output from p2c, the Pascal-to-C translator */
/* From input file "linkmap.p" */
/* This file contains input code shared by modified versions of the */
/* LODSCORE, ILINK,  LINKMAP, and MLINK programs */
/* The modifications are described in the papers: */
/* R. W. Cottingham, Jr., R. M. Idury, and A. A. Schaffer, */
/* Faster Sequential Genetic Linkage Computations */
/* American Journal of Human Genetics, 53(1993), pp. 252-263 */
/* and A. A. Schaffer, S. K. Gupta, K. Shriram, and R. W. Cottingham, Jr. */
/* Avoiding Recomputation in Linkage Analysis*/
/* Human Heredity 44(1994), pp. 225-237.  */
/* Most of the  code  in this file is not changed from LINKAGE 5.1 */
/* The code in this file computes the recombination probabilities*/

#include "commondefs.h"
#if !defined(DOS)
#include "checkpointdefs.h"
#endif
#ifdef LODSCORE
#include "lodefs.h"
#endif

/* Local variables for recombination: */
struct LOC_recombination {
  int here;
  double p1, p2, p3, p4;
} ;

/* Local variables for recombine: */
struct LOC_recombine {
  struct LOC_recombination *LINK;
  double *theta;
  double *segprob;
  int there, nhap;
  boolean thishet[maxlocus];
  hapvector hap1, hap2;
  hapvector thishap1[maxseg], thishap2[maxseg];
} ;

Local Void scramble(LINK)
struct LOC_recombine *LINK;
{
  int whichhap, start, length, i, j, k, nrec;
  double recval, val;
  int FORLIM2;

  start = 0;
  do {
    start++;
  } while (!LINK->thishet[start - 1]);
  length = LINK->there - LINK->LINK->here;
  for (i = 1; i < length; i++) {
    memcpy(LINK->hap1, LINK->thishap1[i], sizeof(hapvector));
    for (j = 1; j <= length; j++) {
      nrec = 0;
      val = 0.5;
      whichhap = 1;
      recval = LINK->theta[start - 1];
#ifdef LODSCORE
      FORLIM2 = nlocus;
#else
      FORLIM2 = mlocus;
#endif
      for (k = start; k < FORLIM2; k++) {
	if (!LINK->thishet[k])
	  recval = recval * (1.0 - LINK->theta[k]) +
		   LINK->theta[k] * (1.0 - recval);
	else {
	  if (whichhap == 1) {
	    if (LINK->thishap1[j - 1][k] == LINK->hap1[k])
	      val *= 1 - recval;
	    else {
	      nrec++;
	      val *= recval;
	      whichhap = 2;
	    }
	  } else {
	    if (LINK->thishap2[j - 1][k] == LINK->hap1[k])
	      val *= 1 - recval;
	    else {
	      nrec++;
	      val *= recval;
	      whichhap = 1;
	    }
	  }
	  recval = LINK->theta[k];
	}
      }
      LINK->there++;
      if (nrec > maxrec)
	LINK->segprob[LINK->there - 1] = 0.0;
      else
	LINK->segprob[LINK->there - 1] = val;
    }
  }
}

Local Void setrec(val, LINK)
double val;
struct LOC_recombine *LINK;
{
  LINK->nhap++;
  memcpy(LINK->thishap1[LINK->nhap - 1], LINK->hap1, sizeof(hapvector));
  memcpy(LINK->thishap2[LINK->nhap - 1], LINK->hap2, sizeof(hapvector));
  LINK->there++;
  LINK->segprob[LINK->there - 1] = val;
}

Local Void dointer(LINK)
struct LOC_recombine *LINK;
{
  int i;
  boolean temphet[3];

  for (i = 0; i < mlocus; i++)
    temphet[i] = LINK->thishet[i];

  if (temphet[0] && temphet[1] && !temphet[2]) {
    setrec(0.5 - 0.5 * LINK->theta[0], LINK);
    setrec(0.5 * LINK->theta[0], LINK);
    setrec(0.5 * LINK->theta[0], LINK);
    setrec(0.5 - 0.5 * LINK->theta[0], LINK);
    return;
  }
  if (temphet[2] && temphet[1] && !temphet[0]) {
    setrec(0.5 - 0.5 * LINK->theta[mlocus - 2], LINK);
    setrec(0.5 * LINK->theta[mlocus - 2], LINK);
    setrec(0.5 * LINK->theta[mlocus - 2], LINK);
    setrec(0.5 - 0.5 * LINK->theta[mlocus - 2], LINK);
    return;
  }
  if (temphet[2] && temphet[0] && !temphet[1]) {
    setrec(0.5 - 0.5 * LINK->theta[mlocus - 1], LINK);
    setrec(0.5 * LINK->theta[mlocus - 1], LINK);
    setrec(0.5 * LINK->theta[mlocus - 1], LINK);
    setrec(0.5 - 0.5 * LINK->theta[mlocus - 1], LINK);
    return;
  }
  if (!(temphet[0] && temphet[1] && temphet[2])) {
    setrec(0.5, LINK);
    return;
  }
  setrec(0.5 * LINK->LINK->p4, LINK);
  setrec(0.5 * LINK->LINK->p2, LINK);
  setrec(0.5 * LINK->LINK->p3, LINK);
  setrec(0.5 * LINK->LINK->p1, LINK);
  setrec(0.5 * LINK->LINK->p2, LINK);
  setrec(0.5 * LINK->LINK->p4, LINK);
  setrec(0.5 * LINK->LINK->p1, LINK);
  setrec(0.5 * LINK->LINK->p3, LINK);
  setrec(0.5 * LINK->LINK->p3, LINK);
  setrec(0.5 * LINK->LINK->p1, LINK);
  setrec(0.5 * LINK->LINK->p4, LINK);
  setrec(0.5 * LINK->LINK->p2, LINK);
  setrec(0.5 * LINK->LINK->p1, LINK);
  setrec(0.5 * LINK->LINK->p3, LINK);
  setrec(0.5 * LINK->LINK->p2, LINK);
  setrec(0.5 * LINK->LINK->p4, LINK);

  /*not informative*/
}

Local Void nexthet(i, val, inphase, LINK)
int i;
double val;
boolInt inphase;
struct LOC_recombine *LINK;
{
  double newval, recval;

  recval = LINK->theta[i - 1];
  do {
    i++;
    LINK->hap1[i - 1] = 0;
    LINK->hap2[i - 1] = 0;
    if (!LINK->thishet[i - 1] && i != mlocus)
      recval = recval * (1 - LINK->theta[i - 1]) +
	       (1 - recval) * LINK->theta[i - 1];
  } while (!(i == mlocus || LINK->thishet[i - 1]));
  if (i != mlocus) {
    if (inphase)
      newval = val * (1 - recval);
    else
      newval = val * recval;
    LINK->hap1[i - 1] = 1;
    LINK->hap2[i - 1] = 2;
    nexthet(i, newval, true, LINK);
    LINK->hap2[i - 1] = 1;
    LINK->hap1[i - 1] = 2;
    if (!inphase)
      newval = val * (1 - recval);
    else
      newval = val * recval;
    nexthet(i, newval, false, LINK);
    return;
  }
  if (!LINK->thishet[i - 1]) {
    setrec(val, LINK);
    return;
  }
  if (inphase)
    newval = val * (1 - recval);
  else
    newval = val * recval;
  LINK->hap1[i - 1] = 1;
  LINK->hap2[i - 1] = 2;
  setrec(newval, LINK);
  if (!inphase)
    newval = val * (1 - recval);
  else
    newval = val * recval;
  LINK->hap2[i - 1] = 1;
  LINK->hap1[i - 1] = 2;
  setrec(newval, LINK);
}


Local Void getrecprob(LINK)
struct LOC_recombine *LINK;
{
  int i;

  LINK->nhap = 0;
  LINK->there = LINK->LINK->here;
  i = 0;
  do {
    i++;
    if (LINK->thishet[i - 1]) {
      LINK->hap1[i - 1] = 1;
      LINK->hap2[i - 1] = 2;
    } else {
      LINK->hap1[i - 1] = 0;
      LINK->hap2[i - 1] = 0;
    }
  } while (!(LINK->thishet[i - 1] || i == mlocus));
  if (i == mlocus)
    setrec(0.5, LINK);
  else if (interfer)
    dointer(LINK);
  else
    nexthet(i, 0.5, true, LINK);
  if (LINK->nhap > 1 && !interfer)
    scramble(LINK);
  LINK->LINK->here = LINK->there;
}

Local Void gethet(system, LINK)
int *system;
struct LOC_recombine *LINK;
{
  int newsystem;

  newsystem = *system + 1;
  LINK->thishet[*system - 1] = false;
  if (*system != mlocus)
    gethet(&newsystem, LINK);
  else
    getrecprob(LINK);
  LINK->thishet[*system - 1] = true;
  if (*system != mlocus)
    gethet(&newsystem, LINK);
  else
    getrecprob(LINK);
}

Local Void recombine(theta_, segprob_, LINK)
double *theta_;
double *segprob_;
struct LOC_recombination *LINK;
{
  struct LOC_recombine V;
  int system;

  V.LINK = LINK;
  V.theta = theta_;
  V.segprob = segprob_;
  LINK->here = 0;
  system = 1;
  gethet(&system, &V);
}

Local Void getfemaletheta(LINK)
struct LOC_recombination *LINK;
{
  double dist;
  int ntheta, i;

  if (interfer)
    ntheta = mlocus;
  else
    ntheta = mlocus - 1;
  for (i = 0; i < ntheta; i++) {
    dist = getdist(&maletheta->theta[i]) * distratio;
    femaletheta->theta[i] = invdist(&dist);
  }
}


Void recombination()
{
  struct LOC_recombination V;
  int i;
  thetarray oldtheta;
  thetavalues *WITH;

  if (interfer) {
    WITH = maletheta;
    if (mapping)
      WITH->theta[mlocus - 1] = mapfunction(WITH->theta[mlocus - 2],
					    WITH->theta[0]);
    memcpy(oldtheta, WITH->theta, sizeof(thetarray));
    if (!mapping && !dolod) {
      for (i = 0; i < mlocus; i++)
	oldtheta[i] = 1 / (1 + exp(oldtheta[i]));
      WITH->theta[0] = oldtheta[0] + oldtheta[mlocus - 1];
      WITH->theta[mlocus - 2] = oldtheta[mlocus - 2] + oldtheta[mlocus - 1];
      WITH->theta[mlocus - 1] = oldtheta[0] + oldtheta[mlocus - 2];
      V.p1 = oldtheta[0];
      V.p2 = oldtheta[mlocus - 2];
      V.p3 = oldtheta[mlocus - 1];
      V.p4 = 1.0 - V.p1 - V.p2 - V.p3;
    } else {
      V.p1 = (WITH->theta[0] +
	      WITH->theta[mlocus - 1] - WITH->theta[mlocus - 2]) / 2.0;
      V.p2 = (WITH->theta[mlocus - 2] +
	      WITH->theta[mlocus - 1] - WITH->theta[0]) / 2.0;
      V.p3 = (WITH->theta[mlocus - 2] +
	      WITH->theta[0] - WITH->theta[mlocus - 1]) / 2.0;
      V.p4 = 1.0 - V.p1 - V.p2 - V.p3;
    }
    recombine(WITH->theta, WITH->segprob, &V);
  } else
    recombine(maletheta->theta, maletheta->segprob, &V);
  if (sexdif) {
    if (!readfemale) {
      if (interfer && !dolod) {
	WITH = maletheta;
	WITH->theta[0] = oldtheta[0] + oldtheta[mlocus - 1];
	WITH->theta[mlocus - 2] = oldtheta[mlocus - 2] + oldtheta[mlocus - 1];
	WITH->theta[mlocus - 1] = oldtheta[0] + oldtheta[mlocus - 2];
      }
      getfemaletheta(&V);
    }
    if (interfer) {
      WITH = femaletheta;
      if (mapping)
	WITH->theta[mlocus - 1] = mapfunction(WITH->theta[mlocus - 2],
					      WITH->theta[0]);
      if (readfemale && !mapping && !dolod) {
	memcpy(oldtheta, WITH->theta, sizeof(thetarray));
	for (i = 0; i < mlocus; i++)
	  oldtheta[i] = 1 / (1 + exp(oldtheta[i]));
	WITH->theta[0] = oldtheta[0] + oldtheta[mlocus - 1];
	WITH->theta[mlocus - 2] = oldtheta[mlocus - 2] + oldtheta[mlocus - 1];
	WITH->theta[mlocus - 1] = oldtheta[0] + oldtheta[mlocus - 2];
	V.p1 = oldtheta[0];
	V.p2 = oldtheta[mlocus - 2];
	V.p3 = oldtheta[mlocus - 1];
	V.p4 = 1.0 - V.p1 - V.p2 - V.p3;
      } else {
	V.p1 = (WITH->theta[0] +
		WITH->theta[mlocus - 1] - WITH->theta[mlocus - 2]) / 2.0;
	V.p2 = (WITH->theta[mlocus - 2] +
		WITH->theta[mlocus - 1] - WITH->theta[0]) / 2.0;
	V.p3 = (WITH->theta[mlocus - 2] +
		WITH->theta[0] - WITH->theta[mlocus - 1]) / 2.0;
	V.p4 = 1.0 - V.p1 - V.p2 - V.p3;
      }
      recombine(WITH->theta, WITH->segprob, &V);
    } else
      recombine(femaletheta->theta, femaletheta->segprob, &V);
  }
#if ((!defined DOS) && (defined(ILINK) || defined(LODSCORE)))
  if (firsttime || (checkpointStatus == normalRun)) {
#else
  if (firsttime) {
#endif
  }
#if PARALLEL
  /* Propagate thetas to slaves if in fun() */
  if(infun) {
    for(i=0;i<nuneed;i++) {
      gmalesegprob[currentthetanum][i] = maletheta->segprob[i];
      gfemalesegprob[currentthetanum][i] = femaletheta->segprob[i];
    }
    if((*slavesPerGroup) > 1) {
#if BARRIER_OUTPUT
      printBarrier("(%d) Group barrier at end of recombination\n",++barnum);
#endif /*BARRIER_OUTPUT*/
#if IS_SHMEM
      BARRIER(partial_barrier[mymaster]->barrier,*slavesPerGroup);
#else
      Tmk_barrier(MakeMask(mymaster,mymaster+(*slavesPerGroup)-1));
#endif /*IS_SHMEM*/
    }
  }
#endif /*PARALLEL*/
}




