/* This file contains some low level probability
   routines used in the MLINK program */

#include "commondefs.h"


/* Local variables for getvect: */
struct LOC_getvect {
  struct LOC_likelihood *LINK;
  thisperson *p;
  hapvector hap1, hap2;
} ;

/*static void getgene (long syste, double val, struct LOC_getvect *LINK);
static void ugetgene (long syste, double val, struct LOC_getvect *LINK);*/

static void getgene ();
static void ugetgene ();

Local double quanfun(phen, thislocus, i, j, mean, LINK)
phenotype *phen;
locusvalues *thislocus;
long i, j;
double *mean;
struct LOC_getvect *LINK;
{
  double val;
  long it, jt, FORLIM, FORLIM1;

  val = 1.0;
  if (phen->missing)
    return val;
  val = 0.0;
  FORLIM = thislocus->UU.U1.ntrait;
  for (it = 0; it < FORLIM; it++) {
    FORLIM1 = thislocus->UU.U1.ntrait;
    for (jt = 0; jt < FORLIM1; jt++) {
      if (i == j)
	val += thislocus->UU.U1.vmat[it]
	       [jt] * (phen->x[jt] - mean[jt]) *
	       (phen->x[it] - mean[it]);
      else
	val += thislocus->UU.U1.conmat * thislocus->UU.U1.vmat[it]
	       [jt] * (phen->x[jt] - mean[jt]) *
	       (phen->x[it] - mean[it]);
    }
  }
  val = thislocus->UU.U1.det * exp(-val * 0.5);
  if (i != j)
    val *= thislocus->UU.U1.contrait;
  return val;
}  /*quanfun*/


Local Void getval(syste, i, j, val, LINK)
long syste, i, j;
double *val;
struct LOC_getvect *LINK;
{
  locusvalues *WITH;
  phenotype *WITH1;

  WITH = thislocus[syste - 1];
  switch (WITH->which) {

  case quantitative:
    *val *= quanfun(LINK->p->phen[syste - 1], thislocus[syste - 1], i, j,
		    WITH->UU.U1.pm[i][j - 1], LINK);
    break;

  case affection:
    WITH1 = LINK->p->phen[syste - 1];
    *val *= WITH->UU.U0.pen[i][j - 1][WITH1->aff]
      [WITH1->liability - 1];
    break;
  /* cgh -- added this for gcc */
  case binary_:
  case null_:
  default:
    break;
  }
}  /*getval*/

#include "comgetvect.c"

