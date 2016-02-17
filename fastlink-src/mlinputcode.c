/* This file contains part of a modified version of the MLINK program */
/* The modifications are described in the papers: */
/* R. W. Cottingham, Jr., R. M. Idury, and A. A. Schaffer, */
/* Faster Sequential Genetic Linkage Computations */
/* American Journal of Human Genetics, 53(1993), pp. 252-263 */
/* and A. A. Schaffer, S. K. Gupta, K. Shriram, and R. W. Cottingham, Jr. */
/* Avoiding Recomputation in Genetic Linkage Analysis, */
/* Human Heredity 44(1994), PP. 225-237 */

#include "commondefs.h"
#include "gemdefs.h"
#include "mldefs.h"

Void readloci(LINK)
struct LOC_inputdata *LINK;
{
  struct LOC_readloci V;
  long i, j, coupling, autosomal, independent, difference;
  locusvalues *WITH;

  V.LINK = LINK;
  lastpriv = 0;
  fscanf(datafile, "%d%d%ld%*[^\n]", &mlocus, &risksys, &autosomal);
  getc(datafile);
  /*Replace the above line with the next when using epistasis*/
  /*readln(datafile,nlocus,risksys,autosomal,lastpriv);*/
  if (mlocus > maxlocus)
    inputerror(0L, mlocus, mlocus, LINK);
  if (mlocus <= 0)
    inputerror(1L, mlocus, mlocus, LINK);
  if (risksys > maxlocus)
    inputerror(37L, risksys, risksys, LINK);
  if (risksys < 0)
    inputerror(38L, risksys, risksys, LINK);
  risk = (risksys != 0);
  sexlink = (autosomal == 1);
#if PARALLEL
  if (Tmk_proc_id == 0) {
#endif
  printf("YOU ARE USING LINKAGE (V%3.1f (1-Feb-1991)) WITH%3d-POINT\n",
	 fVersion, mlocus);
  printf("YOU ARE USING FASTLINK (V%s (6-Oct-1997))",
	 fastversion);
  if (sexlink)
    printf(" SEXLINKED DATA\n");
  else
    printf(" AUTOSOMAL DATA\n");
#if PARALLEL
  }
#endif
  fscanf(datafile, "%d%lf%lf%ld%*[^\n]", &mutsys, &mutmale, &mutfemale,
	 &coupling);
  getc(datafile);
  if (mutsys > maxlocus)
    inputerror(39L, mutsys, mutsys, LINK);
  if (mutsys < 0)
    inputerror(40L, mutsys, mutsys, LINK);
  if (coupling == 1)
    disequi = true;
  else
    disequi = false;
  if (disequi) {
    hapfreq = (thisarray *)Malloc(sizeof(thisarray));
    if (!(hapfreq))
      malloc_err("hapfreq");
  }
  for (i = 1; i <= mlocus; i++) {
    fscanf(datafile, "%ld", &j);
    if (j > mlocus)
      inputerror(2L, i, j, LINK);
    if (j <= 0)
      inputerror(3L, i, j, LINK);
    order[j - 1] = i;
  }
  for (i = 1; i <= mlocus; i++) {
    for (j = 1; j < i; j++) {
      if (order[i - 1] == order[j - 1])
	inputerror(4L, i, j, LINK);
    }
  }
  fscanf(datafile, "%*[^\n]");
  getc(datafile);
  if (mutsys != 0)
    mutsys = order[mutsys - 1];
  if (risksys != 0)
    risksys = order[risksys - 1];
  V.nupriv = 0;
  for (i = 0; i < mlocus; i++)
    getlocus(order[i], &V);
  increment[mlocus - 1] = 1;
  for (i = mlocus - 1; i >= 1; i--)
    increment[i - 1] = increment[i] * thislocus[i]->nallele;
  fgeno = 1;
  for (j = 0; j < mlocus; j++)
    fgeno *= thislocus[j]->nallele;
  mgeno = fgeno;
  nuhap = fgeno;
  for (i = 0; i < mlocus; i++)
    nohom[i] = false;
  if (disequi) {
    allocate_thisarray(hapfreq, mgeno);
    for (i = 0; i < mgeno; i++)
      fscanf(datafile, "%lf", &hapfreq->genarray[i]);
    fscanf(datafile, "%*[^\n]");
    getc(datafile);
  } else {
    for (i = 0; i < mlocus; i++) {
      WITH = thislocus[i];
      if (WITH->which == affection || WITH->which == quantitative) {
	if (WITH->freq[affall - 1] < minfreq)
	  nohom[i] = true;
      }
    }
  }
  fgeno = fgeno * (fgeno + 1) / 2;
  if (!sexlink)
    mgeno = fgeno;
  fscanf(datafile, "%ld", &difference);
  if ((unsigned long)difference > 2) {
    inputwarning(0L, difference, difference, LINK);
    difference = 0;
  }
  sexdif = (difference != 0);
  readfemale = (difference == 2);
  fscanf(datafile, "%ld%*[^\n]", &independent);
  getc(datafile);
  if ((unsigned long)independent > 2) {
    inputwarning(1L, independent, independent, LINK);
    independent = 0;
  }
  interfer = (independent != 0);
  mapping = (independent == 2);
  gettheta(&maletheta, &V);
  if (sexdif)
    gettheta(&femaletheta, &V);
  else
    femaletheta = maletheta;
  if (sexlink && sexdif) {
    inputwarning(2L, difference, difference, LINK);
    sexdif = false;
    readfemale = false;
  }
  fscanf(datafile, "%ld%lf%lf%*[^\n]", &which, &inc, &finish);
  getc(datafile);
  finish += 0.0001;
  if (!sexlink) {
    if (mutsys == 0)
      thispath = auto_;
    else
      thispath = mauto;
  } else if (mutsys == 0)
    thispath = sex;
  else
    thispath = msex;
  segscale = scale;
  for (i = 1; i <= mlocus; i++)
    segscale *= scalemult;
}




