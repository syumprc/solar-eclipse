/* This file contains modified versions of the routines to do single
   likelihood computations for the LINKAGE package. The modifications herein
   are mostly needed to handle loops better and are described in the
   paper: 
   A. A. Schaffer, S. K. Gupta, K. Shriram, R. W. Cottingham, Jr. 
   Avoiding Recomputation in Linkage Analysis
   Human Heredity 44(1994), pp. 225-237. */


/*
   Modified by Dylan in late 1994 to speed up pedigrees with loops.
   #if LOOPSPEED were added at that time.
   Further modified by A. A. Schaffer in mid-1995 to handle loops
   in which loopbreakers have know genotypes. More uses of
   LOOPSPEED added at that time.
*/


#include "commondefs.h"
#ifndef LESSMEMORY
#include "moddefs.h"
#endif

#if defined(LODSCORE)
#include "lodefs.h"
#endif

#if PARALLEL
#if IS_SHMEM

void Tmk_errexit(format, va_alist)
char* format;
va_dcl
{
  va_list args;
  
  va_start(args);
  vfprintf(stderr, format, args);
  va_end(args);
  
  exit(-1);
}

#if defined(GMEM_DEBUG)
extern int totMem;
#endif /* defined(GMEM_DEBUG) */

char *Tmk_malloc(size)
unsigned size;
{
#if defined(GMEM_DEBUG)
  totMem += size;
#endif /* defined(GMEM_DEBUG) */
  return G_MALLOC(size);
}

#endif  /* IS_SHMEM */

unsigned MakeMask(start, end)
int start, end;
{
  unsigned u = 0;
  int i;

  for (i = start; i <= end; i++)
    u |= (1 << i);
  return u;
}

void OneMaster()
{
  iAmAMaster[0] = true;
  whoIsMyMaster[0] = 0;
  thetanumbase[0] = 0;
  thetanumfence[0] = 1;
  *slavesPerGroup = Tmk_nprocs;
  *checkmaster = 0;
  currentthetanum = 0;
  thischild = gMem->gthischild[mymaster];
  pnonzgens = gpnonzgens[currentthetanum];
  qnonzgens = gqnonzgens[currentthetanum];
}

short findemptygen()
{
   short i;
   for (i = 0; i <  maxworkingset; i++)
     if (ggeninuse[currentthetanum][i] == 0) 
       return(i);
   Tmk_errexit("Problem in findemptygen --- the estimate of %d for maxworkingset is too low.\n Specify a larger value with the -w option.\n If this estimate was computed automatically, \n then you have hit a bug that you should report. \n Consult README.p4 or\nREADME.TreadMarks for details.\n", maxworkingset);
   /* cgh -- this can't happen, but it makes gcc happy */
   return -1;
}

#endif /*PARALLEL*/

Void invert(m, n, det)
double (*m)[maxtrait];
int n;
double *det;
{
  covmatrix v;
  double val;
  int i, j, k;

  *det = 1.0;
  for (i = 1; i <= n; i++) {
    val = m[i - 1][i - 1];
    for (k = 0; k <= i - 2; k++)
      val -= v[k][i - 1] * v[k][i - 1];
    *det *= val;
    v[i - 1][i - 1] = sqrt(val);
    for (j = i; j < n; j++) {
      val = m[i - 1][j];
      for (k = 0; k <= i - 2; k++)
	val -= v[k][i - 1] * v[k][j];
      v[i - 1][j] = val / v[i - 1][i - 1];
      v[j][i - 1] = 0.0;
    }
  }
  for (i = 1; i <= n; i++) {
    m[i - 1][i - 1] = 1 / v[i - 1][i - 1];
    for (j = i + 1; j <= n; j++) {
      val = 0.0;
      for (k = i - 1; k <= j - 2; k++)
	val -= m[k][i - 1] * v[k][j - 1];
      m[j - 1][i - 1] = val / v[j - 1][j - 1];
      m[i - 1][j - 1] = 0.0;
    }
  }
  for (i = 0; i < n; i++) {
    for (j = 0; j < i + 1; j++) {
      val = 0.0;
      for (k = j; k < n; k++)
	val += m[k][i] * m[k][j];
      v[i][j] = val;
      v[j][i] = val;
    }
  }
  memcpy(m, v, sizeof(covmatrix));
}

/* static void collapsedown (thisperson *p, struct LOC_likelihood *LINK);*/
static void collapsedown ();
static void precollapsedown();

Void cleanup(p, LINK)
thisperson **p;
struct LOC_likelihood *LINK;
{
  thisperson *WITH;

  WITH = *p;
  if ((!(WITH->loopneeded)) || looplastgen) {   /*changed by A. A. Schaffer*/
#if !PARALLEL
    Free(WITH->gen->sparseflag);
    Free(WITH->gen->genarray);
    Free(WITH->gen);
#endif
    WITH->gen = NULL;
#if PARALLEL
    ggeninuse[currentthetanum][WITH->memindex] = 0; /*Used to free cell in 
                                           genbank for parallel computation*/
#endif
  }
  WITH->newgenexists = false;
}  /*cleanup*/

void allocate_thisarray(location, number)
thisarray *location;
int number;
{
   location->genarray = (double *) malloc(number * sizeof(double));
   if (NULL == location->genarray) 
     malloc_err("genarray field");
   location->sparseflag = (unsigned char *) malloc(number * 
                                                     sizeof(unsigned char));
   if (NULL == location->sparseflag) 
     malloc_err("sparseflag field");
}

Local Void prob(p, LINK)
thisperson **p;
struct LOC_likelihood *LINK;
{
  int i;
  thisperson *WITH;
  unsigned char *tempflag1; /*Next two declarations by R. M. Idury*/
  double *tempwith1; 
#if PARALLEL
  short genindex;  /*Index into genbank for parallel computation*/
#endif

  WITH = *p;
  if (WITH->newgenexists) return;  /*Fixed by A. A. Schaffer*/
  if (WITH->gen == NULL) {
#if !PARALLEL
    WITH->gen = (thisarray *) Malloc(sizeof(thisarray));
    if (WITH->gen == NULL)
      malloc_err("gen field in prob");
    allocate_thisarray(WITH->gen, fgeno);
#else
    /*Next four lines added by A.A. Schaffer for parallel computation*/
    genindex = findemptygen();
    WITH->gen = &(ggenbank[currentthetanum][genindex]);
    WITH->memindex = genindex;
    ggeninuse[currentthetanum][genindex] = 1;
#endif
    WITH->newgenexists = true; /*Inserted by A. A. Schaffer*/
    /*Code to initialze genarray and sparseflag put in by R. M. Idury*/
    tempwith1 = WITH->gen->genarray;
    tempflag1 = WITH->gen->sparseflag;
    for(i=0;i<fgeno;i++) {
      tempflag1[i] = 0;
      tempwith1[i] = 0.0;
    }
    if (WITH->inloop > 0) {
      if (looppers[LINK->thisped - 1][WITH->inloop - 1][0] == *p) { /*G OK*/
	tempwith1[LINK->loopgen[WITH->inloop - 1] - 1] = 
          LINK->holdpoint[WITH->inloop - 1]->
	    genarray[LINK->loopgen[WITH->inloop - 1] - 1];
	tempflag1[LINK->loopgen[WITH->inloop - 1] - 1] = 1;
      }
      else {
	tempwith1[LINK->loopgen[WITH->inloop - 1] - 1] = 1.0;
	tempflag1[LINK->loopgen[WITH->inloop - 1] - 1] = 1;
      }
      return;
    }
    getvect(*p, LINK);
  }
  if ((*p)->pa == NULL && (WITH->inloop==0))
    LINK->nuscale++;
}  /*prob*/


/*pseudoprob does just the scaling function of prob */
/*This routine written by A. A. Schaffer*/
Local Void pseudoprob(p, LINK)
thisperson **p;
struct LOC_likelihood *LINK;
{
  thisperson *WITH;

  WITH = *p;
  if (!(WITH->newgenexists))
    if ((*p)->pa == NULL) {
      LINK->nuscale++;
      if (WITH->inloop > 0) {
	LINK->nuscale--;
      }
    }
}  /*pseudoprob*/

Void initseg(LINK)
struct LOC_seg *LINK;
{
  if ((*LINK->p)->male) {
    LINK->nfirst = mgeno;
    LINK->nsecond = fgeno;
    LINK->firstsex = maletheta;
    LINK->secondsex = femaletheta;
    LINK->pf = mutmale;
    LINK->ps = mutfemale;
  } else {
    LINK->nfirst = fgeno;
    LINK->nsecond = mgeno;
    LINK->firstsex = femaletheta;
    LINK->secondsex = maletheta;
    LINK->pf = mutfemale;
    LINK->ps = mutmale;
  }
  prob(&LINK->father, LINK->LINK);
  prob(&LINK->mother, LINK->LINK);
  LINK->child = LINK->father->foff;
  nchild = 0;
  do {
    prob(&LINK->child, LINK->LINK);
    if (LINK->child->ma == LINK->mother && (!LINK->child->up)) {
      nchild++;
      childarray[nchild - 1] = LINK->child;  /*Line added by A. A. Schaffer*/
      thischild[nchild - 1] = LINK->child->gen;
      malechild[nchild - 1] = LINK->child->male;
    }
    LINK->child = LINK->child->nextpa;
  } while (LINK->child != NULL);  
#if PARALLEL
  gMem->nchild[mymaster] = nchild;
#endif
  if (nchild > maxchild) {
    fprintf(stderr, "\nA nuclear family has more children than maxchild");
    fprintf(stderr, "\nThe program will exit politely to allow you to");
    fprintf(stderr, "\nincrease maxchild in the file commondefs.h");
    exit(EXIT_FAILURE);
  }
} /*initseg*/


/*pseudoseg simulates the calls that initseg would make to prob by
using pseudoprob just to get scaling right*/
/*This routine written by A. A. Schaffer*/
Local Void pseudoseg(LINK)
struct LOC_seg *LINK;
{
  pseudoprob(&LINK->father, LINK->LINK);
  pseudoprob(&LINK->mother, LINK->LINK);
  LINK->child = LINK->father->foff;
  nchild = 0;
  do {
    pseudoprob(&LINK->child, LINK->LINK);
    if (LINK->child->ma == LINK->mother && (!LINK->child->up)) {
      nchild++;
    }
    LINK->child = LINK->child->nextpa;
  } while (LINK->child != NULL);   
  LINK->LINK->nuscale++;  /*does extra scaling skipped in exitseg*/
} /*pseudoseg*/


Void exitseg(LINK)
struct LOC_seg *LINK;
{
  LINK->LINK->nuscale++;
  LINK->child = LINK->father->foff;
  do {
    if (LINK->child->ma == LINK->mother && (!LINK->child->up))
      cleanup(&LINK->child, LINK->LINK);
    LINK->child = LINK->child->nextpa;
  } while (LINK->child != NULL);   /*exitseg*/
}

Local Void getgenetables()
{
  int locfstart, locfend, locfseg;
  int f1, f2, i, j, first;

  currentfence = 1;
  for (first = 0; first < fgeno; first++) {
    base[first] = currentfence;
    locfstart = probstart[first];
    locfend = probend[first];
    locfseg = segstart[first];
    for (i = locfstart-1; i < locfend; i++) {
      f1 = invgenenum1[locfseg - 1];
      haps1[currentfence] = f1;
      f2 = invgenenum2[locfseg - 1];
      haps2[currentfence] = f2;
      hind[currentfence++] = i;
      locfseg++;
    }
    fence[first] = currentfence;
  } 
}  /* getgenetables */

#if ALLELE_SPEED
/*The following procedure recomputes the number of haplotypes and
genotypes for a given pedigree. Written by A. A. Schaffer */
void recompute_haps(pedix, wantlocations)
int pedix;
boolean wantlocations;
{
  int temphap, tempfgeno, tempall;
  int locidx;
  locusvalues **templocus;
#if defined(LODSCORE)
  int locus1, locus2;
#endif

  templocus = thislocus;

  temphap = 1;
  tempfgeno = 1;
  /*need to adjust the increment array for lookups into genenumber*/
  increment[mlocus - 1] = 1;

#if defined(LODSCORE)
  locus1 = locuslist1[iplace - 1] - 1;
  locus2 = locuslist2[jplace - 1] - 1;
  if ((binary_ == templocus[1]->which) &&
      (allformat == templocus[1]->format)) {
    tempall = ped_new_allele_count[pedix - 1][locus2];
    templocus[1]->nallele = tempall;
  }
  else
    tempall = templocus[1]->nallele;
  increment[0] = increment[1] * tempall;
  temphap *= tempall;
  if ((binary_ == templocus[0]->which) &&
      (allformat == templocus[0]->format)) {
    tempall = ped_new_allele_count[pedix - 1][locus1];
    templocus[0]->nallele = tempall;
  }
  else
    tempall = templocus[0]->nallele;
  temphap *= tempall;
#else
  for(locidx = mlocus - 1; locidx >= 0; locidx--) {
    if ((binary_ == templocus[locidx]->which) &&
        (allformat == templocus[locidx]->format)) {
      tempall = ped_new_allele_count[pedix - 1][locidx];
      templocus[locidx]->nallele = tempall;
    }
    else
      tempall = templocus[locidx]->nallele;
    if (locidx > 0)
      increment[locidx - 1] = increment[locidx] * tempall;
    temphap *= tempall;
  }
#endif /*defined(LODSCORE) */

  nuhap = temphap;
  fgeno = nuhap * (nuhap + 1) / 2;
  if (sexlink)
    mgeno = nuhap;
  else
    mgeno = fgeno;
  if (ped_must_change_locations[pedix - 1] && wantlocations) {
    getlocations();
    getgenetables();
  }
}
#endif


/*computenumhaps computes the weighted sum over all genotypes of
how many haplotypes can be passed on to a child.
The following explanation assumes some basic understanding of
combinatorics at the undergraduate level. The total number of
haplotypes that can be passed on is
a weighted sum of the number of genotypes, where
different genotypes get different weights,  We first
classify the genotypes by heterozygosity pattern, where 0
means homozygous and 1 means heterozygous. There are mlocus loci
and therefore 2**mlocus different heterozygosity patterns.
The number of genotypes of a given heterozygosity pattern can be computed
by considering one locus at a time and multiplying the contributions
of all the loci, since they are independent.
A homozygous locus with A alleles contributes a factor of A.
The first heterozygous locus contributes a factor of A *(A-1)/2 if
it has A alleles.
Any other heterozygous locus with B alleles, contributes a factor of B*(B-1).

The weight of a genotype is 1 if it is homozygous; otherwise it
is 2**(h-1), where h is the number of heterozygous loci in the genotype.*/

static int computenumhaps()
{
  int numpatterns;
  hapvector *patternmatrix;
  int p; /*index over heterozygosity patterns*/
  int l; /*index over loci*/
  int value; /*accumulator for return value*/
  int term; /*term for this heterozygosity pattern*/
  boolean hetyet; /*have we seen a heterozygous locus yet*/ 
#if defined(LODSCORE)
  int locus1, locus2;
#endif

  numpatterns = 2;
  for (i=2; i<=mlocus;i++)
  numpatterns*=2;
  patternmatrix = (hapvector *) malloc(numpatterns*sizeof(hapvector));
  if (patternmatrix == NULL)
    malloc_err("patternmatrix");
  for(l=0; l < mlocus; l++)
    patternmatrix[0][l]=0;
  for(p=1; p< numpatterns; p++) {
    l = 0;
    do{
      patternmatrix[p][l] = !patternmatrix[p-1][l];
      l++;
    } while(!patternmatrix[p][l-1]);
    do{
      patternmatrix[p][l] = patternmatrix[p-1][l];
      l++;
    } while(l < mlocus);
  }

  value = 2; /*offset by 1 because 0 index has special meaning*/
  for(p=0; p < numpatterns; p++){
    term = 1;
    hetyet = 0;
    for(l=0; l < mlocus; l++){
      if (!patternmatrix[p][l]) {
        term *= thislocus[l]->nallele;
      }
      else 
        if (!hetyet) {
          hetyet = 1;
          term *= (thislocus[l]->nallele * (thislocus[l]->nallele -1))/2;
        }
        else 
          term *= thislocus[l]->nallele * (thislocus[l]->nallele - 1) * 2;
    }      
    value+=term;
  }
  free(patternmatrix);
  return(value);
}

Void allocategenetables()
{
  int maxfgeno, thisfinalfence, maxfinalfence;
  int ped;
  int i, c;

  maxfemgen = 0;
  maxhaplo = 0;
#if ALLELE_SPEED
  maxfgeno = 0;
  maxfinalfence = 0;
  for(ped = 1; ped <= mymaxped; ped++){
    recompute_haps(ped, false);
    if (fgeno > maxfgeno)
      maxfgeno = fgeno;
    if (nuhap > maxhaplo)
      maxhaplo = nuhap;
    thisfinalfence = computenumhaps();
    if (thisfinalfence > maxfinalfence)
     maxfinalfence = thisfinalfence;
  }
#else
  maxfgeno = fgeno;
  maxhaplo = nuhap;
  maxfinalfence = computenumhaps();
#endif /*ALLELE_SPEED*/

  maxfemgen = maxfgeno;

#ifndef LESSMEMORY
  indpool = (unsigned *) malloc(2 * maxclasssize * maxfgeno * sizeof (unsigned));
  if (indpool == NULL) {
    fprintf(stderr, "\nYou probably need to use the slower version of this program");
    malloc_err("indpool");
  }
  invpool = (unsigned *) malloc(2 * maxclasssize * maxfgeno * sizeof (unsigned));
  if (invpool == NULL) {
    fprintf(stderr, "\nYou probably need to use the slower version of this program");
    malloc_err("invpool");
  }
  nextpool = (unsigned *) malloc(2 * maxclasssize * maxfgeno * sizeof (unsigned));
  if (nextpool == NULL) {
    fprintf(stderr, "\nYou probably need to use the slower version of this program");
    malloc_err("nextpool");
  }
#endif

  genenumber = (int**) malloc(maxhaplo * sizeof(int*));
  if (genenumber == NULL)
    malloc_err("genenumber");
  for(i = 0; i < maxhaplo; i++) {
    genenumber[i] = (int *) malloc(maxhaplo * sizeof(int*));
    if (genenumber[i] == NULL)
      malloc_err("Entry in genenumber");
  }

#if (defined(ILINK) || defined (LODSCORE))
  if (approximate) {
    approxarray = (unschar **) malloc(nuped * sizeof(unschar *));
    if (approxarray == NULL)
      malloc_err("approxarray");
    for(i = 0; i <= nuped; i++) {
      approxarray[i] = (unschar *) malloc(maxfgeno * sizeof(unschar));
      if (approxarray[i] == NULL)
        malloc_err("Entry in approxarray");
    }
  }
#endif

  base = (unsigned *) malloc(maxfgeno * sizeof(unsigned));
  if (base == NULL)
    malloc_err("base");
  fence = (unsigned *) malloc(maxfgeno * sizeof(unsigned));
  if (fence == NULL)
    malloc_err("fence");
  invgenenum1 = (unsigned *) malloc(maxfgeno * sizeof(unsigned));
  if (invgenenum1 == NULL)
    malloc_err("invgenenum1");
  invgenenum2 = (unsigned *) malloc(maxfgeno * sizeof(unsigned));
  if (invgenenum2 == NULL)
    malloc_err("invgenenum2");
  segstart = (int *) malloc(maxfgeno * sizeof(int));
  if (segstart == NULL)
    malloc_err("segstart");
  probstart = (unsigned *) malloc(maxfgeno * sizeof(unsigned));
  if (probstart == NULL)
    malloc_err("probstart");
  probend = (unsigned *) malloc(maxfgeno * sizeof(unsigned));
  if (probend == NULL)
    malloc_err("probend");
  rare = (boolean *) malloc(maxfgeno * sizeof(unsigned));
  if (rare == NULL)
    malloc_err("rare");
  if (risk) {
    risk1 = (boolean *) malloc(maxfgeno * sizeof(boolean));
    if (risk1 == NULL)
      malloc_err("risk1");
    risk2 = (boolean *) malloc(maxfgeno * sizeof(boolean));
    if (risk2 == NULL)
      malloc_err("risk2");
    riskmale = (boolean *) malloc(maxhaplo * sizeof(boolean));
    if (riskmale == NULL)
      malloc_err("riskmale");
  }
  if (mutsys != 0) {
    muthap = (unschar *) malloc(maxhaplo * sizeof(unschar));
    if (muthap == NULL)
      malloc_err("muthap");
  }  

  haps1 = (unsigned short *) malloc(maxfinalfence * sizeof(unsigned short));
  if (haps1 == NULL)
    malloc_err("haps1");
  haps2 = (unsigned short *) malloc(maxfinalfence * sizeof(unsigned short));
  if (haps2 == NULL)
    malloc_err("haps2");
  hind = (unsigned int *) malloc(maxfinalfence * sizeof(unsigned int));
  if (hind == NULL)
    malloc_err("hind");
#if !(PARALLEL)  /*if parallel, used shared arrays*/
  nonzgens = (unsigned int *) malloc(maxfgeno * sizeof(unsigned int));
  if (nonzgens == NULL)
    malloc_err("nonzgens");
#endif
  gene = (double *) malloc(maxfgeno * sizeof(double));
  if (gene == NULL)
   malloc_err("gene");
  flag = (boolean *) malloc(maxfgeno * sizeof(boolean));
  if (flag == NULL)
    malloc_err("flag");
#if !(PARALLEL) /*in parallel, use arrpsumcache and arrqsumcache*/
  psumcache = (double *) malloc(maxhaplo * sizeof(double));
  if (psumcache == NULL)
    malloc_err("psumcache");
  qsumcache = (double *) malloc(maxhaplo * sizeof(double));
  if (qsumcache == NULL)
    malloc_err("qsumcache");
#else
  onechildupqsumcache = (double *) malloc(maxhaplo * sizeof(double));
  if (onechildupqsumcache == NULL)
    malloc_err("onechildupqsumcache");
#endif
#if !defined(LESSMEMORY)
  phapcache1 = (cache *) malloc(maxhaplo * sizeof(cache));
  if (phapcache1 == NULL)
    malloc_err("phapcache1");
#endif
#if PARALLEL
  stripe_pnonzgens = (unsigned int *) malloc(maxfgeno * sizeof(unsigned int));
  if (stripe_pnonzgens == NULL)
    malloc_err("stripe_pnonzgens");
  stripe_qnonzgens = (unsigned int *) malloc(maxfgeno * sizeof(unsigned int));
  if (stripe_qnonzgens == NULL)
    malloc_err("stripe_qnonzgens");
  privatepnonzgens = (unsigned int *) malloc(maxfgeno * sizeof(unsigned int));
  if (privatepnonzgens == NULL)
    malloc_err("privatepnonzgens");
  privateqnonzgens = (unsigned int *) malloc(maxfgeno * sizeof(unsigned int));
  if (privateqnonzgens == NULL)
    malloc_err("privateqnonzgens");
#endif
#if !defined(LESSMEMORY)
#if !PARALLEL
  if (!sexlink) {
    partialprob = (childprob **) malloc(maxclasssize * sizeof(childprob *));
    if (partialprob == NULL) {
      fprintf(stderr, "\nYou probably need to use the slower version of this program");
      malloc_err("partialprob");
    }
    for(c = 0; c < maxclasssize; c++) {
      partialprob[c] = (childprob *) malloc(maxhaplo * sizeof(childprob));
      if (partialprob[c] == NULL) {
	fprintf(stderr, "\nYou probably need to use the slower version of this program");
	malloc_err("Entry in partialprob");
      }
    }
  }
#endif
#endif

}  /* allocgenetables */

Void freegenetables()
{
  int i, c;

#ifndef LESSMEMORY
   free(indpool);
   indpool = NULL;

   free(invpool);
   invpool = NULL;

   free(nextpool);
   nextpool = NULL;
#endif

  for(i = 0; i < maxhaplo; i++) {
    free(genenumber[i]);
    genenumber[i] = NULL;
  }
  free(genenumber);
  genenumber = NULL;

#if (defined(ILINK) || defined (LODSCORE))
  if (approximate) {
    for(i = 0; i <= nuped; i++) {
      free(approxarray[i]);
      approxarray[i] = NULL;
    }
    free(approxarray);
    approxarray = NULL;
  }
#endif

  free(base);
  base = NULL;
  free(fence);
  fence = NULL;
  free(invgenenum1);
  invgenenum1 = NULL;
  free(invgenenum2);
  invgenenum2 = NULL;
  free(segstart);
  segstart = NULL;
  free(probstart);
  probstart = NULL;
  free(probend);
  probend = NULL;
  free(rare);
  rare = NULL;
  if (risk) {
    free(risk1);
    risk1 = NULL;
    free(risk2);
    risk2 = NULL;
    free(riskmale);
    riskmale = NULL;
  }
  if (mutsys != 0) {
    free(muthap);
    muthap = NULL;
  }  

  free(haps1);
  haps1 = NULL;
  free(haps2);
  haps2 = NULL;
  free(hind);
  hind = NULL;

#if !(PARALLEL)  /*if parallel, used shared arrays*/
  free(nonzgens);
  nonzgens = NULL;
#endif
  free(gene);
  gene = NULL;
  free(flag);
  flag = NULL;

#if !(PARALLEL) /*in parallel, use arrpsumcache and arrqsumcache*/
  free(psumcache);
  psumcache = NULL;
  free(qsumcache);
  qsumcache = NULL;
#endif
#if !defined(LESSMEMORY)
  free(phapcache1);
  phapcache1 = NULL;
#endif
#if !defined(LESSMEMORY)
#if !PARALLEL
  if (!sexlink) {
    for(c = 0; c < maxclasssize; c++) {
      free(partialprob[c]);
      partialprob[c] = NULL;
    }
    free(partialprob);
    partialprob = NULL;
  }
#endif
#endif
}  /* freegenetables */


/* getgeneindices is a preprocessing routine that compute the encodings
   of all the haplotypes and genotypes and precomputes for each
   genotype a list of the haplotypes it can pass on in one array
   insted of using indirect addressing by genotype */
Void getgeneindices()
{
  int locfstart, locfend, locfseg;
  int f1, f2, i, j, first;
  int finalfence;

#if defined(LODSCORE)  
  /*find number of recombination probabilities*/  
  nuneed = 7;
  for (i= 3; i<= mlocus; i++)
    nuneed = 5 * nuneed - 3;
#endif  /* defined(LODSCORE) */

#if ((!defined(LESSMEMORY)) && ((!PARALLEL) || (PRECOMPUTE)))
  segprob2 = (double *) malloc(nuneed * nuneed * sizeof(double));
  if (segprob2 == NULL) {
    fprintf(stderr, "\nYou probably need to use the slower version of this program");
    malloc_err("segprob2");
  }
#endif

  /*find size of isozygote class and it's square for joint classes*/
  maxclasssize = 2;
  for (i = 3; i<= mlocus; i++)
    maxclasssize *= 2;
  maxisozygclass = maxclasssize * maxclasssize;

  nuprobclass = 2;
  for(i = 2; i <= mlocus; i++)
    nuprobclass = 3 * nuprobclass - 1;

  segval = (double *) malloc(maxisozygclass * sizeof(double));
  if (segval == NULL)
    malloc_err("segval");
  tempseg = (double *) malloc(maxchild *maxisozygclass * sizeof(double));
  if (tempseg == NULL)
    malloc_err("tempseg");
  segindex = (unsigned *) malloc(4 * maxisozygclass * sizeof(unsigned));
  if (segindex == NULL)
    malloc_err("segindex");
#if (defined(LESSMEMORY) || (defined(PARALLEL) && (0 == PRECOMPUTE)))
  tempseg2 = (double *) malloc(maxchild *maxisozygclass * sizeof(double));
  if (tempseg2 == NULL)
    malloc_err("tempseg2");
#endif

#if ALLELE_SPEED
  allocategenetables();
#endif
#if !ALLELE_SPEED
  getgenetables();
#endif
}  /* getgeneindices */

/*Routine determines in which way to do the next nuclear family genarray
update. Modified by R. M. Idury to distinguish between old versions
of segup (called oldsegup) and new versions of segdown and segup.
Modified by A. A. Schaffer to avoid doing the update in
certain loop traversal pedigree situations */
Void seg(p_, q_, r_, peel, LINK)
thisperson **p_, **q_, **r_;
direction peel;
struct LOC_likelihood *LINK;
{
  struct LOC_seg V;
  boolean phaseunkn;


  V.LINK = LINK;
  V.p = p_;
  V.q = q_;
  V.r = r_;
  phaseunkn = ((*V.q)->pa == NULL && (*V.q)->firstpass &&
	       (*V.q)->inloop == 0 && !disequi);
  if ((*V.p)->male) {
    V.father = *V.p;
    V.mother = *V.q;
  } else {
    V.father = *V.q;
    V.mother = *V.p;
  }
  if (peel == peelup) {
    if (((*V.p)->pa == NULL) && phaseunkn && approximate) {
      if (firstapprox && firsttime) {
	switch (thispath) {

	case auto_:
	  segtop(&V);
	  break;

	case mauto:
	  segtop(&V);
	  break;

	case sex:
          segsextop(&V);
	  break;

	case msex:
	  segsextop(&V);
	  break;
	}
      } else {  /*first approximate not first time*/
	if (firstapprox) {
	  switch (thispath) {

	  case auto_:
	    segctop(&V);
	    break;

	  case mauto:
	    segctop(&V);
	    break;

	  case sex:
	    segsexctop(&V);
	    break;

	  case msex:
	    segsexctop(&V);
	    break;
	  }
	} else {  /*approximate*/
	  switch (thispath) {

	  case auto_:
	    segcapprox(&V);
	    break;

	  case mauto:
	    segcapprox(&V);
	    break;

	  case sex:
	    segsexctop(&V);
	    break;

	  case msex:
	    segsexctop(&V);
	    break;
	  }
	}
      }
    } else {  /*do not approximate*/
      if (phaseunkn) {
	if (firsttime) {
	  switch (thispath) {

	  case auto_:
            if (((*V.p)->loopdepend) || loopfirstgen)  
	      segup(&V);
	    else
	      pseudoseg(&V);
	    break;

	  case mauto:
	    segtop(&V);
	    break;

	  case sex:
            if (((*V.p)->loopdepend) || loopfirstgen)
	      segsexup(&V);
            else
              pseudoseg(&V);
	    break;

	  case msex:
	    segsextop(&V);
	    break;
	  }
	} else {  /*not firsttime*/
	  switch (thispath) {

	  case auto_:
            if (((*V.p)->loopdepend) || loopfirstgen)  
	      segup(&V);
	    else
	      pseudoseg(&V);
	    break;

	  case mauto:
	    segctop(&V);
	    break;

	  case sex:
            if (((*V.p)->loopdepend) || loopfirstgen)
	      segsexup(&V);
            else
              pseudoseg(&V);
	    break;

	  case msex:
	    segsexctop(&V);
	    break;
	  }
	}
      } else {  /*phaseinfo*/
	switch (thispath) {

	case auto_:
          if (((*V.p)->loopdepend) || loopfirstgen)  
	    segup(&V);
          else
            pseudoseg(&V);
	  break;

	case mauto:
	  oldsegup(&V);
	  break;

	case sex:
          if (((*V.p)->loopdepend) || loopfirstgen)
	    segsexup(&V);
          else
            pseudoseg(&V);
	  break;

	case msex:
	  oldsegsexup(&V);
	  break;
	}
      }
    }
  } else {  /*not peelup*/
    switch (thispath) {

    case auto_:
      if (((*V.r)->loopdepend) || loopfirstgen)
        segdown(&V);
      else
        pseudoseg(&V);
      break;

    case mauto:
      msegdown(&V);
      break;

    case sex:
      if (((*V.r)->loopdepend) || loopfirstgen)
        segsexdown(&V);
      else
        pseudoseg(&V);
      break;

    case msex:
      msegsexdown(&V);
      break;
    }
  }
  (*V.q)->firstpass = false;
  (*V.p)->firstpass = false;
}  /*seg*/

Local Void collapseup(p, LINK)
thisperson *p;
struct LOC_likelihood *LINK;
{
  thisperson *q, *child, *nextchild;
  boolean down;
  p->done = true;
  if (p->foff == NULL)
    return;
  down = false;
  child = p->foff;
  while (child != NULL) {
    down = false;
    if (p->male)
      q = child->ma;
    else
      q = child->pa;
    if (!q->done) {
      collapsedown(q, LINK);
      nextchild = child;
      while (nextchild != NULL) {
	if (nextchild->pa == q || nextchild->ma == q) {
	  if (!nextchild->up)
	    collapseup(nextchild, LINK);
	  else
	    down = true;
	}
	if (p->male)
	  nextchild = nextchild->nextpa;
	else
	  nextchild = nextchild->nextma;
      }
      if (q->multi)
	collapseup(q, LINK);
      if (!down)
	seg(&p, &q, &child, peelup, LINK);
      else
	collapsedown(p, LINK);
    }
    if (p->male)
      child = child->nextpa;
    else
      child = child->nextma;
  }
}  /*collapseup*/

Local Void precollapseup(p)
thisperson *p;
{
  thisperson *q, *child, *nextchild, *childonlist;
  boolean down;
  p->done = true;
  if (p->foff == NULL)
    return;
  down = false;
  child = p->foff;
  while (child != NULL) {
    down = false;
    if (p->male)
      q = child->ma;
    else
      q = child->pa;
    if (!q->done) {
      precollapsedown(q);
      nextchild = child;
      while (nextchild != NULL) {
	if (nextchild->pa == q || nextchild->ma == q) {
	  if (!nextchild->up)
	    precollapseup(nextchild);
	  else
	    down = true;
	}
	if (p->male)
	  nextchild = nextchild->nextpa;
	else
	  nextchild = nextchild->nextma;
      }
      if (q->multi)
	precollapseup(q);
      if (!down) {
        if (q->loopdepend) {
	  p->loopdepend = true;
	  p->loopneeded = false;
        }         
        childonlist = p->foff;
        while (childonlist != NULL) {
          if (childonlist->loopdepend) {
            p->loopdepend = true;
            p->loopneeded = false;
          }
	  if (p->male)
	    childonlist = childonlist->nextpa;
	  else
	    childonlist = childonlist->nextma; 
	}
      }
      else
	precollapsedown(p);
    }
    if (p->male)
      child = child->nextpa;
    else
      child = child->nextma;
  }
}  /*precollapseup*/



Local Void collapsedown(p, LINK)
thisperson *p;
struct LOC_likelihood *LINK;
{
  if (p->pa == NULL)
    return;
  p->up = true;
  collapseup(p->pa, LINK);
  seg(&p->pa, &p->ma, &p, peeldown, LINK);
}  /*collapsedown*/


Local Void precollapsedown(p, LINK)
thisperson *p;
struct LOC_likelihood *LINK;
{
  thisperson *childonlist;

  if (p->pa == NULL)
    return;
  p->up = true;
  precollapseup(p->pa);
  if (p->pa->loopdepend || p->ma->loopdepend) {
    p->loopdepend = true;
    p->loopneeded = false;
  }
  childonlist = p->pa->foff;
  while (childonlist != NULL) {
    if (childonlist->loopdepend) {
      p->loopdepend = true;
      p->loopneeded = false;
    }
    if (p->male)
      childonlist = childonlist->nextpa;
    else
      childonlist = childonlist->nextma; 
  }
}  /*precollapsedown*/

Local Void riskcumul(LINK)
struct LOC_likelihood *LINK;
{
  int i;
  thisarray *WITH;

  WITH = LINK->proband->gen;
  if (sexlink && LINK->proband->male) {
    for (i = 0; i < mgeno; i++) {
      if (riskmale[i])
	LINK->hetero += WITH->genarray[i];
    }
    return;
  }
  for (i = 0; i < fgeno; i++) {
    if (risk2[i])
      LINK->homo += WITH->genarray[i];
    else if (risk1[i])
      LINK->hetero += WITH->genarray[i];
  }
}  /*riskcumul*/


Local Void riskcalc(LINK)
struct LOC_likelihood *LINK;
{
  double normal;

  LINK->homo /= like;
  LINK->hetero /= like;
  normal = 1 - LINK->homo - LINK->hetero;
  fprintf(outfile, "RISK FOR PERSON %6d IN PEDIGREE %7d\n",
	  LINK->proband->id, LINK->proband->ped);
  if (!LINK->proband->male || !sexlink)
    fprintf(outfile, "HOMOZYGOTE CARRIER   : %8.5f\n", LINK->homo);
  if (!LINK->proband->male || !sexlink)
    fprintf(outfile, "HETEROZYGOTE CARRIER : %8.5f\n", LINK->hetero);
  else
    fprintf(outfile, "MALE CARRIER         : %8.5f\n", LINK->hetero);
  fprintf(outfile, "NORMAL               : %8.5f\n", normal);
  printf("RISK FOR PERSON %6d IN PEDIGREE %7d\n",
	 LINK->proband->id, LINK->proband->ped);
  if (!LINK->proband->male || !sexlink)
    printf("HOMOZYGOTE CARRIER   : %8.5f\n", LINK->homo);
  if (!LINK->proband->male || !sexlink)
    printf("HETEROZYGOTE CARRIER : %8.5f\n", LINK->hetero);
  else
    printf("MALE CARRIER         : %8.5f\n", LINK->hetero);
  printf("NORMAL               : %8.5f\n", normal);
}  /*riskcalc*/

/*pollutedescendants finds which people descended from startper by
depth first search and ensures that they will be reevaluated for
each genotype of the loopbreaker */
/*This routine written by A. A. Schaffer*/

Local Void pollutedescendants(startper)
thisperson *startper;
{
  thisperson *nextchild;

  startper->loopdepend=true;
  nextchild = startper->foff;
  while(nextchild != NULL) {
    if (nextchild->loopdepend == false) 
      pollutedescendants(nextchild);
    nextchild = nextchild->nextpa;
  }
  nextchild = startper->foff;
  while(nextchild != NULL) {
    if (nextchild->loopdepend == false) 
      pollutedescendants(nextchild);
    nextchild = nextchild->nextma;
  }
}

#if LOOPSPEED
/*
   This procedure takes the loopbreaker, joint-genotype array passed in
   'loopgen' for the pedigree passed in 'ped'.  It initializes the
   global array 'single_locus_vector_num' to hold, for each locus,
   the single-locus, loopbreaker vector number corresponding to 'loopgen'.
   (The loopbreaker vector numbers are indexes into 'loop_vectors'.)

   Written by Dylan in late 1994.  
*/
static void translate_loop_vector(loopgen, ped) 
int loopgen[maxloop];
int ped;
{
  int locus;  /*loop index for loci*/
  int loop;   /*index on loops*/
  int num_alleles; /*number of alleles for a locus*/
  int left, mid, right;  /* for binary search */
  boolean found; /*used in binary search*/
  int breaker; /*loop index for loop_breakers*/
  int geno; /*genotype of loopbreaker*/
  unsigned hap1, hap2; /*haplotypes of loop breaker*/
  unsigned a, b, temp; /*hold haplotypes for swapping*/
  int duplicates; /*duplicate genotypes due to heterozygous symmetry*/
  /* holds at each locus the single locus genotype of each loop breaker */
  int single_locus_geno[maxlocus][maxloop];

  /* if no loops trivial */
  if ( num_loops[ped] == 0 ) {
    for (locus = 0; locus < mlocus; locus++ ) {
      single_locus_vector_num[locus] = 0;
    }
    return;
  }

  /* for each loop breaker */
  for (loop = 0; loop < num_loops[ped]; loop++) {
    /* get the loop breaker joint genotype and joint haplotypes */
    geno = loopgen[loop] - 1;  /*need to decrement since loopgen starts at 1*/
    if (sexlink && (looppers[ped][loop][0]->male)) {
      hap1 = geno;             /*already decremented*/
      hap2 = geno;             
    }
    else {
      hap1 = invgenenum1[geno] - 1;
      hap2 = invgenenum2[geno] - 1;
    }

  
    /* for each locus in reverse order */
    for (locus = (mlocus - 1); locus >= 0; locus--) {
      num_alleles = thislocus[locus]->nallele;
      /* 
	 Calculate the left and right allele for this locus.

	 I'm taking advantage of knowing how the joint haplotypes
	 are computed in setval().
       */
      a = hap1 % num_alleles;
      hap1 = (hap1 - a) / num_alleles;
      b = hap2 % num_alleles;
      hap2 = (hap2 - b) / num_alleles;
      
      /* 
	 Calculate the the single loucs genotype, given the alleles.

	 I'm taking advantage of knowing how the genotypes are
	 calculated in unknown.
      */
      if ( a > b ) {
	temp = a;
	a = b;
	b = temp;
      }
      duplicates = a * (a + 1) / 2;

      /* store the single locus genotype */
      if (sexlink && (looppers[ped][loop][0]->male))
	single_locus_geno[locus][loop] = a;
      else
	single_locus_geno[locus][loop] = (a * num_alleles) - duplicates + b;
    }
  }


  /* for each locus */
  for  (locus = 0; locus < mlocus; locus++) {

    /* 
       binary search in 'loop_vectors[locus]' 
       to find 'single_locus_geno[locus]' 
    */
    left = 0;
    right = num_loop_vectors[locus] - 1;
    mid = (int) right / 2;
    found = false;
     
    do {
      for (breaker = num_loops[ped] - 1; breaker >= 0; breaker--) {
	if ( single_locus_geno[locus][breaker] == 
	    loop_vectors[locus][mid][breaker]) {
	  if ( breaker == 0 ) {
	    found = true;
	  }
	} else 
	  if (single_locus_geno[locus][breaker] < 
	      loop_vectors[locus][mid][breaker]) {
	    right = mid;
	    mid = (int) (right - left) / 2 + left;
	    break;
	  } else {
	      left = mid;
	      mid = (int) (right - left + 1) / 2 + left;
	      break;
	    }
      } /* for breaker */
      if (found == true) {
	break;
      }
    } while (left < right);
    
    if (found == true) {
      single_locus_vector_num[locus] = mid;
    } else {
      printf("Error in translate loop vector.\n"); /* DYLAN error */
      exit(EXIT_FAILURE);
    }

  } /* for each locus */

} /* translate_loop_vector */
#endif

/*identify possible genotypes for one loop breaker and store in the two
  dimensional array loopbreaker_nextgeno, i is index of loop breakers
  and fill in loopmax[i]*/
#if LOOPSPEED
Void find_loopbreaker_genotypes(ped,i, loopmax)
int ped;
int i;
int *loopmax;
{
  int geno, nextgeno; /*loop indices*/
  int temploopmax; /*placeholder for value of loopmax[i]*/
  thisarray *WITH1; /*placeholder for loop breaker genarray*/

  if (i < num_loops[ped -1]) {
    WITH1 = looppers[ped - 1][i][0]->gen;
    for(geno = 0; geno < loopmax[i]; geno++)
      if(!(breaker_poss_genotype[i][geno])) {
	WITH1->genarray[geno] = 0.0;
	WITH1->sparseflag[geno] = false;
      }
  }
  

  /*Initialize sparse array of possible genotypes for loopbreakers*/
    loopbreaker_nextgeno[i] = (int *) malloc(loopmax[i] * sizeof(int));
    if (loopbreaker_nextgeno[i] == NULL)
      malloc_err("loopbreaker_nextgeno entry");
    WITH1 = looppers[ped - 1][i][0]->gen;
    for(geno = 0, nextgeno=0; geno < loopmax[i]; geno++) {
      while ((nextgeno < loopmax[i]) &&  (WITH1->genarray[nextgeno] == 0.0))
	nextgeno++;
      loopbreaker_nextgeno[i][geno] = nextgeno;
      if (nextgeno < loopmax[i])
	temploopmax = nextgeno + 1;
      geno = nextgeno;
      nextgeno++;
    }
    loopmax[i] = temploopmax;
}
#endif /*LOOP_SPEED*/

/*
   likelihood computes the likelihood for thisped_ with proband_ as
   root of the computation.
   Modified by Dylan in late 1994 to use 'unknown_poss'.
   Modified further in mid 1995 by A.A. Schaffer to generate
   possible loop-breaker vectors more effectively
*/
Void likelihood(thisped_, proband_)
int thisped_;
thisperson *proband_;
{
#if !defined(LESSMEMORY)
  int ti, tj; /* put in by R. M. Idury */
#endif  /* !defined(LESSMEMORY) */
  struct LOC_likelihood V;
  int loopmax[maxloop];
  double tmplike;
  int i, j;
  boolean gocalc, alldone;
  thisperson *WITH;
  thisarray *WITH1;
  int lastgen; /*last genotype that first loop breaker can have*/

#if LOOPSPEED
  int locus; /*loop index on loci*/
  boolean is_zero; /*does this loop-breaker vector give a 0.0 likelihood?*/
  int geno, nextgeno; /*loop index, A. A. Schaffer */
  int temploopmax; /*placeholder for maximum genotype of one loopbreaker*/
  int nuscale_count; /*Count how many times we get a non-zero loop-breaker
                       so that scaling will be scale one time through
                       (per_iter_scale) * nuscale_count */
  int per_iter_scale;

#if LOOP_BREAKERS
  int copyIndex;
#endif /*LOOP_BREAKERS*/
#endif /*LOOPSPEED*/

  /* cgh -- added this for gcc warning */
  lastgen = 0;
  
  V.thisped = thisped_;
  V.proband = proband_;
#if ALLELE_SPEED
#if PARALLEL
  sharedCurrentPed[Tmk_proc_id] = thisped_;
#endif /*PARALLEL*/
  recompute_haps(thisped_, true);
  currentped = thisped_;
#endif /*ALLELE_SPEED*/
  if (!informative[V.thisped - 1]) {
    like = 0.0;
    return;
  }
  V.homo = 0.0;
  V.hetero = 0.0;
  tmplike = 0.0;
  alldone = false;
  V.nuscale = 0;
  /* Next two pieces of code added by A. A. Schaffer*/
  /* initialize loop structure; all descendants of place where
     first loop is broken will need to be update for each genotype*/
  for (i = 1; i <= totperson; i++) {
    person[i]->loopdepend = false;
    person[i]->loopneeded = false;
    person[i]->gen = NULL;
  }
  if (looppers[V.thisped-1][0][0] != NULL){ /*G OK*/
#if LOOP_BREAKERS
    for(copyIndex = 0; copyIndex < numCopies[V.thisped-1][0]; copyIndex++)
      pollutedescendants(looppers[V.thisped-1][0][copyIndex]); /*G OK*/
#else
    pollutedescendants(looppers[V.thisped-1][0][0]); /*G OK */
    pollutedescendants(looppers[V.thisped-1][0][1]);  /*G OK*/
#endif /*LOOP_BREAKERS*/
    for (i = 1; i <= totperson; i++) {
      WITH = person[i];
      WITH->newgenexists = false;
      WITH->done = false;
      WITH->up = false;
    }
    precollapseup(V.proband);
    precollapsedown(V.proband);
  }

#if LOOPSPEED
  /*initialization, A.A. Schaffer*/
  nuscale_count =0;
  per_iter_scale = 0;


  infer_genotypes(V.thisped);
#endif
  for (i = 0; i < maxloop; i++) {
    V.loopgen[i] = 1;
    loopmax[i] = 1;
    V.holdpoint[i] = NULL;
    if (looppers[V.thisped - 1][i][0] != NULL) { /*G OK*/
      WITH = looppers[V.thisped - 1][i][0];      /*G OK*/
      WITH->gen = (thisarray *)Malloc(sizeof(thisarray));
      if (WITH->gen == NULL)
	malloc_err("gen field in likelihood");
      allocate_thisarray(WITH->gen, fgeno);
      WITH->newgenexists = true; /*Added by A. A. Schaffer*/
      WITH1 = WITH->gen;
      for (j = 0; j < fgeno; j++)
	WITH1->genarray[j] = 0.0;
      getvect(looppers[V.thisped - 1][i][0], &V);  /*G OK*/
      if (looppers[V.thisped - 1][i][0]->pa == NULL) /*G OK*/
	V.nuscale++;
      V.holdpoint[i] = WITH->gen;
      if (WITH->male)
	loopmax[i] = mgeno;
      else
	loopmax[i] = fgeno;
#if LOOPSPEED
      find_loopbreaker_genotypes(V.thisped, i, loopmax);
      V.loopgen[i] = loopbreaker_nextgeno[i][0] + 1;
#endif
      WITH->gen = NULL;
      WITH->newgenexists = false;  /*Added by A. A. Schaffer*/
    }
  }
  V.loopgen[0] = 0;
  /*find last nonzero genotype in first loop genarray*/
  /*Added by A. A. Schaffer*/
  if (looppers[V.thisped-1][0][0] != NULL){  /*G OK*/
    for (j = loopmax[0]-1; j>=0; j--){
      if (V.holdpoint[0]->genarray[j] > 0.0){
        lastgen =j;
        break;
      }
    }
  }
  /*Next line added by A. A. Schaffer*/
  loopfirstgen = 1;
  do {
    i = 1;
    /*Next line added by A. A. Schaffer*/
    looplastgen = 0;
    do {
      V.loopgen[i - 1]++;
#if LOOPSPEED
        /*If we have a loop get the next possible genotype for first
          loop-breaker or wrap-arraound to its first possible genotype*/
        if (looppers[V.thisped-1][0][0] != NULL){  /*G OK*/
          if (V.loopgen[i - 1] > loopmax[i - 1])
            V.loopgen[i - 1] = loopbreaker_nextgeno[i-1][0] + 1;
          else {
            V.loopgen[i - 1] = loopbreaker_nextgeno[i -1][V.loopgen[i-1] - 1]  + 1;
            i = maxloop;
	  }
        /*Added by A. A. Schaffer*/
        if (i == 1)
          loopfirstgen = 1;
	}
      else {
        loopfirstgen = 1;
        i = maxloop;
      }
#else /*!LOOPSPEED*/
      if (V.loopgen[i - 1] > loopmax[i - 1]) {
	V.loopgen[i - 1] = 1;
        /*Added by A. A. Schaffer*/
        if (i == 1)
          loopfirstgen = 1;
      }
      else
	i = maxloop;
#endif /*LOOPSPEED*/
      i++;
    } while (i <= maxloop);
    gocalc = true;
    for (i = 0; i < maxloop; i++) {
      /*ML change*/
      if (V.holdpoint[i] != NULL) {
	if (V.holdpoint[i]->genarray[V.loopgen[i] - 1] == 0.0)
	  gocalc = false;
        else               /*else clause added by A. A. Schaffer*/
          if ((i == 0) && ((V.loopgen[i] - 1) == lastgen))
            looplastgen = 1;
      }
    }
    if (gocalc) {
#if LOOPSPEED
        nuscale_count++;
        /* get per locus loop vectors from joint genotype loop vector */
        translate_loop_vector(V.loopgen, V.thisped - 1);
        is_zero = false;
        for(locus = 0; locus < mlocus; locus++)
          if (is_zero_breaker[locus]
             [single_locus_vector_num[locus]])
            is_zero = true;
        if (is_zero)
          like = 0.0;
        else {
#endif
      for (i = 1; i <= totperson; i++) {
	WITH = person[i];
        /*next four lines written by A. A. Schaffer*/
        if (loopfirstgen || !(WITH->loopneeded)) {
         if (WITH->gen) 
#if !PARALLEL
	   {
           free(WITH->gen->sparseflag);
           free(WITH->gen->genarray);
           free(WITH->gen);
          }
#else
          ggeninuse[currentthetanum][WITH->memindex] = 0; /*Used to free cell  
                                          in genbank for parallel computation*/
#endif
          WITH->gen = NULL;
        }
        WITH->newgenexists = false;
	WITH->done = false;
	WITH->up = false;
      }

/* accumulate segprob entries */
/* next test written by R. M. Idury*/
#if ((!defined(LESSMEMORY)) && ((!PARALLEL) || (PRECOMPUTE)))
      for(ti = 0; ti < nuneed; ti++)
	for(tj = 0; tj < nuneed; tj++)
	  segprob2[ti*nuneed+tj] = maletheta->segprob[ti] * femaletheta->segprob[tj];
      getprobtable();
#endif
      collapseup(V.proband, &V);
      collapsedown(V.proband, &V);
      like = 0.0;
      WITH1 = V.proband->gen;
      if (V.proband->male) {
	for (i = 0; i < mgeno; i++)
	  like += WITH1->genarray[i];
      } else {
	for (i = 0; i < fgeno; i++)
	  like += WITH1->genarray[i];
      }

      /*Added by Alex*/
      loopfirstgen = 0;
      tmplike += like;
/*Every pass through the pedigree increases the scaling by the
  same amount, so we can assign it to per_iter_scale and mutiply
  by nuscale_count. A. A. Schaffer */
#if LOOPSPEED
      if ((per_iter_scale == 0) && (V.nuscale != 0))
        per_iter_scale = V.nuscale;
#endif
      if (risk && like != 0.0)
	riskcumul(&V);
      for (i = 1; i <= totperson; i++) {
	if (person[i]->gen != NULL)
	  cleanup(&person[i], &V);
      }
#if LOOPSPEED
      }
#endif
    }
    alldone = true;
    for (i = 0; i < maxloop; i++)
      alldone = (alldone && V.loopgen[i] == loopmax[i]);
  } while (!alldone);
  like = tmplike;
  if (risk && like != 0.0)
    riskcalc(&V);
  if (ped_nuscales[V.thisped - 1] == 0)
#if LOOPSPEED
    ped_nuscales[V.thisped - 1] = nuscale_count * per_iter_scale;
#else
    ped_nuscales[V.thisped - 1] = V.nuscale;
#endif /*LOOPSPEED*/

  looplastgen = 1;
  for (i = 1; i <= totperson; i++) {
    if (person[i]->gen != NULL)
      cleanup(&person[i], &V);
  }
  if (like == 0.0)
    like = zerolike;
  else
  /*Next line changed by A. A. Schaffer to avoid recomputation of
    V.nuscale*/
    like = log(like) - ped_nuscales[V.thisped - 1] *log(segscale);
  for (i = 0; i < maxloop; i++) {
    if (V.holdpoint[i] != NULL) {
      Free(V.holdpoint[i]);
      V.holdpoint[i] = NULL;
#if LOOPSPEED
      /*Cleanup the extra array to store loop-breaker genotypes*/
      free(loopbreaker_nextgeno[i]);
      loopbreaker_nextgeno[i]= NULL;
#endif
    }
  }
}  /*likelihood*/

void allocate_loopbreaker_vectors()
{
   int ped; 

   ped_nuscales = (int *) malloc(nuped * sizeof(int *));
   for (ped = 1; ped <= nuped; ped++) {
     ped_nuscales[ped - 1] = 0;
   }
}


#if PARALLEL
void allocthetas()
{
   maletheta = (thetavalues *) malloc(sizeof(thetavalues));
   femaletheta = (thetavalues *) malloc(sizeof(thetavalues));
}

/*allocgen is used by the parallel code to allocate the space for genarrays
  and other shred memory data structures*/
void allocgen()
{
   int i, j, k, w;
   double *tempwith1;
   unsigned char *tempflag1;
   unsigned numparams;       /* number of parallel parameters */
   unsigned maxThetas;       /* max number of thetas */

  if (Tmk_proc_id == 0) {
#if PARALLEL_GCENTRAL
    numparams = maxThetas = 2 * maxn;
#elif (defined(MLINK) || defined(LINKMAP))  /* cgh */
    numparams = calcNumparams();
    maxThetas = numIpeds;

    gZeroTheta = (boolean*) Tmk_malloc(sizeof(boolean) * numIpeds);
    if (!gZeroTheta)
      Tmk_errexit("MALLOC ERROR - gZeroTheta\n");
    Tmk_distribute((char*) &gZeroTheta, sizeof(gZeroTheta));

    absoluteThetaIndex = (int*) Tmk_malloc(sizeof(int) * numIpeds);
    if (!absoluteThetaIndex)
      Tmk_errexit("MALLOC ERROR - absoluteThetaIndex\n");
    Tmk_distribute((char*) &absoluteThetaIndex, sizeof(absoluteThetaIndex));

    firstThetanum = (int*) Tmk_malloc(sizeof(int));
    if (!firstThetanum)
      Tmk_errexit("MALLOC ERROR - firstThetanum\n");
    Tmk_distribute((char*) &firstThetanum, sizeof(firstThetanum));

#else  /* if PARALLEL_GCENTRAL -- cgh */
    numparams = maxThetas = maxn;
#endif  /* if PARALLEL_GCENTRAL */
    gmaletheta = (thetarray *) Tmk_malloc(sizeof(thetarray) * maxThetas);
    if (!gmaletheta)
      Tmk_errexit("MALLOC ERROR - gmaletheta\n");
    Tmk_distribute((char *) &gmaletheta,sizeof(gmaletheta));

    gfemaletheta = (thetarray *) Tmk_malloc(sizeof(thetarray) * maxThetas);
    if (!gfemaletheta)
      Tmk_errexit("MALLOC ERROR - gfemaletheta\n");
    Tmk_distribute((char *) &gfemaletheta,sizeof(gfemaletheta));

    gmalesegprob = (double **) Tmk_malloc(sizeof(double *) * numparams);
    if (!gmalesegprob)
      Tmk_errexit("MALLOC ERROR - gmalesegprob\n");
    Tmk_distribute((char *) &gmalesegprob,sizeof(gmalesegprob));
    for(i = 0; i < numparams; i++) {
      gmalesegprob[i] = (double *) Tmk_malloc(sizeof(double) * nuneed);
      if (!gmalesegprob[i])
	Tmk_errexit("MALLOC ERROR - gmalesegprob entry %d", i);
    } 

    gfemalesegprob = (double **) Tmk_malloc(sizeof(double *) * numparams);
    if (!gfemalesegprob)
      Tmk_errexit("MALLOC ERROR - gfemalesegprob\n");
    Tmk_distribute((char *) &gfemalesegprob,sizeof(gfemalesegprob));
    for(i = 0; i < numparams; i++) {
      gfemalesegprob[i] = (double *) Tmk_malloc(sizeof(double) * nuneed);
      if (!gfemalesegprob[i])
	Tmk_errexit("MALLOC ERROR - gfemalesegprob entry %d", i);
    }

    ggeninuse = (char**) Tmk_malloc(numparams * sizeof(char*));
    if (!ggeninuse)
      Tmk_errexit("MALLOC ERROR - ggeninuse\n");
    Tmk_distribute((char *) (&ggeninuse), sizeof(ggeninuse));

    ggenbank = (thisarray**) Tmk_malloc(numparams * sizeof(thisarray*));
    if (!ggenbank)
      Tmk_errexit("MALLOC ERROR - ggenbank\n");
    Tmk_distribute((char *) (&ggenbank), sizeof(ggenbank));

    for(i = 0; i < numparams; i++) {
      ggeninuse[i] = (char*) Tmk_malloc(maxworkingset * sizeof(char));
      if (!ggeninuse[i])
	Tmk_errexit("MALLOC ERROR - ggeninuse[i]\n");

      ggenbank[i] = (thisarray*) Tmk_malloc(maxworkingset * sizeof(thisarray));
      if (!ggenbank[i])
	Tmk_errexit("MALLOC ERROR - ggenbank[%d]\n", i);

      for (w =0; w < maxworkingset; w++) {

	ggenbank[i][w].sparseflag = (unsigned char *) Tmk_malloc(maxfemgen * 
                                                      sizeof(unsigned char));

        if (!(ggenbank[i][w].sparseflag))
          Tmk_errexit("MALLOC ERROR - ggenbank[%d][%d]\n", i, w);

	ggenbank[i][w].genarray = (double *) Tmk_malloc(maxfemgen * 
                                                        sizeof(double));

        if (!(ggenbank[i][w].genarray))
          Tmk_errexit("MALLOC ERROR - ggenbank[%d][%d]\n", i, w);
      }
    }

     for(k = 0; k < numparams; k++)
       for (i = 0; i < maxworkingset; i++) {
	 tempwith1 = (&(ggenbank[k][i]))->genarray;
	 tempflag1 = (&(ggenbank[k][i]))->sparseflag;
	 for(j = 0; j < fgeno; j++) {
	   tempflag1[j] = 0;
	   tempwith1[j] = 0.0;
	 }
	 ggeninuse[k][i] = 0;
       }
    
    arrgene = (double **) Tmk_malloc(Tmk_nprocs * (sizeof(double *)));
    if (!arrgene)
      Tmk_errexit("MALLOC ERROR - arrgene\n");
    Tmk_distribute((char *) &arrgene, sizeof(arrgene));
    
    for (i=0; i < Tmk_nprocs; i++) {
      arrgene[i] = (double *) Tmk_malloc(maxfemgen * sizeof(double));
      if (!arrgene[i])
	Tmk_errexit("MALLOC ERROR - arrgene[%d]\n", i);
    }
    
    arrqsumcache = (double **) Tmk_malloc(Tmk_nprocs * (sizeof(double *)));
    if (!arrqsumcache)
      Tmk_errexit("MALLOC ERROR - arrqsumcache\n");
    Tmk_distribute((char *) &arrqsumcache, sizeof(arrqsumcache));
    
    for (i=0; i < Tmk_nprocs; i++) {
      arrqsumcache[i] = (double *) Tmk_malloc(sizeof(double)*maxhaplo);
      if (!arrqsumcache[i])
	Tmk_errexit("MALLOC ERROR - arrqsumcache[%d]\n", i);
    }

    arrpsumcache = (double **) Tmk_malloc(Tmk_nprocs * (sizeof(double *)));
    if (!arrpsumcache)
      Tmk_errexit("MALLOC ERROR - arrpsumcache\n");
    Tmk_distribute((char *) &arrpsumcache, sizeof(arrpsumcache));
    
    for (i = 0; i < Tmk_nprocs; i++) {
      arrpsumcache[i] = (double *) Tmk_malloc(sizeof(double)*maxhaplo);
      if (!arrpsumcache[i])
	Tmk_errexit("MALLOC ERROR - arrpsumcache[%d]\n", i);
    }

    gpnonzgens = (unsigned int**) Tmk_malloc(numparams * sizeof(unsigned int*));
    if(!gpnonzgens)
      Tmk_errexit("MALLOC ERROR - gpnonzgens\n");
    Tmk_distribute((char *) &gpnonzgens, sizeof(gpnonzgens));
    
    gqnonzgens = (unsigned int**) Tmk_malloc(numparams * sizeof(unsigned int*));
    if(!gqnonzgens)
      Tmk_errexit("MALLOC ERROR - gqnonzgens\n");
    Tmk_distribute((char *) &gqnonzgens, sizeof(gqnonzgens));
    
    for(i = 0; i < numparams; i++) {
      gpnonzgens[i] = (unsigned int *) Tmk_malloc(maxfemgen * sizeof(unsigned int));
      if(!gpnonzgens[i])
	Tmk_errexit("MALLOC ERROR - gpnonzgens[%d]\n", i);
      
      gqnonzgens[i] = (unsigned int *) Tmk_malloc(maxfemgen * sizeof(unsigned int));
      if(!gqnonzgens[i])
	Tmk_errexit("MALLOC ERROR - gqnonzgens[%d]\n", i);
    } 

#if FUNTIMES
#if PARALLEL_GCENTRAL
    executionTimes = (double *) Tmk_malloc(2 * maxn * sizeof(double));
    if (!executionTimes)
      Tmk_errexit("MALLOC ERROR - executionTimes\n");
    Tmk_distribute((char *) &executionTimes, sizeof(executionTimes)); 
#else /* !PARALLEL_GCENTRAL */
    /* either maxn, or numIpeds */
    executionTimes = (double **) Tmk_malloc(maxThetas * sizeof(double*));
    if (!executionTimes)
      Tmk_errexit("MALLOC ERROR - executionTimes\n");
    Tmk_distribute((char *) &executionTimes, sizeof(executionTimes));

    for(i = 0; i < maxThetas; i++) {
     executionTimes[i] = (double *) Tmk_malloc(2 * sizeof(double));
     if (!executionTimes[i])
       Tmk_errexit("MALLOC ERROR - executionTimes[%d]\n",i);
    } 
#endif /* PARALLEL_GCENTRAL */
#endif /*FUNTIMES*/

#if ALLELE_SPEED
    sharedCurrentPed = (int *) Tmk_malloc(Tmk_nprocs * sizeof(int));
    if (!sharedCurrentPed)
      Tmk_errexit("MALLOC ERROR - sharedCurrentPed\n");
    Tmk_distribute((char *) &sharedCurrentPed, sizeof(sharedCurrentPed));   

    for(i=0; i < Tmk_nprocs; i++)
      sharedCurrentPed[i] = 0;
#endif /*ALLELE_SPEED*/
  }
}

#endif /* PARALLEL */


/* Memory-estimation routines  */

static void mem_collapsedown();

Local Void mem_cleanup(p)
thisperson **p;
{
  thisperson *WITH;

  WITH = *p;
  if (!(WITH->loopneeded)) {
    WITH->newgenexists = false;
    memcount --;
  }
}  /*mem_cleanup*/


Local Void mem_prob(p)
thisperson **p;
{
  thisperson *WITH;


  WITH = *p;
  if (WITH->newgenexists) return;  /*Fixed by A. A. Schaffer*/
  WITH->newgenexists = true;
  memcount++;
  if (memcount > maxmemcount)
    maxmemcount = memcount;
} /* mem_prob */

Local Void mem_initseg(father, mother)
thisperson *father, *mother;
{
  thisperson *child;

  mem_prob(&father);
  mem_prob(&mother);
  child = father->foff;
  nchild = 0;
  do {
    mem_prob(&child);
    if (child->ma == mother && (!child->up)) {
      nchild++;
      childarray[nchild - 1] = child;  /*Line added by A. A. Schaffer*/
    }
    child = child->nextpa;
  } while (child != NULL);  
  if (nchild > maxchild) {
    fprintf(stderr, "\nA nuclear family has more children than maxchild");
    fprintf(stderr, "\nThe program will exit politely to allow you to");
    fprintf(stderr, "\nincrease maxchild in the file commondefs.h");
    exit(EXIT_FAILURE);
  }
} /* mem_initseg*/


Local Void mem_exitseg(father, mother)
thisperson *father, *mother;
{
  thisperson *child;

  child = father->foff;
  do {
    if (child->ma == mother && (!child->up))
      mem_cleanup(&child);
    child = child->nextpa;
  } while (child != NULL); 
}   /*mem_exitseg*/

Local Void mem_segup(p, q, child)
thisperson *p, *q, *child;
{
  thisperson *father, *mother;
  boolean depend;
  int i;

  if (p->male) {
    father = p;
    mother = q;
  }
  else {
    father = q;
    mother = p;
  }
  mem_initseg(father, mother);
  if (!p->loopdepend) {
     depend = false;
     for(i=0; i<nchild; i++)
       if (childarray[i]->loopdepend) {
         depend = true;
         break;
       }
     depend = (depend || (q->loopdepend));
     if (depend) {
       p->loopdepend = depend;
       p->loopneeded = false;
     }
  }
  if (p->loopdepend) {
     if (!(q->loopdepend))
       (q->loopneeded) = true;
     for(i=0; i< nchild; i++)
       if (!(childarray[i]->loopdepend))
         childarray[i]->loopneeded = true;
  }
  mem_cleanup(&q);
  mem_exitseg(father, mother);
} /* mem_segup*/

Local Void mem_segdown(p, q, child)
thisperson *p, *q, *child;
{
    boolean depend;
    int i;

    mem_initseg(p, q);
    if (!(child->loopdepend)) {
      depend = false;
      for(i=0; i<nchild; i++)
        if (childarray[i]->loopdepend) {
          depend = true;
          break;
        }
      depend = (depend || (p->loopdepend) || (q->loopdepend));
      if (depend) {
        child->loopdepend = depend;
        child->loopneeded = false;
      }
    }
    if (child->loopdepend) {
     if (!(p->loopdepend))
       p->loopneeded = true;
     if (!(q->loopdepend))
       q->loopneeded = true;
     for(i=0; i< nchild; i++)
       if (!(childarray[i]->loopdepend))
         childarray[i]->loopneeded = true;
   }
  mem_cleanup(&p);
  mem_cleanup(&q);
  mem_exitseg(p, q);
} /* mem_segdown*/

Local Void mem_collapseup(p)
thisperson *p;
{
  thisperson *q, *child, *nextchild;
  boolean down;
  p->done = true;
  if (p->foff == NULL)
    return;
  down = false;
  child = p->foff;
  while (child != NULL) {
    down = false;
    if (p->male)
      q = child->ma;
    else
      q = child->pa;
    if (!q->done) {
      mem_collapsedown(q);
      nextchild = child;
      while (nextchild != NULL) {
	if (nextchild->pa == q || nextchild->ma == q) {
	  if (!nextchild->up)
	    mem_collapseup(nextchild);
	  else
	    down = true;
	}
	if (p->male)
	  nextchild = nextchild->nextpa;
	else
	  nextchild = nextchild->nextma;
      }
      if (q->multi)
	mem_collapseup(q);
      if (!down)
	mem_segup(p, q, child);
      else
	mem_collapsedown(p);
    }
    if (p->male)
      child = child->nextpa;
    else
      child = child->nextma;
  }
}  /*mem_collapseup*/


Local Void mem_collapsedown(p)
thisperson *p;
{
  if (p->pa == NULL)
    return;
  p->up = true;
  mem_collapseup(p->pa);
  mem_segdown(p->pa, p->ma, p);
}  /*collapsedown*/


/*Written by A. A. Schaffer to estimate the maximum size of the cutset
  so as to know how many thisarray's will be needed*/
int maxw_estimation()

{
  int maxcount;
  int ped, id;
#if LOOP_BREAKERS
  int copyIndex;
#endif  

  maxcount = 0;
  for(ped = 1; ped <=nuped; ped++) {
    memcount = maxmemcount = 0;
    for(id = 1; id <= totperson; id++) {
      person[id]->newgenexists = false;
      person[id]->up = false;
      person[id]->done = false;
      person[id]->loopdepend = false;
      person[id]->loopneeded = false;      
    }
    if(looppers[ped -1][0][0] != NULL) {  /*G OK*/
#if LOOP_BREAKERS
      for(copyIndex = 0; copyIndex < numCopies[ped-1][0]; copyIndex++)
	pollutedescendants(looppers[ped-1][0][copyIndex]); /*G OK*/
#else
      pollutedescendants(looppers[ped-1][0][0]); /*G OK */
      pollutedescendants(looppers[ped-1][0][1]);  /*G OK*/
#endif /*LOOP_BREAKERS*/
      precollapseup(proband[ped - 1]);
      precollapsedown(proband[ped - 1]);
      for(id = 1; id <= totperson; id++) {
        person[id]->up = false;
        person[id]->done = false;
      }
      /*do one extra traversal for loops to simulate second genotype of
        loopbreaker*/
      mem_collapseup(proband[ped - 1]);
      mem_collapsedown(proband[ped - 1]);
      for(id = 1; id <= totperson; id++) {
        person[id]->up = false;
        person[id]->done = false;
        if (person[id]->newgenexists && (!(person[id]->loopneeded)))
          mem_cleanup(&(person[id]));
      }
    }
    mem_collapseup(proband[ped - 1]);
    mem_collapsedown(proband[ped - 1]);
    if (maxmemcount > maxcount)
      maxcount = maxmemcount;
  }
  for(id = 1; id <= totperson; id++) {
    person[id]->newgenexists = false;
    person[id]->up = false;
    person[id]->done = false;
    person[id]->loopdepend = false;
    person[id]->loopneeded = false;      
  }
  return(maxcount);
}

