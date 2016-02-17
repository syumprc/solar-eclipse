/*This file is part of FASTLINK, version 2.2.
The algorithms in FASTLINK are described in: R. W. Cottingham
Jr., R. M. Idury, A. A. Schaffer (1993), Faster sequential genetic
linkage computations, American Journal of Human Genetics 53(1993), pp.
252-263. and
A. A. Schaffer, S. K. Gupta, K. Shriram, and R. W. Cottingham Jr.,
Avoiding Recomputation in Linkage Analysis, Hum. Hered. 44(1994), pp. 225-237.
This file contains low level probability routines.*/

/*
   Modified by Dylan in late 1994 to speed up pedigrees with loops.  Main
   difference is that for unknown people the possible genotypes are stored
   in 'unknown_poss' rather than their 'possible' fields.

   Please read README.loopfile.  See also unknown.c and commondefs.c.
*/


/* Local variables for setval: */
struct LOC_setval {
  struct LOC_getvect *LINK;
  double val;
  int nhap1, nhap2;
} ;

static void prior(LINK)
struct LOC_setval *LINK;
{
  int i;   /*prior*/
  locusvalues *WITH;
  thisarray *WITH1;
  double tempfreq1, tempfreq2; /*frequency of allele*/
  boolean condhap1, condhap2; 
#if defined(LODSCORE)
  int locusindex;
#endif

  condhap1 = condhap2 = false;
  if (disfreqs) {
    for (i = 0; i < mlocus; i++)
      if (affection == thislocus[i]->which) {
     	if (affall == LINK->LINK->hap1[i])
	  condhap1 = true;
        if (((!sexlink ) || (!(LINK->LINK->p->male))) &&
     	    (affall == LINK->LINK->hap2[i]))
	  condhap2 = true;
	break;
      }  
  }
  if (!disequi) {
    if (sexlink && LINK->LINK->p->male) {
      for (i = 0; i < mlocus; i++) {
	WITH = thislocus[i];
#if ALLELE_SPEED
        if ((binary_ == thislocus[i]->which) &&
            (allformat == thislocus[i]->format)) {
#if defined(LODSCORE)
	  if (0 == i)
	    locusindex = locuslist1[iplace - 1] - 1;
	  else
	    locusindex = locuslist2[jplace - 1] - 1;
#endif     /*defined(LODSCORE) */
#if defined(LODSCORE)
	tempfreq1 = ped_loc_all[currentped - 1][locusindex]
		[LINK->LINK->hap1[i]].new_frequency;
        if (disfreqs && condhap1)
	  tempfreq1 = ped_loc_all[currentped - 1][locusindex]
	    [LINK->LINK->hap1[i]].new_frequency_cond;
#else   /*defined(LODSCORE) */
        tempfreq1 =  ped_loc_all[currentped - 1][i]
		[LINK->LINK->hap1[i]].new_frequency;
        if (disfreqs && condhap1)
	  tempfreq1 = ped_loc_all[currentped - 1][i]
	    [LINK->LINK->hap1[i]].new_frequency_cond;
#endif  /*defined(LODSCORE) */
	}
        else
#endif  /*ALLELE_SPEED*/
          if (!disfreqs || (binary_ != thislocus[i]->which) || (!condhap1))
	    tempfreq1 = WITH->freq[LINK->LINK->hap1[i] - 1];
	  else
	    tempfreq1 = WITH->freqcond[LINK->LINK->hap1[i] - 1];
	LINK->val *= tempfreq1;
      }
    } else {
      for (i = 0; i < mlocus; i++) {
	WITH = thislocus[i];
#if ALLELE_SPEED
        if ((binary_ == thislocus[i]->which) &&
            (allformat == thislocus[i]->format)) {
#if defined(LODSCORE)
	  if (0 == i)
	    locusindex = locuslist1[iplace - 1] - 1;
	  else
	    locusindex = locuslist2[jplace - 1] - 1;
#endif     
#if defined(LODSCORE)
          tempfreq1 = 
              ped_loc_all[currentped - 1][locusindex]
		[LINK->LINK->hap1[i]].new_frequency;
          tempfreq2 = 
              ped_loc_all[currentped - 1][locusindex]
		[LINK->LINK->hap2[i]].new_frequency;
          if (disfreqs) {
	    if (condhap1)
	      tempfreq1 = 
		ped_loc_all[currentped - 1][locusindex]
		 [LINK->LINK->hap1[i]].new_frequency_cond;
	    if (condhap2)
	      tempfreq2 = 
		ped_loc_all[currentped - 1][locusindex]
		[LINK->LINK->hap2[i]].new_frequency_cond;
	  }
#else
          tempfreq1 = 
              ped_loc_all[currentped - 1][i]
		[LINK->LINK->hap1[i]].new_frequency;
          tempfreq2 = 
              ped_loc_all[currentped - 1][i]
		[LINK->LINK->hap2[i]].new_frequency;
          if (disfreqs) {
	    if (condhap1)
	      tempfreq1 = 
		ped_loc_all[currentped - 1][i]
		 [LINK->LINK->hap1[i]].new_frequency_cond;
	    if (condhap2)
	      tempfreq2 = 
		ped_loc_all[currentped - 1][i]
		[LINK->LINK->hap2[i]].new_frequency_cond;
	  }
#endif
        }
        else {
#endif
          tempfreq1 = WITH->freq[LINK->LINK->hap1[i] - 1];
          tempfreq2 = WITH->freq[LINK->LINK->hap2[i] - 1];
          if (disfreqs && (binary_ == thislocus[i]->which)) {
	    if (condhap1)
	      tempfreq1 = WITH->freqcond[LINK->LINK->hap1[i] - 1];
	    if (condhap2)
	      tempfreq2 = WITH->freqcond[LINK->LINK->hap2[i] - 1];
          }
#if ALLELE_SPEED
        }
#endif
	LINK->val *= tempfreq1 * tempfreq2;
      }
      if (LINK->nhap1 != LINK->nhap2)
	LINK->val = 2.0 * LINK->val;
    }
  } else {
    WITH1 = hapfreq;
    if (sexlink && LINK->LINK->p->male)
      LINK->val *= WITH1->genarray[LINK->nhap1 - 1];
    else {
      LINK->val *= WITH1->genarray[LINK->nhap1 - 1] *
		   WITH1->genarray[LINK->nhap2 - 1];
      if (LINK->nhap1 != LINK->nhap2)
	LINK->val = 2.0 * LINK->val;
    }
  }
  LINK->val *= segscale;
}  /* prior */

/*getval*/

static void setval(val_, LINK)
double val_;
struct LOC_getvect *LINK;
{
  struct LOC_setval V;
  int here, count, i;
  thisarray *WITH1;

  here = 0;  /*cgh -- to avoid gcc warning*/
  V.LINK = LINK;
  V.val = val_;
  count = 1;
  V.nhap1 = 1;
  V.nhap2 = 1;
  if (sexlink && LINK->p->male) {
    /* 
       WARNING!!!

       If you change the way that the joint haplotype is computed
       you need to change translate_loop_vector().

       Dylan
    */
    for (i = 0; i < mlocus; i++)
      V.nhap1 += increment[i] * (LINK->hap1[i] - 1);
    here = V.nhap1;
  } else {
    for (i = 0; i < mlocus; i++) {
      V.nhap1 += increment[i] * (LINK->hap1[i] - 1);
      V.nhap2 += increment[i] * (LINK->hap2[i] - 1);
      if (LINK->hap1[i] != LINK->hap2[i])
	count *= 2;
      here = genenumber[V.nhap1 - 1][V.nhap2 - 1];
    }
  }
  if (LINK->p->pa == NULL)
    prior(&V);
  WITH1 = LINK->p->gen;
  if (disequi || disfreqs) {
    WITH1->genarray[here - 1] = V.val;
    WITH1->sparseflag[here - 1] = 1; /*R. M. Idury, A. A. Schaffer*/
    return;
  }
  if (count != 1)
    count /= 2;
  for (i = 1; i <= count; i++) {
    WITH1->genarray[here - 1] = V.val;
    WITH1->sparseflag[here - 1] = 1; /*R. M. Idury, A. A. Schaffer*/
    here++;
  }
}  /* setval */

/* Local variables for getgene: */
struct LOC_getgene {
  struct LOC_getvect *LINK;
  int syste;
  double val;
  double newval;
} ;

static void facmale(LINK)
struct LOC_getgene *LINK;
{
  int j;   /*facmale*/
  thisperson *WITH;
  locusvalues *WITH1;
  int FORLIM;

  WITH = LINK->LINK->p;
  WITH1 = thislocus[LINK->syste - 1];
  FORLIM = WITH1->nallele;
  for (j = 1; j <= FORLIM; j++) {
    if (((binformat == WITH1->format) && 
         (WITH->phen[LINK->syste - 1]->phenf == WITH1->UU.allele[j - 1] ||
	  WITH->phen[LINK->syste - 1]->phenf == 0)) ||
        ((allformat == WITH1->format) &&
         (WITH->phen[LINK->syste - 1]->allele1 == j ||
	  WITH->phen[LINK->syste - 1]->allele1 == 0))) {
      LINK->LINK->hap1[LINK->syste - 1] = j;
      if (LINK->syste != mlocus)
	getgene(LINK->syste + 1, LINK->val, LINK->LINK);
      else
	setval(LINK->val, LINK->LINK);
    }
  }
}  /* facmale */

/*facmale*/

static void affmale(LINK)
struct LOC_getgene *LINK;
{
  int j;   /*affmale*/
  locusvalues *WITH;
  int FORLIM;

  WITH = thislocus[LINK->syste - 1];
  FORLIM = WITH->nallele;
  for (j = 1; j <= FORLIM; j++) {
    LINK->newval = LINK->val;
    getval(LINK->syste, 0L, j, &LINK->newval, LINK->LINK);
    LINK->LINK->hap1[LINK->syste - 1] = j;
    if (LINK->newval != 0.0) {
      if (LINK->syste != mlocus)
	getgene(LINK->syste + 1, LINK->newval, LINK->LINK);
      else
	setval(LINK->newval, LINK->LINK);
    }
  }
}  /* affmale */

/*affmale*/

static void quanmale(LINK)
struct LOC_getgene *LINK;
{
  int j;   /*quanmale*/
  thisperson *WITH;
  locusvalues *WITH1;
  int FORLIM;

  WITH = LINK->LINK->p;
  WITH1 = thislocus[LINK->syste - 1];
  if (WITH->phen[LINK->syste - 1]->aff == affall ||
      WITH->phen[LINK->syste - 1]->aff == missaff) {
    LINK->newval = LINK->val;
    LINK->LINK->hap1[LINK->syste - 1] = affall;
    if (LINK->newval != 0.0) {
      if (LINK->syste != mlocus)
	getgene(LINK->syste + 1, LINK->newval, LINK->LINK);
      else
	setval(LINK->newval, LINK->LINK);
    }
  }
  if (WITH->phen[LINK->syste - 1]->aff == affall &&
      WITH->phen[LINK->syste - 1]->aff != missaff)
    return;
  FORLIM = WITH1->nallele;
  for (j = 1; j <= FORLIM; j++) {
    if (j != affall) {
      LINK->newval = LINK->val;
      LINK->LINK->hap1[LINK->syste - 1] = j;
      if (LINK->newval != 0.0) {
	if (LINK->syste != mlocus)
	  getgene(LINK->syste + 1, LINK->newval, LINK->LINK);
	else
	  setval(LINK->newval, LINK->LINK);
      }
    }
  }
}  /* quanmale */

/*quanmale*/

static void fac(LINK)
struct LOC_getgene *LINK;
{
  int i, j;   /*fac*/
  thisperson *WITH;
  locusvalues *WITH1;
  int FORLIM, FORLIM1;

  WITH = LINK->LINK->p;
  WITH1 = thislocus[LINK->syste - 1];
  FORLIM = WITH1->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap1[LINK->syste - 1] = i;
    FORLIM1 = WITH1->nallele;
    for (j = i; j <= FORLIM1; j++) {
      if (((binformat == WITH1->format) &&
          ((WITH->phen[LINK->syste - 1]->phenf ==
	  (WITH1->UU.allele[i - 1] | WITH1->UU.allele[j - 1]) ||
	  WITH->phen[LINK->syste - 1]->phenf == 0))) ||
         ((allformat == WITH1->format) &&
         ((WITH->phen[LINK->syste - 1]->allele1 == i &&
           WITH->phen[LINK->syste - 1]->allele2 == j) ||
           WITH->phen[LINK->syste - 1]->allele1 == 0))) {
	LINK->LINK->hap2[LINK->syste - 1] = j;
	if (LINK->syste != mlocus)
	  getgene(LINK->syste + 1, LINK->val, LINK->LINK);
	else
	  setval(LINK->val, LINK->LINK);
      }
    }
  }
  if ((!disequi) && (!disfreqs))
    return;
  WITH = LINK->LINK->p;
  WITH1 = thislocus[LINK->syste - 1];
  FORLIM = WITH1->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap2[LINK->syste - 1] = i;
    FORLIM1 = WITH1->nallele;
    for (j = i; j <= FORLIM1; j++) {
        if (((binformat == WITH1-> format) &&
            (WITH->phen[LINK->syste - 1]->phenf ==
	     (WITH1->UU.allele[i - 1] | WITH1->UU.allele[j - 1]) ||
	     WITH->phen[LINK->syste - 1]->phenf == 0)) ||
            ((allformat == WITH1->format) &&
	     ((WITH->phen[LINK->syste - 1]->allele1 == i &&
	       WITH->phen[LINK->syste - 1]->allele2 == j) ||
	      WITH->phen[LINK->syste - 1]->allele1 == 0))) {
	LINK->LINK->hap1[LINK->syste - 1] = j;
	if (LINK->syste != mlocus)
	  getgene(LINK->syste + 1, LINK->val, LINK->LINK);
	else
	  setval(LINK->val, LINK->LINK);
      }
    }
  }
}  /* fac */

/*fac*/

static void aff(LINK)
struct LOC_getgene *LINK;
{
  /*Used with an affection status phenotype or when
  thislocus[syste]^which is null*/
  int i, j;
  locusvalues *WITH;
  int FORLIM, FORLIM1;

  WITH = thislocus[LINK->syste - 1];
  FORLIM = WITH->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap1[LINK->syste - 1] = i;
    FORLIM1 = WITH->nallele;
    for (j = i; j <= FORLIM1; j++) {
      if (!nohom[LINK->syste - 1] || i != affall || j != affall) {
	LINK->LINK->hap2[LINK->syste - 1] = j;
	LINK->newval = LINK->val;
	getval(LINK->syste, i, j, &LINK->newval, LINK->LINK);
	if (LINK->newval != 0.0) {
	  if (LINK->syste != mlocus)
	    getgene(LINK->syste + 1, LINK->newval, LINK->LINK);
	  else
	    setval(LINK->newval, LINK->LINK);
	}
      }
    }
  }
  if ((!disequi) && (!disfreqs))
    return;
  WITH = thislocus[LINK->syste - 1];
  FORLIM = WITH->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap2[LINK->syste - 1] = i;
    FORLIM1 = WITH->nallele;
    for (j = i; j <= FORLIM1; j++) {
      if (!nohom[LINK->syste - 1] || i != affall || j != affall) {
	LINK->LINK->hap1[LINK->syste - 1] = j;
	LINK->newval = LINK->val;
	getval(LINK->syste, i, j, &LINK->newval, LINK->LINK);
	if (LINK->newval != 0.0) {
	  if (LINK->syste != mlocus)
	    getgene(LINK->syste + 1, LINK->newval, LINK->LINK);
	  else
	    setval(LINK->newval, LINK->LINK);
	}
      }
    }
  }
}  /* aff */


static void quanval(LINK)
struct LOC_getgene *LINK;
{
  /*Uses this only when thislocus[syste]^.which is not null*/
  int i, j;
  locusvalues *WITH;
  int FORLIM, FORLIM1;

  WITH = thislocus[LINK->syste - 1];
  FORLIM = WITH->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap1[LINK->syste - 1] = i;
    FORLIM1 = WITH->nallele;
    for (j = i; j <= FORLIM1; j++) {
      if (!nohom[LINK->syste - 1] || i != affall || j != affall) {
	LINK->LINK->hap2[LINK->syste - 1] = j;
	LINK->newval = LINK->val;
	getval(LINK->syste, i, j, &LINK->newval, LINK->LINK);
	if (LINK->newval != 0.0) {
	  if (LINK->syste != mlocus)
	    getgene(LINK->syste + 1, LINK->newval, LINK->LINK);
	  else
	    setval(LINK->newval, LINK->LINK);
	}
      }
    }
  }
  if ((!disequi) && (!disfreqs))
    return;
  WITH = thislocus[LINK->syste - 1];
  FORLIM = WITH->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap2[LINK->syste - 1] = i;
    FORLIM1 = WITH->nallele;
    for (j = i; j <= FORLIM1; j++) {
      if (!nohom[LINK->syste - 1] || i != affall || j != affall) {
	LINK->LINK->hap1[LINK->syste - 1] = j;
	LINK->newval = LINK->val;
	getval(LINK->syste, i, j, &LINK->newval, LINK->LINK);
	if (LINK->newval != 0.0) {
	  if (LINK->syste != mlocus)
	    getgene(LINK->syste + 1, LINK->newval, LINK->LINK);
	  else
	    setval(LINK->newval, LINK->LINK);
	}
      }
    }
  }
}  /* quanval */

/*setval*/

static void getgene(syste_, val_, LINK)
int syste_;
double val_;
struct LOC_getvect *LINK;
{
  struct LOC_getgene V;
  locusvalues *WITH;


  V.LINK = LINK;
  V.syste = syste_;
  V.val = val_;
  WITH = thislocus[V.syste - 1];
  if (sexlink && LINK->p->male) {
    switch (WITH->which) {

    case binary_:
      facmale(&V);
      break;

    case affection:
      affmale(&V);
      break;

    case quantitative:
      quanmale(&V);
      break;

    case null_:
      if (WITH->privlocus->which == affection)
	affmale(&V);
      else
	quanmale(&V);
      break;
    }
    return;
  }
  switch (WITH->which) {

  case binary_:
    fac(&V);
    break;

  case affection:
    aff(&V);
    break;

  case quantitative:
    quanval(&V);
    break;

  case null_:
    aff(&V);
    break;
  }
}  /* getgene */

/* Local variables for ugetgene: */
struct LOC_ugetgene {
  struct LOC_getvect *LINK;
  int syste;
  double val;
  double newval;
} ;

/* 
   This procedure finds the possible allele combinations for the
   current person at the current locus.  The person is unknown and
   male the locus is a binary locus.  The run is sexlinked.

   Modified by Dylan in late 1994 to use 'unknown_poss'.  
*/
static void facmale_(LINK)
struct LOC_ugetgene *LINK;
{
  int j;   /*facmale*/
  thisperson *WITH;
#if LOOPSPEED
  int loop_vect_to_check;
#else
  information *WITH1;
#endif
  locusvalues *WITH2;
  int locus;

  locus = LINK->syste - 1;
  WITH = LINK->LINK->p;

#if LOOPSPEED
  /* if known, possibles held in loop vector 0 */
  if ( person[WITH->id]->unknown)
      loop_vect_to_check= single_locus_vector_num[locus];
  else {
    loop_vect_to_check = 0;
  }
#else
  WITH1 = LINK->LINK->p->store;
#endif

  WITH2 = thislocus[locus];

  for (j = 1; j <= WITH2->nallele; j++) {
    if (((binformat == WITH2->format) &&
        (WITH->phen[locus]->phenf == WITH2->UU.allele[j - 1] ||
	WITH->phen[locus]->phenf == 0)) ||
        ((allformat == WITH2->format) &&
         (WITH->phen[LINK->syste - 1]->allele1 == j ||
	  WITH->phen[LINK->syste - 1]->allele1 == 0))) { 
#if LOOPSPEED
      if (unknown_poss[LINK->LINK->p->id][locus][loop_vect_to_check][j - 1] ) {
#else
      if (WITH1->possible[locus][0][j - 1]) {
#endif
	LINK->LINK->hap1[locus] = j;
	if (LINK->syste != mlocus)
	  /* do next locus */
	  ugetgene(locus + 2, LINK->val, LINK->LINK);
	else
	  setval(LINK->val, LINK->LINK);
      }
    }
  }
}  /* facmale_ */

/* 
   This procedure finds the possible allele combinations for the
   current person at the current locus.  The person is unknown and
   male and the locus is an affection locus.  The run is sexlinked.

   Modified by Dylan in late 1994 to use 'unknown_poss'.  
*/
static void affmale_(LINK)
struct LOC_ugetgene *LINK;
{
  int j;
#if LOOPSPEED
  int loop_vect_to_check;
#else
  information *WITH;
#endif
  locusvalues *WITH1;
  int locus;

  locus = LINK->syste - 1;
  WITH1 = thislocus[locus];

#if LOOPSPEED
  /* if known, possibles held in loop vector 0 */
  if ( person[LINK->LINK->p->id]->unknown)
      loop_vect_to_check= single_locus_vector_num[locus];
  else {
    loop_vect_to_check = 0;
  }
#else
  WITH = LINK->LINK->p->store;
#endif

  for (j = 1; j <= WITH1->nallele; j++) {
#if LOOPSPEED
    if (unknown_poss[LINK->LINK->p->id][locus][loop_vect_to_check][j - 1] ) {
#else
    if (WITH->possible[locus][0][j - 1]) {
#endif
      LINK->newval = LINK->val;
      getval(locus + 1, 0L, j, &LINK->newval, LINK->LINK);
      LINK->LINK->hap1[locus] = j;
      if (LINK->newval != 0.0) {
	if (LINK->syste != mlocus)
	  ugetgene(locus + 2, LINK->newval, LINK->LINK);
	else
	  setval(LINK->newval, LINK->LINK);
      }
    }
  }
}  /* affmale_ */

/* 
   This procedure finds the possible allele combinations for the
   current person at the current locus.  The person is unknown and
   male the locus is a quantified locus.  The run is sexlinked.

   Modified by Dylan in late 1994 to use 'unknown_poss'.  
*/
static void quanmale_(LINK)
struct LOC_ugetgene *LINK;
{
  int j;
  thisperson *WITH;
  locusvalues *WITH2;
  int locus;
#if LOOPSPEED
  int loop_vect_to_check;
#else
  information *WITH1;
  information *WITH3;
#endif

  locus = LINK->syste - 1;
  WITH = LINK->LINK->p;

#if LOOPSPEED
  /* if known, possibles held in loop vector 0 */
  if ( person[WITH->id]->unknown)
      loop_vect_to_check= single_locus_vector_num[locus];
  else {
    loop_vect_to_check = 0;
  }
#else
  WITH1 = LINK->LINK->p->store;
#endif

  WITH2 = thislocus[locus];
  if (WITH->phen[locus]->aff == affall ||
      WITH->phen[locus]->aff == missaff) {
#if LOOPSPEED
    if (unknown_poss[LINK->LINK->p->id][locus][loop_vect_to_check][affall-1]) {
#else
    if (WITH1->possible[locus][0][affall - 1]) {
#endif
      LINK->newval = LINK->val;
      LINK->LINK->hap1[locus] = affall;
      if (LINK->newval != 0.0) {
	if (LINK->syste != mlocus)
	  ugetgene(locus + 2, LINK->newval, LINK->LINK);
	else
	  setval(LINK->newval, LINK->LINK);
      }
    }
  }

  if (WITH->phen[locus]->aff == affall && WITH->phen[locus]->aff != missaff)
    return;

#if LOOPSPEED
  /* if known, possibles held in loop vector 0 */
  if ( person[WITH->id]->unknown)
      loop_vect_to_check= single_locus_vector_num[locus];
  else {
    loop_vect_to_check = 0;
  }
#else
  WITH3 = LINK->LINK->p->store;
#endif

  for (j = 1; j <= WITH2->nallele; j++) {
    if (j != affall) {
#if LOOPSPEED
    if (unknown_poss[LINK->LINK->p->id][locus][loop_vect_to_check][j - 1] ) {
#else
      if (WITH3->possible[locus][0][j - 1]) {
#endif
	LINK->newval = LINK->val;
	LINK->LINK->hap1[locus] = j;
	if (LINK->newval != 0.0) {
	  if (LINK->syste != mlocus)
	    ugetgene(locus + 2, LINK->newval, LINK->LINK);
	  else
	    setval(LINK->newval, LINK->LINK);
	}
      }
    }
  }
}  /* quanmale_ */


/* 
   This procedure finds the possible allele combinations for the
   current person at the current locus.  The person is unknown and the
   locus is a binary locus.  The run is not sexlinked and/or the
   person is female.

   Modified by Dylan in late 1994 to use 'unknown_poss'.  
*/
static void fac_(LINK)
struct LOC_ugetgene *LINK;
{
  int i, j;
  thisperson *WITH;
#if LOOPSPEED
  int duplicates;  
  int loop_vect_to_check;
#else
  information *WITH1;
#endif
  locusvalues *WITH2;
  int FORLIM;
  int locus;

  WITH = LINK->LINK->p;
  locus = LINK->syste - 1;

#if LOOPSPEED
  /* if known, possibles held in loop vector 0 */
  if (person[WITH->id]->unknown)
      loop_vect_to_check= single_locus_vector_num[locus];
  else {
    loop_vect_to_check = 0;
  }
#else
  WITH1 = LINK->LINK->p->store;
#endif

  WITH2 = thislocus[locus];
  FORLIM = WITH2->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap1[locus] = i;
    for (j = i; j <= FORLIM; j++) {
      if (((binformat == WITH2->format) &&
         ((WITH->phen[locus]->phenf ==
	   (WITH2->UU.allele[i - 1] | WITH2->UU.allele[j - 1]) ||
	   WITH->phen[locus]->phenf == 0))) ||
         ((allformat == WITH2->format) &&
          ((WITH->phen[LINK->syste - 1]->allele1 == i &&
           WITH->phen[LINK->syste - 1]->allele2 == j) ||
           WITH->phen[LINK->syste - 1]->allele1 == 0))) { 
#if LOOPSPEED
	duplicates = (i - 1) * i / 2;
	if (unknown_poss[LINK->LINK->p->id][locus][loop_vect_to_check]
	    [ (i - 1) * thislocus[locus]->nallele - duplicates + (j - 1)] ) {
#else
        if (WITH1->possible[locus][i - 1][j - 1]) {
#endif
	  LINK->LINK->hap2[locus] = j;
	  if (LINK->syste != mlocus)
	    ugetgene(locus + 2, LINK->val, LINK->LINK);
	  else
	    setval(LINK->val, LINK->LINK);
	}
      }
    }
  }

 if ((!disequi) && (!disfreqs))
    return;

  WITH = LINK->LINK->p;
  locus = LINK->syste - 1;

#if LOOPSPEED
  /* if known, possibles held in loop vector 0 */
  if ( person[WITH->id]->unknown)
      loop_vect_to_check= single_locus_vector_num[locus];
  else {
    loop_vect_to_check = 0;
  }
#else
  WITH1 = LINK->LINK->p->store;
#endif

  WITH2 = thislocus[locus];
  FORLIM = WITH2->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap2[locus] = i;
    for (j = i; j <= FORLIM; j++) {
      if (((binformat == WITH2->format) &&
	   (WITH->phen[locus]->phenf ==
	    (WITH2->UU.allele[i - 1] | WITH2->UU.allele[j - 1]) ||
	    WITH->phen[locus]->phenf == 0)) ||
         ((allformat == WITH2->format) &&
	  ((WITH->phen[LINK->syste - 1]->allele1 == i &&
	    WITH->phen[LINK->syste - 1]->allele2 == j) ||
           WITH->phen[LINK->syste - 1]->allele1 == 0))) {
#if LOOPSPEED
	duplicates = (i - 1) * i / 2;
	if (unknown_poss[LINK->LINK->p->id][locus][loop_vect_to_check]
	    [ (i - 1) * thislocus[locus]->nallele - duplicates + (j - 1)] ) {
#else
	if (WITH1->possible[locus][i - 1][j - 1]) {
#endif
	  LINK->LINK->hap1[locus] = j;
	  if (LINK->syste != mlocus)
	    ugetgene(LINK->syste + 1, LINK->val, LINK->LINK);
	  else
	    setval(LINK->val, LINK->LINK);
	}
      }
    }
  }
}  /* fac_ */

/* 
   This procedure finds the possible allele combinations for the
   current person at the current locus.  The person is unknown and the
   locus is a quantified locus.  The run is not sexlinked and/or the
   person is female.

   Modified by Dylan in late 1994 to compile when LOOPSPEED is defined.
   This was identical to aff_() so I only modified aff_();
*/
#if !LOOPSPEED
static void quanval_(LINK)
struct LOC_ugetgene *LINK;
{
  /*Used with an affection status phenotype or when
  thislocus[syste]^which is null*/
  int i, j;
  thisperson *WITH;
#if !LOOPSPEED
  information *WITH1;
#endif
  locusvalues *WITH2;
  int FORLIM, FORLIM1;

  WITH = LINK->LINK->p;
  WITH1 = LINK->LINK->p->store;
  WITH2 = thislocus[LINK->syste - 1];
  FORLIM = WITH2->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap1[LINK->syste - 1] = i;
    FORLIM1 = WITH2->nallele;
    for (j = i; j <= FORLIM1; j++) {
      if (!nohom[LINK->syste - 1] || i != affall || j != affall) {
	if (WITH1->possible[LINK->syste - 1][i - 1][j - 1]) {
	  LINK->LINK->hap2[LINK->syste - 1] = j;
	  LINK->newval = LINK->val;
	  getval(LINK->syste, i, j, &LINK->newval, LINK->LINK);
	  if (LINK->newval != 0.0) {
	    if (LINK->syste != mlocus)
	      ugetgene(LINK->syste + 1, LINK->newval, LINK->LINK);
	    else
	      setval(LINK->newval, LINK->LINK);
	  }
	}
      }
    }
  }
  if ((!disequi) && (!disfreqs))
    return;
  WITH = LINK->LINK->p;
  WITH1 = LINK->LINK->p->store;
  WITH2 = thislocus[LINK->syste - 1];
  FORLIM = WITH2->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap2[LINK->syste - 1] = i;
    FORLIM1 = WITH2->nallele;
    for (j = i; j <= FORLIM1; j++) {
      if (!nohom[LINK->syste - 1] || i != affall || j != affall) {
	if (WITH1->possible[LINK->syste - 1][i - 1][j - 1]) {
	  LINK->LINK->hap1[LINK->syste - 1] = j;
	  LINK->newval = LINK->val;
	  getval(LINK->syste, i, j, &LINK->newval, LINK->LINK);
	  if (LINK->newval != 0.0) {
	    if (LINK->syste != mlocus)
	      ugetgene(LINK->syste + 1, LINK->newval, LINK->LINK);
	    else
	      setval(LINK->newval, LINK->LINK);
	  }
	}
      }
    }
  }
}  /* quanval_ */
#endif

/* 
   This procedure finds the possible allele combinations for the
   current person at the current locus.  The person is unknown and the
   locus is an affection locus.  The run is not sexlinked and/or the
   person is female.

   Modified by Dylan in late 1994 to use 'unknown_poss'.  When LOOPSPEED
   is defined this is also called when the locus is quantified.
*/
static void aff_(LINK)
struct LOC_ugetgene *LINK;
{
/*   
   aff_ comment:  Used with an affection status phenotype or when
   thislocus[syste]^which is nul
   
   quanval_ comment: Uses this only when thislocus[syste]^.which is 
   not null--quanval_ comment
*/

  int i, j;
  thisperson *WITH;
#if LOOPSPEED
  int duplicates;
  int loop_vect_to_check;
#else
  information *WITH1;
#endif
  locusvalues *WITH2;
  int FORLIM;
  int locus;

  WITH = LINK->LINK->p;
  locus = LINK->syste - 1;

#if LOOPSPEED
  /* if known, possibles held in loop vector 0 */
  if (person[WITH->id]->unknown)
      loop_vect_to_check = single_locus_vector_num[locus];
  else
      loop_vect_to_check = 0;
#else
  WITH1 = LINK->LINK->p->store;
#endif

  WITH2 = thislocus[locus];

  FORLIM = WITH2->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap1[locus] = i;
    for (j = i; j <= FORLIM; j++) {
      if (!nohom[locus] || i != affall || j != affall) {
#if LOOPSPEED
	duplicates = (i - 1) * i / 2;
	if (unknown_poss[WITH->id][locus][loop_vect_to_check]
	    [ (i - 1) * thislocus[locus]->nallele - duplicates + (j - 1)] ) {

#else

	if (WITH1->possible[locus][i - 1][j - 1]) {
#endif
	  LINK->LINK->hap2[locus] = j;
	  LINK->newval = LINK->val;
	  getval(LINK->syste, i, j, &LINK->newval, LINK->LINK);
	  if (LINK->newval != 0.0) {
	    if (LINK->syste != mlocus)
	      ugetgene(locus + 2, LINK->newval, LINK->LINK);
	    else
	      setval(LINK->newval, LINK->LINK);
	  }
	}
      }
    }
  }

  if ((!disequi) && (!disfreqs))
    return;

  WITH = LINK->LINK->p;
  locus = LINK->syste - 1;

#if LOOPSPEED
  /* if known, possibles held in loop vector 0 */
  if ( person[WITH->id]->unknown)
      loop_vect_to_check= single_locus_vector_num[locus];
  else {
    loop_vect_to_check = 0;
  }
#else
  WITH1 = LINK->LINK->p->store;
#endif

  WITH2 = thislocus[locus];
  FORLIM = WITH2->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap2[locus] = i;
    for (j = i; j <= FORLIM; j++) {
      if (!nohom[locus] || i != affall || j != affall) {
#if LOOPSPEED
	duplicates = (i - 1) * i / 2;
	if (unknown_poss[LINK->LINK->p->id][locus][loop_vect_to_check]
	    [ (i - 1) * thislocus[locus]->nallele - duplicates + (j - 1)] ) {
#else
	if (WITH1->possible[locus][i - 1][j - 1]) {
#endif
	  LINK->LINK->hap1[locus] = j;
	  LINK->newval = LINK->val;
	  getval(locus + 1, i, j, &LINK->newval, LINK->LINK);
	  if (LINK->newval != 0.0) {
	    if (LINK->syste != mlocus)
	      ugetgene(locus + 2, LINK->newval, LINK->LINK);
	    else
	      setval(LINK->newval, LINK->LINK);
	  }
	}
      }
    }
  }
}  /* aff_ */


/* 
   This procedure determines the genotype of an unknown individual.

   Modified by Dylan in late 1994 to work with LOOPSPEED defined.
*/
static void ugetgene(syste_, val_, LINK)
int syste_;
double val_;
struct LOC_getvect *LINK;
{
  struct LOC_ugetgene V;
  locusvalues *WITH;


  V.LINK = LINK;
  V.syste = syste_;
  V.val = val_;
  WITH = thislocus[V.syste - 1];
  if (sexlink && LINK->p->male) {
    switch (WITH->which) {

    case binary_:
      facmale_(&V);
      break;

    case affection:
      affmale_(&V);
      break;

    case quantitative:
      quanmale_(&V);
      break;

    case null_:
      if (WITH->privlocus->which == affection)
	affmale_(&V);
      else
	quanmale_(&V);
      break;
    }
    return;
  }
  switch (WITH->which) {

  case binary_:
    fac_(&V);
    break;

  case affection:
    aff_(&V);
    break;

  case quantitative:
#if LOOPSPEED       /* Dylan -- aff_ and quanval_ are identical so I */
    aff_(&V);          /*          only modified aff_                   */
#else
    quanval_(&V);
#endif
    break;

  case null_:
    aff_(&V);
    break;
  }
}  /* ugetgene */


void getvect(p_, LINK)
thisperson *p_;
struct LOC_likelihood *LINK;
{
  struct LOC_getvect V;


  V.LINK = LINK;
  V.p = p_;
  if (V.p->unknown)
    ugetgene(1L, 1.0, &V);
  else
    getgene(1L, 1.0, &V);
}  /* getvect */

/*getvect*/


