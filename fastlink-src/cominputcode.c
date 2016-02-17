/* This file contains input code shared by modified versions of the */
/* ILINK,  LINKMAP, and MLINK programs */
/* The modifications are described in the papers: */
/* R. W. Cottingham, Jr., R. M. Idury, and A. A. Schaffer, */
/* Faster Sequential Genetic Linkage Computations */
/* American Journal of Human Genetics, 53(1993), pp. 252-263 */
/* and A. A. Schaffer, S. K. Gupta, K. Shriram, and R. W. Cottingham, Jr. */
/* Avoiding Recomputation in Linkage Analysis, */
/* Human Heredity 44(1994), pp. 225-237.  */

/*
   Modified by Dylan in late 1994.  Main difference is reads from
   'loopfile' instead of 'speedfile' when LOOPSPEED is set.

   Please read README.loopfile.  See also unknown.c and commondefs.c.
*/


#include "commondefs.h"
#if !defined(DOS)
#include "checkpointdefs.h"
#endif
#include "gemdefs.h"

/* Local variables for getlocations: */
struct LOC_getlocations {
  long ngene, nseg, here, there, start, nhet, thisseg;
  boolean rarepresent, riskhom, riskhet;
  hapvector hap1, hap2;
  boolean thishet[maxlocus];
};

Local boolean checkrare(LINK)
struct LOC_getlocations *LINK;
{
  long i;
  boolean check;
  locusvalues *WITH;

  check = false;
  for (i = 1; i <= mlocus; i++) {
    if (nohom[i - 1]) {
      WITH = thislocus[i - 1];
      if ((WITH->freq[LINK->hap1[i] - 1]< minfreq) ||
	  (WITH->freq[LINK->hap2[i] - 1] < minfreq))
	check = true;
    }
  }
  return check;
}

Local Void checkrisk(riskhet, riskhom, LINK)
boolean *riskhet, *riskhom;
struct LOC_getlocations *LINK;
{
  *riskhet = false;
  *riskhom = false;
  if (LINK->hap1[risksys - 1] == riskall &&
      LINK->hap2[risksys - 1] == riskall)
    *riskhom = true;
  else {
    if ((LINK->hap1[risksys - 1] != riskall &&
	 LINK->hap2[risksys - 1] == riskall) ||
	(LINK->hap2[risksys - 1] != riskall &&
	 LINK->hap1[risksys - 1] == riskall))
      *riskhet = true;
  }
}

Local long gethapn(hap, LINK)
hapvector hap;
struct LOC_getlocations *LINK;
{
  long i, n;

  n = 1;
  for (i = 1; i <= mlocus; i++)
    n += increment[i - 1] * (hap[i - 1] - 1);
  return n;
}

/* Local variables for domalerisk: */
struct LOC_domalerisk {
  struct LOC_getlocations *LINK;
} ;

Local Void setrisk(LINK)
struct LOC_domalerisk *LINK;
{
  long n;

  n = gethapn(LINK->LINK->hap1, LINK->LINK);
  if (LINK->LINK->hap1[risksys - 1] == riskall)
    riskmale[n - 1] = true;
  else
    riskmale[n - 1] = false;
}

Local Void getriskhap(system, LINK)
long system;
struct LOC_domalerisk *LINK;
{
  long i;
  locusvalues *WITH;
  long FORLIM;

  WITH = thislocus[system - 1];
  FORLIM = WITH->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap1[system - 1] = i;
    if (system != mlocus)
      getriskhap(system + 1, LINK);
    else
      setrisk(LINK);
  }
}

Local Void domalerisk(LINK)
struct LOC_getlocations *LINK;
{
  struct LOC_domalerisk V;

  V.LINK = LINK;
  getriskhap(1L, &V);
}

/* Local variables for domutation: */
struct LOC_domutation {
  struct LOC_getlocations *LINK;
} ;

Local Void setmutation(LINK)
struct LOC_domutation *LINK;
{
  long i, n;

  n = gethapn(LINK->LINK->hap1, LINK->LINK);
  if (LINK->LINK->hap1[mutsys - 1] == thislocus[mutsys - 1]->nallele) {
    muthap[n - 1] = n;
    return;
  }
  i = LINK->LINK->hap1[mutsys - 1];
  LINK->LINK->hap1[mutsys - 1] = thislocus[mutsys - 1]->nallele;
  muthap[n - 1] = gethapn(LINK->LINK->hap1, LINK->LINK);
  LINK->LINK->hap1[mutsys - 1] = i;
}

Local Void getmuthap(system, LINK)
long system;
struct LOC_domutation *LINK;
{
  long i;
  locusvalues *WITH;
  long FORLIM;

  WITH = thislocus[system - 1];
  FORLIM = WITH->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap1[system - 1] = i;
    if (system != mlocus)
      getmuthap(system + 1, LINK);
    else
      setmutation(LINK);
  }
}

Local Void domutation(LINK)
struct LOC_getlocations *LINK;
{
  struct LOC_domutation V;

  V.LINK = LINK;
  getmuthap(1L, &V);
}

Local Void setnumbers(LINK)
struct LOC_getlocations *LINK;
{
  long nhap1, nhap2;

  LINK->ngene++;

  segstart[LINK->ngene - 1] = LINK->here + 1;
  probstart[LINK->ngene - 1] = LINK->there + 1;
  probend[LINK->ngene - 1] = LINK->there + LINK->nseg;

  LINK->there += LINK->nseg;

  nhap1 = gethapn(LINK->hap1, LINK);
  nhap2 = gethapn(LINK->hap2, LINK);
  genenumber[nhap1 - 1][nhap2 - 1] = LINK->ngene;
  genenumber[nhap2 - 1][nhap1 - 1] = LINK->ngene;

  if (minfreq != 0.0) {
    if (LINK->rarepresent)
      rare[LINK->ngene - 1] = true;
    else
      rare[LINK->ngene - 1] = false;
  } else
    rare[LINK->ngene - 1] = false;
  if (risk) {
    risk1[LINK->ngene - 1] = LINK->riskhet;
    risk2[LINK->ngene - 1] = LINK->riskhom;
  }

  LINK->thisseg++;
  invgenenum1[LINK->thisseg - 1] = nhap1;
  invgenenum2[LINK->thisseg - 1] = nhap2;
}

Local Void hapscr(system, nscramble, LINK)
long system, nscramble;
struct LOC_getlocations *LINK;
{
  long i, j;

  if (LINK->thishet[system - 1])
    nscramble++;
  if (system != mlocus)
    hapscr(system + 1, nscramble, LINK);
  else
    setnumbers(LINK);
  if (nscramble <= 1)
    return;
  if (LINK->hap1[system - 1] == LINK->hap2[system -1])
    return;
  i = LINK->hap1[system - 1];
  j = LINK->hap2[system - 1];
  LINK->hap1[system - 1] = j;
  LINK->hap2[system - 1] = i;
  if (system != mlocus)
    hapscr(system + 1, nscramble, LINK);
  else
    setnumbers(LINK);
  LINK->hap1[system - 1] = i;
  LINK->hap2[system - 1] = j;
}

Local Void sethap(system, LINK)
long system;
struct LOC_getlocations *LINK;
{
  long i, j;
  locusvalues *WITH;
  long FORLIM, FORLIM1;

  WITH = thislocus[system - 1];
  if (LINK->thishet[system - 1]) {
    FORLIM = WITH->nallele;
    for (i = 1; i < FORLIM; i++) {
      LINK->hap1[system - 1] = i;
      FORLIM1 = WITH->nallele;
      for (j = i + 1; j <= FORLIM1; j++) {
	LINK->hap2[system - 1] = j;
	if (system != mlocus)
	  sethap(system + 1, LINK);
	else {
	  LINK->rarepresent = checkrare(LINK);
	  if (risk)
	    checkrisk(&LINK->riskhet, &LINK->riskhom, LINK);
	  LINK->there = LINK->start;
	  LINK->thisseg = LINK->here;
	  hapscr(1L, 0L, LINK);
	  LINK->here += LINK->nseg;
	}
      }
    }
    return;
  }
  FORLIM = WITH->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->hap1[system - 1] = i;
    LINK->hap2[system - 1] = i;
    if (system != mlocus)
      sethap(system + 1, LINK);
    else {
      LINK->rarepresent = checkrare(LINK);
      if (risk)
	checkrisk(&LINK->riskhet, &LINK->riskhom, LINK);
      LINK->thisseg = LINK->here;
      LINK->there = LINK->start;
      hapscr(1L, 0L, LINK);
      LINK->here += LINK->nseg;
    }
  }
}

Local Void starthap(LINK)
struct LOC_getlocations *LINK;
{
  long i, FORLIM;

  LINK->nseg = 1;
  FORLIM = LINK->nhet;
  for (i = 2; i <= FORLIM; i++)
    LINK->nseg *= 2;
  sethap(1L, LINK);
  LINK->start = LINK->there;
}

Local Void gethet1(system, LINK)
long system;
struct LOC_getlocations *LINK;
{
  LINK->thishet[system - 1] = false;
  if (system != mlocus)
    gethet1(system + 1, LINK);
  else
    starthap(LINK);
  LINK->thishet[system - 1] = true;
  LINK->nhet++;
  if (system != mlocus)
    gethet1(system + 1, LINK);
  else
    starthap(LINK);
  LINK->nhet--;
}


Void getlocations()
{
  struct LOC_getlocations V;

  V.nhet = 0;
  V.here = 0;
  V.there = 0;
  V.ngene = 0;
  V.start = 0;
  gethet1(1L, &V);
  if (mutsys != 0)
    domutation(&V);
  if (sexlink && risk)
    domalerisk(&V);
}


Void inputerror(nerror, par1, par2, LINK)
long nerror, par1, par2;
struct LOC_inputdata *LINK;
{
  printf("Fatal error detected in procedure inputdata\n");
  switch (nerror) {

  case 0:
    printf("Number of loci %2ld exceeds the constant maxlocus\n", par1);
    break;

  case 1:
    printf("Number of loci read %2ld. Less than minimum of 1\n", par1);
    break;

  case 2:
    printf(
      "Error detected reading loci order. Locus number %2ld in position %2ld exceeds number of loci\n",
      par2, par1);
    break;

  case 3:
    printf(
      "Error detected reading loci order. Illegal locus number %2ld in position %2ld\n",
      par2, par1);
    break;

  case 4:
    printf(
      "Error detected reading loci order. Locus number repeated in positions %2ld and %2ld\n",
      par1, par2);
    break;

  case 5:
    printf(
      "Error detected reading locus description. Illegal locus type %2ld for locus %2ld\n",
      par2, par1);
    break;

  case 6:
    printf(
      "Error detected reading locus description for system %2ld. Number of alleles  %2ld exceeds maxall\n",
      par1, par1);
    break;

  case 7:
    printf(
      "Error detected reading locus description for system %2ld. Illegal number of alleles  %2ld\n",
      par1, par2);
    break;

  case 8:
    printf(
      "Error detected reading locus description for system %2ld.\n Number of factors  %2ld exceeds maxfact or length of a long int\n",
      par1, par2);
    break;

  case 9:
    printf(
      "Error detected reading locus description for system %2ld. Illegal number of factors  %2ld\n",
      par1, par2);
    break;

  case 10:
    printf(
      "Error detected reading locus description for system %2ld. Alleles not codominant\n",
      par1);
    break;

  case 11:
    printf("Error detected reading pedigree record %2ld. Illegal code for sex %2ld\n",
	   par1, par2);
    break;

  case 12:
    printf(
      "Error detected reading pedigree record at pedigree%2ld. Maximum number of pedigree records exceeded\n",
      par1);
    break;

  case 13:
    printf(
      "Error detected reading pedigree record %2ld. Maximum number of individuals exceeded\n",
      par1);
    break;

  case 14:
    printf(
      "Error detected reading pedigree record %2ld. Illegal binary factor code %2ld\n",
      par1, par2);
    break;

  case 15:
    printf(
      "Error detected reading pedigree record %2ld. No allelic pair for genotype\n",
      par1);
    break;

  case 16:
    printf(
      "Error detected reading pedigree record %2ld. Allele number %2ld exceeds maxall\n",
      par1, par2);
    break;

  case 17:
    printf(
      "Error detected reading pedigree record %2ld. Illegal allele number %2ld\n",
      par1, par2);
    break;

  case 18:
    printf("Number of systems after factorization (%3ld) exceeds maxsystem\n",
	   par1);
    break;

  case 19:
    printf("Number of systems after factorization (%3ld) less than minimum of 1\n",
	   par1);
    break;

  case 20:
    printf("Number of recombination types (%3ld) exceeds maxrectype\n", par1);
    break;

  case 21:
    printf("Number of recombination types (%3ld) less than minimum of 1\n",
	   par1);
    break;

  case 22:
    printf(
      "End of file detected in tempdat by procedure readthg before all data found\n");
    break;

  case 23:
    printf(
      "Error detected reading iterated locus in datafile. Value (%3ld) greater than nlocus\n",
      par1);
    break;

  case 24:
    printf(
      "Error detected reading iterated locus in datafile. Illegal value (%3ld)\n",
      par1);
    break;

  case 25:
    printf("Number of iterated parameters greater then maxn in gemdefs.h\n");
    break;

  case 26:
    printf(
      "Error detected reading pedigree record %2ld. Liability class (%2ld) exceeds nclass\n",
      par1, par2);
    break;

  case 27:
    printf(
      "Error detected reading pedigree record %2ld. Illegal liability class (%2ld)\n",
      par1, par2);
    break;

  case 28:
    printf(
      "Error detected reading locus description for system%2ld. Liability classes (%3ld) exceed maxliab\n",
      par1, par2);
    break;

  case 29:
    printf(
      "Error detected reading locus description for system%2ld. Illegal number of liability classes (%3ld)\n",
      par1, par2);
    break;

  case 30:
    printf(
      "Error detected reading locus description for system%2ld. Penetrance out of range\n",
      par1);
    break;

  case 31:
    printf(
      "Error detected reading locus description for system%2ld. Number of traits (%3ld) exceeds maxtrait\n",
      par1, par2);
    break;

  case 32:
    printf(
      "Error detected reading locus description for system%2ld. Number of traits out of range (%3ld)\n",
      par1, par2);
    break;

  case 33:
    printf(
      "Error detected reading locus description for system%2ld. Variance must be positive\n",
      par1);
    break;

  case 34:
    printf(
      "Error detected reading locus description for system%2ld. Variance multiplier must be positive\n",
      par1);
    break;

  case 35:
    printf(
      "Error detected reading locus description for system%2ld. Risk allele %3ld) exceeds nallele\n",
      par1, par2);
    break;

  case 36:
    printf(
      "Error detected reading locus description for system%2ld. Illegal risk allele (%3ld)\n",
      par1, par2);
    break;

  case 37:
    printf("Error detected reading datafile. Risk locus %3ld) exceeds nlocus\n",
	   par2);
    break;

  case 38:
    printf("Error detected reading datafile. Illegal value for risk locus %3ld)\n",
	   par2);
    break;

  case 39:
    printf("Error detected reading datafile. Mutation locus %3ld) exceeds nlocus\n",
	   par2);
    break;

  case 40:
    printf(
      "Error detected reading datafile. Illegal value for mutation locus %3ld)\n",
      par2);
    break;

  case 41:
    printf(
      "Error detected reading datafile. Linkage disequilbirium is not allowed with this program\n");
    break;

  case 42:
    printf("Locus %5ld in lod score list exceeds nlocus %5ld\n", par1, par2);
    break;

  case 43:
    printf("Illegal locus number %5ld in lod score list\n", par1);
    break;
  }
  exit(EXIT_FAILURE);
}

Void inputwarning(nwarning, par1, par2, LINK)
long nwarning, par1, par2;
struct LOC_inputdata *LINK;
{
  printf("Warning number from procedure inputdata\n");
  switch (nwarning) {

  case 0:
    printf("Illegal sex difference parameter %2ld Parameter should be 0, 1, or 2\n",
	   par1);
    break;

  case 1:
    printf("Illegal interference parameter %2ld Lack of interference assumed\n",
	   par1);
    break;

  case 2:
    printf(
      "Illegal sex difference parameter %2ld Parameter must be 0 with sex-linked data\n",
      par1);
    break;

  case 3:
    printf(
      "Non-standard affection status%4ld interpreted as normal in pedigree record%5ld\n",
      par2, par1);
    break;
  }
}

#if LOOPSPEED
/*
   Written by Dylan in late 1994.
   Thus routine opens up loopfile.dat
*/
void open_loop_file() 
{
  if (loopfile != NULL) {
    fclose(loopfile);
    loopfile = fopen(LOOPFILE_NAME, "r");
  } else {
    loopfile = fopen(LOOPFILE_NAME, "r");
  }
  if (loopfile == NULL) {
    printf("The file %s was not found.\n", LOOPFILE_NAME);
    printf("Check that LOOPSPEED is defined in unknown.c\n");
    exit(FileNotFound);
  }
}

/*
   Written by Dylan in late 1994.
   This routine closes loopfile.dat
*/
void close_loop_file()
{
  if (loopfile != NULL) {
    fclose(loopfile);
  }
  loopfile = NULL;
}

/*free the memory associated with is_zero_breaker*/
static void free_is_zero_breaker()
{

  int locus;
  for (locus = 0; locus < mlocus; locus++) 
    free(is_zero_breaker[locus]);
  free(is_zero_breaker);
  is_zero_breaker = NULL;
}

/*free the memory associated with breaker_poss_genotype*/
static void free_breaker_poss_genotype()
{
  int loopindex; /*index over loop breakers*/

  for(loopindex = 0; loopindex < maxloop; loopindex++) {
    if (NULL !=  breaker_poss_genotype[loopindex]) {
      free(breaker_poss_genotype[loopindex]);
      breaker_poss_genotype[loopindex] = NULL;
    }
  }
}

/*
   This procedure frees the space used by 'loop_vectors'.  Written by Dylan.
*/
static void free_loop_vectors () 
{
  int locus, vect;

  if (loop_vectors != NULL) {
    for (locus = 0; locus < mlocus; locus++) {
      if ( loop_vectors[locus] != NULL ) {
	for (vect = 0; vect < num_loop_vectors[locus]; vect++) {
	  if ( loop_vectors[locus][vect] != NULL ) {
	    free(loop_vectors[locus][vect]);
	  }
	}
	free(loop_vectors[locus]);
      }
    }
    free(loop_vectors);
    loop_vectors = NULL;
  }
} /* free_loop_vectors */

/* 
   This procedure frees the space used by the data structure
   'unknown_poss'.  Written by Dylan in late 1994.
*/
static void free_unknown_poss() {
  int curr_person, locus, loop_v;

  if (unknown_poss != NULL) {
    for (curr_person = 1; curr_person <= totperson; curr_person++){
      /* only entries for people with unknown genotype */
       if ( person[curr_person]->unknown ) {
	for (locus = 0; locus < mlocus; locus++){
	  if ( person[curr_person]->thisunknown[locus] ) {
	    for (loop_v = 0; loop_v < num_loop_vectors[locus]; loop_v++) {
	      free(unknown_poss[curr_person][locus][loop_v]);
	    } /* for each loop vector */
	  } else {  /* known at this locus so only one loop vector used */
	    free(unknown_poss[curr_person][locus][0]);
	  }
	  free(unknown_poss[curr_person][locus]);
	}  /* for each locus */
	free(unknown_poss[curr_person]);
	person[curr_person]->unknown = false;
      }
    } /* for each person */
    free(unknown_poss);
    unknown_poss = NULL;
  }
} /* free_unknown_poss */

/*
   Written by Dylan in late 1994.
*/
static void read_loop_err() 
{
/*  other errors get written to standard out for some reason
  fprintf(stderr, "Error in reading %s.\n", LOOPFILE_NAME);
*/
  printf("Error in reading %s.\n", LOOPFILE_NAME);

  exit(EXIT_FAILURE);
}

static void malloc_loop_err(message, fewer_vects_size)
char * message;
int fewer_vects_size;
{

  FILE *errorfile;
  time_t secondsNow;

  fprintf(stderr, "\nProblem with malloc, probably not enough space\n");
  fprintf(stderr, "Problem occurred when allocating %s\n", message);
  fprintf(stderr, "while reading %s.\n", LOOPFILE_NAME); 
  fprintf(stderr, "Try reducing max_vectors_considered to %d.\n", fewer_vects_size);

  errorfile = fopen("FASTLINK.err","a");
  if (errorfile) {
    time (&secondsNow);
    fprintf(errorfile,"\n%s",ctime(&secondsNow));
    fprintf(errorfile, "\nProblem with malloc, probably not enough space\n");
    fprintf(errorfile, "Problem occurred when allocating %s\n", message);
    fclose(errorfile);
  }
  exit(EXIT_FAILURE);
}

/*infer possible genotypes for each loop breaker*/

void infer_genotypes(curr_ped) 
int curr_ped;
{
  int loopindex; /*index over loop breakers*/
  boolean **breaker_poss_bylocus[maxloop]; /*possible single locus genotypes
                                              for each loop breaker and locus*/
  int plocus; /*placeholder for locus index*/
  int locus; /*loop index for loci*/
  int genoindex; /*loop index for single locus and multilocus genotypes*/
  int hap1, hap2; /*haplotypes for one person*/
  int num_alleles; /*number of alleles at a locus*/
  int vect; /*loop index for single locus loop vectors over all breakers*/
  int a,b,temp; /*placeholders for alleles*/
  int this_locus_geno; /*genotype for a single person at a single locus*/
  int duplicates; /*used to go from genotype to alleles*/


  for(loopindex = 0; loopindex < num_loops[curr_ped - 1]; loopindex++) {
    breaker_poss_bylocus[loopindex] = (boolean **) malloc (mlocus * sizeof (boolean *));
    if (NULL == breaker_poss_bylocus[loopindex])
      malloc_err("breaker_poss_bylocus middle level");
    for (locus = 0; locus < mlocus; locus++) {
      plocus = order[locus] - 1;
      breaker_poss_bylocus[loopindex][plocus] = (boolean *) malloc (thislocus[plocus]->fgeno);
    if (NULL == breaker_poss_bylocus[loopindex][plocus])
      malloc_err("breaker_poss_bylocus inner level");
      for (genoindex = 0; genoindex < thislocus[plocus]->fgeno; genoindex++)
        breaker_poss_bylocus[loopindex][plocus][genoindex] = false;
    }
  }
  for (locus = 0; locus < mlocus; locus++) {
    plocus = order[locus] - 1;
    for(vect = 0; vect < num_loop_vectors[plocus]; vect++)
      if (!(is_zero_breaker[plocus][vect])) {
	for (loopindex = 0; loopindex < num_loops[curr_ped - 1]; loopindex++)
	  breaker_poss_bylocus[loopindex][plocus]
             [loop_vectors[plocus][vect][loopindex]] = true;
      }
  }
  for(loopindex = 0; loopindex < maxloop; loopindex++)
    breaker_poss_genotype[loopindex] = NULL;
  for(loopindex = 0; loopindex < num_loops[curr_ped - 1]; loopindex++) {
    breaker_poss_genotype[loopindex] = (boolean *) malloc (fgeno * sizeof(boolean));
    if (NULL == breaker_poss_genotype[loopindex]) {
      malloc_err("breaker_poss_genotype");
    }
    for (genoindex = 0; genoindex < fgeno; genoindex++) {
      breaker_poss_genotype[loopindex][genoindex] = true;
      if (sexlink && (looppers[curr_ped - 1][loopindex][0]->male)) {
	hap1 = genoindex;             /*already decremented*/
	hap2 = genoindex;             
      }
      else {
	hap1 = invgenenum1[genoindex] - 1;
	hap2 = invgenenum2[genoindex] - 1;
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
	if (sexlink && (looppers[curr_ped -1][loopindex][0]->male))
	  this_locus_geno = a;
	else
	  this_locus_geno = (a * num_alleles) - duplicates + b;
        if (!(breaker_poss_bylocus[loopindex][locus][this_locus_geno])) {
	  breaker_poss_genotype[loopindex][genoindex] = false;
          break;
	}
      }
    }
  }
  for(loopindex = 0; loopindex < num_loops[curr_ped - 1]; loopindex++) {
    if (NULL != breaker_poss_bylocus[loopindex]) {
      for (locus = 0; locus < mlocus; locus++) 
	free(breaker_poss_bylocus[loopindex][locus]);
      free(breaker_poss_bylocus[loopindex]);
      breaker_poss_bylocus[loopindex] = NULL;
    }
  }
}

/*
   This procedure reads the next pedigree in 'loopfile'.  The information
   is read into 'loop_vectors', 'num_loop_vectors', and 'unknown_poss'.

   Written by Dylan in late 1994.
*/
void read_loop_file(curr_ped) 
int curr_ped;
{
  int ped_read;
  char ch;
  int locus, vect, geno;
  int plocus; /*A. A. Schaffer*/
  int check_locus, check_vect;
  int ind, num_vects;
  int loop;
  int init;
  int fewer_vects_size;
  int b; /*index over loop_breaker vectors*/
  boolean any_possible; /*count of possible genotypes*/
  char extrachar; /*for debugging*/

  /* if only one pedigree, don't reread each time */
  if ( ever_read_loopfile && (nuped == 1) ) {
    return;
  }
  
  /* free data from last pedigree */
  if ( ever_read_loopfile ) {
    free_unknown_poss();
    free_loop_vectors();
    free_is_zero_breaker();
    free_breaker_poss_genotype();
  }

  /* read in pedigree -- assumes pedigrees start at 1 */
  if ( curr_ped == 1 ) {
    if ( !fscanf(loopfile, "Pedigree: %d\n", &ped_read) ) {
      read_loop_err();
    }
  } else {
    if ( !fscanf(loopfile, "edigree: %d\n", &ped_read) ) {
      read_loop_err();
    }
  }

  /* read fewer_vects_size */
  if ( !fscanf(loopfile, "fewer_vects_size: %d\n", &fewer_vects_size) ) {
    read_loop_err();
  }

  /* read in num_loops_considered */
  if (!fscanf(loopfile,"num_loops_considered: %d\n",&(num_loops[ped_read-1]))){
    read_loop_err();
  }

  /* read in num_loop_vectors */
  if ( !fscanf(loopfile, "num_loop_vectors %c\n", &ch) ) {
    read_loop_err();
  }
  for (locus = 0; locus < mlocus; locus++) {
    if ( !fscanf(loopfile, "\t%d :", &check_locus) || (locus != check_locus)){
      read_loop_err();
    }
    fscanf(loopfile, " %d\n", &num_vects);
    num_loop_vectors[order[locus] - 1] = num_vects;
  }

  /*Alex*/
  is_zero_breaker = (boolean **) malloc(mlocus * sizeof(boolean *));
  if (is_zero_breaker == NULL)
    malloc_err("is_zero_breaker");
  for (locus = 0; locus < mlocus; locus++) {
    
    is_zero_breaker[locus] = (boolean *) malloc(num_loop_vectors[locus] * sizeof(boolean));
    if (is_zero_breaker[locus] == NULL)
      malloc_loop_err("is_zero_breaker", locus);
    for(b=0; b < num_loop_vectors[locus]; b++)
      is_zero_breaker[locus][b] = false;
  }

  
  /* read in loop vectors array */
  if ( !fscanf(loopfile, "loop_vectors %c\n", &ch) ) {
    read_loop_err();
  }
  loop_vectors =(vector_for_locus *) malloc(mlocus * sizeof(vector_for_locus));
  if ( loop_vectors == NULL ) {
    malloc_loop_err("loop vectors array", fewer_vects_size);
  }
  for (locus = 0; locus < mlocus; locus++) {
    if (!fscanf(loopfile, "\tL : %d\n",&check_locus) || (locus !=check_locus)){
      read_loop_err();
    }
    plocus = order[locus] -1;  /*Alex*/
    loop_vectors[plocus] = (vector_for_loop *)
      malloc(num_loop_vectors[plocus] * sizeof(vector_for_loop));
    if ( loop_vectors[plocus] == NULL ) {
      malloc_loop_err("loop vectors array", fewer_vects_size);
    }
    for (vect = 0; vect < num_loop_vectors[plocus]; vect++) {
      if ( !fscanf(loopfile, "\t\t%d :", &check_vect) || (vect != check_vect)){
	read_loop_err();
      }
      loop_vectors[plocus][vect] = (int *) 
	malloc (num_loops[curr_ped - 1] * sizeof(int));
      /* DYLAN -- warning EOL */
      loop = 0;
      while ( getc(loopfile) != '\n' ) {
	fscanf(loopfile, "%d", &geno);
	loop_vectors[plocus][vect][loop] = geno;
	loop++;
      }
      if ( loop != num_loops[curr_ped - 1] ) {
	read_loop_err();
      }
    }
  }

  /* read in unknown_poss */
  if ( !fscanf(loopfile, "unknown_poss %c\n", &ch) ) {
    read_loop_err();
  }
  unknown_poss = (geno_for_unknown *) 
    malloc((totperson + 1) * sizeof(geno_for_unknown));
  if ( unknown_poss == NULL ) {
    malloc_loop_err("table of possible genotypes for unknowns", fewer_vects_size);
  }
  while ((extrachar =  getc(loopfile)) == 'i' ) {
    if ( !fscanf(loopfile, "d: %d\n", &ind) ) {
      read_loop_err();
    }
    person[ind]->unknown = true;
    unknown_poss[ind] = (geno_for_locus *) 
      malloc (mlocus * sizeof(geno_for_locus));
    if ( unknown_poss[ind] == NULL ) { 
      malloc_loop_err("table of possible genotypes for unknowns", fewer_vects_size);
    }
    for (locus = 0; locus < mlocus; locus++) {
      if ( !fscanf(loopfile, "\tL: %d\n", &check_locus) || 
	  (locus != check_locus) ) {
	read_loop_err();
      }
      plocus = order[locus] - 1;
      person[ind]->thisunknown[plocus] = true;
      unknown_poss[ind][plocus] = (geno_for_loop_vector *)
	malloc(num_loop_vectors[plocus] * sizeof(geno_for_loop_vector));
      if ( unknown_poss[ind][plocus] == NULL )  {
	malloc_loop_err("table of possible genotypes for unknowns", fewer_vects_size);
      }
      if ( getc(loopfile) == '-' ) {
	person[ind]->thisunknown[plocus] = false;
      }
      for (vect = 0; vect < num_loop_vectors[plocus]; vect++) {
	if ( !fscanf(loopfile, "\t\t%d :", &check_vect) ||
	    (vect != check_vect)) {
	  read_loop_err();
	}
	unknown_poss[ind][plocus][vect] = (boolean *)
	  malloc (thislocus[plocus]->fgeno * sizeof(boolean));
	if ( unknown_poss[ind][plocus][vect] == NULL ) {
	  malloc_loop_err("table of possible genotypes for unknowns", fewer_vects_size);
	}
	for (init = 0; init < thislocus[plocus]->fgeno; init++) {
	  unknown_poss[ind][plocus][vect][init] = false;
	}
        any_possible = false;
	while ( getc(loopfile) != '\n' ) {
	  fscanf(loopfile, "%d", &geno);
	  unknown_poss[ind][plocus][vect][geno] = true;
          any_possible = true;
	}
        if (!any_possible)
          is_zero_breaker[plocus][vect] = true;
      }
    }
  }
  ever_read_loopfile = true;  /* global flag to say we're read this */
}  /* read_loop_file */



static void readsignature()
{

  int allele_speed_signature;  /*used to verify that UNKNOWN also has
                                 ALLELE_SPEED on*/
  int signature;

  signature = 0;
  allele_speed_signature = 0;
  signature = fscanf(speedfile,"%d", &allele_speed_signature);
  if ( ((ALLELE_SPEED &&
        (( signature < 1) || (allele_speed_signature != ALLELE_SPEED_CONSTANT))) ||
       ((!ALLELE_SPEED) && (signature > 0)))) {
    fprintf(stderr,"\nYou appear to be using a version of UNKNOWN");
    fprintf(stderr,"\nthat predates FASTLINK 3.0, or has ALLELE_SPEED set");
    fprintf(stderr,"\ndifferently. ALLELE_SPEED must be set to the same value");
    fprintf(stderr,"\nin unknown.c and commondefs.c.");
    fprintf(stderr,"\nFASTLINK will exit politely to allow you to fix the problem.");
    fprintf(stderr,"\nPlease ensure that your path gives you compatible versions");
    fprintf(stderr,"\nof UNKNOWN and the main programs.\n");
    exit(EXIT_FAILURE);
  }
  return;
}

#else

Void readspeed(LINK)
struct LOC_inputdata *LINK;
{
  long i, a, b, sys;
  char ch;
  information *WITH;
  int allele_speed_signature;  /*used to verify that UNKNOWN also has
                                 ALLELE_SPEED on*/
  int signature;

  signature = 0;
  allele_speed_signature = 0;
  signature = fscanf(speedfile,"%d\n", &allele_speed_signature);
  if ( ((ALLELE_SPEED &&
        ((signature < 1) || (allele_speed_signature != ALLELE_SPEED_CONSTANT))) ||
       ((!ALLELE_SPEED) && (signature > 0)))) {
    fprintf(stderr,"\nYou appear to be using a version of UNKNOWN");
    fprintf(stderr,"\nthat predates FASTLINK 3.0, or has ALLELE_SPEED set");
    fprintf(stderr,"\ndifferently. ALLELE_SPEED must be set to the same value");
    fprintf(stderr,"\nin unknown.c and commondefs.c.");
    fprintf(stderr,"\nFASTLINK will exit politely to allow you to fix the problem.");
    fprintf(stderr,"\nPlease ensure that your path gives you compatible versions");
    fprintf(stderr,"\nof UNKNOWN and the main programs.\n");
    exit(EXIT_FAILURE);
  }

  while (!P_eof(speedfile)) {
    ch = getc(speedfile);
    if (ch == '\n')
      ch = ' ';
    if (ch == 'i' || ch == 'I') {
      fscanf(speedfile, "%c%ld", &ch, &i);
      if (ch == '\n')
	ch = ' ';
      person[i]->unknown = true;
      person[i]->store = (information *)Malloc(sizeof(information));
      if ((person[i]->store) == NULL)
        malloc_err("store field");
      WITH = person[i]->store;
      for (sys = 0; sys < mlocus; sys++) {
	for (a = 0; a < maxall; a++) {
	  for (b = 0; b < maxall; b++)
	    WITH->possible[sys][a][b] = false;
	}
      }
    } else {
      WITH = person[i]->store;
      fscanf(speedfile, "%ld%ld%ld", &sys, &a, &b);
      if (sys <= mlocus)
	WITH->possible[order[sys - 1] - 1][a - 1][b - 1] = true;
    }
    fscanf(speedfile, "%*[^\n]");
    getc(speedfile);
  }
}
#endif

/* Local variables for readped: */
struct LOC_readped {
  struct LOC_inputdata *LINK;
  long sequence;
  long startped[maxped], endped[maxped];
} ;

/* Local variables for getphenotype: */
struct LOC_getphenotype {
  struct LOC_readped *LINK;
  thisperson **p;
  long system;
} ;

Local Void readbin(phen, ourlocus, LINK)
phenotype **phen;
locusvalues *ourlocus;
struct LOC_getphenotype *LINK;
{
  long i, j;
  phenotype *WITH;
  long FORLIM;


  WITH = *phen;
  WITH->which = binary_;
  WITH->phenf = 0;
  FORLIM = LINK->LINK->LINK->nfactor[LINK->system - 1];
  for (i = 1; i <= FORLIM; i++) {
    fscanf(ipedfile, "%ld", &j);
    if (j != 0 && j != 1)
      inputerror(14L, (*LINK->p)->id, j, LINK->LINK->LINK);
    if (j == 1)
      WITH->phenf = ((long)WITH->phenf) | (1 << ((long)i));
  }
}

Local Void readnumber(phen, ourlocus, LINK)
phenotype **phen;
locusvalues *ourlocus;
struct LOC_getphenotype *LINK;
{
  long i, j, l;
  phenotype *WITH;

  /* dwix: begin */
  static int pedidx = 0;
  static int pednum = -1;
  int locidx, allidx;
  /* dwix: end */
  int allfirst, allsecond;

  WITH = *phen;
  WITH->which = binary_;
  WITH->phenf = 0;
  for (i = 1; i <= 2; i++) {
    fscanf(ipedfile, "%ld", &j);
    if (j > maxall)
      inputerror(16L, (*LINK->p)->id, j, LINK->LINK->LINK);
    if ((j < 0) || (j > ourlocus->nallele))
      inputerror(17L, (*LINK->p)->id, j, LINK->LINK->LINK);
    if (j != 0)
#if ALLELE_SPEED
      WITH->phenf = 1;
#else
      WITH->phenf = ((long)WITH->phenf) | (1 << ((long)j));
#endif
    if (1==i)
      allfirst = j;
    else
      allsecond = j;

    /* dwix: begin */
    allidx = j;
    locidx = -1;
    for(l=0; l < mlocus; l++)
      if(ourlocus == thislocus[l])
	locidx = l;
    if (pednum != (*LINK->p)->ped) {
      pedidx++;
      pednum = (*LINK->p)->ped;
    }
    ped_loc_all[pedidx - 1][locidx][allidx].present = true;
    /* dwix: end */
  }
  if (allfirst <= allsecond) {
   WITH->allele1 = allfirst;
   WITH->allele2 = allsecond;
  }
  else {
   WITH->allele1 = allsecond;
   WITH->allele2 = allfirst;
  }
}

Local Void readaff(phen, thislocus, LINK)
phenotype **phen;
locusvalues *thislocus;
struct LOC_getphenotype *LINK;
{
  long thisval;
  phenotype *WITH;

  WITH = *phen;
  WITH->which = affection;
  fscanf(ipedfile, "%ld", &thisval);
  if (thisval == missaff)
    WITH->aff = 0;
  else {
    if (thisval == affval)
      WITH->aff = 2;
    else {
      if (thisval != 1)
	inputwarning(3L, (*LINK->p)->id, thisval, LINK->LINK->LINK);
      WITH->aff = 1;
    }
  }
  if (thislocus->UU.U0.nclass == 1)
    WITH->liability = 1;
  else
    fscanf(ipedfile, "%ld", &(*phen)->liability);
  if (WITH->liability > thislocus->UU.U0.nclass)
    inputerror(26L, (*LINK->p)->id, WITH->liability, LINK->LINK->LINK);
  if (WITH->liability <= 0)
    inputerror(27L, (*LINK->p)->id, WITH->liability, LINK->LINK->LINK);
}

Local Void readquan(phen, thislocus, LINK)
phenotype **phen;
locusvalues *thislocus;
struct LOC_getphenotype *LINK;
{
  long i;
  double xval;
  phenotype *WITH;
  long FORLIM;

  WITH = *phen;
  if (!sexlink || !(*LINK->p)->male) {
    WITH->which = quantitative;
    FORLIM = thislocus->UU.U1.ntrait;
    for (i = 0; i < FORLIM; i++)
      fscanf(ipedfile, "%lf", &(*phen)->x[i]);
    WITH->missing = true;
    FORLIM = thislocus->UU.U1.ntrait;
    for (i = 0; i < FORLIM; i++) {
      if (WITH->x[i] != missval)
	WITH->missing = false;
    }
    return;
  }
  WITH->which = affection;
  fscanf(ipedfile, "%lf", &xval);
  if (xval == missval)
    WITH->aff = missaff;
  else {
    if (xval == affall)
      WITH->aff = affall;
    else
      WITH->aff = -11;
  }
  WITH->liability = 1;
  FORLIM = thislocus->UU.U1.ntrait;
  for (i = 2; i <= FORLIM; i++)
    fscanf(ipedfile, "%lf", &xval);
}

Local Void getphenotype(p_, LINK)
thisperson **p_;
struct LOC_readped *LINK;
{
  struct LOC_getphenotype V;
  long thisread;
  thisperson *WITH;

  V.LINK = LINK;
  V.p = p_;
  WITH = *V.p;
  for (thisread = 0; thisread < mlocus; thisread++) {
    V.system = order[thisread];
    WITH->phen[V.system - 1] = NULL;
    if (thislocus[V.system - 1]->which != null_) {
      WITH->phen[V.system - 1] = (phenotype *)Malloc(sizeof(phenotype));
      if ((WITH->phen[V.system - 1]) == NULL)
        malloc_err("phen field");
    }    
    switch (thislocus[V.system - 1]->which) {

    case quantitative:
      readquan(&WITH->phen[V.system - 1], thislocus[V.system - 1], &V);
      break;

    case affection:
      readaff(&WITH->phen[V.system - 1], thislocus[V.system - 1], &V);
      break;

    case binary_:
      if (thislocus[V.system - 1]->format == allformat)
	readnumber(&WITH->phen[V.system - 1], thislocus[V.system - 1], &V);
      else
	readbin(&WITH->phen[V.system - 1], thislocus[V.system - 1], &V);
      break;
    /* cgh -- added this for gcc warning */
    case null_:
    default:
      break;
    }
    if (lastpriv == V.system) {
      WITH->privphen = (phenotype *)Malloc(sizeof(phenotype));
      if (WITH->privphen == NULL)
        malloc_err("privphen field");
      switch (thislocus[V.system - 1]->privlocus->which) {

      case quantitative:
	readquan(&WITH->privphen, thislocus[V.system - 1]->privlocus, &V);
	break;

      case affection:
	readaff(&WITH->privphen, thislocus[V.system - 1]->privlocus, &V);
	break;
      /* cgh -- added this for gcc warning */
      case binary_:
      case null_:
      default:
	break;
      }
    }
  }
}


Local Void getinformative(LINK)
struct LOC_readped *LINK;
{
  long i, j, k, l, m, count, nchild;
  thisperson *child;
  long FORLIM1, FORLIM3;
  phenotype *WITH;
  locusvalues *WITH1;

  if (fitmodel || risk) {
    for (i = 0; i < nuped; i++)
      informative[i] = true;
    return;
  }
  for (i = 0; i < nuped; i++) {
    informative[i] = false;
    FORLIM1 = LINK->endped[i];
    for (j = LINK->startped[i]; j <= FORLIM1; j++) {
      if (person[j]->foff != NULL) {
	nchild = 0;
	child = person[j]->foff;
	do {
	  nchild++;
	  if (person[j]->male)
	    child = child->nextpa;
	  else
	    child = child->nextma;
	} while (child != NULL);
	count = 0;
	if (nchild > 1 || (nchild == 1 && person[j]->pa != NULL)) {
	  for (k = 0; k < mlocus; k++) {
	    WITH = person[j]->phen[k];
	    WITH1 = thislocus[k];
	    if (WITH1->which != binary_)
	      count++;
	    else {
	      if (WITH->phenf == 0)
		count++;
	      else {
		l = 0;
		FORLIM3 = WITH1->nallele;
                if (binformat == WITH1->format) {
		  for (m = 1; m <= FORLIM3; m++) {
		    if ((unsigned long)m < 32 &&
			((1L << m) & WITH->phenf) != 0)
		      l++;
		  }
                }
                else
                  if((WITH->allele1 > 0) && (WITH->allele2 > 0) &&
                     (WITH->allele1 != WITH->allele2))
                    l = 2;
		if (l > 1)
		  count++;
	      }
	    }
	  }
	}
	if (count > 1)
	  informative[i] = true;
      }
    }
  }
}

Local Void getind(id, LINK)
long *id;
struct LOC_readped *LINK;
{
  fscanf(ipedfile, "%ld", id);
  if (*id == 0)
    return;
  *id += LINK->sequence;
  if (*id > maxind)
    inputerror(13L, *id, *id, LINK->LINK);
  if (person[*id] == NULL) {
    person[*id] = (thisperson *)Malloc(sizeof(thisperson));
    if (person[*id] == NULL)
      malloc_err("person");
  }
}  /*getind*/

Local Void multimarriage(p, LINK)
thisperson **p;
struct LOC_readped *LINK;
{
  thisperson *q, *child, *WITH;

  if ((*p)->foff == NULL) {
    (*p)->multi = false;
    return;
  }
  WITH = *p;
  if (WITH->male)
    q = WITH->foff->ma;
  else
    q = WITH->foff->pa;
  child = WITH->foff;
  (*p)->multi = false;
  do {
    if (WITH->male) {
      WITH->multi = (q == child->ma);
      child = child->nextpa;
    } else {
      WITH->multi = (q == child->pa);
      child = child->nextma;
    }
  } while (!(child == NULL || WITH->multi));
}  /*multimarriage*/


Void readped(LINK)
struct LOC_inputdata *LINK;
{
  struct LOC_readped V;
  /*profield reads column 9 of pedin.dat which indicates proband or
    loopbreaker status*/
  long i, newid, sex_, profield, newped, nuperson, thisone, thisped;
  thisperson *holdloop;
  thisperson *WITH;
  int specialproband; /*seen a proband that is a loopbreaker in this pedigree*/
#if LOOP_BREAKERS
  int copyIndex;
#endif

  V.LINK = LINK;
  for (i = 0; i <= maxind; i++)
    person[i] = NULL;
  V.sequence = 0;
  nuperson = 0;
  nuped = 1;
  specialproband = 0;
#if LOOP_BREAKERS
  for (i = 0; i < maxloop; i++) 
    for(copyIndex = 0; copyIndex < maxloop; copyIndex++)
      looppers[nuped - 1][i][copyIndex] = NULL; /*G OK*/
#else
  for (i = 0; i < maxloop; i++) {
    looppers[nuped - 1][i][0] = NULL;  /*G OK*/
    looppers[nuped - 1][i][1] = NULL;  /*G OK*/
  }
#endif
  proband[nuped - 1] = NULL;
  fscanf(ipedfile, "%ld", &newped);
  thisped = newped;
  V.startped[0] = 1;
  while (!P_eof(ipedfile)) { 
    nuperson++;
    getind(&thisone, &V);
    if (proband[nuped - 1] == NULL)
      proband[nuped - 1] = person[thisone];
    WITH = person[thisone];
    WITH->ped = thisped;
    WITH->id = thisone;
    getind(&newid, &V);
    WITH->pa = person[newid];
    getind(&newid, &V);
    WITH->ma = person[newid];
    getind(&newid, &V);
    WITH->foff = person[newid];
    getind(&newid, &V);
    WITH->nextpa = person[newid];
    getind(&newid, &V);
    WITH->nextma = person[newid];
    fscanf(ipedfile, "%ld", &sex_);
    if (sex_ != 1 && sex_ != 2)
      inputerror(11L, WITH->id, sex_, LINK);
    if (sex_ == 1)
      WITH->male = true;
    else
      WITH->male = false;
    WITH->unknown = false;
    WITH->inloop = 0;
    WITH->loopdepend = false; /*A. A. Schaffer added this line and next*/
    WITH->loopneeded = false;
    fscanf(ipedfile, "%ld", &profield);
    /*Diagnostic added by A. A. Schaffer 7/25/94*/
    if ((profield - 1) > maxloop)
      if (specialproband != nuped)
        specialproband = nuped;
      else {   /*second violation in this pedigree */
	fprintf(stderr, "\nYour pedigree has more loops than allowed by the constant maxloop");
	fprintf(stderr, "\nYou must increase the constant maxloop defined in commondefs.h and recompile");
	fprintf(stderr, "\nYou are encouraged to read the loops.ps document distributed with FASTLINK");
	fprintf(stderr, "\nThe program will exit politely to allow you to correct the problem\n");
	exit(EXIT_FAILURE);
    }
    if (profield == 1)
      proband[nuped - 1] = person[thisone];
#if LOOP_BREAKERS
    else if (profield > (maxloop +1) && proband == NULL) { 
      profield -= maxloop;
      proband[nuped - 1] = person[thisone];
    }
#endif
    else if (profield > 1 && profield - 1 <= maxloop) {
#if LOOP_BREAKERS
      looppers[nuped - 1][profield - 2][numCopies[nuped - 1][profield -2]] =
					person[thisone];
      numCopies[nuped - 1][profield - 2]++;
#else
      if (looppers[nuped - 1][profield - 2][1] == NULL)
	looppers[nuped - 1][profield - 2][1] = person[thisone]; /*G OK*/
      else
	looppers[nuped - 1][profield - 2][0] = person[thisone]; /*G OK*/
#endif
    }
    getphenotype(&person[thisone], &V);
    fscanf(ipedfile, "%*[^\n]");
    getc(ipedfile);
    if (!P_eof(ipedfile))
      fscanf(ipedfile, "%ld", &newped);
    if (newped == thisped)
      continue;
    V.sequence += nuperson;
    V.endped[nuped - 1] = V.sequence;
    nuperson = 0;
    nuped++;
    if (nuped > maxped)
      inputerror(12L, newped, nuped, LINK);
    V.startped[nuped - 1] = V.sequence + 1;
#if LOOP_BREAKERS
  for (i = 0; i < maxloop; i++) 
    for(copyIndex = 0; copyIndex < maxloop; copyIndex++)
      looppers[nuped - 1][i][copyIndex] = NULL; /*G OK*/
#else
    for (i = 0; i < maxloop; i++) {
      looppers[nuped - 1][i][0] = NULL; /*G OK*/
      looppers[nuped - 1][i][1] = NULL; /*G OK*/
    }
#endif
    proband[nuped - 1] = NULL;
    thisped = newped;
  }
  totperson = V.sequence + nuperson;
  V.endped[nuped - 1] = totperson;
  for (newped = 0; newped < nuped; newped++) {
    for (i = 0; i < maxloop; i++) {
      if (looppers[newped][i][0] == NULL) /*G OK*/
	looppers[newped][i][1] = NULL;    /*G OK*/
      else {
#if LOOP_BREAKERS
        for(copyIndex = 0; copyIndex < numCopies[newped][i]; copyIndex++) {
          looppers[newped][i][copyIndex]->inloop = i + 1; /*G OK*/
        }
        holdloop = NULL;
        for(copyIndex = 1; copyIndex < numCopies[newped][i]; copyIndex++) {
          if(looppers[newped][i][copyIndex]->pa != NULL) { /*G OK*/
            holdloop = looppers[newped][i][0];   /*G OK*/
            break;
	  }
	}
	if (holdloop && holdloop->pa == NULL) {
	  looppers[newped][i][0] = looppers[newped][i][copyIndex]; /*G OK*/
	  looppers[newped][i][copyIndex] = holdloop; /*G OK*/
	}
#else
	looppers[newped][i][0]->inloop = i + 1; /*G OK*/
	looppers[newped][i][1]->inloop = i + 1; /*G OK*/
	if (looppers[newped][i][0]->pa == NULL && /*G OK*/
	    looppers[newped][i][1]->pa != NULL) { /*G OK*/
	  holdloop = looppers[newped][i][0];
	  looppers[newped][i][0] = looppers[newped][i][1]; /*G OK*/
	  looppers[newped][i][1] = holdloop; /*G OK*/
	}
#endif
      }
    }
  }
  for (thisone = 1; thisone <= totperson; thisone++)
    multimarriage(&person[thisone], &V);
  getinformative(&V);
}  /*readped*/

/* Local variables for getlocus: */
struct LOC_getlocus {
  struct LOC_readloci *LINK;
  long system;
} ;

Local Void getbin(locus, system, LINK)
locusvalues **locus;
long *system;
struct LOC_getlocus *LINK;
{
  long i, j, k;
  locusvalues *WITH;
  long FORLIM, FORLIM1;

  fscanf(datafile, "%ld%*[^\n]", &LINK->LINK->LINK->nfactor[*system - 1]);
  getc(datafile);
  if ((LINK->LINK->LINK->nfactor[*system - 1] > maxfact) ||
      (LINK->LINK->LINK->nfactor[*system - 1] > (sizeof(long) * 8 - 1)))
    inputerror(8L, *system, LINK->LINK->LINK->nfactor[*system - 1],
	       LINK->LINK->LINK);
  if (LINK->LINK->LINK->nfactor[*system - 1] <= 0)
    inputerror(9L, *system, LINK->LINK->LINK->nfactor[*system - 1],
	       LINK->LINK->LINK);
  WITH = *locus;
  FORLIM = WITH->nallele;
  for (i = 0; i < FORLIM; i++)
    WITH->UU.allele[i] = 0;
  FORLIM = WITH->nallele;
  for (i = 0; i < FORLIM; i++) {
    FORLIM1 = LINK->LINK->LINK->nfactor[*system - 1];
    for (j = 1; j <= FORLIM1; j++) {
      fscanf(datafile, "%ld", &k);
      if (k == 1)
	WITH->UU.allele[i] = ((long)WITH->UU.allele[i]) | (1 << ((long)j));
    }
  }
  fscanf(datafile, "%*[^\n]");
  getc(datafile);
}

Local Void getnumber(locus, system, LINK)
locusvalues **locus;
long *system;
struct LOC_getlocus *LINK;
{
  long i;
  locusvalues *WITH;
  long FORLIM;

  WITH = *locus;
  FORLIM = WITH->nallele;
  for (i = 1; i <= FORLIM; i++)
    WITH->UU.allele[i - 1] = 1;
}


Local Void getpen(locus, LINK)
locusvalues **locus;
struct LOC_getlocus *LINK;
{
  long i, j, k, l;
  locusvalues *WITH;
  long FORLIM, FORLIM1, FORLIM2;

  WITH = *locus;
  fscanf(datafile, "%ld%*[^\n]", &(*locus)->UU.U0.nclass);
  getc(datafile);
  if (WITH->UU.U0.nclass > maxliab)
    inputerror(28L, LINK->system, WITH->UU.U0.nclass, LINK->LINK->LINK);
  if (WITH->UU.U0.nclass <= 0)
    inputerror(29L, LINK->system, WITH->UU.U0.nclass, LINK->LINK->LINK);
  FORLIM = WITH->UU.U0.nclass;
  for (l = 0; l < FORLIM; l++) {
    FORLIM1 = WITH->nallele;
    for (i = 1; i <= FORLIM1; i++) {
      FORLIM2 = WITH->nallele;
      for (j = i - 1; j < FORLIM2; j++) {
	fscanf(datafile, "%lf", &(*locus)->UU.U0.pen[i][j][2][l]);
	if ((unsigned)WITH->UU.U0.pen[i][j][2][l] > 1.0)
	  inputerror(30L, LINK->system, LINK->system, LINK->LINK->LINK);
	WITH->UU.U0.pen[i][j][1][l] = 1 - WITH->UU.U0.pen[i][j][2][l];
	WITH->UU.U0.pen[i][j][0][l] = 1.0;
	for (k = 0; k <= 2; k++)
	  WITH->UU.U0.pen[j + 1][i - 1][k][l] = WITH->UU.U0.pen[i][j][k][l];
      }
    }
    fscanf(datafile, "%*[^\n]");
    getc(datafile);
    FORLIM1 = WITH->nallele;
    for (i = 0; i < FORLIM1; i++)
      WITH->UU.U0.pen[0][i][0][l] = 1.0;
    if (sexlink) {
      FORLIM1 = WITH->nallele;
      for (i = 0; i < FORLIM1; i++) {
	fscanf(datafile, "%lf", &(*locus)->UU.U0.pen[0][i][2][l]);
	if ((unsigned)WITH->UU.U0.pen[0][i][2][l] > 1.0)
	  inputerror(30L, LINK->system, LINK->system, LINK->LINK->LINK);
      }
      FORLIM1 = WITH->nallele;
      for (i = 0; i < FORLIM1; i++)
	WITH->UU.U0.pen[0][i][1][l] = 1.0 - WITH->UU.U0.pen[0][i][2][l];
      fscanf(datafile, "%*[^\n]");
      getc(datafile);
    }
  }
}

Local Void getquan(locus, privelege, LINK)
locusvalues **locus;
boolInt privelege;
struct LOC_getlocus *LINK;
{
  /*Get information of a quantitative locus. Privelege says whether it is
    a priveleged locus or not*/
  long i, j, k;
  locusvalues *WITH;
  long FORLIM, FORLIM1, FORLIM2;

  WITH = *locus;
  fscanf(datafile, "%ld%*[^\n]", &(*locus)->UU.U1.ntrait);
  getc(datafile);
  if (WITH->UU.U1.ntrait > maxtrait)
    inputerror(31L, LINK->system, WITH->UU.U1.ntrait, LINK->LINK->LINK);
  if (WITH->UU.U1.ntrait <= 0)
    inputerror(32L, LINK->system, WITH->UU.U0.nclass, LINK->LINK->LINK);
  FORLIM = WITH->UU.U1.ntrait;
  for (i = 0; i < FORLIM; i++) {
    FORLIM1 = WITH->nallele;
    for (j = 1; j <= FORLIM1; j++) {
      FORLIM2 = WITH->nallele;
      for (k = j; k <= FORLIM2; k++) {
	fscanf(datafile, "%lf", &(*locus)->UU.U1.pm[j][k - 1][i]);
	WITH->UU.U1.pm[k][j - 1][i] = WITH->UU.U1.pm[j][k - 1][i];
      }
    }
  }
  fscanf(datafile, "%*[^\n]");
  getc(datafile);
  if (privelege && LINK->LINK->nupriv != lastpriv)
    return;
  FORLIM = WITH->UU.U1.ntrait;
  for (i = 0; i < FORLIM; i++) {
    FORLIM1 = WITH->UU.U1.ntrait;
    for (j = i; j < FORLIM1; j++) {
      fscanf(datafile, "%lf", &(*locus)->UU.U1.vmat[i][j]);
      if (i + 1 == j + 1 && WITH->UU.U1.vmat[i][j] <= 0.0)
	inputerror(33L, LINK->system, LINK->system, LINK->LINK->LINK);
      WITH->UU.U1.vmat[j][i] = WITH->UU.U1.vmat[i][j];
    }
  }
  fscanf(datafile, "%*[^\n]");
  getc(datafile);
  invert(WITH->UU.U1.vmat, WITH->UU.U1.ntrait, &WITH->UU.U1.det);
  WITH->UU.U1.det = 1 / sqrt(WITH->UU.U1.det);
  fscanf(datafile, "%lf%*[^\n]", &(*locus)->UU.U1.conmat);
  getc(datafile);
  if (WITH->UU.U1.conmat <= 0)
    inputerror(34L, LINK->system, LINK->system, LINK->LINK->LINK);
  WITH->UU.U1.conmat = 1 / WITH->UU.U1.conmat;
  WITH->UU.U1.contrait = 1.0;
  FORLIM = WITH->UU.U1.ntrait;
  for (i = 1; i <= FORLIM; i++)
    WITH->UU.U1.contrait *= WITH->UU.U1.conmat;
  WITH->UU.U1.contrait = sqrt(WITH->UU.U1.contrait);
}

Void getlocus(system_, LINK)
long system_;
struct LOC_readloci *LINK;
{
  struct LOC_getlocus V;
  long j;
  locusvalues *WITH, *WITH1;
  long FORLIM;
  int TEMP;

  V.LINK = LINK;
  V.system = system_;
  thislocus[V.system - 1] = (locusvalues *)Malloc(sizeof(locusvalues));
  if (thislocus[V.system - 1] == NULL)
    malloc_err("entry in thislocus");
  WITH = thislocus[V.system - 1];
  WITH->privlocus = NULL;
  fscanf(datafile, "%ld%d", &LINK->whichtype,
	 &thislocus[V.system - 1]->nallele);
#if LOOPSPEED /* Dylan added */
  /*compute number of genotypes at this locus*/
  thislocus[V.system - 1]->fgeno = thislocus[V.system - 1]->nallele *
    (thislocus[V.system - 1]->nallele + 1) / 2;
#endif
  if (LINK->whichtype < 0 && LINK->whichtype > 4)
    inputerror(5L, V.system, LINK->whichtype, LINK->LINK);
  if (WITH->nallele > maxall)
    inputerror(6L, V.system, WITH->nallele, LINK->LINK);
  if (WITH->nallele <= 0)
    inputerror(7L, V.system, WITH->nallele, LINK->LINK);
  switch (LINK->whichtype) {

  case 0:
    WITH->which = quantitative;
    break;

  case 1:
    WITH->which = affection;
    break;

  case 2:
    WITH->which = binary_;
    break;

  case 3:
    WITH->which = binary_;
    break;
  }
  WITH->format = LINK->whichtype;
  if (lastpriv == 0) {
    fscanf(datafile, "%*[^\n]");
    getc(datafile);
  } else {
    fscanf(datafile, "%ld%*[^\n]", &LINK->whichtype);
    getc(datafile);
    if (LINK->whichtype == 0 || LINK->whichtype == 1) {
      LINK->nupriv++;
      WITH->privlocus = (locusvalues *)Malloc(sizeof(locusvalues));
      if (WITH->privlocus == NULL)
        malloc_err("privlocus field");
      WITH->privlocus->nallele = WITH->nallele;
      WITH1 = WITH->privlocus;
      switch (LINK->whichtype) {

      case 0:
	WITH1->which = quantitative;
	break;

      case 1:
	WITH1->which = affection;
	break;
      }
    }
  }
  if (!disequi) {
    FORLIM = WITH->nallele;
    for (j = 0; j < FORLIM; j++) {
      fscanf(datafile, "%lf", &thislocus[V.system - 1]->freq[j]);
      if (thislocus[V.system - 1]->freq[j] <= minfreq) {
        fprintf(stderr, "\nWarning: one of your allele frequencies is dangerously low\n");
	if (j == 0) {
	  exit(EXIT_FAILURE);
	}
      }
    }
    fscanf(datafile, "%*[^\n]");
    getc(datafile);
    if (disfreqs && (binary_ == WITH->which)) {
      for (j = 0; j < FORLIM; j++) {
	fscanf(datafile, "%lf", &thislocus[V.system - 1]->freqcond[j]);
      }
      fscanf(datafile, "%*[^\n]");
      getc(datafile);
    }
  }
  switch (WITH->which) {

  case binary_:
    if (WITH->format == allformat)
      getnumber(&thislocus[V.system - 1], &V.system, &V);
    else
      getbin(&thislocus[V.system - 1], &V.system, &V);
    break;

  case affection:
    getpen(&thislocus[V.system - 1], &V);
    break;

  case quantitative:
    getquan(&thislocus[V.system - 1], false, &V);
    break;
  /* cgh -- added for gcc */
  case null_:
  default:
    break;
  }
  if (WITH->privlocus != NULL) {
    switch (WITH->privlocus->which) {

    case affection:
      getpen(&WITH->privlocus, &V);
      break;

    case quantitative:
      getquan(&WITH->privlocus, true, &V);
      break;
    /* cgh -- gcc */
    case binary_:
    case null_:
    default:
      break;
    }
  }
  if (LINK->nupriv == lastpriv && lastpriv != 0)
    lastpriv = V.system;
  if (!risk)
    return;
  if (V.system == risksys) {
    fscanf(datafile, "%d%*[^\n]", &TEMP);
    getc(datafile);
    riskall = TEMP;
  }
  if (riskall > thislocus[V.system - 1]->nallele)
    inputerror(35L, V.system, (int)riskall, LINK->LINK);
  if (riskall < 0)
    inputerror(36L, V.system, (int)riskall, LINK->LINK);
}


Void gettheta(sex_, LINK)
thetavalues **sex_;
struct LOC_readloci *LINK;
{
  thetarray oldtheta;
  long i;
  thetavalues *WITH;
  long FORLIM;

#if !PARALLEL
  *sex_ = (thetavalues *)Malloc(sizeof(thetavalues));
  for (i = 0; i < maxlocus; i++)
    (*sex_)->theta[i] = 0.0;
#endif
  nuneed = 7;
  for(i = 2; i < mlocus; i++)
    nuneed = 5 * nuneed - 3;
  if (*sex_ == NULL)
    malloc_err("item of type thetavalues");
   /*Next line added by A. A. Schaffer*/
  (*sex_)->segprob = (double*) malloc(nuneed * sizeof(double));
  if ((*sex_)->segprob == NULL)
    malloc_err("a segprob array in procedure gettheta");
  WITH = *sex_;
  if (*sex_ == maletheta || readfemale) {
    FORLIM = mlocus - 2;
    for (i = 0; i <= FORLIM; i++)
      fscanf(datafile, "%lf", &(*sex_)->theta[i]);
    if (interfer && !mapping)
      fscanf(datafile, "%lf", &(*sex_)->theta[mlocus - 1]);
    fscanf(datafile, "%*[^\n]");
    getc(datafile);
  } else {
    fscanf(datafile, "%lf%*[^\n]", &distratio);
    getc(datafile);
  }
  /*     FOR j:=1 TO maxneed DO segprob[j]:=0.0;*/
  if (!interfer || mapping)
    return;
  memcpy(oldtheta, WITH->theta, sizeof(thetarray));
  WITH->theta[0] = (oldtheta[0] + oldtheta[mlocus - 1] - oldtheta[mlocus - 2]) / 2.0;
  WITH->theta[mlocus - 2] =
    (oldtheta[mlocus - 2] + oldtheta[mlocus - 1] - oldtheta[0]) / 2.0;
  WITH->theta[mlocus - 1] =
    (oldtheta[0] + oldtheta[mlocus - 2] - oldtheta[mlocus - 1]) / 2.0;
  for (i = 0; i < mlocus; i++) {   /*=ln(1/0.0001-1.0)*/
    if (WITH->theta[i] > 0.0) {
      if (WITH->theta[i] < 0.999)
	WITH->theta[i] = log(1.0 / WITH->theta[i] - 1.0);
      else
	WITH->theta[i] = -6.9;   /*=ln(1/0.999-1.0)*/
    } else
      WITH->theta[i] = 9.21;
  }
}


/*
   Moved from ilinputcode.c, mlinputcode.c, and liinputcode.c by
   Dylan in late 1994.
*/
Void inputdata()
{
  struct LOC_inputdata V;

  readloci(&V);
#ifdef ILINK
  setiterate(&V);
#endif
  readped(&V);
#if LOOPSPEED
  readsignature();
#else
  readspeed(&V);
#endif /*LOOPSPEED*/
}



