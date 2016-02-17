/* This file contains a few input/output routines for use with the faster
   versions of LODSCORE, ILINK, LINKMAP, and MLINK described in:
   R. W. Cottingham, Jr., R. M. Idury, and A. A. Schaffer
   Faster Sequential Genetic Linkage Computations
   American Journal of Human Genetics, 53(1993), pp. 252--263
   and A. A. Schaffer, S. K. Gupta, K. Shriram, and R. W. Cottingham, Jr.,
   Avoiding Recomputation in Linkage Analysis
   Human Heredity 44(1994), pp. 225-237. */


/* Two routines taken from 
 * "p2c"  Copyright (C) 1989, 1990, 1991 Free Software Foundation.
 * By Dave Gillespie, daveg@csvax.cs.caltech.edu.  Version --VERSION--.
 * This file may be copied, modified, etc. in any way.  It is not restricted
 * by the licence agreement accompanying p2c itself.
 */


/* Check if at end of file, using Pascal "eof" semantics.  End-of-file for
   stdin is broken; remove the special case for it to be broken in a
   different way. */
#include <stdio.h>
#include "commondefs.h"

#if defined(LODSCORE)
#include "lodefs.h"
#endif
#if defined(ILINK)
#include "ildefs.h"
#endif
#if defined(LINKMAP)
#include "lidefs.h"
#endif
#if defined(MLINK)
#include "mldefs.h"
#endif


int P_eof(f) 
FILE *f;
{
    register int ch;

    if (feof(f))
	return 1;
    if (f == stdin)
	return 0;    /* not safe to look-ahead on the keyboard! */
    ch = getc(f);
    if (ch == EOF)
	return 1;
    ungetc(ch, f);
    return 0;
}


/* Check if at end of line (or end of entire file). */
int P_eoln(f)
FILE *f;
{
    register int ch;

    ch = getc(f);
    if (ch == EOF)
        return 1;
    ungetc(ch, f);
    return (ch == '\n');
}

/*This routine is a simple exit routine, when the sytem is out of memory*/
void malloc_err(message)
char * message;
{
  FILE *errorfile;
  time_t secondsNow;

  fprintf(stderr, "\nProblem with malloc, probably not enough space\n");
  fprintf(stderr, "Problem occurred when allocating %s\n", message);

  errorfile = fopen("FASTLINK.err","a");
  if (errorfile) {
    time (&secondsNow);
    fprintf(errorfile,"\n%s",ctime(&secondsNow));
    fprintf(errorfile, "\nProblem with malloc, probably not enough space\n");
    fprintf(errorfile, "Problem occurred when allocating %s\n", message);
    fclose(errorfile);
  }
#if PARALLEL
  Tmk_errexit("\nProblem with malloc of private memory\n");
#else
  exit(EXIT_FAILURE);
#endif /*PARALLEL*/
}



/*The routine check_consatnts checks to see if some constants
have been set appropriately and reports the maximum number of haplotypes*/
void check_constants()
{
  int prod, l;
#if ALLELE_SPEED
  int pedprod, worst_ped, pedidx;
#endif
  boolean shouldexit, firstprint;
  time_t secondsNow;
  FILE *errorfile = NULL;     /* cgh -- made NULL for gcc */  

  locusvalues **templocus;
#if defined(LODSCORE)
  int locus1, locus2;
#endif

   /* Added by Dylan in late 1994 */  /*Done in 2.3P*/
#if (!LOOPSPEED && !defined(LODSCORE))
  {
    FILE *f;

    f = fopen(LOOPFILE_NAME, "r");
    if ( f != NULL ) {
      fclose(f);
      printf("Warning:  LOOPSPEED is not set for this run, but\n");
      printf("          %s exists, which suggests that\n", LOOPFILE_NAME);
      printf("          LOOPSPEED is set in unknown.c.\n");
    }
  }
#endif

  shouldexit = false;
  firstprint = false;
  time(&secondsNow);

#if defined(LODSCORE)
#if ALLELE_SPEED
  prod = 0;
  locus1 = locuslist1[iplace - 1] - 1;
  locus2 = locuslist2[jplace - 1] - 1;
  for (pedidx = 0; pedidx < mymaxped; pedidx++) {
    pedprod = 1;
      if ((binary_ == holdlocus[locus1]->which) && 
          (allformat == holdlocus[locus1]->format))
        pedprod *= ped_new_allele_count[pedidx][locus1];
      else
        pedprod *= holdlocus[locus1]->nallele;

      if ((binary_ == holdlocus[locus2]->which) && 
          (allformat == holdlocus[locus2]->format))
        pedprod *= ped_new_allele_count[pedidx][locus2];
      else
        pedprod *= holdlocus[locus2]->nallele;
    if (pedprod > prod) {
      prod = pedprod;
      worst_ped = pedidx;
    }
  }
  printf("\nThe maximum number of haplotypes needed is %d, given by a pedigree in which:", prod);
    if ((binary_ == holdlocus[locus1]->which) &&
        (allformat == holdlocus[locus1]->format))
      printf("\n%d alleles are needed at locus %d",ped_new_allele_count[worst_ped][locus1], locus1 + 1);
    else
      printf("\n%d alleles are needed at locus %d",holdlocus[locus1]->nallele, locus1 + 1);
    if ((binary_ == holdlocus[locus2]->which) &&
        (allformat == holdlocus[locus2]->format))
      printf("\n%d alleles are needed at locus %d",ped_new_allele_count[worst_ped][locus2], locus2 + 1);
    else
      printf("\n%d alleles are needed at locus %d",holdlocus[locus2]->nallele, locus2 + 1);

#else
  prod = 1;
  for(l=0; l < mlocus; l++){
    printf("\nNumber of alleles at locus %d is %d",l+1,
	   holdlocus[l]->nallele);
    prod *= holdlocus[l]->nallele;
  }
  printf("\nThe number of haplotypes used is %d", prod);
#endif /*ALLELE_SPEED*/
#else /*defined(LODSCORE)*/
#if ALLELE_SPEED
  prod = 0;
  for (pedidx = 0; pedidx < mymaxped; pedidx++) {
    pedprod = 1;
    for(l = 0; l < mlocus; l++)
      if ((binary_ == thislocus[l]->which) && 
          (allformat == thislocus[l]->format))
        pedprod *= ped_new_allele_count[pedidx][l];
      else
        pedprod *= thislocus[l]->nallele;
    if (pedprod > prod) {
      prod = pedprod;
      worst_ped = pedidx;
    }
  }
  printf("\nThe maximum number of haplotypes needed is %d, given by a pedigree in which:", prod);
  for(l = 0; l < mlocus; l++) 
    if ((binary_ == thislocus[l]->which) &&
        (allformat == thislocus[l]->format))
      printf("\n%d alleles are needed at locus %d",ped_new_allele_count[worst_ped][l], l + 1);
    else
      printf("\n%d alleles are needed at locus %d",thislocus[l]->nallele, l + 1);

#else
  prod = 1;
  for(l=0; l < mlocus; l++){
    printf("\nNumber of alleles at locus %d is %d",l+1,
	   thislocus[l]->nallele);
    prod *= thislocus[l]->nallele;
  }
  printf("\nThe number of haplotypes used is %d", prod);
#endif /*ALLELE_SPEED*/
#endif /*defined(LODSCORE) */
  if (maxseg > BIGMAXSEG){
    printf("\nYour value of maxseg is dangerously high.");
    printf("\nIt is safest to keep maxseg below %d and maxlocus below %d",BIGMAXSEG,BIGMAXLOCUS);
    printf("\nRemember that maxlocus is only a limit on how many loci are in one analysis.");
    printf("\nYou can have many more loci in your pedin.dat file.");
    if (!firstprint)
      errorfile = fopen("FASTLINK.err", "a");
    if (errorfile) {
      fprintf(errorfile,"\n%s", ctime(&secondsNow));
      fprintf(errorfile,"\nYour value of maxseg is dangerously high.");
      fprintf(errorfile,"\nIt is safest to keep maxseg below %d and maxlocus below %d",BIGMAXSEG,BIGMAXLOCUS);
      firstprint = true;
    }
    shouldexit = true;
  }
  printf("\n");
  if (firstprint && errorfile)
     fclose(errorfile);
  if (shouldexit)
#if PARALLEL
    Tmk_errexit("\n Problem with constants");
#else
    exit(EXIT_FAILURE);
#endif /*PARALLEL*/
} 



/* dwix: begin */

/* this struct is a linked list of integers used to store the pedigree */
/* numbers, which will be saved with an association to an index */
struct int_list_
  {
  int ped_number;
  struct int_list_ *next;
  };
typedef struct int_list_ int_list;


/******************************************************************************
function: read_till_blank_line

this function will read from the file given and return the file read
through until a blank line was encountered.  This function returns the
file so the next read will read the line after the blank line, not the
blank line itself.
******************************************************************************/
static int read_till_blank_line(infile)
FILE* infile;  
{
int blank_line = false;                /* is this a blank line */
char intext[DEFAULT_STRING_LENGTH];    /* the line read from the file */
char *tstr;                            /* a temporary char * */

while (!blank_line)   /* while we have not reached a blank line */
  {
  fgets(intext,DEFAULT_STRING_LENGTH,infile);     /* get new line */
  tstr = intext;                                  /* scan whitespace */
  while ((tstr[0]==' ') || (tstr[0]=='\011')) tstr++;
  /* if after scanning through whitespace, we are at then end of tstr, */
  /* then we had a blank line */
  if (tstr[0]=='\0') blank_line=true;
  }
}



/******************************************************************************
function: foundped

this function will scan ipedfile.dat and read the pedigrees in that
file.  It uses a linked list to save all the different pedigree's it
finds.  As it goes through the pedigrees, it counts them, and returns
the count.  The linked list is converted into an array of integers.
This array is the global variable pedidx2num.  It also saves the
calculation in the static variable saved_ped, so there is no
recomputation if the function is called more than once.
******************************************************************************/
static int foundped()
{
static int saved_ped = 0;                /* saved result */
FILE *ipedfile = NULL;                   /* the ipedfile */
char inputline[DEFAULT_STRING_LENGTH];   /* line read from ipedfile */
int ped_count = 0;                       /* pedigrees counted */
int inped;                               /* pedigree read */
int prev_ped = 0;                        /* the last pedigree read */
int_list *ped_list = NULL;               /* linked list to save pedigrees */
int_list *tmp_ped_list = NULL;           /* temp pointer into ped_list */
int count;                               /* counter for loops */

/* if we ran before and saved the result - return it */
if ( saved_ped != 0 ) return(saved_ped);


/* open ipedfile.dat */
ipedfile = fopen("ipedfile.dat","r");
if (ipedfile == NULL) ipedfile = fopen("pedfile.dat","r");
if (ipedfile == NULL)  /* if the open failed - print error and exit */
  {
  /* dwix err msg */
  fprintf(stderr,"Error opening ipedfile.dat and pedfile.dat.\n");
  exit(EXIT_FAILURE);
  }


ped_list = NULL;
/* read all the pedigrees and process */
while (!feof(ipedfile))  /* go through all of file */
  {
  fgets(inputline,DEFAULT_STRING_LENGTH,ipedfile);  /* get a line */
  sscanf(inputline," %d",&inped);                   /* read it's ped */
  if (inped!=prev_ped)                     /* if a new ped number */
    {
    /* make new node in the list */
    tmp_ped_list = NULL;
    tmp_ped_list = (int_list *)malloc(sizeof(int_list));
    if (tmp_ped_list == NULL)  /* if error mallocing exit w/ error */
      malloc_err("tmp_ped_list");

    /* place data into new node */
    tmp_ped_list->ped_number = inped;  /* add ped number */
    tmp_ped_list->next = ped_list;     /* tack list onto end */
				       /* note: this makes the list */
				       /* backward */
    ped_list = tmp_ped_list;           /* make ped_list new leader */
    ped_count++;                       /* increase ped_count */
    prev_ped = inped;                  /* save the previos ped */
    }
  }
fclose(ipedfile);   /* close the file */

/* allocate space for the array to convert index to pedigree number */
pedidx2num = (int *)calloc(ped_count,sizeof(int));
if (pedidx2num==NULL)      /* if error exit w/ message */
  malloc_err("pedidx2num");

/* copy list data into array - go backwards */
tmp_ped_list = ped_list;
for (count=ped_count-1;count >= 0; count--)
  {
  pedidx2num[count] = tmp_ped_list->ped_number;
  tmp_ped_list = tmp_ped_list->next;
  }

/* delete allocated space for the linked list */
while (ped_list!=NULL)  /* while members to delete */
  {
  tmp_ped_list = ped_list;    /* save this member */
  ped_list = ped_list->next;  /* get next member */
  free(tmp_ped_list);         /* delete this member */
  }

saved_ped = ped_count;       /* save the result */
return(ped_count);           /* return the result */
}


/******************************************************************************
function: init_ped_loc_all

this function initializes the variable ped_loc_all after it has been
allocated.  This function initializes all entries in the 3d array to 0
******************************************************************************/
void init_ped_loc_all()
{
int a,b,c;

int pedidx;

locusvalues **templocus;

int loopbound;

#if defined(LODSCORE)
   templocus = holdlocus;
#else
   templocus = thislocus;
#endif


mymaxped = foundped();

if (mymaxped==0)
  {
  /* dwix err msg */
  fprintf(stderr,"foundped() found 0 pedigree's - wrong.\n");
  exit(EXIT_FAILURE);
  }
ped_loc_all = NULL;
ped_loc_all = (loc_all *)malloc(mymaxped*sizeof(loc_all));
if (ped_loc_all==NULL) 
  malloc_err("ped_loc_all");


#if ALLELE_SPEED
  ped_must_change_locations = NULL;
  ped_must_change_locations = (boolean *) malloc((mymaxped + 1)
						 * sizeof(boolean));
  ped_must_change_locations[0] = true;
  for(pedidx = 1; pedidx <= mymaxped; pedidx++)
    ped_must_change_locations[pedidx] = false;
  ped_new_allele_count = NULL;
  ped_new_allele_count = (new_allele_count *) 
     malloc(mymaxped * sizeof(new_allele_count));
#endif

#if defined(LODSCORE)
  loopbound = maxsys;
#else
  loopbound = maxlocus;
#endif

for (a=0;a<mymaxped;a++)
  {
  for (b=0;b<loopbound;b++)
    {
#if ALLELE_SPEED
   ped_new_allele_count[a][b] = 0;
#endif
    ped_loc_all[a][b][0].present = false;
    ped_loc_all[a][b][0].new_allele = 0;
    ped_loc_all[a][b][0].old_allele = 0;
    ped_loc_all[a][b][0].old_frequency = 0.0;
    ped_loc_all[a][b][0].new_frequency = 0.0;
    for (c=1;c<=maxall;c++)
      {
      ped_loc_all[a][b][c].present = false;
      ped_loc_all[a][b][c].new_allele = 0;
      ped_loc_all[a][b][c].old_allele = 0;
      ped_loc_all[a][b][c].old_frequency = 0.0;
      ped_loc_all[a][b][c].new_frequency = 0.0;
      }
    }
  }
return;
}

#if ALLELE_SPEED
/*adjust_alleles renumbers alleles on a pedigree by pedigree basis*/
int adjust_alleles()
{
  int pedidx, locidx, allidx, old_all, new_all;
  double total_freq_sum[maxsys], present_freq_sum;
  /*following declarations are for conditional allele frequencies, Morgan*/
  double cond_total_freq_sum[maxsys], cond_present_freq_sum; 
  int missed;
  
  locusvalues **templocus;

#if defined(LODSCORE)
   templocus = holdlocus;
#else
   templocus = thislocus;
#endif

  /*Compute frequency of alleles sum for each numbered allele locus */
  for(locidx = 0; locidx < mlocus; locidx++) {
    if ((binary_ == templocus[locidx]->which) &&
        (allformat == templocus[locidx]->format)) {
      total_freq_sum[locidx] = 0.0;
      if (disfreqs)
        cond_total_freq_sum[locidx] = 0.0;
      for(allidx = 0; allidx < templocus[locidx]->nallele; allidx++) {
        total_freq_sum[locidx] += templocus[locidx]->freq[allidx];
        if (disfreqs)
	  cond_total_freq_sum[locidx] += templocus[locidx]->freqcond[allidx];
      }
#if PARALLEL
   if (Tmk_proc_id == 0)
#endif
      if ((total_freq_sum[locidx] < (1.0 - epsilon)) ||
          (total_freq_sum[locidx] > (1.0 + epsilon))) {
	fprintf(stderr,"\nALLELE FREQUENCY ALERT: The frequencies at locus index");
        fprintf(stderr,"\n%d sum to %10.7lf rather than 1.0",locidx + 1, total_freq_sum[locidx]);
      }
    }
  }
  for(pedidx = 0; pedidx < nuped; pedidx++)
    for(locidx = 0; locidx < mlocus; locidx++) 
      if ((binary_ == templocus[locidx]->which) &&
          (allformat == templocus[locidx]->format)) {
        for(allidx =1; allidx <= templocus[locidx]->nallele; allidx++) {
          ped_loc_all[pedidx][locidx][allidx].old_frequency =
            templocus[locidx]->freq[allidx - 1];
          if (disfreqs)
	    ped_loc_all[pedidx][locidx][allidx].old_frequency_cond =
	      templocus[locidx]->freqcond[allidx - 1];
	}

        /*map old alleles to new alleles preserving present ones*/
        present_freq_sum = 0.0;
        if (disfreqs)
	  cond_present_freq_sum = 0.0;
        old_all = 1;
        new_all = 1;  
        missed = 0;
        while (old_all <= templocus[locidx]->nallele) {
          if (ped_loc_all[pedidx][locidx][old_all].present) {
             ped_loc_all[pedidx][locidx][old_all].new_allele = new_all;
             ped_loc_all[pedidx][locidx][new_all].old_allele = old_all;
             ped_loc_all[pedidx][locidx][new_all].new_frequency =
               ped_loc_all[pedidx][locidx][old_all].old_frequency;
             present_freq_sum +=
                ped_loc_all[pedidx][locidx][new_all].new_frequency;
             if (disfreqs) {
	       ped_loc_all[pedidx][locidx][new_all].new_frequency_cond =
		 ped_loc_all[pedidx][locidx][old_all].old_frequency_cond;
	       cond_present_freq_sum +=
		 ped_loc_all[pedidx][locidx][new_all].new_frequency_cond;
             }
             new_all++;
	  }
          else
            missed = old_all;
          old_all++;
        }
        
        /*create new allele if not all present as new allele new_all */
        if (new_all < old_all) {  /*missed must have a non-zero value*/
          ped_loc_all[pedidx][locidx][new_all].old_allele = missed;

            ped_loc_all[pedidx][locidx][new_all].new_frequency = 
              total_freq_sum[locidx]   - present_freq_sum;
            if (disfreqs)
	      ped_loc_all[pedidx][locidx][new_all].new_frequency_cond = 
		cond_total_freq_sum[locidx]   - cond_present_freq_sum;
          new_all++;
        }  
        ped_new_allele_count[pedidx][locidx] = new_all - 1;
     /*if this pedigree has different adjusted allele counts
       than the previous one, we must adjust the gene tables before
       working on this pedigree*/
        if ((pedidx > 0) &&
             (ped_new_allele_count[pedidx][locidx] !=
              ped_new_allele_count[pedidx -1][locidx]))
         ped_must_change_locations[pedidx] = true;
      }
}
#endif  /*ALLELE_SPEED*/

/******************************************************************************
function: allele_downcode_check

although named poorly - I can't think of anything else to name it.
This function will use the data in ped_loc_all to determine if the
data can be restructured to decrease the allele product, and thereby
increase the speed of the run.
******************************************************************************/
int allele_downcode_check()
{
int ped_check;       /* if there is a problem with a given pedigree */
int count, ped_count, loc_count, all_count;  /* counters for loops */
int one_missing;
locusvalues **templocus;
boolean firstError = true;
boolean alleleNotFound;

#if defined(LODSCORE)
templocus = holdlocus;
#else
templocus = thislocus;
#endif

for (ped_count = 0; ped_count < nuped; ped_count++) {
  ped_check = false;
  for (loc_count = 0; loc_count < mlocus; loc_count++) {
    /* scan loci */
    if ((binary_ == templocus[loc_count]->which) &&
        (allformat == templocus[loc_count]->format)) {
      /* for each pedigree-locus pair, we start with nothing wrong */
      one_missing = false;
      ped_check = false;
      for (all_count = 1;
	   all_count <= (templocus[loc_count]->nallele); all_count++) {
	/* scan alleles */
        if (!(ped_loc_all[ped_count][loc_count][all_count].present)) {
	  /* allele missing */
          if (one_missing) ped_check = true;
          else one_missing = true;
	}
      }
    
      if (ped_check) {
	/* dwix err msg */
	if (firstError) {
#if PARALLEL
          if (Tmk_proc_id == 0) {
#endif

	  fprintf(stderr,"\nSee README.allele\n");
          fprintf(stderr,"FASTLINK will renumber alleles for you automagically\n");
#if PARALLEL
	}
#endif
	  firstError = false;
	}
#if PARALLEL
          if (Tmk_proc_id == 0) {
#endif
	fprintf(stderr,"ALLELE DIAGNOSIS: pedigree %d does not use all alleles at locus index %d\n", pedidx2num[ped_count], loc_count + 1);
#if PARALLEL
      }
#endif
      }
    }
  }
}


free(pedidx2num);
return(0);  /* if reporting on peds - */
}

/* dwix: end */


#if ALLELE_SPEED
/*The following procedure written by A.A. Schaffer to adjust the
alleles of persons based on renumbering of alleles. Uses the fact that
old allele -> new allele mapping is montonic */

void allele_adjust_persons()
{
  int pedix, locidx, personidx;
  thisperson *next_to_mutate; /*next person to have alleles changed*/
  int prev_ped;
  locusvalues **templocus;

#if defined(LODSCORE)
  templocus = holdlocus;
#else
  templocus = thislocus;
#endif

  personidx = 1;
  pedix = -1;
  prev_ped = -1;
  while ((next_to_mutate = person[personidx]) != NULL) {
    if (prev_ped != (next_to_mutate->ped - 1)) {
      prev_ped = next_to_mutate->ped - 1;
      pedix++;
    }
    for(locidx = 0; locidx < mlocus; locidx++) {
      if ((binary_ == templocus[locidx]->which) &&
          (allformat == templocus[locidx]->format)) {
#if defined(LODSCORE)
        next_to_mutate->holdphen[locidx]->allele1 =
          ped_loc_all[pedix][locidx]
            [next_to_mutate->holdphen[locidx]->allele1].new_allele;
        next_to_mutate->holdphen[locidx]->allele2 =
          ped_loc_all[pedix][locidx]
            [next_to_mutate->holdphen[locidx]->allele2].new_allele;
#else
        next_to_mutate->phen[locidx]->allele1 =
          ped_loc_all[pedix][locidx]
            [next_to_mutate->phen[locidx]->allele1].new_allele;
        next_to_mutate->phen[locidx]->allele2 =
          ped_loc_all[pedix][locidx]
            [next_to_mutate->phen[locidx]->allele2].new_allele;
#endif
      }
    }
    personidx++;
  }
}
#endif

