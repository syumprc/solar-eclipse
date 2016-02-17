/* Using FCLOSE macro for linux compability...files set to null when closed */
/* 25 January 2000 SFBR */

/* Output from p2c, the Pascal-to-C translator */
/* From input file "unknown.p" */

/*Changes at Columbia U. indicated by "change"*/
/*10 July 1993*/

#define FCLOSE(file) (fclose (file), file=0)

/*
   Modified extensively by Dylan Cooper and Alejandro Schaffer
   to produce a file whose name
   is in the macro LOOPFILE_NAME) that is used to restrict the
   genotypes of unknowns in the presence of loops.  Modified in late
   1994 and summer of 1995.
 
   Please read README.loopfile for a description of the file created
   and the concept of loopbreaker vectors.
 
   The basic idea is (for each locus of each pedigree):
       -- find all the permutations (vectors) of single locus genotypes
          that the loop breakers can have
       -- traverse the pedigree from the proband
          for each loopbreaker vector and record the single locus
          genotypes each person may have if the loopbreakers have the
          single locus genotypes specified by the vector
       -- if there are no loops do a faster version of the old alg

  Modified by Alejandro Schaffer in late 1995 to renumber alleles when
  ALLELE_SPEED is set to 1. If 2 or more alleles are unused in a pedigree
  at a locus, all the unused alleles can be combined into 1.

  Allele renumbering affects the output of speedfile.dat and loopfile.dat,
  but not pedfile.dat, datafile.dat, or ipedfile.dat. The reasons for this
  choice are somewhat arbitrary. 


  Further changes were made in April,1996 by Tony Schurtz, in order to
  allow for backward compatibility with the old "speedfile.dat" format.
  Now the application always writes out speedfile.dat in the original
  file format (as though ALLELE_SPEED were zero).  If ALLELE_SPEED is set
  to one, a second speedfile is written, in the new format, to a file called
  "newspeedfile.dat"

  The changes to the source are quite systematic.  Each function which
  depended upon the value of ALLELE_SPEED was simply "cloned", one version
  retaining the ALLELE_SPEED = 0 code, and the other the ALELE_SPEED = 1
  code.  These were named x_old, and x_new, respectively, x being the old
  function name..  main() was cloned in the same way in order to ensure
  that the entire state of the application is reinitialized before
  newspeedfile.dat is generated.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>  /* added by Dylan for use in malloc_err() */

#if defined(vms)
#define unlink delete
#endif

/*
   If LOOPSPEED is defined, the file whose name is defined in the
   macro LOOPFILE_NAME is created.

   If defined here, define in commondefs.h as well.
*/
#if !defined(LOOPSPEED)
#define LOOPSPEED 1
#endif


#define aVersion         "5.20"   /*PRESENT VERSION OF LINKAGE*/
#define fVersion         "3.0P"   /*PRESENT VERSION OF FASTLINK*/

#if LOOPSPEED
/*
   max number of loop vectors that we will consider at any locus
   this is used to bound the number of loops we consider in unknown_poss

   we will always consider at least one loop if there are any
*/
#define max_vectors_considered 20000

int fewer_vects_size;  /* used for diagnostics when above macro is too big */
#endif

/*ALLELE_SPEED is used to set up renumbering of alleles when not all
  alleles are present in a pedigree*/
#if !defined(ALLELE_SPEED)
#define ALLELE_SPEED 1
#endif

#define ALLELE_SPEED_CONSTANT 61192012

#define maxlocus        2    /*MAXIMUM NUMBER OF LOCI*/
#define maxall          50   /*MAX NUMBER OF ALLELES AT A SINGLE LOCUS*/

#define maxgeno         (maxall * (maxall + 1) / 2)
    /*MAX NUMBER OF SINGLE LOCUS GENOTYPES*/
/*Different definition than in analysis programs!*/

#define maxind          10000 /*MAXIMUM NUMBER OF INDIVIDUALS IN A PEDIGREE*/
#define maxmarriage     100   /*MAXIMUM NUMBER OF MARRIAGES FOR ONE MALE*/

#ifndef maxloop
#define maxloop         8
#endif

#define maxfact         maxall
    /*MAXIMUM NUMBER OF BINARY CODES AT A SINGLE LOCUS*/

#define maxtrait        3
    /*MAXIMUM NUMBER OF QUANTITATIVE FACTORS AT A SINGLE LOCUS*/

#define missval         0.0   /*MISSING VALUES FOR QUANTITATIVE TRAITS*/

#define affall          2
    /*DISEASE ALLELE FOR QUANT. TRAITS OR AFFECTION STATUS*/
#define missaff         0   /*MISSING VALUE FOR AFFECTION STATUS*/
#define affval          2   /*CODE FOR AFFECTED INDIVIDUAL*/
#define maxliab         60   /*MAXIMUM NUMBER OF LIABILITY CLASSES*/
#define binformat       2
#define allformat       3

#define FileNotFound   10
#define true            1
#define false           0

#define Malloc          malloc
#define Free            free

#if LOOPSPEED
#define LOOPFILE_NAME "loopfile.dat"
#endif

#if !defined(ONE_ERROR_ONLY)
#define ONE_ERROR_ONLY   0   /*print one sure error, or all possible errors*/
#endif

#if !defined(DEPTH_MULTIPLE)
#define DEPTH_MULTIPLE   3
#endif

#if !defined(DOWN_CHECK)
#define DOWN_CHECK true
#endif

/*from p2c*/

typedef char boolean;

typedef boolean genotype[maxgeno];
typedef enum {
  peelup, peeldown
} direction;
typedef long binset;

typedef binset phenarray[maxall];
typedef boolean possvect[maxall][maxall];
typedef possvect possarray[maxlocus];
typedef enum {
  affection, quantitative, binary_
} locustype;

typedef struct locusvalues {
  long nallele;
  long maxallele;
  locustype which;
  double onefreq;
#if LOOPSPEED
  long fgeno; /* number of genotypes at this locus */
  long mgeno;
#endif
  union {
    long ntrait;
    struct {
      double pen[maxall + 1][maxall][3][maxliab];
      long nclass;
    } U0;
    struct {
      phenarray allele;
      long nfactor, format;
    } U2;
  } UU;
} locusvalues;

typedef struct phenotype {
  locustype which;
  double x[maxtrait];
  boolean missing;
  long aff;
  long  liability;
  binset phenf; /*encodes phenotype/genotype in pedin.dat as a
                    bit vector */
  int allele1;
  int allele2;
} phenotype;

typedef phenotype *indphen[maxlocus];

typedef struct thisarray {
  genotype genarray;
} thisarray;

typedef struct information {
  possarray possible;
} information;

typedef struct thisperson {
  long id, paid, maid, offid, npaid, nmaid, sex, profield, oldped, nseq;
  struct thisperson *pa, *ma, *foff, *nextpa, *nextma;
  thisarray *gen;
  indphen phen;
#if !LOOPSPEED
  information *store;
#endif
  boolean thisunknown[maxlocus];
  boolean unknown, multi, done, up, male;
  boolean downvisit;  /*added by A. A. Schaffer for diagnostic*/
} thisperson;

#if LOOPSPEED
/*
   These typedefs define the table to hold the valid loopbreaker vectors
   for each locus.
*/
typedef int *vector_for_loop;
typedef vector_for_loop *vector_for_locus;
typedef vector_for_locus *loop_vector_array;
#endif

#if LOOPSPEED
/*
   These typedefs define the table to hold the possible genotypes
   for each unknown person at each locus for each valid loopbreaker
   vector.
*/
typedef boolean *geno_for_loop_vector;
typedef geno_for_loop_vector *geno_for_locus;
typedef geno_for_locus *geno_for_unknown;
#endif

#if LOOPSPEED
/* basic linked list stuff */
typedef struct List_Elt {
  struct List_Elt *next;
  int value;
} list_elt;

typedef struct Simple_List {
  list_elt *head;
  list_elt *curr;
} simple_list;
#endif

/* Dylan
         -- changed haplotype from enumerated type of {a,b}
         -- done in mid 1994
*/
typedef unsigned char haplotype;

typedef long subhap[2];

typedef struct new_locus_info {
  int present;  /*is this allele present?*/
  int new_allele; /*what is this allele mapped to?*/
  int old_allele; /*inverse of mapping*/
} new_locus_info;

typedef int new_allele_count[maxlocus];

/*Adjusted allele counts for each pedigree and locus*/
new_allele_count *ped_new_allele_count;

/* Does this pedigree have different adjusted allele counts from
   the previous pedigree? */

boolean *ped_must_change_locations;

/* dwix: begin */

#define DEFAULT_STRING_LENGTH  200

/* type definitions for *allele_downcode_check* */
typedef new_locus_info  loc_all[maxlocus][maxall+1];
loc_all *ped_loc_all;
int currentped, currentlocus;
int *pedidx2num;
int mymaxped;
/* dwix: end */


static subhap seghap[maxgeno];
static locusvalues *thislocus[maxlocus];
static thisperson *person[maxind + 1];
static thisperson *proband, *loop1, *loop2;
static long risksys, mutsys, nsystem;
static long fgeno, mgeno;
static long nsequence, newped, whichsys, totperson;
static boolean sexlink, risk, disequi;
static FILE *speedfile, *datafile, *pedfile, *ipedfile;
static genotype gene;
static double one;   /*changed*/
static boolean makehomozygous;   /* Change - Added 7/8/93 */
static numind;  /*number of individuals in pedigree*/
static depth;  /*depth of recursion*/

/*The following three variables were added by A. A. Schaffer to
 help pinpoint Mendelian incompatibilities*/
static boolean incompat_traversal; /*found incompatibility on this traversal*/
static boolean first_incompat; /*Is this the first incomptibility */
static boolean detected_incompat[maxind+1]; /* array to record whether
                                               an incompatibility with this
                                               person has already been
                                               reported*/


#if LOOPSPEED
static int num_loops_considered; /* Added by Dylan late 1994 */
static int num_loops;
#endif

#if LOOPSPEED
static long genenumber[maxlocus][maxgeno][maxgeno];
#else
static long genenumber[maxgeno][maxgeno];
#endif

#if LOOPSPEED
/* file containing info when pedigree has loops */
static FILE *loopfile;

/* aliases the two loopbreakers for each loop */
static thisperson *looppers[maxloop][2];

/* table holding possible genotypes for unknowns.  Indexed by:
   - person[i]->nseq
   - locus
   - loopbreaker vector #
   - single locus genotype #
*/
static geno_for_unknown *unknown_poss = NULL;

/* array of valid loopbreakers */
static loop_vector_array loop_vectors = NULL;

/* holds number of loopbreaker vectors for each locus */
static int num_loop_vectors[maxlocus];
#endif

/* ----------------------- Procedures ------------------------------- */

/* Dylan -- Copied from iostuff.c */
/*
   This routine is a simple exit routine, when the sytem is out of memory
*/
void malloc_err(message)
char * message;
{
  FILE *errorfile;
  time_t secondsNow;

  fprintf(stderr, "\nProblem with malloc, probably not enough space\n");
  fprintf(stderr, "Problem occurred when allocating %s\n", message);
#if LOOPSPEED
  fprintf(stderr, "Reduce max_vectors_considered to %d.\n", fewer_vects_size);
#endif

  errorfile = fopen("FASTLINK.err","a");
  if (errorfile) {
    time (&secondsNow);
    fprintf(errorfile,"\n%s",ctime(&secondsNow));
    fprintf(errorfile, "\nProblem with malloc, probably not enough space\n");
    fprintf(errorfile, "Problem occurred when allocating %s\n", message);
    FCLOSE(errorfile);
	    }
  exit(EXIT_FAILURE);
  }

/* Two routines taken from */
/* "p2c"  Copyright (C) 1989, 1990, 1991 Free Software Foundation.
 * By Dave Gillespie, daveg@csvax.cs.caltech.edu.  Version --VERSION--.
 * This file may be copied, modified, etc. in any way.  It is not restricted
 * by the licence agreement accompanying p2c itself.
 */


/* Check if at end of file, using Pascal "eof" semantics.  End-of-file for
   stdin is broken; remove the special case for it to be broken in a
   different way. */

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

/* dwix: begin */

/* this struct is a linked list of integers used to store the pedigree */
/* numbers, which will be saved with an association to an index */
struct int_list_
  {
  int ped_number;
  struct int_list_ *next;
  };
typedef struct int_list_ int_list;


/***************************************************************************
function: read_till_blank_line

this function will read from the file given and return the file read
through until a blank line was encountered.  This function returns the
file so the next read will read the line after the blank line, not the
blank line itself.
***************************************************************************/
int read_till_blank_line(infile)
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



/**************************************************************************
function: foundped

this function will scan pedfile.dat and read the pedigrees in that
file.  It uses a linked list to save all the different pedigree's it
finds.  As it goes through the pedigrees, it counts them, and returns
the count.  The linked list is converted into an array of integers.
This array is the global variable pedidx2num.  It also saves the
calculation in the static variable saved_ped, so there is no
recomputation if the function is called more than once.
****************************************************************************/
int foundped()
{
static int saved_ped = 0;                /* saved result */
FILE *pdefile = NULL;                   /* the ipedfile */
char inputline[DEFAULT_STRING_LENGTH];   /* line read from ipedfile */
int ped_count = 0;                       /* pedigrees counted */
int inped;                               /* pedigree read */
int prev_ped = 0;                        /* the last pedigree read */
int_list *ped_list = NULL;               /* linked list to save pedigrees */
int_list *tmp_ped_list = NULL;           /* temp pointer into ped_list */
int count;                               /* counter for loops */

/* if we ran before and saved the result - return it */
if ( saved_ped != 0 ) return(saved_ped);


/* open pedfile.dat */
    pedfile = fopen("pedfile.dat", "r");
if (pedfile == NULL)  /* if the open failed - print error and exit */
  {
  /* dwix err msg */
  fprintf(stderr,"Error opening pedfile.dat in UNKNOWN.\n");
  exit(EXIT_FAILURE);
  }
ped_list = NULL;
/* read all the pedigrees and process */
while (!feof(pedfile))  /* go through all of file */
  {
  fgets(inputline,DEFAULT_STRING_LENGTH,pedfile);  /* get a line */
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
FCLOSE(pedfile);   /* close the file */
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

/***************************************************************************
function: init_ped_loc_all_old

this function initializes the variable ped_loc_all after it has been
allocated.  This function initializes all entries in the 3d array to 0
****************************************************************************/
void init_ped_loc_all_old()
{
int a,b,c;

int pedidx;

locusvalues **templocus;

   templocus = thislocus;

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

for (a=0;a<mymaxped;a++)
  {
  for (b=0;b<maxlocus;b++)
    {
    ped_loc_all[a][b][0].present = false;
    ped_loc_all[a][b][0].new_allele = 0;
    ped_loc_all[a][b][0].old_allele = 0;
    for (c=1;c<=maxall;c++)
      {
      ped_loc_all[a][b][c].present = false;
      ped_loc_all[a][b][c].new_allele = 0;
      ped_loc_all[a][b][c].old_allele = 0;
      }
    }
  }
return;
}

/***************************************************************************
function: init_ped_loc_all

this function initializes the variable ped_loc_all after it has been
allocated.  This function initializes all entries in the 3d array to 0
****************************************************************************/
void init_ped_loc_all_new()
{
int a,b,c;

int pedidx;

locusvalues **templocus;

   templocus = thislocus;

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

  ped_new_allele_count = NULL;
  ped_new_allele_count = (new_allele_count *)
     malloc(mymaxped * sizeof(new_allele_count));

for (a=0;a<mymaxped;a++)
  {
  for (b=0;b<maxlocus;b++)
    {
   ped_new_allele_count[a][b] = 0;
    ped_loc_all[a][b][0].present = false;
    ped_loc_all[a][b][0].new_allele = 0;
    ped_loc_all[a][b][0].old_allele = 0;
    for (c=1;c<=maxall;c++)
      {
      ped_loc_all[a][b][c].present = false;
      ped_loc_all[a][b][c].new_allele = 0;
      ped_loc_all[a][b][c].old_allele = 0;
      }
    }
  }
return;
}

/*adjust_alleles renumbers alleles on a pedigree by pedigree basis*/
int adjust_alleles()
{
  int pedidx, locidx, allidx, old_all, new_all;
  int missed;

  locusvalues **templocus;

   templocus = thislocus;

    pedidx = currentped;
    for(locidx = 0; locidx < nsystem; locidx++)
      if ((binary_ == templocus[locidx]->which) &&
           (allformat == templocus[locidx]->UU.U2.format)) {
        /*map old alleles to new alleles preserving present ones*/
        old_all = 1;
        new_all = 1;
        missed = 0;
        while (old_all <= templocus[locidx]->maxallele) {
          if (ped_loc_all[pedidx][locidx][old_all].present) {
             ped_loc_all[pedidx][locidx][old_all].new_allele = new_all;
             ped_loc_all[pedidx][locidx][new_all].old_allele = old_all;
             new_all++;
	  }
          else
            missed = old_all;
          old_all++;
	}

        /*create new allele if not all present as new allele new_all */
        if (new_all < old_all) {  /*missed must have a non-zero value*/
          ped_loc_all[pedidx][locidx][new_all].old_allele = missed;
          new_all++;
	}
        ped_new_allele_count[pedidx][locidx] = new_all - 1;
      }
}

/*The following procedure written by A.A. Schaffer to adjust the
alleles of persons based on renumbering of alleles. Uses the fact that
old allele -> new allele mapping is monotonic */

void allele_adjust_persons()
{
  int pedix, locidx, personidx;
  thisperson *next_to_mutate; /*next person to have alleles changed*/
  locusvalues **templocus;

  templocus = thislocus;

  personidx = 1;
  pedix = currentped;
  while ((next_to_mutate = person[personidx]) != NULL) {
    for(locidx = 0; locidx < nsystem; locidx++) {
      if ((binary_ == templocus[locidx]->which) &&
          (allformat == templocus[locidx]->UU.U2.format)) {
        next_to_mutate->phen[locidx]->allele1 =
          ped_loc_all[pedix][locidx]
            [next_to_mutate->phen[locidx]->allele1].new_allele;
        next_to_mutate->phen[locidx]->allele2 =
          ped_loc_all[pedix][locidx]
            [next_to_mutate->phen[locidx]->allele2].new_allele;
      }
    }
    personidx++;
  }
}


static void respond()
{
  /*Change - new
  printf("*** Press <Enter> to continue\n");
  scanf("%*[^\n]");
  getchar();*/
  exit(EXIT_FAILURE);
}



static void inputerror(nerror, par1, par2)
long nerror, par1, par2;
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
    printf("Number of iterated parameters greater then maxn\n");
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
      "Error detected reading datafile. Linkage disequilibrium is not allowed with this program\n");
    break;

  case 42:
    printf("Locus %5ld in lod score list exceeds nlocus %5ld\n", par1, par2);
    break;

  case 43:
    printf("Illegal locus number %5ld in lod score list\n", par1);
    break;

  case 44:
    printf("Error detected reading pedigree record %2ld. One 0 allele\n",
	   par1);
    break;
  }
  respond();   /*changed*/
}


static void inputwarning(nwarning, par1, par2)
long nwarning, par1, par2;
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
  respond();   /*changed*/
}

#if LOOPSPEED

/* writes the old speedfile.dat for compatibility.   Written by Dylan. */
static void writespeed_old()
{
  int i, j, a, b;
  thisperson *WITH;
  int allele_count;

  /* for compatibility, when loops exist, speedfile.dat does not eliminate
     anything */
  if (loop1 != NULL || loop2 != NULL) {
    for (i = 1; i <= totperson; i++) {
      WITH = person[i];

      if (WITH->unknown && WITH->foff != NULL) {
        fprintf(speedfile, "id%7ld\n", WITH->nseq);


          for (j = 0; j < nsystem; j++) {
            allele_count = thislocus[j]->nallele;
            for (a = 0; a < allele_count; a++) {
              for (b = 0; b <  allele_count; b++) {
                fprintf(speedfile, "%3ld%3ld%3ld\n", j + 1, a + 1, b + 1);
	      }
	    }
	  }
      }
    }

    return;
  }
 /* else no loops */
  for (i = 1; i <= totperson; i++) {
    WITH = person[i];

    if (WITH->unknown && WITH->foff != NULL) {
      fprintf(speedfile, "id%7ld\n", WITH->nseq);

      /* if male and sexlinked genotype is stored in second strand */
      if (sexlink && WITH->male) {
        for (j = 0; j < nsystem; j++) {
            allele_count = thislocus[j]->nallele;
          for (b = 0; b <  allele_count; b++) {
            if ( unknown_poss[i][j][0][genenumber[j][0][b] - 1]) {
              fprintf(speedfile, "%3ld%3ld%3ld\n", j + 1, 1, b + 1);
	    }
	  }
	}
      }

      /* if not sexlinked and male, need both strands to determine genotype */
      else {
        for (j = 0; j < nsystem; j++) {
            allele_count = thislocus[j]->nallele;
          for (a = 0; a < allele_count; a++) {
            for (b = 0; b <  allele_count; b++) {
              if ( unknown_poss[i][j][0][genenumber[j][a][b] - 1] ) {
                fprintf(speedfile, "%3ld%3ld%3ld\n", j + 1, a + 1, b + 1);
	      }
	    }
	  }
	}
      }
    }
  }
}   /*writespeed_old*/


/* writes the old speedfile.dat for compatibility.   Written by Dylan. */
static void writespeed_new()
{
  int i, j, a, b;
  thisperson *WITH;
  int allele_count;

  /* for compatibility, when loops exist, speedfile.dat does not eliminate
     anything */
   fprintf(speedfile, "%d\n",ALLELE_SPEED_CONSTANT);
  if (loop1 != NULL || loop2 != NULL) {
    for (i = 1; i <= totperson; i++) {
      WITH = person[i];

      if (WITH->unknown && WITH->foff != NULL) {
        fprintf(speedfile, "id%7ld\n", WITH->nseq);

        /* if male and sexlinked genotype is stored in second strand */
        if (sexlink && WITH->male) {
          for (j = 0; j < nsystem; j++) {
            if ((binary_ == thislocus[j]->which) &&
                (allformat == thislocus[j]->UU.U2.format))
	      allele_count = ped_new_allele_count[currentped][j];
            else
              allele_count = thislocus[j]->nallele;
            for (b = 0; b <  allele_count; b++) {
              fprintf(speedfile, "%3ld%3ld%3ld\n", j + 1, 1, b + 1);
	    }
	  }
	} else {
          for (j = 0; j < nsystem; j++) {
            if ((binary_ == thislocus[j]->which) &&
                (allformat == thislocus[j]->UU.U2.format))
	      allele_count = ped_new_allele_count[currentped][j];
            else
              allele_count = thislocus[j]->nallele;
            for (a = 0; a < allele_count; a++) {
              for (b = 0; b <  allele_count; b++) {
                fprintf(speedfile, "%3ld%3ld%3ld\n", j + 1, a + 1, b + 1);
	      }
	    }
	  }
	}
      }
    }

    return;
  }
 /* else no loops */
  for (i = 1; i <= totperson; i++) {
    WITH = person[i];

    if (WITH->unknown && WITH->foff != NULL) {
      fprintf(speedfile, "id%7ld\n", WITH->nseq);

      /* if male and sexlinked genotype is stored in second strand */
      if (sexlink && WITH->male) {
        for (j = 0; j < nsystem; j++) {
            if ((binary_ == thislocus[j]->which) &&
                (allformat == thislocus[j]->UU.U2.format))
	      allele_count = ped_new_allele_count[currentped][j];
            else
              allele_count = thislocus[j]->nallele;
          for (b = 0; b <  allele_count; b++) {
            if ( unknown_poss[i][j][0][genenumber[j][0][b] - 1]) {
              fprintf(speedfile, "%3ld%3ld%3ld\n", j + 1, 1, b + 1);
	    }
	  }
	}
      }

      /* if not sexlinked and male, need both strands to determine genotype */
      else {
        for (j = 0; j < nsystem; j++) {
            if ((binary_ == thislocus[j]->which) &&
                (allformat == thislocus[j]->UU.U2.format))
	      allele_count = ped_new_allele_count[currentped][j];
            else
              allele_count = thislocus[j]->nallele;
          for (a = 0; a < allele_count; a++) {
            for (b = 0; b <  allele_count; b++) {
              if ( unknown_poss[i][j][0][genenumber[j][a][b] - 1] ) {
                fprintf(speedfile, "%3ld%3ld%3ld\n", j + 1, a + 1, b + 1);
	      }
	    }
	  }
	}
      }
    }
  }
}   /*writespeed_new*/


static void write_loopfile_old(ped)
int ped;
{
  int ind, locus, loop, vect, geno;  /* iterators */
  int geno_count;


  fprintf(loopfile, "Pedigree: %d\n", ped);

  fprintf(loopfile, "fewer_vects_size: %d\n", fewer_vects_size);

  fprintf(loopfile, "num_loops_considered: %d\n", num_loops_considered);

  /* print num_loop_vectors */
  fprintf(loopfile, "num_loop_vectors:\n");
  for (locus = 0; locus < nsystem; locus++) {
    fprintf(loopfile, "\t%d : %d\n", locus, num_loop_vectors[locus]);
  }

  /* print loop_vectors */
  fprintf(loopfile, "loop_vectors:\n");
  for (locus = 0; locus < nsystem; locus++) {
    fprintf(loopfile, "\tL : %d\n", locus);
    for (vect = 0; vect < num_loop_vectors[locus]; vect++) {
      fprintf(loopfile, "\t\t%d :", vect);
      for (loop = 0; loop < num_loops_considered; loop++) {
        fprintf(loopfile, " %d", loop_vectors[locus][vect][loop]);
      }
      fprintf(loopfile, "\n");
    }
  }

  /* print unknown_poss */
  fprintf(loopfile, "unknown_poss:\n");

  for (ind = 1; ind <= totperson; ind++) {

    /* DYLAN -- second condition causes us to lose accuracy */
    if (person[ind]->unknown && person[ind]->foff != NULL) {
      fprintf(loopfile, "id: %d\n", person[ind]->nseq);
      for (locus = 0; locus < nsystem; locus++) {
        fprintf(loopfile, "\tL: %d\n", locus);
        if (!(person[ind]->thisunknown[locus]))
          fprintf(loopfile, "-\n"); /* to facilitate reading it back in */
        else
          fprintf(loopfile, "+\n"); /* to facilitate reading it back in */
        for (vect = 0; vect < num_loop_vectors[locus]; vect++) {
          fprintf(loopfile, "\t\t%d :", vect);
            if (sexlink && (person[ind]->male))
	      geno_count = thislocus[locus]->mgeno;
            else
	      geno_count = thislocus[locus]->fgeno;
          for (geno = 0; geno < geno_count; geno++) {
            if (unknown_poss[ind][locus][vect][geno]) {
              fprintf(loopfile, " %d", geno);
	    }
	  }
          fprintf(loopfile, "\n");
	}

      }
    }
  }
} /* write_loopfile_old */

static void write_loopfile_new(ped)
int ped;
{
  int ind, locus, loop, vect, geno;  /* iterators */
  int geno_count;


  fprintf(loopfile, "Pedigree: %d\n", ped);

  fprintf(loopfile, "fewer_vects_size: %d\n", fewer_vects_size);

  fprintf(loopfile, "num_loops_considered: %d\n", num_loops_considered);

  /* print num_loop_vectors */
  fprintf(loopfile, "num_loop_vectors:\n");
  for (locus = 0; locus < nsystem; locus++) {
    fprintf(loopfile, "\t%d : %d\n", locus, num_loop_vectors[locus]);
  }

  /* print loop_vectors */
  fprintf(loopfile, "loop_vectors:\n");
  for (locus = 0; locus < nsystem; locus++) {
    fprintf(loopfile, "\tL : %d\n", locus);
    for (vect = 0; vect < num_loop_vectors[locus]; vect++) {
      fprintf(loopfile, "\t\t%d :", vect);
      for (loop = 0; loop < num_loops_considered; loop++) {
        fprintf(loopfile, " %d", loop_vectors[locus][vect][loop]);
      }
      fprintf(loopfile, "\n");
    }
  }

  /* print unknown_poss */
  fprintf(loopfile, "unknown_poss:\n");

  for (ind = 1; ind <= totperson; ind++) {

    /* DYLAN -- second condition causes us to lose accuracy */
    if (person[ind]->unknown && person[ind]->foff != NULL) {
      fprintf(loopfile, "id: %d\n", person[ind]->nseq);
      for (locus = 0; locus < nsystem; locus++) {
        fprintf(loopfile, "\tL: %d\n", locus);
        if (!(person[ind]->thisunknown[locus]))
          fprintf(loopfile, "-\n"); /* to facilitate reading it back in */
        else
          fprintf(loopfile, "+\n"); /* to facilitate reading it back in */
        for (vect = 0; vect < num_loop_vectors[locus]; vect++) {
          fprintf(loopfile, "\t\t%d :", vect);
            if ((binary_ == thislocus[locus]->which) &&
                (allformat == thislocus[locus]->UU.U2.format))
              if (sexlink && (person[ind]->male))
                geno_count = ped_new_allele_count[currentped][locus];
              else
		geno_count = 
		  ped_new_allele_count[currentped][locus] * 
		    ( 1 +   ped_new_allele_count[currentped][locus])/2;
            else
              if (sexlink && (person[ind]->male))
		geno_count = thislocus[locus]->mgeno;
              else
		geno_count = thislocus[locus]->fgeno;
          for (geno = 0; geno < geno_count; geno++) {
            if (unknown_poss[ind][locus][vect][geno]) {
              fprintf(loopfile, " %d", geno);
	    }
	  }
          fprintf(loopfile, "\n");
	}

      }
    }
  }
} /* write_loopfile_new */

#else

static void writespeed_old()
{
  long i, j, a_, b_;
  thisperson *WITH;
  information *WITH1;
  int allele_count;

  for (i = 1; i <= totperson; i++) {
    WITH = person[i];
    if (WITH->unknown && WITH->foff != NULL) {
      WITH1 = WITH->store;
      fprintf(speedfile, "id%7ld\n", WITH->nseq);
      for (j = 1; j <= nsystem; j++) {
            allele_count = thislocus[j - 1]->nallele;
        for (a_ = 1; a_ <= allele_count; a_++) {
          for (b_ = 1; b_ <= allele_count; b_++) {
            if (WITH1->possible[j - 1][a_ - 1][b_ - 1])
              fprintf(speedfile, "%3ld%3ld%3ld\n", j, a_, b_);
	  }
	}
      }
    }
  }
}  /*writespeed_old*/

static void writespeed_new()
{
  long i, j, a_, b_;
  thisperson *WITH;
  information *WITH1;
  int allele_count;

   fprintf(speedfile, "%d\n",ALLELE_SPEED_CONSTANT);
  for (i = 1; i <= totperson; i++) {
    WITH = person[i];
    if (WITH->unknown && WITH->foff != NULL) {
      WITH1 = WITH->store;
      fprintf(speedfile, "id%7ld\n", WITH->nseq);
      for (j = 1; j <= nsystem; j++) {
            if ((binary_ == thislocus[j - 1]->which) &&
                (allformat == thislocus[j - 1]->UU.U2.format))
	      allele_count = ped_new_allele_count[currentped][j - 1];
            else
              allele_count = thislocus[j - 1]->nallele;
        for (a_ = 1; a_ <= allele_count; a_++) {
          for (b_ = 1; b_ <= allele_count; b_++) {
            if (WITH1->possible[j - 1][a_ - 1][b_ - 1])
              fprintf(speedfile, "%3ld%3ld%3ld\n", j, a_, b_);
	  }
	}
      }
    }
  }
}  /*writespeed_new*/

#endif

static void writeped_old()
{
  long i, j, k, a_, b_;
  thisperson *WITH;
  phenotype *WITH1;
  locusvalues *WITH2;
  long FORLIM2;

  for (i = 1; i <= totperson; i++) {
    WITH = person[i];
    fprintf(ipedfile, "%7ld%5ld%5ld%5ld%5ld%5ld",
	    WITH->oldped, WITH->id, WITH->paid, WITH->maid, WITH->offid,
	    WITH->npaid);
    fprintf(ipedfile, "%5ld%3ld%3ld ", WITH->nmaid, WITH->sex, WITH->profield);
    for (j = 1; j <= nsystem; j++) {
      WITH1 = WITH->phen[j - 1];
      WITH2 = thislocus[j - 1];
      if (WITH2->which == binary_) {
	if (WITH2->UU.U2.format == binformat) {
	  FORLIM2 = WITH2->UU.U2.nfactor;
	  for (k = 1; k <= FORLIM2; k++) {
	    if ((unsigned long)k < (sizeof(long) * 8)
		&& ((1L << k) & WITH1->phenf) != 0)
	      fprintf(ipedfile, " 1");
	    else
	      fprintf(ipedfile, " 0");
	  }
	} else {
	  a_ = 0;
	  b_ = 0;
          a_ = WITH1->allele1;
          b_ = WITH1->allele2;
	  fprintf(ipedfile, "%3ld%3ld", a_, b_);
	}
      } else if (WITH2->which == quantitative) {
	if (!sexlink || !WITH->male) {
	  FORLIM2 = WITH2->UU.ntrait;
	  for (k = 0; k < FORLIM2; k++)
	    fprintf(ipedfile, " %9.4f", WITH1->x[k]);
	} else {
	  FORLIM2 = WITH2->UU.ntrait;
	  for (k = 1; k <= FORLIM2; k++)
	    fprintf(ipedfile, " %9ld", WITH1->aff);
	}
      } else {
	fprintf(ipedfile, "%2ld", WITH1->aff);
	if (WITH2->UU.U0.nclass != 1)
	  fprintf(ipedfile, "%4ld", WITH1->liability);
      }
      if (j != nsystem)
	putc(' ', ipedfile);
    }
    putc('\n', ipedfile);
  }
}  /*writeped_old*/

static void writeped_new()
{
  long i, j, k, a_, b_;
  thisperson *WITH;
  phenotype *WITH1;
  locusvalues *WITH2;
  long FORLIM2;

  for (i = 1; i <= totperson; i++) {
    WITH = person[i];
    fprintf(ipedfile, "%7ld%5ld%5ld%5ld%5ld%5ld",
	    WITH->oldped, WITH->id, WITH->paid, WITH->maid, WITH->offid,
	    WITH->npaid);
    fprintf(ipedfile, "%5ld%3ld%3ld ", WITH->nmaid, WITH->sex, WITH->profield);
    for (j = 1; j <= nsystem; j++) {
      WITH1 = WITH->phen[j - 1];
      WITH2 = thislocus[j - 1];
      if (WITH2->which == binary_) {
	if (WITH2->UU.U2.format == binformat) {
	  FORLIM2 = WITH2->UU.U2.nfactor;
	  for (k = 1; k <= FORLIM2; k++) {
	    if ((unsigned long)k < (sizeof(long) * 8)
		&& ((1L << k) & WITH1->phenf) != 0)
	      fprintf(ipedfile, " 1");
	    else
	      fprintf(ipedfile, " 0");
	  }
	} else {
	  a_ = 0;
	  b_ = 0;
          a_ = ped_loc_all[currentped][j - 1][WITH1->allele1].old_allele;
          b_ = ped_loc_all[currentped][j - 1][WITH1->allele2].old_allele;
	  fprintf(ipedfile, "%3ld%3ld", a_, b_);
	}
      } else if (WITH2->which == quantitative) {
	if (!sexlink || !WITH->male) {
	  FORLIM2 = WITH2->UU.ntrait;
	  for (k = 0; k < FORLIM2; k++)
	    fprintf(ipedfile, " %9.4f", WITH1->x[k]);
	} else {
	  FORLIM2 = WITH2->UU.ntrait;
	  for (k = 1; k <= FORLIM2; k++)
	    fprintf(ipedfile, " %9ld", WITH1->aff);
	}
      } else {
	fprintf(ipedfile, "%2ld", WITH1->aff);
	if (WITH2->UU.U0.nclass != 1)
	  fprintf(ipedfile, "%4ld", WITH1->liability);
      }
      if (j != nsystem)
	putc(' ', ipedfile);
    }
    putc('\n', ipedfile);
  }
}  /*writeped_new*/

/*
   This procedure infers genotypes that can be known and resets
   the 'unknown' flags.

   Modified by Dylan in late 1994 to use the data structure
   'unknown_poss' rather than the 'possible' field of each person.
*/
static void infer_old()
{
  long i, j, k, l, kposs, lposs, count, pacount, macount;
  boolean someknown;
  thisperson *WITH;
  information *WITH1;
  locusvalues *WITH2;
  long FORLIM2, FORLIM3;
  int loop_vect;
  int skip_locus[maxlocus];

  /*Replace by homozygotes if all unknown in a pedigree*/
  if (!makehomozygous)   /*change - added*/
    return;
  for (j = 0; j < nsystem; j++) {
    WITH2 = thislocus[j];
    skip_locus[j] = 0;
    if ((WITH2->which == binary_) && (WITH2->UU.U2.format == allformat))
    {   /*change - 'format=3' added*/
      someknown = false;
      for (i = 1; i <= totperson; i++) {
        if (person[i]->phen[j]->phenf != 0)
          someknown = true;
      }
      if (!someknown) {
        skip_locus[j] = 1;
        if (WITH2->onefreq > 0.0) {
          for (i = 1; i <= totperson; i++) {
	    person[i]->phen[j]->phenf = 1;
	    person[i]->phen[j]->allele1 = 1;
	    person[i]->phen[j]->allele2 = 1;
	  }
	}
      }
    }
  }
  for (i = 1; i <= totperson; i++) {
    if (person[i]->unknown) {
      WITH = person[i];
#if !LOOPSPEED
      WITH1 = WITH->store;
#endif
      for (j = 0; j < nsystem; j++) {
        if ((thislocus[j]->which == binary_) && (!skip_locus[j])) {
          if (WITH->phen[j]->phenf == 0) {
            WITH2 = thislocus[j];
            count = 0;
#if LOOPSPEED
            /* need to deal with case when only one strand carried genotype */
            if (sexlink && WITH->male) {
              for (l = 1; l <= thislocus[j]->nallele; l++) {
                for (loop_vect = 0;
                     loop_vect < num_loop_vectors[j];
                     loop_vect++) {
                  if (unknown_poss[i][j][loop_vect][genenumber[j][0][l-1]-1]){
                    kposs = 1;
                    lposs = l;
                    count++;
}
		      }
		}

	      } else {
#endif
              FORLIM2 = WITH2->nallele;
              for (k = 1; k <= FORLIM2; k++) {
                FORLIM3 = WITH2->nallele;
                for (l = k; l <= FORLIM3; l++) {
#if LOOPSPEED
                  /* check all loopbreaker vectors */
                  for (loop_vect = 0;
                       loop_vect < num_loop_vectors[j];
                       loop_vect++) {
                    if (unknown_poss[i][j][loop_vect][genenumber[j][k-1][l-1] -1]){
#else
                    if (WITH1->possible[j][k - 1][l - 1]) {
#endif
                      kposs = k;
                      lposs = l;
                      count++;
		    }
                  }
                }
#if LOOPSPEED   /* extra { match "} else {" and for statement */
	      }
	    }
#endif

            if (count == 1) {
              if (sexlink && WITH->male) {
                if (allformat == thislocus[j]->UU.U2.format) {
                  if(lposs)
		    WITH->phen[j]->phenf = 1;
                  else
                    WITH->phen[j]->phenf = 0;
		  WITH->phen[j]->allele1 = (WITH->phen[j]->allele2 = lposs);

                }
                else
		  WITH->phen[j]->phenf = WITH2->UU.U2.allele[lposs - 1];
	      }
              else {
                if (allformat == thislocus[j]->UU.U2.format) {
                  WITH->phen[j]->phenf = 
                  (kposs || lposs);
                    WITH->phen[j]->allele1 = kposs;
                    WITH->phen[j]->allele2 = lposs;
		}
                else
                  WITH->phen[j]->phenf = 
                    WITH2->UU.U2.allele[kposs - 1] |
                    WITH2->UU.U2.allele[lposs - 1];
              }
	    }
          }
        }
      } /* for each locus */

      count = 0;
      for (j = 0; j < nsystem; j++) {
        if (thislocus[j]->which != binary_)
          count++;
        else if (WITH->phen[j]->phenf == 0)
          count++;
      }
      WITH->unknown = (count != 0);
    }
  }

  /*Infer children when parents are homozygotes*/
  for (i = 1; i <= totperson; i++) {
    if (person[i]->foff == NULL) {
      WITH = person[i];
      for (j = 0; j < nsystem; j++) {
        WITH2 = thislocus[j];
        if (WITH->phen[j]->which == binary_) {
          if (WITH->phen[j]->phenf == 0) {
            if (WITH->pa != NULL) {
              pacount = 0;
              macount = 0;
              if(binformat == WITH2->UU.U2.format) {
		FORLIM2 = thislocus[j]->nallele;
		for (k = 1; k <= FORLIM2; k++) {
		  if ((WITH2->UU.U2.allele[k - 1] &
		       (~WITH->pa->phen[j]->phenf)) == 0) {
		    kposs = k;
		    pacount++;
		  }
		}
		FORLIM2 = thislocus[j]->nallele;
		for (l = 1; l <= FORLIM2; l++) {
		  if ((WITH2->UU.U2.allele[l - 1] &
		       (~WITH->ma->phen[j]->phenf)) == 0) {
		    lposs = l;
		    macount++;
		  }
		}
              }
              else  {  /*numbered alleles*/
                if (WITH->pa->phen[j]->allele1 > 0) {
                  kposs = WITH->pa->phen[j]->allele1;
                  pacount++;
                }
                if ((WITH->pa->phen[j]->allele1 != WITH->pa->phen[j]->allele2)
                     && (WITH->pa->phen[j]->allele2 > 0)) {
                  kposs = WITH->pa->phen[j]->allele2;
                  pacount++;
                }
                if (WITH->ma->phen[j]->allele1 > 0) {
                  lposs = WITH->ma->phen[j]->allele1;
                  macount++;
                }
                if ((WITH->ma->phen[j]->allele1 != WITH->ma->phen[j]->allele2)
                     && (WITH->ma->phen[j]->allele2 > 0)) {
                  lposs = WITH->ma->phen[j]->allele2;
                  macount++;
                }
              }
              if (macount == 1 && pacount == 1 && !(WITH->male && sexlink))
                if(allformat == WITH2->UU.U2.format) {
		  if (kposs <=lposs) {
		    WITH->phen[j]->allele1 = kposs;
		    WITH->phen[j]->allele2 = lposs;
		  }
		  else {
		    WITH->phen[j]->allele1 = lposs;
		    WITH->phen[j]->allele2 = kposs;
		  }
		  WITH->phen[j]->phenf = ((kposs > 0)  || (lposs > 0));
		}
		else
		  WITH->phen[j]->phenf = WITH2->UU.U2.allele[kposs - 1] |
                                          WITH2->UU.U2.allele[lposs - 1];
              else if (macount == 1 && WITH->male && sexlink)
                if (allformat == WITH2->UU.U2.format) {
                  WITH->phen[j]->allele1 = lposs;
		  WITH->phen[j]->phenf = (lposs > 0);
                }
                else
		  WITH->phen[j]->phenf = WITH2->UU.U2.allele[lposs - 1];

	    }
	  }
	}
      }
    }
  }
}  /*infer_old*/

static void infer_new()
{
  long i, j, k, l, kposs, lposs, count, pacount, macount;
  boolean someknown;
  thisperson *WITH;
  information *WITH1;
  locusvalues *WITH2;
  long FORLIM2, FORLIM3;
  int loop_vect;
  int skip_locus[maxlocus];

  /*Replace by homozygotes if all unknown in a pedigree*/
  if (!makehomozygous)   /*change - added*/
    return;
  for (j = 0; j < nsystem; j++) {
    WITH2 = thislocus[j];
    skip_locus[j] = 0;
    if ((WITH2->which == binary_) && (WITH2->UU.U2.format == allformat))
    {   /*change - 'format=3' added*/
      someknown = false;
      for (i = 1; i <= totperson; i++) {
        if (person[i]->phen[j]->phenf != 0)
          someknown = true;
      }
      if (!someknown) {
        skip_locus[j] = 1;
        if (WITH2->onefreq > 0.0) {
          for (i = 1; i <= totperson; i++) {
	    person[i]->phen[j]->phenf = 1;
	    person[i]->phen[j]->allele1 = 1;
	    person[i]->phen[j]->allele2 = 1;
	  }
	  ped_loc_all[currentped][j][1].old_allele = 1;   
	}
      }
    }
  }
  for (i = 1; i <= totperson; i++) {
    if (person[i]->unknown) {
      WITH = person[i];
#if !LOOPSPEED
      WITH1 = WITH->store;
#endif
      for (j = 0; j < nsystem; j++) {
        if ((thislocus[j]->which == binary_) && (!skip_locus[j])) {
          if (WITH->phen[j]->phenf == 0) {
            WITH2 = thislocus[j];
            count = 0;
#if LOOPSPEED
            /* need to deal with case when only one strand carried genotype */
            if (sexlink && WITH->male) {
              for (l = 1; l <= thislocus[j]->nallele; l++) {
                for (loop_vect = 0;
                     loop_vect < num_loop_vectors[j];
                     loop_vect++) {
                  if (unknown_poss[i][j][loop_vect][genenumber[j][0][l-1]-1]){
                    kposs = 1;
                    lposs = l;
                    count++;
}
		      }
		}

	      } else {
#endif
              FORLIM2 = WITH2->nallele;
              for (k = 1; k <= FORLIM2; k++) {
                FORLIM3 = WITH2->nallele;
                for (l = k; l <= FORLIM3; l++) {
#if LOOPSPEED
                  /* check all loopbreaker vectors */
                  for (loop_vect = 0;
                       loop_vect < num_loop_vectors[j];
                       loop_vect++) {
                    if (unknown_poss[i][j][loop_vect][genenumber[j][k-1][l-1] -1]){
#else
                    if (WITH1->possible[j][k - 1][l - 1]) {
#endif
                      kposs = k;
                      lposs = l;
                      count++;
		    }
                  }
                }
#if LOOPSPEED   /* extra { match "} else {" and for statement */
	      }
	    }
#endif

            if (count == 1) {
              if (sexlink && WITH->male) {
                if (allformat == thislocus[j]->UU.U2.format) {
                  if (ped_loc_all[currentped][currentlocus][lposs].old_allele)
		    WITH->phen[j]->phenf = 1;
                  else
                    WITH->phen[j]->phenf = 0;
		  WITH->phen[j]->allele1 = (WITH->phen[j]->allele2 = lposs);

                }
                else
		  WITH->phen[j]->phenf = WITH2->UU.U2.allele[lposs - 1];
	      }
              else {
                if (allformat == thislocus[j]->UU.U2.format) {
                  WITH->phen[j]->phenf = 
                    (ped_loc_all[currentped][currentlocus][kposs].old_allele ||
                     ped_loc_all[currentped][currentlocus][lposs].old_allele);
                    WITH->phen[j]->allele1 = kposs;
                    WITH->phen[j]->allele2 = lposs;
		}
                else
                  WITH->phen[j]->phenf = 
                    WITH2->UU.U2.allele[kposs - 1] |
                    WITH2->UU.U2.allele[lposs - 1];
              }
	    }
          }
        }
      } /* for each locus */

      count = 0;
      for (j = 0; j < nsystem; j++) {
        if (thislocus[j]->which != binary_)
          count++;
        else if (WITH->phen[j]->phenf == 0)
          count++;
      }
      WITH->unknown = (count != 0);
    }
  }

  /*Infer children when parents are homozygotes*/
  for (i = 1; i <= totperson; i++) {
    if (person[i]->foff == NULL) {
      WITH = person[i];
      for (j = 0; j < nsystem; j++) {
        WITH2 = thislocus[j];
        if (WITH->phen[j]->which == binary_) {
          if (WITH->phen[j]->phenf == 0) {
            if (WITH->pa != NULL) {
              pacount = 0;
              macount = 0;
              if(binformat == WITH2->UU.U2.format) {
		FORLIM2 = thislocus[j]->nallele;
		for (k = 1; k <= FORLIM2; k++) {
		  if ((WITH2->UU.U2.allele[k - 1] &
		       (~WITH->pa->phen[j]->phenf)) == 0) {
		    kposs = k;
		    pacount++;
		  }
		}
		FORLIM2 = thislocus[j]->nallele;
		for (l = 1; l <= FORLIM2; l++) {
		  if ((WITH2->UU.U2.allele[l - 1] &
		       (~WITH->ma->phen[j]->phenf)) == 0) {
		    lposs = l;
		    macount++;
		  }
		}
              }
              else  {  /*numbered alleles*/
                if (WITH->pa->phen[j]->allele1 > 0) {
                  kposs = WITH->pa->phen[j]->allele1;
                  pacount++;
                }
                if ((WITH->pa->phen[j]->allele1 != WITH->pa->phen[j]->allele2)
                     && (WITH->pa->phen[j]->allele2 > 0)) {
                  kposs = WITH->pa->phen[j]->allele2;
                  pacount++;
                }
                if (WITH->ma->phen[j]->allele1 > 0) {
                  lposs = WITH->ma->phen[j]->allele1;
                  macount++;
                }
                if ((WITH->ma->phen[j]->allele1 != WITH->ma->phen[j]->allele2)
                     && (WITH->ma->phen[j]->allele2 > 0)) {
                  lposs = WITH->ma->phen[j]->allele2;
                  macount++;
                }
              }
              if (macount == 1 && pacount == 1 && !(WITH->male && sexlink))
                if(allformat == WITH2->UU.U2.format) {
		  if (kposs <=lposs) {
		    WITH->phen[j]->allele1 = kposs;
		    WITH->phen[j]->allele2 = lposs;
		  }
		  else {
		    WITH->phen[j]->allele1 = lposs;
		    WITH->phen[j]->allele2 = kposs;
		  }
		  WITH->phen[j]->phenf = ((kposs > 0)  || (lposs > 0));
		}
		else
		  WITH->phen[j]->phenf = WITH2->UU.U2.allele[kposs - 1] |
                                          WITH2->UU.U2.allele[lposs - 1];
              else if (macount == 1 && WITH->male && sexlink)
                if (allformat == WITH2->UU.U2.format) {
                  WITH->phen[j]->allele1 = lposs;
		  WITH->phen[j]->phenf = (lposs > 0);
                }
                else
		  WITH->phen[j]->phenf = WITH2->UU.U2.allele[lposs - 1];

	    }
	  }
	}
      }
    }
  }
}  /*infer_new*/

static void getunknown_old()
{

  long i, j, n, ahap, bhap;
  thisperson *WITH;
  information *WITH1;
  locusvalues *WITH2;
  long FORLIM2;

  for (i = 1; i <= totperson; i++) {
    person[i]->unknown = false;
    for (j = 0; j < nsystem; j++)
      person[i]->thisunknown[j] = false;
  }

  /*
     Next bit modified in late 1994 by Dylan.  The if-then-else statements
     used to be nested as if no { were used.  Thus unknown only got set
     true when the locus was binary.   */

  for (i = 1; i <= totperson; i++) {
    WITH = person[i];
    for (j = 0; j < nsystem; j++) {
      if (thislocus[j]->which == binary_) {
        if (WITH->phen[j]->phenf == 0)
          WITH->thisunknown[j] = true;
      }
      else if (thislocus[j]->which == quantitative) {
        if (WITH->phen[j]->x[0] == missval)
          WITH->thisunknown[j] = true;
      }
      else if (WITH->phen[j]->aff == missaff) {
        WITH->thisunknown[j] = true;
      }
      if (WITH->thisunknown[j]) {
        WITH->unknown = true;
      }
    }
  }

#if !LOOPSPEED
  for (i = 1; i <= totperson; i++) {
    WITH = person[i];
    if (WITH->unknown) {
      WITH->store = (information *)Malloc(sizeof(information));
      /* malloc check added by Dylan */
      if (WITH->store == NULL)
        malloc_err("store field");
      WITH1 = WITH->store;
      for (n = 0; n < nsystem; n++) {
        WITH2 = thislocus[n];
        FORLIM2 = WITH2->nallele;
        for (ahap = 0; ahap < FORLIM2; ahap++) {
          for (bhap = 0; bhap < FORLIM2; bhap++)
            WITH1->possible[n][ahap][bhap] = true;
	}
      }
    }
  }
#endif
}  /*getunknown_old*/

static void getunknown_new()
{

  long i, j, n, ahap, bhap;
  thisperson *WITH;
  information *WITH1;
  locusvalues *WITH2;
  long FORLIM2;

  for (i = 1; i <= totperson; i++) {
    person[i]->unknown = false;
    for (j = 0; j < nsystem; j++)
      person[i]->thisunknown[j] = false;
  }

  /*
     Next bit modified in late 1994 by Dylan.  The if-then-else statements
     used to be nested as if no { were used.  Thus unknown only got set
     true when the locus was binary.   */

  for (i = 1; i <= totperson; i++) {
    WITH = person[i];
    for (j = 0; j < nsystem; j++) {
      if (thislocus[j]->which == binary_) {
        if (WITH->phen[j]->phenf == 0)
          WITH->thisunknown[j] = true;
      }
      else if (thislocus[j]->which == quantitative) {
        if (WITH->phen[j]->x[0] == missval)
          WITH->thisunknown[j] = true;
      }
      else if (WITH->phen[j]->aff == missaff) {
        WITH->thisunknown[j] = true;
      }
      if (WITH->thisunknown[j]) {
        WITH->unknown = true;
      }
    }
  }

#if !LOOPSPEED
  for (i = 1; i <= totperson; i++) {
    WITH = person[i];
    if (WITH->unknown) {
      WITH->store = (information *)Malloc(sizeof(information));
      /* malloc check added by Dylan */
      if (WITH->store == NULL)
        malloc_err("store field");
      WITH1 = WITH->store;
      for (n = 0; n < nsystem; n++) {
        WITH2 = thislocus[n];
        if ((binary_ == thislocus[n]->which) &&
            (allformat == thislocus[n]->UU.U2.format))
          FORLIM2 = ped_new_allele_count[currentped][n];
        else
          FORLIM2 = WITH2->nallele;
        for (ahap = 0; ahap < FORLIM2; ahap++) {
          for (bhap = 0; bhap < FORLIM2; bhap++)
            WITH1->possible[n][ahap][bhap] = true;
	}
      }
    }
  }
#endif
}  /*getunknown_new*/

static void getlocation(thislocus)
locusvalues *thislocus;
{
  long ahap, bhap, here, FORLIM, FORLIM1;

  here = 0;
  FORLIM = thislocus->nallele;
  for (ahap = 1; ahap <= FORLIM; ahap++) {
    FORLIM1 = thislocus->nallele;
    for (bhap = ahap; bhap <= FORLIM1; bhap++) {
      here++;
#if LOOPSPEED
      /*
         WARNING!!!

         If you change the way the genotypes are computed from the
         alleles you need to make changes to procedures:
             translate_loop_vector()
             fac_()
             aff_()
             facmale_()
             affmale_()
             quanmale_()
      */
      /* DYLAN -- maybe saving seghaps instead would have been better */
      genenumber[whichsys - 1][ahap - 1][bhap - 1] = here;
      genenumber[whichsys - 1][bhap - 1][ahap - 1] = here;
#else
      genenumber[ahap - 1][bhap - 1] = here;
      genenumber[bhap - 1][ahap - 1] = here;
#endif
      seghap[here - 1][0] = ahap;
      seghap[here - 1][1] = bhap;
    }
  }
}  /*getlocation*/

/*  Local variables for getphenotype: */
struct LOC_getphenotype {
  thisperson **p;
} ;

void readbin(phen, ourlocus, LINK)
phenotype **phen;
locusvalues *ourlocus;
struct LOC_getphenotype *LINK;
{
  long i, j;
  phenotype *WITH1;
  long FORLIM;

  WITH1 = *phen;
  WITH1->which = binary_;
  WITH1->phenf = 0;

  FORLIM = ourlocus->UU.U2.nfactor;
  for (i = 1; i <= FORLIM; i++) {
    fscanf(pedfile, "%ld", &j);
    if (j != 0 && j != 1)
      inputerror(14L, (*LINK->p)->id, j);
    if (j == 1)
      WITH1->phenf = ((long)WITH1->phenf) | (1L << ((int)i));
  }
}


void readnumber_old(phen, LINK)
phenotype **phen;
struct LOC_getphenotype *LINK;
{
  long j, k;
  phenotype *WITH;

  int allidx;

  WITH = *phen;
  WITH->which = binary_;
  WITH->phenf = 0;
  fscanf(pedfile, "%ld%ld", &j, &k);
  if (j > maxall)
    inputerror(16L, (*LINK->p)->id, j);
  if ((j < 0) || (j > thislocus[currentlocus]->maxallele))
    inputerror(17L, (*LINK->p)->id, j);
  if (k > maxall)
    inputerror(16L, (*LINK->p)->id, k);
  if ((k < 0) || (k > thislocus[currentlocus]->maxallele ))
    inputerror(17L, (*LINK->p)->id, k);
  if ((j == 0 || k == 0) && j != k) {
    inputerror(44L, (*LINK->p)->id, j);
    return;
  }
  if (j != 0)
    WITH->phenf = 1;
  if (k != 0)
    WITH->phenf = 1;
  if (j <= k) {
    WITH->allele1 = j;
    WITH->allele2 = k;
  }
  else {
    WITH->allele1 = k;
    WITH->allele2 = j;
  }
}

void readnumber_new(phen, LINK)
phenotype **phen;
struct LOC_getphenotype *LINK;
{
  long j, k;
  phenotype *WITH;

  int allidx;

  WITH = *phen;
  WITH->which = binary_;
  WITH->phenf = 0;
  fscanf(pedfile, "%ld%ld", &j, &k);
  if (j > maxall)
    inputerror(16L, (*LINK->p)->id, j);
  if ((j < 0) || (j > thislocus[currentlocus]->maxallele))
    inputerror(17L, (*LINK->p)->id, j);
  if (k > maxall)
    inputerror(16L, (*LINK->p)->id, k);
  if ((k < 0) || (k > thislocus[currentlocus]->maxallele ))
    inputerror(17L, (*LINK->p)->id, k);
  if ((j == 0 || k == 0) && j != k) {
    inputerror(44L, (*LINK->p)->id, j);
    return;
  }
  if (j != 0)
    WITH->phenf = 1;
  if (k != 0)
    WITH->phenf = 1;
  ped_loc_all[currentped][currentlocus][j].present = true;
  ped_loc_all[currentped][currentlocus][k].present = true;
  if (j <= k) {
    WITH->allele1 = j;
    WITH->allele2 = k;
  }
  else {
    WITH->allele1 = k;
    WITH->allele2 = j;
  }
}

void readaff(phen, thislocus, LINK)
phenotype **phen;
locusvalues *thislocus;
struct LOC_getphenotype *LINK;
{
  long thisval;
  phenotype *WITH;

  WITH = *phen;
  WITH->which = affection;
  fscanf(pedfile, "%ld", &thisval);
  if (thisval == missaff)
    WITH->aff = 0;
  else {
    if (thisval == affval)
      WITH->aff = 2;
    else {
      if (thisval != 1)
	inputwarning(3L, (*LINK->p)->id, thisval);
      WITH->aff = 1;
    }
  }
  if (thislocus->UU.U0.nclass == 1)
    WITH->liability = 1;
  else
    fscanf(pedfile, "%ld", &WITH->liability);
  if (WITH->liability > thislocus->UU.U0.nclass)
    inputerror(26L, (*LINK->p)->id, WITH->liability);
  if (WITH->liability <= 0)
    inputerror(27L, (*LINK->p)->id, WITH->liability);
}

 void readquan(phen, thislocus, LINK)
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
    FORLIM = thislocus->UU.ntrait; 
    for (i = 0; i < FORLIM; i++)
      fscanf(pedfile, "%lg", &WITH->x[i]);
    WITH->missing = true;
    FORLIM = thislocus->UU.ntrait;
    for (i = 0; i < FORLIM; i++) {
      if (WITH->x[i] != missval)
	WITH->missing = false;
    }
    return;
  }
  WITH->which = affection;
  fscanf(pedfile, "%lg", &xval);
  if (xval == missval)
    WITH->aff = missaff;
  else {
    if (xval == affall)
      WITH->aff = affall;
    else
      WITH->aff = -11;
  }
  WITH->liability = 1;
  FORLIM = thislocus->UU.ntrait;
  for (i = 2; i <= FORLIM; i++)
    fscanf(pedfile, "%lg", &xval);
}


 void getphenotype_old(p_)
thisperson **p_;
{
  struct LOC_getphenotype V;
  long thisread, system;
  thisperson *WITH;

  V.p = p_;
  WITH = *V.p;
  for (thisread = 1; thisread <= nsystem; thisread++) {
    system = thisread;
    currentlocus = system - 1;
    WITH->phen[system - 1] = NULL;
    WITH->phen[system - 1] = (phenotype *)Malloc(sizeof(phenotype));
   /* malloc check added by Dylan */
    if (WITH->phen[system - 1] == NULL)
      malloc_err("phenotype field");
    switch (thislocus[system - 1]->which) {

    case quantitative:
      readquan(&WITH->phen[system - 1], thislocus[system - 1], &V);
      break;

    case affection:
      readaff(&WITH->phen[system - 1], thislocus[system - 1], &V);
      break;

    case binary_:
      if (thislocus[system - 1]->UU.U2.format == 3)
	readnumber_old(&WITH->phen[system - 1], &V);
      else
	readbin(&WITH->phen[system - 1], thislocus[system - 1], &V);
      break;
    }
  }
}  /*getphenotype_old*/


 void getphenotype_new(p_)
thisperson **p_;
{
  struct LOC_getphenotype V;
  long thisread, system;
  thisperson *WITH;

  V.p = p_;
  WITH = *V.p;
  for (thisread = 1; thisread <= nsystem; thisread++) {
    system = thisread;
    currentlocus = system - 1;
    WITH->phen[system - 1] = NULL;
    WITH->phen[system - 1] = (phenotype *)Malloc(sizeof(phenotype));
   /* malloc check added by Dylan */
    if (WITH->phen[system - 1] == NULL)
      malloc_err("phenotype field");
    switch (thislocus[system - 1]->which) {

    case quantitative:
      readquan(&WITH->phen[system - 1], thislocus[system - 1], &V);
      break;

    case affection:
      readaff(&WITH->phen[system - 1], thislocus[system - 1], &V);
      break;

    case binary_:
      if (thislocus[system - 1]->UU.U2.format == 3)
	readnumber_new(&WITH->phen[system - 1], &V);
      else
	readbin(&WITH->phen[system - 1], thislocus[system - 1], &V);
      break;
    }
  }
}  /*getphenotype_new*/


 void getind(id, seq)
long *id, *seq;
{
  thisperson *WITH;

  *id = 0;
  fscanf(pedfile, "%ld", seq);
  if (*seq == 0)
    return;
  *id = *seq;
  if (*id > maxind)
    inputerror(13L, *id, *id);
  if (person[*id] != NULL)
    return;
  numind++;
  person[*id] = (thisperson *)Malloc(sizeof(thisperson));
  /* malloc check added by Dylan */
  if (person[*id] == NULL)
    malloc_err("person");
  WITH = person[*id];
  WITH->gen = (thisarray *)Malloc(sizeof(thisarray));
   /* malloc check added by Dylan */
  if (WITH->gen == NULL)
    malloc_err("gen field");
  WITH->nseq = *seq + nsequence;
}  /*getind*/


 void multimarriage(p)
thisperson **p;
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


 void getrest()
{
  char whichchr;

  whichchr = ' ';
  while (!(P_eoln(pedfile) || whichchr == '\f')) {
    whichchr = getc(pedfile);
    if (whichchr == '\n')
      whichchr = ' ';
  }
  fscanf(pedfile, "%*[^\n]");
  getc(pedfile);
}  /*getrest*/


/*
   This procedure reads in one pedigree.

   Modified in late 1994 by Dylan to keep information about the
   loopbreakers in 'looppers'.
*/
static void readped_old()
{
  long i, newid, thisone, thisped, tempid;
  thisperson *WITH;
  thisperson *holdloop;

  for (i = 0; i <= maxind; i++)
    person[i] = NULL;
  totperson = 0;
  loop1 = NULL;
  loop2 = NULL;
  proband = NULL;
  thisped = newped;

#if LOOPSPEED
/* DYLAN added looppers init */
  for (i = 0; i < maxloop; i++) {
    looppers[i][0] = NULL;
    looppers[i][1] = NULL;
  }
  num_loops = 0;  /* added by Dylan late 1994 */
#endif


  printf("Ped. %2ld\n", thisped);   /*change - added*/
  while (!P_eof(pedfile) && thisped == newped) {
    totperson++;
    getind(&thisone, &tempid);
    if (proband == NULL)
      proband = person[thisone];
    WITH = person[thisone];
    WITH->id = tempid;
    WITH->oldped = newped;
    getind(&newid, &WITH->paid);
    WITH->pa = person[newid];
    getind(&newid, &WITH->maid);
    WITH->ma = person[newid];
    getind(&newid, &WITH->offid);
    WITH->foff = person[newid];
    getind(&newid, &WITH->npaid);
    WITH->nextpa = person[newid];
    getind(&newid, &WITH->nmaid);
    WITH->nextma = person[newid];
#if !LOOPSPEED
    WITH->store = NULL;
#endif
    WITH->up = false;
    fscanf(pedfile, "%ld", &WITH->sex);
    if (WITH->sex != 1 && WITH->sex != 2)
      inputerror(11L, WITH->id, WITH->sex);
    WITH->male = (WITH->sex == 1);
    fscanf(pedfile, "%ld", &WITH->profield);

    /* DYLAN, check for more than max loop */
    if ((WITH->profield - 1) > maxloop) {
      fprintf(stderr, "\nUNKNOWN: Your pedigree has more loops than allowed by the constant maxloop");
      fprintf(stderr, "\nYou must increase the constant maxloop defined in unknown.c and recompile");
      fprintf(stderr, "\nYou are encouraged to read the loops.ps document distributed with FASTLINK");
      fprintf(stderr, "\nUNKNOWN will exit politely to allow you to correct the problem\n");
      exit(EXIT_FAILURE);
	      }

    if (WITH->profield == 1)
      proband = person[thisone];

#if LOOPSPEED
    /* DYLAN, init looppers and loopX */
    else if (WITH->profield > 1 && WITH->profield - 1 <= maxloop) {
      if (looppers[WITH->profield - 2][1] == NULL)
        looppers[WITH->profield - 2][1] = person[thisone];
      else
        looppers[WITH->profield - 2][0] = person[thisone];
#endif
      if (WITH->profield == 2) {
        if (loop2 == NULL)
          loop2 = person[thisone];
        else
          loop1 = person[thisone];
      }
#if LOOPSPEED
    }
#endif

    getphenotype_old(&person[thisone]);
    getrest();
    if (!P_eof(pedfile))
      fscanf(pedfile, "%ld", &newped);
  } /* while */

  nsequence += totperson;

#if LOOPSPEED
  /* deal with proband in loop */
  if (loop2 != NULL && loop1 == NULL)
    loop1 = proband;
  if (looppers[0][1] != NULL && looppers[0][0] == NULL)
    looppers[0][0] = proband;

  /* make sure looppers is set up right */
  for (i = 0; i < maxloop; i++) {
    if (looppers[i][0] == NULL)
      looppers[i][1] = NULL;
    else {
      num_loops++;
      if (looppers[i][0]->pa == NULL && looppers[i][1]->pa != NULL) {
        holdloop = looppers[i][0];
        looppers[i][0] = looppers[i][1];
        looppers[i][1] = holdloop;
      }
    }
  }
#endif

  for (thisone = 1; thisone <= totperson; thisone++)
    multimarriage(&person[thisone]);
}  /*readped_old*/

/*
   This procedure reads in one pedigree.

   Modified in late 1994 by Dylan to keep information about the
   loopbreakers in 'looppers'.
*/
static void readped_new()
{
  long i, newid, thisone, thisped, tempid;
  thisperson *WITH;
  thisperson *holdloop;

  for (i = 0; i <= maxind; i++)
    person[i] = NULL;
  totperson = 0;
  loop1 = NULL;
  loop2 = NULL;
  proband = NULL;
  thisped = newped;

#if LOOPSPEED
/* DYLAN added looppers init */
  for (i = 0; i < maxloop; i++) {
    looppers[i][0] = NULL;
    looppers[i][1] = NULL;
  }
  num_loops = 0;  /* added by Dylan late 1994 */
#endif


  printf("Ped. %2ld\n", thisped);   /*change - added*/
  while (!P_eof(pedfile) && thisped == newped) {
    totperson++;
    getind(&thisone, &tempid);
    if (proband == NULL)
      proband = person[thisone];
    WITH = person[thisone];
    WITH->id = tempid;
    WITH->oldped = newped;
    getind(&newid, &WITH->paid);
    WITH->pa = person[newid];
    getind(&newid, &WITH->maid);
    WITH->ma = person[newid];
    getind(&newid, &WITH->offid);
    WITH->foff = person[newid];
    getind(&newid, &WITH->npaid);
    WITH->nextpa = person[newid];
    getind(&newid, &WITH->nmaid);
    WITH->nextma = person[newid];
#if !LOOPSPEED
    WITH->store = NULL;
#endif
    WITH->up = false;
    fscanf(pedfile, "%ld", &WITH->sex);
    if (WITH->sex != 1 && WITH->sex != 2)
      inputerror(11L, WITH->id, WITH->sex);
    WITH->male = (WITH->sex == 1);
    fscanf(pedfile, "%ld", &WITH->profield);

    /* DYLAN, check for more than max loop */
    if ((WITH->profield - 1) > maxloop) {
      fprintf(stderr, "\nUNKNOWN: Your pedigree has more loops than allowed by the constant maxloop");
      fprintf(stderr, "\nYou must increase the constant maxloop defined in unknown.c and recompile");
      fprintf(stderr, "\nYou are encouraged to read the loops.ps document distributed with FASTLINK");
      fprintf(stderr, "\nUNKNOWN will exit politely to allow you to correct the problem\n");
      exit(EXIT_FAILURE);
	      }

    if (WITH->profield == 1)
      proband = person[thisone];

#if LOOPSPEED
    /* DYLAN, init looppers and loopX */
    else if (WITH->profield > 1 && WITH->profield - 1 <= maxloop) {
      if (looppers[WITH->profield - 2][1] == NULL)
        looppers[WITH->profield - 2][1] = person[thisone];
      else
        looppers[WITH->profield - 2][0] = person[thisone];
#endif
      if (WITH->profield == 2) {
        if (loop2 == NULL)
          loop2 = person[thisone];
        else
          loop1 = person[thisone];
      }
#if LOOPSPEED
    }
#endif

    getphenotype_new(&person[thisone]);
    getrest();
    if (!P_eof(pedfile))
      fscanf(pedfile, "%ld", &newped);
  } /* while */

  nsequence += totperson;

#if LOOPSPEED
  /* deal with proband in loop */
  if (loop2 != NULL && loop1 == NULL)
    loop1 = proband;
  if (looppers[0][1] != NULL && looppers[0][0] == NULL)
    looppers[0][0] = proband;

  /* make sure looppers is set up right */
  for (i = 0; i < maxloop; i++) {
    if (looppers[i][0] == NULL)
      looppers[i][1] = NULL;
    else {
      num_loops++;
      if (looppers[i][0]->pa == NULL && looppers[i][1]->pa != NULL) {
        holdloop = looppers[i][0];
        looppers[i][0] = looppers[i][1];
        looppers[i][1] = holdloop;
      }
    }
  }
#endif

  for (thisone = 1; thisone <= totperson; thisone++)
    multimarriage(&person[thisone]);
}  /*readped_new*/

/* Local variables for readloci: */
struct LOC_readloci {
  long whichtype;
} ;

/* Local variables for getlocus: */
struct LOC_getlocus {
  struct LOC_readloci *LINK;
  long system;
} ;

void getquan(locus, LINK)
locusvalues **locus;
struct LOC_getlocus *LINK;
{
  long i;
  locusvalues *WITH;

  WITH = *locus;
  fscanf(datafile, "%ld%*[^\n]", &WITH->UU.ntrait);
  getc(datafile);
  if (WITH->UU.ntrait > maxtrait)
    inputerror(31L, LINK->system, WITH->UU.ntrait);
  if (WITH->UU.ntrait <= 0)
    inputerror(32L, LINK->system, WITH->UU.U0.nclass);
  for (i = 1; i <= 3; i++) {
    fscanf(datafile, "%*[^\n]");
    getc(datafile);
  }
}  /*getquan*/


void getpen(locus, LINK)
locusvalues **locus;
struct LOC_getlocus *LINK;
{
  long i, j, k, l;
  locusvalues *WITH;
  long FORLIM, FORLIM1, FORLIM2;
  WITH = *locus;
  fscanf(datafile, "%ld%*[^\n]", &WITH->UU.U0.nclass);
  getc(datafile);
  if (WITH->UU.U0.nclass > maxliab)
    inputerror(28L, LINK->system, WITH->UU.U0.nclass);
  if (WITH->UU.U0.nclass <= 0)
    inputerror(29L, LINK->system, WITH->UU.U0.nclass);
  FORLIM = WITH->UU.U0.nclass;
  for (l = 0; l < FORLIM; l++) {
    FORLIM1 = WITH->nallele;
    for (i = 1; i <= FORLIM1; i++) {
      FORLIM2 = WITH->nallele;
      for (j = i - 1; j < FORLIM2; j++) {
	fscanf(datafile, "%lg", &WITH->UU.U0.pen[i][j][2][l]);
	if ((unsigned)WITH->UU.U0.pen[i][j][2][l] > one)
	  inputerror(30L, LINK->system, LINK->system);
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
      for (i = 0; i < FORLIM1; i++)
	fscanf(datafile, "%lg", &WITH->UU.U0.pen[0][i][2][l]);
      if ((unsigned)WITH->UU.U0.pen[0][j - 1][2][l] > one)
	inputerror(30L, LINK->system, LINK->system);
      FORLIM1 = WITH->nallele;
      for (i = 0; i < FORLIM1; i++)
	WITH->UU.U0.pen[0][i][1][l] = 1.0 - WITH->UU.U0.pen[0][i][2][l];
      fscanf(datafile, "%*[^\n]");
      getc(datafile);
    }
  }
}  /*getpen*/


void getbin(locus, LINK)
locusvalues **locus;
struct LOC_getlocus *LINK;
{
  long i, j, k;
  locusvalues *WITH;
  long FORLIM, FORLIM1;

  WITH = *locus;
  fscanf(datafile, "%ld%*[^\n]", &WITH->UU.U2.nfactor);
  getc(datafile);
  if ((WITH->UU.U2.nfactor > maxfact) ||
      (WITH->UU.U2.nfactor > (sizeof(long) * 8 - 1)))
    inputerror(8L, LINK->system, WITH->UU.U2.nfactor);
  if (WITH->UU.U2.nfactor <= 0)
    inputerror(9L, LINK->system, WITH->UU.U2.nfactor);
  FORLIM = WITH->nallele;
  for (i = 0; i < FORLIM; i++)
    WITH->UU.U2.allele[i] = 0;
  FORLIM = WITH->nallele;
  for (i = 0; i < FORLIM; i++) {
    FORLIM1 = WITH->UU.U2.nfactor;
    for (j = 1; j <= FORLIM1; j++) {
      fscanf(datafile, "%ld", &k);
      if (k == 1)
	WITH->UU.U2.allele[i] = ((long)WITH->UU.U2.allele[i]) | (1L << ((int)j));
    }
  }
  fscanf(datafile, "%*[^\n]");
  getc(datafile);
}  /*getbin*/


void getnumber(locus, LINK)
locusvalues **locus;
struct LOC_getlocus *LINK;
{
  long i;
  locusvalues *WITH;
  long FORLIM;

  WITH = *locus;
  FORLIM = WITH->nallele;
  for (i = 1; i <= FORLIM; i++)
    WITH->UU.U2.allele[i - 1] = 1;
}  /*getnumber*/


 void getlocus(system_, LINK)
long system_;
struct LOC_readloci *LINK;
{
  struct LOC_getlocus V;
  locusvalues *WITH;


  V.LINK = LINK;
  V.system = system_;
  thislocus[V.system - 1] = (locusvalues *)Malloc(sizeof(locusvalues));
  /* malloc check added by Dylan */
  if (thislocus[V.system - 1] == NULL)
    malloc_err("locus");
  WITH = thislocus[V.system - 1];
  fscanf(datafile, "%ld%ld", &LINK->whichtype, &WITH->nallele);
  WITH->maxallele = WITH->nallele;
#if LOOPSPEED
  WITH->fgeno = (WITH->nallele * (WITH->nallele + 1) / 2);
  if (sexlink)
    WITH->mgeno = WITH->nallele;
  else
    WITH->mgeno = WITH->fgeno;
#endif
  switch (LINK->whichtype) {

  case 0:
    WITH->which = quantitative;
    break;

  case 1:
    WITH->which = affection;
    break;

  case 2:
    WITH->which = binary_;
    WITH->UU.U2.format = binformat;
    break;

  case 3:
    WITH->which = binary_;
    WITH->UU.U2.format = allformat;
    break;
  }
  if (!disequi) {
    fscanf(datafile, "%*[^\n]");    
    fscanf(datafile, "%lf", &thislocus[V.system - 1]->onefreq);
  }
  else {
  fscanf(datafile, "%*[^\n]");
  getc(datafile);
  }
  if (!disequi) {
    fscanf(datafile, "%*[^\n]");
    getc(datafile);
  }
  switch (WITH->which) {

  case quantitative:
    getquan(&thislocus[V.system - 1], &V);
    break;

  case affection:
    getpen(&thislocus[V.system - 1], &V);
    break;

  case binary_:
    if (WITH->UU.U2.format == binformat)
      getbin(&thislocus[V.system - 1], &V);
    else
      getnumber(&thislocus[V.system - 1], &V);
    break;
  }
  if (risk && V.system == risksys) {
    fscanf(datafile, "%*[^\n]");
    getc(datafile);
  }
}  /*getlocus*/



static void readloci()
{
  struct LOC_readloci V;
  long coupling, autosomal;
  double mu;
  int i, j;
#if LOOPSPEED
  int order[maxlocus];
#endif

  fscanf(datafile, "%ld%ld%ld%*[^\n]", &nsystem, &risksys, &autosomal);
  getc(datafile);
  if (nsystem > maxlocus)
    inputerror(0L, nsystem, nsystem);
  if (nsystem <= 0)
    inputerror(1L, nsystem, nsystem);
  risk = (risksys != 0);
  sexlink = (autosomal == 1);
  printf("YOU ARE USING LINKAGE (V%s) WITH%3ld-POINT\n", aVersion, nsystem);
  printf("YOU ARE USING FASTLINK (V%s)\n", fVersion);
  if (sexlink)
    printf(" SEXLINKED DATA\n");
  else
    printf(" AUTOSOMAL DATA\n");
  fscanf(datafile, "%ld%lf%lf%ld%*[^\n]", &mutsys, &mu, &mu, &coupling);
  getc(datafile);
  disequi = (coupling == 1);

  fscanf(datafile, "%*[^\n]");
  getc(datafile);
  for (i = 1; i <= nsystem; i++) 
    getlocus(i, &V);
  
  if(disequi && ALLELE_SPEED) {

    fprintf(stderr,"\nYou cannot use the disequilibrium model with ALLELE_SPEED set to 1.");
    fprintf(stderr,"\nChange ALLELE_SPEED to 0 in unknown.c and commendefs.h.");
    fprintf(stderr,"\nRecompile unknown  and the main programs.\n");
    exit(EXIT_FAILURE);
  }


}  /*readloci*/

#if !LOOPSPEED
static void cleanup(p)
thisperson **p;
{
  long i, j;
  thisperson *WITH;
  information *WITH1;
  thisarray *WITH2;

  WITH = *p;
  if (!WITH->unknown)
    return;
  WITH1 = WITH->store;
  WITH2 = WITH->gen;
  if (sexlink && WITH->male) {
    for (i = 0; i < mgeno; i++) {
      if (!WITH2->genarray[i])
	WITH1->possible[whichsys - 1][0][i] = false;
    }
    for (i = 1; i < mgeno; i++) {
      for (j = 0; j < mgeno; j++)
	WITH1->possible[whichsys - 1][i][j] = false;
    }
    return;
  }
  for (i = 0; i < fgeno; i++) {
    if (!WITH2->genarray[i])
      WITH1->possible[whichsys - 1][seghap[i][0] - 1]
	[seghap[i][1] - 1] = false;
    WITH1->possible[whichsys - 1][seghap[i][1] - 1]
      [seghap[i][0] - 1] = WITH1->possible[whichsys - 1][seghap[i][0] - 1]
      [seghap[i][1] - 1];
  }
}  /*cleanup*/
#endif

static void getgene(system, p, phen)
long system;
thisperson *p;
phenotype **phen;
{
  long here, i, j;
  locusvalues *WITH;
  long FORLIM;
  phenotype *WITH1;
  long FORLIM1;

  here = 0;
  WITH = thislocus[system - 1];
  if (sexlink && p->male) {
    FORLIM = WITH->nallele;
    for (i = 1; i <= FORLIM; i++) {
      here++;
      switch (WITH->which) {

      case quantitative:
	WITH1 = phen[system - 1];
	if (i == affall)
	  p->gen->genarray[here - 1] = (WITH1->aff == affall ||
					WITH1->aff == missaff);
	else
	  p->gen->genarray[here - 1] = (WITH1->aff != affall ||
					WITH1->aff == missaff);
	break;

      case affection:
	WITH1 = phen[system - 1];
	p->gen->genarray[here - 1] = (WITH->UU.U0.pen[0][i - 1]
				      [WITH1->aff]
				      [WITH1->liability - 1] > 0.0);
	break;

      case binary_:
	WITH1 = phen[system - 1];
        if (allformat == WITH->UU.U2.format)
          if ((0 == WITH1->allele1) ||
              (i == WITH1->allele1))
	    p->gen->genarray[here - 1] = 1;
          else
	    p->gen->genarray[here - 1] = 0;
        else
	  p->gen->genarray[here - 1] =
	    (WITH1->phenf == WITH->UU.U2.allele[i - 1] ||
	     WITH1->phenf == 0);
	break;
      }
    }
    return;
  }
  FORLIM = WITH->nallele;
  for (i = 1; i <= FORLIM; i++) {
    FORLIM1 = WITH->nallele;
    for (j = i - 1; j < FORLIM1; j++) {
      here++;
      switch (WITH->which) {

      case quantitative:
	p->gen->genarray[here - 1] = true;
	break;

      case affection:
	WITH1 = phen[system - 1];
	p->gen->genarray[here - 1] = (WITH->UU.U0.pen[i][j][WITH1->aff]
				      [WITH1->liability - 1] > 0.0);
	break;

      case binary_:
        WITH1 = phen[system - 1];
        if (allformat == WITH->UU.U2.format)
          if ((0 == WITH1->allele2) ||
              ((i == WITH1->allele1) && ((j + 1) == WITH1->allele2)))
            p->gen->genarray[here - 1] = 1;
          else
            p->gen->genarray[here - 1] = 0;
        else 
	  p->gen->genarray[here - 1] = (WITH1->phenf ==
	    (WITH->UU.U2.allele[i - 1] | WITH->UU.U2.allele[j]) ||
	    WITH1->phenf == 0);
	break;
      }
    }
  }
}  /*getgene*/

#if LOOPSPEED
/*
   This procedure determines the number of loops that will be considered
   on a per vector basis.  It also sets the global variable
   'fewer_vects_size'.

   Written by Dylan in late 1994.
*/
static void set_num_loops_considered()
{
  int num_geno[maxlocus][maxloop];
  int num_vects[maxlocus];
  int max_vects;
  int locus, loop, geno;
  int max_loop_vectors;
  int last_fewer_vects_size;
  int geno_bound;

  /* count the number of genotypes each loop breaker may have */
  for (locus = 0; locus < nsystem; locus++) {
    for (loop = 0; loop < num_loops; loop++) {
      getgene(locus + 1, looppers[loop][0], looppers[loop][0]->phen);
      num_geno[locus][loop] = 0;
      if ((!sexlink) || (!looppers[loop][0]->male))
        geno_bound = thislocus[locus]->fgeno;
      else
        geno_bound = thislocus[locus]->mgeno;
      for (geno = 0; geno < geno_bound; geno++) {
        if (looppers[loop][0]->gen->genarray[geno] == true ) {
          num_geno[locus][loop]++;
	}
      }
    }
  }
  /* find the number of loops to consider */
  num_loops_considered = 1;
  for (locus = 0; locus < nsystem; locus++) {
    num_vects[locus] = num_geno[locus][0];
  }
  fewer_vects_size = 1;
  last_fewer_vects_size = 1;
  for (loop = 1; loop < num_loops; loop++) {
    for (locus = 0; locus < nsystem; locus++) {
      num_vects[locus] *= num_geno[locus][loop];
    }
    max_vects = num_vects[0];
    for (locus = 1; locus < nsystem; locus++) {
      if ( num_vects[locus] > max_vects ) {
        max_vects = num_vects[locus];
      }
    }
    fewer_vects_size = max_vects - 1;
    if (max_vects <= max_vectors_considered ) {
      num_loops_considered++;
      last_fewer_vects_size = fewer_vects_size;
    } else {
      fewer_vects_size = last_fewer_vects_size;
      break;
    }
  }

}  /* set_num_loops_considered */
#endif


static void collapsedown ();

/* Local variables for seg: */
struct LOC_seg {
  thisperson **p, **q, **r, *child, *father, *mother;
  subhap firsthap, secondhap;
  long nfirst, nsecond;
} ;


/*
   Dylan -- Unrolled loops and removed unnecessary computations
         -- Done in mid 1994
*/
 boolean segfun(child, first, second, LINK)
thisperson **child;
long first, second;
struct LOC_seg *LINK;
{
  boolean temp;
  haplotype thishap, thathap;

  if (!sexlink) {
#if LOOPSPEED
    temp = ((*child)->gen->genarray[genenumber[whichsys - 1 ][LINK->secondhap[0]  - 1][LINK->firsthap[0] - 1] - 1] ||
            (*child)->gen->genarray[genenumber[whichsys - 1][LINK->secondhap[0] - 1][LINK->firsthap[1] - 1] - 1] ||
            (*child)->gen->genarray[genenumber[whichsys - 1][LINK->secondhap[1]  - 1][LINK->firsthap[0] - 1] - 1] ||
            (*child)->gen->genarray[genenumber[whichsys - 1][LINK->secondhap[1] - 1][LINK->firsthap[1] - 1] - 1]);
    return temp;
#else

    temp = ((*child)->gen->genarray[genenumber[LINK->secondhap[0] - 1][LINK->firsthap[0] - 1] - 1] ||
            (*child)->gen->genarray[genenumber[LINK->secondhap[0] - 1][LINK->firsthap[1] - 1] - 1] ||
            (*child)->gen->genarray[genenumber[LINK->secondhap[1] - 1][LINK->firsthap[0] - 1] - 1] ||
            (*child)->gen->genarray[genenumber[LINK->secondhap[1] - 1][LINK->firsthap[1] - 1] - 1]);
    return temp;
#endif

  }

  if ((*child)->male) {
    if ((*LINK->p)->male) {
      temp = ((*child)->gen->genarray[LINK->secondhap[0] - 1] ||
              (*child)->gen->genarray[LINK->secondhap[1] - 1]);
    } else {
      temp = ((*child)->gen->genarray[LINK->firsthap[0] - 1] ||
              (*child)->gen->genarray[LINK->firsthap[1] - 1]);
    }
    return temp;
  }

  if ((*LINK->p)->male) {

#if LOOPSPEED
    temp = ((*child)->gen->genarray[genenumber[whichsys - 1][LINK->secondhap[0] - 1][first] - 1] ||
            (*child)->gen->genarray[genenumber[whichsys - 1][LINK->secondhap[1]  - 1][first] - 1]);
  } else {
    temp = ((*child)->gen->genarray[genenumber[whichsys - 1][LINK->firsthap[0] - 1][second] - 1] ||
            (*child)->gen->genarray[genenumber[whichsys - 1][LINK->firsthap[1] - 1][second] - 1]);

#else
    temp = ((*child)->gen->genarray[genenumber[LINK->secondhap[0] - 1][first] - 1] ||
            (*child)->gen->genarray[genenumber[LINK->secondhap[1] - 1][first] - 1]);
  } else {
    temp = ((*child)->gen->genarray[genenumber[LINK->firsthap[0] - 1][second] - 1] ||
            (*child)->gen->genarray[genenumber[LINK->firsthap[1] - 1][second] -  1]);
#endif

  }
  return temp;
}  /*segfun*/



/* Dylan -- modified to cache non-zero entries in q's genarray
         -- rewrote loops to start at 0
         -- eliminated FORMLIM and WITH variables
         -- done in mid 1994
*/
void segup(LINK)
struct LOC_seg *LINK;
{
  long first;
  boolean segval, val;
  int nonzindex, nonzcount;
  boolean compat;
  int j;

  /* array of indexes into genarray of non-zero elements */
  unsigned int nonzgens[maxgeno];

  if ((*LINK->p)->male) {
    LINK->nfirst = mgeno;
    LINK->nsecond = fgeno;
  } else {
    LINK->nfirst = fgeno;
    LINK->nsecond = mgeno;
  }
  /*find nonzero entries in q's genarray and make a list of them
    stored in nonzgens */
  nonzcount = 0;
  { int i;

    for(i = 0; i < LINK->nsecond; i++) {
      if( (*LINK->q)->gen->genarray[i] ) {
        nonzgens[nonzcount++] = i;
      }
    }
  }

  for (first = 0; first < LINK->nfirst; first++) {
    if ( (*LINK->p)->gen->genarray[first] ) {
      segval = false;
      memcpy(LINK->firsthap, seghap[first], sizeof(subhap));
      for (nonzindex = 0; nonzindex < nonzcount; nonzindex++) {
        memcpy(LINK->secondhap, seghap[nonzgens[nonzindex]], sizeof(subhap));
        val = true;
        LINK->child = LINK->father->foff;
        do {
          if (LINK->mother == LINK->child->ma)
            val = segfun(&LINK->child, first, nonzgens[nonzindex], LINK);
          LINK->child = LINK->child->nextpa;
	} while (val && LINK->child != NULL);
        segval = (val || segval);
      }
      (*LINK->p)->gen->genarray[first] = segval;
    }
  }

#if LOOPSPEED
  /*do family by family incompatibility test porvided that LOOPSPEED is
    off or there are no loops*/
  if (0 == num_loops) {
#endif
  if (!incompat_traversal) {
    compat = false;
    for (j = 0; j < LINK->nfirst; j++)
       compat = (compat || ((*LINK->p)->gen->genarray[j]));
    if (!compat && (!detected_incompat[(*LINK->p)->id]) &&
	(!detected_incompat[(*LINK->q)->id])) {
      printf("\n One incompatibility involves the family in which person %d is a parent",(*LINK->p)->id);
      printf("\n The person number refers to the second column in the pedigree file input to UNKNOWN");
      incompat_traversal = true;
      detected_incompat[(*LINK->p)->id] = true;
      detected_incompat[(*LINK->q)->id] = true;
      if (ONE_ERROR_ONLY)
	exit(EXIT_FAILURE);
      if (first_incompat)
	respond();
    }
  }
#if LOOPSPEED
}
#endif

#if !LOOPSPEED
  /*
     cleanup modifies a data structure 'possible' not used in
     LOOPSPEED case
  */
  cleanup(LINK->q);
  LINK->child = LINK->father->foff;
  do {
    if (LINK->child->ma == LINK->mother)
      cleanup(&LINK->child);
    LINK->child = LINK->child->nextpa;
  } while (LINK->child != NULL);
#endif

}  /*segup*/


/* 
   Dylan -- modified to cache non-zero entries in q's genarray
         -- rewrote loops to start at 0
         -- eliminated FORMLIM and WITH variables
         -- Done in mid 1994
*/
 void segdown(LINK)
struct LOC_seg *LINK;
{
  long first, second, here;
  boolean val;
  haplotype thishap, thathap;
  int nonzindex, nonzcount;
  unsigned int nonzgens[maxgeno];
  boolean compat;
  int j;

  for (first = 0; first < fgeno; first++)
    gene[first] = false;

  /*find nonzero entries in q's genarray and make a list of them
    stored in nonzgens */
  nonzcount = 0;
  { int i;

    for(i = 0; i < fgeno; i++) {
      if( (*LINK->q)->gen->genarray[i] ) {
        nonzgens[nonzcount++] = i;
      }
    }
  }
  for (first = 0; first < mgeno; first++) {
    if ((*LINK->p)->gen->genarray[first]) {
      memcpy(LINK->firsthap, seghap[first], sizeof(subhap));
      for (nonzindex = 0; nonzindex < nonzcount; nonzindex++) {
        memcpy(LINK->secondhap, seghap[nonzgens[nonzindex]], sizeof(subhap));

        /* Dylan -- changed from an or expression that was always true */
        val = true; 

        LINK->child = LINK->father->foff;
        do {
          if (LINK->child->ma == LINK->mother) {
            if (!LINK->child->up)
              val = segfun(&LINK->child, first, nonzgens[nonzindex], LINK);
	  }
          LINK->child = LINK->child->nextpa;
	} while (val && LINK->child != NULL);
        if (val) {
          if (!sexlink) {
            for (thishap = 0;
                 thishap <= 1;
                 thishap++) {
              for (thathap = 0;
                   thathap <= 1;
                   thathap++) {
#if LOOPSPEED
                here = genenumber[whichsys - 1][LINK->secondhap[thishap] - 1]
                  [LINK->firsthap[thathap] - 1] - 1;

#else
                here = genenumber[LINK->secondhap[thishap] - 1]
                  [LINK->firsthap[thathap] - 1] - 1;

#endif

                gene[here] = (gene[here] || val);
	      }
	    }
	  } else if ((*LINK->r)->male) {
            for (thathap = 0;
                 thathap <= 1;
                 thathap++) {
              here = LINK->secondhap[thathap] - 1;
              gene[here] = (gene[here] || val);
	    }
	  } else {
            for (thathap = 0;
                 thathap <= 1;
                 thathap++) {

#if LOOPSPEED
              here = genenumber[whichsys - 1][LINK->secondhap[thathap] - 1]
                [first] - 1;

#else

              here = genenumber[LINK->secondhap[thathap] - 1][first] - 1;
#endif
              gene[here] = (gene[here] || val);
	    }
	  }
	}
      }
    }
  }
  for (first = 0; first < fgeno; first++)
    (*LINK->r)->gen->genarray[first] = ((*LINK->r)->gen->genarray[first] && 
                                        gene[first]);

#if LOOPSPEED
  /*do family by family incompatibility test porvided that LOOPSPEED is
    off or there are no loops*/
  if (0 == num_loops) {
#endif
/*compatibility test*/
  if (!incompat_traversal) {
    compat = false;
    for (j = 0; j < fgeno; j++)
       compat = (compat || ((*LINK->r)->gen->genarray[j]));
    if (!compat && (!detected_incompat[(*LINK->p)->id]) &&
	(!detected_incompat[(*LINK->q)->id])) {
      printf("\n One incompatibility involves the family in which person %d is a child",(*LINK->r)->id);
      printf("\n The person number refers to the second column in the pedigree file input to UNKNOWN");
      incompat_traversal = true;
      detected_incompat[(*LINK->p)->id] = true;
      detected_incompat[(*LINK->q)->id] = true;
      if (ONE_ERROR_ONLY)
	exit(EXIT_FAILURE);
      if (first_incompat)
	respond();
    }
  }
#if LOOPSPEED
}
#endif

#if !LOOPSPEED
  /* cleanup modifies a data structure not used in loopspeed case */
  cleanup(LINK->p);
  cleanup(LINK->q);
  LINK->child = LINK->father->foff;
  do {
    if (LINK->child->ma == LINK->mother) {
      if (!LINK->child->up)
        cleanup(&LINK->child);
    }
    LINK->child = LINK->child->nextpa;
  } while (LINK->child != NULL);
#endif

} /*segdown*/



 void seg(p_, q_, r_, peel)
thisperson **p_, **q_, **r_;
direction peel;
{
  struct LOC_seg V;


  V.p = p_;
  V.q = q_;
  V.r = r_;
  if ((*V.p)->male) {
    V.father = *V.p;
    V.mother = *V.q;
  } else {
    V.father = *V.q;
    V.mother = *V.p;
  }
  if (peel == peelup)
    segup(&V);
  else
    segdown(&V);
}  /*seg*/


 void collapseup(p)
thisperson *p;
{
  thisperson *q, *child, *nextchild;
  boolean down;

  depth++;
  if (depth > (DEPTH_MULTIPLE * numind)) {
    printf("The next pedigree appears to have an unbroken loop\n");
    printf("If you do not think so, increase DEPTH_MULTIPLE in unknown.c\n");
    exit(EXIT_FAILURE);
  }
  p->done = true;
  if (p->foff == NULL) {
    depth--;
    return;
  }
  down = false;
  child = p->foff;
  while (child != NULL) {
    down = false;
    if (p->male)
      q = child->ma;
    else
      q = child->pa;
    if (!q->done) {
      collapsedown(q);
      nextchild = child;
      while (nextchild != NULL) {
	if (nextchild->pa == q || nextchild->ma == q) {
	  if (!nextchild->up)
	    collapseup(nextchild);
	  else
	    down = true;
	}
	if (p->male)
	  nextchild = nextchild->nextpa;
	else
	  nextchild = nextchild->nextma;
      }
      if (q->multi)
	collapseup(q);
      if (!down)
	seg(&p, &q, &child, peelup);
      else
	collapsedown(p);
    }
    if (p->male)
      child = child->nextpa;
    else
      child = child->nextma;
  }
  depth--;
}  /*collapseup*/


static void collapsedown(p)
thisperson *p;
{
  depth++;
  if (depth > (DEPTH_MULTIPLE * numind)) {
    printf("The next pedigree appears to have an unbroken loop\n");
    printf("If you do not think so, increase DEPTH_MULTIPLE in unknown.c\n");
    exit(EXIT_FAILURE);
  }
  if (DOWN_CHECK && p->downvisit && (p != proband)) {
    printf("The next pedigree appears to have an unbroken loop\n");
    printf("If you do not think so, change DOWN_CHECK to false in unknown.c\n");
    exit(EXIT_FAILURE);
  }
  else
    p->downvisit = true;
  if (p->pa == NULL) {
    depth--;
    return;
  }
  p->up = true;
  collapseup(p->pa);
  seg(&p->pa, &p->ma, &p, peeldown);
  depth--;
}  /*collapsedown*/


/*
   This procedure mallocs the space needed at this time for the
   data structure 'unknown_poss'.  Written by Dylan in late 1994.
*/
#if LOOPSPEED
static void malloc_unknown_poss(curr_person, init_value)
int curr_person;
boolean init_value;
{
  int i;
  int geno;

  if (unknown_poss == NULL) {
    unknown_poss =
      (geno_for_unknown *) malloc((totperson + 1) * sizeof(geno_for_unknown));
    if (unknown_poss == NULL)
      malloc_err("unknown possibilities table");
    for (i = 1; i <= totperson; i++){
      unknown_poss[i] = NULL;
    }
  }
  if (unknown_poss[curr_person] == NULL) {
    unknown_poss[curr_person] = (geno_for_locus *) 
      malloc(nsystem * sizeof(geno_for_locus));
    if (unknown_poss[curr_person] == NULL)
      malloc_err("unknown possibilities table");
    for (i = 0; i < nsystem; i++){
      unknown_poss[curr_person][i] = NULL;
    }
  }
  if (unknown_poss[curr_person][whichsys - 1] == NULL) {
    unknown_poss[curr_person][whichsys - 1] = (geno_for_loop_vector *)
      malloc(num_loop_vectors[whichsys - 1] * sizeof(geno_for_loop_vector));
    if (unknown_poss[curr_person][whichsys - 1] == NULL)
      malloc_err("unknown possibilities table");
    for (i = 0; i < num_loop_vectors[whichsys - 1]; i++) {
      unknown_poss[curr_person][whichsys - 1][i] = (boolean *)
        malloc (thislocus[whichsys - 1]->fgeno * sizeof(boolean));
      if (unknown_poss[curr_person][whichsys - 1][i] == NULL)
        malloc_err("unknown possibilities table");
      for (geno = 0; geno < thislocus[whichsys - 1]->fgeno; geno++) {
        unknown_poss[curr_person][whichsys - 1][i][geno] = init_value;
      }
    }
  }
} /* malloc_unknown_poss */
#endif

/* 
   This procedure frees the space used by the data structure
   'unknown_poss'.  Written by Dylan in late 1994.
*/
#if LOOPSPEED
static void free_unknown_poss() {
  int curr_person, locus, loop_v;

  if (unknown_poss != NULL) {
    for (curr_person = 1; curr_person <= totperson; curr_person++){
      if (unknown_poss[curr_person] != NULL) {
        for (locus = 0; locus < nsystem; locus++){
          if (unknown_poss[curr_person][locus] != NULL) {
            for (loop_v = 0; loop_v < num_loop_vectors[locus]; loop_v++) {
              if (unknown_poss[curr_person][locus][loop_v] != NULL) {
                free(unknown_poss[curr_person][locus][loop_v]);
	      }
	    } /* for each loop vector */
            free(unknown_poss[curr_person][locus]);
	  }
	}  /* for each locus */
        free(unknown_poss[curr_person]);
      }
    } /* for each person */
    free(unknown_poss);
    unknown_poss = NULL;
  }
} /* free_unknown_poss */
#endif

#if LOOPSPEED
/*
   This procedure traverses the pedigree, using the proband declared
   in pedin.dat, and determines for each valid loopbreaker vector at
   this locus what genotypes the proband may have.  This information
   is stored in the table 'unknown_poss'.

   It also compiles info for compatibility testing.

   This procedure was rewritten by Dylan in late 1994 to consider each
   valid loopbreaker vector individually.  
*/
void likelihood()
{
  int loopmax[maxloop];
  int loopgen[maxloop];
  int geno;
  int FORLIM;
  int i, j;
  int inner_vect;  /* current loopbreaker vector being considered */
  boolean valid_vector;
  thisperson *WITH;
  thisarray *WITH1;
  genotype compat_test_genarray; /* temp. to acculate geno in compat test */
  boolean alldone;
  int nextgeno, temploopmax; /*AAS*/
  int *loopbreaker_nextgeno[maxloop];
  boolean first;

  /* set loopmax */
  for (i = 0; i < num_loops; i++) {
      if (looppers[i][0]->male)
        loopmax[i] = mgeno;
      else
        loopmax[i] = fgeno;
  }

  /* set compat_test_array if testing compatitibility and
     proband is known at this locus */
  if (proband->thisunknown[whichsys - 1] == false) {
    for (i = 0; i < fgeno; i++) {
      compat_test_genarray[i] = false;
    }
  }



  first = false;
  /* init variables to get genotype vector for loopbreakers */
  for (i = num_loops_considered; i < num_loops; i++) {
    getgene(whichsys, looppers[i][0], looppers[i][0]->phen);
    loopbreaker_nextgeno[i] = (int *) malloc(loopmax[i] * sizeof(int));
    if (NULL == loopbreaker_nextgeno[i])
      malloc_err("loopbreaker_nextgeno entry");
    for(geno = 0; geno < loopmax[i]; geno++) 
      loopbreaker_nextgeno[i][geno] = 0;
    for(geno = 0, nextgeno = 0; geno < loopmax[i]; geno++) {
      while((nextgeno < loopmax[i]) && (looppers[i][0]->gen->genarray[nextgeno] == 0))
	nextgeno++;
      loopbreaker_nextgeno[i][geno] = nextgeno;
      if (nextgeno < loopmax[i])
	temploopmax = nextgeno + 1;
      geno = nextgeno;
      nextgeno++;
    }
    loopmax[i] = temploopmax;
    loopgen[i] = loopbreaker_nextgeno[i][0] + 1;
    first = true;
  } 

  /* iterate over all possible genotype vectors for loopbreakers */
  do {

    /* get next genotype vector for loopbreakers */
      i = num_loops_considered;
      if (num_loops > num_loops_considered) { 
	if (!first)
	  do {
	    if (loopgen[i] >= loopmax[i]) 
	      loopgen[i] = loopbreaker_nextgeno[i][0] + 1;
	    else {
	      loopgen[i] = loopbreaker_nextgeno[i][loopgen[i]] + 1;
	      i = num_loops;
	    }
	    i++;
	  } while (i <= num_loops);
	else
	  first = false;
      }

    /* check if this loopbreaker genotype vector is possible */
    valid_vector = true;
    for (i = num_loops_considered; i < num_loops; i++) {
        if (loopgen[i] > loopmax[i]) {
          valid_vector = false;
          break;
	}
    }

    /* if outer loop breaker genotype vector is possible */
    if ( valid_vector == true ) {

      /* iterate over each possible inner loop breaker vector */
      for (inner_vect = 0;
           inner_vect < num_loop_vectors[whichsys - 1];
           inner_vect++){

        /* init each person */
        /* DYLAN -- really only need to get everyone if a loopbreaker can
           have multiple values i.e. it is unknown or a locus is "affection" */
        for (i = 1; i <= totperson; i++) {
          getgene(whichsys, person[i], person[i]->phen);
          WITH = person[i];
          WITH->done = false;
          WITH->up = false;
          WITH->downvisit = false;
	}

        /* set loopbreak genotypes to match loopbreaker genotype vector */
        for (i=0; i < num_loops; i++) {
          for (j=0; j < loopmax[i]; j++) {
            looppers[i][0]->gen->genarray[j] = false;
            looppers[i][1]->gen->genarray[j] = false;
	  }
	}
        for (i = 0; i < num_loops_considered; i++) {
          looppers[i][0]->gen->genarray
            [loop_vectors[whichsys - 1][inner_vect][i]] = true;
          looppers[i][1]->gen->genarray
            [loop_vectors[whichsys - 1][inner_vect][i]] = true;
	}
        for (i = num_loops_considered; i < num_loops; i++) {
          looppers[i][0]->gen->genarray[loopgen[i] - 1] = true;
          looppers[i][1]->gen->genarray[loopgen[i] - 1] = true;
	}

        /* traverse pedigree */
        collapseup(proband);
        collapsedown(proband);
        /* record possible genotypes */
        for (i = 1; i <= totperson; i++) {
 /*The following line changed by A. A. Schaffer to add inference capability
  in the case where person[i] is unknown only at some, but
  not all loci*/
/*        if ( person[i]->thisunknown[whichsys - 1] == true ) {*/
          if ( person[i]->unknown) {
            malloc_unknown_poss(i, false);
/*	    if (!person[i]->thisunknown[whichsys - 1] == true ) 
	      getgene(whichsys, person[i], person[i]->phen);*/
            if (person[i]->male)
              FORLIM = mgeno;
            else
              FORLIM = fgeno;
            for (geno =0; geno < FORLIM; geno++) {
              unknown_poss[i][whichsys - 1][inner_vect][geno] =
                (person[i]->gen->genarray[geno] ||
                 unknown_poss[i][whichsys - 1][inner_vect][geno]);
	    }
	  }
	}

        /* DYLAN -- really should rule out some loopbreaker vectors here */
        /* accumulate info for compat test */
        if ( proband->thisunknown[whichsys - 1] == false ) {
          if (proband->male)
            FORLIM = mgeno;
          else
            FORLIM = fgeno;
          for (i = 0; i < FORLIM; i++) {
            compat_test_genarray[i] = (compat_test_genarray[i] ||
                                       proband->gen->genarray[i]);
	  }
	}

      } /* for */
    } /* if valid vector */

    alldone = true;
    for (i = num_loops_considered; i < num_loops; i++)
      alldone = (alldone && loopgen[i] == loopmax[i]);
  } while (!alldone);


  /* copy out compatibility test array for use in check in iterpeds() */
  if (proband->thisunknown[whichsys - 1] == false) {
    if (proband->male)
      FORLIM = mgeno;
    else
      FORLIM = fgeno;
    for (i = 0; i < FORLIM; i++) {
      proband->gen->genarray[i] = compat_test_genarray[i];
    }
  }

  for (i = num_loops_considered; i < num_loops; i++)
    free(loopbreaker_nextgeno[i]);


}  /*likelihood*/

/*
   This procedure is called when the pedigree has no loops.  Just do the
   traversal.  'unknown_poss' is set later.

   Written by Dylan in late 1994.
*/
static void no_loop_likelihood(proband_index)
int proband_index;
{
  int i;

  proband = person[proband_index];

  /* traverse pedigree */
  collapsedown(proband);
  collapseup(proband);

}  /* no_loop_likelihood*/

#else

static void likelihood(proband)
thisperson **proband;
{
  long i, j;
  thisperson *WITH;
  information *WITH1;
  thisarray *WITH2;
  locusvalues *WITH3;
  long FORLIM, FORLIM1;


  collapsedown(*proband);
  collapseup(*proband);
  if (!(*proband)->thisunknown[whichsys - 1])
    return;
  WITH = *proband;
  WITH1 = WITH->store;
  WITH2 = WITH->gen;
  WITH3 = thislocus[whichsys - 1];
  if (sexlink && WITH->male) {
    FORLIM = WITH3->nallele;
    for (j = 0; j < FORLIM; j++)
      WITH1->possible[whichsys - 1][0][j] = WITH2->genarray[j];
  } else {
    FORLIM = WITH3->nallele;
    for (i = 0; i < FORLIM; i++) {
      FORLIM1 = WITH3->nallele;
      for (j = i; j < FORLIM1; j++) {
	WITH1->possible[whichsys - 1][i]
	  [j] = WITH2->genarray[genenumber[i][j] - 1];
	WITH1->possible[whichsys - 1][j][i] = WITH1->possible[whichsys - 1][i]
	  [j];
      }
    }
  }
  cleanup(proband);
}  /*likelihood*/
#endif

#if LOOPSPEED
/*
   This procedure mallocs enough space in 'loop_vectors' for the row for
   the current locus.  Written by Dylan in late 1994.
*/
static void malloc_loop_vectors()
{
  int i;

  if (loop_vectors == NULL) {
    loop_vectors =
      (vector_for_locus *) malloc(nsystem * sizeof(vector_for_locus));
    if (loop_vectors == NULL)
      malloc_err("loop vector table");
    for (i = 0; i < nsystem; i++){
      loop_vectors[i] = NULL;
    }
  }
  loop_vectors[whichsys - 1] = (vector_for_loop *)
    malloc(num_loop_vectors[whichsys - 1] * sizeof(vector_for_loop));
  if (loop_vectors[whichsys - 1] == NULL)
    malloc_err("loop vector table");
  for (i = 0; i < num_loop_vectors[whichsys - 1]; i++){
    loop_vectors[whichsys - 1][i] = (int *) 
      malloc (num_loops_considered * sizeof(int));
    if (loop_vectors[whichsys - 1][i] == NULL)
      malloc_err("loop vector table");
						    }
} /* malloc_loop_vectors */

/*
   This procedure frees the space used by 'loop_vectors'.  Written by Dylan.
*/
static void free_loop_vectors ()
{
  int locus, vect;

  if (loop_vectors != NULL) {
    for (locus = 0; locus < nsystem; locus++) {
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
#endif

#if LOOPSPEED
/*
   Frees space used by simple linked list.  Written by Dylan in late 1994.
*/
static void free_simple_list(l)
simple_list l;
{
  list_elt *curr, *temp;

  for (curr = l.head; curr; ) {
    temp = curr;
    curr = curr->next;
    free(temp);
  }
} /* free_simple_list */
#endif


/*
   This procedure records all the valid loopbreaker vectors for the
   current locus.  In particular...

   This procedure fills in the 'loop_vectors' row corresponding to the
   current locus, 'whichsys', with the valid loopbreaker vectors.  It
   also records the number of valid loopbreaker vectors in
   'num_loop_vectors'.  Written by Dylan in late 1994.
*/
#if LOOPSPEED
static void get_loopbreaker_vectors() {
  int i, j;
  boolean valid_vector;
  boolean alldone;
  int loopgen[maxloop];
  int loopmax[maxloop];
  int nextgeno, geno, temploopmax; /*AAS*/
  int *loopbreaker_nextgeno[maxloop];
  simple_list temp_list;
  int loop_vect;

  /* if no loops, little to do */
  if ( num_loops_considered == 0 ) {
    num_loop_vectors[whichsys - 1] = 1; /* the null loop vector */
    return;
  }
  /* init */
  num_loop_vectors[whichsys - 1] = 0;
  temp_list.head = NULL;
  temp_list.curr = NULL;

  /* init variables to get genotype vector for loopbreakers */
  for (i = 0; i < num_loops_considered; i++) {
    loopgen[i] = 1;
    loopmax[i] = 1;
    if (looppers[i][0] != NULL) {
      if (looppers[i][0]->male)
        loopmax[i] = mgeno;
      else
        loopmax[i] = fgeno;
      getgene(whichsys, looppers[i][0], looppers[i][0]->phen);
    }
    loopbreaker_nextgeno[i] = (int *) malloc(loopmax[i] * sizeof(int));
    if (NULL == loopbreaker_nextgeno[i])
      malloc_err("loopbreaker_nextgeno entry");
    for(geno = 0; geno < loopmax[i]; geno++) 
      loopbreaker_nextgeno[i][geno] = 0;
    for(geno = 0, nextgeno = 0; geno < loopmax[i]; geno++) {
      while((nextgeno < loopmax[i]) && (looppers[i][0]->gen->genarray[nextgeno] == 0))
        nextgeno++;
      loopbreaker_nextgeno[i][geno] = nextgeno;
      if (nextgeno < loopmax[i])
        temploopmax = nextgeno + 1;
      geno = nextgeno;
      nextgeno++;
    }
    loopmax[i] = temploopmax;
    loopgen[i] = loopbreaker_nextgeno[i][0] + 1;
  } 


  do {

    /* accumulate possible values (to make vectors) into a list */
      if (temp_list.curr == NULL) {
        temp_list.curr = (list_elt *) malloc (sizeof(list_elt));
        if (temp_list.curr == NULL)
          malloc_err("loop breaker vector");
      } else {
        temp_list.curr->next = (list_elt *) malloc (sizeof(list_elt));
        if (temp_list.curr == NULL)
          malloc_err("loop breaker vector");
        temp_list.curr= temp_list.curr->next;
      }
      temp_list.curr->value = loopgen[0] - 1;
      temp_list.curr->next = NULL;
      if (temp_list.head == NULL) {
        temp_list.head = temp_list.curr;
      }
      for (i = 1; i < num_loops_considered; i++) {
        temp_list.curr->next = (list_elt *) malloc (sizeof(list_elt));
        if (temp_list.curr->next == NULL)
          malloc_err("loop breaker vector");
        temp_list.curr = temp_list.curr->next;
        temp_list.curr->next = NULL;
        temp_list.curr->value = loopgen[i] - 1;
      }
      num_loop_vectors[whichsys - 1]++;
    alldone = true;
    for (i = 0; i < num_loops_considered; i++)
      alldone = (alldone && loopgen[i] == loopmax[i]);

    if (!alldone) {
      /* get next genotype vector for loopbreakers */
      i = 0;
      do {
	if (loopgen[i] >= loopmax[i]) 
	  loopgen[i] = loopbreaker_nextgeno[i][0] + 1;
	else {
	  loopgen[i] = loopbreaker_nextgeno[i][loopgen[i]] + 1;
	  i = num_loops_considered;
	}
	i++;
      } while (i <= num_loops_considered);
    }
  } while (!alldone);

  for (i = 0; i < num_loops_considered; i++)
    free(loopbreaker_nextgeno[i]);

  /* move vectors from list to table now that size is known */
  malloc_loop_vectors();
  loop_vect = 0;
  temp_list.curr = temp_list.head;
  for (j = 0; j < num_loop_vectors[whichsys - 1]; j++) {
    for (i = 0; i < num_loops_considered; i++) {
      loop_vectors[whichsys-1][loop_vect][i] = temp_list.curr->value;
      temp_list.curr = temp_list.curr->next;
    }
    loop_vect++;
  }
  free_simple_list(temp_list);

} /* get_loopbreaker_vectors */
#endif

#if LOOPSPEED
/*
   This procedure calls likelihood for each person who is unknown at
   the current locus.  It also calls likelihood on the first person
   in the pedigree to check that the pedigree is compatitible.

   This procedure was rewritten by Dylan in late 1994.
*/

static void iterpeds()
{
  int i, j;
  int geno;
  int FORLIM;
  int loop_vect, loop_vect2;
  boolean compat;

  compat = false;

  /* find all valid loopbreaker vectors at this locus */
  get_loopbreaker_vectors();

  likelihood();

  if ( proband->thisunknown[whichsys - 1] == true ) {
    for (j = 1; j <= totperson; j++) {
      if ( proband == person[j] ) {
        i = j;
        break;
      }
    }
    if (proband->male)
      FORLIM = thislocus[whichsys - 1]->mgeno;
    else
      FORLIM = thislocus[whichsys - 1]->fgeno;
    for (loop_vect = 0;
         loop_vect < num_loop_vectors[whichsys - 1];
         loop_vect++)
      for (j = 0; j < FORLIM; j++)
        compat = (compat || unknown_poss[i][whichsys-1][loop_vect][j]);
  } else {
    if (proband->male)
      FORLIM = thislocus[whichsys - 1]->mgeno;
    else
      FORLIM = thislocus[whichsys - 1]->fgeno;
    for (j = 0; j < FORLIM; j++)
      compat = (compat || proband->gen->genarray[j]);
  }
  if (!compat) {
    printf("ERROR: Incompatibility detected in this family for locus %12ld\n",
           whichsys);
    respond();
}
  /* set unknown_poss entry to known genotype for this locus */
  /*A. A. Schaffer modified the code to do this in likelihood
    instead*/
/*  if ( compat ) {
    for (i = 1; i <= totperson; i++) {
      if (person[i]->unknown && !person[i]->thisunknown[whichsys - 1] ) {
        malloc_unknown_poss(i, false);
        getgene(whichsys, person[i], person[i]->phen);
  Note possible problem for SEXLINK here, but commented out
        for (geno = 0; geno < thislocus[whichsys - 1]->fgeno; geno++) {
          for (loop_vect2 = 0;
               loop_vect2 < num_loop_vectors[whichsys - 1];
               loop_vect2++)
          unknown_poss[i][whichsys - 1][loop_vect2][geno] =
            person[i]->gen->genarray[geno];
	}
      }
    }
}*/

}  /*iterpeds*/

/*
   This procedure is used when there are no loops in the pedigree.  It
   sets 'unknown_poss' to hold the same values as 'possible' would have
   with LOOPSPEED undefined, but only updates 'unknown_poss' at the end.

   Note that we only getgene() for each person once, allowing their
   genotypes to be narrowed on each traversal.

   Written by Dylan in late 1994.
   Modified by A. A. Schaffer in summer 1995 to allow for incompatibility
   location within nuclear families
*/
static void no_loop_iterpeds()
{
  int i, j;
  int FORLIM;
  int geno;
  int loop_vect;
  boolean compattest, compat;

  first_incompat = false;
  compattest = false;  /* have not done a comptability test yet */
  compat = false;

  /* trivially init */
  get_loopbreaker_vectors();

  for (i= 1; i <= totperson; i++) {
    getgene(whichsys, person[i], person[i]->phen);
  }

  for (i = 1; i <= totperson; i++) {
    if (!compattest || person[i]->thisunknown[whichsys - 1]) {
      incompat_traversal = false;
      /* init for segup() and segdown() */
      for (j = 1; j <= totperson; j++) {
        person[j]->done = false;
        person[j]->up = false;
        person[j]->downvisit = false;
      }

      /* traverse for each person */
      no_loop_likelihood(i);

      /* test compatibilty */
      if (!compattest) {
        compattest = true;
      if (proband->male)
        FORLIM = thislocus[whichsys - 1]->mgeno;
      else
        FORLIM = thislocus[whichsys - 1]->fgeno;
        for (j = 0; j < FORLIM; j++) {
          compat = (compat || person[i]->gen->genarray[j]);
	}
        if (!compat) {
          printf("ERROR: Incompatibility detected in this family for locus %12ld\n",
                 whichsys);
          respond();
	  first_incompat = true;
        }
      }
    }
  }

  /* set unknown_poss from the genarrays */
  for (i = 1; i <= totperson; i++) {
    if ( person[i]->unknown ) {
      malloc_unknown_poss(i, false);
      if (person[i]->male)
        FORLIM = thislocus[whichsys - 1]->mgeno;
      else
        FORLIM = thislocus[whichsys - 1]->fgeno;
      for (j = 0; j < FORLIM; j++) {
        unknown_poss[i][whichsys - 1][0][j] = person[i]->gen->genarray[j];
      }
    }
  }

}  /* no_loop_iterpeds*/

#else

/*old iterpeds*/
static void iterpeds()
{
  long i, j;
  boolean compattest, compat;
  long FORLIM1;

  /* This means that this part of unknown is not active for pedigrees with loops! */
  if (loop1 != NULL || loop2 != NULL)
    return;
  for (i = 1; i <= totperson; i++)
    getgene(whichsys, person[i], person[i]->phen);
  first_incompat = false;
  compattest = false;
  compat = false;
  for (i = 1; i <= totperson; i++) {
    if (!compattest || person[i]->thisunknown[whichsys - 1]) {
      incompat_traversal=false;
      for (j = 1; j <= totperson; j++) {
	person[j]->done = false;
	person[j]->up = false;
        person[j]->downvisit = false;
      }
      likelihood(&person[i]);
      if (!compattest) { /*Only do the overall compatibility test once per
                           locus-pedigree pair*/
	compattest = true;
        if (person[i]->male)
          FORLIM1 = mgeno;
        else
	  FORLIM1 = fgeno;
	for (j = 0; j < FORLIM1; j++)
	  compat = (compat || person[i]->gen->genarray[j]);
	if (!compat) {
	  printf("ERROR: Incompatibility detected in this family for locus %12ld\n",
		 whichsys);
	  respond();
	  first_incompat = true;
	}
      }
    }
  }
  for (i = 1; i <= totperson; i++) {
    if (person[i]->unknown)
      cleanup(&person[i]);
  }
}  /*iterpeds*/
#endif

static void reinit()
{
  long i, j;

  for (i = 1; i <= totperson; i++) {
    for (j = 0; j < nsystem; j++)
      Free(person[i]->phen[j]);
  }

#if !LOOPSPEED
  for (i = 1; i <= totperson; i++) {
    if (person[i]->store != NULL)
      Free(person[i]->store);
  }
#endif

  for (i = 1; i <= totperson; i++) {
    Free(person[i]->gen);
    Free(person[i]);
    person[i] = NULL;
  }

#if LOOPSPEED
  free_unknown_poss();
  free_loop_vectors();
#endif

}  /*reinit*/


 boolean testhets()
{
  /*Change: Function added by Joe Terwilliger 7/8/93*/
  double a_, b_, c;
  long prog, d, numl, lc, nall, sexl, nqv, i, j, sexd, int_;
  boolean tmp;
  char fff;

  tmp = true;
  fscanf(datafile, "%ld%lg%ld%ld%*[^\n]", &numl, &b_, &sexl, &prog);
  getc(datafile);
  fscanf(datafile, "%lg%lg%lg%ld%*[^\n]", &a_, &b_, &c, &d);
  getc(datafile);
  if (d == 1) {
    tmp = false;
    goto _L10;
  }
  if (prog != 1 && prog != 3)
    goto _L10;
  fscanf(datafile, "%lg%*[^\n]", &a_);
  getc(datafile);
  for (j = 1; j <= numl; j++) {
    fscanf(datafile, "%ld", &lc);
    switch (lc) {

    case 0:
      fscanf(datafile, "%ld%*[^\n]", &nall);
      getc(datafile);
      fscanf(datafile, "%*[^\n]");
      getc(datafile);
      fscanf(datafile, "%ld%*[^\n]", &nqv);
      getc(datafile);
      for (i = 1; i <= nqv; i++) {
	fscanf(datafile, "%lg%*[^\n]", &a_);
	getc(datafile);
      }
      fscanf(datafile, "%lg%*[^\n]", &a_);
      getc(datafile);
      fscanf(datafile, "%lg%*[^\n]", &a_);
      getc(datafile);
      break;

    case 1:
      fscanf(datafile, "%ld%*[^\n]", &nall);
      getc(datafile);
      fscanf(datafile, "%lg%*[^\n]", &a_);
      getc(datafile);
      fscanf(datafile, "%ld%*[^\n]", &nall);
      getc(datafile);
      if (sexl == 0) {
	for (i = 1; i <= nall; i++) {
	  fscanf(datafile, "%lg%*[^\n]", &a_);
	  getc(datafile);
	}
      } else {
	for (i = 1; i <= nall + nall; i++) {
	  fscanf(datafile, "%lg%*[^\n]", &a_);
	  getc(datafile);
	}
      }
      break;

    case 2:
      fscanf(datafile, "%ld%*[^\n]", &nall);
      getc(datafile);
      for (i = 1; i <= nall + 2; i++) {
	fscanf(datafile, "%lg%*[^\n]", &a_);
	getc(datafile);
      }
      break;

    case 3:
      fscanf(datafile, "%ld%*[^\n]", &nall);
      getc(datafile);
      fscanf(datafile, "%lg%*[^\n]", &a_);
      getc(datafile);
      break;
    }
  }
  fscanf(datafile, "%ld%ld%*[^\n]", &sexd, &int_);
  getc(datafile);
  if (sexd != 0) {
    fscanf(datafile, "%lg%*[^\n]", &a_);
    getc(datafile);
  }
  fscanf(datafile, "%lg%*[^\n]", &a_);
  getc(datafile);
  numl--;
  if (numl == 2 && int_ == 1)
    numl = 3;
  if (sexd == 1)
    numl++;
  if (sexd == 2)
    numl += numl;
  fscanf(datafile, "%lg%*[^\n]", &a_);
  getc(datafile);
  for (i = 1; i <= numl; i++)
    fscanf(datafile, "%ld", &d);
  fff = ' ';
  while (!P_eoln(datafile) && fff != '1') {
    fff = getc(datafile);
    if (fff == '\n')
      fff = ' ';
  }
  if (fff == '1')
    tmp = false;
_L10:
  return tmp;
}  /*testhets*/


static void initunknown_old()
{
  printf("Program UNKNOWN version %s\n", aVersion);
  printf("The following maximum values are in effect:\n");
  printf("%8ld loci\n", (long)maxlocus);
  printf("%8ld single locus genotypes\n", maxgeno);
  printf("%8ld alleles at a single locus\n", (long)maxall);
  printf("%8ld individuals in one pedigree\n", (long)maxind);
  printf("%8ld marriage(s) for one male\n", (long)maxmarriage);
  printf("%8ld quantitative factor(s) at a single locus\n", (long)maxtrait);
  printf("%8ld liability classes\n", (long)maxliab);
  printf("%8ld binary codes at a single locus\n", maxfact);
  printf("%8ld maximum number of loops\n", maxloop);
  one = 1.00001;   /*change*/

  printf("Opening DATAFILE.DAT\n");
  if (datafile != NULL)
    datafile = freopen("datafile.dat", "r", datafile);
  else
    datafile = fopen("datafile.dat", "r");
  if (datafile == NULL)
    exit(FileNotFound);
  if (P_eof(datafile)) {
    printf(
      "ERROR: File empty or nonexistent. Press <Enter> to continue or <Ctrl-C> to abort\n");
    scanf("%*[^\n]");
    getchar();
  }
  makehomozygous = testhets();   /*Change - 2 lines added*/
  rewind(datafile);


  if (ipedfile != NULL)
    ipedfile = freopen("ipedfile.dat", "w", ipedfile);
  else
    ipedfile = fopen("ipedfile.dat", "w");
  if (ipedfile == NULL)
    exit(FileNotFound);

#if LOOPSPEED
  if (loopfile != NULL)
    loopfile = freopen(LOOPFILE_NAME, "w", loopfile);
  else
    loopfile = fopen(LOOPFILE_NAME, "w");
  if (loopfile == NULL)
    exit(FileNotFound);
#endif

  if (speedfile != NULL)
    speedfile = freopen("speedfile.dat", "w", speedfile);
  else
    speedfile = fopen("speedfile.dat", "w");
  if (speedfile == NULL)
    exit(FileNotFound);

  /*changed from speedfileName*/
}  /*initunknown_old*/


static void initunknown_new()
{
  printf("Reopening DATAFILE.DAT\n");
  if (datafile != NULL)
    datafile = freopen("datafile.dat", "r", datafile);
  else
    datafile = fopen("datafile.dat", "r");
  if (datafile == NULL)
    exit(FileNotFound);
  if (P_eof(datafile)) {
    printf(
      "ERROR: File empty or nonexistent. Press <Enter> to continue or <Ctrl-C> to abort\n");
    scanf("%*[^\n]");
    getchar();
  }
  makehomozygous = testhets();   /*Change - 2 lines added*/
  rewind(datafile);


  if (ipedfile != NULL)
    ipedfile = freopen("ipedfile.dat", "w", ipedfile);
  else
    ipedfile = fopen("ipedfile.dat", "w");
  if (ipedfile == NULL)
    exit(FileNotFound);

#if LOOPSPEED
  if (loopfile != NULL)
    loopfile = freopen(LOOPFILE_NAME, "w", loopfile);
  else
    loopfile = fopen(LOOPFILE_NAME, "w");
  if (loopfile == NULL)
    exit(FileNotFound);
#endif

  if (speedfile != NULL)
    speedfile = freopen("newspeedfile.dat", "w", speedfile);
  else
    speedfile = fopen("newspeedfile.dat", "w");
  if (speedfile == NULL)
    exit(FileNotFound);

  /*changed from speedfileName*/
}  /*initunknown_new*/


main_old(argc, argv) /* generate "speedfile.dat" in old format - Tony */
int argc;
char *argv[];
{
  locusvalues *WITH;
  int ind, ped;
  int num_alleles;

  ipedfile = NULL;
  pedfile = NULL;
  datafile = NULL;
  unlink("newspeedfile.dat");
  speedfile = NULL;
#if LOOPSPEED
  unlink(LOOPFILE_NAME);
  loopfile = NULL;
#endif
  initunknown_old();
  readloci();
  printf("Opening PEDFILE.DAT\n");
  if (pedfile != NULL)
    pedfile = freopen("pedfile.dat", "r", pedfile);
  else
    pedfile = fopen("pedfile.dat", "r");
  if (pedfile == NULL)
    exit(FileNotFound);
  if (P_eof(pedfile)) {
    printf(
      "ERROR: File empty or nonexistent. Press <Enter> to continue or <Ctrl-C> to abort\n");
    scanf("%*[^\n]");
    getchar();
  }
  /*CLOSE(datafile);*/
  nsequence = 0;
  if (!P_eof(pedfile))
    fscanf(pedfile, "%ld", &newped);

  ped = 0;
  while (!P_eof(pedfile)) {
    numind = 0;
    depth = 0;
    currentped = ped;
    ped++;
    readped_old();
    getunknown_old();
#if LOOPSPEED
    if ( num_loops > 0 ) {
      set_num_loops_considered();
    } else {
      num_loops_considered = 0;
    }
#endif
    for (whichsys = 1; whichsys <= nsystem; whichsys++) {
      if (mutsys != whichsys) {
        for(ind = 0; ind <=maxind; ind++)
	  detected_incompat[ind] = false;
	WITH = thislocus[whichsys - 1];
	fgeno = WITH->nallele * (WITH->nallele + 1) / 2;
	if (sexlink)
	  mgeno = WITH->nallele;
	else
	  mgeno = fgeno;
	getlocation(thislocus[whichsys - 1]);

#if LOOPSPEED
        WITH->mgeno = mgeno;
        WITH->fgeno = fgeno;
        if ( num_loops == 0 ) {
          no_loop_iterpeds();
	} else {
          iterpeds();
	}
#else
        iterpeds();
#endif

      }
      else {
#if LOOPSPEED
        fprintf(stderr, "\nYou cannot have LOOPSPEED set to 1 and use the mutation model.");
        fprintf(stderr,"\nIf you did not intend to use the mutation model,");
        fprintf(stderr,"\nchange the second line of the locus file to be");
        fprintf(stderr,"\n0 0.0 0.0 0");
        fprintf(stderr,"\nand redo your run.");
        fprintf(stderr,"\nIf you really intended to use the mutation model,");
        fprintf(stderr,"\nthen recompile UNKNOWN and the main programs");
        fprintf(stderr,"\nwith LOOPSPEED set to 0, and then rerun.");
        exit(EXIT_FAILURE);
#else
       ;
#endif
      }
    }
    infer_old();
    writeped_old();
    writespeed_old();
#if LOOPSPEED
    write_loopfile_old(ped);
#endif
    reinit();
  }
  /*CLOSE(pedfile);
  CLOSE(speedfile);
  CLOSE(ipedfile);*/
  if (datafile != NULL)
    FCLOSE(datafile);
  if (pedfile != NULL)
    FCLOSE(pedfile);
  if (ipedfile != NULL)
    FCLOSE(ipedfile);
#if LOOPSPEED
  if (loopfile != NULL)
    FCLOSE(loopfile);
#endif
  if (speedfile != NULL)
    FCLOSE(speedfile);
} /* main_old */

main_new(argc, argv) /* generate "newspeedfile.dat" - Tony */
int argc;
char *argv[];
{
  locusvalues *WITH;
  int ind, ped;
  int num_alleles;

  ipedfile = NULL;
  pedfile = NULL;
  datafile = NULL;
  speedfile = NULL;
#if LOOPSPEED
  unlink(LOOPFILE_NAME);
  loopfile = NULL;
#endif
  initunknown_new();
  readloci();
  init_ped_loc_all_new();
  printf("Reopening PEDFILE.DAT\n");
  if (pedfile != NULL)
    pedfile = freopen("pedfile.dat", "r", pedfile);
  else
    pedfile = fopen("pedfile.dat", "r");
  if (pedfile == NULL)
    exit(FileNotFound);
  if (P_eof(pedfile)) {
    printf(
      "ERROR: File empty or nonexistent. Press <Enter> to continue or <Ctrl-C> to abort\n");
    scanf("%*[^\n]");
    getchar();
  }
  /*CLOSE(datafile);*/
  nsequence = 0;
  if (!P_eof(pedfile))
    fscanf(pedfile, "%ld", &newped);

  ped = 0;
  while (!P_eof(pedfile)) {
    numind = 0;
    depth = 0;
    currentped = ped;
    ped++;
    readped_new();
    adjust_alleles();
    allele_adjust_persons();
    getunknown_new();
#if LOOPSPEED
    if ( num_loops > 0 ) {
      set_num_loops_considered();
    } else {
      num_loops_considered = 0;
    }
#endif
    for (whichsys = 1; whichsys <= nsystem; whichsys++) {
      if (mutsys != whichsys) {
        for(ind = 0; ind <=maxind; ind++)
	  detected_incompat[ind] = false;
	WITH = thislocus[whichsys - 1];
        if ((binary_ == WITH->which) &&
            (allformat == WITH->UU.U2.format))
          num_alleles = ped_new_allele_count[currentped][whichsys -1];
        else
          num_alleles = WITH->nallele;
	fgeno = num_alleles * (num_alleles + 1) / 2;
        WITH->nallele = num_alleles;
	if (sexlink)
	  mgeno = num_alleles;
	else
	  mgeno = fgeno;
	getlocation(thislocus[whichsys - 1]);

#if LOOPSPEED
        WITH->mgeno = mgeno;
        WITH->fgeno = fgeno;
        if ( num_loops == 0 ) {
          no_loop_iterpeds();
	} else {
          iterpeds();
	}
#else
        iterpeds();
#endif

      }
      else {
#if LOOPSPEED
        fprintf(stderr, "\nYou cannot have LOOPSPEED set to 1 and use the mutation model.");
        fprintf(stderr,"\nIf you did not intend to use the mutation model,");
        fprintf(stderr,"\nchange the second line of the locus file to be");
        fprintf(stderr,"\n0 0.0 0.0 0");
        fprintf(stderr,"\nand redo your run.");
        fprintf(stderr,"\nIf you really intended to use the mutation model,");
        fprintf(stderr,"\nthen recompile UNKNOWN and the main programs");
        fprintf(stderr,"\nwith LOOPSPEED set to 0, and then rerun.");
        exit(EXIT_FAILURE);
#else
       ;
#endif
      }
    }
    infer_new();
    writeped_new();
    writespeed_new();
#if LOOPSPEED
    write_loopfile_new(ped);
#endif
    reinit();
  }
  /*CLOSE(pedfile);
  CLOSE(speedfile);
  CLOSE(ipedfile);*/
  if (datafile != NULL)
    FCLOSE(datafile);
  if (pedfile != NULL)
    FCLOSE(pedfile);
  if (ipedfile != NULL)
    FCLOSE(ipedfile);
#if LOOPSPEED
  if (loopfile != NULL)
    FCLOSE(loopfile);
#endif
  if (speedfile != NULL)
    FCLOSE(speedfile);
} /* main_new */

main(argc, argv)
int argc;
char *argv[];
{
   main_old(argc,argv);
  if (ALLELE_SPEED == 1)
    main_new(argc,argv);
  exit(EXIT_SUCCESS);
}
