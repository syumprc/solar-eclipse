#ifndef  _COMMONDEFS_H

#define  _COMMONDEFS_H  1

/* Output from p2c, the Pascal-to-C translator */
/* From input file "ilink.p" */
/* This file contains definitions for a modified version of the ILINK program*/
/* The modifications are described in the papers: */
/* R. W. Cottingham, Jr., R. M. Idury, and A. A. Schaffer */
/* Faster Sequential Genetic Linkage Computations */
/* American Journal of Human Genetics, 53(1993), pp. 252-263*/
/* and A. A. Schaffer, S. K. Gupta, K. Shriram, and R. W. Cottingham, Jr., */
/* Avoiding Recomputation in Linkage Analysis, */
/* Human Heredity 44(1994), pp. 225-237. */


/* SOLAR modified minimally 2012 to compile */

#include <string.h>

/* VMS: MAY NEED TO CHANGE --
   possibly uncomment the following 2 lines */
/* #include unixio */
/* #include file */

#include <stdio.h>

/* Shriram: begin */
#include <errno.h>
#include <math.h>

/* VMS: MAY NEED TO CHANGE --
   possibly remove "sys/" from the next 3 lines */
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
/* Shriram: end */

/* cgh */
#if !defined(vms)
/* #include <malloc.h> SOLAR */
#endif

/* VMS: MAY NEED TO CHANGE --
   comment out any of the next 3 lines for files your system can't find */
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
/* end cgh */

/* Define LOOPSPEED to be 1 to run with loop speedups added by Dylan and
   by A. A. Schaffer; cannot work with LODSCORE because lcp-produced
   scripts for LODSCORE do not use unknown*/

#if  defined(LODSCORE)
#define LOOPSPEED 0
#endif

#if !defined(LOOPSPEED)
#define LOOPSPEED 1
#endif /*LOOPSPEED*/

#if !defined(LOOP_BREAKERS)
#if LOOPSPEED
#define LOOP_BREAKERS 1
#else 
#define LOOP_BREAKERS 0
#endif /*LOOP_BREAKERS*/
#endif /*defined(LOOP_BREAKERS)*/


/*ALLELE_SPEED is used to set up renumbering of alleles when not all
  alleles are present in a pedigree*/
#if !defined(ALLELE_SPEED)
#define ALLELE_SPEED 1
#endif
#define ALLELE_SPEED_CONSTANT 61192012


#define true 1
#define false 0
#ifndef NULL
#define NULL 0
#endif

/*From p2c.h*/
#ifndef FileNotFound
#define FileNotFound 10 
#endif

/* Kimmo Kallio */
#if defined(vms)
#define unlink delete
#endif  /* defined(vms) */
/* end Kimmo */


#ifndef EXIT_FAILURE
#ifdef vms
#define EXIT_FAILURE (02000000000L)
#else
#ifndef SOLARIS
#define EXIT_FAILURE (-1)
#endif   /* SOLARIS */
#endif   
/*from p2c.h*/
#endif

#ifndef EXIT_SUCCESS
#ifdef vms
#define EXIT_SUCCESS 1
#else
#define EXIT_SUCCESS 0
#endif   
/*from p2c.h*/
#endif

#ifndef Void
#define Void void
#endif

#ifndef Local
#define Local static
#endif

#ifndef Free
#define Free free
#endif

#ifndef Malloc
#define Malloc malloc
#endif


#ifndef LOOPFILE_NAME  /*Dylan added*/
#define LOOPFILE_NAME "loopfile.dat"
#endif



#define fVersion          5.1     /*PRESENT VERSION OF LINKAGE*/
#define fastversion     "4.1P"   /*PRESENT VERSION OF FASTLINK*/
/*A. A. Schaffer*/
#define DIAGNOSTIC       true /*should we test settings of constants*/

    /*MAXIMUM NUMBER OF RECOMBINATION PROBABILITIES*/
/*THE FOLLOWING SHOULD BE LARGER THAN MININT*/

#ifndef maxcensor
#define maxcensor       50000L   /*MAXIMUM FOR CENSORING ARRAY*/
#endif

#ifndef maxsys
#define maxsys          70   /*MAXIMUM NUMBER OF SYSTEMS*/
#endif

#ifndef maxlocus
#define maxlocus        8   /*MAXIMUM NUMBER OF LOCI */
#endif

#define maxrec          maxlocus   /*MAXIMUM POSSIBLE NUMBER OF RECOMB/MEI*/

#ifndef maxseg
#define maxseg          128   /*    (maxlocus-1)        */
/* = 2              I.E. 2 TO THE POWER maxlocus-1 */
#endif

#ifndef BIGMAXSEG
#define BIGMAXSEG       16384
#endif

#ifndef BIGMAXLOCUS
#define BIGMAXLOCUS     15
#endif

#ifndef maxall
#define maxall          100   /*MAX NUMBER OF ALLELES AT A SINGLE LOCUS*/
#endif

#ifndef maxind
#define maxind          10000   /*MAXIMUM NUMBER OF INDIVIDUALS*/
#endif 

#ifndef maxped
#define maxped          5000   /*MAXIMUM NUMBER OF PEDIGREES*/
#endif

#ifndef maxchild
#define maxchild        500   /*MAXIMUM NUMBER OF FULLSIBS IN A SIBSHIP*/
#endif

#ifndef maxloop
#define maxloop         8  /*MAXIMUM NUMBER OF LOOPS PER PEDIGREE*/
#endif

#define minfreq         0.0

#define affall          2
/*DISEASE ALLELE FOR QUANTITATIVE TRAITS
                               OR AFFECTION STATUS*/
/* QUANTITATIVE TRAIT */
#define maxtrait        3
    /*MAXIMUM NUMBER OF QUANTITATIVE FACTORS AT A SINGLE LOCUS*/

#define missval         0.0   /*MISSING VALUES FOR QUANTITATIVE TRAITS */
/* AFFECTION STATUS */

#define missaff         0   /*MISSING VALUE FOR AFFECTION STATUS */
#define affval          2   /*CODE FOR AFFECTED INDIVIDUAL*/

#ifndef maxliab
#define maxliab        60   /*MAXIMUM NUMBER OF LIABILITY CLASSES */
#endif

/*Do not change the following two declarations*/
#define binformat       2    /*DO NOT CHANGE THIS !*/
#define allformat       3    /*DO NOT CHANGE THIS, AAS*/

/* BINARY (FACTOR UNION) SYSTEM */
#ifndef maxfact
#define maxfact         maxall
#endif
    /*MAXIMUM NUMBER OF BINARY CODES AT A SINGLE LOCUS*/
/* OTHERS */


#define scalemult       2.0   /*SCALE WEIGHT FOR EACH LOCUS INCLUDED*/

    /*TRUE IF ESTIMATING PARAMETERS OTHER THAN REC*/
#ifndef fitmodel
#define fitmodel       true
#endif

/*Note: It is probably a bad idea to change dostream to false if using 
  LINKMAP*/
#ifndef dostream
#define dostream        true   /*STREAM FILE OUTPUT*/
#endif

#define zerolike        (-1.0e20)
    /*FOR INCONSISTENT DATA OR RECOMBINATION */

/* Shriram: begin                       Shriram's additions commented out
#ifdef M_LN10
#define log10_          M_LN10
#else
*/
#define log10_          2.30259
/*
#endif
   Shriram: end */

#define minint          (-32767)   /*MINIMUM ALLOWED INTEGER*/
/*GRADIENT APPROXIMATIONS*/

#ifndef approximate
#define approximate     false
#endif

#ifndef PARALLEL
#define PARALLEL  0
#endif

#if PARALLEL
int maxworkingset;
#endif /*PARALLEL*/


#if PARALLEL
#ifndef maxn
#define maxn 8
#endif /*maxn*/

#ifndef MAXNUMFAM
#define MAXNUMFAM 1000
#endif

#endif /*PARALLEL*/

/* dwix */
#define DEFAULT_STRING_LENGTH 200
/* end dwix */

#if PARALLEL

void Tmk_errexit();

#if defined(IS_PARMACS) || defined(IS_P4)
#define IS_SHMEM (IS_PARMACS || IS_P4)

char* Tmk_malloc(unsigned);

#if IS_P4

#include "p4.h"
#define G_MALLOC(size)    p4_shmalloc(size)
#define CREATE(proc)      p4_create(proc);
#define WAIT_FOR_END(n)   p4_wait_for_end();
#define Tmk_free(p)       p4_shfree(p)

#endif  /* IS_P4 */

#else  /* defined(IS_PARMACS) || defined(IS_P4) */
#define IS_SHMEM 0
#endif  /* defined(IS_PARMACS) || defined(IS_P4) */

#ifndef BARRIER_OUTPUT
#define BARRIER_OUTPUT  0
#endif
#ifndef LOADBALANCE_OUTPUT
#define LOADBALANCE_OUTPUT 0
#endif
#ifndef FUNTIMES
#define FUNTIMES   1
#endif
#ifndef PARALLEL_GCENTRAL
#define PARALLEL_GCENTRAL 0
#endif
#ifndef DO_LOADBALANCE
#define DO_LOADBALANCE    0
#endif
#ifndef DO_PARALLELTHETAS
#define DO_PARALLELTHETAS  1
#endif
#define  maxprocs  8 /*maximum number of processors*/
#if IS_SHMEM
unsigned Tmk_proc_id;
unsigned Tmk_nprocs;
#endif /*IS_SHMEM*/

/*flags to determine which kind of nuclear family visit is coming next*/
#define ONEDOWN  0
#define MANYDOWN 1
#define ONEUP    2
#define MANYUP   3

#ifndef PRECOMPUTE
#define PRECOMPUTE 0
#endif  /*PRECOMPUTE*/
#endif

#if PARALLEL 
extern char *optarg;
#if !IS_SHMEM
#include "Tmk.h"
#endif /* !IS_SHMEM */
#endif /* PARALLEL */

/* cgh -- as per suggestion from alc: we use this type as paramater to
   functions which take boolean argument to guarantee automagic
   promotion to int for K&R function declarations with ANSI
   prototypes. */
typedef int boolInt;

/*next two typedefs come from p2c*/
typedef char boolean;
typedef short unschar;

unschar **approxarray;

typedef struct censorrec {
  boolean censor[maxcensor - minint];
} censorrec;

typedef char hapvector[maxlocus];
typedef double covmatrix[maxtrait][maxtrait];
typedef int mutarray[3][2];
typedef double thesemeans[maxtrait];
typedef thesemeans means[maxall + 1][maxall];
typedef enum {
  auto_,    /* autosomal, regular inheritance */
  mauto,    /* autosomal, mutation model */
  sex,      /* sexlinked, regular inheritance */
  msex      /* sexlinked, mutation model */
} pathway;
typedef enum {
  peelup, peeldown
} direction;

typedef int binset;

typedef binset phenarray[maxall];


typedef enum {
  affection, quantitative, binary_, null_
} locustype;

typedef struct locusvalues {
  int nallele, format;
#if LOOPSPEED /* Dylan added */
  int fgeno;  /*holds the number of single-locus genos possible*/
#endif
  double freq[maxall];
  double freqcond[maxall]; /*conditional frequencies for carriers*/
  struct locusvalues *privlocus;
  locustype which;
  union {
    struct {
      double pen[maxall + 1][maxall][3][maxliab];
      long nclass;
    } U0;
    struct {
      long ntrait;
      means pm;
      covmatrix vmat;
      double det, contrait, conmat;
    } U1;
    phenarray allele;
  } UU;
} locusvalues;

typedef struct phenotype {
  locustype which;
  binset phenf; /*encodes phenotype/genotype in pedin.dat as a
                  bit vector*/
/*Next two fields added by AAS to avoid hassles in dealing with phenf*/
  int allele1;
  int allele2;
  double x[maxtrait];
  boolean missing;
  int aff;
  long liability;
} phenotype;

typedef phenotype *hindphen[maxsys];
typedef phenotype *indphen[maxlocus];


/* a structure of type 'thisarray' stores the conditional
genotype probabilities for an individual. The
probabilities are stored in the array genarray.
We expect this array to be sparse. The field sparseflag
is a boolean array such that sparseflag[i] is nonzero
if and only if genarray[i] is nonzero. */

typedef struct thisarray {
  unsigned char *sparseflag;
  double *genarray;
} thisarray;

#if !LOOPSPEED
typedef boolean possvect[maxall][maxall];
typedef possvect possarray[maxlocus];

typedef struct information {
  possarray possible;
} information;
#endif

/* a record of type thisperson stores information about one person
   in one pedigree */

typedef struct thisperson {
  int id, ped, inloop;
  struct thisperson *pa;      /* father */
  struct thisperson *ma;      /* mother */
  struct thisperson *foff;    /* first offspring */
  struct thisperson *nextpa;  /* next (half)-sibling with the same father */
  struct thisperson *nextma;  /* next (half)-sibling with the same mother */
  thisarray *gen;             /* genotype information */
  hindphen holdphen;
  indphen phen;
  phenotype *privphen;
  boolean unknown, multi, done, up, male, firstpass;
  /* next three added by A. A. Schaffer */
  boolean newgenexists; 
  boolean loopdepend;  /* true if genotype of loopbreaker depends on this */
  boolean loopneeded;  /* true if genotype of this depends on loopbreaker */
#if LOOPSPEED /* Dylan added */
  boolean thisunknown[maxlocus];  /*true if unknown at this locus*/
#else
  information *store;
#endif
#if PARALLEL
  short memindex;    /*added to keep track of where gen is in genbank*/
#endif
} thisperson;



typedef int subhap[maxseg];
typedef double thetarray[maxlocus];
/* typedef double happrob[maxneed]; */
unsigned nuneed; /* Introduced by R. M. Idury, actual size of segprob
                    arrays */
unsigned nuprobclass; /* Introduced by A. A. Schaffer, actual number
			 of probclasses */

typedef struct thetavalues {
  thetarray theta;
  double *segprob;
} thetavalues;

thetavalues *maletheta, *femaletheta;
thetarray *gmaletheta, *gfemaletheta;

typedef struct new_locus_info {
  int present;  /*is this allele present?*/
  int new_allele; /*what is this allele mapped to?*/
  int old_allele; /*inverse of mapping*/
  double old_frequency;
  double new_frequency;
  double new_frequency_cond;
  double old_frequency_cond; /*for conditional allele frequencies, Morgan */
} new_locus_info;

#if ALLELE_SPEED
#if defined(LODSCORE)
typedef int new_allele_count[maxsys]; 
#else
typedef int new_allele_count[maxlocus]; 
#endif

/*Adjusted allele counts for each pedigree and locus*/
new_allele_count *ped_new_allele_count; 

/* Does this pedigree have different adjusted allele counts from
   the previous pedigree? */

boolean *ped_must_change_locations;
#endif
  
/* dwix: begin */
/* type definitions for allele_downcode_check */
#if defined(LODSCORE)
typedef new_locus_info  loc_all[maxsys][maxall+1];
#else
typedef new_locus_info  loc_all[maxlocus][maxall+1];
#endif

loc_all *ped_loc_all;
int currentped;
int *pedidx2num;
int mymaxped;
/* dwix: end */


int maxhaplo, maxfemgen; /*A. A. Schaffer*/
int maxclasssize, maxisozygclass;

/*The following declarations help store a partial correspondence
between genotypes and haplotypes,
haps1 between base[i] and fence[i] stores the left haplotypes
that genotype i can pass on to a child; haps2 stores the right
haplotypes; hind stores the index into an array of recombination probabilities
indicating the probability of this haplotype getting passed on from
genotype i. currentfence is used as a counter to fill
haps1, haps2, hind, base, and fence. */
unsigned currentfence;           /* R. M. Idury */
unsigned int *base, *fence; /*R. M. Idury*/
unsigned short *haps1, *haps2;
unsigned int *hind;
/* The arrays invgenenum1 and invgenenum2 store an effective inverse
to the genenumber mapping above. They convert an index into
genenumber into two haplotypes */
unsigned *invgenenum1, *invgenenum2; /* R. M. Idury */

typedef struct cache{
 unsigned first;
 unsigned last;
} cache;

 pathway thispath;
 boolean informative[maxped];
 boolean *rare, *risk1,  *risk2;
 boolean *riskmale;
 censorrec *censorstruct;
 int thisc;
 boolean malechild[maxchild];
#if !PARALLEL
 thisarray *thischild[maxchild];  /*if parallel must be allocated dynamically*/
#endif
 thisperson *childarray[maxchild]; /*A. A. Schaffer*/
 int nchild;
 int *segstart;
 unsigned int *probstart, *probend;
 boolean nohom[maxlocus];
 locusvalues *thislocus[maxlocus];
 int increment[maxlocus], order[maxlocus];
 unsigned int *nonzgens; /*used to hold nonzero genotypes in pedigree
                            traversal, A.A. Schaffer */
 boolean *flag;  /*R. M. Idury*/  
 double *gene; /*used for local computations in pedigree traversal routines*/

/* The arrays psumcache and qsumcache store conditional probabilities of
   different haplotypes being passed on from p and q respectively.
   They are used in some of the pedigree traversal routines. */
 double *psumcache, *qsumcache;

#if PARALLEL
 double *onechildupqsumcache;
#endif

/*Used in segup to keep track of entries in indpool, invpool,
  and nextpool that correspond to different haplotypes of a child*/
 cache *phapcache1;

/*PEOPLE*/
 thisperson *person[maxind + 1];
 thisperson *proband[maxped];
 thisperson *looppers[maxped][maxloop][maxloop];
 int numCopies[maxped][maxloop];  /*number of copies of loop breaker*/

/*MUTATION */
 unschar *muthap;
 int **genenumber;
/*RECOMBINATION*/
 unsigned *segindex;
 double *tempseg, *tempseg2;
 double *segval;
/*FREQUENCIES */
 thisarray *hapfreq;
/*RISK*/
 int riskall;
/*OTHERS*/
 int risksys, mutsys, mlocus, lastpriv, i;
 int nuhap;
 int fgeno, mgeno;
 int nuped, totperson;
 double segscale, mutmale, mutfemale, like, alike, distratio;
 boolean interfer, disequi, sexlink, risk, sexdif, readfemale,
	        dolod, firstapprox, lasttime;
 boolean firsttime; /* true if in or before first likelihood evaluation */
 boolean disfreqs; /*AAS for K. Morgan, tells if we are getting two
                     sets of allele frequencies*/
 FILE *outfile, *ipedfile, *datafile, *stream, *speedfile;

 FILE *loopfile;   /*holds information for loop speedup*/
#if LOOPSPEED /*Dylan and A. A. Schaffer  added */
 boolean ever_read_loopfile;  /*set true if loopfile has ever been read*/
 int *loopbreaker_nextgeno[maxloop]; /*sparse array of possible genotypes
                                       for loopbreakers*/
#endif
/*ILINK*/
 FILE *final;
 boolean mapping;
/*Added by A. A. Schaffer to keep track of which traversal we are on
 for handling a loop*/
 char loopfirstgen; /* false if doing a traveral for the loopbreaker
		       and it is not the first genotype of the innermost
		       loop breaker */
 char looplastgen;

#if LOOPSPEED  /* Dylan added */
 /*
   These typedefs define the table to hold the valid loopbreaker vectors
   for each locus.
 */
 typedef int *vector_for_loop;
 typedef vector_for_loop *vector_for_locus;
 typedef vector_for_locus *loop_vector_array;

 /*
   These typedefs define the table to hold the possible genotypes
   for each unknown person at each locus for each valid loopbreaker
   vector.
 */
 typedef boolean *geno_for_loop_vector;
 typedef geno_for_loop_vector *geno_for_locus;
 typedef geno_for_locus *geno_for_unknown;

 geno_for_unknown *unknown_poss;
 loop_vector_array loop_vectors;
 int num_loop_vectors[maxlocus];
 int num_loops[maxped];            /* number of loops in pedigree */
 int single_locus_vector_num[maxlocus];

 boolean **is_zero_breaker; /*AAS, identifies loop breaker vectors with
                              0.0 likelihood*/
 boolean *breaker_poss_genotype[maxloop];

#endif
 int *ped_nuscales;        /*AAS, keeps track of scaling factors used
                             in each pedigree */

 int memcount, maxmemcount; /*AAS, used to count how many genarrays needed for
                    memory estimation*/

#if PARALLEL

#if IS_SHMEM
#if IS_PARMACS
typedef struct Barrier
{
  BARDEC(barrier)
} Barrier;
#else  /* IS_PARMACS */

/* cgh -- p4 barrier syntax */
typedef struct Barrier
{
  p4_barrier_monitor_t barrier;
} Barrier;

#define BARINIT(bar)      p4_barrier_init(&(bar))
#define BARRIER(bar, np)  p4_barrier(&(bar), (np))

#endif  /* IS_PARMACS */
Barrier **partial_barrier;
#endif  /* IS_SHMEM */

typedef struct GlobalMemory
{
  int finished;

  /* cgh && schaffer -- got rid of distinction between
     PARALLEL_GCENTRAL and not, and changed maxn to maxprocs.
     Now indexed by mymaster. */
  int nchild[maxprocs];
  int trav_flag[maxprocs];   	/* indicates whether in segup or segdown */
  int pmale[maxprocs];       /* indicates whether parent p is male */

  			/* indicates whether there's one or many children */
#if IS_PARMACS
  BARDEC(full_barrier)
#elif IS_P4
  Barrier full_barrier;
#else
  int seg_barrier;
  int sum_barrier;
#endif
  int pprocbase[maxprocs];
  int pprocfence[maxprocs];
  int qprocbase[maxprocs];
  int qprocfence[maxprocs];

  /* cgh && schaffer -- got rid of distinction between
     PARALLEL_GCENTRAL and not, and changed maxn to maxprocs.
     Now indexed by mymaster. */
  thisarray *pgenptr[maxprocs]; /*first parent */
  thisarray *qgenptr[maxprocs]; /*sceond parent */
  thisarray *rgenptr[maxprocs]; /*child*/
  thisarray *gthischild[maxprocs][maxchild];
} GlobalMemory;


GlobalMemory *gMem;

/* Global variables added by Sandeep for Phase II: Parallel Thetas */
int *checkmaster;     /* Boolean indicating if iAmAMaster should be checked */
int *usecentral;      /*Are we using central differences*/
int *slavesPerGroup;  /* number of slaves per parallel theta group */
int *thetasPerGroup;  /* number of thetas that are calculated by each group */
int *nextTheta;       /* next theta number that has not been assigned */
int *numGroups;       /* number of parallel theta groups */
int *thetanumbase;    /* element i contains number of first job to be done */
                      /* by the group led by processor i */
int *thetanumfence;   /* element i contains 1 + number of last job to be done */
                      /* by the group led by processor i */
int currentthetanum;  /* number of the theta for the current function eval */
int mymaster;         /* processor number of job's master */
int parallelThetas;   /* bollean indicating if doing fun calls in parallel */
int infun;            /* bollean indicating if we are in fun() */
int *funfinished;     /* array of booleans indicating if fun is finished */
#if ALLELE_SPEED
int *sharedCurrentPed; /* array  of current pedigree numbers*/
#endif

double **gmalesegprob;
double **gfemalesegprob;
double *gf;           /* global f */
#if PARALLEL_GCENTRAL
double **gcentralf;   /* f's calculated during gcentral */
#endif /* PARALLEL_GCENTRAL */
int barnum;           /* Number of barriers (used in debugging) */
#if FUNTIMES
#if PARALLEL_GCENTRAL
double *executionTimes; /* Records times of executions during parallelThetas */
#else /* !PARALLEL_GCENTRAL */
double **executionTimes; /* Records times of executions during parallelThetas */
int fe;      /* first or second fe during gforward, gcentral */
#endif /* PARALLEL_GCENTRAL */
#endif /* FUNTIMES */
double accumfamilytime;


/* Global variables needed by Phase I: Load Balance */
double *rowtime;      /* row times calculated by processor 0 on the */
double *qrowtime;     /* first round of parallel thetas */

typedef int fencearray[MAXNUMFAM][maxprocs];

fencearray *rowfence; /* Stores fence values (numnuc x numproc) */
fencearray *qfence;   /* Stores fence values (numnuc x numproc) */
fencearray *temprowfence; /* Temporary storage for when fences are calculated */
fencearray *tempqfence;   /* Temporary storage for when fences are calculated */
int *timeExecutions;  /* boolean that indicates if execution times should be
                         stored in [q]rowtime.  Only used if mymaster = 0 */
int *numNuclear;      /* The number of nuclear families */
int *firstfe;           /* boolean that indicates if first function evaluation */
int *secondfe;        /* boolean that indicates if second function evaluation */
int onlymaster;       /* indicates that only master is doing family */

/*The array genbank keeps a set of pointers to gen arrays that can be
used dynamically; geninuse will keep track of which 
entries in genbank are in use*/

thisarray *genbank;
thisarray **ggenbank;  /* cgh - used to be [maxn] or [2*maxn] */

thisarray **thischild;

char *geninuse;
char **ggeninuse;  /* cgh - used to be [maxn] or [2*maxn] */

double **arrgene;
double **arrqsumcache; /* added to provide individual global qsumcaches that
			* can be summed - Sandhya
			*/
double **arrpsumcache; /* added to provide individual global psumcaches that
			* can be summed - Sandhya
			*/

unsigned int *pnonzgens, *qnonzgens, *stripe_pnonzgens, *stripe_qnonzgens,
    *privatepnonzgens, *privateqnonzgens;
unsigned int **gpnonzgens, **gqnonzgens;  /* cgh - used to be [maxn] or [2*maxn] */
#endif /*PARALLEL*/



/* Local variables for seg: */
struct LOC_seg {
  struct LOC_likelihood *LINK;  /* information for pedigree */
  thisperson **p;           /* in segdown(): father
			       in segup(): parent connecting traversal */
  thisperson **q;           /* mate of p */
  thisperson **r;           /* first child of p and q */
  thisperson *child, *father, *mother;  /* same people as p, q, and r */
  int fseg, sseg, sstart, send, fstart, fend, firstseg, secondseg;
  int nfirst;               /* number of genotypes of gender of p */
  int nsecond;              /* number of genotypes of gender of q */
  double pf, ps;
  thetavalues *firstsex;    /* theta info for gender of p */
  thetavalues *secondsex;   /* theta info for gender of q */
} ;

/*Local variables for likelihood*/
struct LOC_likelihood {
  int thisped;
  thisperson *proband;
  int loopgen[maxloop];
  double homo, hetero;
  int nuscale;
  thisarray *holdpoint[maxloop];
} ;

/*Local variables for inputdata*/
#if defined(LODSCORE)
struct LOC_inputdata {
  long nfactor[maxsys];
};
#else
struct LOC_inputdata {
  long nfactor[maxlocus];
};
#endif

/*Local variables for readloci*/
struct LOC_readloci {
  struct LOC_inputdata *LINK;
  long i, whichtype, nupriv;
};

/*The following definitions are needed to do checkpointing*/
/* K. Shriram: begin */

#define DateTimeStampStringLength       27      /* see ctime()'s man page */

#define SystemCallStringLength  1024    /* length of arg for system() */
char  systemCallString [ SystemCallStringLength ] ;

#define ensureWrite(f)          { fflush ( (f) ) ; fsync ( fileno ( f ) ) ; }

FILE  * checkpointDatafile ;

/* K. Shriram: end */

/* begin: cgh */

/* program name */
#if defined(MLINK)
#define PROGRAM "MLINK"
#elif defined(LINKMAP)
#define PROGRAM "LINKMAP"
#elif defined(ILINK)
#define PROGRAM "ILINK"
#else
#define PROGRAM ""
#endif  /* defined(MLINK) */

/* min, max macros */
#define Max(a, b) ((a) > (b)) ? (a) : (b)
#define Min(a, b) ((a) < (b)) ? (a) : (b)

/* used to avoid round-off errors with floats */
#define ROUNDOFF_INC  0.0001

/* no reason to call sprintf 35 times for this */
#define LINE "-----------------------------------"


/* Parallel definitions */
#if PARALLEL

#if IS_SHMEM
#define Tmk_distribute(a,b);
#endif  /* IS_SHMEM */

/* cgh -- new variables for MLINK and LINKMAP */

/* Global private variable for MLINK and LINKMAP indicating which
   _absolute_ theta vector the current processor is working on.  This
   is in contrast to currentthetanum which indicates which theta vector
   _relative_ to the current group of parallel evaluations the current
   processor is working on. */
int absoluteThetanum;

/* Global shared array indicating reordering of theta vectors so
   all nonzero vectors appear first.  Entry i = k means that for the
   ith evaluation, absoluteThetanum is k. */
int* absoluteThetaIndex;

/* Global shared variable indicating which theta is the first among
   those being worked on in parallel */
int* firstThetanum;

/* Global private variable for processor 0 indicating which relative
   theta vector should be considered next in AssignThetas() */
int nextThetaToAssign;

/* Total number of calls to iterpeds() calculated in simIpedLoop() */
int numIpeds;

/* Global variable indicating the number of theta vectors that
   do contain zero entries -- valid only for processor 0 */
int numZeroThetas;

/* Global variable indicating the number of theta vectors that
   contain no zero entries -- valid only for processor 0 */
int numNonzeroThetas;

/* array of booleans indicating whether or not the corresponding entry
   in g{male,female}theta has a zero entry */
boolean* gZeroTheta;

/* Array of booleans indicating if processor is master */
boolean* iAmAMaster;

/* element i contains master of processor i */
int* whoIsMyMaster;  

/* Unlinked likelihood evaluation --
   corresponds to scorevalue in MLINK,
   second variable stores score for each family 
   for lodbyfamily */
double* unlinkedLike;
double* unlinkedLikeByPed;

/* Status of iterpeds() loop simulation */
typedef enum { countIpeds, computeThetas } loopStatus;

#include "strbuff.h"

#define TEMPBUFF_SIZE 256    /* size of temporary output buffer */

/* outBuff structure
   An array of strBuff objects.  There is one for each iteration of
   iterpeds(). */
typedef struct outBuff {
  int numIpeds;             /* number of calls to iterpeds() */
  strBuff** ipeds;         /* array of all iterpeds() output */ 
} outBuff;

/* global outBuff objects for stdout, outfile, and stream */
outBuff* stdoutBuff;
outBuff* outfileBuff;
outBuff* streamBuff;

/* flag which dictates whether or not to report memory usage stats */
boolean reportMemStats;

#endif  /* PARALLEL */

#if 0   
/* Unfortunately, the following breaks cc.  We'll put this off until
   we can figure out a fix. -- cgh */
/* The VAX requires an `l' with an `e' in the print format string;
   however, on other platforms, gcc complains about this use.  We use
   this macro instead, and rely on compile-time string concatenation
   to handle it for us. */
#if defined(vms)
#define pt5e ".5le"
#else 
#define pt5e ".5e"
#endif  /* if defined(vms) */
#endif  /* if 0 */


/* function prototypes */
#if !defined(KNR_PROTO)

#if PARALLEL

void writeBuff(outBuff*, FILE*);
void bufferPedOutput(int, double);
void bufferTotalsOutput(double, double);
void bufferLodOutput(double, double, boolInt);
void printBarrier(char*, int);

#if defined(MLINK)
void bufferLikeOutput(double, double);
#else  /* if defined(MLINK) */
void bufferLikeOutput(double);
#endif  /* if defined(MLINK) */

#endif  /* PARALLEL */

void malloc_err(char*);
int P_eof(FILE*);
int P_eoln(FILE*);
unsigned MakeMask(int start, int end);
void getlocus(long, struct LOC_readloci*);
void inputerror(long, long, long, struct LOC_inputdata*);
void inputwarning(long, long, long, struct LOC_inputdata*);
void readped(struct LOC_inputdata*);

#if (defined(ILINK) || defined(LINKMAP) || defined(MLINK))
void gettheta(thetavalues**, struct LOC_readloci*);
#endif  /* (defined(ILINK) || defined(LINKMAP) || defined(MLINK)) */

#endif  /* !defined(KNR_PROTO) */

#if defined(KNR_PROTO)
void seqStartup();
#else
void seqStartup(int, char**);
#endif

#if PARALLEL

outBuff* newOutBuff();
void initOutBuffs();
void writeOutBuffs();
void preLikeBufferOutput();

#if (defined(MLINK) || defined(LINKMAP))
int calcNumparams();
#endif /* (defined(MLINK) || defined(LINKMAP)) */

#endif  /* PARALLEL */

void checkzero();
void  getCkptTuple();
void allocgen();
void allocprob();
void allocthetas();
void manychilddown();
void manychildup();
void onechilddown();
void onechildup();

/* cgh - routines in commoncode.c */
#if !defined(KNR_PROTO)
double mapfunction(double, double);
double getdist(double*);
double invdist(double*);
#endif  /* !defined(KNR_PROTO) */
void initialize();

/* cgh -- routines in commoncode.c */
void openFiles();
void ckptInit();
void ckptEnsureFileContents();
void closeInputFiles();
void miscInit();
void ckptPreIpedRecover();
void closeOutputFiles();
void ckptCleanup();
void initParams();

/* prototypes for SunOS 4.1.x gcc */
#if defined(SUNOS_4_1_X_GCC)
int _filbuf();
int _flsbuf();
int close();
int fclose();
int fflush();
int fprintf();
int fread();
int fscanf();
int fseek();
int fsync();
int fwrite();
int getopt();
int gettimeofday();
int printf();
int puts();
int rename();
void rewind();
void setbuf();
time_t time();
int ungetc();
int vfprintf();
#endif  /* defined(SUNOS_4_1_X_GCC) */


#if defined(KNR_PROTO)

#if PARALLEL
void writeBuff();
void bufferPedOutput();
void bufferTotalsOutput();
void bufferLodOutput();
void bufferLikeOutput();
void printBarrier();
#endif  /* PARALLEL */

#if (defined(ILINK) || defined(LINKMAP) || defined(MLINK))
void gettheta();
#endif  /* (defined(ILINK) || defined(LINKMAP) || defined(MLINK)) */

double mapfunction();
double getdist();
double invdist();

#endif  /* defined(KNR_PROTO) */


#endif  /* _COMMONDEFS_H */

extern void getprobtable();
extern void segup();
extern void segdown();
extern void segsexup();
extern void segsexdown();
extern void performCheckpoint();
extern void recoverCheckpoint ();
extern void getvect();
extern void recombination();
extern void getlocations();
extern void readspeed();
extern void invert();
extern void cleanup();
extern void initseg();
extern void exitseg();
extern void allocategenetables();
extern void getgeneindices();
extern void freegenetables();
#if ALLELE_SPEED
extern void recompute_haps();
#endif
extern void seg();
extern void likelihood();
extern void allocate_loopbreaker_vectors();
extern int maxw_estimation();
extern void segsexctop();
extern void segsextop();
extern void segctop();
extern void segtop();
extern void segcapprox();
extern void msegsexdown();
extern void msegdown();
extern double mapfunction();
extern double getdist();
extern double invdist();
extern void setparam();
extern void inputdata();
extern void oldsegsexup();
extern void oldsegup();
extern void malloc_err();
extern int P_eof();
extern int P_eoln();
extern void check_constants();
extern void printErrMesg();
extern void init_ped_loc_all();
#if ALLELE_SPEED
extern int adjust_alleles();
extern void allele_adjust_persons();
#endif
extern int allele_downcode_check();
extern unsigned MakeMask();
extern void getlocus();
extern void inputerror();
extern void inputwarning();
extern void open_loop_file();
extern void close_loop_file();
extern void read_loop_file();
extern void readped();
extern void allocate_thisarray();
