/*This file contains definitions for new versions of the routines */
/*segup, segdown, segsexup, segsexdown, and some auxiliary routines */
/* for use with LODSCORE, ILINK, LINKMAP, and MLINK. */
/* The modifications are described in the papers: */
/* R. W. Cottingham, Jr., R. M. Idury, and A. A. Schaffer, */
/* Faster Sequential Genetic Linkage Computations */
/* American Journal of Human Genetics, 53(1993) pp. 252-263*/
/* and A. A. Schaffer, S. K. Gupta, K. Shriram, and R. W. Cottingham, Jr. */
/* Avoiding Recomputation in Linkage Analysis, */
/* Human Heredity 44(1994), pp. 225-237. */
/* The versions using this file use a lot of memory */
/* Definitions for somewhat slower versions that use less memory */
/* are in slowmoddefs.h */

#ifndef  _MODDEFS_H

#define  _MODDEFS_H  1



#ifndef AUTOSOMAL_RUN
#define AUTOSOMAL_RUN  1  /*1 is safe; can be set to 0 to make
                            sexlinked runs space efficient*/
#endif

#ifndef SEXDIF_RUN
#define SEXDIF_RUN     1  /*1 is safe; can be set to 0 on sexlinked runs or
                            runs in which maletheta and femaletheta will be
                            the same*/
#endif

/*segprob2 stores products of recombination or
nonrecombination probabilities; one probability
comes from each parent */
double *segprob2;



/*used in the autosomal case when sexdif=1*/
#if PARALLEL
double **probtabledif;
double **probtable;
unsigned **probtableindex;
#else
double *probtabledif;
double *probtable;
unsigned *probtableindex;
#endif

unsigned *classsize;
unsigned *classbase;

typedef double childprob[maxchild]; 
childprob **partialprob;

unsigned *invpool, *nextpool, *indpool;


/* _MODDEFS_H */

#endif
