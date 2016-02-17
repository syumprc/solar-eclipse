/* This file contains new versions of the routines */
/* segsexup, segsexdown, and some auxiliary routines for use with */
/* LODSCORE, ILINK, LINKMAP, and MLINK */
/* The modifications are described in the papers: */
/* R. W. Cottingham, Jr., R. M. Idury, and A. A. Schaffer */
/* Faster Sequential Genetic Linkage Computations */
/* American Journal of Human Genetics, 53(1993), pp. 252-263*/
/* and A. A. Schaffer, S. K. Gupta, K. Shriram, and R. W. Cottingham, Jr. */
/* Avoiding Recomputation in Linkage Analaysis, */
/* Human Heredity 44(1994), pp. 225-237. */
/* The versions in this file use a lot of memory */
/* Somewhat slower versions that use less memory */
/* are in slowsexmodified.c */

/* 2-8-94; fixed two bugs AAS*/

/*Most of the code in this file was written by R. M. Idury */
#include "commondefs.h"
#include "moddefs.h"

unsigned lsegsexfun2(first,second)
unsigned first, second;
{
  int g1, g2; /*two gene numbers*/
  int j, k; /* loop indices*/
  int s1, s2; /*haplotype numbers*/
  int FORLIM1, FORLIM2; /*loop bounds*/
  int *TEMPGENE1; /* store pointers into genenumber*/
  unsigned char *tempflag3; /*stores sparsity pattern for child's genarray*/

  FORLIM1 = base[second];/*find beginning and end of isozygote class
			   for second*/
  FORLIM2 = fence[second];
  TEMPGENE1 = genenumber[first];/*get pointer into genenumber for
					 this haplotype*/
  /*try to find a non-zero value for each child*/
  for (k = 0; k < nchild; k++) {

    /*retrieve sparsity pattern for child k*/
    tempflag3 = thischild[k]->sparseflag;
  
    /*iterate over all the recombined isozygotes of second*/
    for (j = FORLIM1; j < FORLIM2; j++) {
      s1 = haps1[j]; /*get haplotypes of this genotype*/
      s2 = haps2[j];
      
      if(malechild[k]) {
	if(tempflag3[s1 - 1] != 0 || tempflag3[s2 - 1] != 0)
	  goto notzero;
	else
	  continue;
      }
      g1 = TEMPGENE1[s1 - 1];
      g2 = TEMPGENE1[s2 - 1];
      
      /*if any of the flags is TRUE, then this child can have this genotype,
	and we move immediately to testing the next child*/
      if(tempflag3[g1 - 1] != 0 || tempflag3[g2 - 1] != 0) goto notzero;
    }
    return 0; /*this child failed to have any of the genotypes for
		the isozygote class, so this isozygote class is not possible*/
  notzero:
    continue;
  }
  return 1; /*all children passed the test and have a possible joint genotype*/
}

/* segsexsum2 is used in segsexup to compute some common subexpressions in
the probability updates.
first and second are the joint genotypes of the two parents.
fslength is the number of probabilities needed for the combined
isozygote class of both parental genotypes; fs stands for the product
of first and second. LINK is used to pass the genetic information for
parents and children */

Local void segsexsum2(first,second,fslength)
unsigned first,second, fslength;
{
  int g1, g2; /*indices to store gene numbers*/
  int f1, s1, s2; /*indices to store haplotype numbers*/
  int index; /*counter for combined isozygote equivalence class,
                will range from 0 to fslength-1*/
  int i, j, k, l; /*loop indices*/
  int FORLIM1, FORLIM2; /*lower and upper limits on a for loop*/
  int *TEMPGENE1; /*temporary pointers into genenumber table*/
  double *tempwith3; /*temporarily stores genarray for current child*/

  FORLIM1 = base[second]; /*start of the isozygote class for second*/
  FORLIM2 = fence[second]; /*end of the isozygote class for second*/
  index = 0; /*start with first isozygote in the class*/
  
  i = base[first];
  f1 = haps1[i]; /*retrieve first haplotype of joint genotype i*/
  TEMPGENE1 = genenumber[f1 - 1]; /*lookup partial information
					    in genenumber table*/
  /*iterate over isozygotes of j*/
  for (j = FORLIM1; j < FORLIM2; j++) {
    s1 = haps1[j]; /*retrieve first haplotype of joint genotype j*/
    s2 = haps2[j]; /*retrieve second haplotype*/
   /*lookup the two ways to combine one haplotype from i and one from j*/
    g1 = TEMPGENE1[s1 - 1];
    g2 = TEMPGENE1[s2 - 1];

    /*iterate over children; update partial computations of
      new probabilities and store in array tempseg
      note that tempseg has nchild*fslength relevant entries
      The fslength entries for child 1 come first, then
      the fslength entries, for child 2, etc. This is why
      the increment on l is fslength for each change of child*/

      for (l = index, k = 0 ; k < nchild; l += fslength, k++) {
	tempwith3 = thischild[k]->genarray; /*retrieve genarray*/

       /* sum the probabilities for the two joint genotypes that
        the child can get from the current parental joint genotypes*/

	if(malechild[k])
	  tempseg[l] = (tempwith3[s1 - 1] + tempwith3[s2 - 1]);
	else
	  tempseg[l] = (tempwith3[g1 - 1] + tempwith3[g2 - 1]);
      }
      index++; /*increment isozygote class counter*/
    }
}

/*segsexsumdown2 is used in segsexdown to compute some common subexpressions in
the probability updates.
first and second are the joint genotypes of the two parents.
fslength is the number of probabilities needed for the combined
isozygote class for both parental genotypes; fs stands for the
product of first and second. */

Local void segsexsumdown2(first,second,fslength,male)
unsigned first,second,fslength;
unsigned char male;
{
  int g1, g2; /*indices to store gene numbers*/
  int f1, s1, s2; /*indices to store haplotype numbers*/
  int index; /*counter for combined isozygote equivalence class,
                will range from 0 to fslength-1*/
  int index2; /*used as index for isozygote class*/
  int i, j, k, l; /*loop indices*/
  int FORLIM1, FORLIM2; /*lower and upper limits on a for loop*/
  int *TEMPGENE1; /*temporary pointers into genenumber table*/
  double *tempwith3; /*temporarily stores genarray for current child*/

  FORLIM1 = base[second]; /*start of the isozygote class for second*/
  FORLIM2 = fence[second]; /*end of the isozygote class for second*/
  index = 0; /*start with first isozygote in the class*/
  index2 = 0;
  /*iterate over isozygotes of i*/
  i = base[first];
  f1 = haps1[i]; /*retrieve first haplotype of joint genotype i*/
  TEMPGENE1 = genenumber[f1 - 1]; /*lookup partial information
					    in genenumber table*/
  /*iterate over isozygotes of j*/
  for (j = FORLIM1; j < FORLIM2; j++) {
    s1 = haps1[j]; /*retrieve first haplotype of joint genotype j*/
    s2 = haps2[j]; /*retrieve second haplotype*/
    /*lookup the two ways to combine one haplotype from i and one from j*/
    g1 = TEMPGENE1[s1 - 1];
    g2 = TEMPGENE1[s2 - 1];
    /*store these gene numbers for later use; the
      key point is that we will want these numbers
      consecutively later, so we store them consecutively
      in segindex*/
    if(male) {
      segindex[index2++] = s1;
      segindex[index2++] = s2;
    }
    else {
      segindex[index2++] = g1;
      segindex[index2++] = g2;
    }

    /*iterate over children; update partial computations of
      new probabilities and store in array tempseg
      note that tempseg has nchild*fslength relevant entries
      The fslength entries for child 1 come first, then
      the fslength entries, for child 2, etc. This is why
      the increment on l is fslength for each change of child*/

      for (l = index, k = 0 ; k < nchild; l += fslength, k++) {
	tempwith3 = thischild[k]->genarray; /*retrieve genarray*/

       /* sum the probabilities for the two joint genotypes that
        the child can get from the current parental joint genotypes*/

	if(malechild[k])
	  tempseg[l] = (tempwith3[s1 - 1] + tempwith3[s2 - 1]);
	else
	  tempseg[l] = (tempwith3[g1 - 1] + tempwith3[g2 - 1]);
      }
      index++; /*increment isozygote class counter*/
    }
}

/*segsexup updates the genotype probabilities for a parent (p) based
on the probabilities of children. The parameter LINK includes
all the genetic information about p, p's spouse q, and their
children*/

Void segsexup(LINK)
struct LOC_seg *LINK;
{
  int  FORLIM; /*loop bound*/
  int nonzindex, nonzcount;    /*loop index and counter for nonzero values*/
  unsigned int  step1, step2; /*size of isozygote classes*/
  static unsigned current; /*index for haplotype pools*/
  double valtemp; /*intermediate value in probability updates*/
  double val, temp1, tempr; /*temporaries to store intermediate values*/
  unsigned int f1, f2; /*two haplotypes from parent*/
  unsigned int here, i, j, first, second; /*genotype indices*/
  int jointisoindex; /*index to work within isozygote classes*/
  unsigned k; /*index to loop over children*/
  unsigned char skip; /*used to skip iterations of a loop*/
  thisarray *WITH2; /*stores genetic information about p*/
  thisarray *WITH3; /*stores gnetic information about q*/
  double *newwith2, *newwith3, *newwithr; /*store genarrays for
						       p,q, and children*/
  unsigned int c1, c2; /*haplotypes*/
  boolean depend; /*used to handle loops*/
  unsigned char *newflag2, *newflag3, *newflagr; /*store sparsity patterns for
					p and q genarrays*/


  /* newsegprob and newsegprob2 are  used to hold segprob
    arrays, which contain the probabilities for various patterns of
    recombination events*/
  double *newsegprob, *newsegprob2; 
  double *tempprob; /*temporary holder for probability array*/
  int ind;      /*used to store offset for probability array*/

  initseg(LINK); /*get data about this p,q,children triple*/

  /*get sparsity patterns for p and q genarrays*/
  newflag2 = (*LINK->p)->gen->sparseflag; 
  newflag3 = (*LINK->q)->gen->sparseflag;

  WITH2 = (*LINK->p)->gen; /*get genetic data for p*/
  WITH3 = (*LINK->q)->gen; /*get genetic data for q*/
  newwith2=WITH2->genarray; /*get genarray for p*/
  newwith3=WITH3->genarray; /*get genarray for q*/
  newwithr=thischild[0]->genarray; /*get genarray for first child*/

  if((*LINK->p)->male) {
/*The case of 1 child is handled specially because the subcomputations
are much simpler. In particular we do not need to multiply the
probabilities among the different children. In a typical pedigree
many pairs will have only one child about which information is known
(in essence this is the child that keeps the pedigree connected), so
this is an important and common special case. */

if(nchild == 1) {
/*initialize cache data structures*/
  for(i = 0;i < maxhaplo;i++) {
    flag[i] = 0;
    phapcache1[i].first = 0;
    qsumcache[i] = 0.0;
  }

/*initialize gene array and set up flag array*/
  if(malechild[0]) {
    for(i = 0; i < mgeno; i++) {
      gene[i] = 0.0;
      if(newwithr[i] == 0.0) continue;
      flag[i] = 1;
    }
  }
  else {
    for(i = 0; i < fgeno; i++) {
      gene[i] = 0.0;
      if(newwithr[i] == 0.0) continue;
      flag[invgenenum1[i]-1] = 1;
      flag[invgenenum2[i]-1] = 1;
    }
  }

  /*This section of the code precomputes for each haplotype the the probability
  that the child will inherit this haplotype from q. Each genotype has
  two haplotypes, but can produce others by recombination. Therefore,
  for each genotype we must sum over the different haplotypes that can
  be produced by its isozygote class. The contributions for each haplotype
  are stored in qsumcache.*/

  newsegprob = femaletheta->segprob; /*get probabilities to work with q*/

 /*iterate over all joint genotypes*/
  for (first = 0; first < fgeno; first++ ) {
    if(newflag3[first] == 0) continue; /*check if joint genotype is possible*/
    valtemp = newwith3[first]; /*get cond. prob. that q has genotype first*/
    FORLIM = fence[first]; /*find bounds for the isozygote class of first*/

 /*iterate over all the recombined genotypes of the isozygote class*/
    for (i = base[first]; i < FORLIM; i++) {
      f1 = haps1[i];
      f2 = haps2[i];
      if((flag[f1-1] != 0) || (flag[f2-1] != 0)) {
/*conditional probability of first as a genotype multiplied by probability
  of passing on f1 (alternatively f2) as a haplotype*/
        val = valtemp*newsegprob[hind[i]];

       /*store in qsumcache*/
        if(flag[f1-1] != 0) {
	  qsumcache[f1 - 1] += val;
        }
        if(flag[f2-1] != 0) {
	  qsumcache[f2 - 1] += val;
        }
      }
    }
  }


 /*In this section of the code we update the probabilities for
  the parent p based on the probabilities for the child*/  
 /*Iterate over all joint genotypes of the child*/
  if(malechild[0]) {
  newflagr=thischild[0]->sparseflag;
  temp1 = 0.0;
  for(here = 0; here < mgeno; here++) {
    if(newflagr[here] == 0) continue;
    tempr = newwithr[here];
    temp1 += tempr * qsumcache[here]; /*Fixed bug here 2-10-94, AAS*/
  }

  for(here = 0; here < mgeno; here++)
    gene[here] = temp1;
  }
  else {
  newflagr=thischild[0]->sparseflag;
  for(here = 0; here < fgeno; here++) {
    if(newflagr[here] == 0) continue;
    tempr = newwithr[here];
    c1 = invgenenum1[here];
    c2 = invgenenum2[here];
    gene[c1 - 1] += tempr * qsumcache[c2 - 1];
    if(c1 != c2)
      gene[c2 - 1] += tempr * qsumcache[c1 - 1];
  }
  }

/*set up new genarray for p; it is gene scaled by segscale*/      
  for (first = 0; first < mgeno; first++) {
    if(newflag2[first] == 0.0) continue;
    newwith2[first] *= segscale*gene[first];
  }

}


else  {  /*nchild is bigger than 1*/

/*find nonzero entries in q's genarray and make a list of them
  stored in nonzgens; just get one per isozygote class*/
  nonzcount = 0;
  
  /*iterate over genotypes for q*/
  for(i = 0; i < fgeno; i += step2) {
    /*number of distinct probabilties needed for i's isoz. class*/
    step2 = probend[i] - probstart[i] + 1;
    for(j = i; j < i+step2; j++)
      if(newflag3[j] != 0) {
        nonzgens[nonzcount++] = i; /*store index of nonzero value*/
	break;                        /*go to next isozygote class*/
      }
  }

  newsegprob2 = maletheta->segprob; /*get p probabilties*/

  /*iterate over genotypes for p*/
  for (first = 0; first < mgeno; first++) {

    if(newflag2[first] == 0) continue;

    /*initialize update multiple for each isozygote in class*/
    segval[0] = 0.0;

    /*iterate over the genotypes representing isozygote classes that
     q may have*/
    for (nonzindex = 0; nonzindex < nonzcount; nonzindex++) {
      second = nonzgens[nonzindex];
      /*check if this first, second pair yield a nonzero value
        among children's probabilities*/
      if(lsegsexfun2(first,second) == 0) continue;

      /*number of distinct probabilties needed for second's isoz. class*/
      step2 = probend[second] - probstart[second] + 1;
  
   /*call segsum2 to compute the part of the conditional probability
     update that is common to all the members of the combined isozygote
     class defined by first and second */
      segsexsum2(first,second,step2);

      /*now specialize update for each member of first's class*/
        /*further specialize update for each member of second's isozygote
          class*/
	for(j = 0; j < step2; j++) {
        /*skip if this isozygote not possible*/
	  if(newflag3[second+j] == 0) continue; 

          /*get offset into newsegprob2; the offset depends on
            the isozygote class size and index of each parent
            note that fisozygoteindex gets incremented by the size
            of the joint class each time, so it represents a
            sum of all the numbers of distinct probabilities
            needed for all joint iso. classes  considered before the
            current p isozygote class*/
	  tempprob = newsegprob2 + probstart[second+j] - 1;

        /*combine for all children*/
        /*due to the arrangement in segsexsum all the probability contributions
          for a given child are contiguous in the tempseg array.
          the number of contributions is fslength which is the number of
          probabilities needed for the joint isozygote class of the parents.
          We get the contribution of the first child (index 0) and
          then loop over the rest*/
	  val = tempprob[0] * tempseg[0];
         /*iterate over joint isozygote class*/
	  for(jointisoindex = 1; jointisoindex < step2; jointisoindex++)
	    val += tempprob[jointisoindex] * tempseg[jointisoindex];
	  newwithr = tempseg + step2; /*move the base of newwithr*/

         /*iterate over all the remaining children; notice that the
           loop index k is not used inside the loop this is because the
           last loop statement which moves the offset of newwithr
           has the effect of incrementing the child*/
	  for(k = 1; k < nchild; k++) {
	    temp1 = tempprob[0] * newwithr[0];
	    for(jointisoindex = 1; jointisoindex < step2; jointisoindex++)
	      temp1 += tempprob[jointisoindex] * newwithr[jointisoindex];
	    val *= temp1;
	    newwithr += step2;
	  }
         /*update segval entry for this isozygote of first*/
	  segval[0] += newwith3[second+j] * val;
	}
    }
    /*update p's genarray for each isozygote in this class*/
    newwith2[first] *= segval[0] * segscale;
  }
}

/* If any of the nonzero entries in p's genarray became 0,
   we want to set them to zero to avoid computations on subsequent
   calls*/

  for(i = 0; i < mgeno; i++)
    if((newflag2[i]==1) && (newwith2[i] == 0.0)) newflag2[i] = 0;
  }
  else { /* p is female */

if(nchild == 1) {
/*initialize cache data structures*/
  for(i = 0;i < maxhaplo;i++) {
    flag[i] = 0;
    phapcache1[i].first = 0;
    qsumcache[i] = 0.0;
  }

/*initialize gene array and set up flag array*/
    for(i = 0; i < fgeno; i++)   /*Fixed bug here 2-8-94, AAS*/
      gene[i] = 0.0;
  if(malechild[0]) {
    for(i = 0; i < mgeno; i++) {
      if(newwithr[i] == 0.0) continue;
      flag[i] = 1;
    }
  }
  else {
    for(i = 0; i < fgeno; i++) {
      if(newwithr[i] == 0.0) continue;
      flag[invgenenum1[i]-1] = 1;
      flag[invgenenum2[i]-1] = 1;
    }
  }

  current = 1; /* current should start with 1, since 0 means a null ptr */

  for (first = 0; first < fgeno; first++) {
    if(newflag2[first] == 0) continue; /*use only possible genotypes*/
    FORLIM = fence[first]; /*find end of isozygote class of first*/

  /*iterate over all members of first's isozygote calss*/
    for (i = base[first]; i < FORLIM; i++) {
      f1 = haps1[i]; /*get haplotypes*/
      f2 = haps2[i];
      ind = hind[i];  /*get probability offset for i*/

   /*if f1 is a possible haplotype for child add the
     probability offset and value of first to the cache*/
      if(flag[f1-1] != 0) {
	indpool[current] = ind; /*store values in cache*/
	invpool[current] = first;
	nextpool[current] = 0;
        
       /* increment cache indices*/
	if(phapcache1[f1-1].first == 0) {  
	  phapcache1[f1-1].first = current;
	  phapcache1[f1-1].last = current;
	}
	else {
	  nextpool[phapcache1[f1-1].last] = current;
	  phapcache1[f1-1].last = current;
	}
	current++;
      }
 
      /*similar to f1, if f2 is a possible haplotype for
        a child, cache the ind and first values for it*/
      if(flag[f2-1] != 0) {
	indpool[current] = ind;
	invpool[current] = first;
	nextpool[current] = 0;
	if(phapcache1[f2-1].first == 0) {
	  phapcache1[f2-1].first = current;
	  phapcache1[f2-1].last = current;
	}
	else {
	  nextpool[phapcache1[f2-1].last] = current;
	  phapcache1[f2-1].last = current;
	}
	current++;
      }
    }
  }

  for(i = 0; i < mgeno; i++) {
    qsumcache[0] += newwith3[i];
  }

 /*In this section of the code we update the probabilities for
  the parent p based on the probabilities for the child*/  
    newsegprob = femaletheta->segprob; /*use female probabilities to work with p*/
 /*Iterate over all joint genotypes of the child*/
  newflagr=thischild[0]->sparseflag;
  if(malechild[0]) {
  for(here = 0; here < mgeno; here++) {
    if(newflagr[here] == 0) continue;
    tempr = qsumcache[0] * newwithr[here];

    for(i=phapcache1[here].first; i != 0; i = nextpool[i]) {
      first = invpool[i];
      gene[first] += tempr*newsegprob[indpool[i]];
    }

  }
  }
  else {
  for(here = 0; here < fgeno; here++) {
    if(newflagr[here] == 0) continue;
    tempr = newwithr[here];
    c1 = invgenenum1[here];
    c2 = invgenenum2[here];

  /*now find all the values of first that as genotypes of the
 father could have passed on the other haplotype c1; multiply
  by their probabilties, which are stored in phapcache, and
  add to the entry gene[first], which will be used to form
 the new genarray for p*/

    temp1 = tempr * newwith3[c2 - 1];
    for(i=phapcache1[c1-1].first; i != 0; i = nextpool[i]) {
      first = invpool[i];
      gene[first] += temp1*newsegprob[indpool[i]];
    }

  /*if c1 is distinct from c2 we need to do the same computation reversing
    the roles of the two haplotypes*/
    if(c1 != c2) {
      temp1 = tempr * newwith3[c1 - 1];
      for(i=phapcache1[c2-1].first; i != 0; i = nextpool[i]) {
        first = invpool[i];
	gene[first] += temp1*newsegprob[indpool[i]];
      }
    }
  }
  }

/*set up new genarray for p; it is gene scaled by segscale*/      
  for (first = 0; first < fgeno; first++) {
    if(newflag2[first] == 0.0) continue;
    newwith2[first] *= segscale*gene[first];
  }

}


else  {  /*nchild is bigger than 1*/


/*find nonzero entries in q's genarray and make a list of them
  stored in nonzgens; just get one per isozygote class*/
  nonzcount = 0;

  /*iterate over genotypes for q*/
  for(i = 0; i < mgeno; i++)
    if(newflag3[i] != 0)
      nonzgens[nonzcount++] = i; /*store index of nonzero value*/

  newsegprob2 = maletheta->segprob; /*get male probabilties for q*/

  /*iterate over genotypes for p*/
  for (first = 0; first < fgeno; first += step1) {
      /*number of distinct probabilties needed for first's isoz. class*/
    step1 = probend[first] - probstart[first] + 1;
    skip = 1;

    /*work only on those isozygotes that are possible*/
    for(i = first; i < first+step1; i++)
      if(newflag2[i] != 0) {
	skip = 0;
	break; /*go to next isozygote in class*/
      }
    if(skip) continue;

    /*initialize update multiple for each isozygote in class*/
    for(i = 0; i < step1; i++) segval[i] = 0.0;

    /*iterate over the genotypes representing isozygote classes that
     q may have*/
    for (nonzindex = 0; nonzindex < nonzcount; nonzindex++) {
      second = nonzgens[nonzindex];
      /*check if this first, second pair yield a nonzero value
        among children's probabilities*/
      if(lsegsexfun2(second,first) == 0) continue;

      /*number of distinct probabilties needed for second's isoz. class*/

  
   /*call segsexsum2 to compute the part of the conditional probability
     update that is common to all the members of the combined isozygote
     class defined by first and second */
      segsexsum2(second,first,step1);

      /*now specialize update for each member of first's class*/
      for(i = 0; i < step1; i++) {
	if(newflag2[first+i] == 0) {
	  continue; /*skip if this isozygote is not possible*/
	}
        /*further specialize update for each member of second's isozygote
          class*/
        /*skip if this isozygote not possible*/
	  if(newflag3[second] == 0) continue; 

          /*get offset into newsegprob2; the offset depends on
            the isozygote class size and index of each parent
            note that fisozygoteindex gets incremented by the size
            of the joint class each time, so it represents a
            sum of all the numbers of distinct probabilities
            needed for all joint iso. classes  considered before the
            current p isozygote class*/
	  tempprob = newsegprob2 + probstart[first+i] - 1;

        /*combine for all children*/
        /*due to the arrangement in segsexsum2
          all the probability contributions
          for a given child are contiguous in the newwithr array.
          the number of contributions is fslength which is the number of
          probabilities needed for the joint isozygote class of the parents.
          We get the contribution of the first child (index 0) and
          then loop over the rest*/
	  val = tempprob[0] * tempseg[0];
         /*iterate over joint isozygote class*/
	  for(jointisoindex = 1; jointisoindex < step1; jointisoindex++)
	    val += tempprob[jointisoindex] * tempseg[jointisoindex];
	  newwithr = tempseg + step1; /*move the base of newwithr*/

         /*iterate over all the remaining children; notice that the
           loop index k is not used inside the loop. this is because the
           last loop statement which moves the offset of newwithr
           has the effect of incrementing the child*/
	  for(k = 1; k < nchild; k++) {
	    temp1 = tempprob[0] * newwithr[0];
	    for(jointisoindex = 1; jointisoindex < step1; jointisoindex++)
	      temp1 += tempprob[jointisoindex] * newwithr[jointisoindex];
	    val *= temp1;
	    newwithr += step1;
	  }
         /*update segval entry for this isozygote of first*/
	  segval[i] += newwith3[second] * val;
      }
    }
    /*update p's genarray for each isozygote in this class*/
    for(i = 0; i < step1; i++) newwith2[first+i] *= segval[i] * segscale;
  }
}

/* If any of the nonzero entries in p's genarray became 0,
   we want to set them to zero to avoid computations on subsequent
   calls*/

  for(i = 0; i < fgeno; i++)
    if((newflag2[i]==1) && (newwith2[i] == 0.0)) newflag2[i] = 0;
}

/*Added by A. A. Schaffer to handle loops */
  if (loopfirstgen && (!(*LINK->p)->loopdepend)) {
    depend = false;
    for(i=0; i<nchild; i++)
      if (childarray[i]->loopdepend) {
        depend = true;
        break;
      }
    depend = (depend || ((*LINK->q)->loopdepend));
    if (depend) {
      (*LINK->p)->loopdepend = depend;
      (*LINK->p)->loopneeded = false;
    }
  }
  if ((*LINK->p)->loopdepend) {
    if (!((*LINK->q)->loopdepend))
      (*LINK->q)->loopneeded = true;
    for(i=0; i< nchild; i++)
      if (!(childarray[i]->loopdepend))
        childarray[i]->loopneeded = true;
  }

  cleanup(LINK->q, LINK->LINK);
  exitseg(LINK);
}  /*segsexup*/

/* segsexdown updates the genotype probabilities for a child  based
on the probabilities of parents and other children. The parameter LINK includes
all the genetic information about p, p's spouse q, and their
children */

Void segsexdown(LINK)
struct LOC_seg *LINK;
{
  int  FORLIM; /*loop bound*/
  int nonzindex, nonzcount;    /*loop index and counter for nonzero values*/
  unsigned cgeno;
  unsigned int  step2; /*size of isozygote class*/
  double valtemp; /*intermediate values in probability updates*/
  double val, temp1; /*temporaries to store intermediate values*/
  unsigned int f1, f2; /*four haplotypes from parents*/
  unsigned int here, i, j, first, second; /*genotype indices*/
  int jointisoindex; /*indices to work within isozygote classes*/
  unsigned currentindex; /*index to update genarray within isozygote class*/
  unsigned k; /*index to loop over children*/
  thisarray *WITH2; /*stores genetic information about p*/
  thisarray *WITH3; /*stores gnetic information about q*/
  double *newwith2, *newwith3, *newwithr, *newwithc; /*store genarrays for
						       p,q, and children*/
  unsigned int c1, c2; /*haplotypes*/
  unsigned char *newflag2, *newflag3, *newflagr; /*store sparsity patterns for
					p and q genarrays*/
  unsigned char male;


  /* newsegprob and newsegprob2 are  used to hold segprob
    arrays, which contain the probabilities for various patterns of
    recombination events*/
  double *newsegprob, *newsegprob2; 
  double *tempprob; /*temporary holder for probability array*/
  boolean depend; /*used to handle loops*/

  initseg(LINK); /*get data about this p,q,children triple*/

  /*get sparsity patterns for p, q, and child genarrays*/
  newflag2 = (*LINK->p)->gen->sparseflag; 
  newflag3 = (*LINK->q)->gen->sparseflag;
  newflagr = (*LINK->r)->gen->sparseflag;

  /*In segsexdown p is always male and q is always female;
    This is not true for segsexup*/
  WITH2 = (*LINK->p)->gen; /*get genetic data for p*/
  WITH3 = (*LINK->q)->gen; /*get genetic data for q*/
  newwith2 = WITH2->genarray; /*get genarray for p*/
  newwith3 = WITH3->genarray; /*get genarray for q*/
  newwithr = (*LINK->r)->gen->genarray; /*get genarray for first child*/
  male = (*LINK->r)->male;
  if(male) cgeno = mgeno;
  else cgeno = fgeno;

/*The case of 1 child (nchild==0) is handled specially because the
subcomputations are much simpler. In particular we do not need to multiply the
probabilities among the different children. In a typical pedigree
many pairs will have only one child about which information is known
(in essence this is the child that keeps the pedigree connected), so
this is an important and common special case. */

if(nchild == 0) {
/*initialize cache data structures*/
  for(i = 0;i < maxhaplo;i++) {
    flag[i] = 0;
    psumcache[i] = 0.0;
    qsumcache[i] = 0.0;
  }

/*initialize gene array and set up flag array*/
    for(i = 0; i < cgeno; i++) {
      gene[i] = 0.0;
      flag[invgenenum1[i]-1] = 1;
      flag[invgenenum2[i]-1] = 1;
    }

  newsegprob = femaletheta->segprob; /*use female probabilities to work with q*/
 /*iterate over all joint genotypes*/
  for (first = 0; first < fgeno; first++ ) {
    if(newflag3[first] == 0) continue; /*check if joint genotype is possible*/
    valtemp = newwith3[first]; /*get cond. prob. that p has genotype first*/
    FORLIM = fence[first]; /*find bounds for the isozygote class of first*/

 /*iterate over all the recombined genotypes of the isozygote class*/
    for (i = base[first]; i < FORLIM; i++) {
      f1 = haps1[i];
      f2 = haps2[i];

      if((flag[f1-1] !=0) || (flag[f2-1] !=0)) {

/*condition probability of first as a genotype multiplied by probability
  of passing on f1 (alternatively f2) as a haplotype*/
        /*get probability offset for i*/
        val = valtemp*newsegprob[hind[i]];
     /*store in qsumcache*/
        if(flag[f1-1] != 0) {
	  qsumcache[f1 - 1] += val;
        }
        if(flag[f2-1] != 0) {
	  qsumcache[f2 - 1] += val;
        }
      }
    }
  }


 /*In this section of the code we update the probabilities for
  the child based on the probabilities for the parents*/  

 /*Iterate over all joint genotypes of the child*/

  if(male) {
    double temp1 = 0.0;
    for(here = 0; here < mgeno; here++) {
      temp1 += newwith2[here];
    }
    for(here = 0; here < mgeno; here++) {
      if(newflagr[here] == 0) continue;
  
      gene[here] =  temp1 * qsumcache[here];
  
    }
  }
  else {
    for(here = 0; here < fgeno; here++) {
      if(newflagr[here] == 0) continue;
      c1 = invgenenum1[here];
      c2 = invgenenum2[here];
  
     /*probability of child getting genotype here as a
       result of p passing on c1 and q passing on c2 is
       summed to gene[here] */
  
      gene[here] += newwith2[c1-1] * qsumcache[c2-1];
  
    /*if c1 is distinct from c2 we need to do the same computation reversing
      the roles of the two haplotypes*/
      if(c1 != c2) {
        gene[here] += newwith2[c2-1] * qsumcache[c1-1];
      }
    }
  }

/*set up new genarray for r; it is gene scaled by segscale*/      
    for (first = 0; first < cgeno; first++) {
      if(newflagr[first] == 0) continue;
      if(gene[first] == 0.0) newflagr[first] = 0;
      newwithr[first] *= segscale*gene[first];
    }

}


else  {  /*nchild is bigger than 0*/

 /*initialize genarray entries for child to 0*/
    for(i = 0; i < cgeno; i++) {
      gene[i] = 0.0;
    }

/*find nonzero entries in q's genarray and make a list of them
  stored in nonzgens; just get one per isozygote class*/
  nonzcount = 0;

  /*iterate over genotypes for q*/
  for(i = 0; i < fgeno; i += step2) {
    /*number of distinct probabilties needed for i's isoz. class*/
    step2 = probend[i] - probstart[i] + 1;
    for(j = i; j < i+step2; j++)
      if(newflag3[j] != 0) {
        nonzgens[nonzcount++] = i; /*store index of nonzero value*/
	break;                        /*go to next isozygote class*/
      }
  }

 newsegprob2 = maletheta->segprob; /*get p probabilties*/

  /*iterate over genotypes for p*/
  for (first = 0; first < mgeno; first++) {

    if(newflag2[first] == 0) continue;

    /*iterate over the genotypes representing isozygote classes that
     q may have*/
    for (nonzindex = 0; nonzindex < nonzcount; nonzindex++) {
      second = nonzgens[nonzindex];
      /*check if this first, second pair yield a nonzero value
        among children's probabilities*/
      if(lsegsexfun2(first,second) == 0) continue;

      /*number of distinct probabilties needed for second's isoz. class*/
      step2 = probend[second] - probstart[second] + 1;
  
   /*call segsexsumdown2 to compute the part of the conditional probability
     update that is common to all the members of the combined isozygote
     class defined by first and second */
      segsexsumdown2(first,second,step2,male);

      for(jointisoindex = 0; jointisoindex < step2; jointisoindex++)
	segval[jointisoindex] = 0.0;

      /*now specialize update for each member of first's class*/
	if(newflag2[first] == 0) continue;

        /*further specialize update for each member of second's isozygote
          class*/
	for(j = 0; j < step2; j++) {
        /*skip if this isozygote not possible*/
	  if(newflag3[second+j] == 0) {
	    continue; 
	  }

	  tempprob = newsegprob2 + probstart[second+j] - 1;

        /*combine for all children*/
        /*due to the arrangement in segsexsumdown2
          all the probability contributions
          for a given child are contiguous in the newwithr array.
          the number of contributions is fslength which is the number of
          probabilties needed for the joint isozygote class for the parents.
          We get the contribution of the first child (index 0) and
          then loop over the rest*/
	  val = tempprob[0] * tempseg[0];
         /*iterate over joint isozygote class*/
	  for(jointisoindex = 1; jointisoindex < step2; jointisoindex++)
	    val += tempprob[jointisoindex] * tempseg[jointisoindex];
	  newwithc = tempseg + step2; /*move the base of newwithc*/

         /*iterate over all the remaining children; notice that the
           loop index k is not used inside the loop. this is because the
           last loop statement which moves the offset of newwithc
           has the effect of incrementing the child*/
	  for(k = 1; k < nchild; k++) {
	    temp1 = tempprob[0] * newwithc[0];
	    for(jointisoindex = 1; jointisoindex < step2; jointisoindex++)
	      temp1 += tempprob[jointisoindex] * newwithc[jointisoindex];
	    val *= temp1;
	    newwithc += step2;
	  }
          /*probability of this combination of parent genotypes*/
	  val *= newwith2[first] * newwith3[second+j];
	  for(jointisoindex = 0; jointisoindex < step2; jointisoindex++) {
          /*probability of this recombination pattern (based on other children)*/
	    segval[jointisoindex] += tempprob[jointisoindex] * val;
	  }
	}
      /*update the probabilities of two joint genotypes the
       child might get; each different choice of recombination
       pattern will lead to a different set of two genotypes*/
      currentindex = 0;
      for(jointisoindex = 0; jointisoindex < step2; jointisoindex++) {
	temp1 = segval[jointisoindex];
	gene[segindex[currentindex++]-1] += temp1;
	gene[segindex[currentindex++]-1] += temp1;
      }
    }
  }
  /*finally update child's real genarray by coppy gene multiplied
    by scale factor segscale*/
  for(i = 0; i < cgeno; i++) {
    if(newflagr[i] == 0) continue;
    if(gene[i] == 0.0) newflagr[i] = 0; /*if probability changes from nonzero
					  to 0.0 change flag to 0*/
    newwithr[i] *= segscale * gene[i];
  }
}

/*Added by A. A. Schaffer to handle loops */
  if (loopfirstgen && (!((*LINK->r)->loopdepend))) {
    depend = false;
    for(i=0; i<nchild; i++)
      if (childarray[i]->loopdepend) {
        depend = true;
        break;
      }
    depend = (depend || ((*LINK->p)->loopdepend) || ((*LINK->q)->loopdepend));
    if (depend) {
     (*LINK->r)->loopdepend = depend;
     (*LINK->r)->loopneeded = false;
    }
  }
  if ((*LINK->r)->loopdepend) {
    if (!((*LINK->p)->loopdepend))
      (*LINK->p)->loopneeded = true;
    if (!((*LINK->q)->loopdepend))
      (*LINK->q)->loopneeded = true;
    for(i=0; i< nchild; i++)
      if (!(childarray[i]->loopdepend))
        childarray[i]->loopneeded = true;
  }
  cleanup(LINK->q, LINK->LINK);
  exitseg(LINK);
}  /*segsexdown*/






