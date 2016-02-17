/* This file contains a modified version of the MLINK program */
/* The modifications are described in the papers: */
/* R. W. Cottingham, Jr., R. M. Idury, and A. A. Schaffer */
/* Faster Sequential Genetic Linkage Computations */
/* American Journal of Human Genetics, 53(1993), pp. 252-263 */
/* and A. A. Schaffer, S. K. Gupta, K. Shriram, and R. W. Cottingham, Jr. */
/* Avoiding Recomputation in Linkage Analysis, Human Heredity */
/* 44(1994), pp. 225-237. */

#include "commondefs.h"
#include "gemdefs.h"
#include "mldefs.h"
#if !defined(DOS)
#include "checkpointdefs.h"
#endif
#ifndef LESSMEMORY
#include "moddefs.h"
#endif


#if PARALLEL  /* cgh */
#include "compar.h"           /* parallel support code */
#endif  /* defined(PARALLEL) -- cgh */


extern void checkzero();
#if !defined(KNR_PROTO)
extern boolean zerotest(double);
#else
extern boolean zerotest();
#endif  /* !defined(KNR_PROTO) */


#if PARALLEL  /* cgh */

/* Calculate the number of calls to iterpeds() by simulating the
   nested loop in main().  We utilize temporary copies of all global
   variables to ensure nothing external is affected. */
void simIpedLoop(status)
loopStatus status;  
{
  int i, j;
  int currentIter;        /* current number of simulated iteration */
  long datafilePos;       /* current position in datafile */
  
  /* initialize temporaries */
  thetarray tMaletheta;
  int tWhich = which;
  double tInc = inc;
  double tFinish = finish;
  for (i = 0; i < mlocus - 1; i++)
    tMaletheta[i] = maletheta->theta[i];

  /* remember where we are in datafile.dat */
  datafilePos = ftell(datafile);

  if (status == countIpeds)
    /* initialize zero theta counter */
    numZeroThetas = 0;
  
  /* initialize the loop counter */
  if (score && (!risk)) {
    /* we will do an initial iteration of iterpeds() with 0.5 theta */ 
    currentIter = 1;
    if (status == computeThetas) {
      /* we have an initial call to iterpeds() here with theta = 0.5 */
      for (i = 0; i < mlocus - 1; i++)
	gmaletheta[0][i] = tMaletheta[i];
      gmaletheta[0][tWhich - 1] = 0.5;
      /* see if this thetarray has a 0 */
      gZeroTheta[0] = false;
      for (i = 1; i < mlocus - 1; i++)
	if (zerotest(tMaletheta[i]))
	  gZeroTheta[0] = true;
    } else {
      /* we are counting iterations -- but we also want to count
	 the number of zero thetas */
      for (i = 1; i < mlocus - 1; i++)
	if (zerotest(tMaletheta[i])) {
	  numZeroThetas++;
	  break;
	}
    }
  } else {
    /* no initial call, so the first iteration is 0 */
    currentIter = 0;
  }

  /* here we simulate the nested loop in main that calls iterpeds() */
  do {
    while (tMaletheta[tWhich - 1] <= tFinish) {
      /* calculate thetas, if appropriate */
      if (status == computeThetas) {
	gZeroTheta[currentIter] = false;
	for (i = 0; i < mlocus - 1; i++) {
	  gmaletheta[currentIter][i] = tMaletheta[i];
	  if (zerotest(tMaletheta[i]))
	    /* if this thetarray has a 0, record it */
	    gZeroTheta[currentIter] = true;
	}
      } else {
	for (i = 0; i < mlocus - 1; i++) {
	  if (zerotest(tMaletheta[i])) {
	    /* count the number of zero thetas */
	    numZeroThetas++;
	    break;
	  }
	}
      }
      
      /* this is where iterpeds() is called */
      currentIter++;

      tMaletheta[tWhich - 1] += tInc;
      if (zerotest(tMaletheta[tWhich - 1]))
	tMaletheta[tWhich - 1] = 0.0;
    }
    if (!P_eof(datafile)) {
      fscanf(datafile, "%lf%lf%lf%*[^\n]",
	     &tMaletheta[tWhich - 1], &tInc, &tFinish);
      getc(datafile);
    } else
      tWhich = 0;
    /* this increment is used to avoid round-off errors when checking
       for the finishing value */ 
    tFinish += ROUNDOFF_INC;
  } while (tWhich != 0);

  if (datafile == NULL)
    Tmk_errexit("simIpedLoop() error - datafile is NULL\n");
  /* reset datafile to its previous position */
  if (fseek(datafile,datafilePos,0) != 0)
      Tmk_errexit("simIpedLoop() error - fseek failed\n");
  /* we probably hit feof, so clear it */
  clearerr(datafile);

  /* if it was our job to count the calls to iterpeds(), set the
     global numIpeds variable */
  if (status == countIpeds) {
    numIpeds = currentIter;
  } else {
    /* calculate the number of nonzero theta vectors */
    for (i = 0, numNonzeroThetas = 0; i < numIpeds; i++)
      if (gZeroTheta[i] == false)
	/* increment count of nonzero thetas, and save ordering
	   of nozero thetas */
	absoluteThetaIndex[numNonzeroThetas++] = i;
    /* save ordering of theta vectors that contain zeros */
    for (i = 0, j = numNonzeroThetas; i < numIpeds; i++)
      if (gZeroTheta[i] == true)
	absoluteThetaIndex[j++] = i; 
  }
  return;
}

#else  /* if PARALLEL -- cgh */


/*The following routine iterates over the different pedigrees and handles 
  output. Carol Haynes suggested adding the printing of lod scores for
  each family*/
void iterpeds()
{
  int i;
  int thisped;
  int FORLIM;

  static double eachlod[maxped]; /*C. Haynes */
  int II=0;   /*C. Haynes*/

  tlike = 0.0;
  alike = 0.0;
  for (i = 1; i <= totperson; i++) {
    person[i]->done = false;
    person[i]->firstpass = true;
  }
  thisc = minint;
  recombination();
  checkzero();
  
  if (normalRun == checkpointStatus) {
    for (i = 1; i <= 35; i++)
      putc('-', outfile);
    putc('\n', outfile);
    for (i = 1; i <= 35; i++)
      putc('-', outfile);
    putc('\n', outfile);
    if (sexdif)
      fprintf(outfile, "MALE THETAS   ");
    else
      fprintf(outfile, "THETAS ");
  }
  FORLIM = mlocus - 2;
  if (normalRun == checkpointStatus) {
    for (i = 0; i <= FORLIM; i++)
      fprintf(outfile, "%6.3f", maletheta->theta[i]);
    if (interfer)
      fprintf(outfile, "%6.3f", maletheta->theta[mlocus - 1]);
    putc('\n', outfile);
    if (dostream) {
      FORLIM = mlocus - 2;
    for (i = 0; i <= FORLIM; i++)
      fprintf(stream, "%6.3f\n", maletheta->theta[i]);
    }
    if (interfer && dostream)
      fprintf(stream, "%6.3f\n", maletheta->theta[mlocus - 1]);
    if (sexdif) {
      fprintf(outfile, "FEMALE THETAS ");
      FORLIM = mlocus - 2;
      for (i = 0; i <= FORLIM; i++)
	fprintf(outfile, "%6.3f", femaletheta->theta[i]);
      if (interfer)
	fprintf(outfile, "%6.3f", femaletheta->theta[mlocus - 1]);
      putc('\n', outfile);
      if (dostream) {
	FORLIM = mlocus - 2;
	for (i = 0; i <= FORLIM; i++)
	  fprintf(stream, "%6.3f\n", femaletheta->theta[i]);
      }
      if (interfer && dostream)
	fprintf(stream, "%6.3f\n", femaletheta->theta[mlocus - 1]);
    }
    for (i = 1; i <= 35; i++)
      putc('-', outfile);
    fprintf(outfile, "\nPEDIGREE |  LN LIKE  | LOG 10 LIKE\n");
    for (i = 1; i <= 35; i++)
      putc('-', outfile);
    putc('\n', outfile);

    for (i = 1; i <= 35; i++)
      putchar('-');
    putchar('\n');
    for (i = 1; i <= 35; i++)
      putchar('-');
    putchar('\n');
    if (sexdif)
      printf("MALE THETAS   ");
    else
      printf("THETAS ");
    FORLIM = mlocus - 2;
    for (i = 0; i <= FORLIM; i++)
      printf("%6.3f", maletheta->theta[i]);
    if (interfer)
      printf("%6.3f", maletheta->theta[mlocus - 1]);
    putchar('\n');
    if (sexdif) {
      printf("FEMALE THETAS ");
      FORLIM = mlocus - 2;
      for (i = 0; i <= FORLIM; i++)
	printf("%6.3f", femaletheta->theta[i]);
      if (interfer)
	printf("%6.3f", femaletheta->theta[mlocus - 1]);
      putchar('\n');
    }
    for (i = 1; i <= 35; i++)
      putchar('-');
    printf("\nPEDIGREE |  LN LIKE  | LOG 10 LIKE\n");
    for (i = 1; i <= 35; i++)
      putchar('-');
    putchar('\n');
  }
#if LOOPSPEED
  open_loop_file();
#endif
  for (thisped = 0; thisped < nuped; thisped++) {
#if LOOPSPEED
    read_loop_file(thisped + 1);
#endif
    likelihood((thisped + 1), proband[thisped]);
    if (byfamily && (normalRun == checkpointStatus))
      fprintf(outfile, "%9d %12.6f ", proband[thisped]->ped, like);
    if (dostream && (normalRun == checkpointStatus))
      fprintf(stream, "%12d % .5e", proband[thisped]->ped, like);
    if (normalRun == checkpointStatus)
      printf("%9d %12.6f ", proband[thisped]->ped, like);
    alike += like;
    like /= log10_;

    /* C. Haynes*/   
    if (maletheta->theta[which - 1] == 0.5)
      eachlod[II] = like;

    if (byfamily && (normalRun == checkpointStatus))
      if (lodbyfamily)
	fprintf(outfile, "%12.6f LOD= %12.6f\n", like, 
	      ((like-eachlod[II]) > BIGNEGATIVE) ? like-eachlod[II] : -999.999999);
      else
	fprintf(outfile, "%12.6f\n", like);
    if (dostream && (normalRun == checkpointStatus))
      fprintf(stream, " %.5e\n", like); /* ERouillard -- "large" likelihoods and LRP */
    if (normalRun == checkpointStatus)
      if (lodbyfamily)
	printf("%12.6f LOD= %12.6f\n", like,
	      ((like - eachlod[II]) > BIGNEGATIVE) ? like-eachlod[II] : -999.999999);
      else
	printf("%12.6f\n", like);
    tlike += like;
    II++;  /*C. Haynes*/
  }
#if LOOPSPEED
  close_loop_file();
#endif
  if (normalRun == checkpointStatus) {
    for (i = 1; i <= 35; i++)
      putc('-', outfile);
    fprintf(outfile, "\nTOTALS    %12.6f %12.6f\n", alike, tlike);
    for (i = 1; i <= 35; i++)
      putchar('-');
    printf("\nTOTALS    %12.6f %12.6f\n", alike, tlike);
    if (dostream)
      fprintf(stream, "% .5e  % .5e\n", alike, tlike);
  }
  alike = -2 * alike;
  if (normalRun == checkpointStatus) {
    fprintf(outfile, "-2 LN(LIKE) = % .5e", alike);
    printf("-2 LN(LIKE) = % .5e", alike);
  }
  if (score && (!risk)) {
    if (mlocus == 2) {
      if (maletheta->theta[which - 1] == 0.5)
	scorevalue = tlike;
      alike = tlike - scorevalue;
      if (normalRun == checkpointStatus) {
	fprintf(outfile, " LOD SCORE = %12.6f", alike);
	printf(" LOD SCORE = %12.6f", alike);
      }
    } else {
      if (maletheta->theta[which - 1] == 0.5)
	scorevalue = alike;
      alike = scorevalue - alike;
      if (normalRun == checkpointStatus) {
	fprintf(outfile, " LOG LIKE DIFFERENCE = %12.6f", alike);
	printf(" LOG LIKE DIFFERENCE = %12.6f", alike);
      }
    }
  }
  if (score && (!risk) && dostream)
    if (normalRun == checkpointStatus) 
      fprintf(stream, "% .5e\n", tlike - scorevalue);
  if (normalRun == checkpointStatus) {
    putc('\n', outfile);
    putchar('\n');
  }
  if (firsttime) {
    if (thisc < maxcensor)
      printf("Maxcensor can be reduced to %12d\n", thisc);
    else {
      if (thisc > maxcensor)
	printf("you may gain efficiency by increasing maxcensor\n");
    }
  }
  firsttime = false;
}
#endif  /* if PARALLEL -- cgh */


void preIpedLoop()
{
  if (dostream && (checkpointedRun != checkpointStatus)) {
    fprintf(stream, "@\n");
    fprintf(stream, "MLINK\n");
    fprintf(stream, "%5d%5ld\n", mlocus, which);
    if (score && (!risk))
      fprintf(stream, " 1\n");
    else
      fprintf(stream, " 0\n");
    if (interfer)
      fprintf(stream, " 1\n");
    else
      fprintf(stream, " 0\n");
    if (sexdif && readfemale)
      fprintf(stream, "2 ");
    else {
      if (sexdif)
	fprintf(stream, "1\n");
      else
	fprintf(stream, "0\n");
    }
    fprintf(stream, "%12d\n", nuped);
    for (i = 1; i <= mlocus; i++) {
      j = 0;
      do {
	j++;
      } while (j != order[i - 1]);
      fprintf(stream, "%12ld\n", j);
    }
    putc('\n', stream);

    /* This call to recombination() is needed to make sure that the
       proper theta values are printed to stream.  recombination()
       is called again before liklihood() is ever called from within
       iterpeds() (in sequential MLINK), or multiPedLike() (in
       parallel MLINK).  Therefore, if dostream is false, this
       call can be safely omitted.
       In parallel MLINK, we don't want the barrier inside of
       recombination() to get called, and we don't want
       g{male,female}segprob assigned to.  This code in
       recombination() is only executed when the boolean flag infun
       is true, so we set it to false before the call. -- cgh */
#if PARALLEL  /* cgh */
    infun = false;
#endif  /* if PARALLEL -- cgh */    
    recombination();
    
    for (i = 1; i < mlocus; i++)
      fprintf(stream, "% .5e\n", maletheta->theta[i - 1]);
    if (interfer)
      fprintf(stream, "%6.3f\n", maletheta->theta[mlocus - 1]);
    if (sexdif) {
      for (i = 1; i < mlocus; i++)
	fprintf(stream, "% .5e\n", femaletheta->theta[i - 1]);
    }
    if (interfer && sexdif)
      fprintf(stream, "%6.3f\n", femaletheta->theta[mlocus - 1]);
  }

  fprintf(outfile, "LINKAGE (V%3.1f) WITH%3d-POINT", fVersion, mlocus);
  if (sexlink)
    fprintf(outfile, " SEXLINKED DATA\n");
  else
    fprintf(outfile, " AUTOSOMAL DATA\n");
  fprintf(outfile, " ORDER OF LOCI: ");
  for (i = 1; i <= mlocus; i++) {
    thissystem = 1;
    while (order[thissystem - 1] != i)
      thissystem++;
    fprintf(outfile, "%3ld", thissystem);
  }
  putc('\n', outfile);
  scorevalue = 0.0;
  
#if !PARALLEL  /* cgh */
  if (score && (!risk)) {
    holdtheta = maletheta->theta[which - 1];
    maletheta->theta[which - 1] = 0.5;
    if (checkpointedRun == checkpointStatus) {
      printf("\nTo get LODSCOREs correct the program will redo the");
      printf("\nlikelihood evaluation with moving locus unlinked");
      printf("\nOutput will be suppressed\n\n");
      fflush(stdout);
    }
    iterpeds();
    if (normalRun == checkpointStatus) {
      ensureWrite(outfile);
      ensureWrite(stream);
    }
    maletheta->theta[which - 1] = holdtheta;
  }
#endif  /* if !PARALLEL -- cgh */
}


#if !PARALLEL
void ipedLoop()
{
#if !defined(DOS)  
  iterpeds_counter = 0;  /*set up for checkpointing*/
#endif  /* !defined(DOS) */
  do {
    while (maletheta->theta[which - 1] <= finish) {
#if !defined(DOS)
      if ((0 == checkpoint_counter) &&
	  (iterped_call_before == checkpoint_place))
	checkpointStatus = normalRun;
      if (normalRun == checkpointStatus) {
	performCheckpoint(functionLocation, iterpeds_counter, iterped_call_before);
	copyFile("outfile.dat", MainMLINKOutfileDat);
        if (dostream)
          copyFile("stream.dat", MainMLINKStreamDat);
      }
      if (( 0 == checkpoint_counter) &&
	  (iterped_call_before == checkpoint_place)) {
#endif  /* if !defined(DOS) */
	iterpeds();
	ensureWrite(outfile);
	ensureWrite(stream);
#if !defined(DOS)
	performCheckpoint(functionLocation, iterpeds_counter, iterped_call_after);
      }
      if ((checkpoint_counter == 0) && (iterped_call_after == checkpoint_place)){
	checkpoint_place = iterped_call_before;
	checkpointStatus = normalRun;
      }
      if (normalRun == checkpointStatus) {
	copyFile("outfile.dat", MainMLINKOutfileDat);
	if (dostream)
	  copyFile("stream.dat", MainMLINKStreamDat);
      }
#endif  /* if !defined(DOS) */
      maletheta->theta[which - 1] += inc;
      if (zerotest(maletheta->theta[which - 1]))
	maletheta->theta[which - 1] = 0.0;
#if !defined(DOS)      
      if (checkpoint_counter > 0)
	checkpoint_counter--;
      iterpeds_counter++;
#endif  /* if !defined(DOS) */
    }
    if (!P_eof(datafile)) {
      fscanf(datafile, "%lf%lf%lf%*[^\n]", &maletheta->theta[which - 1], &inc,
	     &finish);
      getc(datafile);
    } else
      which = 0;
    /* cgh -- this increment is used to avoid round-off errors when
       checking for the finishing value */ 
    finish += ROUNDOFF_INC;
  } while (which != 0);
}
#endif  /* !PARALLEL */



/* End. */
