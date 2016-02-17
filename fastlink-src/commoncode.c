/* This file contains common routines abstracted from main() in all
   versions of the FASTLINK programs ILINK, LINKMAP, and MLINK.
   Sequential FASTLINK is an improved version of LINKAGE.  Improvements
   are described in described in: R. W. Cottingham, Jr., R. M. Idury, and
   A. A. Schaffer Faster Sequential Genetic Linkage Computations American
   Journal of Human Genetics, 53(1993), pp. 252--263 and A. A. Schaffer,
   S. K. Gupta, K. Shriram, and R. W. Cottingham, Jr., Avoiding
   Recomputation in Linkage Analysis Human Heredity 44(1994),
   pp. 225-237.  The parallel implementations of ILINK are described in:
   S. Dwarkadas, A. A. Schaffer, R. W. Cottingham Jr., A. L. Cox,
   P. Keleher, and W. Zwaenepoel, Parallelization of General Linkage
   Analysis Problems, Human Heredity 44(1994), pp. 127-141 and
   S. K. Gupta, A. A. Schaffer, A. L. Cox, S. Dwarkadas, and
   W. Zwaenepoel, Integerating Parallelization Strategies for Linkage
   Analysis, Computers and Biomedical Research, to appear.
   The code in this file was written by Chris Hyams. */
   
#include "commondefs.h"
#include "gemdefs.h"

#if !defined(LESSMEMORY)
#include "moddefs.h"
#endif  /* if !defined(LESSMEMORY) */

#if !defined(DOS)
#include "checkpointdefs.h"
#endif  /* if !defined(DOS) */

#if defined(MLINK)
#include "mldefs.h"
#endif  /* if defined(MLINK) */

#if defined(LINKMAP)
#include "lidefs.h"
#endif  /* if defined(LINKMAP) */

#if defined(ILINK)
#include "ildefs.h"
#endif  /* if defined(ILINK) */

#if PARALLEL
#include "compar.h"           /* parallel support code */
#endif  /* defined(PARALLEL) */

#if (defined(MLINK) || defined(LINKMAP))
extern void preIpedLoop();

#if PARALLEL

#if !defined(KNR_PROTO)
extern void simIpedLoop(loopStatus);
#else  /* !defined(KNR_PROTO) */
extern void simIpedLoop();
#endif   /* !defined(KNR_PROTO) */

#else  /* if PARALLEL */

extern void ipedLoop();

#endif  /* if PARALLEL */
#endif   /* (defined(MLINK) || defined(LINKMAP) */

#if defined(GMEM_DEBUG)
int totMem = 0;
#endif  /* defined(GMEM_DEBUG) */


double mapfunction(theta1, theta2)
double theta1, theta2;
{
  /*User defined function giving recombination between
    flanking markers as a function of recombination
    between adjacent markers*/
  return ((theta1 + theta2) / (1 + 4 * theta1 * theta2));
}


double getdist(theta)
double *theta;
{
  if (*theta < 0.5)
    return (log(1.0 - 2.0 * *theta) / -2.0);
  else
    return 10.0;
}


double invdist(dist)
double *dist;
{
#if defined(ILINK)
  if (*dist != 10.0)
#else
  if (*dist < 10.0)
#endif  /* defined(ILINK) */
    return ((1 - exp(-2 * *dist)) / 2.0);
  else
    return 0.5;
}


static void printVersion()
{
  printf("\nProgram %s version%6.2f (1-Feb-1991)\n\n", PROGRAM, fVersion);
  printf("\nFASTLINK ");
#if defined(LESSMEMORY)
  printf("(slow) ");
#endif  
  printf("version %s (30-Jun-1999)", fastversion);
#if PARALLEL
  printf(" (PARALLEL)");
#endif  
  printf("\n\n");
}


void initialize()
{
#if PRALLEL
  if (Tmk_proc_id == 0) {
#endif
  printVersion();
  printf("The program constants are set to the following maxima:\n");
#if defined(ILINK)
  printf("%6d loci in mapping problem\n", (int)maxlocus);
#else
  printf("%6d loci in mapping problem (maxlocus)\n", (int)maxlocus);
#endif  /* defined(ILINK) */
  printf("%6d alleles at a single locus (maxall)\n", (int)maxall);
  printf("%6ld maximum of censoring array (maxcensor)\n", maxcensor);
#if defined(ILINK)
  printf("%6d individuals in all pedigrees combined\n", (int)maxind);
#else
  printf("%6d individuals in all pedigrees combined (maxind)\n",
	 (int)maxind);
#endif  /* defined(ILINK) */  
  printf("%6d pedigrees (maxped)\n", (int)maxped);
#if !defined(ILINK)
  printf("%6d binary codes at a single locus (maxfact)\n", (int)maxfact);
#endif  /* !defined(ILINK) */  
  printf("%6d quantitative factor(s) at a single locus\n", (int)maxtrait);
  printf("%6d liability classes\n", (int)maxliab);
  printf("%6d binary codes at a single locus\n", (int)maxfact);
  printf("%8.2f base scaling factor for likelihood (scale)\n", scale);
  printf("%8.2f scale multiplier for each locus (scalemult)\n", scalemult);
  printf("%8.5f frequency for elimination of heterozygotes (minfreq)\n",
	 minfreq);
  if (minfreq != 0.0) {
    printf("IMPORTANT : RECOMPILE THIS PROGRAM WITH MINFREQ=0.0\n");
    printf("FOR THE ANALYSIS OF RECESSIVE TRAITS\n");
  }
  putchar('\n');
#if PRALLEL
}
#endif
}  /*initialize*/



void printInfo()
{
  printVersion();
  printf("\n%s has been compiled with the following options:\n\n",
	 PROGRAM);

#if PARALLEL
  printf("       PARALLEL computation (PARALLEL defined)\n");
#endif  
  
  printf("       CHECKPOINTING is ");
#if defined(DOS)
  printf("disabled (DOS defined)\n");
#else
  printf("enabled (DOS not defined)\n");
#endif  /* defined(DOS) */

#if defined(LESSMEMORY)  
  printf("       SLOW version (LESSMEMORY defined)\n");
#else
  printf("       FAST version (LESSMEMORY not defined)\n");  
#endif  /* defined(LESSMEMORY) */

#if PARALLEL
  printf("       PRECOMPUTE = ");
#if PRECOMPUTE
  printf("1\n");
#else
  printf("0\n");
#endif  
#endif  /* PARALLEL */
  
  printf("\nProgram constants are set to the following maxima:\n\n");
  printf("%6d maximum number of loci (maxlocus)\n", (int)maxlocus);
  printf("%6d maximum number of alleles at a single locus (maxall)\n",
	 (int)maxall);
  printf("%6d maximum number of individuals in a pedigree (maxind)\n",
	 (int)maxind);
  printf("%6d maximum number of loops (maxloop)\n", (int)maxloop);
  printf("%6d maximum number of children of a parent (maxchild)\n",
	 (int)maxchild);
  printf("\n");
  exit(0);
}



void openFiles()
{
  final = NULL;
  speedfile = NULL;
  stream = NULL;
  datafile = NULL;
  ipedfile = NULL;
  outfile = NULL;

  datafile = fopen("datafile.dat", "r");
  if (datafile == NULL)
    exit(FileNotFound);
  ipedfile = fopen("ipedfile.dat", "r");
  if (ipedfile == NULL)
    exit(FileNotFound);
#if ALLELE_SPEED
  speedfile = fopen("newspeedfile.dat", "r");
#else
  speedfile = fopen("speedfile.dat", "r");
#endif
  if (speedfile == NULL)
    exit(FileNotFound);

#if PARALLEL
#if defined(LESSMEMORY)
  if (PRECOMPUTE)
    Tmk_errexit("Cannot have LESSMEMORY defined and PRECOMPUTE set to 1\n");
#endif
  
/* If you change things in connexion with either final or stream, be
   sure to make the appropriate modifications to the checkpointing
   code near the beginning of outf().  --Shriram */
     if (Tmk_proc_id ==0) {
#endif /* PARALLEL */
#ifdef vms
    outfile = fopen("outfile.dat", "w", "ctx=rec","shr=get,put,upd" ) ;
#else  /* ifdef vms */
    outfile = fopen("outfile.dat", "w");
#endif  /* ifdef vms */
    if (outfile == NULL)
      exit(FileNotFound);
    if (dostream) {
#ifdef vms
      stream = fopen("stream.dat", "w", "ctx=rec","shr=get,put,upd" ) ;
#else  /* ifdef vms */
      stream = fopen("stream.dat", "w");
#endif  /* ifdef vms */
      if (stream == NULL)
	exit(FileNotFound);
    }
#if defined(ILINK)
#ifdef vms
    final = fopen("final.dat", "w", "ctx=rec","shr=get,put,upd" ) ;
#else  /* ifdef vms */
    final = fopen("final.dat", "w");      
#endif  /* ifdef vms */
    if (final == NULL)
      exit(FileNotFound);
#endif  /* defined(ILINK) */
#if PARALLEL
  }
#endif  /* if PARALLEL */
}


void ckptInit()
{
#if PARALLEL  /* cgh */
  if (Tmk_proc_id == 0) {
#endif /* if PARALLEL -- cgh */

#if !defined(DOS)
  /* Shriram: begin */
  /*Perform the script level checkpointing here*/
#ifdef vms  
  scriptCheckpoint = fopen ( ScriptCheckpointFilename ,
			     "r", "ctx=rec","shr=get,put,upd" ) ;
#else  /* ifdef vms */
  scriptCheckpoint = fopen ( ScriptCheckpointFilename , "r" );
#endif  /* ifdef vms */
  if ( NULL != scriptCheckpoint )
  {
    fscanf ( scriptCheckpoint , "%d %d" ,
	     & scriptRun , & scriptToSleep ) ;
    fclose ( scriptCheckpoint ) ;
    if ( 0 != scriptToSleep )
    {
#ifdef vms
      scriptCheckpoint = fopen ( ScriptCheckpointFilename ,
				 "w", "ctx=rec","shr=get,put,upd" ) ;
#else  /* ifdef vms */
      scriptCheckpoint = fopen ( ScriptCheckpointFilename , "w" ) ;
#endif  /* ifdef vms */
      scriptRun ++ ;  scriptToSleep -- ;
      fprintf ( scriptCheckpoint , "%d %d\n" ,
		scriptRun , scriptToSleep ) ;
      fclose ( scriptCheckpoint ) ;
      printf ( "Recovering to checkpoint: %d more run(s)\n" ,
               scriptToSleep + 1 ) ;
#if defined(ILINK)
      if ( NULL != final )
        fclose ( final ) ;
#endif  /* defined(ILINK) */
#if defined(MLINK) || defined(LINKMAP)
      if ( NULL != outfile )
        fclose ( outfile ) ;
#endif  /* defined(MLINK) || defined(LINKMAP) */
      if ( NULL != stream )
        fclose ( stream ) ;
      
      if ( 0 == scriptToSleep )
      {
#ifdef vms
        fileTester = fopen ( ScriptFinalOut ,
			     "r" , "ctx=rec","shr=get,put,upd" ) ;
#else  /* ifdef vms */
        fileTester = fopen ( ScriptFinalOut , "r" ) ;
#endif  /* ifdef vms */
        if ( NULL != fileTester )
        {
          fclose ( fileTester ) ;
	  copyFile ( ScriptFinalOut , "final.out" ) ;
        }
        if ( dostream )
        {
#ifdef vms
          fileTester = fopen ( ScriptStreamOut ,
			       "r", "ctx=rec","shr=get,put,upd" ) ;
#else  /* ifdef vms */
          fileTester = fopen ( ScriptStreamOut , "r" ) ;
#endif  /* ifdef vms */
          if ( NULL != fileTester )
          {
            fclose ( fileTester ) ;
	    copyFile ( ScriptStreamOut , "stream.out" ) ;
          }
        }
      }
      exit ( EXIT_FAILURE ) ;
    }
  }
#endif  /* if !defined(DOS) */
  
#if PARALLEL  /* cgh */
  }    /* if (Tmk_proc_id == 0) */
#endif /* if PARALLEL -- cgh */
}


#if !defined(DOS)
void ckptEnsureFileContents()
{
  /* Shriram: begin */
  
  /* We will use stat() to determine the size of the checkpointing
     file.  If stat() reports any of the following: the file cannot
     be accessed; there is a problem with the filename or with the
     buffer; the file does not exist; or the file is shorter than the
     length of the date/time stamp we put in it, then we assume that
     no checkpointing has taken place.
     If we don't ensure the file has contents now, we would be doing
     this later, and might need to return here, etc. */
  
  statBuffer.st_size = 0 ;
  statReturnCode = stat(CheckpointFilename, &statBuffer);
  if (statReturnCode != STAT_SUCCESS)
    checkpointStatus = normalRun ;
  else {
#ifdef vms
    checkpointDatafile =
      fopen(CheckpointFilename , "r", "ctx=rec","shr=get,put,upd");
#else  /* vms */
    checkpointDatafile = fopen ( CheckpointFilename , "r" ) ;
#endif  /* vms */
    checkpointStatus = checkpointedRun ;
    puts("NOTE: attempting to continue previous (unfinished) run");
    fgets(dateTimeStamp ,
	  DateTimeStampStringLength , checkpointDatafile ) ;
    printf ( "      from %s" , dateTimeStamp ) ;
    getCkptTuple ( ) ;
  }
  /* Shriram: end */
}
#endif  /* if !defined(DOS) */


void closeInputFiles()
{
#if (defined(LINKMAP) || defined(ILINK))
  if (datafile != NULL)
    fclose(datafile);
  datafile = NULL;
#endif  /* defined(LINKMAP) || defined(ILINK) */
  if (ipedfile != NULL)
    fclose(ipedfile);
  ipedfile = NULL;
  if (speedfile != NULL)
    fclose(speedfile);
  speedfile = NULL;
}


void miscInit()
{
#if defined(ILINK)
#if !defined(DOS)
  if ( checkpointStatus != checkpointedRun )	/* Shriram */
  {
#endif  /* !defined(DOS) */
    h = sqrt(exp(-nbit * log(2.0)));
    tol = tolconst;
    /* tol:=tolconst*sqrt(n);*/
    trupb = sqrt(h);
    maxit = n * iterationMultiple;
    firsttime = true;
#if !defined(DOS)
  }
  firsttime = ( checkpointStatus == normalRun ) ;    /* K. Shriram */
#endif  /* !defined(DOS) */
#elif defined(MLINK)  /* defined(ILINK) */
  firsttime = true;
#else  /* defined(ILINK) */
  firsttime = ( checkpointStatus != checkpointedRun ); /*K. Shriram*/
#endif  /* defined(ILINK) */
  
  lasttime = false;
  dolod = false;
  censorstruct = (censorrec *)Malloc(sizeof(censorrec));
  if (censorstruct == NULL)
    malloc_err("censorstruct");

  
  firstapprox = true;
  if (DIAGNOSTIC)  /* A. A. Schaffer*/
    {
      allele_downcode_check(); /* dwix */
#if ALLELE_SPEED
      adjust_alleles();
#endif
#if PARALLEL
    if (Tmk_proc_id == 0) 
#endif /*PARALLEL*/
      check_constants();
#if ALLELE_SPEED
      allele_adjust_persons();
#endif
    }
}


#if (defined(MLINK) || defined(LINKMAP))
#if !defined(DOS)
void ckptPreIpedRecover()
{
  if (checkpointedRun == checkpointStatus) {
    recoverCheckpoint(&checkpoint_place);
    checkpoint_counter = ckptTuple.ckptAttribute;    
    fclose(outfile);
    copyFile(MainOutfileDat, "outfile.dat");
#ifdef vms
    outfile = fopen("outfile.dat", "a", "ctx=rec","shr=get,put,upd");
#else  /* ifdef vms */
    outfile = fopen( "outfile.dat", "a");
#endif  /* ifdef vms */
    if (dostream) {
      fclose(stream);
      copyFile(MainStreamDat, "stream.dat");
#ifdef vms
      stream = fopen("stream.dat", "a", "ctx=rec","shr=get,put,upd");
#else  /* ifdef vms */
      stream = fopen( "stream.dat", "a");
#endif  /* ifdef vms */
    }
  }
  else{
    checkpoint_counter = 0;
    checkpoint_place = iterped_call_before;
  }
}
#endif  /* !defined(DOS) */
#endif  /* (defined(MLINK) || defined(LINKMAP)) */


void closeOutputFiles()
{
#if PARALLEL
  if (Tmk_proc_id == 0) {
#endif  /* if PARALLEL */
    if (outfile != NULL)
      fclose(outfile);
    outfile = NULL;
    if (dostream) {
#if !defined(ILINK)
      fprintf(stream, "~\n");
#endif  /* !defind(ILINK) */
      if (stream != NULL)
	fclose(stream);
      stream = NULL;
    }
#if defined(ILINK)
    if (final != NULL)
      fclose(final);
    final = NULL;
#endif  /* defined(ILINK) */
#if defined(MLINK)
    if (datafile != NULL)
      fclose(datafile);
    datafile = NULL;
#endif  /* defined(MLINK) */
#if PARALLEL
  }
#endif  /* if PARALLEL */  
}

#if !defined(DOS)
void ckptCleanup()
{
  /* Shriram: begin */

  /* This is the time at which to update the script-level
     checkpointing routine, since we have closed final and stream;
     any disasters that occur after this are "safe". */

#ifdef vms
  scriptCheckpoint = fopen(ScriptCheckpointFilename,
	  "r+", "ctx=rec","shr=get,put,upd" );
#else  /* ifdef vms */
  scriptCheckpoint = fopen ( ScriptCheckpointFilename , "r+" ) ;
#endif  /* ifdef vms */
  if ( NULL != scriptCheckpoint )
  {
    rewind ( scriptCheckpoint ) ;
    scriptRun ++ ;
    fprintf ( scriptCheckpoint , "%d 0\n" , scriptRun ) ;
    fclose ( scriptCheckpoint ) ;

    copyFile ( "final.out" , ScriptFinalOut ) ;
    appendFile ( "outfile.dat" , ScriptFinalOut ) ;
    
    if ( dostream )
    {
      copyFile ( "stream.out" , ScriptStreamOut ) ; 
      appendFile ( "stream.dat" , ScriptStreamOut ) ; 
    }
  }
  
  if ( NULL != checkpointDatafile )
    fclose ( checkpointDatafile ) ;
  
  unlink ( CheckpointFileBackup ) ;
  unlink ( CheckpointFilename ) ;
  unlink ( MainOutfileDat );
  unlink ( MainStreamDat  );
  /* Shriram: end */
}
#endif  /* if !defined(DOS) */


void initParams()
{
  /*find number of recombination probabilities*/  
  nuneed = 7;
  for (i= 3; i<= mlocus; i++)
    nuneed = 5 * nuneed - 3;

#if !defined(LESSMEMORY)  
  /*find size of isozygote class and it's square for joint classes*/
  maxclasssize = 2;
  for (i = 3; i<= mlocus; i++)
    maxclasssize *= 2;
  maxisozygclass = maxclasssize * maxclasssize;
#endif
  
  nuprobclass = 2;
  for(i = 2; i <= mlocus; i++)
    nuprobclass = 3 * nuprobclass - 1;

#if defined(ILINK)
  initialize();
  ihess = 0;
  ibnd = 1;
  icall = 1;
  ivar = 0;
  ihx = 1;
  inconsistent = false;
#endif  /* defined(ILINK) */
}


/* sequential startup code */
void seqStartup(argc, argv)
int argc;
char** argv;
{
  int c;

  disfreqs = false;
#ifdef vms
  ;
#else
  while ((c = getopt(argc, argv, "ic")) != -1)
    switch (c) {
    case 'i':
      printInfo();
      break;
    case 'c':
      disfreqs = true;
      break;
    case '?':
      /* fprintf(stderr, "Unrecognized option\n"); */
      exit(-1);
    default:
      break;
    }
#endif
}



#if (defined(MLINK) || defined(LINKMAP))

void checkzero()
{
  int i;  /* cgh */
  
  if (!firsttime) {
    for (i = 1; i < mlocus; i++) {
      if (maletheta->theta[i - 1] != 0.0 && zeromale[i - 1])
	firsttime = true;
    }
    for (i = 1; i < mlocus; i++) {
      if (femaletheta->theta[i - 1] != 0.0 && zerofemale[i - 1])
	firsttime = true;
    }
  }
#if defined(MLINK)
  if (maletheta->theta[which - 1] == 0.0)
    firsttime = true;
  if (femaletheta->theta[which - 1] == 0.0)
    firsttime = true;
#elif defined(LINKMAP)
  if (maletheta->theta[whichvary - 1] == 0.0)
    firsttime = true;
  if (femaletheta->theta[whichvary - 1] == 0.0)
    firsttime = true;
#endif  /* defined(MLINK) */
  if (!firsttime)
    return;
  for (i = 1; i < mlocus; i++)
    zeromale[i - 1] = (maletheta->theta[i - 1] == 0.0);
  for (i = 1; i < mlocus; i++)
    zerofemale[i - 1] = (femaletheta->theta[i - 1] == 0.0);
}


/*The following procedure was introduced by A. A. Schaffer to
catch components of theta vectors that should be 0.0 but aren't*/
boolean zerotest(thetacomponent)
double thetacomponent;
{
  return((thetacomponent < epsilon) && (thetacomponent >= 0.0));
}


/* This main() serves for both MLINK and LINKMAP */
main(argc, argv)
int argc;
char *argv[];
{
  int numgens;

#if PARALLEL  /* cgh */
  parStartup(argc, argv);
  if (Tmk_proc_id == 0)
#else
  seqStartup(argc, argv);
#endif  /* if PARALLEL -- cgh */
    initialize();
  
#if !defined(DOS)
  checkpointStatus = normalRun;
#endif  /* if !defined(DOS) */

  init_ped_loc_all(); /* dwix */
  openFiles();

#if !defined(DOS)   
  ckptInit();
#endif  /* if !defined(DOS) */

#if PARALLEL
#if !IS_SHMEM
  if (reportMemStats == false)
    gMemAndBarrierAlloc();
#endif  /* !IS_SHMEM */
  allocthetas();
#endif  /* if PARALLEL */
  
  inputdata();
  if (interfer) {
    fprintf(stderr,"Interference is not supported in mlink or linkmap\n");
    exit(EXIT_FAILURE);
  }
  initParams();
  
#if PARALLEL  /* cgh */
  if (Tmk_proc_id == 0)
    checkNotImplemented();
#endif  /* if PARALLEL -- cgh */
  
#if !defined(DOS)
  ckptEnsureFileContents();
#endif  /* if !defined(DOS) */

  miscInit();
#if !ALLELE_SPEED
  allocategenetables();
  getlocations();
#endif
  getgeneindices();  /*R. M. Idury*/
  allocate_loopbreaker_vectors(); /*A. A. Schaffer*/
  closeInputFiles();
#if PARALLEL  
  if (0 == maxworkingset) {   /*number was not input by user*/
    maxworkingset = maxw_estimation();
    if (Tmk_proc_id == 0)
      printf("\nEstimating %d for maxworkingset\n", maxworkingset);
  }
  else 
    printf("\nOverriding computed value of maxworkingset with input value %d\n",maxworkingset);
#endif

#if PARALLEL  /* cgh */
  if (Tmk_proc_id == 0) {
#endif  /* if PARALLEL -- cgh */

  preIpedLoop();
  
#if PARALLEL  /* cgh */
  /* simulate loop to calculate the number of calls to iterpeds() */
  simIpedLoop(countIpeds);
  initParLib(argc, argv);

#if IS_SHMEM  
  gMemAndBarrierAlloc();
#endif  /* IS_SHMEM */

  initOutBuffs();        /* initialize output buffers */
  allocgen();            /*Added by Alex*/
#if PRECOMPUTE
  allocprob();
#endif  /* if PRECOMPUTE */
  parThetaSetup();

  /* simulate loop iterpeds() loop again, this time buffering output
     to stream, as well as computing and storing theta values */
  simIpedLoop(computeThetas);
  }  /* if (Tmk_proc_id == 0) */
#endif  /* if PARALLEL -- cgh */

#if !defined(DOS)
  ckptPreIpedRecover();
#endif  /* if !defined(DOS) */

#if !PARALLEL  /* cgh */
  ipedLoop();
#endif  /* if !PARALLEL -- cgh */
  
#if PARALLEL  /* cgh */
  childLoop();
  if (Tmk_proc_id == 0) {
    writeOutBuffs();  /* done -- write buffered output to files */
#endif  /* if PARALLEL -- cgh */
    closeOutputFiles();
#if PARALLEL
  }  /* if (Tmk_proc_id == 0) */
#endif  /* PARALLEL */

#if !defined(DOS)
  ckptCleanup();
#endif  /* if !defined(DOS) */
  
#if PARALLEL
  parFinish();
#endif  /* if PARALLEL */
   
  exit(EXIT_SUCCESS);
}

#endif  /* defined(MLINK) || defined(LINKMAP) */

