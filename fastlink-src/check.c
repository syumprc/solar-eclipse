/* This file contains the various routines used to perform checkpointing
   in ILINK, LODSCORE, LINKMAP, or MLINK.
   When compiling it, provide either ILINK. LODSCORE, LINKMAP, or MLINK
   LODSCORE as a #define'd value; this can be done on the command line
   for standard Unix C compilers with

   % cc -DILINK <rest of command line>

   Some of the checkpointing code is also to be found in ilink.c,
   lodscore.c, linkmap.c, and mlink.c
   within the specific functions that deal with checkpointing. */

/* The code in this file was written by K. Shriram */
/* The checkpointing scheme is described in the paper: */
/* A. A. Schaffer, S. K. Gupta, K. Shriram, and R. W. Cottingham, Jr. */
/* Avoiding Recomputation in Linkage Analysis,*/
/*   Hum. Hered. 44(1994), pp. 225-237*/

#include "commondefs.h"
#include "checkpointdefs.h"
#include "gemdefs.h"
#include "ildefs.h"

/* While copying or appending files, if we have an error, we refer it
   to fileErrorHelp(). */

void
fileErrorHelp ( fileName , problemType )
  char  * fileName ;
  int  problemType ;
{
  char  theLine [ 256 ] ;
  
  puts ( "A problem was encountered with the file" ) ;
  printf ( "    %s\n" , fileName ) ;
  printf ( "while attempting to " ) ;
  switch ( problemType )
  {
  case 1:
    puts ( "open it for reading." ) ;
    break ;
  case 2:
    puts ( "open it for writing to." ) ;
    break ;
  case 3:
    puts ( "append to it." ) ;
    break ;
  case 4:
    puts ( "write data to it." ) ;
    break ;
  }
  printf ( "\nThis will possibly cause further errors, and will also\n" ) ;
  puts ( "prevent correct recovery from a checkpoint.  Hence, it is" ) ;
  puts ( "recommended that you break the execution of this program at" ) ;
  puts ( "this time and re-start computation." ) ;

  gets ( theLine ) ;
}

/* In the process of checkpointing, we need to make copies of files;
   sometimes we also need to append one file to another, to mimic the
   operations that take place in the scripts.  ourCopyAppendFile() does
   that; it is front-ended by the macro definitions copyFile() and
   appendFile() in checkpointdefs.h. */

void
ourCopyAppendFile ( fromName , toName , operationType )
  char  * fromName ;
  char  * toName ;
  int  operationType ;
{
  int  flags ;
  int  fromFile ;
  int  toFile ;
  int  bytesRead ;
  char  theBuffer [ CopyBlockSize ] ;

  if ( CopyOperation == operationType )
    flags = O_WRONLY | O_CREAT ;
  else	/* AppendOperation == operationType */
    flags = O_WRONLY | O_APPEND | O_CREAT ;

#ifdef vms
  fromFile = open ( fromName , O_RDONLY, 0, "ctx=rec", "shr=get,put,upd" ) ;
#else
  fromFile = open ( fromName , O_RDONLY ) ;
#endif
  if ( -1 == fromFile )
    fileErrorHelp ( fromName , 1 ) ;
  toFile = open ( toName , flags , CopyAppendPerms ) ;
  if ( -1 == toFile )
    fileErrorHelp ( toName ,
		    ( CopyOperation == operationType ) ? 2 : 3 ) ;

  do
  {
    bytesRead = read ( fromFile , theBuffer , CopyBlockSize ) ;
    if ( 0 < bytesRead )
      if ( bytesRead != write ( toFile , theBuffer , bytesRead ) )
	fileErrorHelp ( toName , 4 ) ;
  }
  while ( 0 < bytesRead ) ;

  close ( fromFile ) ;
  close ( toFile ) ;
}

/* Another routine to deal with errors; this time if the checkpoint
   file being read in isn't complete (ie, the end marker is not found
   at the appropriate point).  There are two options; either a backup
   file is found or it isn't.  The authenticity of the backup is not
   examined. */

void  recoveryErrorHelp ( )
{
  char  theLine [ 256 ] ;
  FILE  * helpTextFile ;
  FILE  * checkpointBackup ;

#ifdef ILINK
  checkpointBackup = fopen ( CheckpointILINKFileBackup , "r" ) ;
#elif LODSCORE
  checkpointBackup = fopen ( CheckpointLODSCOREFileBackup , "r" ) ;
#elif LINKMAP
  checkpointBackup = fopen ( CheckpointLINKMAPFileBackup , "r" ) ;
#elif MLINK
  checkpointBackup = fopen ( CheckpointMLINKFileBackup , "r" ) ;
#endif
  if ( NULL == checkpointBackup ) {
    fprintf(stderr,"\nIn the process of recovering, an error has been detected.  The most");
    fprintf(stderr,"\nlikely reason is that a system error such as a machine crash took");
    fprintf(stderr,"\nplace while the program was saving certain important information");
    fprintf(stderr,"\nnecessary for the recovery process.");
    fprintf(stderr,"\n");
    fprintf(stderr,"\nThis program attempts to create back-ups of such information from time");
    fprintf(stderr,"\nto time.  However, no such backup is to be found.  The most likely");
    fprintf(stderr,"\nreasons for this are:");
    fprintf(stderr,"\n");
    fprintf(stderr,"\n(a) The program had not been executing for long enough when the crash");
    fprintf(stderr,"\noccurred.");
    fprintf(stderr,"\n(b) The file was deleted during the crash.");
    fprintf(stderr,"\n(c) This run is taking place in a different directory than the");
    fprintf(stderr,"\nprevious one did.");
    fprintf(stderr,"\n");
    fprintf(stderr,"\nThe program will now await your input.  The best course of action");
    fprintf(stderr,"\nmight be to 'break' this run, delete the checkpoint file and resume.");
    fprintf(stderr,"\nShould this error occur repeatedly, please report it to the distributors.");
    fprintf(stderr,"\nTo continue, hit <return> or <enter>.");
    fprintf(stderr,"\n");
  }
  else {
  fprintf(stderr,"\nIn the process of recovering, an error has been detected.  The most");
  fprintf(stderr,"\nlikely reason is that a system error (such as a machine crash) took");
  fprintf(stderr,"\nplace while the program was saving certain important information");
  fprintf(stderr,"\nnecessary for the recovery process.");
  fprintf(stderr,"\n");
  fprintf(stderr,"\nA backup from an earlier save has been located.  The program will");
  fprintf(stderr,"\nsuspend after displaying this message and will await your input; the");
  fprintf(stderr,"\nfollowing course of action is suggested:");
  fprintf(stderr,"\n");
  fprintf(stderr,"\n(a) Stop execution by 'breaking' the program.");
  fprintf(stderr,"\n(b) Copy the backup file to the primary file.");
  fprintf(stderr,"\n(c) Run the program again, as before.");
  fprintf(stderr,"\n");
  fprintf(stderr,"\nHopefully, the data recovery will now proceed smoothly.  Should you");
  fprintf(stderr,"\nstill get this error, however, please contact the distributors of this");
  fprintf(stderr,"\nprogram so that the problem may be evaluated.");
  fprintf(stderr,"\n");
  fprintf(stderr,"\nIf you wish to continue execution with the existing recovered data,");
  fprintf(stderr,"\nhit <return> or <enter>.");
  fprintf(stderr,"\n");
  }
  fclose ( checkpointBackup ) ;
  checkpointBackup = NULL;
  gets ( theLine ) ;
}

/* Puts the information tuple at the head of the checkpoint file. */

void  putCkptTuple ( )
{
  fprintf ( checkpointDatafile , "%d %d %d %d\n" ,
	   ckptTuple . ckptLocation , ckptTuple . ckptAttribute ,
	   ckptTuple . iplace , ckptTuple . jplace ) ;
}

/* Reads the information tuple. */

void  getCkptTuple ( )
{
  int  newlineDummy ;

  /* cgh -- cast &ckptTuple.ckptLocation to int* for gcc warning */
  fscanf ( checkpointDatafile , "%d %d %d %d" ,
	  (int*) & ckptTuple . ckptLocation , & ckptTuple . ckptAttribute ,
	  & ckptTuple . iplace , & ckptTuple . jplace ) ;
  newlineDummy = getc ( checkpointDatafile ) ;
  if ( '\n' != newlineDummy )
    recoveryErrorHelp ( ) ;
}

/* Macro front-ends to putCkptNumber(). */

#define  putCkptFloat(f)	( putCkptNumber ( & (f) , sizeof ( float ) ) )
#define  putCkptDouble(d)	( putCkptNumber ( & (d) , sizeof ( double ) ) )
#define	 putCkptInt(i)		( putCkptNumber ( & (i) , sizeof ( int ) ) )
#define  putCkptLong(l)		( putCkptNumber ( & (l) , sizeof ( long ) ) )
#define  putCkptBoolean(b)	( putCkptNumber ( & (b) , sizeof ( boolean ) ) )

/* General routine for putting the internal representation of a number. */

void  putCkptNumber ( numberPtr , numberSize )
void  * numberPtr ;
int   numberSize ;
{
  fwrite ( numberPtr , numberSize , 1 , checkpointDatafile ) ;
}

/* Routine for putting an object of type vector. */

void    putCkptVector ( theVector )
vector  theVector ;
{
  int  vectorRef ;

  for ( vectorRef = 0 ; vectorRef < maxn ; vectorRef ++ )
    putCkptDouble ( theVector [ vectorRef ] ) ;
}

/* Code to put a matrix, vector-by-vector. */

void    putCkptMatrix ( theMatrix )
matrix  theMatrix ;
{
  int  matrixRef ;

  for ( matrixRef = 0 ; matrixRef < maxn ; matrixRef ++ )
    putCkptVector ( theMatrix [ matrixRef ] ) ;
}

/* Code to put an iter type. */

void      putCkptIter ( theIter )
itertype  theIter ;
{
  int  iterRef ;

  for ( iterRef = 0 ; iterRef < maxn ; iterRef ++ )
    putCkptInt ( theIter [ iterRef ] ) ;
}

/* Routine that makes calls to write out all the variables we are going
   to checkpoint.  Large and boring. */

void  writeCheckpointInformation ()
{
  int  rowCtr , colCtr ;

  /* from MIN1 */
  putCkptMatrix ( tl ) ;
  putCkptVector ( d ) ;  putCkptVector ( g ) ;  putCkptVector ( gsave ) ;
  putCkptVector ( y ) ;  putCkptVector ( p ) ;
  
  /* from MIN2 */
  putCkptInt ( nit ) ;  putCkptInt ( nfe ) ;  putCkptInt ( idg ) ;
  putCkptInt ( idif ) ;  putCkptInt ( isw ) ;  putCkptInt ( iret ) ;
  putCkptInt ( ibnd ) ;  putCkptInt ( ihess ) ;  putCkptInt ( ivar ) ;
  putCkptInt ( ihx ) ;  putCkptInt ( maxit ) ;
  putCkptDouble ( tol ) ;  putCkptDouble ( tmin ) ;  putCkptDouble ( h ) ;
  putCkptDouble ( trupb ) ;  putCkptDouble ( ptg ) ;

  /* from GEMINI */
  putCkptVector ( xall ) ;
  putCkptVector ( x ) ;
  putCkptVector ( v ) ;     putCkptVector ( se ) ;
  putCkptIter ( itp ) ;
  for ( rowCtr = 0 ; rowCtr < maxn ; rowCtr ++ )
    for ( colCtr = 0 ; colCtr < 2 ; colCtr ++ )
      putCkptDouble ( bnd [ rowCtr ] [ colCtr ] ) ;
  putCkptInt ( nall ) ;  putCkptInt ( n ) ;  putCkptInt ( icall ) ;
  putCkptDouble ( tbnd ) ;  putCkptDouble ( f ) ;  putCkptDouble ( fsmf ) ;
  putCkptDouble ( fsav2 ) ;  putCkptDouble ( t ) ;  putCkptDouble ( hx ) ;
  putCkptDouble ( xsave ) ;  putCkptDouble ( fxph ) ;  putCkptDouble ( fxmh ) ;
  putCkptDouble ( xp ) ;  putCkptDouble ( xpm ) ;  putCkptDouble ( ytp ) ;
  
  /* from ALEX */
  putCkptVector ( outsavex ) ;
  putCkptDouble ( outsavefvalue ) ;
  for ( rowCtr = 0 ; rowCtr < maxlocus ; rowCtr ++ )
    putCkptDouble ( savedf [ rowCtr ] ) ;

  /* from UPDATE */
  putCkptVector ( wtil ) ;  putCkptVector ( ztil ) ;  putCkptVector ( w ) ;
  putCkptVector ( z ) ;  putCkptVector ( wtjp1 ) ;  putCkptVector ( ztjp1 ) ;
  putCkptVector ( s ) ;  putCkptVector ( dp ) ;
  putCkptDouble ( nu ) ;  putCkptDouble ( muj ) ;  putCkptDouble ( lambj ) ;
  putCkptDouble ( lambj2 ) ;  putCkptDouble ( sb ) ;  putCkptDouble ( sc ) ;
  putCkptDouble ( sd ) ;  putCkptDouble ( fbcd ) ;  putCkptDouble ( alpha ) ;
  putCkptDouble ( sa ) ;  putCkptDouble ( thet1 ) ;  putCkptDouble ( thet2 ) ;
  putCkptDouble ( aa ) ; putCkptDouble ( bb ) ;  putCkptDouble ( cc ) ;
  putCkptDouble ( del2 ) ;  putCkptDouble ( alph1 ) ;  putCkptDouble ( alph2 ) ;
  putCkptDouble ( rdiv ) ;  putCkptDouble ( eta ) ;  putCkptDouble ( aj ) ;
  putCkptDouble ( thj ) ;  putCkptDouble ( bj ) ;  putCkptDouble ( gamlj ) ;
  putCkptDouble ( betlj ) ;  putCkptDouble ( del ) ;
  putCkptInt ( iswup ) ;

  /* from STEP */
  putCkptVector ( xt ) ;
  putCkptDouble ( fsave ) ;  putCkptDouble ( sumt ) ;  putCkptDouble ( twot ) ;
  putCkptDouble ( ft ) ;  putCkptDouble ( f2t ) ;  putCkptDouble ( ft2 ) ;
  putCkptDouble ( scal ) ;
  
  /* from CHKBND */
  putCkptDouble ( clb ) ;  putCkptDouble ( check ) ;  putCkptDouble ( eh ) ;
  putCkptDouble ( teq ) ;
  
  /* from OTHERS */
  putCkptInt ( itsys ) ;

  /* from INIB */
  putCkptMatrix ( bmat ) ;

  /* Miscellaneous */

  putCkptInt ( continue_ ) ;
  putCkptBoolean ( firstapprox ) ;
  putCkptBoolean ( firsttime ) ;
  putCkptBoolean ( lasttime ) ;
  putCkptBoolean ( inconsistent ) ;
  putCkptBoolean ( dolod ) ;
  putCkptInt ( thisc ) ;
  putCkptDouble ( penlike ) ;
  putCkptDouble ( like ) ;
  for ( rowCtr = 0 ; rowCtr < maxped ; rowCtr ++ )
  {
    putCkptDouble ( likebyped [ rowCtr ] ) ;
    putCkptDouble ( outsavelike [ rowCtr ] ) ;
  }
  for ( rowCtr = 0 ; rowCtr < maxlocus ; rowCtr ++ )
  {
    putCkptDouble ( maletheta -> theta [ rowCtr ] ) ;
    putCkptDouble ( femaletheta -> theta [ rowCtr ] ) ;
  }
  for ( rowCtr = 0 ; rowCtr < nuneed ; rowCtr ++ )
  {
    putCkptDouble ( maletheta -> segprob [ rowCtr ] ) ;
    putCkptDouble ( femaletheta -> segprob [ rowCtr ] ) ;
  }
}

/* The actual routine that gets called when we want to checkpoint.  It
   takes arguments telling it how where the checkpoint is being made,
   and what sequence of calls got us there.  The extra variable takes
   some extra information about local variables; whatever value is
   sent is always written, and it is for the call to recoverCheckpoint()
   to ignore it or not. */

void                performCheckpoint ( locationType , locationDatum , extra )
checkpointLocation  locationType ;
int                 locationDatum ;
int		    extra ;
{
  time_t  secondsNow ;

  if ( iterationLocation == locationType )
    return ;
 
/*  printf ( "Checkpointing: type %d location %d\n" , locationType , locationDatum ) ; */

#ifdef ILINK
  rename ( CheckpointILINKFilename , CheckpointILINKFileBackup );
#elif LODSCORE
  rename ( CheckpointLODSCOREFilename , CheckpointLODSCOREFileBackup ) ;
#elif LINKMAP
  rename ( CheckpointLINKMAPFilename , CheckpointLINKMAPFileBackup ) ;
#elif MLINK
  rename ( CheckpointMLINKFilename , CheckpointMLINKFileBackup ) ;
#endif

#ifdef ILINK
  checkpointDatafile = fopen ( CheckpointILINKFilename , "w" ) ;
#elif LODSCORE
  checkpointDatafile = fopen ( CheckpointLODSCOREFilename , "w" ) ;
#elif LINKMAP
  checkpointDatafile = fopen ( CheckpointLINKMAPFilename , "w" ) ;
#elif MLINK
  checkpointDatafile = fopen ( CheckpointMLINKFilename , "w" ) ;
#endif
  if ( NULL != checkpointDatafile )
  {
    time ( & secondsNow ) ;
    fprintf ( checkpointDatafile , "%s" , ctime ( & secondsNow ) ) ;
    /* ctime() automatically adds a newline at the end of the string */
/*    fprintf ( checkpointDatafile , "%d    %d\n" ,
	     locationType , locationDatum ) ;
	     */
    ckptTuple . ckptLocation = locationType ;
    ckptTuple . ckptAttribute = locationDatum ;
    putCkptTuple ( ) ;
#ifdef ILINK 
    writeCheckpointInformation ( ) ;
#elif LODSCORE 
    writeCheckpointInformation ( ) ;
#endif
    fprintf ( checkpointDatafile , "\n%d\n" , extra ) ;
    fprintf ( checkpointDatafile , "%s\n" , EndOfDataString ) ;
    ensureWrite ( checkpointDatafile ) ;
    fclose ( checkpointDatafile ) ;
    checkpointDatafile = NULL;
  }
}

/* Front-ends to retrieve numbers. */

#define  getCkptFloat(f)	( getCkptNumber ( & (f) , sizeof ( float ) ) )
#define  getCkptDouble(d)	( getCkptNumber ( & (d) , sizeof ( double ) ) )
#define	 getCkptInt(i)		( getCkptNumber ( & (i) , sizeof ( int ) ) )
#define  getCkptLong(l)		( getCkptNumber ( & (l) , sizeof ( long ) ) )
#define  getCkptBoolean(b)	( getCkptNumber ( & (b) , sizeof ( boolean ) ) )

/* Get a number in machine-representation form. */

void  getCkptNumber ( numberPtr , numberSize )
void  * numberPtr ;
int   numberSize ;
{
  fread ( numberPtr , numberSize , 1 , checkpointDatafile ) ;
}

/* Get a vector. */

void    getCkptVector ( theVector )
vector  theVector ;
{
  int    vectorRef ;

  for ( vectorRef = 0 ; vectorRef < maxn ; vectorRef ++ )
    getCkptDouble ( theVector [ vectorRef ] ) ;
}

/* Get a matrix by vector. */

void    getCkptMatrix ( theMatrix )
matrix  theMatrix ;
{
  int  matrixRef ;

  for ( matrixRef = 0 ; matrixRef < maxn ; matrixRef ++ )
    getCkptVector ( theMatrix [ matrixRef ] ) ;
}

/* Get an iter type. */

void      getCkptIter ( theIter )
itertype  theIter ;
{
  int   iterRef ;

  for ( iterRef = 0 ; iterRef < maxn ; iterRef ++ )
    getCkptInt ( theIter [ iterRef ] ) ;
}

/* Read back all the variables we wrote. */
/* cgh -- made this void */
void readCheckpointInformation ()
{
  int    rowCtr , colCtr ;

  /* from MIN1 */
  getCkptMatrix ( tl ) ;
  getCkptVector ( d ) ;  getCkptVector ( g ) ;  getCkptVector ( gsave ) ;
  getCkptVector ( y ) ;  getCkptVector ( p ) ;
  
  /* from MIN2 */
  getCkptInt ( nit ) ;  getCkptInt ( nfe ) ;  getCkptInt ( idg ) ;
  getCkptInt ( idif ) ;  getCkptInt ( isw ) ;  getCkptInt ( iret ) ;
  getCkptInt ( ibnd ) ;  getCkptInt ( ihess ) ;  getCkptInt ( ivar ) ;
  getCkptInt ( ihx ) ;  getCkptInt ( maxit ) ;
  getCkptDouble ( tol ) ;  getCkptDouble ( tmin ) ;  getCkptDouble ( h ) ;
  getCkptDouble ( trupb ) ;  getCkptDouble ( ptg ) ;

  /* from GEMINI */
  getCkptVector ( xall ) ;
  getCkptVector ( x ) ;
  getCkptVector ( v ) ;     getCkptVector ( se ) ;
  getCkptIter ( itp ) ;
  for ( rowCtr = 0 ; rowCtr < maxn ; rowCtr ++ )
    for ( colCtr = 0 ; colCtr < 2 ; colCtr ++ )
      getCkptDouble ( bnd [ rowCtr ] [ colCtr ] ) ;
  getCkptInt ( nall ) ;  getCkptInt ( n ) ;  getCkptInt ( icall ) ;
  getCkptDouble ( tbnd ) ;  getCkptDouble ( f ) ;  getCkptDouble ( fsmf ) ;
  getCkptDouble ( fsav2 ) ;  getCkptDouble ( t ) ;  getCkptDouble ( hx ) ;
  getCkptDouble ( xsave ) ;  getCkptDouble ( fxph ) ;  getCkptDouble ( fxmh ) ;
  getCkptDouble ( xp ) ;  getCkptDouble ( xpm ) ;  getCkptDouble ( ytp ) ;
  
  /* from ALEX */
  getCkptVector ( outsavex ) ;
  getCkptDouble ( outsavefvalue ) ;
  for ( rowCtr = 0 ; rowCtr < maxlocus ; rowCtr ++ )
    getCkptDouble ( savedf [ rowCtr ] ) ;

  /* from UPDATE */
  getCkptVector ( wtil ) ;  getCkptVector ( ztil ) ;  getCkptVector ( w ) ;
  getCkptVector ( z ) ;  getCkptVector ( wtjp1 ) ;  getCkptVector ( ztjp1 ) ;
  getCkptVector ( s ) ;  getCkptVector ( dp ) ;
  getCkptDouble ( nu ) ;  getCkptDouble ( muj ) ;  getCkptDouble ( lambj ) ;
  getCkptDouble ( lambj2 ) ;  getCkptDouble ( sb ) ;  getCkptDouble ( sc ) ;
  getCkptDouble ( sd ) ;  getCkptDouble ( fbcd ) ;  getCkptDouble ( alpha ) ;
  getCkptDouble ( sa ) ;  getCkptDouble ( thet1 ) ;  getCkptDouble ( thet2 ) ;
  getCkptDouble ( aa ) ; getCkptDouble ( bb ) ;  getCkptDouble ( cc ) ;
  getCkptDouble ( del2 ) ;  getCkptDouble ( alph1 ) ;  getCkptDouble ( alph2 ) ;
  getCkptDouble ( rdiv ) ;  getCkptDouble ( eta ) ;  getCkptDouble ( aj ) ;
  getCkptDouble ( thj ) ;  getCkptDouble ( bj ) ;  getCkptDouble ( gamlj ) ;
  getCkptDouble ( betlj ) ;  getCkptDouble ( del ) ;
  getCkptInt ( iswup ) ;

  /* from STEP */
  getCkptVector ( xt ) ;
  getCkptDouble ( fsave ) ;  getCkptDouble ( sumt ) ;  getCkptDouble ( twot ) ;
  getCkptDouble ( ft ) ;  getCkptDouble ( f2t ) ;  getCkptDouble ( ft2 ) ;
  getCkptDouble ( scal ) ;
  
  /* from CHKBND */
  getCkptDouble ( clb ) ;  getCkptDouble ( check ) ;  getCkptDouble ( eh ) ;
  getCkptDouble ( teq ) ;
  
  /* from OTHERS */
  getCkptInt ( itsys ) ;

  /* from INIB */
  getCkptMatrix ( bmat ) ;
  
  /* Miscellaneous */
  getCkptInt ( continue_ ) ;
  getCkptBoolean ( firstapprox ) ;
  getCkptBoolean ( firsttime ) ;
  getCkptBoolean ( lasttime ) ;
  getCkptBoolean ( inconsistent ) ;
  getCkptBoolean ( dolod ) ;
  getCkptInt ( thisc ) ;
  getCkptDouble ( penlike ) ;
  getCkptDouble ( like ) ;
  for ( rowCtr = 0 ; rowCtr < maxped ; rowCtr ++ )
  {
    getCkptDouble ( likebyped [ rowCtr ] ) ;
    getCkptDouble ( outsavelike [ rowCtr ] ) ;
  }
  for ( rowCtr = 0 ; rowCtr < maxlocus ; rowCtr ++ )
  {
    getCkptDouble ( maletheta -> theta [ rowCtr ] ) ;
    getCkptDouble ( femaletheta -> theta [ rowCtr ] ) ;
  }
  for ( rowCtr = 0 ; rowCtr < nuneed ; rowCtr ++ )
  {
    getCkptDouble ( maletheta -> segprob [ rowCtr ] ) ;
    getCkptDouble ( femaletheta -> segprob [ rowCtr ] ) ;
  }
}

/* Function called when recovery is desired.  The only argument is a
   pointer to the extra data which (the pointer), if null, ignores the
   value recovered.  Note the check being made for the end of data marker. */

void  recoverCheckpoint ( extra )
int  * extra ;
{
  char   	      dateTimeStamp [ DateTimeStampStringLength ] ;
  char		      theLine [ 80 ] ;
/* cgh -- unused variables
   checkpointLocation  locationType ;
   int                 locationDatum ; */
  int                 dummy ;

#ifdef ILINK
  checkpointDatafile = fopen ( CheckpointILINKFilename , "r" ) ;
#elif LODSCORE
  checkpointDatafile = fopen ( CheckpointLODSCOREFilename , "r" ) ;
#elif LINKMAP
  checkpointDatafile = fopen ( CheckpointLINKMAPFilename , "r" ) ;
#elif MLINK
  checkpointDatafile = fopen ( CheckpointMLINKFilename , "r" ) ;
#endif

  fgets ( dateTimeStamp , DateTimeStampStringLength , checkpointDatafile ) ;
  getCkptTuple ( ) ;
#ifdef ILINK 
  readCheckpointInformation ( ) ;
#elif LODSCORE
  readCheckpointInformation ( ) ;
#endif
  if ( NULL != extra )
    fscanf ( checkpointDatafile , "%d" , extra ) ;
  else
    fscanf ( checkpointDatafile , "%d" , & dummy ) ;
  fscanf ( checkpointDatafile , "%s" , theLine ) ;
  if ( ! ( strcmp ( theLine , EndOfDataString ) ) )
    puts ( "Data recovered" ) ;
  else
    recoveryErrorHelp ( ) ;
  fclose ( checkpointDatafile ) ;
  checkpointDatafile = NULL;
#ifdef ILINK
  checkpointStatus = normalRun ;
#elif LODSCORE
  checkpointStatus = normalRun ;
#endif
}
