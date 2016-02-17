/* This file contains the definitions that are used in checkpoint.c, and
   in ilink.c or lodscore.c (as the case may be).  Unlike checkpoint.c,
   this file does not use the ILINK or LODSCORE #define's; it could, but
   it likely isn't worth the bother.  Note that distinct names are set up
   for each of ILINK and LODSCORE. */
/* The checkpointing code was written by K. Shriram */
/* The checkpointing process is described in the paper: */
/* A. A. Schaffer, S. K. Gupta, K. Shriram, and R. W. Cottingham, Jr. */
/* Avoiding Recomputation in Linkage Analysis,*/
/* Human Heredity 44(1994), pp. 225-237. */

#ifndef  _CHECKPOINTDEFS_H
#define  _CHECKPOINTDEFS_H  1


/* When recovering and the primary checkpoint file is corrupted ... */

#define  HelpTextFoundBackup		"recoveryFoundText"
#define  HelpTextNotFoundBackup		"recoveryNotFoundText"

/* String marker to check for corruption of the primary checkpoint */

#define  EndOfDataString 		"END_OF_DATA"

/* The grand old checkpoint files, which contain the interation- and
   function-checkpoints */

#define CheckpointILINKFilename		"./checkpoint.ILINK"
#define CheckpointILINKFileBackup	"./checkpoint.ILINK.bak"
#define CheckpointLODSCOREFilename	"./checkpoint.LODSCORE"
#define	CheckpointLODSCOREFileBackup	"./checkpoint.LODSCORE.bak"
#define CheckpointLINKMAPFilename	"./checkpoint.LINKMAP"
#define	CheckpointLINKMAPFileBackup	"./checkpoint.LINKMAP.bak"
#define CheckpointMLINKFilename	        "./checkpoint.MLINK"
#define	CheckpointMLINKFileBackup	"./checkpoint.MLINK.bak"

/* We need to store information about how many runs we have made in the
   current script and how many have been skipped over during recovery */

#define ScriptILINKCheckpointFilename	 "./script.checkpoint.ILINK"
#define ScriptLODSCORECheckpointFilename "./script.checkpoint.LODSCORE"
#define ScriptLINKMAPCheckpointFilename  "./script.checkpoint.LINKMAP"
#define ScriptMLINKCheckpointFilename    "./script.checkpoint.MLINK"

/* We maintain copies of final.out and stream.out */

#define ScriptILINKFinalOut		"script.ILINK.final.out"
#define ScriptILINKStreamOut		"script.ILINK.stream.out"
#define ScriptLODSCOREFinalOut		"script.LODSCORE.final.out"
#define ScriptLODSCOREStreamOut		"script.LODSCORE.stream.out"
#define ScriptLINKMAPFinalOut		"script.LINKMAP.final.out"
#define ScriptLINKMAPStreamOut		"script.LINKMAP.stream.out"
#define ScriptMLINKFinalOut		"script.MLINK.final.out"
#define ScriptMLINKStreamOut		"script.MLINK.stream.out"

/* outf() stores *it's* own copies of final, stream and recfile */

#define OutfILINKFinalDat		"outf.ILINK.final.dat"
#define	OutfILINKStreamDat		"outf.ILINK.stream.dat"
#define OutfLODSCORERecfileDat		"outf.LODSCORE.recfile.dat"
#define OutfLODSCOREStreamDat		"outf.LODSCORE.stream.dat"

/* main() stores copies of recfile and stream in LODSCORE */

#define MainLODSCORERecfileDat		"main.LODSCORE.recfile.dat"
#define MainLODSCOREStreamDat		"main.LODSCORE.stream.dat"
#define MainLINKMAPStreamDat		"main.LINKMAP.stream.dat"
#define MainLINKMAPOutfileDat		"main.LINKMAP.outfile.dat"
#define MainMLINKStreamDat		"main.MLINK.stream.dat"
#define MainMLINKOutfileDat		"main.MLINK.outfile.dat"

/* cgh -- abstracted filenames for sharing code */

#if defined(MLINK)
#define ScriptCheckpointFilename    ScriptMLINKCheckpointFilename
#define ScriptFinalOut              ScriptMLINKFinalOut
#define ScriptStreamOut             ScriptMLINKStreamOut
#define CheckpointFilename          CheckpointMLINKFilename
#define CheckpointFileBackup        CheckpointMLINKFileBackup
#define MainOutfileDat              MainMLINKOutfileDat
#define MainStreamDat               MainMLINKStreamDat
#endif  /* if defined(MLINK) */

#if defined(LINKMAP)
#define ScriptCheckpointFilename    ScriptLINKMAPCheckpointFilename
#define ScriptFinalOut              ScriptLINKMAPFinalOut
#define ScriptStreamOut             ScriptLINKMAPStreamOut
#define CheckpointFilename          CheckpointLINKMAPFilename
#define CheckpointFileBackup        CheckpointLINKMAPFileBackup
#define MainOutfileDat              MainLINKMAPOutfileDat
#define MainStreamDat               MainLINKMAPStreamDat
#endif  /* if defined(LINKMAP) */

#if defined(ILINK)
#define ScriptCheckpointFilename    ScriptILINKCheckpointFilename
#define ScriptFinalOut              ScriptILINKFinalOut
#define ScriptStreamOut             ScriptILINKStreamOut
#define CheckpointFilename          CheckpointILINKFilename
#define CheckpointFileBackup        CheckpointILINKFileBackup
#define MainOutfileDat              OutfILINKFinalDat
#define MainStreamDat               OutfILINKStreamDat
#endif  /* if defined(ILINK) */

/* end cgh */


/* We will be copying files block-by-block, so we use a constant for
   the block size.  Set this to the optimal value for your system. */

#define  CopyBlockSize  512

/* And now we define the two forms of copying: copyFile and appendFile.
   In each case, we are going to punt to ourCopyAppendFile, with the 
   appropriate flags being passed. */

#define  CopyOperation  0
#define  AppendOperation  1

#define  copyFile(f,t)		ourCopyAppendFile((f),(t),CopyOperation)
#define  appendFile(f,t)	ourCopyAppendFile((f),(t),AppendOperation)

#define  CopyAppendPerms  0666

/* bugfix for checkpoint code */
#ifndef STAT_SUCCESS
#define STAT_SUCCESS 0
#endif


/* We define some enumerated types that enable us to keep track of
   the type of checkpointing that is taking/took place, and of the
   calling sequence that lead to it */

typedef enum
{
  normalRun ,
  checkpointedRun
}
checkpointType ;

checkpointType  checkpointStatus ;

typedef enum
{
  iterationLocation ,
  functionLocation
}
checkpointLocation ;

struct ckptTupleType
{
  checkpointLocation  ckptLocation ;
  int                 ckptAttribute ;
  int                 iplace ;
  int                 jplace ;
} ;

struct ckptTupleType  ckptTuple ;

typedef enum
{
  iterped_call_before   = 0,         
  iterped_call_after    = 1,
  fCP_outf              = 2 ,
  fCP_gem_init          = 3 ,
  fCP_gem_iter_st_first = 4 ,
  fCP_gem_iter_st_decr  = 5 ,
  fCP_gem_iter_st_incr  = 6 ,
  fCP_gem_iter_gcen1    = 7 ,
  fCP_gem_iter_gfor     = 8 ,
  fCP_gem_gfor          = 9 ,
  fCP_gem_iter_gcen2    = 10 ,
  fCP_gem_ ,
  fCP_gem_iter_ ,
  fCP_gem_iter_gcen1_ ,
  fCP_gem_iter_gcen2_
}
funCallPathType ;

funCallPathType  funCallPath ;

/* This needs to be prototyped to get rid of warnings -- cgh */
#if !defined(KNR_PROTO)
void ourCopyAppendFile(char*, char*, int);
#else
void ourCopyAppendFile();
#endif  /* !defined(KNR_PROTO) */

/* These used to be in {ilink,mlink,linkmap}.c -- cgh */
int          statReturnCode ;
struct stat  statBuffer ;
char         dateTimeStamp [ DateTimeStampStringLength ] ;
FILE         * scriptCheckpoint ;
int          scriptRun , scriptToSleep ;
FILE         * fileTester;
int          iterpeds_counter;
int          checkpoint_counter;
int          checkpoint_place;

#endif  /* _CHECKPOINTDEFS_H */
