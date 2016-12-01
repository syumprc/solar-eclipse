/* 
 * solar.cc
 * Written by Charles P. Peterson
 * Copyright (c) 1997 Southwest Foundation for Biomedical Research
 *
 * This is the main routine for solar, a Tcl-interpreter based application
 *   for genetic analysis: Statistical Oligogenic Linkage Analysis Routines.
 *
 * Solar models consist of Solar/Tcl commands (i.e., they are scripts).
 *
 * Includes code derived from tclAppInit.c --
 * Copyright (c) 1993 The Regents of the University of California.
 * Copyright (c) 1994-1995 Sun Microsystems, Inc.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "solar.h"
#include "tclInt.h"

// globals from solarmain.c

extern "C" {
    int Batch_Mode;
    int Argc;
    char** Argv;
    extern char* StartingFile;
}


extern "C" void mp_setup ();

const char *SolarVersion ()
{
    return "8.1.4";
}

/*
 * Startup script run before first prompt appears
 */
const char *Internal_Startup = "\
history keep 200 \n\
set tcl_prompt1 {puts -nonewline \"solar> \"} \n\
proc quit {} {exit} \n\
set tcl_precision 17 \n\
set SOLAR_constraint_tol 0.0001 \n\
if {-1 == [lsearch [array names env] SOLAR_LIB]} { \n\
    set DEFAULT_LIB /usr/local/lib \n\
    puts stderr \"Warning.  SOLAR_LIB not defined.  Defaulting to $DEFAULT_LIB.\" \n\
} else { \n\
    set DEFAULT_LIB $env(SOLAR_LIB) \n\
} \n\
global auto_path \n\
if {[key -is-key-on-command-line?]} { \n\
    set auto_path [list . $DEFAULT_LIB $DEFAULT_LIB/tcl8.4] \n\
} else { \n\
    set auto_path [list . ~/lib $DEFAULT_LIB $DEFAULT_LIB/tcl8.4] \n\
} \n\
set found 0 \n\
# Messages are now appended to a global and printed later \n\
global Solar_Gotcl \n\
set Solar_Gotcl \"\" \n\
# Directory containing active solar.tcl is moved to head of the auto_path \n\
# this may change in future release \n\
# puts \"pre-sort auto_path is $auto_path\" \n\
foreach d $auto_path { \n\
  if {[file exists $d/solar.tcl]} { \n\
    set found 1 \n\
    if {\".\" != $d} { \n\
      set epos [lsearch -exact $auto_path $d] \n\
      set auto_path [lreplace $auto_path $epos $epos] \n\
      set auto_path [linsert $auto_path 0 $d] \n\
    } \n\
    if {![string compare $d .]} { \n\
      set Solar_Gotcl \"working directory\" \n\
    } elseif {![string compare $d ~/lib]} { \n\
      set Solar_Gotcl \"$env(HOME)/lib\" \n\
    } elseif {![file exists $d/tclIndex]} { \n\
      puts stderr \"Missing tclIndex for solar.tcl in $d\" \n\
      puts stderr \"Run solar as owner in that directory to build\" \n\
      puts stderr \"Or copy solar.tcl to current directory or ~/lib\" \n\
      exit \n\
    } \n\
    break \n\
  } \n\
} \n\
# puts \"final auto_path is now $auto_path\" \n\
if {!$found} { \n\
  puts stderr \"ERROR: solar.tcl not found in path: $auto_path\" \n\
  exit \n\
} \n\
# newtcl always looks in solar.tcl first, then in other files \n\
# this may change in future release \n\
proc newtcl {} { \n\
  if {[llength [glob -nocomplain ./*.tcl]]} { \n\
#   puts stderr \"Indexing tcl scripts in working directory\" \n\
    auto_mkindex . solar.tcl *.tcl \n\
  } \n\
  if {![key -is-key-on-command-line?] && \
       [llength [glob -nocomplain ~/lib/*.tcl]]} { \n\
#   puts stderr \"Indexing tcl scripts in ~/lib\" \n\
    auto_mkindex ~/lib solar.tcl *.tcl \n\
  } \n\
  auto_reset \n\
  return {} \n\
} \n\
newtcl \n\
set first_flag [llength [info globals finished_first_time]] \n\
global DEFAULT_LIB \n\
if [file exists .solar] { \n\
if {!$first_flag} {puts stderr \"Reading from .solar in working directory\"} \n\
source .solar \n\
} elseif [file exists ~/.solar] { \n\
if {!$first_flag} {puts stderr \"Reading from .solar in HOME directory\n\"} \n\
source ~/.solar \n\
} elseif [file exists $DEFAULT_LIB/.solar] { \n\
if {!$first_flag} {puts stderr \"Reading from .solar in $DEFAULT_LIB\n\"} \n\
source $DEFAULT_LIB/.solar \n\
} \n\
";

/*
 * Script run at startup and whenever a model is loaded
 */

const char *Model_Startup = "\
set finished_first_time {} \
";

// Use macro to declare command procedures (they all have same arguments)

#define DECL(name) \
extern "C" int name \
(ClientData clientData, Tcl_Interp* interp, int argc, const char *argv[]);

DECL(DefinitionCmd)
DECL(SolarFileCmd)
DECL(TableFileCmd)
DECL(ChiCmd)
DECL(NormalCmd)
DECL(KeyCmd)
DECL(HelpCmd)
DECL(FieldCmd)
DECL(PhenotypesCmd)
DECL(ParameterCmd)
DECL(ScaleCmd)
DECL(TraitCmd)
DECL(CovariateCmd)
DECL(OmegaCmd)
DECL(MatrixCmd)
DECL(MaximizeCmd)
DECL(LoglikeCmd)
DECL(QuadraticCmd)
DECL(ConstraintCmd)
DECL(OptionCmd)
DECL(MuCmd)
DECL(ModelCmd)
DECL(VerbosityCmd)
DECL(PedigreeCmd)
DECL(MarkerCmd)
DECL(FreqCmd)
DECL(MapCmd)
DECL(SnpCmd)
DECL(IbsCmd)
DECL(IbdCmd)
DECL(IbdOptCmd)
DECL(MibdCmd)
DECL(SimqtlCmd)
DECL(DrandCmd)
DECL(SolarBinaryVersionCmd)
DECL(SolarCompiledDateCmd)
DECL(SolarCompiledTimeCmd)
DECL(TclgrCmd)
DECL(AlnormCmd)
DECL(EVDCmd)
DECL(MaskCmd)
DECL(VoxelCmd)
DECL(ImoutCmd)
DECL(MathMatrixCmd)
DECL(TransposeCmd)
DECL(TimesCmd)
DECL(PlusCmd)
DECL(MinusCmd)
DECL(OlsCmd)
DECL(InverseCmd)
DECL(DinverseCmd)
DECL(EigenCmd)
DECL(MeanSumCmd)
DECL(SolveCmd)
DECL(MathMatrixShowCmd)
DECL(RowCmd);
DECL(ColCmd);
DECL(DiagonalCmd);
DECL(InsertCmd);
DECL(ShuffleCmd);
DECL(IdentityCmd);
DECL(MatrixPowerCmd);
DECL(ConcatenateCmd);
DECL(MathMatrixOutputCmd);
DECL(PermuteFCmd);
DECL(PermuteYCmd);
DECL(RunCreateFakePedigreeCmd);
DECL(RunSplitPhenoFileCmd);
DECL(RunPlinkConverter);
DECL(RunreorderphenoCmd);
int validate_solar (Tcl_Interp *interp);

extern void setup_functions ();

void add_solar_command (const char *name, Tcl_CmdProc proc, Tcl_Interp *interp)
{
    Solar_CreateCommand (interp, name, proc, (ClientData) NULL,
			 (Tcl_CmdDeleteProc *) NULL);
}

extern "C" int Solar_Init (Tcl_Interp *interp)
{
    TclRenameCommand (interp, (char*) "load", (char*) "loadbinary");
    add_solar_command ("define", DefinitionCmd, interp);
    add_solar_command ("solar_binary_version", SolarBinaryVersionCmd, interp);
    add_solar_command ("solar_compiled_date", SolarCompiledDateCmd, interp);
    add_solar_command ("solar_compiled_time", SolarCompiledTimeCmd, interp);
    add_solar_command ("tablefile", TableFileCmd, interp);
    add_solar_command ("solarfile", SolarFileCmd, interp);
    add_solar_command ("tclgr", TclgrCmd, interp);
    add_solar_command ("chi", ChiCmd, interp);
    add_solar_command ("normal", NormalCmd, interp);
    add_solar_command ("key", KeyCmd, interp);
    add_solar_command ("verbosity", VerbosityCmd, interp);
    add_solar_command ("help", HelpCmd, interp);
    add_solar_command ("field", FieldCmd, interp);
    add_solar_command ("cmibd", MibdCmd, interp);
    add_solar_command ("ibdoption", IbdOptCmd, interp);
    add_solar_command ("cibd", IbdCmd, interp);
    add_solar_command ("cibs", IbsCmd, interp);
    add_solar_command ("drand", DrandCmd, interp);
    add_solar_command ("csimqtl", SimqtlCmd, interp);
    add_solar_command ("csnp", SnpCmd, interp);
    add_solar_command ("cmap", MapCmd, interp);
    add_solar_command ("cfreq", FreqCmd, interp);
    add_solar_command ("cmarker", MarkerCmd, interp);
    add_solar_command ("cpedigree", PedigreeCmd, interp);
    add_solar_command ("loglike", LoglikeCmd, interp);
    add_solar_command ("quadratic", QuadraticCmd, interp);
    add_solar_command ("ccmaximize", MaximizeCmd, interp);
    add_solar_command ("omega", OmegaCmd, interp);
    add_solar_command ("matrix", MatrixCmd, interp);
    add_solar_command ("ctranspose",TransposeCmd, interp);
    add_solar_command ("times", TimesCmd, interp);
    add_solar_command ("plus", PlusCmd, interp);
    add_solar_command ("minus", MinusCmd, interp);
    add_solar_command ("inverse", InverseCmd, interp);
    add_solar_command ("dinverse", DinverseCmd, interp);
    add_solar_command ("parameter", ParameterCmd, interp);
    add_solar_command ("scale", ScaleCmd, interp);
    add_solar_command ("covariate", CovariateCmd, interp);
    add_solar_command ("constraint", ConstraintCmd, interp);
    add_solar_command ("trait", TraitCmd, interp);
    add_solar_command ("phenotypes", PhenotypesCmd, interp);
    add_solar_command ("option", OptionCmd, interp);
    add_solar_command ("mu", MuCmd, interp);
//    add_solar_command ("save", LoadSaveCmd, interp);
//    add_solar_command ("load", LoadSaveCmd, interp);
    add_solar_command ("model", ModelCmd, interp);
    add_solar_command ("alnorm", AlnormCmd, interp);
    add_solar_command ("evd", EVDCmd, interp);
    add_solar_command ("voxel", VoxelCmd, interp);
    add_solar_command ("mask", MaskCmd, interp);
    add_solar_command ("imout", ImoutCmd, interp);
    add_solar_command ("mathmatrix",MathMatrixCmd, interp);
    add_solar_command ("rows",MathMatrixCmd, interp);
    add_solar_command ("cols",MathMatrixCmd, interp);
    add_solar_command ("peek",MathMatrixCmd, interp);
    add_solar_command ("ols",OlsCmd, interp);
    add_solar_command ("evalues",EigenCmd, interp);
    add_solar_command ("evectors",EigenCmd, interp);
    add_solar_command ("mean",MeanSumCmd, interp);
    add_solar_command ("min",MeanSumCmd, interp);
    add_solar_command ("max",MeanSumCmd, interp);
    add_solar_command ("solve",SolveCmd, interp);
    add_solar_command ("show",MathMatrixShowCmd, interp);
    add_solar_command ("row",RowCmd, interp);
    add_solar_command ("col",ColCmd, interp);
    add_solar_command ("diagonal",DiagonalCmd, interp);
    add_solar_command ("insert",InsertCmd, interp);
    add_solar_command ("shuffle",ShuffleCmd, interp);
    add_solar_command ("identity",IdentityCmd, interp);
    add_solar_command ("matrixpower",MatrixPowerCmd, interp);
    add_solar_command ("concatenate",ConcatenateCmd, interp);
    add_solar_command ("output", MathMatrixOutputCmd, interp);
    add_solar_command ("permutef", PermuteFCmd, interp);
    add_solar_command ("permutey", PermuteYCmd, interp);
    add_solar_command ("create_fake_pedigree", RunCreateFakePedigreeCmd, interp);
    add_solar_command ("split_class_file", RunSplitPhenoFileCmd, interp);
    add_solar_command ("plink_converter", RunPlinkConverter, interp);
    add_solar_command ("reorder_phenotype", RunreorderphenoCmd, interp);

//  internal initializations (not subject to error)

    setup_functions ();
    Option::setup ();

//  Copyright message

    if (Batch_Mode)
    {
	Solar_Eval (interp, "global Solar_Batch ; set Solar_Batch 1");
    }
    else
    {
	Solar_Eval (interp, "global Solar_Batch ; set Solar_Batch 0");
    }

//  Initializations (subject to error)

    if ( TCL_OK != Solar_Eval (interp, Internal_Startup))
    {
	if (StartingFile) unlink (StartingFile);
	fprintf (stderr, "Error in tcl script file:\n");
	fprintf (stderr, "%s\n", Tcl_GetStringResult (interp));
	exit (EXIT_FAILURE);
    }
    if (StartingFile) unlink (StartingFile);

    if ( TCL_OK != Solar_Eval (interp, "solar_tcl_startup"))
    {
	fprintf (stderr, "Error in solar_tcl_startup:\n");
	fprintf (stderr, "%s\n", Tcl_GetStringResult (interp));
	exit (EXIT_FAILURE);
    }

    if ( TCL_OK != Solar_Eval (interp, "make_solar_aliases"))
    {
	fprintf (stderr, "Error making default aliases: \n");
	fprintf (stderr, "%s\n", Tcl_GetStringResult (interp));
	exit (EXIT_FAILURE);
    }

    Phenotypes::start();  // Loads phenotypes.info
    Field::Start();      // Loads field.info

    Solar_Eval (interp, "Start_Dirs");

    

    if ( TCL_OK != Model::renew(interp))
    {
	fprintf (stderr, "%s\n", Tcl_GetStringResult (interp));
	exit (EXIT_FAILURE);
    }

    if ( TCL_OK != Solar_Eval (interp, "fatal_error_checks") )
    {
	if (!Strcmp (Tcl_GetStringResult (interp), "invalid command name \"fatal_error_checks\""))
	{
	    fprintf (stderr, 
"Obsolete version of solar.tcl in path: . ~/lib $SOLAR_LIB\n");
	}
	else
	{
	    fprintf (stderr, "%s\n", Tcl_GetStringResult (interp));
	}
	exit (EXIT_FAILURE);
    }

    Solar_Eval (interp, "check_os");
    Tcl_ResetResult (interp);
	
    validate_solar (interp);

    if (Argc > 1)
    {
	char *command_string = Tcl_Concat (Argc-1, &Argv[1]);
//	fprintf (stderr, "SOLAR executing command:  %s\n", command_string);
	int code = Solar_Eval (interp, command_string);
	Tcl_Free (command_string);
	if (*interp->result != 0)
	{
	    printf ("%s\n", interp->result);
	}
	if (code != TCL_OK)
	{
	    exit (EXIT_FAILURE);
	}
	exit (EXIT_SUCCESS);
    }

    return TCL_OK;
}

int SolarBinaryVersionCmd (ClientData clientData, Tcl_Interp *interp,
		     int argc, const char *argv[])
{
	
    char buf[256];
    sprintf (buf, "%s", SolarVersion());
    RESULT_BUF (buf);
    return TCL_OK;
}

int SolarCompiledDateCmd (ClientData clientData, Tcl_Interp *interp,
		     int argc, const char *argv[])
{
    RESULT_LIT (__DATE__);
    return TCL_OK;
}

int SolarCompiledTimeCmd (ClientData clientData, Tcl_Interp *interp,
		     int argc, const char *argv[])
{
    RESULT_LIT (__TIME__);
    return TCL_OK;
}


extern "C" void error (const char *message)
{
    fprintf (stderr, "solar: %s\n", message);
    exit (EXIT_FAILURE);
}



