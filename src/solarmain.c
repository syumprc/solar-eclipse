/* 
 * solar.cc
 * Written by Charles P. Peterson
 * Copyright (c) 1999 Southwest Foundation for Biomedical Research
 *
 * This is the (now nearly vestigial) main routine for solar, which needs
 * to be compiled as C (not C++) for full compatibility with Tcl (now that
 # "C" is a compiler type.
 *
 * Includes code derived from tclAppInit.c --
 * Copyright (c) 1993 The Regents of the University of California.
 * Copyright (c) 1994-1995 Sun Microsystems, Inc.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "tcl.h"

/* This is really the only thing in SOLAR referenced here */
int Solar_Init (Tcl_Interp *interp);


/*
 * The following variable is a special hack that is needed in order for
 * Sun shared libraries to be used for Tcl (copied from Tcl source).
 */
#ifdef DUMATHPTR
extern int matherr();
int *tclDummyMathPtr = (int *) matherr;
#endif


#ifdef BREAK_ON_NANS
#include <sunmath.h>
#endif

int Batch_Mode;
int Condor;
int NisKey;
char Key[10];

int Argc;
char** Argv;
char* StartingFile = 0;

int main(int argc, char *argv[])
{
    char **argvi = argv;
    NisKey = 0;
/*
 * version 6.4.3 starts with .SOLARstarting filename that
 * must be deleted asap in solar.cc
 * then all following arguments, if any, are moved forwards
 */
    if (argc > 1 && !strncmp(argvi[1],".SOLARstarting", 14))
    {
	StartingFile = argvi[1];
	argvi++;
	argc--;
    }

/*
 * version 6.4.4 and greater allow -noce argument to suppress command line
 * editing.  That is handled in startup script, here it must simply be
 * ignored.
 */

    if (argc > 1 && !strcmp (argvi[1],"-noce"))
    {
	argvi++;
	argc--;
    }

/*
 * version 6.5.1 and greater permit -condor argument which enables user
 * identification on Condor parallel systems
 */
    if (argc > 1 && !strcmp (argvi[1],"-condor"))
    {
	Condor = 1;
	argvi++;
	argc--;
    }
    else
    {
	Condor = 0;
    }

#ifdef BREAK_ON_NANS
    ieee_handler ("set", "common", SIGFPE_ABORT);
#endif

    Key[0] = '\0';
    if (argc > 2 && !strcmp (argvi[1],"-key"))
    {
	strncpy (Key, argvi[2], 8);
	Key[9] = '\0';
	argc -= 2;
	argvi = argv+2;
    }
    else if (argc > 1 && !strcmp (argvi[1],"-niskey"))
    {
	if (argc==2)
	{
	    fprintf (stderr, "-niskey requires argument\n");
	    return 1;
	}
	strncpy (Key, argvi[2], 8);
	Key[9] = '\0';
	argc -= 2;
	argvi = argv+2;
	NisKey = 1;
    }

    Argc = argc;
    Argv = argvi;

    if (argc > 1) {
	Batch_Mode = 1;
    } else {
	Batch_Mode = 0;
    }
    Tcl_Main(1, argvi, Tcl_AppInit);
    return 0;			/* Needed only to prevent compiler warning. */
}

int Tcl_AppInit (Tcl_Interp *interp)
{
    if (Tcl_Init(interp) == TCL_ERROR) {
	return TCL_ERROR;
    }
    return Solar_Init (interp);
}

