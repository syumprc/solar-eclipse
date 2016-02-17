/*
 * loadsave.cc implements the load and save commands
 * Written by Charles Peterson beginning on January 19, 1998
 * Copyright (c) 1997 Southwest Foundation for Biomedical Research
 */

// As of Version 8.0.1, this is now obsolete.  It is no longer added to
//   the interpreter so this routine is no longer called.
//
// It was originally created as a "fast" way of making it possible
//   for the "load" and "save" operators of many commands to come first,
//   for an imperative voice that seems natural.  Or in other words "load"
//   is a virtual command overloaded for each type of thing that can be loaded.
//
// However until 8.0.1 I didn't realize or entirely understand the serious flaw
//   this has.  The following C++ code is only reached after Tcl itself has
//   "evaluated" the command and arguments.  Then, from here, we evaluate the
//   command and arguments again.
//
// The consequence of this double evaluation is that arguments that are lists
//   or have spaces inside them (such as filenames or directory names with
//   spaces in them) get divided as separate arguments.
//
// The easiest way to solve this, and I now judge it to have immaterial
//    consequences for speed, is to implement load and save in Tcl.  Within
//    Tcl the list of remaining arguments can be a variable expanded in the
//    second evaluation, which thereby becomes harmless.  Duplicating this
//    behavior in C++ is not so easy.
//
//  Wrt speed, the only consequence really is having to load the "load" command
//    from solar.tcl, and this is only done once.  The command itself is
//    compiled into Tcl intermediate code and executes about as fast as this
//    C++ would if it were redesigned to handle the list of arguments correctly.

#include <stdlib.h>
#include <stdio.h>
// tcl.h from solar.h
#include "solar.h"

extern "C" int LoadSaveCmd (ClientData clientData, Tcl_Interp *interp,
		  int argc, char *argv[])
{
    if (argc <= 2)
    {
	if (!(argc == 2 && !StringCmp (argv[1], "help", case_ins)))
	{
	    fprintf (stderr, "Invalid %s command\n\n", argv[0]);
	}
	printf (
"Purpose:  The load command is an alias to make some other commands nicer.\n");
	printf (
"          For example, the command \"load pedigree\" is actually an alias for\n");
	printf (
"          \"pedigree load.\"  For more information about a particular \"load\"\n");
	printf (
"          command, see the help for the underlying command, for example:\n");
	printf (
"          \"help pedigree\".  Here is a summary of some load commands:\n");
	printf ("\n");
	printf (
"              load pedigree <filename>\n");
	printf (
"              load phenotypes <filename>\n");
	printf (
"              load matrix <filename> <name1> [<name2>]\n");
	printf (
"              load model <filename>\n");
	printf (
"              load freq <filename>\n");
	printf (
"              load marker <filename>\n");
	printf (
"              load map <filename>\n");


	if (argc == 2 && !StringCmp (argv[1], "help", case_ins)) return TCL_OK;
	return TCL_ERROR;
    }

    char buf[1024];
    char *arguments = Tcl_Concat (argc-2, argv+2);
    sprintf (buf, "%s %s %s", argv[1], argv[0], arguments);
    Tcl_Free (arguments);
    return Solar_Eval (interp, buf);
}
