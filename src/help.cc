/*
 * omega.cc implements the help command
 * Written by Charles Peterson beginning on October 27, 1997
 * Copyright (c) 1997 Southwest Foundation for Biomedical Research
 */

#include <string.h>
#include <stdlib.h>
// tcl.h from solar.h
#include "solar.h"

extern "C" int HelpCmd (ClientData clientData, Tcl_Interp *interp,
		  int argc, char *argv[])
{
    if (argc == 1 || argc == 3 && (!Strcmp (argv[1], "-output") ||
	   !Strcmp (argv[1], "-out") || !Strcmp (argv[1], "-o")))
    {
	if (TCL_OK == Solar_Eval (interp, "help-all-commands"))
	{
	    char buf[256];
	    char filename[256];
	    strncpy (filename, Tcl_GetStringResult (interp), 256);
	    if (argc == 1)
	    {
		sprintf (buf, "exec >&@stdout <@stdin more %s", filename);
		Solar_Eval (interp, buf);
		unlink (filename);
	    }
	    else
	    {
		sprintf (buf, "mv %s %s", Tcl_GetStringResult (interp), argv[2]);
		system (buf);
		Tcl_ResetResult (interp);
	    }
	    return TCL_OK;
	}
	return TCL_ERROR;
    }
    if (argc == 2 || argc == 4 && (!Strcmp (argv[2], "-output") ||
	   !Strcmp (argv[2], "-output") || !Strcmp (argv[2], "-output")))
    {
	char buf[1024];
	sprintf (buf, "help-on-command %s", argv[1]);
	if (TCL_OK == Solar_Eval (interp, buf)) {
	    char filename[256];
	    strncpy (filename, Tcl_GetStringResult (interp), 256);
	    sprintf (buf, "exec >&@stdout <@stdin more %s", filename);
	    if (argc == 2)
	    {
		sprintf (buf, "exec >&@stdout <@stdin more %s", filename);
		Solar_Eval (interp, buf);
		unlink (filename);
	    }
	    else
	    {
		sprintf (buf, "mv %s %s", Tcl_GetStringResult (interp), argv[3]);
		system (buf);
		Tcl_ResetResult (interp);
	    }
	    return TCL_OK;
	}
	return TCL_ERROR;
    }
    RESULT_LIT ("Invalid help command\n");
    return TCL_ERROR;
}

 

    




