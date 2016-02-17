/*
 * tclgr.c
 * tcl command interface for ACE/gr plotting package
 * Utilizes acegr_np (named pipe interface) provided with ACE/gr 4.1.1
 *
 * Written by Charles Peterson, 11 May 1998
 * Copyright (C) 1998 Southwest Foundation for Biomedical Research
 *
 * acegr_np is Copyright (C) 1997 Henrik Seidel
 * acegr is Copyright (c) 1991-95 Paul J Turner, Portland, OR
 * AND Copyright (c) 1996-97 ACE/gr Development Team
 *
 */

#include <string.h>
#include "tcl.h"
#include "plotpipe.h"

char gr_command[1024] = "xmgr";
static int terminated_by_user = 0;

int TclgrCmd (ClientData clientData, Tcl_Interp *interp, int argc, 
	      const char *argv[]);

int Tclgr_Init (Tcl_Interp *interp)
{
    Tcl_CreateCommand (interp, "tclgr", TclgrCmd, (ClientData) 0,
		       (Tcl_CmdDeleteProc *) 0);
    return TCL_OK;
}

int TclgrCmd (ClientData clientData, Tcl_Interp *interp, int argc, 
	      const char *argv[])
{
    if (argc == 1 || (argc == 2 && (!strcmp (argv[1], "help") ||
	!strcmp (argv[1], "usage"))))
    {
	return Tcl_Eval (interp, "help tclgr");
    }

    if (!strcmp (argv[1], "syscommand"))
    {
	if (argc < 3)
	{
	    interp->result = gr_command;
	    return TCL_OK;
	}
	else
	{
	    strcpy (gr_command,argv[2]);
	    return TCL_OK;
	}
    }

    if (!strcmp (argv[1], "open"))
    {
	int status;
	int bufsize = 1000;
	int argindex;

	for (argindex = 2; argindex < argc; argindex++)
	{
	    if (!strcmp (argv[argindex], "-buffersize"))
	    {
		argindex++;
		if (argindex < argc)
		{
		    char dummy[512];
		    int count = sscanf (argv[argindex], 
					"%d %s", &bufsize, &dummy);
		    if (count == 1) 
		    {
			continue;
		    }
		}
		interp->result = "Ill-formed -buffersize argument";
		return TCL_ERROR;
	    }
	    else
	    {
		interp->result = "Invalid tclgr open arguments";
		return TCL_ERROR;
	    }
	}
	if (terminated_by_user)
	{
	    interp->result = 
    "ACE/gr session apparently terminated by user; use 'tclgr close' to reset";
	    return TCL_ERROR;
	}
	status = ACEgrOpen (bufsize, gr_command);
	if (status)
	{
	    if (status == -2)
	    {
		interp->result = 
		    "tclgr session already opened from this solar session";
	    }
	    else
	    {
		interp->result = "Error returned from ACEgrOpen";
	    }
	    return TCL_ERROR;
	}
	return TCL_OK;
    }
    if (!strcmp (argv[1], "close"))
    {
	int status = ACEgrClose (terminated_by_user);
	terminated_by_user = 0;
	if (status)
	{
	    interp->result = "Error returned from ACEgrClose";
	    return TCL_ERROR;
	}
	return TCL_OK;
    }
    if (terminated_by_user)
    {
	interp->result = 
    "ACE/gr session apparently terminated by user; use 'tclgr close' to reset";
	return TCL_ERROR;
    }
    if (!strcmp (argv[1], "buffer"))
    {
	char *command = Tcl_Concat (argc-2, &argv[2]);
	int status = ACEgrCommand (command);
	if (status)
	{
	    if (status == -2)
	    {
		terminated_by_user = 1;
		interp->result = 
    "ACE/gr session apparently terminated by user; use 'tclgr close' to reset";
		return TCL_ERROR;
	    }
	    interp->result = "Error returned from ACEgrCommand";
	    return TCL_ERROR;
	}
	return TCL_OK;
    }

    if (!strcmp (argv[1], "send"))
    {
	char *command = Tcl_Concat (argc-2, &argv[2]);
	int status = ACEgrCommand (command);
	if (status)
	{
	    interp->result = "Error returned from ACEgrCommand";
	    return TCL_ERROR;
	}
	status = ACEgrFlush ();
	if (status)
	{
	    if (status == -2)
	    {
		terminated_by_user = 1;
		interp->result = 
    "ACE/gr session apparently terminated by user; use 'tclgr close' to reset";
		return TCL_ERROR;
	    }
	    interp->result = "Error returned from ACEgrFlush";
	    return TCL_ERROR;
	}
	return TCL_OK;
    }
	
    if (!strcmp (argv[1], "flush"))
    {
	int status = ACEgrFlush ();
	if (status)
	{
	    if (status == -2)
	    {
		terminated_by_user = 1;
		interp->result = 
    "ACE/gr session apparently terminated by user; use 'tclgr close' to reset";
		return TCL_ERROR;
	    }
	    interp->result = "Error returned from ACEgrFlush";
	    return TCL_ERROR;
	}
	return TCL_OK;
    }

    interp->result = "Unrecognized tclgr command; see tclgr help";
    return TCL_ERROR;
}



