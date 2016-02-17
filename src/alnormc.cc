/*
 * alnormc.cc provides a Tcl interface to the alnorm subroutine
 * Written by Charles Peterson beginning on August 10, 2000
 * Copyright (c) 2000 Southwest Foundation for Biomedical Research
 */

#include "solar.h"

extern "C" int AlnormCmd (ClientData clientData, Tcl_Interp *interp,
			  int argc, char *argv[])
{
    double z, alnorm_return, alnorm_ (double*,int*);
    int upper;

    if (argc == 2 && !Strcmp ("help", argv[1]))
    {
	return Solar_Eval (interp, "help alnorm");
    }

    if (argc != 3)
    {
	RESULT_LIT ("Wrong number of arguments to alnorm");
	return TCL_ERROR;
    }

    char *dbuf = (char*) Malloc (strlen (argv[1]) + 1);
    if (1 != sscanf (argv[1], "%lf%s", &z, dbuf))
    {
	free (dbuf);
	RESULT_LIT ("First argument to alnorm must be a number");
	return TCL_ERROR;
    }
    free (dbuf);

    if (!Strcmp (argv[2], "t") || !Strcmp (argv[2], "1"))
    {
	upper = 1;
    }
    else if (!Strcmp (argv[2], "f") || !Strcmp (argv[2], "0"))
    {
	upper = 0;
    }
    else
    {
	RESULT_LIT ("Second argument to alnorm must be t or f");
	return TCL_ERROR;
    }

    alnorm_return = alnorm_ (&z, &upper);

    char buf[128];
    sprintf (buf, "%-10g", alnorm_return);
    RESULT_BUF (buf);
    return TCL_OK;
}

