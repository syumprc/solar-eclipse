/*
 * normal.cc implements the normal command
 * interface to cdfnor.f
 *
 * Written by Charles Peterson beginning on May 17, 2005
 * Copyright (c) 2005 Southwest Foundation for Biomedical Research
 */

#include "solar.h"

extern "C" void cdfnor_ (int* which, double* p, double* q, double* x, double* mean, double* sd, int* status, double* bound);

extern "C" int NormalCmd (ClientData clientData, Tcl_Interp* interp,
			  int argc, char* argv[])
{
    char buf[1024];
/*
 * normal -inverse
 */
    if (argc==3 && (!Strcmp ("-inverse", argv[1]) ||
		    !Strcmp ("-i", argv[1]) ||
		    !Strcmp ("-in", argv[1]) ||
		    !Strcmp ("-inv", argv[1]) ||
		    !Strcmp ("-inve", argv[1]) ||
		    !Strcmp ("-inver", argv[1]) ||
		    !Strcmp ("-invers", argv[1])))
    {
	double pct;

	int con = sscanf (argv[2], "%lg %1c", &pct, buf);
	if (con != 1)
	{
	    printf ("con is %d\n", con);
	    sprintf (buf, "normal: invalid number %s\n", argv[2]);
	    RESULT_BUF (buf);
	    return TCL_ERROR;
	}

	int which = 2;
	double q = 1.0 - pct;
	double mean = 0;
	double sd = 1;
	int status = 0;
	double bound = 0;
	double x = 0;

//	printf ("Calling cdfnor\n");
	cdfnor_ (&which, &pct, &q, &x, &mean, &sd, &status, &bound);
//	printf ("Called cdfnor\n");
	if (status != 0)
	{
	    sprintf (buf, "normal: Inverse range error for %g", pct);
	    RESULT_BUF (buf);
	    return TCL_ERROR;
	}

	sprintf (buf, "%14.9f", x);
	double ftest = 0;
	sscanf (buf, "%g", &ftest);
	if (ftest == 0)
	{
	    sprintf (buf, "%14.9g", x);
	}
	RESULT_BUF (buf);
	return TCL_OK;
    }
    RESULT_LIT ("normal: Invalid arguments");
    return TCL_ERROR;
}

