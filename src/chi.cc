/*
 * chi.cc implements the Chi command
 * Written by Charles Peterson beginning on March 10, 1998
 * Copyright (c) 1997 Southwest Foundation for Biomedical Research
 */

#include "safelib.h"
// tcl.h from solar.h
#include "solar.h"

extern "C" void cdfchi_ (int*, double*, double*, double*, double*,
			 int*, double*);

extern "C" int ChiCmd (ClientData clientData, Tcl_Interp *interp,
	    int argc, char *argv[])
{
    if (argc == 2 && !StringCmp ("help", argv[1], case_ins))
    {
	return Solar_Eval (interp, "help chi");
    }
    if (!Strcmp (argv[1], "-inverse"))
    {
	if (argc != 4)
	{
	    RESULT_LIT (
		"chi -inverse option requires 2 arguments: p and d.f.");
	    return TCL_ERROR;
	}
	double p, q, chival, df, bound;
	int status = 0;
	int which = 2;
	char dummy[512];
	if (1 == sscanf (argv[2], "%lf %s", &q, &dummy))
	{
	    p = 1 - q;
	    if (1 == sscanf (argv[3], "%lf %s", &df, &dummy))
	    {
		cdfchi_ (&which, &p, &q, &chival, &df, &status, &bound);
		if (status == 0)
		{
		    char buf[256];
		    sprintf (buf, "%f", chival);
		    RESULT_BUF (buf);
		    return TCL_OK;
		} 
		else if (status == -1)
		{
		    RESULT_LIT ("Internal error (invalid operation)");
		    return TCL_ERROR;
		}
		else if (status == -2)
		{
		    RESULT_LIT ("The p-value must be in the range [0,1]");
		    return TCL_ERROR;
		}
		else if (status == -3)
		{
		    RESULT_LIT ("The p-value must be in the range [0,1]");
		    return TCL_ERROR;
		}
		else if (status == -4)
		{
		    RESULT_LIT ("The chi-square value must be >= 0");
		    return TCL_ERROR;
		}
		else if (status == -5)
		{
		    RESULT_LIT ("The degrees of freedom must be > 0");
		    return TCL_ERROR;
		}
		else if (status == 1)
		{
		    RESULT_LIT \
			("Internal error (answer lower than search bound)");
		    return TCL_ERROR;
		}
		else if (status == 2)
		{
		    RESULT_LIT \
			("Internal error (answer higher than search bound)");
		    return TCL_ERROR;
		}
		else if (status == 3)
		{
		    RESULT_LIT ("Internal error (p + q != 1)");
		    return TCL_ERROR;
		}
		else if (status == 10)
		{
		    RESULT_LIT ("Internal error (in cdfgam routine)");
		    return TCL_ERROR;
		}
		RESULT_LIT ("Error calculating chi-square from p-value");
		return TCL_ERROR;
	    }
	}
	RESULT_LIT ( "Invalid chi command");
	return TCL_ERROR;
    }
    if (argc == 3 || argc == 4)
    {
	int cvi = 1;
	int dfi = 2;
	char prefix[6];
	strcpy (prefix, "p = ");
	if (argc == 4)
	{
	    if (Strcmp (argv[1], "-number"))
	    {
		RESULT_LIT ("Invalid chi command argument");
		return TCL_ERROR;
	    }
	    cvi = 2;
	    dfi = 3;
	    strcpy (prefix, "");
	}
	double chival;
	double df;
	char dummy[512];
	if (1 == sscanf (argv[cvi], "%lf %s", &chival, &dummy))
	{
	    if (1 == sscanf (argv[dfi], "%lf %s", &df, &dummy))
	    {
		double p, q, bound;
		int status = 0;
		int which = 1;
		cdfchi_ (&which, &p, &q, &chival, &df, &status, &bound);
		if (status == 0)
		{
		    double prob = q;
/*		    if (prob >= 0.0000001) */
		    if (prob >= 0.00001)
		    {
			char buf[256];
			sprintf (buf, "%s%-9.7f", prefix, prob);
			RESULT_BUF (buf);
		    }
		    else if (prob >= 1e-307)
		    {
			char buf[256];
			sprintf (buf, "%s%.8g", prefix, prob);
			RESULT_BUF (buf);
		    }
		    else
		    {
			if (strlen (prefix))
			{
			    strcpy (prefix, "p < ");
			}
			char buf[256];
			sprintf (buf, "%s1e-307", prefix);
			RESULT_BUF (buf);
		    }
		    return TCL_OK;
		}
		else if (status == -1)
		{
		    RESULT_LIT ("Internal error (invalid operation)");
		    return TCL_ERROR;
		}
		else if (status == -2)
		{
		    RESULT_LIT ("The p-value must be in the range [0,1]");
		    return TCL_ERROR;
		}
		else if (status == -3)
		{
		    RESULT_LIT ("The p-value must be in the range [0,1]");
		    return TCL_ERROR;
		}
		else if (status == -4)
		{
		    RESULT_LIT ("The chi-square value must be >= 0");
		    return TCL_ERROR;
		}
		else if (status == -5)
		{
		    RESULT_LIT ("The degrees of freedom must be > 0");
		    return TCL_ERROR;
		}
		else if (status = 1)
		{
		    RESULT_LIT \
			("Internal error (answer lower than search bound)");
		    return TCL_ERROR;
		}
		else if (status = 2)
		{
		    RESULT_LIT \
			("Internal error (answer higher than search bound)");
		    return TCL_ERROR;
		}
		else if (status = 3)
		{
		    RESULT_LIT ("Internal error (p + q != 1)");
		    return TCL_ERROR;
		}
		else if (status = 10)
		{
		    RESULT_LIT ("Internal error (in cdfgam routine)");
		    return TCL_ERROR;
		}
		RESULT_LIT ("Error calculating p-value for chi-square");
		return TCL_ERROR;
	    }
	}
    }
    RESULT_LIT ( "Invalid chi command");
    return TCL_ERROR;
}


