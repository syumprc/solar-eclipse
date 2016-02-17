/*
 * loglike.cc implements the loglike command
 * Written by Charles Peterson beginning on December 3, 1997
 * Copyright (c) 1997 Southwest Foundation for Biomedical Research
 */

#include <stdio.h>
#include "solar.h"
// tcl.h from solar.h

model_status Loglike::_status = unoptimized;
double Loglike::_loglike = 0.0;
double Loglike::_quadratic = 0.0;

extern "C" int LoglikeCmd (ClientData clientData, Tcl_Interp *interp,
			   int argc, char *argv[])
{
    if (argc == 3 && !StringCmp (argv[1], "set", case_ins))
    {
	char junk[256];
	double newloglike;
	if (1 != sscanf (argv[2], "%lf %s", &newloglike, &junk))
	{
	    RESULT_LIT ("Invalid number for loglike");
	    return TCL_ERROR;
	}
	Loglike::loglike (newloglike);
	return TCL_OK;
    }

    if (argc != 1)
    {
	Solar_Eval (interp, "help loglike");
	if (!StringCmp (argv[1], "help", case_ins)) return TCL_ERROR;
	return TCL_OK;
    }
    try
    {
	char buf[128];
	sprintf (buf, loglike_format, Loglike::loglike() );
	RESULT_BUF (buf);
	return TCL_OK;
    }
    catch (Safe_Error_Return& ser)
    {
	RESULT_BUF (ser.message ());
	return TCL_ERROR;
    }
}


extern "C" int QuadraticCmd (ClientData clientData, Tcl_Interp *interp,
			     int argc, char *argv[])
{
    if (argc == 3 && !StringCmp (argv[1], "set", case_ins))
    {
	char junk[256];
	double newquadratic;
	if (1 != sscanf (argv[2], "%lf %s", &newquadratic, &junk))
	{
	    RESULT_LIT ("Invalid number for quadratic");
	    return TCL_ERROR;
	}
	Loglike::quadratic (newquadratic);
	return TCL_OK;
    }

    if (argc != 1)
    {
	Solar_Eval (interp, "help quadratic");
	if (!Strcmp (argv[1], "help")) return TCL_ERROR;
	return TCL_OK;
    }
    try
    {
	char buf[128];
	sprintf (buf, loglike_format, Loglike::quadratic() );
	RESULT_BUF (buf);
	return TCL_OK;
    }
    catch (Safe_Error_Return& ser)
    {
	RESULT_BUF (ser.message ());
	return TCL_ERROR;
    }
}


char *Loglike::get_string (char *buf)
{
    sprintf (buf, loglike_format, loglike() );
    return buf;
}

char *Loglike::get_quadratic_string (char *buf)
{
    sprintf (buf, loglike_format, quadratic() );
    return buf;
}

void Loglike::check_status ()
{
    if (_status == valid)
    {
	return;
    }
    else if (_status == unoptimized)
    {
	throw Safe_Error_Return ("ERROR: Model must be maximized first");
    }

// _status == failed

    throw Safe_Error_Return ("ERROR: Last model maximization failed");
}

