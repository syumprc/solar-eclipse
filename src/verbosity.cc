/*
 * verbosity.cc implements the verbosity command
 * Written by Charles Peterson beginning on January 16, 1998
 * Copyright (c) 1997 Southwest Foundation for Biomedical Research
 */

#include "solar.h"

unsigned int Verbosity::_verbosity = 255;

extern "C" int VerbosityCmd (ClientData clientData, Tcl_Interp *interp,
		  int argc, char *argv[])
{
    if (argc == 2 && !StringCmp ("help", argv[1], case_ins))
    {
	return Solar_Eval (interp, "help verbosity");
    }

    if (argc == 1)
    {
	if (Verbosity::max())
	{
	    RESULT_LIT ("verbosity max");
	    return TCL_OK;
	}
	if (Verbosity::plus())
	{
	    RESULT_LIT ("verbosity plus");
	    return TCL_OK;
	}
	if (Verbosity::def())
	{
	    RESULT_LIT ("verbosity default");
	    return TCL_OK;
	}
	if (Verbosity::min())
	{
	    RESULT_LIT ("verbosity min");
	    return TCL_OK;
	}
	char buf[256];
	sprintf (buf, "verbosity 0x%0x", Verbosity::get());
	RESULT_BUF (buf);
	return TCL_OK;
    }

    if (argc == 2)
    {
	if (!StringCmp (argv[1], "-number", case_ins))
	{
	    int vlevel = Verbosity::get ();
	    char buf[256];
	    sprintf (buf, "0x%0x", vlevel);
	    RESULT_BUF (buf);
	    return TCL_OK;
	}
	if (!StringCmp (argv[1], "default", case_ins))
	{
	    Verbosity::set (default_verbosity);
	    return TCL_OK;
	}
	if (!StringCmp (argv[1], "plus", case_ins))
	{
	    Verbosity::set (plus_verbosity);
	    return TCL_OK;
	}
	if (!StringCmp (argv[1], "min", case_ins))
	{
	    Verbosity::set (min_verbosity);
	    return TCL_OK;
	}
	if (!StringCmp (argv[1], "max", case_ins))
	{
	    Verbosity::set (max_verbosity);
	    return TCL_OK;
	}
	char dummy[256];
	int new_verbosity = 0;
	if (1 != sscanf (argv[1], "%i %s", &new_verbosity, &dummy) ||
	    new_verbosity < 0)
	{
	    RESULT_LIT ("Invalid verbosity specification");
	    return TCL_ERROR;
	}
	Verbosity::set (new_verbosity);
	return TCL_OK;
    }
    RESULT_LIT ("Invalid verbosity command");
    return TCL_ERROR;
}

extern "C" int verbose_ (char *keyword, int len_keyword_)
{
    char buf[32];
    int len = (len_keyword_ < 31) ? len_keyword_ : 31;
    strncpy (buf, keyword, len);
    buf[len] = '\0';
    return verbose (buf);
}

int verbose (const char *keyword)
{
    unsigned int flags = Verbosity::get();

    if (!StringCmp (keyword, "LOGLIKE", case_ins))
    {
	return (flags & 1) ? 1 : 0;
    }
    if (!StringCmp (keyword, "RESULTS", case_ins))
    {
	return (flags & 2) ? 1 : 0;
    }
    if (!StringCmp (keyword, "ITERATE", case_ins))
    {
	return (flags & 4) ? 1 : 0;
    }
    if (!StringCmp (keyword, "SAMPLE", case_ins))
    {
	return (flags & 8) ? 1 : 0;
    }
    if (!StringCmp (keyword, "PHENOTYPES_LOAD", case_ins))
    {
	return (flags & 16) ? 1 : 0;
    }
    if (!StringCmp (keyword, "NUMERICAL", case_ins))
    {
	return (flags & 2048) ? 1 : 0;
    }
    if (!StringCmp (keyword, "ARTTOL", case_ins))
    {
	return (flags & 4096) ? 1 : 0;
    }
    if (!StringCmp (keyword, "SUBITER", case_ins))
    {
	return (flags & 0x40000) ? 1 : 0;
    }
    if (!Strcmp (keyword, "PREMAX"))
    {
	return (flags & 0x80000) ? 1 : 0;
    }
    if (!StringCmp (keyword, "LOGLIKE", case_ins))
    {
	return (flags & 1) ? 1 : 0;
    }
    if (!StringCmp (keyword, "LOGLIKE", case_ins))
    {
	return (flags & 1) ? 1 : 0;
    }
    if (!StringCmp (keyword, "MAXIMUM", case_ins))
    {
	return (flags == max_verbosity) ? 1 : 0;
    }
    return 0;
}




