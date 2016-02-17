/*
 * option.cc implements the option command
 * Written by Charles Peterson beginning on October 27, 1997
 * Copyright (c) 1997 Southwest Foundation for Biomedical Research
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "solar.h"
// tcl.h from solar.h
#include "safelib.h"

Option *Option::Options[] = {0};
int Option::count = 0;

extern "C" int OptionCmd (ClientData clientData, Tcl_Interp *interp,
	       int argc, char *argv[])
{
    if (argc == 2 && !StringCmp ("help", argv[1], case_ins))
    {
	return Solar_Eval (interp, "help option");
    }

    if (argc == 1)
    {
	Option::show_all ();
	return TCL_OK;
    }

    if (argc == 2)
    {
	if (!Strcmp (argv[1], "return"))
	{
	    return Option::return_commands( interp );
	}
	return Option::query (argv[1], interp);
    }

    if ((argc == 4 && !Strcmp(argv[2],"+")) || (argc==3 && '+'== argv[2][0]))
    {
	char* addin;
	if (argc==4) {
	    addin = argv[3];
	} else {
	    addin = &argv[2][1];
	    if (!strlen(addin))
	    {
		RESULT_LIT ("Invalid option syntax: + not followed by token");
		return TCL_ERROR;
	    }
	}
		
	return Option::set (argv[1], addin, true, interp);
    }
    if (argc == 3)
    {
	return Option::set (argv[1], argv[2], false, interp);
    }

    RESULT_LIT ("Invalid option command");
    return TCL_ERROR;
}


int Option::get_int (const char *name)
{
    for (int i = 0; i < Option::count; i++)
    {
	if (!StringCmp (name, Option::index(i)->name, case_ins))
	{
	    char dummy[256];
	    int integer;
	    if (1 == sscanf (Option::index(i)->value, "%d %s", &integer,
			     &dummy))
	    {
		return integer;
	    }
	    char message[256];
	    sprintf (message, "Option %s must be integer", name);
	    throw Safe_Error_Return (message);
	}
    }
    char message[256];
    sprintf (message, "Missing option %s", name);
    throw Safe_Error_Return (message);
}

// Fortran interface to options, returns double

extern "C" void ioption_ (char *fstring, int *value, int flen)
{
    char cstring[MAX_OPTION_NAME_LEN+1];
    int uselen = (flen < MAX_OPTION_NAME_LEN) ? flen : MAX_OPTION_NAME_LEN;
    strncpy (cstring, fstring, uselen);
    cstring[uselen] = '\0';
// Apply lnblank
    for (int i = uselen-1; i >= 0; i++)
    {
	if (!isspace (cstring[i])) break;
	cstring[i] = '\0';
    }
    *value = Option::get_int (cstring);
}

extern "C" void doption_ (char *fstring, double *value, int flen)
{
    char cstring[MAX_OPTION_NAME_LEN+1];
    int uselen = (flen < MAX_OPTION_NAME_LEN) ? flen : MAX_OPTION_NAME_LEN;
    strncpy (cstring, fstring, uselen);
    cstring[uselen] = '\0';
// Apply lnblank
    for (int i = uselen-1; i >= 0; i++)
    {
	if (!isspace (cstring[i])) break;
	cstring[i] = '\0';
    }
    *value = Option::get_double (cstring);
}

extern "C" void soption_ (char* fstring, char* value, int flen, int vlen)
{
    char cstring[MAX_OPTION_NAME_LEN+1];
    int uselen = (flen < MAX_OPTION_NAME_LEN) ? flen : MAX_OPTION_NAME_LEN;
    strncpy (cstring, fstring, uselen);
    cstring[uselen] = '\0';
// Apply lnblank
    for (int i = uselen-1; i >= 0; i++)
    {
	if (!isspace (cstring[i])) break;
	cstring[i] = '\0';
    }
    const char* vstring = Option::get_string (cstring);
    if (!vstring) vstring = "";
    strncpy (value, vstring, vlen);
    fortspad (value, vlen);
    value[0] = toupper (value[0]);
    for (int i = 1; i < vlen; i++)
    {
	value[i] = tolower (value[i]);
    }
}


const char* Option::get_string (const char* name)
{
    for (int i = 0; i < Option::count; i++)
    {
	if (!Strcmp (name, Option::index(i)->name))
	{
	    return Option::index(i)->value;
	}
    }
    char message[256];
    sprintf (message, "Missing option %s", name);
    throw Safe_Error_Return (message);
}    


double Option::get_double (const char *name)
{
    for (int i = 0; i < Option::count; i++)
    {
	if (!StringCmp (name, Option::index(i)->name, case_ins))
	{
	    char dummy[256];
	    double d;
	    if (1 == sscanf (Option::index(i)->value, "%lg %s", &d,
			     &dummy))
	    {
		return d;
	    }
	    char message[256];
	    sprintf (message, "Option %s must be double", name);
	    throw Safe_Error_Return (message);
	}
    }
    char message[256];
    sprintf (message, "Missing option %s", name);
    throw Safe_Error_Return (message);
}

int Option::query (const char *name, Tcl_Interp *interp)
{
    for (int i = 0; i < Option::count; i++)
    {
	if (!StringCmp (name, Option::index(i)->name, case_ins))
	{
	    Option* op = Option::index(i);
	    char *outstring = Strdup ("");

	    string__append (&outstring, op->value);
	    if (op->list_type)
	    {
		for (int j = 0; j < op->list_len; j++)
		{
		    string__append (&outstring, " \n");
		    string__append (&outstring, op->list[j]);
		}
	    }
	    RESULT_BUF (outstring);
	    free (outstring);
	    return TCL_OK;
	}
    }
    RESULT_LIT ("No such option");
    return TCL_ERROR;
}

int Option::set (const char *name, const char *value, bool append, 
		 Tcl_Interp *interp)
{
    for (int i = 0; i < Option::count; i++)
    {
	if (!StringCmp (name, Option::index(i)->name, case_ins))
	{
	    Option *opt = Option::index(i);
	    if (append)
	    {
		if (!opt->list_type)
		{
		    RESULT_LIT ("this option does not allow +");
		    return TCL_ERROR;
		}
		if (!Strcmp (value,opt->initial_value))
		{
		    RESULT_LIT ("+ cannot be used with null option value");
		    return TCL_ERROR;
		}
		opt->append (value);
	    }
	    else
	    {
		if (opt->list_type && opt->list_len > 0)
		{
		    opt->free_list();
		}
		if (opt->value != opt->initial_value) 
		    free ((char*) opt->value);
		opt->value = Strdup (value);
	    }
	    return TCL_OK;
	}
    }
    RESULT_LIT ("No such option");
    return TCL_ERROR;
}

void Option::free_list ()
{
    for (int j = 0; j < list_len; j++)
    {
	free (list[j]);
    }
    free (list);
    list = 0;
    list_len = 0;
}


void Option::show_all ()
{
    for (int i = 0; i < Option::count; i++)
    {
	Option *opt = Option::index (i);
	if (opt->visible)
	{
	    printf ("option %s %s\n", opt->name, opt->value);
	    if (opt->list_type)
	    {
		for (int j = 0; j < opt->list_len; j++)
		{
		    printf ("option %s + %s\n", opt->name, opt->list[j]);
		}
	    }
	}
    }
}

Option::Option (const char *name_, const char *value_)
{
    name = Strdup (name_);
    value = Strdup (value_);
    initial_value = value;
    visible = 1;
    writeable = 1;
    list_type = false;
    session = false;
    list = 0;
    list_len = 0;
}

Option::~Option ()
{
    if (Strcmp (value, initial_value))
    {
	free (initial_value);
    }
    free (name);
    free (value);
    if (list_len) free_list();
}

void Option::add (const char *name, const char *value)
{
    Options[count++] = new Option (name, value);
}

// add invisible options
void Option::addi (const char *name, const char *value)
{
    Options[count++] = new Option (name, value);
    Options[count-1]->visible = 0;
}

// add options that don't get written to models
void Option::addnw (const char *name, const char *value)
{
    Options[count++] = new Option (name, value);
    Options[count-1]->writeable = 0;
}

// add session options that stay active for an entire session until canceled
void Option::addses (const char *name, const char *value)
{
    Options[count++] = new Option (name, value);
    Options[count-1]->writeable = 0;
    Options[count-1]->session = true;
}

// add list-type options
void Option::addl (const char *name, const char *value)
{
    Options[count++] = new Option (name, value);
    Options[count-1]->list_type = true;
    Options[count-1]->list = 0;
    Options[count-1]->list_len = 0;
}


// The way lists and append works requires a little description:
//
// Synopsis: the list is an extension to, not a replacement of, the scalar
// base value.  (Hmmn, this is beginning to sound like lisp car and cdr.)
//
// The scalar base element, often zero, is at the base of the list
//   The base element has "value" and "initial value".
// If set without +, the base element only is changed, and a pre-existing
//   list is discarded.
// If set with +, there are four possibilities...
//    If the there is no list, and the scalar base value is the "initial value"
//      the scalar base value is changed, and the list remains empty.
//    If there is no list, but the scalar base value is NOT the initial value,
//      the added element goes at the beginning of the extension list.
//    If the there is a list, then the new values goes at the end of it.
//    If an attempt is made to + the initial value, that is an error; the
//      initial value is assumed to be some kind of null/off/default for which
//      + is semantically wrong.  That is true even though we nicely allow
//      + when the current value is initial value, in which case actually
//      a substitution is being performed.
//
//    From the starting state, option can specify + or not.  (Not would be
//       more logical, but + may be more convenient in a loop.)
//
//    Now when reading the list, we can always start with the scalar value.
//      the complexity having been moved into the add operation.
//    


void Option::append (const char* newvalue)
{
    if (list_len == 0 && !Strcmp (value, initial_value))
    {
	value = Strdup (newvalue);
    }
    else if (list_len == 0)
    {
	list = (char**) Calloc (1, sizeof(char*));
	list[0] = Strdup (newvalue);
	list_len = 1;
    }
    else
    {
	list = (char**) realloc (list, ++list_len*sizeof(char**));
	list[list_len-1] = Strdup (newvalue);
    }
}


void Option::setup ()
{
    add ("ModelType","Default");
    add ("CMDiagonal", "0");
    add ("StandErr", "1");
    addi ("RobustSE", "0");
    add ("StandLogLike", "0");
    add ("AutoCovarBound", "1.25");
    add ("Grid", "0");
    add ("GridPoints", "1");
    add ("MaxIter", "1000");
    add ("Outlier", "0");
    add ("CutPeople", "0.05");
    add ("CutPed", "1.0");
    add ("TDist", "0");
    add ("Conv", "1.0E-6");
    add ("Conv(Discrete)", "1.0E-4");
    add ("NConv", "4");
    add ("Tol", "1.0E-8");
    add ("MaxStep", "5");
    add ("BCliff", "0.1");
    add ("MaxCliffs","15");
    add ("ScoreOnlyIndex", "-1");
    add ("MergeHousePeds", "1");
    add ("MergeAllPeds", "0");
    add ("RobustEst", "0");
    add ("Tune", "3");
     addl ("PedSelect","0");
    add ("HouseGroupShow","0");
    add ("EnableDiscrete", "1");
    add ("DiscreteOrder", "1");
    add ("DiscreteMethod","1");
    add ("UnbalancedTraits","1");
    add ("CorrectDeltas","0");
    add ("EnforceBounds","1");
    add ("BounDiff","0");
    add ("EnforceConstraints","0");
    add ("AbsVarianceParms","1");
    add ("SingularTrait","0");
    add ("zscore", "0");
    add ("zmean1", "0.0");
    add ("zmean2", "0.0");
    add ("zsd1", "0.0");
    add ("zsd2", "0.0");
    add ("ParameterFormat","16");
    add ("MatrixNumberFormat","15");
    add ("CheckPhenIDs","1");
    add ("RicVolOffset","1");
    add ("PolyClasses","");
    addnw ("SampleSameTrustMe", "0");
    addnw ("EvdPhase", "0");
    addnw ("DontAllowSampleChange","0");
    addnw ("EvdOptimizations", "0");
    add   ("Eigenvectors", "0");
    add   ("EVDmat", "0");
    add   ("EVDcovs", "0");
    add   ("FPHIMethod","1");
    add   ("ResetRandom", "0");
    addses ("CsvBufSize","0");
    addses ("ExpNotation","0");
    addi ("PedLike","0");
    addi ("LenInteger", "64000");
    addi ("LenReal", "64000");
    addi ("LenExtra", "2");
    add  ("ShuffleReseeding", "1");
}


void Option::write_commands (FILE *file)
{
    Option *opt;
    for (int i = 0; opt = Option::index (i); i++)
    {
	if (!strcmp (opt->value, opt->initial_value)) continue;
	if (opt->writeable)
	{
	    fprintf (file, "option %s %s\n", opt->name, opt->value);
	    if (opt->list_type)
	    {
		for (int j = 0; j < opt->list_len; j++)
		{
		    fprintf (file, "option %s + %s\n", opt->name,
			     opt->list[j]);
		}
	    }
	}
    }
}

int Option::return_commands (Tcl_Interp* interp)
{
    Option *opt;
    std::string rstring;
    bool exist = false;
    bool first = true;
    for (int i = 0; opt = Option::index (i); i++)
    {
	if (!strcmp (opt->value, opt->initial_value)) continue;

	exist = true;
	if (!first) rstring += "\n";
	rstring += "option ";
	rstring += opt->name;
	rstring += " ";
	rstring += opt->value;
	first = false;

	if (opt->list_type)
	{
	    for (int j = 0; j < opt->list_len; j++)
	    {
		rstring += "\noption ";
		rstring += opt->name;
		rstring += " + ";
		rstring += opt->value;
		if (Option::index(i+1)) rstring += "\n";
	    }
	}
    }
    if (exist)
    {
	RESULT_BUF (rstring.c_str());
    }
    else
    {
	RESULT_LIT ("");
    }
    return TCL_OK;
}


void Option::reset()
{
    Option *opt;
    for (int i=0; opt=Option::index(i); i++)
    {
	if (opt->session) continue;

	opt->value = opt->initial_value;
	if (opt->list_type)
	{
	    opt->free_list();
	}
    }
}


extern "C" void pedselect_ (int* iped, int* iselect)
{
    *iselect = 1;
    int start = Option::get_int ("PedSelect");
    if (start==0)
    {
	return;
    }

    if (Option::check_list ("PedSelect", *iped))
    {
	return;
    }
    *iselect = 0;
    return;
}

bool Option::check_list (const char* name, int ivalue)
{
    int integer;
    char dummy[1024];

    for (int i = 0; i < Option::count; i++)
    {
	if (!Strcmp (name, Option::index(i)->name))
	{
	    Option* opt = Option::index(i);
	    if (1 == sscanf (opt->value, "%d %s", &integer, dummy))
	    {
		if (integer == ivalue)
		{
		    return true;
		}
	    }
	    for (int j = 0; j < opt->list_len; j++)
	    {
		if (1 == sscanf (opt->list[j], "%d %s", &integer, dummy))
		{
		    if (integer == ivalue)
		    {
			return true;
		    }
		}
	    }
	    return false;
	}
    }
    return false;
}
