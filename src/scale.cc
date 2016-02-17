/*
 * scale.cc implements the scale command using the Scale class
 * Written by Charles Peterson beginning in 2001
 * Copyright (c) 2001 Southwest Foundation for Biomedical Research
 */

#include <stdio.h>
#include "solar.h"

Scale *Scale::Scales[] = {0};
int Scale::Total_Count = 0;
bool Scale::_all = false;

static Scale *SNP = new Scale ("SNP",0.0);

extern "C" int ScaleCmd (ClientData clientData, Tcl_Interp *interp,
			 int argc, char *argv[])
{
// No argument, display all

    if (argc == 1)
    {
	Scale::Write_Commands (stdout);
	return TCL_OK;
    }

// Help

    if (argc == 2 && !Strcmp (argv[1], "help"))
    {
	return Solar_Eval (interp, "help scale");
    }
    
// Default all

    if (argc == 2 && !Strcmp (argv[1], "default_all"))
    {
	Scale::Reset();
	return TCL_OK;
    }

// Display one

    if (argc == 2)
    {
	Scale* s = Scale::Find (argv[1]);
	char buf[256];
	if (s == 0)
	{
	    sprintf (buf, "scale %s default", argv[1]);
	    RESULT_BUF (buf);
	    return TCL_ERROR;
	}
	sprintf (buf, "scale %s %-.10g", argv[1], s->adj());
	RESULT_BUF(buf);
	return TCL_OK;
    }

    if (argc == 3)
    {

// Delete scale

	if (!Strcmp ("default", argv[2]))
	{
	    Scale::Remove (argv[1]);
	    return TCL_OK;
	}
	if (!Strcmp ("default", argv[1]))
	{
	    Scale::Remove (argv[2]);
	    return TCL_OK;
	}


// Scale all variables to 0 (noscale all)
	if (!Strcmp ("all", argv[1]))
	{
	    if (0.0 != atof (argv[2]))
	    {
		RESULT_LIT ("Only 0 is allowed for scale all");
		return TCL_ERROR;
	    }
	    Scale::All();
	    return TCL_OK;
	}

// Add or modify scale

	double newadj;
	char buf[256];
	if (1 != sscanf (argv[2], "%lf %s", &newadj, buf))
	{
	    char outbuf[256];
	    sprintf (outbuf, "Scale adjustment %s is not numeric", buf);
	    RESULT_BUF (outbuf);
	    return TCL_ERROR;
	}
	Scale* s = Scale::Find (argv[1]);

// Modify scale

	if (s != 0)
	{
	    s->adj (newadj);
	}
	else
	{
	    Scale* s = new Scale (argv[1], newadj);
	    if (!(s->add()))
	    {
		RESULT_LIT ("Too many scales");
		return TCL_ERROR;
	    }
	}
	return TCL_OK;
    }
    RESULT_LIT ("Invalid scale command");
    return TCL_ERROR;
}

Scale* Scale::add()
{
    if (Total_Count < MAX_PARAMETERS)
    {
	Scales[Total_Count++] = this;
	return this;
    }
    return 0;
}

Scale* Scale::Find (const char *testname)
{
    int i;
    char testchars[5];
    strncpy (testchars, testname, 4);
    testchars[4] = '\0';
    if (_all) return Scales[0];

    for (i = 0; i < Total_Count; i++)
    {
	if (!Strcmp (testname, Scales[i]->_name))
	{
	    return Scales[i];
	}
    }
    return 0;
}

void Scale::Reset ()
{
    for (int i = 0; i < Total_Count; i++)
    {
	delete Scales[i];
    }
    Total_Count = 0;
    _all = false;
}

void Scale::All ()
{
    Scale::Reset();
    Scale* s = new Scale ("all", 0.0);
    s->add();
    _all = true;
}

void Scale::Remove (const char *testname)
{
    if (!Strcmp (testname, "all"))
    {
	Scale::Reset();
    }
    else
    {
	int i;
	for (i = 0; i < Total_Count; i++)
	{
	    if (!Strcmp (testname, Scales[i]->_name))
	    {
		delete Scales[i];
		--Total_Count;
		int j;
		for (j = i; j < Total_Count; j++)
		{
		    Scales[j] = Scales[j+1];
		}
		Scales[Total_Count] = 0;
	    }
	}
    }
}

void Scale::Write_Commands (FILE *fptr)
{
    int i;
    for (i = 0; i < Total_Count; i++)
    {
	fprintf (fptr, "scale %s %-.10g\n", Scales[i]->_name, Scales[i]->_adj);
    }
}

