/*
 * maximize.cc implements the maximize command
 * Written by Charles Peterson beginning on October 29, 1997
 * Copyright (c) 1997 Southwest Foundation for Biomedical Research
 */

#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "solar.h"
// tcl.h from solar.h
#include "safelib.h"

int maximize (Tcl_Interp *interp);
double ccsearch (const char *outfilename, int get_who, bool sampledata);
int validate_solar (Tcl_Interp *interp);

static const int DEFAULT_COV_BND_RETRIES = 10;
static const int DEFAULT_COV_INCREASE_FACTOR = 5;
static const int STARTING = 4646;
static const int BOUNDARY_ERROR = 4747;
static int Who = 0;  // If 1, we only want to know who is included
static bool Sampledata = false;  // If true, we only want to write data
bool Maximize_InitPar = false;  // If true, we only want to initialize params
// At this time, who and sampledata are mutually exclusive
static bool Premax = false;    // Verbosity premax
static int Omega_Type = 0; // 0=unknown 1=sporadic, 2=polygenic

// Fortran subroutines to open and close FORTRAN units

extern "C" void openphen_ (char *phenfilename, int *status, int filenamelen);
extern "C" void closephen_ ();
extern "C" void openout_ (int *status);
extern "C" void closeout_ ();

// Save interp for when needed upstream in callbacks

Tcl_Interp* Save_Maximize_Interp = 0;

extern "C" int MaximizeCmd (ClientData clientData, Tcl_Interp *interp,
		  int argc, char *argv[])
{
    try {

    if (argc == 2 && !StringCmp (argv[1], "help", case_ins))
    {
	return Solar_Eval (interp, "help maximize");
    }
    Premax  = verbose ("PREMAX");
    if (Premax)
    {
	fprintf (stderr, "    **  Entering cmaximize\n");
    }
    if (argc == 1 || argc == 2 || argc == 3)
    {
	int status = STARTING;

	int total_retries = DEFAULT_COV_BND_RETRIES;
	if (!Solar_Eval (interp, "covboundretries"))
	{
	    total_retries = atoi (Tcl_GetStringResult (interp));
	}
	int retries_left = total_retries;

	double increase_factor = DEFAULT_COV_INCREASE_FACTOR;
	if (!Solar_Eval (interp, "covboundincr"))
	{
	    increase_factor = atof (Tcl_GetStringResult (interp));
	}

	while (status == STARTING || 
	       (status == BOUNDARY_ERROR && retries_left-- > 0))
	{
	    if (status == BOUNDARY_ERROR)
	    {
		char pnames[1024];
		Covariate::boundary_check (increase_factor, pnames, 1024);
		if (verbose ("ITERATE"))
		{
		    fprintf (stderr, 
		 "\n    *** Retry with new boundaries for: %s\n\n", 
			     pnames);
		}
	    }

	    const char *outfilename = "solar.out";
	    Who = 0;
	    Sampledata = false;
	    Maximize_InitPar = false;
	    if (argc==2)
	    {
		if (!Strcmp(argv[1],"-who"))
		{
		    Who = 1;
		}
		else if (!Strcmp(argv[1],"-runwho"))
		{
		    Who = 2;
		}
		else if (!Strcmp(argv[1],"-sampledata"))
		{
		    Sampledata = true;
		}
		else
		{
		    outfilename = argv[1];
		}
	    }
	    else if (argc==3)
	    {
		if (!Strcmp (argv[1],"-who"))
		{
		    Who = 1;
		    outfilename = argv[2];
		}
		else if (!Strcmp(argv[1],"-runwho"))
		{
		    Who = 2;
		    outfilename = argv[2];
		}
		else if (!Strcmp (argv[2],"-who"))
		{
		    Who = 1;
		    outfilename = argv[1];
		}
		else if (!Strcmp(argv[2],"-runwho"))
		{
		    Who = 2;
		    outfilename = argv[1];
		}
		else if (!Strcmp (argv[1],"-sampledata"))
		{
		    Sampledata = true;
		    outfilename = argv[2];
		}
		else if (!Strcmp (argv[2],"-sampledata"))
		{
		    Sampledata = true;
		    outfilename = argv[1];
		}
		else if (!Strcmp (argv[1],"-initpar"))
		{
		    Maximize_InitPar = true;
		    outfilename = argv[2];
		}
		else if (!Strcmp (argv[2],"-initpar"))
		{
		    Maximize_InitPar = true;
		    outfilename = argv[1];
		}
		else
		{
		    RESULT_LIT ("Invalid arguments to maximize");
		    return TCL_ERROR;
		}
	    }
	    status = Loglike::maximize (interp, outfilename);
	}
	if (status == BOUNDARY_ERROR)
	{
	    char buf[256];
	    sprintf (buf, "Covariate boundary errors after %d retries",
		     total_retries);
	    RESULT_BUF (buf);
	    return TCL_ERROR;
	}
	return status;
    }
    RESULT_LIT ("Invalid maximize command.  Use help maximize.\n");
    return TCL_ERROR;

    }
    catch (Safe_Error_Return& ser)
    {
	fprintf (stderr, "%s\n", ser.message());
	exit (EXIT_FAILURE);
    }
}


int Loglike::maximize (Tcl_Interp *interp, const char *outfilename)
{
    Save_Maximize_Interp = interp;

    Loglike::_status = failed;
    EVD::Reset ();

    if (Premax)	fprintf (stderr, "    **  Validating key\n");
    if (!validate_solar (interp)) return TCL_ERROR; // Validate registration

// Free Equation Variable storage (avoid memory leak)

    EqVar::bind (interp);

// check Trait

    Zscore::Reset();
    if (Premax) fprintf (stderr, "    **  Binding trait\n");
    if (Trait::Bind (interp)) return TCL_ERROR;

// Check Phenotypes

    if (Premax) fprintf (stderr, "    **  Binding phenotypes\n");
    if (Phenotypes::bind (interp)) return TCL_ERROR;

// Check FisherPedigree

    if (Premax) fprintf (stderr, 
			 "    **  Joining pedigree and phenotypes data\n");
    if (FisherPedigree::bind (interp)) return TCL_ERROR;


    if (Omega::check_if_defined (interp)) return TCL_ERROR;

    char* omega = strdup (Omega::expression()->string());
    int ptr=0;
    while (omega[ptr] != '\0')
    {
	omega[ptr] = tolower (omega[ptr]);
	ptr++;
    }
    int found_delta7 = !(!strstr (omega, "delta7"));

//  fprintf (stderr, "Omega is %s\n", omega);
//  fprintf (stderr, "Bool: %d\n", found_delta7);

// Check omega for sporadic/polygenic model type

	Omega_Type = 0;
	squeeze (omega);
	if (!Strcmp (omega,"pvar*(phi2*h2r+i*e2)"))
	{
	    if (Constraint::Check_by_name ("h2r"))
	    {
		Omega_Type = 1;
	    }
	    else
	    {
		Omega_Type = 2;
	    }
	}
	free (omega);

// If using discrete/mixed code, be sure phi2 matrix is loaded
//   and delta7 is loaded only if needed.
// Also constrain SD to 1.0...if actually discrete

    if (Strcmp("Default",Option::get_string("ModelType")) ||
	(Trait::Maximization_Type() == discrete && 
	 Option::get_int ("EnableDiscrete")))
    {
//	if (Trait::Number_Of() > 1) {
//	    RESULT_LIT (\
//	      "Mixed Discrete and Quantitative Bivariate not yet supported");
//	    return TCL_ERROR;
//	}
	if (Premax) fprintf (stderr, "    **  Binding mixed/discrete code stuff\n");
	const char *message = 0;
	if (!Matrix::find ("phi2"))
	{
	    if (found_delta7)
	    {
		message = Matrix::setup (0,"phi2","phi2","delta7");
	    }
	    else
	    {
		message = Matrix::setup (0,"phi2","phi2");
	    }
	}
	else  // phi2 already loaded...but load delta7 if needed
	{
	    if (found_delta7)
	    {
		if (!Matrix::find ("delta7"))
		{
		    message = Matrix::setup (0, "phi2","phi2","delta7");
		}
	    }

// though well intentioned, unloading matrix causes problem if memory short
// already.  This may be fixed in next update
#if 0
	    else // unload delta7 if not needed
	    {
		if (Matrix::find ("delta7"))
		{
		    message = Matrix::setup (0, "phi2","phi2");
		}
	    }
#endif

	}
	if (message)
	{
	    fprintf (stderr, "%s\n", message);
	    RESULT_LIT (
		"Error loading phi2 matrix required for discrete traits");
	    return TCL_ERROR;
	}
	for (int itrait = 0; itrait < Trait::Number_Of(); itrait++)
	{
	    if (discrete == Trait::Type(itrait))
	    {
// If polyclasses, cycle through them, otherwise null prefix
		char polysuffix[32];
		const char* polyclasses = Option::get_string ("polyclasses");
		while (1)
		{
		    polysuffix[0] = 0;
		    if (strlen (polyclasses))
		    {
			int cindex = 0;
			polysuffix[0] = 0;
			strcpy (polysuffix,"_c");
			while (polyclasses[cindex] != ',' && 
			       polyclasses[cindex] != 0)
			{
			    polysuffix[cindex+2] = polyclasses[cindex];
			    cindex++;
			}
			polysuffix[cindex+2] = 0;
			if (polyclasses[cindex] == 0)
			{
			    polyclasses = polyclasses + cindex;
			} 
			else
			{
			    polyclasses = polyclasses + cindex + 1;
			}
		    }
		    char command[1024];
		    if (Trait::Number_Of() > 1)
		    {
			sprintf (command, "constraint <sd%s(%s)> = 1", 
				 polysuffix, Trait::Name(itrait));
		    }
		    else 
		    {
			sprintf (command, "constraint sd%s = 1",
				 polysuffix);
		    }
		    if (Solar_Eval (interp, command) == TCL_ERROR)
		    {
			RESULT_LIT 
           ("Error setting constraint sd = 1 as required for discrete traits");
			return TCL_ERROR;
		    }
		    char mean_basename[1024];
		    sprintf (mean_basename, "mean%s", polysuffix);
		    int pindex = Trait::Find_Parameter (itrait, mean_basename);
		    if (pindex > -1)
		    {
			Parameter *p = Parameter::index(pindex);
			if (p)
			{
			    if (p->lower == 0.0 && p->upper == 0.0)
			    {
				p->lower = -8.0;
				p->upper = 8.0;
			    }
			}
		    }
		    if (!strlen(polyclasses)) break;
		}
	    }
	}
    }

// Continue with remaining tests

    int retval;
    retval = Loglike::max_over_this_phen (interp, outfilename);
    return retval;
}

// Middle-inner part, assuming Phenotypes and FisherPedigree are OK

int Loglike::max_over_this_phen (Tcl_Interp *interp, const char *outfilename)
{
// ensure matrices are happy (ibdid present)

    if (Premax) fprintf (stderr, "    **  Binding matrices\n");
    if (Matrix::bind (interp)) return TCL_ERROR;

// bind all variables in Omega equation

    if (Premax) fprintf (stderr, "    **  Binding omega\n");
    if (Omega::bind (interp)) return TCL_ERROR;

// bind all variables in Mu equation

    if (Premax) fprintf (stderr, "    **  Binding mu\n");
    if (Mu::bind (interp)) return TCL_ERROR;

// bind Covariates

    if (Premax) fprintf (stderr, "    **  Binding covariates\n");
    if (Covariate::bind (interp)) return TCL_ERROR;

// Save "start" as "last" values

    if (Premax) fprintf (stderr, "    **  Binding parameters\n");
    if (Parameter::bind (interp)) return TCL_ERROR;

// bind Constraints

    if (Premax) fprintf (stderr, "    **  Binding constraints\n");
    if (Constraint::bind (interp)) return TCL_ERROR;

// Change to output directory

    char dirname[1024];
    char filename[1024];
    strcpy (dirname, "");

    if (Premax) fprintf (stderr, "    **  Moving to output directory\n");
    char *fullfilename = append_extension (outfilename, ".out");

    char buf[1024];
    sprintf (buf, "file tail %s", fullfilename);
    if (TCL_OK != Solar_Eval (interp, buf))
    {
	free (fullfilename);
	return TCL_ERROR;
    }
    strcpy (filename, Tcl_GetStringResult (interp));
    Tcl_ResetResult (interp);

    sprintf (buf, "file dirname %s", fullfilename);
    if (TCL_OK != Solar_Eval (interp, buf))
    {
	free (fullfilename);
	return TCL_ERROR;
    }
    free (fullfilename); // done with this now

    strcpy (dirname, Tcl_GetStringResult (interp));
    Tcl_ResetResult (interp);
    if (!strcmp (dirname, "."))
    {
	strcpy (dirname, "");
    }
    else
    {
	char olddir[1024];
	if (!getcwd (olddir, 1024))
	{
	    RESULT_LIT ("Unable to save current directory");
	    return TCL_ERROR;
	}
	if (chdir (dirname))
	{
	    sprintf (buf, "Unable to access directory %s", dirname);
	    RESULT_BUF (buf);
	    return TCL_ERROR;
	}
	strcpy (dirname, olddir);
    }
    int retval = Loglike::max_inside_outdir (interp, filename);
    if (strlen (dirname))
    {
	if (chdir (dirname))
	{
	    fprintf (stderr, "Unable to restore working directory");
	    retval = TCL_ERROR;
	}
    }
    return retval;
}


// inner-inner part assuming we are in output directory

int Loglike::max_inside_outdir (Tcl_Interp *interp, const char *outfilename)
{
// Save current model as "last.mod"
    if (Premax) fprintf (stderr, "    **  Saving model\n");
    FILE *lastmod = fopen ("last.mod", "w");
    if (lastmod)
    {
	Model::write (lastmod);
	Fclose (lastmod);
    }
    else
    {
	RESULT_LIT ("Error opening last.mod in output directory");
	return TCL_ERROR;
    }
	
// Open FORTRAN units

    if (Premax) fprintf (stderr, "    **  Opening FORTRAN units\n");
    int STATUS;
    openout_ (&STATUS);
    if (STATUS)
    {
	
        if (STATUS == 1 || STATUS == 2)
        {
	    RESULT_LIT ("Error opening FORTRAN scratch file");
        }
        else if (STATUS == 3)
        {
	    RESULT_LIT ("Error opening output file");
        }
        else // STATUS == 4
        {
	    RESULT_LIT ("Error opening simple output file");
        }
	return TCL_ERROR;
    }

    int retval = TCL_ERROR;
    try
    {
	if (Premax) fprintf (stderr, "    **  Entering search wrapper\n");
	double loglike = ccsearch (outfilename,Who,Sampledata);
	char pnames[1024];
	if (Covariate::boundary_check(0,pnames,1024))
	{
	    retval = BOUNDARY_ERROR;
	    loglike = 0;
	}
	if (loglike != 0.0 && !isnan (loglike))
	{
	    Loglike::loglike (loglike);
	    Loglike::_status = valid;
	    retval = TCL_OK;
	}
	else if (Who==1 || Sampledata || Maximize_InitPar)
	{
	    retval = TCL_OK;
	}
	else
	{
	    RESULT_LIT ("Convergence failure");  // Convergence error
	}
    }
    catch (Safe_Error_Return& ser)  //errors requiring user intervention
    {
	char Error_Message[1024];
	strcpy (Error_Message, "\n");
	strncpy (&Error_Message[1], ser.message(), 1023);
	RESULT_BUF (Error_Message);
    }
    closeout_ ();
    char combuf[1024];
    FILE *hgm = 0;
    if (FisherPedigree::HouseholdGroupMessage[0] != '\0' && 
	(hgm = fopen ("houseinfo.tmp", "w")))
    {
	cfprintf (hgm, "\n\n%s", FisherPedigree::HouseholdGroupMessage);
	fclose (hgm);
	sprintf (combuf, "cat solar.smp houseinfo.tmp solar.hst >\"%s\"", 
		 outfilename);
	system (combuf);
	system ("rm -f solar.smp houseinfo.tmp solar.hst");
    }
    else
    {
	sprintf (combuf, "cat solar.smp solar.hst >\"%s\"", outfilename);
	system (combuf);
	system ("rm -f solar.smp solar.hst");
    }

#ifdef MATRIX_DEBUG
    extern FILE *MDfile;
    if (MDfile) fclose (MDfile);
#endif

    return retval;
}

extern "C" void converr_ (double *f, double* par)
{
    bool looks_ok = false;
    if (*f != 0 && *f > -1e36 && *f < 1e36)
    {
	looks_ok = true;
	for (int i = 0; i < Parameter::count(); i++)
	{
	    Parameter *p = Parameter::index(i);
	    double result = par[i];
	    if (fabs(result) < -1e36 || fabs(result) > 1e36)
	    {
		looks_ok = false;
		result = 0.0;
	    }
	    else if (fabs(result) < MIN_VISIBLE)
	    {
		result = 0.0;
	    }
	    Parameter::set (i, result, p->lower, p->upper, p->se, p->score);
	}
    }
    
//    if (verbose ("ITERATE"))
//    {
//	fprintf (stderr, "\n    **  Unable to solve quadratic\n");
//    }
    if (looks_ok)
    {
	throw Safe_Error_Return ("Convergence failure (Restartable) (Unable to Solve Quadratic)");
    } else {
	throw Safe_Error_Return ("Convergence failure (Unsolvable Quadratic)");
    }
	
}

// FORTRAN interface to Omega_Type

extern "C" int omegatyp_ ()
{
    return Omega_Type;
}

Tcl_Interp* get_maximize_interp ()
{
    return Save_Maximize_Interp;
}
