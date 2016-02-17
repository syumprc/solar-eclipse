/*
 * mu.cc implements the Mu class and command
 * Written by Charles Peterson beginning on December 30, 1997
 * Copyright (c) 1997 Southwest Foundation for Biomedical Research
 */

#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "solar.h"
#include "expression.h"

extern int VirtualMultivar;


Expression *Mu::expression = 0;
double *Mu::vmeans = 0;
double *Mu::vmins = 0;
double *Mu::vmaxs = 0;
Context *Mu::context = 0;
char *Mu::before_string = 0;
char *Mu::after_string = 0;
char *Mu::combined_string = Strdup ("0");
bool Mu::_default_part = true;

// The following local statics are used during mu function evaluation
// (iff there is a non-defaulted mu)

static double *Var;
static double *Par;
static int *Male;
static int Trait;
static int Mean[MAX_TRAITS];


extern "C" int MuCmd (ClientData clientData, Tcl_Interp *interp,
		  int argc, char *argv[])
{
    if (argc == 2 && !StringCmp ("help", argv[1], case_ins))
    {
	return Solar_Eval (interp, "help mu");
    }

    if (argc == 1)
    {
	char *string = Mu::free_this_string();
	printf ("mu = %s\n", string);
	free (string);
	return TCL_OK;
    }

    if (argc == 2 && !StringCmp ("command", argv[1], case_ins))
    {
	char *substring = Mu::free_this_string();
	char* outstring = Strdup ("mu = ");
	StringAppend (&outstring, substring);
	RESULT_BUF (outstring);
	free (outstring);
	free (substring);
	return TCL_OK;
    }

    if (argc == 2 && !StringCmp ("reset", argv[1], case_ins))
    {
	Mu::reset ();
	return TCL_OK;
    }

    if (argc > 1 && argv[1][0] == '=')
    {
    // Set a new Mu equation
	char *eq_string = Tcl_Concat (argc-1, &argv[1]);
	char *test_string = Strdup (&eq_string[1]);
	char *expr_string = 0;

	squeeze (test_string);

	try  // Syntax exceptions might be thrown
	{
	    int mulen = strlen (test_string);
	    int li = mulen - 1;
	    bool expr_added = false;
	    if (strlen(test_string) > 3)
	    {
		if ((tolower(test_string[0])=='m') &&
		    (tolower(test_string[1])=='u') &&
		    test_string[2]=='+' || test_string[2]=='-')
		{
		    expr_string = &test_string[2];  // include leading op +/-
		    Mu::add_term_after (expr_string);
		    expr_added = true;
		}
		else if
		    ((tolower(test_string[li])=='u') &&
		     (tolower(test_string[li-1])=='m') &&
		     (test_string[li-2]=='+'))
		{
		    test_string[li-1] = '\0';  // include trailing +
		    expr_string = test_string;
		    Mu::add_term_before (expr_string);
		    expr_added = true;
		}
		else if (test_string[0]=='{' && strstr(test_string,"}"))
		{
		    expr_string = 1+strstr(test_string,"}");
		    Mu::reset();
		    Mu::add_term_after(expr_string);
		    expr_added = true;
		}
	    }
	    if (!expr_added)
	    {
		expr_string = test_string;
		Mu::new_expression (expr_string);
	    }
	}
	catch (Syntax_Error)
	{
	    RESULT_LIT ("Syntax Error");
	    return TCL_ERROR;
	}
	catch (Undefined_Function)
	{
	    RESULT_LIT ("Undefined Function");
	    return TCL_ERROR;
	}
	return TCL_OK;
    }

    RESULT_LIT ("Invalid mu command");
    return TCL_ERROR;
}

void Mu::add_term_before (const char *term) // trailing + required
{
    char* new_combined_string = Strdup (term);
    StringAppend (&new_combined_string, combined_string);
    Expression* new_expression = new Expression (new_combined_string);

// Success!

    if (expression) delete expression;
    if (combined_string) free (combined_string);
    expression = new_expression;
    combined_string = new_combined_string;

    char* new_before_string = Strdup (term);
    StringAppend (&new_before_string, before_string);
    if (before_string) free (before_string);
    before_string = new_before_string;
}

void Mu::add_term_after (const char *term)  // leading +/- required
{
    char* new_combined_string = Strdup (combined_string);
    StringAppend (&new_combined_string, term);
    Expression *new_expression = new Expression (new_combined_string);

// Success!

    if (expression) delete expression;
    if (combined_string) free (combined_string);
    expression = new_expression;
    combined_string = new_combined_string;

    StringAppend (&after_string, term);
}    
	
void Mu::new_expression (const char *estring)
{
    reset ();
    expression = new Expression (estring);

// Success !!!

    _default_part = false;
    combined_string = Strdup (estring);
}

// The MUCC function is an abbreviated mean function entry point
// used by discrete trait routine fun_mehd.f and also by mean_ below.

extern "C"
double mucc_ (int* male, double* var, double* par, int* trait)
{
    Var = var;
    Par = par;
    Male = male;
    if (VirtualMultivar)
    {
	Trait = ((int) Var[0]) - 1;
    }
    else
    {
	Trait = *trait - 1;
    }

    double mean = 0.;

// WARNING!  Covariate::eval must not be called if no default mu !!!
// Covariate binding is never done if no default mu !!!

    if (Mu::default_part_used()) 
    {
	mean = Covariate::eval (male, var, par, Trait);
    }

// Custom mu, in whole or part
    mean += Mu::eval ();
    return mean;
}

//
// The "MEAN" function used in FORTRAN is now defined here
//

extern "C" 
double mean_ (double *extra, double *par, double *var,
	      int *nextra, int *npar, int *ntrait, int *nvar,
	      int *ped, int *per, int *trait, int *male)
{
    return mucc_ (male, var, par, trait);
}

static double get_parameter_value (int pindex_0)
{
    return Par[pindex_0];
}

static double get_mean (int pindex_0)
{
    return Par[Mean[Trait]];
}

static double get_t1 (int pindex_0)
{
    return (double) (Trait==0);
}

static double get_t2 (int pindex_0)
{
    return (double) (Trait==1);
}

static double get_female (int pindex_0)
{
    return 1 - *Male;
}

static double get_male (int pindex_0)
{
    return *Male;
}
    
static double get_variable_value (int vindex)
{
    int pi = EqVar::get_packed_index (vindex);
    return Var[pi];
}

static double get_variable_mean (int vindex)
{
    int pi = EqVar::get_packed_index (vindex);
    return Mu::vmeans[pi];
}

static double get_variable_min (int vindex)
{
    int pi = EqVar::get_packed_index (vindex);
    return Mu::vmins[pi];
}

static double get_variable_max (int vindex)
{
    int pi = EqVar::get_packed_index (vindex);
    return Mu::vmaxs[pi];
}


int Mu::bind (Tcl_Interp *interp)
{
    int i;
    char buf[1024];

// If no user-defined expression, covariate.cc does everything

    if (!expression) return TCL_OK;

// First, build the Mu context
    if (context) delete context;
    context = new Context;

// Indexes to needed trait-dependent parameters
    int found_means = 0;
    for (i = 0; i < Trait::Number_Of(); i++)
    {
	if (-1 != (Mean[i] = Trait::Find_Parameter (i, "mean"))) found_means++;
    }

// The built-in stuff
    
    if (found_means == Trait::Number_Of())
    {
	context->add ("Mean", get_mean, 0);
    }
    context->add ("Male", get_male, 0);
    context->add ("Female", get_female, 0);
    context->add ("t1", get_t1, 0);
    context->add ("t2", get_t2, 0);

// The parameters

    Parameter *p;
    for (i = 0; p = Parameter::index (i); i++)
    {
	context->add (p->name(), get_parameter_value, i);
    }

// Try bind the Mu expression to the Mu context
//   Each undefined name should be a new variable needed


    bool ok = false;
    while (!ok)
    {
	try
	{
	    expression->bind (context);
	    ok = true;
	}
	catch (Undefined_Name& un)
	{
	    try
	    {
		if (!Strcmp (un.name, "Sex"))
		{
		    RESULT_LIT 
		("Use 'Male' or 'Female' instead of 'Sex' in Mu equation");
		    return TCL_ERROR;
		}
		if (strlen (un.name) > strlen ("x_"))
		{
		    int nindex = strlen ("x_");
		    char *tname = Strdup (un.name);
		    char tsave = tname[nindex];
		    tname[nindex] = '\0';
		    if (!StringCmp (tname, "x_", case_ins))
		    {
			tname[nindex] = tsave;
			char *vname = &tname[nindex];
			if (!Phenotypes::available (vname))
			    throw Undefined_Name (vname);
			int vindex = EqVar::add (vname);
			context->add (tname, get_variable_mean, vindex);
			free (tname);
			continue;
		    }
		    free (tname);
		}
		if (strlen (un.name) > strlen ("min_"))
		{
		    int nindex = strlen ("min_");
		    char *tname = Strdup (un.name);
		    char tsave = tname[nindex];
		    tname[nindex] = '\0';
		    if (!StringCmp (tname, "min_", case_ins))
		    {
			tname[nindex] = tsave;
			char *vname = &tname[nindex];
			if (!Phenotypes::available (vname))
			    throw Undefined_Name (vname);
			int vindex = EqVar::add (vname);
			context->add (tname, get_variable_min, vindex);
			free (tname);
			continue;
		    }
		    free (tname);
		}
		if (strlen (un.name) > strlen ("max_"))
		{
		    int nindex = strlen ("max_");
		    char *tname = Strdup (un.name);
		    char tsave = tname[nindex];
		    tname[nindex] = '\0';
		    if (!StringCmp (tname, "max_", case_ins))
		    {
			tname[nindex] = tsave;
			char *vname = &tname[nindex];
			if (!Phenotypes::available (vname))
			    throw Undefined_Name (vname);
			int vindex = EqVar::add (vname);
			context->add (tname, get_variable_max, vindex);
			free (tname);
			continue;
		    }
		    free (tname);
		}
		if (!Phenotypes::available (un.name))
		    throw Undefined_Name (un.name);
		int vindex = EqVar::add (un.name);
		context->add (un.name, get_variable_value, vindex);
	    }
	    catch (Undefined_Name& un)
	    {
		sprintf (buf, 
			 "Undefined variable `%s' in Mu equation", un.name);
		RESULT_BUF (buf);
		return TCL_ERROR;
	    }
	}
    }
    return TCL_OK;
}

void insert_string (char* string, const char* into)
{
    int move = strlen (into);
    int start = strlen (string);  // moves null also
    for (;start >= 0; start--)
    {
	string[start+move] = string[start];
    }
    for (;move > 0; move--)
    {
	string[move-1] = into[move-1];
    }
}


// Only used for default comment string
// Must be freed after use
// Line Wraps automatically
// uses comment character

char* Mu::comment_string ()
{
    if (0==Trait::Number_Of())
    {
	return Strdup ("# Mu requires trait to be specified");
    }
    char *mustring = Mu::free_this_string();

// Copy mu into buf large enough for line-wrapping

    int buflen = 3 * strlen (mustring);
    char* buf = (char*) Calloc (buflen, 1);
    strcpy (buf, "# mu = ");
    strcat (buf, mustring);
    free (mustring);

// While last line is greater than 76 characters, wrap line back to last +
// or simply force at 76

    int lastptr = 0;

// Must be sure lastptr does not go beyond string!!!

    while (76 < strlen (&buf[lastptr]))
    {
	bool wrapped = false;
	for (int i = 75; i >= 1; i--)
	{
	    if (buf[lastptr+i] == '+')
	    {
		insert_string (&buf[lastptr+i], "\n# ");
		lastptr = lastptr+i+3;
		wrapped = true;
		break;
	    }
	}
	if (!wrapped)
	{
	    insert_string (&buf[lastptr+75], "\n# ");
	    lastptr = lastptr+75+3;
	}
    }
    char* retbuf = Strdup (buf);
    free (buf);
    return retbuf;
}

char* Mu::free_this_string ()
{
    if (default_part_used())
    {
	char* ptr = Strdup ("");;
	if (before_string) StringAppend (&ptr, before_string);
	StringAppend (&ptr, "\\{");

// Single trait
	
	Covariate *c;
	if (1 ==Trait::Number_Of())
	{
	    StringAppend (&ptr, "Mean");
	    for (int i = 0; c = Covariate::index (i); i++)
	    {
		if (c->active())
		{
		    const char *must = c->mu_string (0);
		    if (must && strlen(must))
		    {
			StringAppend (&ptr, must);
		    }
		}
	    }
	} else {

// Multiple traits

	    for (int t = 0; t < Trait::Number_Of(); t++)
	    {
		char tnbuf[32];
		sprintf (tnbuf, "t%d*(<Mean(", t+1);
		StringAppend (&ptr, tnbuf);
		StringAppend (&ptr, Trait::Name(t));
		StringAppend (&ptr, ")>");
		for (int i = 0; c = Covariate::index (i); i++)
		{
		    if (c->active())
		    {
			const char *must = c->mu_string (t);
			if (must && strlen(must))
			{
			    StringAppend (&ptr, must);
			}
		    }
		}
		if (t+1 == Trait::Number_Of())
		{
		    StringAppend (&ptr, ")");
		} else {
		    StringAppend (&ptr, ") + ");
		}
	    }
	}
	StringAppend (&ptr, "\\}");
	if (after_string) StringAppend(&ptr, after_string);
	return ptr;
    }
    else
    {
	return Strdup (combined_string);
    }
}

void Mu::reset ()
{
    if (before_string) free (before_string), before_string = 0;
    if (after_string) free (after_string), after_string = 0;
    if (combined_string) free (combined_string), combined_string = 
			     Strdup ("0");
    if (expression) delete expression, expression = 0;
    _default_part = true;
}

    
void Mu::write_commands (FILE *file)
{
    if (default_part_used () && (!before_string && !after_string))
    {
	char* free_string = comment_string();
	fprintf (file, "%s\n", free_string);
	free (free_string);
    }
    else
    {
	char *mu_string = Mu::free_this_string ();
	fprintf (file, "mu = %s\n", mu_string);
	free (mu_string);
    }
}

