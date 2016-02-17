/*
 * covariate.cc implements the Covariate class and covariate command
 * Written by Charles Peterson beginning in 1997
 * Copyright (c) 1997 Southwest Foundation for Biomedical Research
 */

#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "safelib.h"
// tcl.h from solar.h
#include "solar.h"

/*
 * For better or worse, the current implementation uses pointers (i.e.
 * hard links) to parameters.  This is enforced by not allowing user
 * deletion of parameters, and requires certain ordering of reset
 * operations (covariates deleted before parameters).  Soft-linking by
 * name at "bind" time might be better.
 */

/*
 * "Generic" or "unqualified" covariates apply to all trait(s).
 * Covariates qualified by a trait name (suffixed in parentheses)
 * apply while the trait is active, otherwise become inactive.
 * Inactive covariates still force a variable to be included in the
 * analysis (just as with suspended covariates).  Because of this,
 * you can also have covariates simply for the purpose of forcing
 * variable inclusion.  To make that intuitive, it is possible to
 * to have a covariate qualified by the "null" trait.
 *
 * covar sex        sex covariate applied to any/all traits ("generic")
 * covar age(q4)    age covariate only applied to q4 ("specific")
 * covar weight()   weight variable required for inclusion in analysis
 *
 * When trait(s) change, the "new_traits" method be called must clear out
 * beta parameters which are no longer applicable.
 */

/*
 * There are two orthogonal dimensions of effectiveness:
 *   non-suspended vs. suspended
 *   applicable vs. not-applicable
 * 
 * Active covariates are those which are applicable (to a current trait)
 * and not suspended.
 *
 * Covariate suspension is a temporary way of removing covariates from a
 * model for hypothesis testing.  The beta parameter, and its estimates
 * and boundaries remain, but they are ignored, effectively excised from
 * the mu equation.  The underlying variable (phenotype) continues to be
 * required in the sample.
 *
 * Applicability is related to the genericity or specificity of a covariate.
 * Generic covariates are applicable to all traits.  Specific covariates
 * are specific to particular traits, or to the null trait.  If a covariate
 * is specific to a trait not in the current model, it is inapplicable (and
 * also not effective).  Covariates specific to the null trait are never
 * applicable (or effective) but still require their variable.
 */


/*
 * ABOUT COVARIATE SUSPENSION...an unrequested feature which causes some
 *   misunderstandings in its current form and might be removed or changed
 *   in future...I've thought about retaining feature internally but
 *   making it appear as a virtual constraint.
 *
 * For hypothesis testing, the user should set up all covariates and
 *   then suspend them one by one.
 *
 * Historically, there was no such thing as "suspension" and the user
 *   would constrain them to zero instead.  During such constraining,
 *   the previous optimized value would be lost.  During "suspension"
 *   we save the current value, which is re-loaded during "restore."
 *
 * Suspension is functionally equivalent to zero constraining, but more
 *   efficient (mainly because the term is removed from MU evaluation
 *   evaluation, which is one of the biggest users of processing time).
 *   As with constraining, the requirement that the variables
 *   be present in the user sample is maintained.
 */

Context* Covariate::_Context = 0;

// Save trait mean indexes at bind time

static int Mean[MAX_TRAITS];
int Covariate::last_number_of_traits;

CovariateTerm::CovariateTerm (const char *vname, int exp, 
			      const char* qualifier)
{
    name = Strdup(vname); 
    if (qualifier)
    {
	trait = Strdup(qualifier);
    }
    else
    {
	trait = 0;
    }
    exponent=exp; 
    sex=(Strcmp(name,"sex"))?false:true;
    packed_number=0; 
    next=0; 
    adj_value_set=false;
    expression = 0;
}

void CovariateTerm::add (const char* newname, int exp, const char* qualifier)
{
    CovariateTerm **e=end();
    *e=new CovariateTerm (newname, exp, qualifier);
}


/*
 * CovariateList is used to expand covariates entered
 * in shorthand forms.  It is only used during user operations
 * on covariates to present a nice interface.  It should not be
 * confused with CovariateTerm which is an integral part of Covariate.
 */
class CovariateList
{
    friend class Covariate;

    CovariateTerm *terms;
    CovariateList *next;

    CovariateList (CovariateList *cl, const char* var, int exp, 
		   const char* qualifier) {
	terms = new CovariateTerm (var, exp, qualifier);
	next = 0;
	if (cl) cl->next = this;
    }
    void add (const char *var, int exp, const char* qualifier) 
	{terms->add (var, exp, qualifier);}
    const char *string (char *buf);
    static void parse_or_throw (CovariateList **cl, char* user_name_string,
				char** user_qualifier_string_ref);
public:
    static CovariateList* parse (const char *user_string);
    void cleanup () {delete terms;}
    ~CovariateList () {
	if (next) delete next;
    }
};

double Covariate::_minrange = 1.0;
enum covariate_action {CREATE, SUSPEND, RESTORE, DELETE};
Covariate *Covariate::Covariates[MAX_PARAMETERS] = {0};
int Covariate::total_count = 0;
static int VirtualTraits = 0;


extern "C" int CovariateCmd (ClientData clientData, Tcl_Interp *interp,
		  int argc, char *argv[])
{
    if (argc == 2 && !Strcmp (argv[1], "-betanames"))
    {
	Covariate::show_betanames (interp);
	return TCL_OK;
    }
    if (argc == 2 && !Strcmp (argv[1], "-applicable"))
    {
	Covariate::show_applicable (interp, false);
	return TCL_OK;
    }
    if (argc == 2 && !Strcmp (argv[1], "-active"))
    {
	Covariate::show_applicable (interp, true);
	return TCL_OK;
    }

    if (argc == 2 && !StringCmp (argv[1], "help", case_ins))
    {
	return Solar_Eval (interp, "help covariate");
    }

    if (argc == 1)
    {
	Covariate::show_all (interp);
	return TCL_OK;
    }

    if (argc == 2 && !Strcmp (argv[1], "delete_all"))
    {
	Covariate::reset();
	return TCL_OK;
    }

    int first_arg = 1;

    covariate_action action = CREATE;
    if (argc >= 3 && !Strcmp (argv[1], "suspend"))
    {
	action = SUSPEND; 
	first_arg = 2;
    }
    else if (argc >= 3 && !Strcmp (argv[1], "restore"))
    {
	action = RESTORE;
	first_arg = 2;
    }
    else if (argc >= 3 && (!Strcmp (argv[1], "delete") ||
	                   !Strcmp (argv[1], "delete_fully")))
    {
	action = DELETE;
	first_arg = 2;
    }

// Operation also allowed to follow identifier if 3 arguments

    else if (argc == 3 && !Strcmp (argv[2], "suspend"))
    {
	action = SUSPEND; 
	first_arg = 1;
	argc = 2;
    }
    else if (argc == 3 && !Strcmp (argv[2], "restore"))
    {
	action = RESTORE;
	first_arg = 1;
	argc = 2;
    }
    else if (argc == 3 && (!Strcmp (argv[2], "delete") ||
	                   !Strcmp (argv[2], "delete_fully")))
    {
	action = DELETE;
	first_arg = 1;
	argc = 2;
    }

// If no explicit operator, it's "create"

    else
    {
	action = CREATE;
	first_arg = 1;
    }


// Perform action on list of covariate(s)

    int i;
    for (i = first_arg; i < argc; i++)
    {
	CovariateList *cl = 0;
	try
	{
	    cl = CovariateList::parse (argv[i]);
	    switch (action)
	    {
	    case CREATE:
		Covariate::add (cl);
		break;
	    case SUSPEND:
		Covariate::suspend (cl);
		break;
	    case RESTORE:
		Covariate::restore (cl);
		break;
	    case DELETE:
		Covariate::remove (cl);
	    }
	}
	catch (Safe_Error_Return& ser)
	{
	    if (cl) 
	    {
		if (action != CREATE) cl->cleanup();
		delete cl;
	    }
	    char *message = ser.message();
	    RESULT_BUF (message);
	    return TCL_ERROR;
	}
	if (cl) 
	{
	    if (action != CREATE) cl->cleanup();
	    delete cl;
	}
    }
    return TCL_OK;	
}

static char *make_betaname (char *bname, CovariateTerm *terms)
{
    strcpy (bname, "b");
    terms->fullname (&bname[1]);
    return bname;
}

Covariate::Covariate (CovariateTerm *terms)
{
    _terms = terms;
    _number = total_count++;
    _active = true;
    _mu_string = 0;
    for (int t = 0; t < MAX_TRAITS; t++)
    {
	_beta[t] = 0; 
	_save_beta[t] = 0.;
    }
     Covariates[_number] = this;
}


Covariate::~Covariate()
{
    for (int t = 0; t < Trait::Number_Of(); t++)
    {
	if (_beta[t])
	{	
	    _beta[t]->covariate = 0;
	    _beta[t]->remove ();
	}
    }
    for (int i = _number+1; i < total_count; i++)
    {
	Covariates[i-1] = Covariates[i];
	Covariates[i-1]->_number--;
    }
    Covariates[--total_count] = 0;
    if (_mu_string) free (_mu_string);
    delete _terms;
}

void Covariate::new_traits()
{
// Delete old beta parameters
    int i,t;
    for (i = 0; i < total_count; i++)
    {
	Covariate* c = Covariates[i];
	for (t = 0; t < last_number_of_traits; t++)
	{
	    if (c->_beta[t])
	    {
		c->_beta[t]->covariate = 0;
		c->_beta[t]->remove ();
	    }
	}

// Check for existing new beta parameters to be used

	for (t = 0; t < Trait::Number_Of(); t++)
	{

// A trait qualifier attached to any term affects them all
// (and the first one overrides...)

	    char* usetrait = 0;
	    CovariateTerm* tt = c->_terms;
	    while (tt)
	    {
		if (tt->trait)
		{
		    usetrait = tt->trait;
		    break;
		}
		tt = tt->next;
	    }
	    if (usetrait && Strcmp (usetrait, Trait::Name(t)))
	    {
		c->_beta[t] = 0;
		continue;
	    }

	    char bname[1024];
	    make_betaname (bname, c->_terms);
	    if (Trait::Number_Of() > 1 && !usetrait)
	    {
		strcat (bname, "(");
		strcat (bname, Trait::Name(t));
		strcat (bname, ")");
	    }
	    c->_beta[t] = Parameter::find (bname);
	    if (!c->_beta[t])
	    {

// Make new beta parameters

		c->_beta[t] = Parameter::add (bname);
	    }
	    c->_beta[t]->covariate = c;
	    c->_beta[t]->trait = t;
	}
    }
    last_number_of_traits = Trait::Number_Of();
}


Covariate *Covariate::find (CovariateTerm *terms)
{
    Covariate *c;
    int i;
    for (i = 0; c = Covariate::index(i); i++)
    {
	if (CovariateTerm::compare (c->_terms, terms))
	{
	    return c;
	}
    }
    return 0;
}


void Covariate::reset ()
{
    while (total_count)
    {
	Covariate *c = Covariate::index (0);
	delete c;
    }
}

void Covariate::remove (CovariateList *cl)
{
    do
    {
	Covariate *c = Covariate::find (cl->terms);
	if (c==0) 
	{
	    char buf[1024];
	    char nbuf[512];
	    sprintf (buf, "Covariate %s not found", 
		     cl->terms->fullname (nbuf));
	    throw Safe_Error_Return (buf);
	}
	delete c;

    } while (0 != (cl = cl->next));
}

void Covariate::zero_all ()
{
    int i, t;
    Covariate *c;
    for (i = 0; c = Covariate::index (i); i++)
    {
	for (t = 0; t < Trait::Number_Of(); t++)
	{
	    if (c->_beta[t]) c->_beta[t]->zero();
	}
    }
}


// Original implementation of fullname had memory bug, fixed in 4.0.4.
// Intent is to put covariate trait in any term at the end, with
// the first having precedence.  Most straightforward implementation
// is recursive, with qualification (traitname) passing through second
// argument hidden from outside caller (see declaration in solar.h).
//
// parse_or_throw now ensures only one trait can be specified anyway

const char *CovariateTerm::fullname_imp (char *buf, char* traitname)
{
    if (!traitname) traitname = trait;

    strcpy (buf, name);
    if (exponent != 1)
    {
	sprintf (&buf[strlen(buf)], "^%d", exponent);
    }
    if (next)
    {
	strcat (buf, "*");
	next->fullname_imp(&buf[strlen(buf)], traitname);
    } 
    else if (traitname)
    {
	sprintf (&buf[strlen(buf)], "(%s)", traitname);
    }
    return buf;
}


CovariateList *CovariateList::parse (const char *string)
{
    CovariateList *cl = 0;
    char* string_copy = Strdup (string);
    char* qualifier_copy = Strdup ("");
    try
    {
	CovariateList::parse_or_throw (&cl, string_copy, &qualifier_copy);
    }
    catch (Safe_Error_Return& ser)
    {
	free (string_copy);
	free (qualifier_copy);
	if (cl)
	{
	    cl->cleanup ();
	    delete cl;
	}
	throw ser;
    }
    free (string_copy);
    free (qualifier_copy);
    return cl;
}

void CovariateList::parse_or_throw (CovariateList** first_cl, 
				    char* string, char** qualifier_ref)
{

// Find qualifier, if any, now
// If found, stash it, and remove from remainder of string
// Assuming validity, we look for first and last parens
//  (May be more parens inside, if contains trait expression)

// (This is a radical change from how it used to be done, made
//  necessary because of expression traits, but come to think of it,
//  it should always have been done like this...much more intuitive,
//  not that any of this procedure is intuitive.)

    char* qualifier = qualifier_ref[0];

    char* lparen;
    if (0 != (lparen = strchr (string, '(')))
    {
	char* ptr = lparen+1;
	char* rparen = 0;
	int lpcount = 1;
	while (*ptr != 0)
	{
	    if (*ptr == ')')
	    {
		lpcount--;
		if (!lpcount)
		{
		    rparen = ptr;
		    break;
		}
	    }
	    else if (*ptr == '(')
	    {
		lpcount++;
	    }
	    ptr++;
	}
	if (!rparen)
	{
	    throw Safe_Error_Return
		("Missing right parenthesis");
	}
	*lparen = *rparen = '\0';
	qualifier = StringAppend (&qualifier, &lparen[1]);
	*qualifier_ref = qualifier;

// Now that qualifier has been copied, cut it out

	while (*++ptr != '\0')
	{
	    *lparen++ = *ptr;
	}
	*lparen = *ptr;

// Now ensure that this is the only qualifier

	if (strchr (string, '('))
	{
	    throw Safe_Error_Return
		("More than one trait specified for covariate");
	}
    }
    else
    {
	qualifier = 0;
    }


// Repeat until no more interactors

    CovariateList *cl = 0;
	
    char *variable1 = string;
    char *pound1p = 0;
    int exponent1 = 1;
    int degree1 = 1;

    for (;;)
    {

// Set pointers to delimiters

	char *variable = string;
	char *interactor = 0;
	int exponent = 1;
	int degree = 1;

	char *poundp = strchr (string, '#');
	char *starp = strchr (string, '*');

// Find beginning of next term, if any
// Terminate this term there

	if (starp || poundp)
	{
	    if (starp && poundp)
	    {
		throw Safe_Error_Return 
		    ("Pound and Star are not allowed in same covariate term");
	    }
	    else if (starp)
	    {
		interactor = starp+1;
		*starp = '\0';
	    }
	    else
	    {
		interactor = poundp+1;
		*poundp = '\0';
	    }
	}

// Find exponent(s), if any, attached to this term
//  (and qualifier might follow exponent)

	char *caratp = strchr (string, '^');
	char *commap = strchr (string, ',');
	if (commap && !caratp)
	{
	    throw Safe_Error_Return
		("Invalid exponent expression or use of comma");
	}

	if (caratp)
	{

// Handle ^1,2 shortcut

	    if (commap)
	    {
		degree = 0;
		if (!strncmp (caratp, "^1,2", strlen ("^1,2")))
		{
		    degree = 2;
		}
		if (!strncmp (caratp, "^1,2,3", strlen ("^1,2,3")))
		{
		    degree = 3;
		}
		if (!strncmp (caratp, "^1,2,3,4", strlen ("^1,2,3,4")))
		{
		    degree = 4;
		}
	        if ('\0' != caratp[degree*2] || degree < 0)
		{
		    throw Safe_Error_Return 
			("Invalid exponent expression or use of comma");
		}
		*caratp = '\0';  // terminates variable name at ^
	    }
	    else
	    {
		char junk[1024];
		bool invalid = false;
		int count = sscanf (caratp+1, "%d %s", &exponent, junk);
		if (count==1)
		{
		    *caratp = '\0';
		}
		else
		{
		    throw Safe_Error_Return ("Invalid exponent expression");
		}
	    }  // end else (!commap)
	} // end if (caratp)

	if (!cl)
	{
	    if (0==strlen(variable)) throw Safe_Error_Return (
		"Missing covariate variable name");
	    cl = *first_cl = new CovariateList (0, variable, exponent, 
						qualifier);
	    variable1 = variable;
	    exponent1 = exponent;
	    degree1 = degree;
	    pound1p = poundp;
	    if (!interactor && degree>1) // variable^1,2....
	    {
		int i;
		for (i = 2; i <= degree; i++)
		{
		    cl = new CovariateList (cl, variable, i, qualifier);
		}
	    }
	}
	else
	{

// If shortcuts, handle them here (exiting afterwards)
// Shortcuts are only allowed for 1 or 2 variables
// Check that there is not another "interactor"

	    if (0==strlen(variable)) throw Safe_Error_Return (
		"Missing covariate variable name");
	    if (commap)
	    {
		throw Safe_Error_Return (
		   "Covariate comma shortcut only allowed for first variable");
	    }
	    if (interactor)
	    {
		if (degree1>1 || pound1p || poundp)
		{
		    throw Safe_Error_Return (
			"Covariate shortcuts only allowed with 2 variables");
		}
	    }

	    if (pound1p && degree1==1)  // variable1[^n1]#variable2[^n2]
	    {
		cl = new CovariateList (cl, variable, exponent,
					qualifier);  // var2
		cl = new CovariateList (cl, variable1, exponent1,
					qualifier); // var1*var2
		cl->add (variable, exponent, qualifier);
	    }
	    else if (degree1>1 && !pound1p)  // variable1^1,2,...*variable2^n2
	    {
		cl->add (variable, exponent, qualifier);
		for (int i = 2; i <= degree1; i++)
		{
		    cl = new CovariateList (cl, variable1, i, qualifier);
		    cl->add (variable, exponent, qualifier);
		}
	    }
	    else if (degree1>1 && pound1p)  // variable1^1,2,...#variable2^n2
	    {
		cl = new CovariateList (cl, variable, exponent,
					qualifier);  //var2
		cl = new CovariateList (cl, variable1, exponent1,
					qualifier); //var1*var2
		cl->add (variable, exponent, qualifier);

		for (int i = 2; i <= degree1; i++)
		{
		    cl = new CovariateList (cl, variable1, i,
					    qualifier);
		    cl = new CovariateList (cl, variable1, i,
					    qualifier);
		    cl->add (variable, exponent, qualifier);
		}
	    }
	    else  // The non-shortcut case of a*b
	    {
		cl->add (variable, exponent, qualifier);  // Now a*b is completed
	    }
	}

// We repeat if there is another "interactor"
	
	string = interactor;
	if (!string) return;
    }
}

	
void Covariate::add (CovariateList *cl)
{
    int t;
    int ntraits = Trait::Number_Of ();
    Parameter* bp[MAX_TRAITS];

    do
    {
	CovariateTerm *terms = cl->terms;

	if (Covariate::find (terms)) continue;

// Make (or find) beta parameter for each trait
// And check for use by another covariate

	for (t=0; t < ntraits; t++)
	{

// A trait qualifier attached to any term affects them all
// (and the first one overrides...)

	    char* usetrait = 0;
	    CovariateTerm* tt = terms;
	    while (tt)
	    {
		if (tt->trait)
		{
		    usetrait = tt->trait;
		    break;
		}
		tt = tt->next;
	    }
	    if (usetrait && Strcmp (usetrait, Trait::Name(t)))
	    {
		bp[t] = 0;
		continue;
	    }

	    char bname[1024];
	    make_betaname (bname, terms);
	    if (Trait::Number_Of() > 1 && !usetrait)
	    {
		strcat (bname, "(");
		strcat (bname, Trait::Name(t));
		strcat (bname, ")");
	    }
	    bp[t] = Parameter::find (bname);

// Look for old-style parameter name with :

	    if (!bp[t] && strchr (bname, '*'))
	    {
		char bname2[1024];
		strcpy (bname2, bname);
		char *ptr;
		while (ptr = strchr (bname2, '*'))
		{
		    *ptr = ':';
		}
		bp[t] = Parameter::find (bname2);
	    }
	
// If we didn't find old or new, make a new one

	    if (!bp[t])
	    {
		bp[t] = Parameter::add (bname);
	    }

	    if (bp[t]->covariate)
	    {
		char buf[1024];
		sprintf (buf, 
           "Unable to use parameter '%s' for two different covariates", 
			 bname);
		cl->cleanup ();  // Only cleanup unused terms;
		throw Safe_Error_Return (buf);
	    }
	}

// We've created beta parameters for each trait.  
// Now create covariate and bind to beta parameter.

	Covariate *c = new Covariate (terms);
	for (t = 0; t < ntraits; t++)
	{
	    c->_beta[t] = bp[t];
	    if (bp[t])
	    {
		bp[t]->covariate = c;
		bp[t]->trait = t;
	    }
	}

    } while (0 != (cl = cl->next));
    last_number_of_traits = Trait::Number_Of();
}

// The mu-string is now generated only when required to allow for scaling.
// However, the string is saved for deletion when covariate itself goes away
//    (Very lazy garbage collection)

const char* Covariate::mu_string (int trait)
{
// Free last mu string (if any)
    if (_mu_string) free (_mu_string);

    const char* bname;
    bool found = false;
    bool some_not_found = false;
    if (!_beta[trait])
    {
	_mu_string = Strdup ("");
	return _mu_string;
    } else {
	bname = _beta[trait]->name();
    }

    char addstring[1024];
    addstring[0] = '\0';
    
// Surround betaname with <> if betaname contains * or + or ()

    if (0!=strchr(bname,'*') || 0!=strchr(bname,'^') || 
	0!=strchr(bname,')'))
    {
	sprintf (addstring, "%s+<%s>", addstring, bname);
    }
    else
    {
	sprintf (addstring, "%s+%s", addstring, bname);
    }


// Write each variable term, subtracting appropriate adjustment value
    CovariateTerm* term = _terms;
    for ( ; term; term = term->next)
    {

	if (term->sex)
	{
	    sprintf (addstring, "%s*Female", addstring);
	}
	else
	{
	    char expstr[32];
	    strcpy (expstr, "");
	    if (term->exponent > 1)
	    {
		sprintf (expstr, "^%d", term->exponent);
	    }
	    Scale* s = Scale::Find (term->name);
	    if (s)
	    {
		if (s->adj() == 0.0)
		{
		    sprintf (addstring, "%s*%s%s", addstring,
			     term->name, expstr);
		} else {
		    sprintf (addstring, "%s*(%s-%-.10g)%s", addstring,
			     term->name, s->adj(), expstr);
		}
	    } else {
		char testchars[5];
		strncpy (testchars, term->name, 4);
		testchars[4] = '\0';
		if (!Strcmp ("snp_",testchars) || !Strcmp("hap_",testchars))
		{
		    sprintf (addstring, "%s*%s%s", addstring,
			     term->name, expstr);
		}
		else if (term->adj_value_set)
		{
		    if (term->adj_value == 0.0)
		    {
			sprintf (addstring, "%s*%s%s", addstring,
				 term->name, expstr);
		    } else {
			sprintf (addstring, "%s*(%s-%-.10g)%s", addstring,
				 term->name, term->adj_value, expstr);
		    }
		} else {
		    sprintf (addstring, "%s*(%s-x_%s)%s", addstring,
			     term->name, term->name, expstr);
		}
	    }
	}
    }

// Add final bracket if trait-qualified covariate

// Save as part of covariate and return
    _mu_string = Strdup (addstring);
    return _mu_string;
}
	

void Covariate::suspend (CovariateList *cl)
{
    do
    {
	CovariateTerm *term = cl->terms;

	Covariate *c = Covariate::find (term);
	if (!c) 
	{
	    char buf[1024];
	    char nbuf[512];
	    sprintf (buf, "Covariate %s not found", cl->terms->fullname(nbuf));
	    throw Safe_Error_Return (buf);
	}
	c->_active = false;
	for (int i = 0; i < Trait::Number_Of(); i++)
	{
	    if (c->_beta[i])
	    {
		c->_save_beta[i] = c->_beta[i]->start;
		c->_beta[i]->start = 0.;
	    }
	}
    }
    while (0 != (cl = cl->next));
}


void Covariate::restore (CovariateList *cl)
{
    do
    {
	CovariateTerm *term = cl->terms;

	Covariate *c = Covariate::find (term);
	if (!c) 
	{
	    char buf[1024];
	    char nbuf[512];
	    sprintf (buf, "Covariate %s not found", cl->terms->fullname(nbuf));
	    throw Safe_Error_Return (buf);
	}
	c->_active = true;
	for (int t = 0; t < Trait::Number_Of(); t++)
	{
	    if (c->_beta[t])
	    {
		c->_beta[t]->start = c->_save_beta[t];
	    }
	}
    }
    while (0 != (cl = cl->next));
}


void Covariate::write_commands (FILE *file)
{
    Covariate *c;
    for (int i = 0; c = Covariate::index (i); i++)
    {
	char nbuf[512];
	fprintf (file, "covariate %s\n", c->fullname (nbuf));
	if (!c->active())
	{
	    fprintf (file, "covariate suspend %s\n", c->fullname (nbuf));
	}
    }
}

void Covariate::show_all (Tcl_Interp *interp)
{
    Covariate *c;
    for (int i = 0; (c = Covariate::index(i)); i++)
    {
	char buf[512];
	size_t first_pos = 0;
	if (!c->active())
	{
	    sprintf (buf, "Suspended[");
	    first_pos = strlen (buf);
	}
	c->fullname (&buf[first_pos]);
	if (!c->active())
	{
	    strcat (buf, "]");
	}
	Solar_AppendElement (interp, buf);
    }
}

// Returns list of covariates which are applicable to current trait(s)
// optionally only those active also

void Covariate::show_applicable (Tcl_Interp *interp, bool only_active)
{
    Covariate *c;
    bool betanames = false;
    for (int i = 0; c = Covariate::index(i); i++)
    {
	char buf[1024];
	if (!only_active || c->active())
	{
	    bool found_beta = false;
	    for (int t = 0; t < Trait::Number_Of(); t++)
	    {
		if (c->beta(t))
		{
		    if (betanames)
		    {
			strcpy (buf, c->beta(t)->name());
		    }
		    else if (!found_beta)
		    {
			size_t first_pos = 0;
			if (!c->active())
			{
			    sprintf (buf, "Suspended[");
			    first_pos = strlen (buf);
			}
			c->fullname (&buf[first_pos]);
			if (!c->active())
			{
			    strcat (buf, "]");
			}
		    }
		    Solar_AppendElement (interp, buf);
		    found_beta = true;
		}
	    }
	}
    }
}


void Covariate::show_betanames (Tcl_Interp *interp)
{
    Covariate *c;
    for (int i = 0; (c = Covariate::index(i)); i++)
    {
	char buf[512];
	for (int t = 0; t < Trait::Number_Of(); t++)
	{
	    Parameter *beta = c->beta(t);
	    if (!beta) continue;
	    Solar_AppendElement (interp, (char*) beta->name());
	}
    }
}

double eval_blank (int eqindex)
{
    return MISSING_PHENOTYPE;
}

int Covariate::bind (Tcl_Interp *interp)
{
    if (_Context) delete _Context;
    _Context = 0;

// Check for VirtualTraits

    VirtualTraits = 0;
    if (-1 != Option::get_int ("UnbalancedTraits") && Trait::Number_Of() > 1)
    {
	VirtualTraits = Trait::Number_Of();
    }

// Set up indexes to each mean
// iff default part of mu is used (and Covariate::eval will be called)

    int i;
    if (Mu::default_part_used())
    {
	for (i = 0; i < Trait::Number_Of(); i++)
	{
	    if (-1 == (Mean[i] = Trait::Find_Parameter (i, "mean")))
	    {
		char buf[1024];
		if (Trait::Number_Of() > 1)
		{
		    sprintf (buf, "Missing parameter mean(%s)", 
			     Trait::Name(i));
		}
		else
		{
		    sprintf (buf, "Missing parameter mean");
		}
		RESULT_BUF (buf);
		return TCL_ERROR;
	    }
	}
    }

// Make sure every variable required is present

    Covariate *c = 0;
    for (i = 0; (c = Covariate::index(i)); i++)
    {
	CovariateTerm *term;
	for (term = c->_terms; term; term = term->next)
	{
	    if (term->expression) delete term->expression;
	    if (!term->sex)
	    {
		bool use_phen = Phenotypes::available (term->name);
		Definition* definition = Definition::Find (term->name);
		if (!use_phen && !definition)
		{
		    char buf[256];
		    sprintf (buf, "Missing covariate variable %s",
			     term->name);
		    RESULT_BUF (buf);
		    return TCL_ERROR;
		}
		if (use_phen && definition)
		{
		    char buf[256];
		    sprintf (buf, 
		      "Covariate variable '%s' both phenotype and definition",
			     term->name);
		    RESULT_BUF (buf);
		    return TCL_ERROR;
		}

// Defined traits must have context established now
// which will require adding additional phenotypes and possibly an inormal

		if (definition)
		{
		    try
		    {
			term->expression = new Expression 
			    (definition->definition());
		    }
		    catch (...)
		    {
			char* errmes = Strdup 
	         ("Expression error in define command for covariate term ");
			StringAppend (&errmes, term->name);
			RESULT_BUF (errmes);
			free (errmes);
			return TCL_ERROR;
		    }
		    if (!_Context) _Context = new Context;
		    _Context->add ("blank", eval_blank, 0);
		    bool ok = false;
		    char* inor_match_class (char* string, bool *ifclass, int* pclass);
		    char* zscor_match (char* string);
		    int pclass;

		    while (!ok)
		    {
			bool ifclass = true; // check for class
			try
			{
			    term->expression->bind (_Context);
			    ok = true;
			}
			catch (Undefined_Name& un)
			{
			    char *pname;
			    try
			    {
				if (Phenotypes::available (un.name))
				{
				    int vindex = -1;
				    if (term->trait && 
					Strcmp (term->trait, ""))
				    {
					int ttt;
					for (ttt = 0; 
					     ttt < Trait::Number_Of();
					     ttt++)
					{
					    if (!Strcmp (Trait::Name(ttt),
							 term->trait))
					    {
						vindex = EqVar::add_for_trait 
						    (un.name, ttt);
					    }
					}
				    }
				    else
				    {
					vindex = EqVar::add (un.name);
				    }

//        Assumption is that if covariate is trait specific, and that
//        trait is not present, the covariate will never be accessed anyway
//        so the fact that pretend to we set it up here, without actually
//        having set up the EqVar for it is unimportant because the
//        eval_phenotype will never get called.
//
//        This assumption is checked for in eval_phenotype.  If violated,
//        program termination is forced and error messages advises
//        contacting solar@txbiomedgenetics.org

				    _Context->add (un.name, eval_phenotype, 
						   vindex);

// Incomplete...this might not work if proband indicator

				}
				else if (!Strcmp ("sex", un.name))
				{
				    _Context->add (un.name, eval_sex, 0);
				}
				else if (0 != (pname = zscor_match (un.name)))
				{
				    if (!Phenotypes::available (pname))
				    {
					throw Undefined_Name (pname);
				    }

// zscorexp will either return pre-stored Mean and SD
// or it will return error code indicating no stored Mean,SD for this trait

				    bool meanvalid = false;
				    double mean, sd;
				    char cbuf[1024];
				    char junk[1024];
				    sprintf (cbuf, "zscorexp get %s mean",pname);
				    int errstatus = Solar_Eval (interp, cbuf);
				    if (errstatus)
				    {
					RESULT_LIT ("");
				    }
				    else
				    {
					int scount = sscanf (interp->result,
							     "%lg %s", &mean, junk);
					if (scount == 1)
					{
					    sprintf (cbuf, "zscorexp get %s sd",pname);
					    errstatus = Solar_Eval (interp, cbuf);
					    if (errstatus)
					    {
						RESULT_LIT ("");
					    }
					    else
					    {
						scount = sscanf (interp->result,
								 "%lg %s", &sd, junk);
						if (scount == 1)
						{
						    meanvalid = true;
						}
					    }
					}
				    }
			    // vindex is index to underlying phenotype
			    // zindex is index to zscore object
				    int vindex = EqVar::add_for_trait (pname, i);
				    int zindex;
				    if (!meanvalid)
				    {
					zindex = Zscore::Setup (pname, vindex);
				    }
				    else
				    {
					zindex = Zscore::SetupMean (pname, mean, sd,
								    vindex);
				    }
				    _Context->add (un.name, get_zscore, zindex);
				}
			        else if (0 != (pname = inor_match_class (un.name,
								 &ifclass, &pclass)))
				{
				    if (!Phenotypes::available (pname))
				    {
					throw Undefined_Name (pname);
				    }
			    
				    char cbuf[1024];
				    if (ifclass) {
					sprintf (cbuf, 
					 "inormal -trait %s -phenfile -class %d",
						 pname, pclass);
				    } else {
					sprintf (cbuf, 
						 "inormal -trait %s -phenfile", 
						 pname);
				    }

				    int errstatus = Solar_Eval (interp, cbuf);
				    if (errstatus)
				    {
					printf (
		  "maximize: error applying inormal to covarite %s\n", pname);
                                 // Error message from inormal in TCL_Result
					return TCL_ERROR;   
				    }
			    
				    if (!Phenotypes::INCount)
				    {
					Phenotypes::INData = (double*)
					    Calloc (2, sizeof (double));
					Phenotypes::INNames = (char**)
					    Calloc (2, sizeof (char*));
					Phenotypes::INIfclass = (bool*) Calloc (2, sizeof (int));
					Phenotypes::INClass = (int*) Calloc (2, sizeof (bool));
				    }
				    else
				    {
					int newsize = Phenotypes::INCount + 2;
					Phenotypes::INData = (double*) Realloc
					    (Phenotypes::INData, 
					     newsize*sizeof (double));
					Phenotypes::INNames = (char**) Realloc
					    (Phenotypes::INNames, 
					     newsize*sizeof (char*));
					Phenotypes::INIfclass = (bool*) Realloc (Phenotypes::INIfclass, newsize*sizeof(bool));
					Phenotypes::INClass = (int*) Realloc (
					    Phenotypes::INClass, newsize*sizeof(int));
				    }
				    Phenotypes::INNames[Phenotypes::INCount]
					= Strdup (pname);
				    Phenotypes::INIfclass[Phenotypes::INCount] = ifclass;
				    Phenotypes::INClass[Phenotypes::INCount] = pclass;
				    double get_inormal (int inindex);
				    _Context->add (un.name, get_inormal, Phenotypes::INCount);
				    Phenotypes::INCount++;
				} else {
				    throw Undefined_Name (un.name);
				}
			    }
			    catch (Undefined_Name& un)
			    {
				char buf[1024];
				sprintf (buf, 
		   "Undefined variable `%s' in covariate definition", un.name);
				RESULT_BUF (buf);
				return TCL_ERROR;
			    }
			}
		    }
		}
	    }
	}
    }
    return TCL_OK;
}

/*
 * FUNCTION INITCOV is called by SUBROUTINE INITAL to initialize covariates
 * It sets up the par (start), parmin (lower) and parmax (upper) for
 * this parameter.
 */
extern "C" void initcov_ (int *PNUM, double *PAR, double *parmin,
			 double *parmax,
			 double *vvar, double *vmin, double *vmax)
{
    Covariate::set_boundaries (PNUM, PAR, parmin,
			       parmax, vvar, 
			       vmin, vmax);
}

int Covariate::set_boundaries (int *PNUM, double *PAR, double *parmin,
			       double *parmax,
			       double *vvar, double *vmin, double *vmax)
{
    double spread_factor = Option::get_double ("AutoCovarBound");  // default 2

// Estimate trait factors here based on assumption of 1 trait
// reset later to correct trait values if trait known

    double deltaT = vmax[0] - vmin[0];
    double mulT = spread_factor * deltaT;
//  printf ("mulT is %g\n",mulT);

    int pnum = *PNUM - 1;
    Parameter *par = Parameter::index (pnum);
    Covariate *cov = par->covariate;

    if (!cov)
    {
// Default determined by variable with minimum range
	*parmax = mulT / Covariate::_minrange;
    } 
    else
    {
	int trait = par->trait;
	deltaT = vmax[trait] - vmin[trait];
	mulT = spread_factor * deltaT;
	if (cov->_terms->sex==true && cov->_terms->next==0)
	{
	    *parmax = mulT;
	}
	else
	{
	    double sum_of_squares = 0.0L;
	    for (CovariateTerm *term = cov->_terms; term; term = term->next)
	    {
		if (term->sex) continue;
		int vindex = term->packed_number;
		double deltaV = vmax[vindex] - vmin[vindex];
		if (deltaV == 0.0) continue;
		double vterm = deltaV;
		for (int i = 2; i <= term->exponent; i++)
		{
		    vterm *= deltaV;
		}
		sum_of_squares += vterm * vterm;
	    }
//	    printf ("sum_of_squares is %g\n",sum_of_squares);
	    double divisor = sqrt (sum_of_squares);
	    if (divisor == 0.0L) divisor = 1.0L;
//	    printf ("mulT/divisor = %g / %g\n",mulT, divisor);
	    *parmax = mulT / divisor;
//	    printf ("parmax is %g\n",parmax);
	}
    }
    if (*parmax < 2.0E-4) *parmax = 2.0E-4;
    *parmin = - *parmax;
    if (!cov) return false;
    return true;
}

//
// Fisher never bothered to save means for each variable
// Since we need this info, we do it here
//
extern "C" int savemeans_ (int *nvar, double *vmean, double *vmin, 
			    double *vmax, int *output_unit)
{
    return Covariate::savemeans (nvar, vmean, vmin, vmax, output_unit);
}

extern "C" void _center_string (char*, char*);

extern int *Vbinary;

int Covariate::savemeans  (int *nvar, double *vmean, double *vmin, 
			    double *vmax, int *output_unit)
{

// This is where we set adj_value for each covariate variable

    int i;
    Covariate *cov;
    for (i = 0; (cov = Covariate::index (i)); i++)
    {
	if (!cov->active ()) continue;
	CovariateTerm *term;
	for (term = cov->_terms; term; term = term->next)
	{
	    if (term->sex) continue;  // Don't need adj value for sex
	    Scale* s = Scale::Find (term->name);
	    if (s)
	    {
		term->adj_value = s->adj ();
	    }
	    else
	    {
		char testchars[5];
		strncpy (testchars, term->name, 4);
		testchars[4] = '\0';
		if (!Strcmp ("snp_",testchars) || !Strcmp("hap_",testchars))
		{
		    term->adj_value = 0;
		}
		else 
		{
		    int vindex = term->packed_number;
		    term->adj_value = vmean[vindex];
		    if (Vbinary[vindex])
		    {
			// binaries adj'd to lower val
			term->adj_value = vmin[vindex]; 
		    }
		}
	    }
	    term->adj_value_set = true;
	}
    }

// Now, find the variable with the smallest non-zero range
//   (This is used to scale beta values w/o known variable associations)

    Covariate::_minrange = 1; // Could be 'sex' which has range 1

    for (i = 1; i < *nvar; i++)
    {
	double delta = fabs (vmax[i] - vmin[i]);
	if (delta > 0 && delta < Covariate::_minrange) 
	    Covariate::_minrange = delta;
    }
    return 0;
}

static int evalcount = 0;

double Covariate::eval (int *malep, double *var, double *par, int trait)
{
    double meancov;
    if (VirtualTraits)
    {
	trait = (int) var[0] - 1;
    }
    meancov = par[Mean[trait]];


#if 0
    int id = var[2];
    fprintf (stderr,"ID:%d  Trait:%d  Mean:%f\n", id, trait, meancov);
    if (++evalcount > 120) exit (0);
#endif

    Covariate *c;
    for (int i = 0; c = Covariate::index(i); i++)
    {
	if (!c->_active) continue;

	CovariateTerm *term;
	double factor = 1.0;
	bool not_this_trait = false;
	for (term = c->_terms; term; term = term->next)
	{
	    if (c->_beta[trait] == 0)
	    {
		not_this_trait = true;
		break;
	    }
	    if (term->sex)
	    {
		if (*malep!= 0)
		{
		    factor = 0.0;
		    break;
		}
	    }
	    else
	    {
		for (int j = 0; j < term->exponent; j++)
		{
		    factor *= var[term->packed_number] - term->adj_value;
		}
	    }
	}
	if (!not_this_trait)
	{
	    meancov += factor * par[c->_beta[trait]->number()];
	}
    }
    return meancov;
}

int Covariate::boundary_check (double factor, char *pnames, int plen)
{
    *pnames = '\0';

// Check bounds for each active covariate.
// Adjust if factor != 0

    if (factor>1.0) factor -= 1;

    bool adjusted = false;
    Covariate *c;
    for (int i = 0; (c = Covariate::index(i)); i++)
    {
	if (!c->active ()) continue;
	for (int t = 0; t < Trait::Number_Of(); t++)
	{
	    Parameter *beta = c->beta(t);
	    if (!beta) continue;

	    double current = beta->start;
	    double upper = beta->upper;
	    double lower = beta->lower;
	    bool fixup = beta->fix_upper;
	    bool fixlow = beta->fix_lower;


	    double uldelta = upper - lower;
	    double dupper = fabs (upper - current);
	    double dlower = fabs (current - lower);
	    
	    double ddelta = 0.001;
	    if (ddelta > uldelta/1000.)
	    {
		ddelta = uldelta/1000.;
	    }

// If value is within 0.001 (or scaled value) of boundary, adjust

	    if (dupper < ddelta && !fixup)
	    {
		beta->upper += factor * uldelta;
		adjusted = true;
		char buf[256];
		if (Verbosity::max() && factor > 0.0)
		{
		    printf ("    *** Adjusting %s upper bound by %g to %g\n", 
			    c->fullname (buf), factor+1.0, beta->upper);
		}
		if (strlen(c->fullname(buf))+1+1+strlen(pnames) < plen)
		{
		    strcat (pnames, c->fullname(buf));
		    strcat (pnames, " ");
		}
	    }
	    if (dlower < ddelta && !fixlow)
	    {
		beta->lower -= factor * uldelta;
		adjusted = true;
		char buf[256];
		if (Verbosity::max() && factor > 0.0)
		{
		    printf ("    *** Adjusting %s lower bound by %g to %g\n", 
			    c->fullname (buf), factor+1.0, beta->lower);
		}
		if (strlen(c->fullname(buf))+1+1+strlen(pnames) < plen)
		{
		    strcat (pnames, c->fullname(buf));
		    strcat (pnames, " ");
		}
	    }
	}
    }
    return adjusted?1:0;
}

CovariateTerm::~CovariateTerm ()
{
    if(next) delete next; 
    free (name); 
    if (trait) free (trait);
}

