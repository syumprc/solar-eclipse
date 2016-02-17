/*
 * constraint.cc implements the Constraint class and command
 * Written by Charles Peterson beginning in November, 1997
 * Copyright (c) 1997 Southwest Foundation for Biomedical Research
 */

#include <stdio.h>
#include <ctype.h>
// tcl.h from solar.h
#include "solar.h"

Constraint *Constraint::Constraints[] = {0};
int Constraint::_count = 0;

extern bool Maximize_InitPar;

extern "C" int ConstraintCmd (ClientData clientData, Tcl_Interp *interp,
		  int argc, char *argv[])
{
    if (argc == 1)
    {
	Constraint::display_all ();
	return TCL_OK;
    }

// Now we have 2 or more arguments

    if (argc == 2 && !StringCmp (argv[1], "command", case_ins))
    {
	Constraint::return_commands (interp);
	return TCL_OK;
    }

    if (argc == 2 && !StringCmp (argv[1], "help", case_ins))
    {
	return Solar_Eval (interp, "help constraint");
    }

    const char* message = 0;
    if (argc >= 3 && (!Strcmp (argv[1], "delete") ||
                      !Strcmp (argv[1], "delete_fully") ||
                      !Strcmp (argv[1], "del") ||
	              !Strcmp (argv[2], "delete") ||
	              !Strcmp (argv[2], "del") ||
	              !Strcmp (argv[2], "delete_fully")))
    {
	char *cnumstr = argv[2];
	if (!Strcmp (argv[2], "delete") ||
	    !Strcmp (argv[2], "del") ||
	    !Strcmp (argv[2], "delete_fully"))
	{
	    cnumstr = argv[1];
	}
	int number;
	char extra[128];
	int count = sscanf (cnumstr, "%d %s", &number, &extra);
	if (count != 1 || argc > 3)
	{
	    if (argc > 3)
	    {
		char* con_string = Tcl_Concat (argc-2, &argv[2]);
		message = Constraint::remove_spec (con_string);
		Tcl_Free (con_string);
	    }
	    else
	    {
		message = Constraint::remove_spec (cnumstr);
	    }
	}
	else
	{
	    message = Constraint::remove (number);

	}
	if (message)
	{
	    RESULT_LIT (message);
	    return TCL_ERROR;
	}
	return TCL_OK;
    }

    if (argc == 2 && !StringCmp (argv[1], "delete_all", case_ins))
    {
	Constraint::reset();
	return TCL_OK;
    }

// This must be a new constraint (or error)

    char *con_string = Tcl_Concat (argc-1, &argv[1]);
    message = Constraint::parse_and_add (con_string);
    Tcl_Free (con_string);
    if (message)
    {
	RESULT_LIT (message);
	return TCL_ERROR;
    }
    return TCL_OK;
}

const char* Constraint::parse_and_add (const char* cstring)
{
    const char* message;
    Constraint* c = Constraint::parse (cstring, &message);
    if (!c)
    {
	return message;
    }
    return c->add ();
}

void Constraint::display_all ()
{
    for (int i = 0; i < _count; i++)
    {
	printf ("[%d] %s\n", i+1, Constraints[i]->string);
    }
}


int Constraint::return_commands (Tcl_Interp *interp)
{
    for (int i = 0; i < _count; i++)
    {
	char buf[1024];
	sprintf (buf, "constraint %s", Constraints[i]->string);
	Solar_AppendElement (interp, buf);
    }
    return TCL_OK;
}

void Constraint::reset()
{
    while (_count)
    {
	remove ( _count );
    }
}


const char* Constraint::remove (int n)
{
    n--;  // make 0-based
    if (n < 0 || n >= _count) return "No such constraint";
    delete Constraints[n];
    for (n += 1; n < _count; n++)
    {
	Constraints[n-1] = Constraints[n];
    }
    _count--;
    return 0;
}

void Constraint::remove_if_needs_par (const char* pname)
{
    int i = 0;
    while (i < _count)
    {
	bool remove_this = false;
	Constraint* c = Constraints[i];
	for (int j = 0; j < c->Nlterms; j++)
	{
	    Term* t = c->lterms[j];
	    const char* tname = t->name();
	    if (tname && !Strcmp(tname, pname))
	    {
		remove_this = true;
		break;
	    }
	}
	Term* t = c->rterm;
	const char* tname = t->name();
	if (tname && !Strcmp(tname, pname))
	{
	    remove_this = true;
	}
	if (remove_this)
	{
	    remove (i+1);
	}
	else
	{
	    i++;
	}
    }
}
	    
const char* Constraint::remove_spec (const char *cstr)
{

// If string does not contain =, add "=1"
    const char* ptr = cstr;
    char* freestring = 0;
    if (!strstr (cstr, "="))
    {
	freestring = (char*) malloc (strlen(cstr)+4);
	strcpy (freestring, cstr);
	strcat (freestring, "=1");
	ptr = freestring;
    }
    const char* message = 0;
    Constraint* c = parse (ptr, &message);
    if (freestring) free (freestring);
    if (!c) 
    {
	if (message)
	{
	    return message;
	}
	else
	{
	    return \
                "Did not understand, please try full constraint specification";
	}
    }
    int oldindex = c->find_matching_left_side();
    if (oldindex > 0)
    {
	delete c;
	return remove (oldindex);
    }
    delete c;
    return "No matching constraint to delete";
}





Constraint::~Constraint ()
{
    free (string);
    if (rterm) delete rterm;
    for (int i = 0; i < Nlterms; i++)
    {
	if (lterms[i]) delete lterms[i];
    }
}

Constraint* Constraint::parse (const char *str, const char** message)
{
    Constraint *c = new Constraint (str);
    char *freestring = Strdup (str);
    char *fsptr = freestring;
    squeeze (fsptr);
    for (;;)
    {
	Term *t = Term::scan ((const char**) &fsptr, message);
	if (!t)
	{
	    if (fsptr[0] == '=') break;
	    free (freestring);
	    delete c;
	    return 0;
	}
	c->lterms[c->Nlterms++] = t;
    }
    fsptr++;  // Skip over equals sign
    Term *t = Term::scan ((const char**) &fsptr, message);
    if (!t)
    {
	free (freestring);
	delete c;
	return 0;
    }
    c->rterm = t;
    if (TCL_OK != c->transpose ())
    {
	*message = "Invalid constraint semantics";
	free (freestring);
	delete c;
	return 0;
    }
    free (freestring);
    return c;
}

// Recognize principal constraint e2 + ... = 1
// and return corresponding trait

char* Constraint::principal ()
{
    char* traitname = 0;

    if (rterm->name() || rterm->factor() != 1.0) return 0;

    for (int i = 0; i < Nlterms; i++)
    {
	
	Term* t = lterms[i];
	const char* tname = t->name();
	if (!Strcmp (tname, "e2"))
	{
	    traitname = Strdup ("");
	    break;
	}
	if (strlen (tname) > 4 && strlen (tname) < 1024)
	{
	    char buf[1024];
	    strncpy (buf, tname, 1024);
	    buf[3] = '\0';
	    if (!Strcmp (buf, "e2(") &&	')' == buf[strlen(tname)-1])
	    {
		buf[strlen(tname)-1] = '\0';
		traitname = Strdup (&tname[3]);
		break;
	    }
	}
    }
    return traitname;
}
		
const char* Constraint::add ()
{
// Check for duplicates, replace old with new

    int dupnum = find_matching_left_side();
    if (dupnum > 0)
    {
	delete Constraints[dupnum-1];
	Constraints[dupnum-1] = this;
	return 0;
    }

// If this is a principal variance constraint,
//   Remove previous principal variance constraint from model

    char* traitname;
    char* testname;
    if (0 != (traitname = principal()))
    {
	for (int i = 0; i < _count; i++)
	{
	    char *testname;  // remember to free()
	    if (0 != (testname = Constraints[i]->principal()))
	    {
		if (!Strcmp (traitname, testname))
		{
		    free (traitname);
		    free (testname);
		    delete Constraints[i];
		    Constraints[i] = this;
		    return 0;
		}
		free (testname);
	    }
	}
    }
    free (traitname);
	
// This is a new ordinary constraint

    if (_count+1 < MAX_CONSTRAINTS)
    {
	Constraints[_count++] = this;
	return 0;
    }
    return "Constraint limit exceeded";
}


int Constraint::find_matching_left_side ()
{
    int i, j, k;
    for (i = 0; i < _count; i++)
    {
	Constraint *c = Constraints[i];
	if (Nlterms == c->Nlterms || (c->transposed &&
		Nlterms == (c->Nlterms - 1)))
	{

// Terms may be in a different order, so we compare using a quasi-set
// of terms

	    bool lterms_match = true;
	    Term *tset[MAX_TERMS];
	    for (j = 0; j < Nlterms; j++)
	    {
		tset[j] = c->lterms[j];
	    }
	    for (j = 0; j < Nlterms; j++)
	    {
		Term *t = lterms[j];
		int matching_term = -1;
		for (k = j; k < Nlterms; k++)
		{
		    if (t->cmp(tset[k]))
		    {
			matching_term = k;
			break;
		    }
		}
		if (matching_term >= 0)
		{
		    tset[k] = tset[j];  // if k>j, keep j's term active
		}
		else
		{
		    lterms_match = false;
		    break;
		}
	    }
	    if (lterms_match)
	    {
		return i+1;
	    }
	}
    }
    return 0;
}


bool is_allowable_char_in_name (char c)
{
    if (isspace(c)) return false;
    if (c == '+' || c == '-' || c == '/' || c == '*' ||
	c == '\0' ||
	c == '<' || c == '>' || c == '=')
    {
//	fprintf (stderr, "Not allowing char: %c\n", c);
	return false;
    }
//  printf ("allowing %c\n", c);
    return true;
}

// Assumes pre-squeezed string
Term *Term::scan (const char **string, const char **error_message)  
{
    *error_message = 0;
    double number = 1.0;
    char name[128];
    name[0] = '\0';
    double sign = 1.0;
    bool found_number = false;
    bool found_name = false;

    const char *s = *string;
    
    if (*s == '+' || *s == '-')
    {
	if (*s == '-')
	{
	    sign = -1.0;
	}
	s++;
    }

    bool found_bracket = false;
    if (isdigit (*s) || *s == '.' || *s == '+' || *s == '-')
    {
	char nbuf[128];
	int index=0;
	while ((!isalpha(*s) || *s=='e' || *s=='E') && *s != '\0')
	{
	    nbuf[index++] = *s++;
	    if (nbuf[index-1] == '<') break;
	}
	if (nbuf[index-1] == '<')
	{
	    index--;
	    found_bracket = true;
	}
	if (nbuf[index-1] == '*') index--;
	nbuf[index] = '\0';
	char junk[128];
	int count = sscanf (nbuf, "%lg %s", &number, &junk);
	if (count != 1)
	{
	    *error_message = "Invalid factor or sign in constraint";
	    return 0;
	}
	found_number = true;
    }
    int i;
    if (found_bracket || *s == '<')
    {
	if (*s == '<') s++;
	for (i = 0; *s != '\0' && *s != '>'; i++)
	{
	    name[i] = *s++;
	}
	if (*s != '>')
	{
	    *error_message = "Unbalanced angle brackets <>";
	    return 0;
	}
	s++;
    } else {
	for (i = 0; is_allowable_char_in_name (*s); i++)
	{
	    name[i] = *s++;
	}
    }
    name[i] = '\0';
    char *nname = name;
    if (strlen (nname) == 0)
    {
	nname = 0;
    }
    else
    {
	if (!Parameter::find (name))
	{
	    fprintf (stderr, "No parameter %s\n", name);
	    *error_message = "Undefined parameter in constraint";
	    return 0;
	}
	found_name = true;
    }
    if (!found_number && !found_name) 
    {
	*error_message = "Invalid term in constraint";
	return 0;
    }
    number *= sign;
    *string = s;
    return new Term (number, nname);
}


int Term::bind (Tcl_Interp *interp)
{
    if (name())
    {
	if (!(_parameter = Parameter::find (name())))
	{
	    char buf[256];
	    sprintf (buf, "Missing parameter named %s", name() );
	    RESULT_BUF (buf);
	    return TCL_ERROR;
	}
    }
    return TCL_OK;
}

int Constraint::bind (Tcl_Interp *interp)
{
    for (int i = 0; i < _count; i++)
    {
	Constraint *c = Constraints[i];
	for (int j = 0; j < c->Nlterms; j++)
	{
	    if (TCL_OK != c->lterms[j]->bind (interp))
	    {
		return TCL_ERROR;
	    }
	}
    }
    return TCL_OK;
}

int Constraint::transpose ()
{
    if (rterm->name())
    {
	if (Nlterms+1 > MAX_TERMS)
	{
	    return TCL_ERROR;
	}
	rterm->negate ();
	lterms[Nlterms++] = rterm;
	rterm = new Term (0.0, 0);
	transposed = true;
    }
    double accumulator = 0.0;
    int i = 0;
    while (i < Nlterms)
    {
	if (!lterms[i]->name())
	{
	    accumulator += lterms[i]->factor();
	    delete lterms[i];
	    for (int j = i+1; j < Nlterms; j++)
	    {
		lterms[j-1] = lterms[j];
	    }
	    Nlterms--;
	}
	else
	{
	    i++;
	}
    }
    rterm->add_to (accumulator);
    return TCL_OK;
}
	
extern "C" void getcons_ (double *cnstr, double *cvalue) 
{
    Constraint::apply (cnstr, cvalue);
}

// Return count of real and virtual covariates

int Constraint::count()
{
    int nvcons = 0;
    Covariate *c;
    for (int i = 0; c = Covariate::index(i); i++)
    {
	if (!c->active())
	{
	    for (int t = 0; t < Trait::Number_Of(); t++)
	    {
		if (c->beta(t))
		{
		    nvcons++;
		}
	    }
	}
    }		
    return _count + nvcons;
}


void Constraint::apply (double *cnstr, double *cvalue)
{
    int ncons = Constraint::_count;
    int nvcons = Constraint::count();
    int i;
    for (i = 0; i < ncons; i++)
    {
	Constraint *c = Constraints[i];
	for (int j = 0; j < c->Nlterms; j++)
	{
	    Parameter *p = c->lterms[j]->parameter();
	    cnstr[nvcons * p->number() + i] = c->lterms[j]->factor();
	}
	cvalue[i] = c->rterm->factor();
    }

// Generate virtual constraints for suspended covariates

    int k = ncons;
    Covariate *c;
    for (i = 0; c = Covariate::index(i); i++)
    {
	if (!c->active())
	{
	    for (int t = 0; t < Trait::Number_Of(); t++)
	    {
		Parameter *p = c->beta(t);
		if (p)
		{
		    cnstr[nvcons * p->number() + k] = 1.0;
		    cvalue[k] = 0.0;
		    k++;
		}
	    }
	}
    }
    if (Maximize_InitPar)
    {
	throw Safe_Error_Return ("InitPar only");
    }
}

void Constraint::write_commands (FILE *file)
{
    for (int i = 0; i < _count; i++)
    {
	char buf[1024];
	fprintf (file, "constraint %s\n", Constraints[i]->string);
    }
}

// Check if a parameter is constrained to a constant value
// If so, make sure it is set to that value, and that bounds are appropriate

extern "C" void concheck_ (int *PNUM, double *START, double *LOWER,
			   double *UPPER)
{
    int pnum = *PNUM - 1;
    Constraint::check_parameter (pnum, START, LOWER, UPPER);
}

const double Fudge = 0.001;

//
// Check, by name, if a parameter is part of
// constant constraint (e.g. h2r = 0)
//
bool Constraint::Check_by_name (const char* name)
{
    bool constrained = false;
    for (int i = 0; i < _count; i++)
    {
	Constraint *c = Constraint::Constraints[i];
	if (c->Nlterms > 1) continue;
	if (!Strcmp (name,c->lterms[0]->name()))
	{
	    constrained = true;
	    break;
	}
    }
    return constrained;
}

int Constraint::check_parameter (int pnum, double *START, double *LOWER,
				  double *UPPER)
{
    Parameter *par = Parameter::index (pnum);
    if (!par) return 0;
    
    for (int i = 0; i < _count; i++)
    {
	Constraint *c = Constraint::Constraints[i];
	if (c->Nlterms > 1) continue;
	if (c->lterms[0]->parameter() == par)
	{
	    double fac = c->lterms[0]->factor();
	    double rvalue = c->rterm->factor();
	    if (fac == 0.) return 0;           // Let FORTRAN deal with it
	    double value = rvalue / fac;
	    *START = value;
	    if (*LOWER >= *START) *LOWER = *START - Fudge;
	    if (*UPPER <= *START) *UPPER = *START + Fudge;
	    return 1;
	}
    }
    return 0;
}
	    
extern "C" void conotsat_ (int* I)
{
    Constraint::Not_Satisfied ((*I) - 1);
}

void Constraint::Not_Satisfied (int i)
{
    if (i < _count)
    {
	char buf[10000];
	sprintf (buf, "Constraint not satisfied: %s", Constraints[i]->string);
	throw Safe_Error_Return (buf);
    } else {
	throw Safe_Error_Return ("Implicit model constraint not satisfied");
    }
}

