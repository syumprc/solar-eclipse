/*
 * trait.cc implements the trait command
 * Written by Charles Peterson beginning on October 27, 1997
 * Copyright (c) 1997 Southwest Foundation for Biomedical Research
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "solar.h"
// tcl.h from solar.h
#include "safelib.h"

#define INVALID_TRAIT_NAME "Must_give_command_model_new"

Trait*     Trait::Traits[MAX_TRAITS] = {0};
int        Trait::Number_Of_Traits = 0;
int        Trait::_Number_Of_Expressions = 0;
bool       Trait::Valid = true;
Context*   Trait::_Context = 0;

Trait::Trait (char* trait_name, char* stdev_name)
{
    _name = Strdup (trait_name);
    _stdev = stdev_name;
    if (stdev_name)
    {
	_stdev = Strdup (stdev_name);
    }
    _type = unknown;
    _expression = 0;
}

void Trait::Set_Quantitative (int i)
{
    Traits[i]->_type = quantitative;
}

void Trait::Set_Discrete (int i)
{
    Traits[i]->_type = discrete;
}


Trait::~Trait ()
{
    free (_name);
    if (_stdev) free (_stdev);
    if (_expression) delete _expression;
}

int Trait::One (char* trait_name, char* stdev_name)
{
    Reset ();  // Ensure no other traits
    Trait *t = new Trait (trait_name, stdev_name);
    Traits[0] = t;
    Number_Of_Traits = 1;
    return 0;
}

int Trait::Add (char* trait_name, char* stdev_name)
{
    if (Number_Of_Traits >= MAX_TRAITS)
    {
	return 1;
    }	
    Trait *t = new Trait (trait_name, stdev_name);
    Traits[Number_Of_Traits++] = t;
    return 0;
}

// Update parameters is called when trait changes
// and when -noparm is not in effect

void Trait::Update_Parameters (Tcl_Interp* interp, char* oldtrait1, char* oldtrait2)
{
    Parameter *p;

// One trait to one trait is fairly simple
//   zero out mean and sd
//   reset boundaries (not starting points) of variance components

    if (1==Trait::Number_Of())
    {
	p = Parameter::find("mean");
	if (p) p->zero();
	p = Parameter::find("sd");
	if (p) p->zero();
	
// Reset variance component boundaries to default values

	p = Parameter::find("e2");
	if (p) 
	{
	    p->lower = (p->start <= 0.01) ? 0.0 : 0.01;
	    p->upper = 1;
	}
	p = Parameter::find("h2r");
	if (p) {p->lower = 0; p->upper = 1;}
        int hindex = 1;
	for (;;hindex++)
	{
	    char hname[128];
	    sprintf (hname, "h2q%i", hindex);
	    p = Parameter::find (hname);
	    if (p) {
		p->lower = 0; p->upper = 1;
	    } else {
		break;
	    }
	}

// Two traits to two traits is more complicated
// Old parameters need to be renamed as well as zeroed or rebounded
// renaming has to propagate to omega, constraints, and mu

// This is now obsolete.  Model new must be given before changing traits

    } else if (0) {
	
	for (int trno = 0; trno < 2; trno++)
	{
	    char* tname = (trno==0) ? oldtrait1 : oldtrait2;
	    char pname[1024];
	    sprintf (pname, "mean(%s)", tname);
	    p = Parameter::find (pname);
	    if (p)
	    {
		sprintf (pname, "mean(%s)", Trait::Name(trno));
		p->rename (pname);
		p->zero();
	    }
	    sprintf (pname, "sd(%s)", tname);
	    p = Parameter::find (pname);
	    if (p)
	    {
		sprintf (pname, "sd(%s)", Trait::Name(trno));
		p->rename (pname);
		p->zero();
	    }
	    sprintf (pname, "e2(%s)", tname);
	    p = Parameter::find (pname);
	    if (p)
	    {
		sprintf (pname, "e2(%s)", Trait::Name(trno));
		p->rename (pname);
		p->lower = (p->start <= 0.01) ? 0.0 : 0.01;
		p->upper = 1.0;
	    }
	    sprintf (pname, "h2r(%s)", tname);
	    p = Parameter::find (pname);
	    if (p)
	    {
		sprintf (pname, "sd(%s)", Trait::Name(trno));
		p->rename (pname);
		p->zero();
	    }
	    for (int hindex = 1; 1; hindex++)
	    {
		sprintf (pname, "h2q%i(%s)", hindex, tname);
		p = Parameter::find (pname);
		if (p) {
		    sprintf (pname, "h2q%i(%s)", hindex, Trait::Name(trno));
		    p->rename (pname);
		    p->lower = 0; 
		    p->upper = 1.0;
		} else {
		    break;
		}
	    }
	    char command[1024];
	    sprintf (command, "omegasubst %s %s", tname, Trait::Name(trno));
	    int status = Solar_Eval (interp, command);
	    sprintf (command, "constsubst %s %s", tname, Trait::Name(trno));
	    status = Solar_Eval (interp, command);
	    sprintf (command, "musubst %s %s", tname, Trait::Name(trno));
	    status = Solar_Eval (interp, command);
	}
	p = Parameter::find ("rhoe");
	p->lower = (p->start <= -0.9) ? -1.0 : -0.9;
	p->upper = (p->start >= 0.9)  ?  1.0 : 0.9;
	p = Parameter::find ("rhog");
	p->lower = (p->start <= -0.9) ? -1.0 : -0.9;
	p->upper = (p->start >= 0.9)  ?  1.0 : 0.9;
	for (int hindex = 1; 1; hindex++)
	{
	    char hname[128];
	    sprintf (hname, "rhoq%i", hindex);
	    p = Parameter::find (hname);
	    if (p)
	    {
		p->lower = (p->start <= -0.9) ? -1.0 : -0.9;
		p->upper = (p->start >= 0.9)  ?  1.0 : 0.9;
	    }
	}
    }
}

void Trait::Reset()
{
    for (int i = 0; i < Trait::Number_Of(); i++)
    {
	delete Traits[i];
	Traits[i] = 0;
    }
    Number_Of_Traits = 0;
}

int Trait::Find_Parameter (int index, const char* name)
{
    if (index >= Trait::Number_Of())
    {
	error ("Trait index invalid; contact solar@txbiomedgenetics.org");
    }

    Parameter *p;
    char buf[1024];
    if (Trait::Number_Of() < 2)
    {
	strcpy (buf,name);
    }
    else
    {
	sprintf (buf, "%s(%s)", name, Trait::Name(index));
    }
    p = Parameter::find (buf);

    if (!p) return -1;

    return p->number();
}


const char* Trait::Name (int i)
{
    if (i > Number_Of() || i < 0)
    {
	fprintf (stderr, 
		 "Trait index out of bounds; contact solar@txbiomedgenetics.org");
	exit (0);
    }
    return Traits[i]->_name;
}

Expression* Trait::expression (int i)
{
    if (i > Number_Of() || i < 0)
    {
	fprintf (stderr, 
		 "Trait index out of bounds; contact solar@txbiomedgenetics.org");
	exit (0);
    }
    return Traits[i]->_expression;
}

const char* Trait::Name_Stdev (int i)
{
    if (i > Number_Of())
    {
	fprintf (stderr, 
		 "Trait index out of bounds; contact solar@txbiomedgenetics.org");
	exit (0);
    }
    return Traits[i]->_stdev;
}


trait_type Trait::Type (int i)
{
    if (i < Number_Of_Traits && i >= 0)
    {
	return Traits[i]->_type;
    }
    return unknown;
}


//  If any trait is unknown, the mt is unknown, then,
//    if any trait is discrete, the mt is discrete, otherwise
//      the mt is quantitative
trait_type Trait::Maximization_Type ()
{
    trait_type mt = quantitative;
    for (int i = 0; i < Number_Of_Traits; i++)
    {
	if (unknown == Traits[i]->_type)
	{
	    mt = unknown;
	    break;
	}
	else if (discrete == Traits[i]->_type)
	{
	    mt = discrete;
	}
    }
    return mt;
}
    

// Check for inor m a l _ and return phenotype name

char* inor_match_class (char* string, bool* ifclass, int* pclass)
{
    char* rval = 0;
    *pclass = 0;

    if (!strncmp (string, "inor_", 5))
	rval = string+5;
    else if (!strncmp (string, "inorm_", 6))
	rval = string+6;
    else if (!strncmp (string, "inorma_", 7))
	rval = string+7;
    else if (!strncmp (string, "inormal_", 8))
	rval = string+8;
    else if (*ifclass && !strncmp (string, "inormalc_", 9))
    {
	char* classp = string+9;
	char* endptr;
	errno = 0;
	long longc = strtol (classp, &endptr, 0);
	if (errno || longc > INT_MAX || longc < INT_MIN)
	{
	    printf ("Invalid class specified in %s\n",string);
	    rval = string + strlen(string);
	    return rval;
	}
	*pclass = (int) longc;
	*ifclass = true;
	rval = endptr+1;
	return rval;
    }
    *ifclass = false;
    return rval;
}

char* inor_match (char* string)
{
    bool ifclass = false;
    int ignore;
    return inor_match_class (string, &ifclass, &ignore);
}

char* zscor_match (char* string)
{
    char* rval = 0;

    if (!strncmp (string, "zscor_", 6))
	rval = string+6;
    else if (!strncmp (string, "zscore_", 7))
	rval = string+7;

    return rval;
}



int Trait::Bind (Tcl_Interp *interp)
{
    if (!Valid)
    {
	RESULT_LIT ("Invalid trait specification");
	return TCL_ERROR;
    }

    if (Number_Of_Traits == 0)
    {
	RESULT_LIT ("No trait has been specified");
	return TCL_ERROR;
    }
    if (_Context) delete _Context;
    _Context = 0;
    _Number_Of_Expressions = 0;

// cleanup inverse normals from last time

    if (Phenotypes::INNames)
    {
	int iin;
	for (iin = 0; iin < Phenotypes::INCount; iin++)
	{
	    if (Phenotypes::INNames[iin]) free (Phenotypes::INNames[iin]);
	}
	free (Phenotypes::INNames);
	Phenotypes::INNames = 0;
    }
    Phenotypes::INCount = 0;
    if (Phenotypes::INData)
    {
	free (Phenotypes::INData);
	Phenotypes::INData = 0;
    }
    if (Phenotypes::INIfclass)
    {
	free (Phenotypes::INIfclass);
	Phenotypes::INIfclass = 0;
    }
    if (Phenotypes::INClass)
    {
	free (Phenotypes::INClass);
	Phenotypes::INClass = 0;
    }
	    
    if (Solar_Eval (interp, "inormal -reset"))
    {
	RESULT_LIT ("Attempt to reset inormal data failed");
	return TCL_ERROR;
    }

// Setup each trait

    for (int i = 0; i < Number_Of_Traits; i++)
    {
	Trait *t = Traits[i];
	bool use_phen = Phenotypes::available (t->_name);
	Definition* use_def = Definition::Find(t->_name);
	if (!use_phen && !use_def)
	{
	    char buf[256];
	    sprintf (buf, "Trait variable '%s' not found",  t->_name);
	    RESULT_BUF (buf);
	    return TCL_ERROR;
	}
	if (use_phen && use_def)
	{
	    char buf[256];
	    sprintf (buf, 
		     "Trait variable '%s' both phenotype and definition",  
		     t->_name);
	    RESULT_BUF (buf);
	    return TCL_ERROR;
	}
	if (t->_expression)
	{
	    delete t->_expression;
	    t->_expression = 0;
	}
	if (use_def)
	{
	    _Number_Of_Expressions++;
	    try
	    {
		t->_expression = new Expression (use_def->definition());
	    }
	    catch (...)
	    {
		char* errmes = Strdup ("Expression error in define command for trait ");
		StringAppend (&errmes, t->_name);
		RESULT_BUF (errmes);
		free (errmes);
		return TCL_ERROR;
	    }
	    Expression* e = t->_expression;
	    if (!_Context)
	    {
		_Context = new Context;
	    }

// Establish context, adding required variables and normals until everything defined

	    bool ok = false;
	    int pclass;

	    while (!ok)
	    {
		bool ifclass = true;
		try
		{
		    e->bind (_Context);
		    ok = true;
		}
		catch (Undefined_Name& un)
		{
		    char *pname;
		    try
		    {
			if (Phenotypes::available (un.name))
			{
			    int vindex = EqVar::add_for_trait (un.name, i);
			    _Context->add (un.name, eval_phenotype, vindex);
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
			    if (ifclass)
			    {
				sprintf (cbuf, "inormal -trait %s -phenfile -class %d",
					 pname, pclass);
			    } else {
				sprintf (cbuf, "inormal -trait %s -phenfile", pname);
			    }				
			    int errstatus = Solar_Eval (interp, cbuf);
			    if (errstatus)
			    {
				printf ("maximize: error applying inormal to trait %s\n", pname);
				return TCL_ERROR;   // Error message from inormal in TCL_Result
			    }
			    
			    if (!Phenotypes::INCount)
			    {
				Phenotypes::INData = (double*) Calloc (2, sizeof (double));
				Phenotypes::INNames = (char**) Calloc (2, sizeof (char*));
				Phenotypes::INIfclass = (bool*) Calloc (2, sizeof (int));
				Phenotypes::INClass = (int*) Calloc (2, sizeof (bool));
			    }
			    else
			    {
				int newsize = Phenotypes::INCount + 2;
				Phenotypes::INData = (double*) Realloc (Phenotypes::INData, newsize*sizeof (double));
				Phenotypes::INNames = (char**) Realloc (Phenotypes::INNames, newsize*sizeof (char*));
				Phenotypes::INIfclass = (bool*) Realloc (Phenotypes::INIfclass, newsize*sizeof(bool));
				Phenotypes::INClass = (int*) Realloc (
				    Phenotypes::INClass, newsize*sizeof(int));
			    }
			    Phenotypes::INNames[Phenotypes::INCount] = Strdup (pname);
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
				 "Undefined variable `%s' in trait expression", un.name);
			RESULT_BUF (buf);
			return TCL_ERROR;
		    }
		}
	    }
	}
	if (t->_stdev)
	{
	    char* sname = t->_stdev;
	    if (!Phenotypes::available (sname))
	    {
		char buf[256];
		sprintf (buf, "Trait stdev variable '%s' not found", sname);
		RESULT_BUF (buf);
		return TCL_ERROR;
	    }
	}
    }
    return TCL_OK;
}

    


extern "C" int TraitCmd (ClientData clientData, Tcl_Interp *interp,
		  int argc, char *argv[])
{
    bool noparm = false;

    if (argc == 2 && !StringCmp ("help", argv[1], case_ins))
    {
	return Solar_Eval (interp, "help trait");
    }

    if (argc == 2 && !Strcmp ("-is-trait-discrete?", argv[1]))
    {
	if (Trait::Number_Of() > 1)
	{
	    RESULT_LIT ("0");
	    return TCL_OK;
	}
	if (Trait::Type (0) == discrete)
	{
	    RESULT_LIT ("1");
	    return TCL_OK;
	} 
	else if (Trait::Type(0) == quantitative)
	{
	    RESULT_LIT ("0");
	    return TCL_OK;
	}
	RESULT_LIT ("Model has not been maximized");
	return TCL_ERROR;
    }

    if (argc >= 3 && !Strcmp ("-noparm", argv[1]))
    {
	noparm = true;
	int i;
	for (i = 1; i < argc-1; i++)
	{
	    argv[i] = argv[i+1];
	}
	argc--;
    }
	
// Describe trait(s)
    if (argc == 1)
    {
	if (!Trait::Number_Of())
	{
	    RESULT_LIT ("No trait has been specified");
	    return TCL_ERROR;
	}
	char* trait_description = Trait::Describe_All();
	RESULT_BUF (trait_description);
	free (trait_description);
	return TCL_OK;
    }

// Setup new trait(s)

// Check if model new would be needed

    int oldtraitcount = Trait::Number_Of();
    if (Parameter::count() > 0 &&
	(oldtraitcount!=0 && (oldtraitcount >= 2 || argc > 2)))
    {
	Trait::Valid = false;
	RESULT_BUF (
	    "Must give command 'model new' before changing traits if more than one trait");
	return TCL_ERROR;
    }

// Start by saving old traits for later if required

    char oldtrait1[1024] = {0};
    char oldtrait2[1024] = {0};
    if (oldtraitcount > 0)	strcpy (oldtrait1, Trait::Name(0));
    if (oldtraitcount > 1)	strcpy (oldtrait1, Trait::Name(1));
    
// Reset old trait(s), especially so failure results in null trait

    Trait::Valid = false;
    int argindex = 1;
    bool first_trait_set = false;
    bool invalid_trait_name = false;
    while (argindex < argc)
    {
	char* trait_name = argv[argindex++];
	char *stdev_name = 0;

	if (!Strcmp (trait_name, INVALID_TRAIT_NAME))
	{
	    invalid_trait_name = true;
	}
	else if (strlen(trait_name) > 80) {
	    RESULT_LIT ("trait name exceeds 80 character limit");
	    return TCL_ERROR;
	}
	else if (!Phenotypes::available (trait_name) && !Definition::Find(trait_name))
	{
	    char buf[256];
	    sprintf (buf,
        "Selected trait '%s' not available; have you loaded phenotypes file?", 
		     trait_name);
	    RESULT_BUF (buf);
	    return TCL_ERROR;
	} 
	if (argindex+2 < argc)
	{
	    if (!Strcmp (argv[argindex+1], "-stdev"))
	    {
		stdev_name = argv[argindex+2];
		if (!Phenotypes::available (stdev_name))
		{
		    char buf[256];
		    sprintf (buf,
			 "Selected stdev variable %s not available",
			     stdev_name);
		    RESULT_BUF (buf);
		    return TCL_ERROR;
		}
		argindex += 2;
	    }
	}
	if (!first_trait_set)
	{
	    Trait::One (trait_name, stdev_name);
	    first_trait_set = true;
	}
	else
	{
	    if (Trait::Add (trait_name, stdev_name))
	    {
		RESULT_LIT ("Only two traits are permitted");
		return TCL_ERROR;
	    }
	}
    }

/*
 * do stuff we need to do for new trait
 */
    if (oldtraitcount!=0 && !noparm)
    {
	Trait::Update_Parameters (interp, oldtrait1, oldtrait2);
    }
 
    Covariate::new_traits ();

    if (!invalid_trait_name)
    {
	Trait::Valid = true;
    }

    int status = TCL_OK;
    if (Option::get_int("zscore")) {
	status = Solar_Eval (interp, "zscore -off");
	RESULT_LIT ("Zscore turned off because new trait");
    }
    return status;
}

void Trait::Write_Commands (FILE *file)
{
    if (Number_Of_Traits)
    {
	char* trait_description = Trait::Describe_All();
	fprintf (file, "trait %s\n", trait_description);
	free (trait_description);
    }
}

char *Trait::Describe_All ()
{
    char* trait_description = 0;
    if (!Valid)
    {
	StringAppend (&trait_description, INVALID_TRAIT_NAME);
	return trait_description;
    }
    for (int i = 0; i < Number_Of_Traits; i++)
    {
	if (i > 0) StringAppend (&trait_description, " ");
	StringAppend (&trait_description, Traits[i]->_name);
	if (Traits[i]->_stdev)
	{
	    StringAppend (&trait_description, " -stdev ");
	    StringAppend (&trait_description, Traits[i]->_stdev);
	}
    }
    return trait_description;
}

// Fortran interface to get "actual" number of traits


extern "C" int ttraits_ ()
{
    return Trait::Number_Of();
}


