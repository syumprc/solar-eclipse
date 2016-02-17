/*
 * field.cc implements the field class and command
 * Written by Charles Peterson beginning on June 3, 1998
 * Copyright (c) 1998 Southwest Foundation for Biomedical Research
 *
 * Field maps a user field name to an internal (default) name
 */

#include "solar.h"

Field *Field::Fields = Fields_Setup();

Field *Field::Fields_Setup ()
{
    Field *start = 
	new Field ("ID","<Ego ID> (Individual Permanent ID) [mandatory]",
		   mandatory);
    start->add ("FA","<Father ID> (Father's Permanent ID) [mandatory]", 
		mandatory);
    start->add ("MO","<Mother ID> (Mother's Permanent ID) [mandatory]", 
		mandatory);
    start->add ("SEX", "<Sex> (M/F, m/f, or 1/2) [mandatory]", 
		mandatory);
    start->add ("PROBND","<Proband Status> [optional]", 
		optional);
    start->add ("MZTWIN","<Monozygotic Twin Group> [optional]", 
		optional);
    start->add ("FAMID","<Family ID> [optional]",
		optional);
    start->add ("HHID","<Household ID> [optional]", 
		optional);
    return start;
}

void Field::Start()
{
    FILE *finfo = fopen ("field.info", "r");
    if (finfo)
    {
	char fline[1024];
	bool errors = false;
	while (fgets (fline, 1024, finfo))
	{
	    if (strlen (fline) && fline[0]!='#')
	    {
		char field[1024];
		char map[1024];
		char junk[1024];
		if (2 != sscanf (fline, "%s %s %s", field, map, junk))
		{
		    fprintf  (stderr, "\nInvalid line in field.info:  %s",
			      fline);
		    errors = true;
		}
		else
		{
		    const char* message = Field::Change (field, map);
		    if (message)
		    {
			fprintf (stderr, "\nError in field.info: %s:  %s",
				message, fline);
			errors = true;
		    }
		}
	    }
	}
	if (errors) fprintf (stderr, "\n");
	Fclose (finfo);
    }
}

			
Field *Field::add (const char *name, const char *description, mandate man)
{
    Field *f = this;
    while (f->_next) f = f->_next;
    f->_next = new Field (name, description, man);
    return this;
}

Field *Field::Find (const char *int_name)
{
    Field *f = Fields;
    do
    {
	if (!Strcmp (int_name, f->_internal_name))
	{
	    return f;
	}
    }
    while (f=f->_next);
    return 0;
}

const char *Field::Map (const char *int_name)
{
    Field *f = Find (int_name);
    if (!f)
    {
	char buf[256];
	sprintf (buf, "Internal name %s has not been defined", int_name);
	throw Safe_Error_Return (buf);
    }
    return f->_user_name;
}

bool Field::Changed (const char *int_name)
{
    Field *f = Find (int_name);
    if (!f)
    {
	char buf[256];
	sprintf (buf, "Internal name %s has not been defined", int_name);
	throw Safe_Error_Return (buf);
    }
    return f->_changed;
}



const char *Field::Change (const char *int_name, const char *user_name)
{
// Find int_name
    Field *f = Find (int_name);
    if (!f) return "No such internal field name is defined";
    if (!Strcmp (user_name, "-none"))
    {
	if (f->_mandate == mandatory)
	{
	    return "Field is mandatory";
	}
    }
    if (!Strcmp (user_name, "-default"))
    {
	f->_user_name = Strdup (f->_internal_name);
	f->_changed = false;
    } else {
	f->_user_name = Strdup (user_name);
	f->_changed = true;
    }
    Phenotypes::start();  // May have changed interpretation
    return 0;
}
    


extern "C" int FieldCmd (ClientData clientData, Tcl_Interp *interp,
	      int argc, char *argv[])
{
    if (argc == 2 && !Strcmp ("help", argv[1]))
    {
	return Solar_Eval (interp, "help field");
    }

    if (argc == 1)
    {
	Field *f = Field::All();
	const char *name;
	const char *user_name;
	const char *description;
#define FIELD_FORMAT "%-9s %-15s %s\n"

	printf (FIELD_FORMAT, "Default", "Current", "Meaning");
	printf (FIELD_FORMAT, "-------", "-------", "-------");
	while (Field::Info (&f, &name, &user_name, &description))
	{
	    printf (FIELD_FORMAT, name, user_name, description);
	}
	return TCL_OK;
    }

    if (argc == 2)
    {
	const char *user_name = 0;
	try
	{
	    RESULT_BUF (Field::Map (argv[1]));
	    return TCL_OK;
	}
	catch (Safe_Error_Return)
	{
	    RESULT_LIT ("No such internal name is defined for mapping");
	}
	return TCL_ERROR;
    }

    if (argc == 3)
    {
	const char *message = Field::Change (argv[1], argv[2]);
	if (message)
	{
	    RESULT_BUF (message);
	    return TCL_ERROR;
	}
	if (strlen(argv[1]) < 256 && strlen(argv[2]) < 256)
	{
	    char command[1024];
	    sprintf (command, "field_info_update %s %s", argv[1], argv[2]);
	    return Solar_Eval (interp, command);
	}
    }

    RESULT_LIT ("Invalid field command");
    return TCL_OK;
}

	
