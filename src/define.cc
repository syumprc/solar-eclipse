/*
 * definition.cc implements the "define" command
 * Written by Charles Peterson beginning on Feb 3, 2005
 * Copyright (c) 2005 Southwest Foundation for Biomedical Research
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "solar.h"
// tcl.h from solar.h
#include "safelib.h"

Definition** Definition::Definitions = 0;
int Definition::End = -1;

Definition* Definition::Find (const char* name)
{
    int index = find_index (name);
    if (index < 0) return 0;
    return Definitions[index];
}

extern "C" int DefinitionCmd (ClientData clientData, Tcl_Interp *interp,
		  int argc, char *argv[])
{
    if (argc == 1)
    {
	return Definition::Show_All (interp);
    }
    if (argc == 2 && !Strcmp ("help", argv[1]))
    {
	return Solar_Eval (interp, "help trait");
    }

    if (argc == 2 && !Strcmp ("Delete_All", argv[1]))
    {
	Definition::Delete_All();
	return TCL_OK;
    }

    if (argc == 2 && !Strcmp ("new", argv[1]))
    {
	Definition::Delete_All();
	return TCL_OK;
    }

    if (argc == 2 && !Strcmp ("names", argv[1]))
    {
	return Definition::Show_Names(interp);
    }

    if (argc == 3 && !Strcmp ("delete", argv[1]))
    {
	return Definition::Delete (argv[2], interp);
    }

    if (argc == 3 && !Strcmp ("delete", argv[2]))
    {
	return Definition::Delete (argv[1], interp);
    }

    if (argc == 4 && !Strcmp ("rename", argv[1]))
    {
	return Definition::Rename (argv[2], argv[3], interp);
    }

    
    if (argc == 2 || (argc == 3 && !Strcmp (argv[2], "=")))
    {
	return Definition::Show (argv[1], interp);
    }

    if (argc == 4 && !Strcmp ("=", argv[2]))
    {
	return Definition::Assign (argv[1], argv[3], interp);
    }

    if (argc > 4 && !Strcmp ("=", argv[2]))
    {
	char* def_string = Tcl_Concat (argc-3, &argv[3]);
	int status = Definition::Assign (argv[1], def_string, interp);
	Tcl_Free (def_string);
	return status;
    }

    RESULT_LIT ("Invalid Define Command");
    return TCL_ERROR;
}

int Definition::find_index (const char* name)
{
    if (Definitions)
    {
	for (int i = 0; i <= End; i++)
	{
	    if (!Strcmp (name, Definitions[i]->_name))
	    {
		return i;
	    }
	}
    }
    return -1;
}

Definition::~Definition ()
{
    if (_name) free (_name);
   if (_definition) free (_definition);
}

int Definition::Rename (const char* name1, const char* name2,
			   Tcl_Interp* interp)
{
// See if name1 exists (Don't save index now, it might change.)
    if (find_index (name1) < 0)
    {
	char* errmes = Strdup ("Definition not found: ");
	StringAppend (&errmes, name1);
	RESULT_BUF (errmes);
	free (errmes);
	return TCL_ERROR;
    }

// See if first and second names are exactly the same
    if (!strcmp (name1, name2))
    {
	RESULT_LIT ("define rename: old and new names are identical");
	return TCL_ERROR;
    }

// See if first and second names are the same except for case
    bool change_case = false;
    if (!Strcmp (name1, name2))
    {
	change_case = true;
    }

// See if second name exists already.  If it does, it's deleted,
// (unless case changes)

    if (!change_case && find_index(name2) >= 0)
    {
	Delete (name2, interp);
    }

// OK, now do the rename
    int index = find_index(name1);
    if (index < 0)
    {
// This can't happen.  We already checked for it...
	RESULT_BUF ("Definition rename failure 9975.  Contact SOLAR developers.");
	return TCL_ERROR;
    }

    free (Definitions[index]->_name);
    Definitions[index]->_name = Strdup (name2);
    return TCL_OK;
}

int Definition::Show (const char* name, Tcl_Interp* interp)
{
    Definition* exp = Find (name);
    if (!exp)
    {
	char* errmes = Strdup ("Definition not found: ");
	StringAppend (&errmes, name);
	RESULT_BUF (errmes);
	free (errmes);
	return TCL_ERROR;
    }
    char* result = Strdup (name);
    StringAppend (&result, " = ");
    StringAppend (&result, exp->_definition);
    RESULT_BUF (result);
    free (result);
    return TCL_OK;
}

int Definition::Show_Names (Tcl_Interp* interp)
{
    if (Definitions)
    {
	char* name_list = Strdup ("");

	for (int i = 0; i <= End; i++)
	{
	    if (i > 0) StringAppend (&name_list, " ");
	    StringAppend (&name_list, Definitions[i]->_name);
	}
	RESULT_BUF (name_list);
	free (name_list);
    }
    return TCL_OK;
}

void Definition::Write_Commands (FILE *file)
{
    if (Definitions)
    {
	for (int i = 0; i <= End; i++)
	{
	    fprintf (file, "define %s = %s\n", Definitions[i]->_name,
		     Definitions[i]->_definition);
	}
    }
}

int Definition::Show_All (Tcl_Interp* interp)
{
    if (Definitions)
    {
	char* result = Strdup ("");
	for (int i = 0; i <= End; i++)
	{
	    StringAppend (&result, Definitions[i]->_name);
	    StringAppend (&result, " = ");
	    StringAppend (&result, Definitions[i]->_definition);
	    if (i != End)
	    {
		StringAppend (&result, " \n");
	    }
	}
	RESULT_BUF (result);
	free (result);
    }
    return TCL_OK;
}
	    

int Definition::Delete (const char* name, Tcl_Interp* interp)
{
    int index = find_index (name);
    if (index < 0)
    {
	char* errmes = Strdup ("Definition not found: ");
	StringAppend (&errmes, name);
	RESULT_BUF (errmes);
	free (errmes);
	return TCL_ERROR;
    }
    delete Definitions[index];
    int i = index;
    for (; i < End; i++)
    {
	Definitions[i] = Definitions[i+1];
    }
    Definitions[i] = 0;
    --End;
    if (End < 0)
    {
	free (Definitions);
	Definitions = 0;
    }
    else
    {
	Definitions = (Definition **) Realloc (Definitions, (End+1)*
						 sizeof (Definition*));
    }
    return TCL_OK;
}
	

void Definition::Delete_All ()
{
    if (Definitions)
    {
	for (int i = 0; i <= End; i++)
	{
	    delete Definitions[i];
	}
    }
    free (Definitions);
    Definitions = 0;
    End = -1;
}

int Definition::Assign (const char* name, const char* estring, Tcl_Interp* interp)
{
    try
    {
	Expression* val = new Expression (estring);
	delete val;
    }
    catch (...)
    {
	char* errmes = Strdup ("Invalid expression: ");
	StringAppend (&errmes, estring);
	RESULT_BUF (errmes);
	free (errmes);
	return TCL_ERROR;
    }

// Check if name already in use, if so, it will be reused
    Definition* newexp;
    Definition* oldexp = Definition::Find (name);
    if (oldexp)
    {
	newexp = oldexp;
	free (newexp->_definition);
    }
    else
    {

// Add this to the collection
	newexp = new Definition;
	newexp->_name = Strdup (name);
	End++;
	if (!Definitions)
	{
	    Definitions = (Definition**) malloc (sizeof (Definition*));
	}
	else
	{
	    Definitions = (Definition**) realloc (Definitions, (End+1)*sizeof(Definition*));
	}
	Definitions[End] = newexp;
    }
    newexp->_definition = Strdup (estring);
    return TCL_OK;
}



