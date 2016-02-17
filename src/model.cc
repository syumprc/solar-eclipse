/* 
 * model.cc
 * Written by Charles P. Peterson
 * Starting on December 2, 1997
 * Copyright (c) 1997 Southwest Foundation for Biomedical Research
 *
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "solar.h"

bool Model::_loading = {false};

extern "C" int ModelCmd (ClientData clientData, Tcl_Interp *interp,
		  int argc, char *argv[])
{
    if (argc == 2 && !StringCmp (argv[1], "help", case_ins))
    {
	return Solar_Eval (interp, "help model");
    }
    
    if (argc == 1)
    {
	Model::write (stdout);
	return TCL_OK;
    }

    if (argc == 2 && !StringCmp (argv[1], "new", case_ins))
    {
	return Model::renew (interp);
    }

    if (argc == 3)
    {
	if (!StringCmp (argv[1], "save", case_ins))
	{
	    char *given_filename = argv[2];
	    char *use_filename = append_extension (given_filename, ".mod");
	    FILE *file = fopen (use_filename, "w");
	    if (!file)
	    {
		char buf[1024];
		sprintf (buf, "File %s could not be opened", use_filename);
		free (use_filename);
		RESULT_BUF (buf);
		return TCL_ERROR;
	    }
	    free (use_filename);
	    Model::write (file);
	    Fclose (file);
	    return TCL_OK;
	}
	if (!StringCmp (argv[1], "load", case_ins))
	{
	    char *given_filename = argv[2];
	    char *use_filename = append_extension (given_filename, ".mod");
	    int status = Model::load (use_filename, interp);
	    free (use_filename);
	    return status;
	}
    }

    RESULT_LIT ("Invalid model command");
    return TCL_ERROR;
}

int Model::write (FILE *file)
{
    fprintf (file, "solarmodel %s\n", SolarVersion ());
    Definition::Write_Commands (file);
    Matrix::write_commands (file);
    Trait::Write_Commands (file);
    Parameter::write_commands (file);
    Covariate::write_commands (file);
    Scale::Write_Commands (file);
    Constraint::write_commands (file);
    Omega::write_commands (file);
    Mu::write_commands (file);
    Option::write_commands (file);
    if (Loglike::status())
    {
	char buf[128];
	fprintf (file, "loglike set %s\n", Loglike::get_string (buf));
    }
    return TCL_OK;
}

int Model::load (const char *filename, Tcl_Interp *interp)
{
    Model::reset ();
    FILE *modf = fopen (filename, "r");
    char buf[1024];
    if (!modf)
    {
	sprintf (buf, "Couldn't read file \"%s\": no such file", filename);
	RESULT_BUF (buf);
	return TCL_ERROR;
    }
    if (!fgets(buf, 1024, modf))
    {
	sprintf (buf, "Model file \"%s\" is empty", filename);
	RESULT_BUF (buf);
	return TCL_ERROR;
    }
    if (strncmp (buf, "solarmodel", strlen ("solarmodel")))
    {
	sprintf (buf, 
"Must use upgrade command to upgrade this model: upgrade %s", filename);
	RESULT_BUF (buf);
	return TCL_ERROR;
    }
    fclose (modf);

// Set loading status and reset afterwards

    int success = TCL_ERROR;
    _loading = true;
    try
    {
	success = Solar_EvalFile (interp, filename);
    }
    catch (Safe_Error_Return& ser)
    {
	RESULT_BUF (ser.message());
    }
    _loading = false;
    return success;
}    

int Model::renew (Tcl_Interp *interp)
{
    Model::reset ();
    if ( TCL_OK != Solar_Eval (interp, "Solar_Model_Startup"))
    {
	Model::reset ();
	return TCL_ERROR;
    }
    return TCL_OK;
}

void Model::reset()
{
    Loglike::reset();
    Option::reset();
    Mu::reset();
    Omega::reset();
    Constraint::reset();
    Covariate::reset();
    Scale::Reset();
    Parameter::reset();
    Trait::Reset();
    Matrix::reset();
//  Definition::Reset();
}


char *append_extension (const char *given_filename, const char *extension)
{
    int elen = strlen (extension);
    int glen = strlen (given_filename);
    if (glen >= elen)
    {
	if (!strcmp (&given_filename[glen-elen], extension))
	{
	    return Strdup (given_filename);
	}
    }
    char *filename = Strdup (given_filename);
    return StringAppend (&filename, extension);
}
