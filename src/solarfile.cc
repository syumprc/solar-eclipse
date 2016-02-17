//
// Filename: solarfile.cc
// Purpose:  implements class SolarFile which hides TableFile hierarchy
//           and adds field name mapping.
//
// Written:  February 5, 1999  Charles Peterson
// Copyright (c) 1999, Southwest Foundation for Biomedical Research
//
// Note: As with TableFile, any error causes lock-up.
//       No further action can be taken, except rewind (which clears
//       EOF, and EOF only) or close.


#include "solar.h"

SolarFile* SolarFile::open (const char *file_desc, const char *filename, 
			    const char **errmsg)
{
    SolarFile *sf = new SolarFile (file_desc);
    sf->tf = TableFile::open (filename, errmsg);
    if (sf->tf) {
	return sf;
    }
    delete sf;
    return 0;
}

int SolarFile::setup (const char *generic_name, const char **errmsg)
{
    if (0 != (*errmsg = _errmsg)) return 0;

    const char *use_name = establish_name (generic_name, errmsg);
    if (*errmsg) {
	_errmsg = *errmsg;
	return 0;
    }

    int spos = tf->setup (use_name, errmsg);
    if (*errmsg) {
	_errmsg = *errmsg;
	return 0;
    }
    return spos;
}


// test_name does not set error if name not found; it simply returns false
// However, if there is already an error, it always returns false.

bool SolarFile::test_name (const char *generic_name, const char **errmsg)
{
    if (0 != (*errmsg = _errmsg)) return false;

    const char *local_errmsg = 0;
    const char *use_name = establish_name (generic_name, &local_errmsg);

    if (!local_errmsg)
    {
	bool found = tf->test_name (use_name, errmsg);
	_errmsg = *errmsg;
	return found;
    }
    return false;
}

const char *SolarFile::establish_name (const char *generic_name, 
				       const char **errmsg)
{
    const char *local_errmsg;

    if (!Strcmp (generic_name, "id") ||
	!Strcmp (generic_name, "fa") ||
	!Strcmp (generic_name, "mo") ||
	!Strcmp (generic_name, "famid") ||
	!Strcmp (generic_name, "mztwin") ||
	!Strcmp (generic_name, "hhid") ||
	!Strcmp (generic_name, "sex"))
    {

//      if (!Strcmp (generic_name, "sex")) {

	const char *description = "Sex";
	const char *secondary_name = "gender";
	const char *both_names = "Sex or Gender";

	if (!Strcmp (generic_name, "id"))
	{
	    description = "Ego ID";
	    secondary_name = "ego";
	    both_names = "ego or id";

	} else if (!Strcmp (generic_name, "fa")) {

	    description = "Father ID";
	    secondary_name = "sire";
	    both_names = "fa or sire";

	} else if (!Strcmp (generic_name, "mo")) {

	    description = "Mother ID";
	    secondary_name = "dam";
	    both_names = "mo or dam";

	} else if (!Strcmp (generic_name, "famid")) {

	    description = "Family ID";
	    secondary_name = "";
	    both_names = "famid";
	    
	} else if (!Strcmp (generic_name, "mztwin")) {

	    description = "MZ Twin";
	    secondary_name = "";
	    both_names = "mztwin";

	} else if (!Strcmp (generic_name, "hhid")) {

	    description = "Household ID";
	    secondary_name = "";
	    both_names = "hhid";

	}

// First see if a field command is in effect, if so, use that name only

	const char *established_name = Field::Map (generic_name);
	if (Field::Changed (generic_name))
	{
	    if (!tf->test_name (established_name, &local_errmsg))
	    {
		sprintf (error_message, "\n\
Field not found in %s file for <%s>\n\
Expected field name '%s'\n\
field command mapping <%s> to '%s' is in effect\n\
If you did not enter a field command, it might be from .solar file or script\n\
You can change this with another field command, e.g.:\n\
    solar> field %s newname"
			 ,file_description, description, established_name,
			 description, established_name, generic_name);
		*errmsg = error_message;
		return 0;
	    }
	    return established_name;
	}

// Next try unmapped primary and secondary names

	bool primary_found = tf->test_name (generic_name, &local_errmsg);
	bool secondary_found = false;
	if (secondary_name[0] != '\0')
	{
	    secondary_found = tf->test_name (secondary_name, &local_errmsg);
	}
	if (primary_found && secondary_found)
	{
	    sprintf (error_message, "\n\
Multiple matches for <%s> are present in %s file:\n\
    '%s' \n\
You may use 'field' command to map <%s> to a particular name, e.g:\n\
    solar> field %s newname",
		 description, file_description, both_names, 
		 description, generic_name);
	    *errmsg = error_message;
	    return 0;
	}
	if (!primary_found && !secondary_found)
	{
	    sprintf (error_message, "\n\
Field not found in %s file for <%s>\n\
Expected default field name '%s' \n\
You may use 'field' command to map <%s> to another name, e.g.:\n\
    solar> field %s newname",
		 file_description, description, both_names,
		 description, generic_name);
	    *errmsg = error_message;
	    return 0;
	}
	if (primary_found) return generic_name;
	return secondary_name;
    }
    return generic_name;
}

