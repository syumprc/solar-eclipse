/*
 * phenotypes.cc implementes the phenotypes class and command
 * Written by Charles Peterson beginning on June 3, 1998
 * Copyright (c) 1998 Southwest Foundation for Biomedical Research
 *
 * In 2004 this was modified to allow multiple phenotypes files.  Basically,
 * the multi-file capability is added through an additional level of
 * indirection in required variables, which is referenced through a zero
 * based file index.
 *
 * The data vector pointer, Current_Person_Phenotypes, is now an array
 * of such pointers, one for each phenotypes tablefile.  To translate
 * the parameter index "pindex" into file_index and file_pos, there
 * are File_Index and File_Pos vectors.
 *
 * Look below for "Multifile" notes.
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "solar.h"
// tcl.h from solar.h
#include "safelib.h"


/*
 * For each phenotypes file, there is an index of pointers to each record
 */

class PhenoPointers
{
public:
    char *_id;
    char *_famid;
    int _fposition;
};

class PhenoFileIndex
{
public:
    PhenoFileIndex() {Pointers=0; Count=0;}
    PhenoPointers* Pointers;
    int Count;
    void add (const char *id, const char *famid, int fposition);
    void reset ();
    bool established () {return (0!=Pointers) ? true : false;}

    static void Reset_All ();
};

PhenoFileIndex* Pheno_Index[MAX_PHENOTYPE_FILES] = {0};


int Phenotypes::Filecount = 0;
const char **Phenotypes::PhenoNames[] = {0};
int Phenotypes::PhenoCount[] = {0};
SolarFile *Phenotypes::Sfile[] = {0};
famid_status Phenotypes::_found_famid = not_known;
int Phenotypes::_proband_pindex = -1;
int Phenotypes::File_Index[] = {0};
int Phenotypes::File_Pos[] = {0};
int Phenotypes::File_Useage[] = {0};
int Phenotypes::Setup_Count = 0;
time_t Phenotypes::Last_Modified[] = {0};
Expression* Phenotypes::Expressions[] = {0};
int Phenotypes::Current_Person_Sex = 0;

char Required_By_Trait[MAX_PHENOTYPES][MAX_TRAITS];

char **Phenotypes::Current_Person_Phenotypes[] = {0};
StringArray *Phenotypes::Setup_Names  = 0;
Tcl_Interp* Phenotypes::Interpp = 0;

// Inverse normals (some setup done in trait.cc also)

double* Phenotypes::INData = 0;
char**    Phenotypes::INNames = 0;
int      Phenotypes::INCount = 0;
bool* Phenotypes::INIfclass = 0;
int* Phenotypes::INClass = 0;

// Main command routine

extern "C" int PhenotypesCmd (ClientData clientData, Tcl_Interp *interp,
		  int argc, char *argv[])
{
    if (argc == 2 && !StringCmp ("help", argv[1], case_ins))
    {
	return Solar_Eval (interp, "help phenotypes");
    }

    if (argc == 1)
    {
	return Phenotypes::describe (interp, true);
    }

    if (argc == 2 && !Strcmp ("-files", argv[1]))
    {
	return Phenotypes::showfiles (interp);
    }

// THE "load phenotypes" command

    if (argc == 2 || (argc >= 3 && !StringCmp ("load", argv[1], case_ins)))
    {
    // Check for missing famid error and missing file error
	int start_phen = 1 + !Strcmp ("load", argv[1]);
	int file_count = 0;
	char* freebuf = Strdup ("check_phenotypes");
	for (int i = start_phen; i < argc; i++)
	{
	    StringAppend (&freebuf, " ");
	    StringAppend (&freebuf, argv[i]);
	    file_count++;
	}
	int errstatus = Solar_Eval (interp, freebuf);
	free (freebuf);
	if (errstatus)
	{
		Phenotypes::reset (hard);
		return TCL_ERROR;
	}

    // _found_famid is initially set here (based on tcl check_phenotypes)
    // It gets unset if file modification detected or reset
    //   (That forces re-check during build_index)

	Phenotypes::_found_famid = not_found;
	if (!Strcmp (interp->result,"1"))
	{
	    Phenotypes::_found_famid = found;
	}
	RESULT_LIT ("");

// Setup a new phenotypes file
	const char *errmsg = 0;
	Phenotypes::load_files (&argv[start_phen], file_count, &errmsg);
	if (errmsg)
	{
	    RESULT_BUF (errmsg);
	    return TCL_ERROR;
	}

// Display information
	if (verbose ("PHENOTYPES_LOAD"))
	{
	    return Phenotypes::describe(interp, false);
	}
	return TCL_OK;
    }
    RESULT_LIT ("Invalid phenotypes command");
    return TCL_ERROR;
}

void Phenotypes::reset(ResetMode rm)
{
    for (int i = 0; i < MAX_PHENOTYPE_FILES; i++)
    {
	if (Sfile[i]) delete Sfile[i];
	Sfile[i] = 0;
	PhenoNames[i] = 0;
	PhenoCount[i] = 0;
	Current_Person_Phenotypes[i] = 0;
    }
    Filecount = 0;
    PhenoFileIndex::Reset_All();
    _proband_pindex = -1;
    Setup_Names = 0;
    _found_famid = not_known;
    if (rm == hard)
    {
	unlink ("phenotypes.info");
    }
}

// load_files manages .info file
// actual loading done by call to reload_files

void Phenotypes::load_files (char** filenames, int filecount,
			    const char **errmsg)
{
    Phenotypes::reset (soft);
    reload_files (filenames, filecount, errmsg);
    if (*errmsg)
    {
	Phenotypes::reset (hard);
    }
    else
    {
	bool must_update = true;
	FILE *pnf = fopen ("phenotypes.info", "r");
	if (pnf)
	{
	    must_update = false;
	    char old_filename[1024];
	    int counted = 0;
	    while (fgets (old_filename, 1024, pnf) && strlen (old_filename))
	    {
		if (old_filename[strlen(old_filename)-1] == '\n')
		    old_filename[strlen(old_filename)-1] = '\0';
		if ((counted >= filecount) ||
		    strcmp (filenames[counted], old_filename))
		{
		    must_update = true;
		    break;
		}
		counted++;
	    }
	    if (counted < filecount) must_update = true;
	    fclose (pnf);
	}
	if (must_update)
	{
	    pnf = fopen ("phenotypes.info", "w");
	    if (pnf)
	    {
		for (int i = 0; i < filecount; i++)
		{
		    fprintf (pnf, "%s\n", filenames[i]);
		}
		fclose (pnf);
	    }
	}
    }
}
	
// reload_files does not delete or rewrite .info file
// used by start and load_files
// this fixes a bug in parallel operation from the same working directory

void Phenotypes::reload_files (char** filenames, int nfiles, 
			      const char **errmsg)
{
    *errmsg = 0;
    Phenotypes::reset (soft);
    FisherPedigree::Changing_Pedigree();

    for (int i = 0; i < nfiles; i++)
    {
	Sfile[i] = SolarFile::open ("phenotypes", filenames[i], errmsg);
	if (*errmsg)
	{
	    Phenotypes::reset (soft);
	    return;
	}
	struct stat statbuf;
	if (stat (filenames[i], &statbuf))
	{
	    Phenotypes::reset (soft);
	    *errmsg = "Phenotypes file disappeared disappeared after loading";
	    return;
	}
	Last_Modified[i] = statbuf.st_mtime;
	PhenoNames[i] = Sfile[i]->names (&PhenoCount[i], errmsg);
	if (*errmsg)
	{
	    Phenotypes::reset (soft);
	    return;
	}
	Filecount++;
    }
}

void Phenotypes::start ()
{
    FILE *pnf = fopen ("phenotypes.info", "r");
    const char *errmsg = 0;
    if (pnf)
    {
	Filecount = 0;
	char filename[1024];
	char** filenames = (char**) Malloc (1);
	filename[0] = '\0';
	while (fgets (filename, 1024, pnf) && strlen (filename))
	{
	    if (filename[strlen(filename) - 1] == '\n')
		filename[strlen(filename) - 1] = '\0';

	    filenames = (char**) Realloc (filenames, ++Filecount *
		sizeof (char**));
	    filenames[Filecount-1] = Strdup (filename);
	}
	fclose (pnf);
	Phenotypes::reload_files (filenames, Filecount, &errmsg);
	if (errmsg)
	{
	    fprintf (stderr, "\nError reading phenotypes file %s:\n%s\n\n",
		     filename, errmsg);
	}
	for (int i = 0; i < Filecount; i++)
	{
	    free (filenames[i]);
	}
	free (filenames);
    }
}
	

int Phenotypes::showfiles (Tcl_Interp *interp)
{
    char* buf = Strdup ("");
    for (int i = 0; i < Filecount; i++)
    {
	if (i > 0) string__append (&buf, "  ");
	string__append (&buf, Sfile[i]->filename());
    }
    RESULT_BUF (buf);
    free (buf);
    return TCL_OK;
}
	

int Phenotypes::describe (Tcl_Interp* interp, bool showall)
{
    if (Filecount==0)
    {
	RESULT_LIT ("No phenotypes file is loaded");
	return TCL_ERROR;
    }
    char* buf = Strdup ("");
    for (int i = 0; i < Filecount; i++)
    {
	int totalchars = 0;
	if (i > 0) string__append (&buf, "\n\n");
	string__append (&buf, Sfile[i]->filename());
	string__append (&buf, ": ");
	for (int j = 0; j < PhenoCount[i]; j++) 
	{
	    string__append (&buf, PhenoNames[i][j]);
	    string__append (&buf, " ");
	    if (!showall)
	    {
		totalchars += strlen (PhenoNames[i][j]);
		if (totalchars > 150)
		{
		    string__append (&buf, "...\n");
		    string__append (&buf, "[you may show all phenotypes using the \'phenotypes\' command]\n");
		    break;
		}
	    }
	}
    }
    RESULT_BUF (buf);
    free (buf);
    return TCL_OK;
}

void Phenotypes::write_commands (FILE *file)
{
    if (Sfile[0])
    {
	fprintf (file, "load phenotypes");
	for (int i = 0; i < Filecount; i++)
	{
	    fprintf (file, " %s", Sfile[i]->filename());
	}
	fprintf (file, "\n");
    }
}

bool Phenotypes::available (const char *name)
{
    if (Sfile[0])
    {
	const char *errmsg = 0;
	bool found = false;
	for (int i = 0; i < Filecount; i++)
	{
	    found = Sfile[i]->test_name (name, &errmsg);
	    if (errmsg) return false;
	    if (found) return found;
	}
    }
    return false;
}

int Phenotypes::bind (Tcl_Interp *interp)
{
    if (!Sfile || !PhenoNames)
    {
	RESULT_LIT ("Phenotypes file has not been loaded");
	return TCL_ERROR;
    }
    for (int i = 0; i < MAX_PHENOTYPES; i++)
    {
	for (int j = 0; j < Trait::Number_Of(); j++)
	{
	    Required_By_Trait[i][j] = 0;
	}
    }
    return build_index (interp);
}

int Phenotypes::build_index (Tcl_Interp *interp)
{

// build_index builds pointers to each record of each phenotypes
// file, AND checks the discrete trait coding.

// (Assuming the coding is correct, tests in preped/pinput actually
//  do the effective determination.)

// See if index(es) already built and files haven't changed since then

    int ifile;
    bool index_invalid = false;
    bool needs_reload = false;
    for (ifile = 0; ifile < Filecount; ifile++)
    {
	struct stat statbuf;
	if (stat (Sfile[ifile]->filename(), &statbuf) ||
	    statbuf.st_mtime != Last_Modified[ifile])
	{
	    index_invalid = true;
	    needs_reload = true;
	    _found_famid = not_known;
	} 
	else if (!Pheno_Index[ifile])
	{
	    index_invalid = true;
	}
    }

// If reload needed due to file modification, do it now

    if (needs_reload)
    {
	char** filenames = (char**) Calloc (Filecount+1, sizeof (char*));
	int fin;
	for (fin = 0; fin < Filecount; fin++)
	{
	    filenames[fin] = Strdup (Sfile[fin]->filename());
	}
	const char* errmsg = 0;
	reload_files (filenames, Filecount, &errmsg);
	if (errmsg)
	{
	    RESULT_BUF (errmsg);
	    return TCL_ERROR;
	}
	for (int fin = 0; fin < Filecount; fin++)
	{
	    free (filenames[fin]);
	}
	free (filenames);
	if (errmsg)
	{
	    RESULT_BUF (errmsg);
	    return TCL_ERROR;
	}
    }


// check famid status if unchecked or if files changed

    if (_found_famid == not_known)
    {
	char* freebuf = Strdup ("check_phenotypes");
	for (int i = 0; i < Filecount; i++)
	{
	    StringAppend (&freebuf, " ");
	    StringAppend (&freebuf, Sfile[i]->filename());
	}
	int errstatus = Solar_Eval (interp, freebuf);
	free (freebuf);
	if (errstatus)
	{
	    return TCL_ERROR;
	}
	_found_famid = not_found;
	if (!Strcmp (interp->result,"1"))
	{
	    _found_famid = found;
	}
	RESULT_LIT ("");
    }



    if (index_invalid || Trait::Maximization_Type()==unknown)
    {

	const char *errmsg = 0;
	bool message_shown = false;
	int ntraits = Trait::Number_Of();
	bool trait_found[MAX_TRAITS];
	bool trait_found_anywhere[MAX_TRAITS];
	int itrait;

	int values_seen[MAX_TRAITS];
	double values[MAX_TRAITS][2];

	for (itrait = 0; itrait < ntraits; itrait++)
	{
	    trait_found_anywhere[itrait] = false;
	    values_seen[itrait] = 0;
	}

	for (ifile = 0; ifile < Filecount; ifile++)
	{
	    if (!message_shown && verbose ("ITERATE"))
	    {
		cfprintf (stdout, 
		  "\nBuilding index for phenotypes file...\n\n","");
		message_shown = true;
	    }
	    
	    if (Pheno_Index[ifile])
	    {
		Pheno_Index[ifile]->reset ();
	    }
	    else
	    {
		Pheno_Index[ifile] = new PhenoFileIndex;
	    }
		
	    int famid_index = 1;

	    Sfile[ifile]->start_setup (&errmsg);
	    Sfile[ifile]->setup ("id", &errmsg);  // Index 0
	    for (int i = 0; i < ntraits; i++)
	    {
		if (!Trait::expression(i) &&
		    Sfile[ifile]->test_name (Trait::Name(i), &errmsg))
		{
		    Sfile[ifile]->setup (Trait::Name(i), &errmsg);
		    if (errmsg)
		    {
			RESULT_BUF(errmsg);
			return TCL_ERROR;
		    }
		    trait_found[i] = true;
		    trait_found_anywhere[i] = true;
		    famid_index++;
		}
		else
		{
		    trait_found[i] = false;
		}
	    }
	    if (Sfile[ifile]->test_name ("FAMID", &errmsg) && found_famid() )
	    {
		Sfile[ifile]->setup ("FAMID", &errmsg); // Index 1 + ntraits
	    }
	    else
	    {
		if ( found_famid() ) 
		{
		    Phenotypes::start ();  // Reset to starting state
		    RESULT_LIT (
		    "FAMID mapping changed; check that 'field famid' matches phenotypes file");
		    return TCL_ERROR;
		}
	    }

	    Sfile[ifile]->rewind (&errmsg);

	    if (errmsg)  // If _not_ EOF
	    {
		char message[1024];
		{
		    sprintf (message, 
			     "Error reading phenotypes file %s: %s\n", 
			     Sfile[ifile]->filename(), errmsg);
		}
		RESULT_BUF (message);
		Phenotypes::start ();  // Reset to starting state
		return TCL_ERROR;
	    }

	    int offset;
	    char **data;
	    int reccnt = 0;

// Scan the Phenotypes file for two separate reasons:
//    Get file pointers to the beginning of each record
//    Check to see whether there are actually more than two
//      values for trait.  This helps check for improper
//      coding.
//
// We're not setting up covariates or even traits now.  That's
// done later in ccsearch using Phenotypes::start_setup

	    while (0 != (data = Sfile[ifile]->get (&errmsg)))
	    {

// If we've seen more than 2 values, trait must be quantitative
		int idata = -1;
		for (itrait = 0; itrait < ntraits; itrait++)
		{
		    if (trait_found[itrait])
		    {
			idata++;
			if (unknown == Trait::Type (itrait))
			{
			    if (strlen (data[idata+1]))
			    {
				double trait = atof (data[idata+1]);
				if (values_seen[itrait] == 2)
				{
				    if (trait != values[itrait][0] && 
					trait != values[itrait][1])
				    {
					Trait::Set_Quantitative (itrait);
				    }
				}
				else if (values_seen[itrait] == 1)
				{
				    if (values[itrait][0] != trait)
				    {
					values[itrait][1] = trait;
					values_seen[itrait] = 2;
				    }
				}
				else if (values_seen[itrait] == 0)
				{
				    values[itrait][0] = trait;
				    values_seen[itrait] = 1;
				}
			    }
			}
		    }
		}

// Now actually add this record to index
		int position = Sfile[ifile]->get_position ();
		if (found_famid())
		{
		    Pheno_Index[ifile]->add 
			(data[0], data[famid_index], position);
		}
		else
		{
		    Pheno_Index[ifile]->add (data[0], "", position);
		}
		reccnt++;
	    }

	    Sfile[ifile]->rewind (&errmsg);
	    if (errmsg)  // If _not_ EOF
	    {
		char message[1024];
		{
		    sprintf (message, 
	            "Error reading file %s (data record %d):\n  %s\n", 
			 Sfile[ifile]->filename(), reccnt+1, errmsg);
		}
		RESULT_BUF (message);
		Phenotypes::start ();  // Reset to starting state
		return TCL_ERROR;
	    }
	}

// Now that we've seen all the values, determine if each trait is binary
// and if coded correctly

	for (itrait = 0; itrait < ntraits; itrait++)
	{
	    if (trait_found_anywhere[itrait])
	    {
		if (unknown == Trait::Type(itrait))
		{
		    if (values_seen[itrait] > 1)
		    {
			if (fabs(values[itrait][0] - values[itrait][1]) == 1.0
			    && (ceil(values[itrait][0]) == values[itrait][0]))
			{
			    Trait::Set_Discrete (itrait);
			}
			else
			{
			    printf (" ** Warning. Only two trait %d values but not coded as discrete\n\n",itrait+1);
			    printf (" ** Trait %d will be considered quantitative\n",itrait+1);
			    Trait::Set_Quantitative (itrait);
			}
		    }
		    else
		    {			
			int soption = Option::get_int("SingularTrait");
			if (soption == 0)
			{
			    char buf[1024];
			    sprintf (buf,"Only one non-blank value for trait %d found in phenotypes file!\nIf this is intended, set option SingularTrait to permit this.\n   1 for discrete, 2 for quantitative.\nSee 'help option' for more details.",itrait+1);
			    RESULT_BUF(buf);
			    return TCL_ERROR;
			}
			else if (soption == 1)
			{
			    Trait::Set_Discrete (itrait);
			}
			else
			{
			    Trait::Set_Quantitative (itrait);
			}
		    }
		}
	    }
	}

// Ensure that trait was seen in one file or another
// This should only be caught here if phenotypes file changed after
// the trait command was given.

	if (message_shown)
	{
	    for (itrait = 0; itrait < ntraits; itrait++)
	    {
		if (!Trait::expression(itrait) &&
		    !trait_found_anywhere[itrait])
		{
		    char message[1024];
		    sprintf (message, 
			     "Trait %d not found in any phenotypes file\n");
		    RESULT_BUF(message);
		    Phenotypes::start();
		    return TCL_ERROR;
		}
	    }
	}

// For expression traits, we don't test all the values, that would be very
// complex since multiple files would have to be indexed simultaneously.  
// However, we check to see if the top level operator is boolean,
// and let that determine whether the trait is considered discrete or not.
// User can simply add a top level inequality (if not using one already)
// to satisfy this requirement.

	for (itrait = 0; itrait < ntraits; itrait++)
	{
	    Expression *exp = Trait::expression(itrait);
	    if (exp)
	    {
		if (exp->boolean())
		{
		    Trait::Set_Discrete(itrait);
		}
		else
		{
		    Trait::Set_Quantitative(itrait);
		}
	    }
	}
    }

    Phenotypes::Interpp = interp;
    
    return TCL_OK;
}


extern "C" double getpheno_ (int *pindex)
{
    const char *phenotype = Phenotypes::get_indexed_phenotype (*pindex-1);
    if (!phenotype)
    {
	throw Safe_Error_Return ("Error accessing indexed phenotype");
    }
    double value;
    if (!strlen (phenotype))
    {
	value = MISSING_PHENOTYPE;
    }
    else
    {
	char junk[128];
	int count;
	count = sscanf (phenotype, "%lg %32c", &value, junk);
	if (count != 1)
	{
	    char buf[512];
	    char phen[128];
	    strncpy (phen, phenotype, 80);
// convert fortran D exponential format
	    char* dpos = strchr (phen, 'D');
	    if (dpos)
	    {
		*dpos = 'e';
		count = sscanf (phen, "%lg %32c", &value, junk);
		if (count != 1)
		{
		    phen[79] = 0;
		    sprintf (buf, "Non-numeric data found in phenotype field:  %s\n(Note: blank must be used to specify missing data; see help file-phenotypes)",
		     phen);
		    throw Safe_Error_Return (buf);
		}
	    }
	}
    }
    return value;
}

const char *Phenotypes::test_proband ()
{
    if (_proband_pindex < 0) return "0";
    int findex = File_Index[_proband_pindex];
    int pindex = File_Pos[_proband_pindex];
    if (!Current_Person_Phenotypes[findex])
    {
	return "0";  // If no vars available, assume not a proband
    }
    const char *probvar = Current_Person_Phenotypes[findex][pindex];
    if (!strlen (probvar))
    {
	return "0";
    }
    double testval = 0;
    int converted = sscanf (probvar, "%lg", &testval);
    if (converted == 1 && testval == 0.)
    {
	return "0";
    }
    return "1";
}


// The following methods used by expression traits

double eval_phenotype (int eqindex)
{

// The following error should never occur.
// However if it does, it is probably related to a trait specific
// covariate being accessed when that trait is not present.  See
// definition context binding section in Covariate::bind where
// an assumption is violated.

    if (eqindex < 0) {
	error ("Internal inconsistency in accessing defined covariate value\n\
Please report this error to solar@txbiomedgenetics.org");
    }

    int pindex = EqVar::get_packed_index (eqindex) + 1; // Fortranize
    double test = getpheno_ (&pindex);
    if (test == MISSING_PHENOTYPE) {
	throw Missing_Data();
    }
    return test;
}

double eval_sex (int eqindex)
{
    return Phenotypes::Current_Person_Sex;
}


// The following method used to perform inormal_ transformation

double get_inormal (int inindex)
{
    double value = Phenotypes::INData[inindex];
    if (value == MISSING_PHENOTYPE) {
	throw Missing_Data();
    }
    return value;
}

double get_zscore (int zindex)
{
    Zscore* z = Zscore::Get(zindex);
    double phen = eval_phenotype (z->vindex);
    if (!z->meanvalid)
    {
	return phen;
    }
    double zscore = (phen - z->mean) / z->sd;
    return zscore;
}

static int etraits = 0;  // See below for lengthy rationale

const char *Phenotypes::get_indexed_phenotype (int pindex)
{
    int ntrait = Trait::Number_Of();
    if (pindex < ntrait)
    {

// The way this is going to be used by caller pinput is we are going to
// always start with trait 0, then trait 1, and so on.  We're never going to
// jump to trait N.  Some of the traits may actually be expressions, in which
// case they aren't actually "variables."  For any given trait I, we need to
// know the number of expression traits which preceded it in order to set the
// index into phenotypes vector accordingly.
//
// If at some point in the future it becomes necessary to permit jumping
// straight to trait I, we'll either have to explicitly loop through all
// preceding traits or store a separate information vector of some kind.
//
// The virtual indexes are the ones maintained by Setup_Names.  See note
// before setup_virtual_position.

	if (pindex == 0) etraits = 0;

	Expression* exp = Trait::expression(pindex);
	if (exp)
	{
	    etraits++;
	    try
	    {
		double result = exp->eval();
		static char buf[256];
		sprintf (buf, "%18.12e", result);
		return buf;
	    }
	    catch (...)
	    {
		return MISSING_PHENOTYPE_STRING;
	    }
	}
	else
	{
// Now virtual traits are included in the virtual index

	    int findex = File_Index[pindex];
	    if (!Current_Person_Phenotypes[findex])
	    {
		return MISSING_PHENOTYPE_STRING;
	    }
	    return Current_Person_Phenotypes [findex]
	    [File_Pos[pindex]];
	}
    }
    else if (pindex == ntrait) // IBDID
    {
	throw Safe_Error_Return ("Attempted to access IBDID pseudovariable");
    }
    else if (pindex == ntrait+1) // Proband test
    {
	return test_proband ();
    }
    else
    {

// The following adjustments are now obsolete due to changes in
// setup_virtual_position
//
//	int cpindex = pindex - Trait::Number_Of_Expressions();
// Adjust for old IBDID placeholder
//	cpindex--;
// Adjust for absence of proband
//	if (_proband_pindex < 0)
//	{
//	    cpindex--;
//	}
//


// Get Covariates that are defined as expressions

	if (Expressions[pindex])
	{
	    try
	    {
		double result = Expressions[pindex]->eval();
		static char buf[256];
		sprintf (buf, "%18.12e", result);
		return buf;
	    }
	    catch (...)
	    {
		return MISSING_PHENOTYPE_STRING;
	    }
	}		

// Now retrieve data from tablefile vector

	int findex = File_Index[pindex];
	if (!Current_Person_Phenotypes[findex])
	{
	    return MISSING_PHENOTYPE_STRING;
	}
	return Current_Person_Phenotypes [findex]
	    [File_Pos[pindex]];
    }
}

void Phenotypes::seek (const char *id, const char *famid)
{
    const char* errmsg = 0;
    int j;
    for (j = 0; j < Filecount; j++)
    {
	Current_Person_Phenotypes[j] = 0;
	PhenoFileIndex* pf = Pheno_Index[j];
	for (int i = 0; i < pf->Count; i++)
	{
	    if (strcmp (pf->Pointers[i]._id, id)) continue;
	    if ( found_famid() )
	    {
		if (strcmp (pf->Pointers[i]._famid, famid)) continue;
	    }

// Found position, now "set" it

	    Sfile[j]->set_position (pf->Pointers[i]._fposition, &errmsg);
	    Current_Person_Phenotypes[j] = Sfile[j]->get (&errmsg);
	    if (errmsg)
	    {
		char buf [256];
		sprintf (buf, "Phenotype file error: %s", errmsg);
		throw Safe_Error_Return (errmsg);
	    }
	}
    }

// Load this individual's inverse normals if applicable

    if (Phenotypes::INData)
    {
	free (Phenotypes::INData);
	Phenotypes::INData = 0;
    }
    if (Phenotypes::INCount)
    {
	Phenotypes::INData = (double*) Calloc (Phenotypes::INCount+1, sizeof (double));
    }	
    for (j = 0; j < Phenotypes::INCount; j++)
    {
	char cbuf[2048];
	if (strlen (famid))
	{
	    sprintf (cbuf, "inormal -trait %s -id \"%s\" -famid \"%s\"", Phenotypes::INNames[j], id, famid);
	}
	else
	{
	    sprintf (cbuf, "inormal -trait %s -id \"%s\"", Phenotypes::INNames[j], id);
	}
	if (Phenotypes::INIfclass[j])
	{
	    sprintf (&cbuf[strlen(cbuf)], " -class %d\n", Phenotypes::INClass[j]);
	}

	int tclerr = Solar_Eval (Phenotypes::Interpp, cbuf);
	if (tclerr)
	{
	    sprintf (cbuf, "\nmaximize: Error getting inverse normal for individual %s\n", id);
	    printf (cbuf);  // Safe_Error_Return doesn't always work here
	    throw Safe_Error_Return (cbuf);
	}
	const char* instring = Tcl_GetStringResult (Phenotypes::Interpp);
	double invalue = 0;
	if (strlen(instring) == 0)
	{
	    Phenotypes::INData[j] = MISSING_PHENOTYPE;
	}
	else
	{
	    int converted = sscanf (instring, "%lf %c", &invalue, cbuf);
	    
	    Tcl_SetResult (Phenotypes::Interpp, (char*) "", TCL_STATIC);

	    if (converted != 1)
	    {
		sprintf (cbuf, "\nmaximize: Error reading inverse normal for individual %s\n", id);
		throw Safe_Error_Return (cbuf);
	    }
//	    printf ("For ID %s we got %g\n", id, invalue);

	    Phenotypes::INData[j] = invalue;
	}
    }
}

void Phenotypes::start_setup ()
{
    int i;
    if (Filecount==0) throw Safe_Error_Return ("Phenotype file not opened");
    const char *errmsg = 0;
    for (i=0; i < Filecount; i++)
    {
	Sfile[i]->start_setup (&errmsg);
	if (errmsg)
	{
	    char buf [256];
	    sprintf (buf, "Phenotype file error: %s", errmsg);
	    throw Safe_Error_Return (errmsg);
	}
	File_Useage[i] = 0;
    }
    for (i=0; i < MAX_PHENOTYPES; i++)
    {
	File_Index[i] = 0;
	File_Pos[i] = 0;
    }
    Setup_Count = 0;
    if (Setup_Names) delete Setup_Names;
    Setup_Names = new StringArray;
    _proband_pindex = -1;
}


void Phenotypes::setup_proband (const char *name)
{
    setup (name);  // setup as a "real" phenotype
    _proband_pindex = Setup_Count-1; // Setup_Count was incremented by setup
}

/*
 * Multifile note:
 *   When setting up a phenotype, we must do the following things:
 *     1) Find the file(s) including this phenotype
 *     2) If this name is found in multiple files, throw a fatal error
 *     3) Setup the phenotype name in specific tablefile
 *     4) Add to index vectors File_Index and File_Position so that
 *        pindex can later be used to access this field
 *     5) Update File_Useage and Setup_Count to allow for later vector
 *        updates.
 *     6) Add to "Setup_Names" list of setup phenotypes.  This list does
 *        not account for the specific file phenotype came from, but reflects
 *        which phenotypes have been setup.
 */

void Phenotypes::setup (const char *name)
{
    if (Filecount==0) throw Safe_Error_Return ("Phenotype file not opened");
    const char *errmsg = 0;
    int found = 0;
    int i;
    char* justname = Strdup (name);
    char* colonpos = strchr (justname, ':');
    if (colonpos) *colonpos = '\0';

    for (i = 0; i < Filecount; i++)
    {
        if (Sfile[i]->test_name (justname, &errmsg) && !errmsg)
	{
	    if (found)
	    {
	        char buf[1024];
		sprintf (buf, "Phenotype named %s found in several files",
			 justname);
		free (justname);
	        throw Safe_Error_Return (buf);
	    }
	    Sfile[i]->setup (justname, &errmsg);
	    found = i+1;
	}
	if (errmsg)
	{
	    char buf [1024];
	    sprintf (buf, "Phenotype file error: %s", errmsg);
	    free (justname);
	    throw Safe_Error_Return (buf);
	}
    }
    if (!found)
    {
	char buf[1024];
	sprintf (buf, "Phenotype name %s not found", justname);
	free (justname);
	throw Safe_Error_Return (buf);
    }
    free (justname);
    File_Index[Setup_Count] = found - 1;
    File_Pos[Setup_Count] = File_Useage[found - 1];
    File_Useage[found-1]++;
    setup_virtual_position (name);
}

void Phenotypes::Covar_Of_Trait (int covindex, int traitindex)
{
    Required_By_Trait[covindex][traitindex] = 1;
}

extern "C" int traitneeds_ (int* traitindex, int* covindex)
{
    return (int) Required_By_Trait[(*covindex)-1][(*traitindex)-1];
}


// One thing that has become tricky since the advent of expression traits
// is that there are two indexes used.  One, represented by the "Setup_Count"
// represents the next index into File_Index and File_Pos to get to the actual
// data.  The other, represented by Setup_Names, is the virtual index used
// by pindex/getpheno_/get_indexed_phenotypes to access the all phenotypes,
// including virtual ones.  This is actually managed by Setup_Names.  Names
// are (potentially) found in this list by scanning through them.
//
// Therefore, all names, even virtual ones, are added to the Setup_Names
// with setup_virtual_position.
//
// HOWEVER, the expression trait is added to this list as a blank in case
// the trait is set up as a unary expression of a single phenotype.  For
// example, define a = q4.  In that case, we don't want to "find" q4 on this
// list when actually adding q4 the phenotype later.  Perhaps confusingly, the
// phenotype is then reported (in maximization output) twice, once as the
// expression and once as the actual phenotype.  But at least it "works" even
// in this obscure test case.

void Phenotypes::setup_virtual_position (const char* name)
{
    Setup_Names->add (name);
    Setup_Count++;
}

const char* Phenotypes::get_var_name (int index)
{
    return Setup_Names->index (index);
}

int Phenotypes::find_virtual_position (const char *name)
{
    return Setup_Names->find (name);
}

const char* Phenotypes::filenames ()
{
    static char filenamebuf[MAX_PHEN_FILE_CHARS];
    strcpy (filenamebuf, "");
    if (Sfile[0])
    {
	for (int i = 0; i < Filecount; i++)
	{
	    if (i) strncat (filenamebuf, "  ", MAX_PHEN_FILE_CHARS);
	    strncat (filenamebuf, Sfile[i]->filename(), MAX_PHEN_FILE_CHARS);
	}
    }
    return filenamebuf;
}

void PhenoFileIndex::add (const char *id, const char *famid, int fposition)
{
    if (0==Pointers)
    {
	Pointers = (PhenoPointers*) Calloc (1,sizeof(PhenoPointers));
	Count = 1;
    } else {
	int newsize = (++Count) * sizeof (PhenoPointers);
//	fprintf (stderr, "\nNew size is %d\n", newsize);
	Pointers = (PhenoPointers*) Realloc (Pointers, newsize);
    }
    Pointers[Count-1]._id = Strdup (id);
    Pointers[Count-1]._famid = Strdup (famid);
    Pointers[Count-1]._fposition = fposition;
}


void PhenoFileIndex::reset ()
{
    if (Pointers)
    {
	for (int i = 0; i < Count; i++)
	{
	    
	    free (Pointers[i]._id);
	    free (Pointers[i]._famid);
	}
	free (Pointers);
	Pointers = 0;
    }
    Count = 0;
}

void PhenoFileIndex::Reset_All ()
{
    for (int i = 0; i < MAX_PHENOTYPE_FILES; i++)
    {
	if (Pheno_Index[i])
	{
	    Pheno_Index[i]->reset();
	    Pheno_Index[i] = 0;
	}
    }
}

void Phenotypes::Initialize_Expressions ()
{
    for (int i = 0; i < MAX_PHENOTYPES; i++)
    {
	Expressions[i] = 0;
    }
}

void Phenotypes::Setup_Expression (Expression* expression, int nvar)
{
    if (nvar >= MAX_PHENOTYPES) {
	fprintf (stderr, "More than %d phenotypes\n", MAX_PHENOTYPES);
	exit (10);
    }
    Expressions[nvar] = expression;
}

void StringArray::add (const char* string)
{
    if (++_count>_size)
    {
	_array=(char**)Realloc(_array,sizeof(char*)*(++_size));
    }
    _array[_count-1]=Strdup(string);
}

