// tablefile.cc...implements table file classes
// See tablefile.h for more information

#include <strings.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>

#include "solar.h"
#include "tablefile.h"

#ifndef NO_RICVOLUMESET
#define RICVOLUMESET
#include "RicVolumeSet.h"
#endif


// ** Binary parameters from command line not yet implemented

// Binary parameters can be saved for every setup position
//   They are only actually applied to binary datafield objects (>filename),
//   and not tested otherwise.  If the parameters are missing for a
//   binary datafield object, and error is raised.
//
// Bparam is actually an internally managed list of binary parameter objects.
// It is not allocated until binary parameters are actually created (by
// the setup_binary method).
//
// Assumption is made that there are relatively few (<1000) binary objects
// included as variables, in fact, typically, there will only be one.
//
// Not much optimization is done now.  They could be sorted each
// time a new one is added, reducing the need to scan the list for each lookup.

#if 0
class Bparam
{
    BparamE* bb;
    BparamE* be;
public:
    Bparam() {be=bb=0;};
    add(pos,i,j,k,d) {be

    Bparam()int ix,int iy,int iz, int it=0) {x=iz;y=iy;z=iz;t=it;next=0};
    Bparam
    int setpos;
    int x, y, z, t;
    Bparam* next;
};
#endif


// FISHER2 no longer supported

bool Fisher_2_Format_Allowed;

// Local classes

class BeginWidth
{
public:
    int begin;
    int width;
};

class CodeRecord  // Matches PEDSYS code file
{
public:
    char width_field[2];
    char sep1;
    char label[21];
    char sep2;
    char m1[6];
    char sep3;
    char m2[6];
    char sep4;
    char m3[6];
    char sep5;
    char m4[6];
    char sep6;
    char data_type;
    char nl;
    int width;
};

const char *Standard_Pedsys_Mnemonics[] = {
"GENO  ",
"PHENO ",
"G'TYPE",
"P'TYPE",
"INFERD",
0};

class PedsysFile : public TableFile
{
    friend TableFile *TableFile::open (const char*, const char**);
protected:
// Dynamic arrays
    BeginWidth *user_setup;
    CodeRecord **code_recs;
    virtual void load_names (const char **errmsg);
    virtual const char *read_code_file (FILE *code_file);
    int *_nominal_widths;
    char user_buffer[MAX_REC_BYTES];
    char *user_array[MAX_REC_FIELDS+1];
public:
    PedsysFile (const char *fname) : TableFile () 
	{user_setup=0; code_recs=0; _filename=Strdup(fname); _nominal_widths=0;}
    virtual ~PedsysFile ();

    virtual int *widths (int *count, const char **errmsg);
    virtual int *nominal_widths (int *count, const char **errmsg);
    void start_setup (const char **errmsg);
    virtual int setup (const char *name, const char **errmsg);
    virtual char **get(const char **errmsg);
    virtual void rewind (const char **errmsg);
};

class CommaDelimitedFile : public TableFile
{
    int *user_indexes;
    int highest_user_index;
    int data_starting_position;
    void load_names(const char **errmsg);
    const char *read_header ();
    friend TableFile *TableFile::open (const char*,const char**);
    char* prestring;
    bool getwidths;
    char** user_array;
    int *table_pointers;       // indexes into buffer
    int table_pointer_count;
    char **record_buffers;       // array of buffers
    int record_buffer_count;
     int suggested_bufsize;
    char** freelist;
    size_t freelistcnt;
    short* _types;
public:
    CommaDelimitedFile (const char *fname) : TableFile () 
      {user_indexes=0; highest_user_index=-1; 
	  _filename=Strdup(fname); data_starting_position=0;
          prestring=0; getwidths=false; user_array=0;
	  table_pointers=0;table_pointer_count=0;record_buffers=0;
	  record_buffer_count=0;suggested_bufsize=10000;
          freelist=0;freelistcnt=0;_types=0;}
    ~CommaDelimitedFile ();
    int *widths (int *count, const char **errmsg);
    void start_setup (const char **errmsg);
    int setup (const char *name, const char **errmsg);
    char **get(const char **errmsg);
    void rewind (const char **errmsg);
#ifdef RICVOLUMESET
    static RicVolumeSet* GlobalBin;
    static char* GlobalBinFilename;
#endif
};

#ifdef RICVOLUMESET
char* CommaDelimitedFile::GlobalBinFilename = 0;
RicVolumeSet* CommaDelimitedFile::GlobalBin = 0;
#endif

// C function stubs

extern "C" void *tab_open (const char *filename, const char **errmsg)
{
    return (void *) TableFile::open (filename, errmsg);
}

extern "C" const char **tab_names (void *tabp, int *count, const char **errmsg)
{
    TableFile *ft = (TableFile*) tabp;
    return ft->names (count, errmsg);
}

extern "C" const char **tab_short_names (void *tabp, int *count, const char **errmsg)
{
    TableFile *ft = (TableFile*) tabp;
    return ft->short_names (count, errmsg);
}

extern "C" int tab_test_name (void *tabp, const char *name, const char **errmsg)
{
    TableFile *ft = (TableFile*) tabp;
    return ft->test_name (name, errmsg);
}

extern "C" int *tab_widths (void *tabp, int *count, const char **errmsg)
{
    TableFile *ft = (TableFile*) tabp;
    return ft->widths (count, errmsg);
}

extern "C" void tab_start_setup (void *tabp, const char **errmsg)
{
    TableFile *ft = (TableFile*) tabp;
    ft->start_setup (errmsg);
}

extern "C" int tab_setup (void *tabp, const char *name, const char **errmsg)
{
    TableFile *ft = (TableFile*) tabp;
    return ft->setup (name, errmsg);
}

extern "C" char **tab_get (void *tabp, const char **errmsg)
{
    TableFile *ft = (TableFile*) tabp;
    return ft->get (errmsg);
}

extern "C" void tab_rewind (void *tabp, const char **errmsg)
{
    TableFile *ft = (TableFile*) tabp;
    ft->rewind (errmsg);
}

extern "C" void tab_close (void *tabp)
{
    TableFile *ft = (TableFile*) tabp;
    delete ft;
}

extern "C" const char *tab_filename (void *tabp)
{
    TableFile *ft = (TableFile*) tabp;
    return ft->filename ();
}

// extern "C" int tab_get_position (void *tabp)
// {
//     TableFile *ft = (TableFile*) tabp;
//     return ft->get_position ();
// }

extern "C" void tab_set_position (void *tabp, int pos, const char **errmsg)
{
    TableFile *ft = (TableFile*) tabp;
    ft->set_position (pos, errmsg);
    return;
}

// Utility function (replaces "fgets" where commenting is allowed)


// 7.5.0 paranoia, fptr changed to ffptr to be sure it doesn't change class
// fptr; this should have no effect as this is not a class member

char *fgets_skip_comments (char *buf, int count, FILE *ffptr)
{
    char *ret;
    for (;;)
    {
	ret = fgets (buf, count, ffptr);
	if (!ret) break;
	if (buf[0] == '\0' || buf[0] == '#' || buf[0] == '\n' ||
	    buf[0] == '\r')
	{
	    continue;
	}
	break;
    }
    return ret;
}

// Method definitions

TableFile::~TableFile ()
{
    int i;
    for (i = 0; i < field_count; i++)
    {
	if (_names) free (_names[i]);
	if (_short_names) free (_short_names[i]);
    }
    if (_names) free (_names);
    if (_short_names) free (_short_names);
    if (_widths) free (_widths);
    if (_filename) free (_filename);
    if (fptr) fclose (fptr);
}

void TableFile::set_position (int pos, const char **errmsg)
{
    if (_errmsg && !Strcmp (_errmsg, "EOF")) _errmsg = 0;
    if (0 != (*errmsg = _errmsg)) return;
    if (fseek (fptr, pos, 0))
    {
	*errmsg = _errmsg = "Error setting position in data file";
    }
    return;
}

PedsysFile::~PedsysFile () 
{
    if (code_recs)
    {
	for (int i = 0; i < field_count; i++)
	{
	    if (code_recs[i]) free (code_recs[i]);
	}
	free (code_recs);
    }
    if (user_setup) free (user_setup);
    if (_nominal_widths) free (_nominal_widths);
    
    code_recs = 0;
    user_setup = 0;
    _nominal_widths = 0;
}


TableFile *TableFile::open (const char *input_filename, const char **errmsg)
{
//  printf ("Tablefile opening %s\n",input_filename);
    *errmsg = 0;
    TableFile *ft = 0;
    if (strlen(input_filename) > 1024)
    {
	*errmsg = "Filename too long, greater than 1024 characters";
	return 0;
    }

// Try to open code file.
// If code file exists with name X.cde, we assume pedsys file
// (For this test, we do NOT remove last filename extension)

    char code_filename[1030];
    FILE *code_file = 0;

// If last char is ".", add "cde"
    if (input_filename[strlen(input_filename)-1] == '.')
    {
	strcpy (code_filename, input_filename);
	strcat (code_filename, "cde");
	code_file = fopen (code_filename, "r");
    }
    else
    {

// If last char is not ".", add ".cde"
	strcpy (code_filename, input_filename);
	strcat (code_filename, ".cde");
	code_file = fopen (code_filename, "r");
    }
    
// If that didn't work, try CDE

    if (!code_file)
    {
	if (input_filename[strlen(input_filename)-1] == '.')
	{
	    strcpy (code_filename, input_filename);
	    strcat (code_filename, "CDE");
	    code_file = fopen (code_filename, "r");
	}
	else
	{

// If last char is not ".", add ".cde"
	    strcpy (code_filename, input_filename);
	    strcat (code_filename, ".CDE");
	    code_file = fopen (code_filename, "r");
	}
    }
    
// If we didn't find unambiguous code file, try looking for commas in first rec

    if (!code_file)
    {

// Assumption: most speed critical data is comma delimited files
// So, rather than doing yet another test here, which would slow
// down comma delimited files, we try to make CommaDelimitedFile object
//   If that fails with no comma error, we do Pedsys last ditch test

	ft = new CommaDelimitedFile (input_filename);
	CommaDelimitedFile *cdf = (CommaDelimitedFile*) ft;
	if ((*errmsg = cdf->read_header ()))
	{
	    delete ft;
	    if (!strcmp (*errmsg, "File not Found"))
	    {
		return 0;
	    }
	    if (strcmp (*errmsg, "No commas"))
	    {
		return 0;
	    }
	}
	else
	{
	    return ft;
	}

// File didn't have commas or other distinguishing characteristics
// Make a last ditch effort to find code file by removing last extension
// (This is typical, though undesireable in our view, for pedsys files.)

// scan back to previous dot, remove postfix and add "cde"

	strcpy (code_filename, input_filename);
	if (strlen (code_filename) > 1)
	{
	    char *lastp = code_filename + strlen (code_filename) - 2;
	    while (lastp >= code_filename && *lastp != '.') lastp--;
	    if (*lastp == '.')
	    {
		*++lastp = '\0';
		strcat (code_filename, "cde");
		code_file = fopen (code_filename, "r");
	    }
	}
	if (!code_file)  // if not, try CDE
	{
	    strcpy (code_filename, input_filename);
	    if (strlen (code_filename) > 1)
	    {
		char *lastp = code_filename + strlen (code_filename) - 2;
		while (lastp >= code_filename && *lastp != '.') lastp--;
		if (*lastp == '.')
		{
		    *++lastp = '\0';
		    strcat (code_filename, "CDE");
		    code_file = fopen (code_filename, "r");
		}
	    }
	}
	if (!code_file)
	{
	    *errmsg = "Missing Pedsys code file or unsupported data format";
	    return 0;
	}
    }

// Code file found
    
    ft = new PedsysFile (input_filename);
    PedsysFile *pf = (PedsysFile*) ft;
    if ((*errmsg = pf->read_code_file (code_file)))
    {
	fclose (code_file);
	delete ft;  // Closes tfile
	return 0;
    }
    fclose (code_file);
    ft->_last_position = ftell (ft->fptr);
    return ft;
}

const char *PedsysFile::read_code_file (FILE *code_file)
{
// 7.5.0 paranoia change
    fptr = 0;
    fptr  = fopen (_filename, "r");
    if (!fptr)
    {
	return "File not Found";
    }

    code_recs = (CodeRecord**) Calloc (1, sizeof(CodeRecord*));
    code_recs[0] = 0;
    field_count = 0;

    char codebuf[64];
    char forline[1024];
    fgets (forline, 1024, code_file);
    while (fgets (codebuf, 64, code_file))
    {
#ifdef FISHER2
	if (!Strcmp (codebuf, "!EndCode\n")) return 0;  // Fisher2File
#endif
	if (strlen (codebuf) != 55)
	{
	    return "Invalid code file record length";
	}
	CodeRecord *cfr = (CodeRecord*) Calloc (1, sizeof (CodeRecord));
	strncpy ((char *) cfr, codebuf, 55);
	codebuf[2] = '\0';
	char junk[4];
	if (1 != sscanf (codebuf, "%d %s", &cfr->width, &junk))
	{
	    free (cfr);
	    return "Invalid code file record: width";
	}
	switch (cfr->data_type)
	{
	case ' ':
	    cfr->data_type = 'C';
	    break;
	case 'I':
	case 'C':
	case 'R':
	case 'D':
	    break;
	default:
	    free (cfr);
	    return "Invalid code file record: data type";
	}
	cfr->sep1 = cfr->sep2 = cfr->sep5 = cfr->sep6 = '\0';

// sep 3,4,5,6

	code_recs = (CodeRecord **) Realloc ((void*) code_recs, 
				     (1+(++field_count))*sizeof(CodeRecord*));
	code_recs[field_count-1] = cfr;
	code_recs[field_count] = 0;
    }
    
    return 0;
}


const char **TableFile::names (int *count, const char **errmsg)
{
     return TableFile::names_length (count, false, errmsg);
}


const char **TableFile::short_names (int *count, const char **errmsg)
{
    return TableFile::names_length (count, true, errmsg);
}

bool TableFile::test_name (const char *name, const char **errmsg)
{
    load_names (errmsg);
    if (*errmsg) return false;

    const char **test_names = (const char**) 
	((short_names_switch) ? _short_names : _names);

    int i;
    for (i = 0; i < field_count; i++)
    {
	if (!Strcmp (test_names[i], name)) return true;
    }
    return false;
}


const char **TableFile::names_length (int *count, bool short_length, const char **errmsg)
{
    short_names_switch = short_length;

    *count = 0;
    load_names (errmsg);
    if (*errmsg) return 0;

    *count = field_count;
    if (short_length)
    {
	return (const char**) _short_names;
    }
    else
    {
	return (const char**) _names;
    }
}


void CommaDelimitedFile::load_names (const char **errmsg) 
{*errmsg=_errmsg;}  // Loaded when opened

bool check_std_mnemonic (char *string)
{
    int i;
    const char **mntest = Standard_Pedsys_Mnemonics;
    char test[7];
    strncpy (test, string, 6);
    test[6] = '\0';

// Force user string into upper case

    for (i = 0; i < 6; i++)
    {
	test[i] = toupper (test[i]);
    }

// Compare with standard mnemonics

    for (i = 0; mntest[i]; i++)
    {
	if (!strncmp (test, mntest[i], 6))
	{
	    return true;
	}
    }
    return false;
}


void PedsysFile::load_names (const char **errmsg)
{
    if (0 != (*errmsg = _errmsg)) return;

    if (_names) return;
    _names = (char **) Calloc (field_count+1, sizeof (char*));
    _names[field_count] = 0;
    _short_names = (char **) Calloc (field_count+1, sizeof (char*));
    _short_names[field_count] = 0;

    for (int i=0; i < field_count; i++)
    {
	char mnemonics[22];
	
	code_recs[i]->sep5 = 0;
	strcpy (mnemonics, code_recs[i]->m1);
	if (check_std_mnemonic (&mnemonics[7]))
	{
	    mnemonics[7] = '\0';
	} 
	else if (check_std_mnemonic (&mnemonics[14]))
	{
	    mnemonics[14] = '\0';
	}
	else
	{
	    mnemonics[20] = '\0';
	}
	squeeze (mnemonics);    // Squeeze out whitespace
	mnemonics[18] = '\0';   // Truncate to 18 characters (limit)
	int len = strlen (mnemonics);
	int j;
#if 0
	for (j = 0; j < len; j++)
	{
	    if (mnemonics[j] == ',' || mnemonics[j] == '/')
	    {
		_errmsg = "Field name contains invalid char: , or /";
		break;
	    }
	}
#endif
	_names[i] = Strdup (mnemonics);

// Short name is made from first mnemonic only

	strcpy (mnemonics, code_recs[i]->m1);
	squeeze (mnemonics);
	_short_names[i] = Strdup (mnemonics);
    }
    if (_errmsg)
    {
	free (_names);
	_names = 0;
	free (_short_names);
	_short_names = 0;
	*errmsg = _errmsg;
    }
}


// This now returns "real" widths (as actually used, ignoring leading and
//   trailing spaces), not the code file widths (see nominal_widths()).

int *PedsysFile::widths (int *count, const char **errmsg)
{
    *count = 0;
    if (0 != (*errmsg = _errmsg)) return 0;

// If widths already found, return them

    *count = field_count;
    if (_widths) return _widths;

// Get nominal widths from code file

    int ccount;
    if (!_nominal_widths) nominal_widths (&ccount, errmsg);
    if (*errmsg) return 0;

// Save file position and rewind

    int save_position = ftell (fptr);
    if (save_position < 0)
    {
	*errmsg = _errmsg = "Save position error during field width scanning";
	*count  = 0;
	return 0;
    }
    PedsysFile:rewind (errmsg);
    if (*errmsg) return 0;

// Allocate widths

    _widths = (int*) Calloc (field_count+1, sizeof (int));
    _widths[field_count] = 0;

// Scan through each record in file

    char buf[MAX_REC_BYTES+1];
    while (fgets (buf, MAX_REC_BYTES, fptr))
    {
	int i;
	int cpos = 0;
	for (i = 0; i < field_count; i++)
	{
	    char save_char = buf[cpos+_nominal_widths[i]];
	    buf[cpos+_nominal_widths[i]] = '\0';
	    trim_blank_sides (&buf[cpos]);
	    int len = strlen (&buf[cpos]);
	    if (len > _widths[i]) _widths[i] = len;
	    buf[cpos+_nominal_widths[i]] = save_char;
	    cpos += _nominal_widths[i];
	}
    }

// Restore file position

    if (fseek (fptr, save_position, 0))
    {
	*errmsg = _errmsg = 
	    "Restore position error during field width scanning";
	free (_widths);
	_widths = 0;
	*count  = 0;
    }
    return _widths;
}

int *PedsysFile::nominal_widths (int *count, const char **errmsg)
{
    *count = 0;
    if (0 != (*errmsg = _errmsg)) return 0;

    *count = field_count;
    if (_nominal_widths) return _nominal_widths;

    _nominal_widths = (int*) Calloc (field_count+1, sizeof (int));
    _nominal_widths[field_count] = 0;

    for (int i=0; i < field_count; i++)
    {
	_nominal_widths[i] = code_recs[i]->width;
    }
    return _nominal_widths;
}


void PedsysFile::start_setup (const char **errmsg)
{
    if (0 != (*errmsg = _errmsg)) return;

    if (user_setup) free (user_setup);
    user_setup = 0;
    user_field_count = 0;
    return;
}

int PedsysFile::setup (const char *name, const char **errmsg)
{
    if (0 != (*errmsg = _errmsg)) return 0;

    int count;
    if (!_names) names (&count, errmsg);
    if (*errmsg) return 0;
    if (!_nominal_widths) nominal_widths (&count, errmsg);
    if (*errmsg) return 0;

    int running_position = 0;
    int found_width = 0;
    bool found = false;
    int found_position = 0;
    int found_index = 0;
    char **test_names = (short_names_switch) ? _short_names : _names;

    for (int i=0; i < field_count; i++)
    {
	if (!Strcmp (name, test_names[i]))
	{
	    if (found)
	    {
		sprintf (error_message, "Ambiguous field name: %s", name);
		*errmsg = _errmsg = error_message;
		return 0;
	    }
	    found = true;
	    found_width = _nominal_widths[i];
	    found_position = running_position;
	    found_index = i;
	}
	running_position += _nominal_widths[i];
    }
    if (!found)
    {
	sprintf (error_message, "Field not found: %s", name);
	*errmsg = _errmsg = error_message;
	return 0;
    }

    user_field_count++;
    if (!user_setup)
    {
	user_setup = (BeginWidth*) Calloc (user_field_count,
					   sizeof (BeginWidth));
    }
    else
    {
	user_setup = (BeginWidth*) Realloc (user_setup, user_field_count * 
					    sizeof (BeginWidth));
    }
    user_setup[user_field_count-1].begin = found_position;
    user_setup[user_field_count-1].width = found_width;
    return found_index;
}


char **PedsysFile::get (const char **errmsg)
{
    if (0 != (*errmsg = _errmsg)) return 0;

    char input_buffer[MAX_REC_BYTES];
    _last_position = ftell (fptr);
    if (!fgets (input_buffer, MAX_REC_BYTES, fptr))
    {
	*errmsg = _errmsg = "EOF";
	return 0;
    }
    int user_position = 0;
    if (user_field_count > MAX_REC_FIELDS)
    {
	*errmsg = _errmsg = "More than 40000 fields exceeds limit";
	return 0;
    }
    for (int i=0; i < user_field_count; i++)
    {
	int width = user_setup[i].width;
	if (user_position+width+1 > MAX_REC_BYTES)
	{
	    *errmsg = _errmsg = "More than 800000 bytes exceeds record limit";
	    return 0;
	}
	 bcopy (&input_buffer[user_setup[i].begin], 
		&user_buffer[user_position], width);
	 user_buffer[user_position+width] = '\0';
	 user_array[i] = &user_buffer[user_position];
	 trim_blank_sides (user_array[i]);
	 user_position += 1 + strlen (user_array[i]);
    }
    user_array[user_field_count] = 0;
    return user_array;
}


void PedsysFile::rewind (const char **errmsg)
{
    if (_errmsg && !Strcmp (_errmsg,"EOF")) _errmsg = 0;
    if (0 != (*errmsg = _errmsg)) return;

    errno = 0;
    ::rewind (fptr);
    if (errno)
    {
	*errmsg = _errmsg = "Unable to rewind data file";
    }
    _last_position = 0;
}

// in the current design, a CommaDelimitedFile must always
// invoke read_header during open, and that sets data
// starting position. Otherwise it is initialized to zero,
// which would work also, but it is never used that way.

void CommaDelimitedFile::rewind (const char **errmsg)
{
    if (_errmsg && !Strcmp (_errmsg,"EOF")) _errmsg = 0;
    if (0 != (*errmsg = _errmsg)) return;

    fseek (fptr, data_starting_position, SEEK_SET);
}

const char *CommaDelimitedFile::read_header ()
{
    fptr  = fopen (_filename, "r");
    if (!fptr)
    {
	return "File not Found";
    }

// Get user bufsize

    int user_bufsize = Option::get_int ("CsvBufSize");

// Initial allocation of _names and _short_names

    field_count = 0;
    _names = (char**) Calloc (1, sizeof (char*));
    _short_names = (char **) Calloc (1, sizeof (char*));
    _types = 0;

// Process header one buffer at a time to allow for infinite length
//   if more characters remain at end of record, copy to buf beginning and
//   read next record with an offset, which is then cleared.

    char* ret;
    int offset = 0;
    bool first = true;
    char* pointer;
    char* end_pointer;
    char echar; // end character

    if (user_bufsize) suggested_bufsize = user_bufsize;
    int header_bufsize = suggested_bufsize;

    char buf[header_bufsize];
    int buffers_used = 0;

    for (;;)
    {
	bool got_a_field = false;
	ret = fgets (&buf[offset], header_bufsize-offset, fptr);
	if (!ret)
	{
	    return "Unable to read comma delimited header";
	}
	offset = 0;
	pointer = buf;
	if (buffers_used++ > 200000.0 / header_bufsize)
	{
	    if (!user_bufsize) suggested_bufsize = 200000;
	}

// Check for commas (first buffer only)

	if (first)
	{
	    pointer = buf;
	    int ccount = 0;
	    int tested = 0;
	    while (*pointer != '\0')
	    {
		if (*pointer==',') ccount++;
		pointer++;
		tested++;
	    }
	    if (ccount < 1)
	    {
		return "No commas";
	    }
	    first = false;
	}

// Process header

	pointer = buf;
	bool is_type = false;
	bool next_type = false;
	while (1)
	{

// get next name up to terminator, ignoring type

	    end_pointer = pointer;
	    while (1)
	    {
		if (*end_pointer != ','  && *end_pointer != '\0' && 
		       *end_pointer != '\n' &&
		       *end_pointer != '\r')
		{
		    end_pointer++;
		}
		else
		{
		    if (*end_pointer == ':')  // we forcibly i
		    {
			*end_pointer = 0;
			end_pointer++;
		    }
		    else
		    {
			break;
		    }
		}
	    }
	    echar = *end_pointer;
	    *end_pointer = '\0';

// now make sure name is complete, that buffer ending didn't truncate name
//   truncation occurred if the last character is null, because if record
//   terminates normally, a newline comes first and will be detected above.

// on the other hand, if we ended with \n or \r, this is the last name, but
// instead we test for comma at the bottom of the loop and return if not comma

	    if (echar == '\0')
	    {
		if (!got_a_field)
		{
		    return "One data field exceeds bufsize in csv file";
		}

		int i;
		for (i = 0; pointer[i] != '\0'; i++)
		{
		    buf[i] = *(pointer+i);

		}
		buf[i] = '\0';
		offset = i;
		break;
	    }

// found end of name, save to arrays of names

	    got_a_field = true;
	    *end_pointer = '\0';
	    _names = (char**) Realloc ((void*) _names, (1+(++field_count)) *
				   sizeof (char*));
	    _short_names = (char**) Realloc ((void*) _short_names, 
					 (1+field_count) * sizeof (char*));
	    _types = (short*) Realloc ((void*) _types,
					 (1+field_count) * sizeof (short*));

	    trim_blank_sides (pointer);

            char* colonpos = 0;
	    if (0 != (colonpos = strchr (pointer, ':')))
	    {
//		printf ("found colon in field %d\n",field_count-1);
		*colonpos = 0;
		char* typepos = colonpos+1;
		if (!Strcmp (typepos,"nifti"))
		{
//		    printf ("Setting field %d = 11\n",field_count-1);
		    _types[field_count-1] = 11;
		}
	    }
	    else
	    {
		_types[field_count-1] = 0;
	    }

	    _names[field_count-1] = Strdup (pointer);
	    if (strlen(pointer) > SHORT_NAME_LENGTH)
	    {
		pointer[SHORT_NAME_LENGTH] = '\0';
	    }
	    _short_names[field_count-1] = Strdup (pointer);


// If this name not terminated by comma, this is the end, return here
//   But check if this line is terminated by \r, since fgets does not
//   support that old Mac format (still used by many Mac programs!)
//   but if /n follows, that's MSDOS format, and it's OK

	    if (echar != ',')
	    {
		if (echar == '\r') {
		    if (*(end_pointer+1) != '\n') {
			fprintf (stderr, 
    "File %s has unsupported text line terminators\n",_filename);
			return (
    "Use retext command to fix file before using");
		    }
		}
		_names[field_count] = 0;
		_short_names[field_count] = 0;
		data_starting_position = ftell (fptr);
		return 0;
	    }
	    pointer = end_pointer+1;

	} // continue parsing record

    } // continue reading into buffer

// above loop has no exit can't get here

}


int *CommaDelimitedFile::widths (int *count, const char **errmsg)
{
    if (0 != (*errmsg = _errmsg)) return 0;

// If widths already found, return them

    *count = field_count;
    if (_widths)
    {
	return _widths;
    }

// Save file position

    int save_position = ftell (fptr);
    if (save_position < 0)
    {
	*errmsg = _errmsg = "Save position error during field width scanning";
	*count  = 0;
	return 0;
    }

    CommaDelimitedFile::rewind (errmsg);
    if (*errmsg) return 0;
    _widths = (int *) Calloc ((1+field_count), sizeof (int));

    getwidths = true;

    while (!_errmsg)
    {
	CommaDelimitedFile::get (errmsg);
    }

    getwidths = false;

    if (_errmsg && strcmp (_errmsg, "EOF"))
    {
	return 0;
    }
    *errmsg = _errmsg = 0;

// Restore file position

    if (fseek (fptr, save_position, 0))
    {
	*errmsg = _errmsg = 
	    "Restore position error during field width scanning";
	free (_widths);
	_widths = 0;
	*count  = 0;
    }

    return _widths;
}


    
void CommaDelimitedFile::start_setup (const char **errmsg)
{
    if (0 != (*errmsg = _errmsg)) return; 
    
    if (user_indexes) free (user_indexes);
    user_indexes = (int *) Calloc (1, sizeof (int));
    user_field_count = 0;
    highest_user_index = -1;
}

int CommaDelimitedFile::setup (const char *name, const char **errmsg)
{
    if (0 != (*errmsg = _errmsg)) return 0; 

    for (int i = 0; i < field_count; i++)
    {
	char **test_names = (short_names_switch) ? _short_names : _names;

	if (!Strcmp (name, test_names[i]))
	{
	    user_field_count++;
	    user_indexes = (int *) Realloc (user_indexes, 
					    user_field_count*sizeof(int));
	    user_indexes[user_field_count-1] = i;
	    user_array = (char **) Realloc (user_array,
				       (1+user_field_count)*sizeof(char*));
	    if (i > highest_user_index) 
	    {
		highest_user_index = i;
	    }
	    return i;
	}
    }
    sprintf (error_message, "Name not found: %s", name);
    *errmsg = _errmsg = error_message;
    return 0;
}



#define DEBUGGET 1

char **CommaDelimitedFile::get (const char **errmsg)
{
    if (0 != (*errmsg = (const char*) _errmsg)) return 0; 
    _last_position = ftell (fptr);


// This function may also determine field widths if called by the widths
// method and the getwidths flag is set.  In that case, it does not actually
// return data, instead it reads EVERY field and fills the _widths
// for that field if the current number is lower than the current width

    int current_record_buffer = 0;
    int offset = 0;
    int iscan = 0;
    int icopy = 0;
    int copypos = 0;
    bool more_to_read = false;
    int firstindex = 0;
    int lastindex = 0;
    int scanp = 0;
    int last_incomplete = 0;
    char testchar = '\0';
    bool reuse_buffer = true;
    char* buf = 0;
    int width_count = 0;

// record_buffers is an array of buffer pointers to buffers into which
// data is read by fgets.  Only the buffers actually needed are retained,
// the others are overwritten on the next buf read if not needed. In the end,
// the buffers in record_buffers must be retained until the next call to get.
// If not deleted here, the record_buffers and the buffers it it points to
// are deleted when the object is deleted.

// table_pointers point to the beginning of every field, needed or not
//  initial scanning finds table pointers but only to highest index
//  table_pointers are freed here normally, but if not then by
//  tablefile object when it is is deleted

// user_indexes points to the particular fields requested by caller and was
// previously set up by the setup function.

// user_array is the array of char* pointers returned to the caller.  It is
// also allocated by the setup function for efficiency (one allocation for
// many uses).

    int highest_index = highest_user_index;
    if (getwidths)
    {
	highest_index = field_count-1;
    }

    for (;;)
    {
	if (!record_buffers)
	{
//	    printf ("allocating new record buffers\n");
	    record_buffer_count = current_record_buffer+1;
	    record_buffers = (char**) malloc (record_buffer_count+2 
						  * sizeof (char*));
	    buf = (char*) malloc (suggested_bufsize);
	    record_buffers[current_record_buffer] = buf;
	}
	else if (!reuse_buffer)
	{
	    current_record_buffer++;
	    if (current_record_buffer+1 > record_buffer_count)
	    {
//		printf ("allocating new buffer\n");
		record_buffer_count = current_record_buffer+1;
		record_buffers = (char**) realloc (record_buffers, 
				       (record_buffer_count+2)*sizeof(char*));
		buf = (char*) malloc (suggested_bufsize);
		record_buffers[current_record_buffer] = buf;
	    }
	    else
	    {
//		printf ("re-using buf[%d]\n",current_record_buffer);
		buf = record_buffers[current_record_buffer];
	    }
	}
	else
	{
//	    printf ("simply reusing buf[%d]\n",current_record_buffer);
	    buf = record_buffers[current_record_buffer];
	}
	reuse_buffer = true;

	if (offset)
	{
	    strcpy (buf,prestring);
	}

	if (!(fgets_skip_comments (&buf[offset], suggested_bufsize-offset, 
				   fptr)))
	{

// We expected to read something (either in first buffer or later)
// but read nothing because of EOF (likely) or some error (unlikely).
// Previously tablefile handled all such as EOF, not fully safe

	    *errmsg = _errmsg = "EOF";
	    if (!feof(fptr))
	    {
		printf ("ERRNO reading file is %d\n",errno);
		*errmsg = _errmsg = "Error during file read";
	    }
	    return 0;
	}

//	printf ("Read buffer in commadelimitedfile.get: %s\n",buf);

	offset = 0;
	scanp = 0;

// Set pointers to the records in this buffer
//   iscan is the index (zero based) of the buffer field we are working on
//   firstindex is index of first field in buffer
//   lastindex is index of last full field in buffer
//   userpointer[iscan] is indexed pointer to fields


//	printf ("Scanning from %d to %d\n",iscan,highest_index);
	lastindex = firstindex = iscan;
	for (; iscan <= highest_index; iscan++)
	{
//	    printf ("Scanning index %d\n",iscan);
	    if (!table_pointers || table_pointer_count < iscan+1)
	    {
		table_pointer_count = iscan+1;
		table_pointers = (int*) Realloc (table_pointers,
						   (table_pointer_count+1) * 
						   sizeof(char*));
	    }
	    last_incomplete = table_pointers[iscan] = scanp;
	    lastindex = iscan;
	    testchar = buf[scanp++];
	    while (testchar != ',' && testchar != '\0' && testchar != '\n'
		   && testchar != '\r')
	    {
		testchar = buf[scanp++];
	    }
	    
// null terminate this datum

	    if (testchar != '\0') buf[scanp-1] = '\0';

// If this datum ends with null, it is incomplete, because at the true EOL
// a newline comes first and would break the scanning loop above.

	    if (testchar == '\0')
	    {
		lastindex--;
		break;
	    }

// if this ends with linefeed, it is end of record so break

	    if (testchar == '\n' || testchar == '\r')
	    {
		break;
	    }
	}
    
// If scanning ended by finding all requested fields, we haven't actually
// tested whether there are more buffers for this record which need to be read
// and discarded.  So we check for that now
// If testchar is null, we need to continue reading
// If testchar is comma, we need to scan to end of current record
// If testchar is newline, we are done

	if (testchar == '\n' || testchar == '\r')
	{
	    more_to_read = false;
	    if (iscan < highest_index)
	    {
		char showbuf[40];
		fseek (fptr,_last_position,0);
		fgets (showbuf,40,fptr);
		fprintf (stderr,"Short record: %s ...\n",showbuf);
		_errmsg = *errmsg = "Short record in input file";
		return 0;
	    }
	}
	else if (testchar == '\0')
	{
	    more_to_read = true;
	}
	else
	{
	    more_to_read = true;
	    char tchar;
	    int tscanp = scanp;
	    while (0 != (tchar = buf[tscanp++]))
	    {
		if (tchar == '\n' || tchar == '\r')
		{
		    more_to_read = false;
		    break;
		}
	    }
	    if (more_to_read)
	    {
//		printf ("Didn't find EOL in post scan so must read more.\n");
	    }
	}

// Now copy these data fields
//   only if complete within this segment

	int last_field = user_field_count-1;
	int first_field = 0;
	if (getwidths)
	{
	    first_field = firstindex;
	    last_field = lastindex;
	}
//	printf ("last_field is %d\n",last_field);
//	printf ("first_field is %d\n",first_field);
	for (icopy = first_field; icopy <= last_field; icopy++)
	{
//	    printf ("icopy is %d\n",icopy);
// user_indexes point to the fields requested for a get
// user_indexes is not needed or available for getwidths, which
//   gets the widths of all fields

	    int needindex;
	    if (!getwidths)
	    {
		needindex = user_indexes[icopy];
	    }
	    else
	    {
		needindex = icopy;
	    }
	    if (needindex >= firstindex &&
		needindex <= lastindex)
	    {
		if (!getwidths) reuse_buffer = false;

// null terminated data is no longer copied (big change)!
// Now we simply point to null terminated data within the initial buffer,
// which is retained until next record or file is closed anyway.
//

//		printf ("pointer is %d\n",table_pointers[needindex]);
//		printf ("Copying >%s<\n",&buf[table_pointers[needindex]]);

		if (!getwidths)
		{
		    short data_type = _types[needindex];

		    if (data_type < 11)
		    {

// ordinary numeric data...just copy pointer to datum

			user_array[icopy] = &buf[table_pointers[needindex]];
			trim_blank_sides (&buf[table_pointers[needindex]]);
		    }
		    else
		    {

// This field is not raw numeric data.
// Instead, it is the filename of a binary file.
// We obtain the data by employing the required method on the file
//   Currently, RicVolumeSet is the only supported filetype so we don't check
//   Otherwise, the binary filetype would be specified by the field type

// We need to pull out arguments from the data string here:
//    the NIFTI filename, and the Volume number, like this:
//
//    images.gz:5
//
// These are combined with the current voxel to get the actual datum.

#ifndef RICVOLUMESET
			_errmsg = *errmsg = "RicVolumeSet not yet supported on this system";
			return 0;
#else

			char* data_string = Strdup 
			    (&buf[table_pointers[needindex]]);
			
			int volume_number = 0;
			char* valstring;
			char* colonpos;

			char* scanpos = data_string;
			if (0 != (colonpos = strchr (scanpos,':')))
			{
			    *colonpos = '\0';
			    scanpos = valstring = colonpos+1;
			}
			else
			{
			    printf ("Binary filename but no volume spec: %s\n"
				    ,data_string);
			    _errmsg = *errmsg= "Missing volume spec";
			    free (data_string);
                            return 0;
			}
			char* endptr;
			errno = 0;
			long larg = strtol (valstring, &endptr, 0);
			if (errno || larg > INT_MAX || larg < INT_MIN)
			{
			    *colonpos = ':';
			    printf (
				"Invalid volume spec %s in file %s\n",
					valstring, _filename);
			    _errmsg = *errmsg = 
				"Invalid binary volume specified";
			    free (data_string);
			    return 0;
			}
			volume_number = (int) larg;
			int volume_adjust = Option::get_int ("RicVolOffset");
			volume_number -= volume_adjust;

// Get the current voxel value

			if (!Voxel::Valid())
			{
			    printf ("Voxel not defined.\n");
			    _errmsg = *errmsg =
				"Voxel not defined";
			    free (data_string);
			    return 0;
			}

// Get value from RicVolumeSet

			if (0==GlobalBinFilename ||
			    strcmp(GlobalBinFilename,data_string))
			{
			    if (GlobalBinFilename)
			    {
				free (GlobalBinFilename);
				GlobalBinFilename = 0;
				delete GlobalBin;
				GlobalBin = 0;
			    }
			    printf ("Opening data RicVolumeSet %s...\n",
				data_string);
			    try
			    {
 				GlobalBin = new RicVolumeSet (data_string);
			    }
			    catch (...)
			    {
				printf ("Unable to read binary file\n");
				_errmsg = *errmsg = "Unknown binary file type";
				free (data_string);
				return 0;
			    }
			    printf ("Got new RicVolumeSet: %ld\n",(long) GlobalBin);
			    printf ("This image has %d volumes\n",
				    GlobalBin->nvol);

			    GlobalBinFilename = Strdup (data_string);
			}
			free (data_string);
			float val = 1e-10;

			if (GlobalBin->nvol <= volume_number)
			{
			    printf ("Invalid volume number %d\n",volume_number);
			    _errmsg = *errmsg = "Invalid volume number";
			    return 0;
			}
			try 
			{
//			    printf ("pointer is %ld\n", (long) &GlobalBin->VolSet[volume_number]);
			    int x = Voxel::X;
			    int y = Voxel::Y;
			    int z = Voxel::Z;
			    val = 
				GlobalBin->VolSet[volume_number].vox[x][y][z];
			}
			catch (...)
			{
			    printf ("Got to catch statement\n");
			    return 0;
			}

			char fbuf[64];
			sprintf (fbuf,"%14.8g", (double) val);
			
			char* sptr = Strdup (fbuf);
			user_array[icopy] = sptr;
			freelist = (char**) Realloc (
			    freelist,(freelistcnt+2)*sizeof(char**));
			freelist[freelistcnt++] = sptr;
#endif
		    }
		}

		if (getwidths)
		{
		    trim_blank_sides (&buf[table_pointers[needindex]]);
		    int newwidth = strlen(&buf[table_pointers[needindex]]);
//		    printf ("width[%d] = %d\n",width_count,newwidth);
		    if (newwidth > _widths[width_count++])
		    {
			_widths[width_count-1] = newwidth;
		    }
		}
	    }
	}

// If no more buffers, exit here

	if (lastindex >= highest_index)
	{
//	    printf ("last_index, highest_index: %d, %d\n",lastindex,highest_index);
	    if (!getwidths) user_array[user_field_count] = 0;

	    if (more_to_read)
	    {
//		printf ("Reading to end of buffer\n");
		char eolbuf[suggested_bufsize];
		bool end_of_record = false;
		for (;;)
		{
		    if (!(fgets(eolbuf,suggested_bufsize,fptr)))
		    {
			if (feof(fptr))
			{
			    end_of_record = true;
			    break;
			}
			*errmsg = _errmsg = "Error during file read";
			return 0;
		    }
		    for (int j = 0; eolbuf[j] != '\0'; j++)
		    {
			if (eolbuf[j] == '\r' || eolbuf[j] == '\n')
			{
			    end_of_record = true;
			    break;
			}
		    }
		    if (end_of_record)
		    {
			break;
		    }
		}
	    }
	    else
	    {
//		printf ("no more to read\n");
	    }
	    if (!getwidths)
	    {
		return user_array;
	    }
	    else
	    {
		return 0; // No return needed for getwidths
	    }
	}

// Must read another buffer...
// Copy start of last datum to beginning of buffer and set offset

	scanp = last_incomplete;
	offset = strlen (&buf[scanp]);
	if (prestring) free (prestring);
	prestring = (char*) malloc (offset+1);
	strcpy (prestring,&buf[scanp]);
//	printf ("Copied >%s<\n",prestring);


    } // loop multiple buffer reads if needed

}

CommaDelimitedFile::~CommaDelimitedFile()
{
    int i;
    for (i = 0; i < freelistcnt; i++)
    {
	free (freelist[i]);
    }
    if (freelist) free (freelist);
    freelist = 0;
	
    if (user_indexes) free (user_indexes);
    user_indexes = 0;

    if (user_array) free (user_array);
    user_array = 0;

    if (record_buffers)
    {
	for (i = 0; i < record_buffer_count; i++)
	{
	    if (record_buffers[i]) free (record_buffers[i]);
	}
	free (record_buffers);
	record_buffers = 0;
	record_buffer_count = 0;
    }

// table_pointers is only used internally in "get"
// however, because of all the 0 returns from "get", it's hard to free it
// there

    if (table_pointers)
    {
	free (table_pointers);
    }
    table_pointers = 0;

    if (prestring) free (prestring);
    prestring = 0;
    if (_types) free (_types);
}


void squeeze (char *s)  // Squeeze out whitespace
{
    char *s2 = s;
    while (*s2 != '\0')
    {
	if (*s2 != ' ' && *s2 != '\t')
	{
	    *s++ = *s2++;
	}
	else
	{
	    *s2++;
	}
    }
    *s = '\0';
}

