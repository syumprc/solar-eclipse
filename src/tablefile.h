#ifndef TABLEFILE_H
#define TABLEFILE_H

#define MAX_REC_BYTES 800000
#define MAX_REC_FIELDS 40000
#define SHORT_NAME_LENGTH 6


/*
 * Tablefile: A library for accessing data files ("tables") which may
 *   physically exist in different formats.  The library figures out
 *   the correct format automatically.
 *
 * Currently supported formats are PEDSYS and Comma Delimited.
 * Comma Delimited files must have first record having field names.
 * (Fisher2 format files are also supported, though inefficiently.)
 *
 * Charles Peterson, 20 May 1998.
 * Copyright (c) 1998, Southwest Foundation for Biomedical Research
 */


/*
 * C Interface: Be sure to read notes which follow.
 *              An example is also shown below.
 *
 * (The C interface may also, of course, be used in C++.
 *  In fact, this is currently recommended because the C++ interface is
 *  subject to future refinement since it has more potential.)
 */

#ifdef __cplusplus
extern "C" 
{
#endif
    typedef int status;
    typedef void tab;    /* table object...do not modify this from C! */

/* Open table file for reading */
    tab *tab_open (const char *filename, const char **errmsg);  

/* Get names (and set "short status" if calling tab_short_names) */
    const char **tab_names (tab *tabp, int *count, const char **errmsg);
    const char **tab_short_names (tab *tabp, int *count, const char **errmsg);

/* Test for existance of a name */
    int tab_test_name (tab *tabp, const char *name, const char **errmsg);

/* Get data widths */
    int *tab_widths (tab *tabp, int *count, const char **errmsg);

/* Setup desired data record...call tab_setup once for each field desired */
    void tab_start_setup (tab *tabp, const char **errmsg);
    int tab_setup (tab *tabp, const char *name, const char **errmsg);
/* tab_setup returns index into withds and names records */

/* Get data in user record */
    char **tab_get (tab *tabp, const char **errmsg);

/* Rewind and Close */
    void tab_rewind (tab *tabp, const char **errmsg);
    void tab_close (tab *tabp);

/* Remind me what filename this is */
    const char *tab_filename (tab *tabp);

/* Save current file position; move to previous file position */
//    int tab_get_position (tab *tabp);
    void tab_set_position (tab *tabp, int position, const char **errmsg);

#ifdef __cplusplus
};
#endif

/*
 * Notes:
 * (1) All errmsg returns should be tested.  (NULL means no error.)  The
 *     string "EOF" is returned when tab_get reaches EOF.  Any other string
 *     indicates a serious error.  EOF is cleared by rewind.
 *
 *     It is particularly important to test tab_setup returns even if you
 *     are using verified names.  Redundant names are detected at this
 *     point.
 *
 * (2) Although not necessarily recommended,
 *     in some cases (e.g. multiple tab_setup's) you may do a
 *     series of operations and then test the last errmsg.*  Once an error
 *     occurs, error state is latched and no further operations are done.
 *     However, beware that if an error has occurred, functions returning
 *     arrays will return NULL, and if you attempt to use the value you
 *     will SEGV.  For functions that return an object or array, you can
 *     test for 0 return instead of errmsg.
 *
 * (3) The user must copy all data returned in arrays (including names and
 *     widths.  All storage is re-used by the table object on the next call, 
 *     and then freed when the table object is closed.
 *
 * (4) Requesting short names will set short name status.  That will affect
 *     how names are tested in the test or setup functions.
 *
 * (5) Name comparisons are done with case insensitivity.
 *
 * Example:
 *
 * #include "tablefile.h"
 * 
 * void errout (char *message)
 * {
 *     fprintf (stderr, "%s\n", message);
 *     exit (1);
 * }
 * 
 * int main (int argc, char **argv)
 * {
 *     char *errmsg;
 *     tab *table;
 *     int count;
 *     char **names;
 *     int *widths;
 *     char **data;
 *     int i;
 * 
 *     table = tab_open (argv[1], &errmsg);
 *     if (errmsg) errout (errmsg);
 * 
 *     names = tab_names (table, &count, &errmsg);
 *     if (errmsg) errout (errmsg);
 * 
 *     widths = tab_widths (table, &count, &errmsg);
 *     if (errmsg) errout (errmsg);
 * 
 *     tab_start_setup (table, &errmsg);
 *     if (errmsg) errout (errmsg);
 * 
 *     for (i = 2; i < argc; i++)
 *     {
 *         tab_setup (table, argv[i], &errmsg);
 *         if (errmsg) errout (errmsg);
 *     }
 * 
 *     while ((data = tab_get (table, &errmsg)))
 *     {
 *         for (i = 2; i < argc; i++)
 * 	   {
 * 	       if (i == 2)
 * 	       {
 * 		   printf ("%s", data[i-2]);
 * 	       }
 * 	       else
 * 	       {
 * 		   printf (",%s", data[i-2]);
 * 	       }
 * 	   }
 * 	   printf ("\n");
 *     }
 *     if (errmsg && !strcmp (errmsg, "EOF")) errmsg = 0;
 *     if (errmsg) errout (errmsg);
 *     return 0;
 * }
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#ifdef USE_SAFELIB
#include "safelib.h"
#else
#define Malloc malloc
#define Calloc calloc
#define Realloc realloc
#define Strcmp strcmp
#endif

/* 
 * C++ interface
 */

#ifdef __cplusplus

// bool type not defined in some compilers
// but it's a reserved word (not a definition) in ANSI C++ compilers
// Makefile must define NEEDS_BOOL for pre-ANSI compilers

#ifdef NEEDS_BOOL
enum bool {false, true};
#undef NEEDS_BOOL
#endif

void squeeze (char *s);  // Squeeze out whitespace

class TableFile
{
protected:
    const char *_errmsg;
    FILE *fptr;
    int field_count;

    int user_field_count;
    bool short_names_switch;
    char **_names;
    char **_short_names;
    int *_widths;

    char error_message[256];
    virtual void load_names(const char **errmsg) = 0;
    const char **names_length (int *count, bool short_length, 
			       const char **errmsg);
    long  _last_position;
public:
    char *_filename;
    TableFile () {_errmsg=0; fptr = 0; field_count=0; 
                            user_field_count=0; 
                            short_names_switch=false;
                            _names=0; _short_names=0; _widths=0; _filename=0;
                            }
    virtual ~TableFile();
    static TableFile *open (const char *filename, const char **errmsg);

    const char **names (int *count, const char **errmsg);
    const char **short_names (int *count, const char **errmsg);
    bool test_name (const char *name, const char **errmsg);
    virtual int *widths (int *count, const char **errmsg) = 0;

    virtual void start_setup (const char **errmsg) = 0;
    virtual int setup (const char *name, const char **errmsg) = 0; // ret index
    virtual char **get (const char **errmsg) = 0;
    virtual void rewind (const char **errmsg) = 0;
    const char *filename () {return (const char*) _filename;}
    long get_position () {return _last_position;}
    virtual void set_position (int pos, const char **errmsg);
};

void inline trim_blank_sides (char *buffer)
{
    int j;
    int k;
    bool blank = true;
    for (j=0; '\0' != buffer[j]; j++)
    {
	if (buffer[j] != ' ' && buffer[j] != '\t')
	{
	    blank = false;
	    int last_non_blank=j;
	    if (j > 0)
	    {
	    // Trim leading space
		for (k=0; ;k++,j++)
		{
		    char copied = buffer[k] = buffer[j];
		    if (copied == '\0') break;
		    if (copied != ' ' && copied != '\t') last_non_blank = k;
		}
	    }
	    else
	    {
	    // No leading space, but find last_non_blank
		for (k=0; '\0' != buffer[k]; k++)
		{
		    if (buffer[k] != ' ' && buffer[k] != '\t') 
			last_non_blank = k;
		}
	    }
	    // Trim following space
	    buffer[last_non_blank+1] = '\0';
	    break;
	}
    }
    if (blank)
    {
	buffer[0] = '\0';
    }
}



#endif
#endif
