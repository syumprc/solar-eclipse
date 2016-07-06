/* safelib.c
 * safe memory allocation and other functions
 *   (returns usually don't have to be checked; they are checked here)
 *
 * comments: Requires that there be a function named 'error' to which a
 *   string is passed.  Also writes directly to standard error stream (2)
 *   to avoid additional dynamic memory usage if memory is low.
 *
 */

#include <errno.h>
#include <termios.h>
#include <sys/ioctl.h>
#include "safelib.h"

void error (char *error_message); /* defined by client program */

char *MMessage = "\nFatal error: Not enough memory.\n";
char *FMessage = "Warning.  File didn't close.\n";


void *smalloc (size_t size)
{
    void *mem = malloc (size);
    if (!mem)
    {
	write (2, MMessage, strlen (MMessage)); /* requires little memory */
	exit (EXIT_FAILURE);
    }
    return mem;
}

void *scalloc (size_t nelem, size_t elsize)
{
    void *mem = calloc (nelem, elsize);
    if (!mem)
    {
	write (2, MMessage, strlen (MMessage)); /* requires little memory */
	exit (EXIT_FAILURE);
    }
    return mem;
}

void *srealloc (void *ptr, size_t size)
{
    void *mem = realloc (ptr, size);
    if (!mem)
    {
	write (2, MMessage, strlen (MMessage)); /* requires little memory */
	exit (EXIT_FAILURE);
    }
    return mem;
}

char *sstrdup (const char *str)
{
    char *newstr = NULL;

    if (str)
    {
	newstr = smalloc (1 + strlen (str));
	strcpy (newstr, str);
    }
    return newstr;
}

void sclose (int fdes)
{
    if (close (fdes))
    {
	write (2, FMessage, strlen (FMessage));
    }
}

void sfclose (FILE *fptr)
{
    if (fclose (fptr))
    {
	write (2, FMessage, strlen (FMessage));
    }
}

void ssigset (int sig, void (*disp)(int))
{
    if (SIG_ERR == sigset (sig, disp))
    {
	error ("Error setting up signal handler!");
    }
}

void ssighold (int sig)
{
    if (sighold (sig))
    {
	perror ("Error holding signal!");
    }
}

void ssigrelse (int sig)
{
    if (sigrelse (sig))
    {
	perror ("Error holding signal!");
    }
}

/*
 * Fairly safe dynamic string handling
 * stringpp is a pointer pointer to a dynamically allocated string (ONLY!!!)
 * either *stringpp or append (or both) may be NULL
 * if either/both is/are non-null, either/both MUST be null terminated
 * string is resized (or created) to fit current data
 * if no data, a 'zero' length null-terminated string will be created anyway
 */

char *string__append (char **stringpp, const char *appendum)
{
    int oldlen = (*stringpp) ? strlen (*stringpp) : 0;
    int applen = (appendum) ? strlen (appendum) : 0;

    *stringpp = srealloc (*stringpp, oldlen+applen+1);
    if (!oldlen) **stringpp = '\0';

    if (applen)
    {
	strcat (*stringpp, appendum);
    }
    return *stringpp;
}

/*
 * spadf pads out a c string with blanks to the specified length for use 
 * by FORTRAN.  FORTRAN strings are not terminated...they are blank padded.
 * Note: the length had better be correct or a buss error will occur!
 * The string should be null terminated (as for a C string) on input, but
 * will probably not be null terminated on output.
 */
void spadf (char *string, int length)
{
    int slen = strlen (string);
    while (slen < length)
    {
	string[slen++] = ' ';
    }
}

/*
 * The main virtue of sfgets is that it
 * removes the trailing \n from each line, which is a pain in the neck
 * to have to do each time.
 */
char *sfgets (char *buffer, size_t bufsiz, FILE *file)
{
    char *rchar = fgets (buffer, bufsiz, file);
    int len;

    if (rchar)
    {
	len = strlen (buffer);
	if (buffer[len-1] == '\n')
	{
	    buffer[len-1] = '\0';
	}
    }
    else
    {
	buffer[0] = '\0';
    }
    return rchar;
}

/*
 * The main virtue of sgets is that it checks the return for input error.
 * Generally, input error does not occur, except in wierd cases
 * (such as running in background).  If input error occurs, the best
 * thing may be to print a message and exit, otherwise we might get
 * caught in an infinite loop.
 */
char *sgets (char *inbuf)
{
    if (!gets (inbuf))
    {
	if (errno != EINTR)
	{
	    if (errno == 5)
	    {
		errno = 0;
		error ("Background job is unable to read from terminal");
	    }
	    else
	    {
		error ("Terminal input error");
	    }
	    exit (EXIT_FAILURE);
	}
    }
    return inbuf;
}

char *scan_r (char *ibuf, char *delim, char **next, int get_remainder)
/*
 * This should be called through either the NextToken or
 * RestOfTokens interfaces, which set the 'get_remainder' switch
 * appropriately.
 *
 * This replaces strtok_r, but doesn't mess up the original
 * buffer.  Instead, it returns a copy of the next token.
 * You can free it when you are done, unless NULL is returned.
 *
 * Null is returned if there is no token (nothing or only delimiters left).
 *
 */
{
    char savech;
    char *token = NULL;
    size_t tbegin, tlen;
    size_t len = strlen (ibuf);

    char *buf = ibuf;         /* Save in case it's same var as next */
    if (next) *next = NULL;

    tbegin = strspn (buf, delim);
    if (tbegin == len) return NULL;

    if (!get_remainder)
    {
	tlen = strcspn (&buf[tbegin], delim);
    }
    else
    {
    /*
     * Remove trailing whitespace instead, then return rest of line
     */
	int tend = len-1;
	while (tend > tbegin)
	{
	    if (buf[tend] == ' ' || buf[tend] == '\t')
	    {
		tend--;
	    }
	    else
	    {
		break;
	    }
	}
	tlen = 1 + tend - tbegin;
    }
    if (tlen == 0) return NULL;
    savech = buf[tbegin+tlen];
    buf[tbegin+tlen] = '\0';
    token = Strdup (&buf[tbegin]);
    buf[tbegin+tlen] = savech;

    if (next) *next = &buf[tbegin+tlen];
    return token;
}

char *NextToken (char *buf, char *delim, char **next)
{
    return scan_r (buf, delim, next, 0);
}

char *RestOfTokens (char *buf, char *delim, char **next)
{
    return scan_r (buf, delim, next, 1);
}


char *NextArg (char *string, char delimiter, char **lastp)
{
/* This is a scanner not entirely unlike strtok... or strchr.
 *
 * Intended for argument lists like (1,,2)
 *
 * See NextToken for a more general scanning interface
 *
 * But, Instead of passing over multiple delimiters, it considers that
 * an adjacent pair of delimiters marks a 'NULL' argument.
 * If the first character is a delimiter, the first argument is a NULL
 * argument.  NextArg returns NULL for each NULL argument.  (Do not
 * assume otherwise!)  Also, any arguments retrieved after the end of
 * the string are also NULL arguments.
 *
 * Also unlike strtok,  the delimiter is a single character (which makes
 * the most sense in this application anyway, besides being simpler).
 *
 * Also unlike strtok, this IS THREAD SAFE.  A pointer, lastp, is returned
 * for use on the next scan.  It should be used as the next 'string' 
 * argument.  'string' and 'lastp' may be the same variable.
 *
 *
 * Finally, the strings returned are Malloc'd strings, which can (and
 * should be if possible) freed after use.
 */
    char *nextarg = NULL;
    char *endp = string;
    char save;

    while (*endp != delimiter && *endp != '\0') endp++;

    if (endp != string)
    {
	save = *endp;
	*endp = '\0';
	nextarg = Strdup (string);
	*endp = save;
    }
    if (*endp == delimiter)
    {
	*lastp = endp+1; /* Advance past delimiter */
    }
    else
    {
	*lastp = endp;   /* Don't go past terminating \0 */
    }
    return nextarg;
}


/*
 * Currently, str_lower is internal use by StringCmp only.
 *   Needs a 'safe' shell for external use.
 */

char *str_lower (char *string)
{
    char *ptr = string;
    while (*ptr != '\0')
    {
	*ptr = tolower (*ptr);
	ptr++;
    }
    return string;
}

/* StringCmp ()
 * A 'safe' string comparitor with case sensitivity switch
 *
 * Unlike strcmp, this permits either or both of its arguments to be
 * null pointers.  Such a NULL string will match either another
 * null string or an empty string.  If the other string has characters,
 * it will be greater.
 *
 * Also unlike strcmp, there is a case insensitivity switch.  Effectively
 * this causes both strings to be downcased before comparison, which
 * may also affect +/- values of strings which don't match.
 *
 */
int StringCmp (const char *string1, const char *string2, int case_sensitive)
{
    int cmp = 0;  /* This forms the default for many NULL cases */

    if (!string1 || !string2)    /* If either is NULL */
    {
	if (string1 || string2)  /* If one is NULL but not other */
	{
	    if (!string1)        /* If string1 is NULL */
	    {
		if (string2[0] != '\0')
		{
		    cmp = -1;    /* string1 is NULL and string2 is non-empty */
		}
	    }
	    else
	    {
		if (string1[0] != '\0')
		{
		    cmp = 1;     /* string2 is NULL and string1 is non-empty */
		}
	    }
	}
    }
    else if (case_sensitive)
    {
	cmp = strcmp (string1, string2);
    }
    else
    {
#if 0
/* My original inefficient version */
	char *lower1 = str_lower (Strdup (string1));
	char *lower2 = str_lower (Strdup (string2));
	cmp = strcmp (lower1, lower2);
	free (lower1);
	free (lower2);
#elif defined  __SUNPRO_C
	cmp = strcasecmp (string1, string2);
#else
/* Generic version */
/* Thanks to Ben Pfaff and Chris Torek and others from usenet postings */
	for (;;)
	{
	    if (!*string1 ||
		tolower ((unsigned char) *string1) != 
		tolower ((unsigned char) *string2))
	    {
		cmp = (tolower ((unsigned char) *string1) -
		       tolower ((unsigned char) *string2));
		break;
	    }
	    string1++;
	    string2++;
	}
#endif
    }
    return cmp;
}

/*
 * A non-proprietary basename function
 */
char *Basename (char *pathname)
{
    int filename_start = 0;
    int last_filename_start = 0;
    int i;

    if (!pathname || !strlen(pathname))
    {
	return (Strdup ("."));
    }

    for (i = 0; pathname[i] != '\0'; i++)
    {
	if (pathname[i] == '/')
	{
	    last_filename_start = filename_start;
	    filename_start = i+1;
	}
    }
    if (filename_start == i)
    {
	filename_start = last_filename_start;
    }
    return Strdup (&pathname[filename_start]);
}

/*****************************************************************/
/*   FUNCTION CFPRINTF - print centered string              */
/*   Warning...takes only one argument for argument list         */
/*****************************************************************/

#define LINE_WIDTH 80

char *EmptyString = "";

void _center_string (char *outstring, char *instring);

void cfprintf (FILE *file, const char *string, const void *arg)
{
    char buf[BUFSIZ];
    char cbuf[BUFSIZ];

    if (!arg) arg="";  /* BIG source of error removed */
    if (!string) string = "";
    if (!file) file = stderr;
    sprintf (buf, string, arg);
    _center_string (cbuf, buf);
    fprintf (file, cbuf);
}

void _center_string (char *outstring, char *instring)
{
    int i, j;
    int prenewlinecount = 0;
    int postnewlinecount = 0;
    int spaces;

    for (i = 0; instring[i] == ' ' || instring[i] == '\t'
	 || instring[i] == '\n'; i++)
    {
	if (instring[i] == '\n')
	{
	    prenewlinecount++;
	}
    }

    for (j = strlen(instring) - 1; j >= 0; j--)
    {
	if (instring[j] != '\n') break;
	postnewlinecount++;
    }

    spaces = (LINE_WIDTH - (int) (strlen (&instring[i]) - postnewlinecount))
	      / 2;

    for (j = 0; j < prenewlinecount; j++) outstring[j] = '\n';
    for (j = prenewlinecount; j < spaces+prenewlinecount; j++)
    {
	outstring[j] = ' ';
    }
    outstring[j] = '\0';

    strcat (outstring, &instring[i]);
}

/*
 * Returns 0 on success, 1 on error
 */
int TermSize (int *lines, int *columns)
{
    struct winsize ws;

    if (ioctl (0, TIOCGWINSZ, &ws))
    {
	return 1;
    }
    *lines = ws.ws_row;
    *columns = ws.ws_col;
    return 0;
}
    
