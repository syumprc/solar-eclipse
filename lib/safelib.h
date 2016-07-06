/* safelib.h
 * safe memory allocation and other generally useful functions
 *   (returns don't have to be checked; they are checked here)
 *
 * comments: Requires that there be a function named 'error' to which a
 *   string is passed.  Also writes directly to standard error stream (2)
 *   to avoid additional dynamic memory usage if there is a memory
 *   allocation error.
 */

#ifndef SAFELIB_H
#define SAFELIB_H

#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <signal.h>

#ifdef __cplusplus
extern "C" {
    enum case_sensitivity {case_ins, case_sens};
#endif

extern char *EmptyString;
void *smalloc (size_t size);
void *scalloc (size_t nelem, size_t elsize);
void *srealloc (void *ptr, size_t size);
char *sstrdup (const char *str);
void sclose (int fdes);
void sfclose (FILE *fptr);
void ssigset (int sig, void (*disp)(int));
void ssighold (int sig);
void ssigrelse (int sig);
char *string__append (char **stringpp, const char *appendum);
void spadf (char *cstring, int padlen);
char *sfgets (char *buffer, size_t buflen, FILE *file);
char *sgets (char *buffer);
char *NextArg (char *string, char delimiter, char **lastp);
char *NextToken (char *buf, char *delim, char **next);
char *RestOfTokens (char *buf, char *delim, char **next);
int StringCmp (const char *string1, const char *string2, int case_sensitive);
char *Basename (char *pathname);
void cfprintf (FILE *file, const char *string, const void *arg);
int TermSize (int *lines, int *columns);

#ifdef __cplusplus
}
#endif


/*
 * A useful macro (if not entirely safe) for FORTRAN interfacing
 * 'max' for two arguments is already provided, but not for 3.
 * We could do 4 args, etc., but 3 probably hits almost all
 * cases beyond 2.
 */

#define max3(a,b,c) \
((a) > (b)) ? (((a)>(c))?(a):(c)) : (((b)>(c))?(b):(c))
    
/*
 * Define nice names for most of these functions (that don't already
 * have nice names...)
 */

#define Malloc smalloc
#define Calloc scalloc
#define Realloc srealloc
#define Strdup sstrdup
#define Close sclose
#define Fclose sfclose
#define Sigset ssigset
#define Sighold ssighold
#define Sigrelse ssigrelse
#define StringAppend string__append
#define Fgets sfgets
#define Gets sgets
#define Centered_Fprintf cfprintf
#define Strcmp(a,b) StringCmp ((a), (b), case_ins)
#endif



