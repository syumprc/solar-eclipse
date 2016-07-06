/* error.c
 * default error function for libsafe.a
 *
 * (Now that we're using a "real" library archive,
 *  we can have a default.)
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>

void error (char *error_message)
{
    fprintf (stderr, "%s\n", error_message);
    exit (EXIT_FAILURE);
}
