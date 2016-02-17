/* error.c
 * default error function for libsafe.a
 *
 * (Now that we're using a "real" library archive,
 *  we can have a default.)
 *
 *      Copyright (C) 1997 Southwest Foundation for Biomedical Research
 *                          All rights reserved.
 *                 Absolutely no warranty, express or implied.
 *
 * author:  Charles P. Peterson
 *   date:  November 4, 1997
 *
 */

#include <stdio.h>
#include <stdlib.h>

void error (char *error_message)
{
    fprintf (stderr, "%s\n", error_message);
    exit (EXIT_FAILURE);
}
