/* pipeback.c
 *   pipeback...return fd piped to a program's stdout
 *   add_argv...build variable list of arguments
 *
 *      Copyright (C) 1995 Southwest Foundation for Biomedical Research
 *                          All rights reserved.
 *                 Absolutely no warranty, express or implied.
 *
 * author:  Charles P. Peterson
 *   date:  June 23, 1995
 *
 */

#ifndef PIPEBACK_H
#define PIPEBACK_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

int pipeback (const char *filename, const char *argv[]);
FILE *pipeback_shell_open (const char *filename, const char *argv[]);
void pipeback_shell_close (FILE *file);
void argv_add (const char **argv[], const char *argp);
void error (const char *message);
int pipeback_shared (const char *filename, const char *argv[], int exec_type);
int getcrc (const char* filename, unsigned int* crc, unsigned int* noct);

#ifdef __cplusplus
}
#endif

#endif
