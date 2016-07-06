/* pipeback.c
 *   pipeback...return fd piped to a program's stdout
 *   add_argv...build variable list of arguments
 *
 */

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>

#include "safelib.h"
#include "pipeback.h"

#include <signal.h>

int local_pipe[2];

int PipebackChildPID = 0;

static int Usage = 0;

FILE *pipeback_shell_open (const char *filename, const char *argv[])
{
    int exec_type = 1;
    int fdesc = pipeback_shared (filename, argv, exec_type);
    FILE *pfile = fdopen (fdesc, "r");
    if (!pfile)
    {
/*	close (fdesc);  Not sure if this is safe */
	return 0;
    }
    Usage++;
    return pfile;
}

void pipeback_shell_close (FILE *file)
{
#ifndef NOCLDWAIT
    int wait_status = 0;
#endif
    int cresult = fclose (file);
    if (cresult) fprintf (stderr, "Warning.  Error closing pipeback file.\n");
    if (0 == --Usage)
    {
	if (SIG_ERR == sigset (SIGCHLD, SIG_DFL))
	{
	    fprintf (stderr, 
		 "Warning.  Error resetting SIGCHLD action after pipeback.\n");
	}
#ifndef NOCLDWAIT
    wait (&wait_status);
#endif
    }
}

int pipeback (const char *filename, const char *argv[])
{
    int exec_type = 0;
    return pipeback_shared (filename, argv, exec_type);
}

int pipeback_shared (const char *filename, const char *argv[], int exec_type)
{
    struct sigaction *sigact;
    char mbuf[256];

    sigact = scalloc (sizeof (struct sigaction), 1);
    sigact->sa_handler = SIG_IGN;
#ifdef NOCLDWAIT
    sigact->sa_flags = SA_NOCLDWAIT;
#endif
    if (sigaction (SIGCHLD, sigact, NULL))
    {
	error ("Error protecting against zombies!");
	exit (EXIT_FAILURE);
    }

    if (pipe (local_pipe))
    {
	error ("Error opening pipeback pipe!");
	exit (EXIT_FAILURE);
    }

    if (PipebackChildPID = fork())
    {
	sclose (local_pipe[1]);
	return local_pipe[0];
    }
    else
    {
	sclose (local_pipe[0]);
	sclose (1);
	dup (local_pipe[1]);
	sclose (2);
	dup (local_pipe[1]);
	sclose (local_pipe[1]);
	if (0 == exec_type)
	{
	    execv (filename, (char**) argv);
	}
	else
	{
	    execvp (filename, (char**) argv);
	}
	sprintf (mbuf, "Pipeback was unable to load program %s", filename);
	error (mbuf);
    }
}

void argv_add (const char **argv[], const char *argp)
{
    int length = 0;

    if (!(*argv))
    {
	*argv = (const char **) smalloc (2 * sizeof (char*));
	(*argv)[0] = argp;
	(*argv)[1] = NULL;
    } 
    else
    {
	int i = 0;

	while ((*argv)[i]) i++;
	(*argv)[i] = argp;
	*argv = (const char **) srealloc (*argv, (i+2) * sizeof (char*));
	(*argv)[++i] = NULL;
    }
}

int getcrc (const char* filename, unsigned int* crc, unsigned int* noct)
{
    int error = 1;  // set to 0 on success

    FILE* cfile;
    const char* argv[3];
    argv[0] = "cksum";
    argv[1] = filename;
    argv[2] = 0;

    cfile = pipeback_shell_open ("cksum", argv);

    if (!cfile)
    {
	fprintf (stderr, "\nError opening %s for checksum\n", filename);
    }
    else
    {
	char buf[256];
	if (!fgets (buf, 256, cfile))
	{
	    fprintf (stderr, "\nError computing checksum on %s\n", filename);
	}
	else
	{
	    int count = sscanf (buf, "%u %u", crc, noct);
	    if (count != 2)
	    {
		fprintf (stderr, "\nError getting checkum on %s\n", filename);
	    }
	    else
	    {
		error = 0;
	    }
	}
	pipeback_shell_close (cfile);
    }
    return error;
}
