/*
 * Filename: plotpipe.c
 * Purpose:  drive xmgr or compatible plot program through named pipe
 * 
 * Derived from acegr_np.c
 * Modified by Charles Peterson, SFBR
 */


/*************************************************************************/
/*  (C) 04.08.1997 Henrik Seidel (HS)                                    */
/*  <henrik@itb.biologie.hu-berlin.de>                                   */
/*************************************************************************/

#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <unistd.h>
#include <signal.h>
#include <errno.h>
#include <sys/param.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef HAVE_FCNTL_H
#  include <fcntl.h>
#endif
#ifdef HAVE_SYS_WAIT_H
#  include <sys/wait.h>
#endif

#include "plotpipe.h"

/* static global variables */
static char* buf = NULL;               /* global write buffer */
static int bufsize;                    /* size of the global write buffer */
static int bufsizeforce;               /* threshold for forcing a flush */
static char fifo[MAXPATHLEN] = "";     /* name of the fifo */
static int fd_fifo = -1;               /* file descriptor of the fifo */
static pid_t pid = (pid_t) -1;         /* pid of acegr */
static int installed_unlink_fifo = 0;      /* allow multiple sequential sessions */

/* Delete named pipe if it exists */
void unlink_fifo_void (void)
{
    if (*fifo != 0) {
        unlink (fifo);
        *fifo = 0;
    }
    if (pid > 0)
    {
	void (*last_sigfunc)(int);
	last_sigfunc = signal (SIGCHLD, SIG_IGN);
	kill (pid, SIGKILL);
	signal (SIGCHLD, last_sigfunc);
/*	fprintf (stderr, "Killed plot process"); */
	pid = -1;
    }
}

/* Interface used by atexit or on_exit */
static void
#if (defined(HAVE_ON_EXIT))
unlink_fifo (int status, void* arg)
#else
unlink_fifo (void)
#endif
{
    unlink_fifo_void ();
}

static void
handler (const int signum)
{
    /* Call the exit procedure to ensure proper deleting of the named */
    /* pipe (we registered unlink_fifo with atexit) */
    exit (EXIT_FAILURE);
}

static void
handler2 (const int signum)
{
    unlink_fifo_void ();
}


int
ACEgrOpen (const int arg, char *xmgr_name)
{
    int ratexit = 0;

    /* Only one simultaneous session allowed */
    if (fd_fifo != -1) {
	return -2;
    }

    /* Set the buffer sizes according to arg */
    if (arg < 64) {
        fprintf (stderr, "The buffer size in ACEgrOpen should be >= 64\n");
        return -1;
    }
    bufsize = arg;
    bufsizeforce = arg / 2;
    
    /* Don't exit on SIGPIPE */
    signal (SIGPIPE, SIG_IGN);
    
    /* Make sure that fifo is removed on exit and interruptions */
    if (!installed_unlink_fifo) {
#if (defined(HAVE_ON_EXIT))
	on_exit (unlink_fifo, NULL);
	installed_unlink_fifo = 1;
#else
	ratexit = atexit (unlink_fifo);
	if (ratexit) {
	    fprintf (stderr, 
"Warning: unable to build handler to clean up named pipe on exit\n");
	} else {
	    installed_unlink_fifo = 1;
	}
#endif
    }

    if (signal (SIGINT, handler))  {} /* fprintf (stderr,
"Warning: non-default SIGINT handler being superceded for ACE/gr\n"); */
    if (signal (SIGTERM, handler)) fprintf (stderr,
"Warning: non-default SIGTERM handler being superceded for ACE/gr\n");
    if (signal (SIGQUIT, handler)) fprintf (stderr,
"Warning: non-default SIGQUIT handler being superceded for ACE/gr\n");
/*  signal (SIGHUP, handler2);
    signal (SIGCHLD, handler); */

    /* Make the named pipe */
    if (mkfifo (tmpnam(fifo), 0600)) {
        perror ("ACEgrOpen");
        return (-1);
    }

/*  fprintf (stderr, "fifo is named %s\n", fifo); */

    /* Fork a subprocess for starting acegr */
    pid = vfork ();
    if (pid == (pid_t) (-1)) {
        perror ("ACEgrOpen");
        return (-1);
    }

    /* If we are the child, replace ourselves with acegr */
    if (pid == (pid_t) 0) {
        /* we use a timer delay of 0.9 seconds -- it should not be an
           integer number since this might produce race conditions
           with sleep() */
	int retval;

	int null_desc;
	null_desc = open ("/dev/null", O_WRONLY);
	close (2);
	dup (null_desc);

	{
	    int appleX11 = 0;

#ifdef __APPLE__

	    struct stat astat;
	    
	    if (stat ("/usr/X11R6/bin/Xquartz", &astat))
	    {
		if (stat ("/usr/X11R6/NoSolarWarning", &astat))
		{
		    printf (
"\nApple X11 not installed; /usr/X11R6/bin/Xquartz not found\n");
		    printf (
"You must install Apple X11 from install disk or start alternate X11 first\n");
		    printf (
"If X11 from XFree86 is running, a mouse click may be required to position\n");
		    printf (
"and open xmgr window now.  Read README.Mac for more details.\n");
		    printf (
"Create file named /usr/X11R6/NoSolarWarning to disable this warning\n\n");
		}
            }
	    else
	    {
		appleX11 = 1;
		if (stat ("/usr/X11R6/NoSolarWarning", &astat))
		{
		    printf (
"Please be patient...this may take 30 seconds first time\n");
		    printf (
"If X11 creates an XTerm when starting, you may delete it\n");
		}
	    }

#endif
/* end if __APPLE__ */

	    if (appleX11)
	    {
		retval = execlp ("open-x11a", "open-x11a", xmgr_name, "-noask",
				 "-npipe", fifo, "-timer", "900", "-maxgraph",
				 "30", (char *) NULL);
	    }
	    else
	    {
		retval = execlp (xmgr_name, xmgr_name, "-noask", "-npipe",
				 fifo, "-timer", "900", "-maxgraph", "30",
				 (char *) NULL);
	    }
	}
        if (retval == -1) {
	    fprintf (stdout, "Could not start plot program: %s\n", xmgr_name);
            exit (-1);
        }
    }

    /* We are the parent -> open the fifo for writing and allocate the */
    /* write buffer */
    fd_fifo = open (fifo, O_WRONLY);
    if (fd_fifo == -1) {
        perror ("ACEgrOpen");
        return -1;
    }
    buf = (char *) malloc ((size_t) bufsize);
    if (buf == NULL) {
        fprintf (stderr, "ACEgrOpen: Not enough memory\n");
	unlink_fifo_void ();
        return -1;
    }
    *buf = 0;
    return 0;
}

int
ACEgrClose (int terminated_by_user)
{
    void (*last_sigfunc)(int);
    int status = 0;

    /* Verify that pipe is open */
    if (fd_fifo == -1) {
	fprintf (stderr, "ACE/gr not open\n");
	return -1;
    }

    /* Don't exit when our child "acegr" exits */
    last_sigfunc = signal (SIGCHLD, SIG_IGN);
    /* Tell acegr to exit and wait until it did so */
    if (!terminated_by_user)
    {
	ACEgrCommand ("exit");
	ACEgrFlush ();
    }
    waitpid (pid, NULL, 0);
    pid = -1;

    /* Free the buffer, close the fifo, and delete the fifo */
    if (buf) free (buf);
    buf = NULL;
    if (fd_fifo != -1) close (fd_fifo);
    fd_fifo = -1;
    unlink_fifo_void();

    signal (SIGCHLD, last_sigfunc);
    signal (SIGINT, SIG_DFL);
    signal (SIGTERM, SIG_DFL);
    signal (SIGQUIT, SIG_DFL);
    return (0);
}

int
ACEgrFlush (void)
{
    int loop = 0;
    int len, res;

    /* Verify that pipe is open */
    if (fd_fifo == -1) {
	fprintf (stderr, "ACE/gr not open\n");
	return -1;
    }

    len = strlen (buf);
    while (loop < 10) {
        res = write (fd_fifo, buf, len);
        if (res == len) {
            /* We successfully wrote everything -> clear the buffer
               and return */
            *buf = 0;
            return (0);
        } else {
            /* There was an error, we could not write the buffer */
            if (errno == EPIPE) {
                /* Wait a second, since every 1.8 seconds there is a 0.9 s
                   time window for writing */
                sleep (1);
            } else {
                /* There was an error we cannot handle */
                perror ("ACEgrFlush");
                return -1;
            }
        }
        loop++;
    }
/*    fprintf (stderr, "ACEgrFlush: ran into eternal loop\n"); */
    return -2;
}

int
ACEgrCommand (const char* cmd)
{
    int len, res;
    
    /* Verify that pipe is open */
    if (fd_fifo == -1) {
	fprintf (stderr, "ACE/gr not open\n");
	return -1;
    }

    /* Append the new string to the global write buffer */
    strncat (buf, cmd, bufsize);
    strncat (buf, "\n", bufsize);
    len = strlen (buf);
    
    /* Try to send the global write buffer to acegr */
    res = write (fd_fifo, buf, len);
    if (res < len) {
        /* We could not write the buffer */
        if (errno == EPIPE) {
            /* The reason is that we got a SIGPIPE, we can handle
               this. If the filling of the global write buffer exceeds
               some threshold, flush the buffer. This involves waiting
               for a time window to write to acegr. */
            if (strlen (buf) >= bufsizeforce) {
                if (ACEgrFlush () == EOF) {
                    fprintf (stderr,
                             "ACEgrCommand: Can't flush write buffer\n");
                    return (-1);
                }
            }
        } else {
            /* There was some error we cannot handle */
            perror ("ACEgrCommand");
            return (-1);
        }
    } else {
        /* We could write the entire buffer -> clear it */
        *buf = 0;
    }
    return (0);
}



