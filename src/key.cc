/*
 * key.cc implements the key command and validate function
 * Written by Charles Peterson beginning on March 4, 1998
 * Copyright (c) 1997 Southwest Foundation for Biomedical Research
 *
 */

#include <stdio.h>
#include <time.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
// tcl.h from solar.h
#include "solar.h"
#include "pipeback.h"

static const int ROTOR_SIZE = 32;
static const int COMPOSEED = 29;
static const int MAXID = 256;
static const int KEYLEN = 8;
static const int MIN_DEPUTY_KEYLEN = 10;

static char *make_key (char *username, unsigned seed);
static void shuffle (char *dsttable, char *srctable, int *random, int number);
static char* validate_key_from_deputy (char* username, char* key);
static char* deputy_make_key (char* username, unsigned seed);
static char* make_deputy_key (char* username, unsigned seed);
int validate_solar (Tcl_Interp *interp);
static int key_validated = 0;



// Damn GCC compiler bug workaround

extern "C" char* ccuserid (char*);

// Condor parallel system flag

extern "C" int Condor;
extern "C" int NisKey;

// Command line Key

extern char Key[10];

static void keyhelp () 
{
    printf (
"Usage:    key make username\n");
    printf ("\n");
}

// This table contains EASY to recognize letters, number, and special chars
//   Hard to recognize characters such as 0/letter O are excluded

const int ezchar_len = 32;  // not including last fake char
static char ezchar[33] = {
'a',
'b',
'c',
'd',
'e',
'f',
'g',
'h',
'i',
'j',
'k',
'm',
'n',
'p',
'q',
'r',
's',
't',
'u',
'v',
'w',
'x',
'y',
'z',
'2',
'3',
'4',
'5',
'6',
'7',
'8',
'9',
'?'};


// This table contains random numbers 1-32
static int ezrandom[32] = {
18,
14,
10,
15,
24,
23,
12,
9,
28,
32,
5,
20,
29,
22,
25,
3,
26,
11,
27,
7,
13,
2,
21,
4,
8,
30,
17,
1,
31,
6,
19,
16,
};

static unsigned int masterlock[8] =
{
    21,
    19,
    4,
    28,
    108,
    23,
    31,
    80
};

static unsigned int mastermasterlock[8] =
{
    27,
    108,
    37,
    15,
    29,
    23,
    3,
    111
};

extern "C" int KeyCmd (ClientData clientData, Tcl_Interp *interp,  int argc, char *argv[])
{
	
	RESULT_LIT("1");
	return TCL_OK;
    if (argc == 2 && !Strcmp (argv[1], "help"))
    {
	keyhelp ();
	return TCL_OK;
    }
    if (argc == 2 && !Strcmp (argv[1], "-whoami"))
    {
	char* username = ccuserid(NULL);
	RESULT_BUF (username);
	return TCL_OK;
    }

    if (argc == 2 && !Strcmp (argv[1], "-is-key-on-command-line?"))
    {
	if (strlen (Key))
	{
	    RESULT_LIT ("1");
	} 
	else
	{
	    RESULT_LIT ("0");
	}
	return TCL_OK;
    }

// Tcl-command key validation comes in two flavors.  The first, 
// key -check, is a check for whether a key was detected at startup
// or the last maximization.  It is used by the "register" command
// before creating the key file so that the user doesn't get all
// sorts of error messages during the normal execution of the register
// command.

// The key -validate command does the full validation.  This is used
// by the register command AFTER a key file has been created to be sure
// it was validated correctly.  Previous validation status is cleared
// to be sure the key is re-read.

// The startup and maximize routines call the function
// validate_solar() directly.  It accepts previous successful validation
// status so that the key does not have to be re-read for each maximization.

    if (argc == 2 && !Strcmp (argv[1], "-check"))
    {
	if (key_validated == 1)
	{
	    RESULT_LIT ("1");
	}
	else
	{
	    RESULT_LIT ("0");
	}
	return TCL_OK;
    }

    if (argc == 2 && !Strcmp (argv[1], "-validate"))
    {
	key_validated = 0;
	int status = validate_solar (interp);
	if (status == 1)
	{
	    RESULT_LIT ("1");
	}
	else
	{
	    RESULT_LIT ("0");
	}
	return TCL_OK;
    }


    char masterkey[32];
    int i;
    for (i=0; i<8; i++)
    {
	if (masterlock[i] > 31)
	{
	    masterkey[i] = (int) masterlock[i];
	}
	else
	{
	    masterkey[i] = ezchar[ezrandom[(int) masterlock[i]]];
	}
    }
    masterkey[i] = '\0';

    if (argc == 3)
    {
	if (!StringCmp (argv[1], "make", case_ins))
	{
	    for (i=0; i<8; i++)
	    {
		if (mastermasterlock[i] > 31)
		{
		    masterkey[i] = (int) mastermasterlock[i];
		}
		else
		{
		    masterkey[i] = ezchar[ezrandom[(int) mastermasterlock[i]]];
		}
	    }
	    masterkey[i] = '\0';
	    system ("date");
	    char *password = getpass ("Enter master password> ");
	    time_t rawtime = time (NULL);
	    struct tm *tmp = localtime (&rawtime);
	    int date = tmp->tm_mday;
	    int code = date + 11;
	    char target[32];
	    sprintf (target, masterkey, code);
	    if (!strcmp (password, target))
	    {
		unsigned seed = time (NULL) % ROTOR_SIZE;
		char *newkey = make_key (argv[2], seed);
		char buf[256];
		sprintf (buf, "%s", newkey);
		RESULT_BUF (buf);
		return TCL_OK;
	    }
	}
	else if (!StringCmp (argv[1], masterkey, case_sens))
	{
	    unsigned seed = time (NULL) % ROTOR_SIZE;
	    char *newkey = make_key (argv[2], seed);
	    char buf[256];
	    sprintf (buf, "%s", newkey);
	    RESULT_BUF (buf);
	    return TCL_OK;
	}
	fprintf (stderr, "Don't do this again.\n");
	exit (EXIT_FAILURE);
    }

// to make key FOR deputy:

    if (argc >= 3 && !strcmp (argv[1], masterkey) &&
	!Strcmp (argv[3], "deputy"))
    {
	    unsigned seed = time (NULL) % ROTOR_SIZE;
	    char *newkey = make_deputy_key (argv[2], seed);
	    RESULT_BUF (newkey);
	    free (newkey);
	    return TCL_OK;
    }

// Handle deputy make key command

    if (argc >= 4 && !Strcmp (argv[1], "deputy") &&
	!Strcmp (argv[2], "7solar7"))
    {

	unsigned seed = time (NULL) % ROTOR_SIZE;
	char *newkey = deputy_make_key (argv[3], seed);
	if (!newkey)
	{
	    RESULT_LIT ("Error making key");
	    return TCL_ERROR;
	}
	RESULT_BUF (newkey);
	free (newkey);
	return TCL_OK;
    }

    return TCL_ERROR;
}



static void shuffle (char *dsttable, char *srctable, int *random, int number)
{
    for (int i = 0; i < number; i++)
    {
	dsttable[i] = srctable[random[i]-1];
    }
}

static char *make_key (char *username, unsigned seed)
{
    int i;
    char cryptable[ROTOR_SIZE+1];
    shuffle (cryptable, ezchar, ezrandom, ROTOR_SIZE);
    cryptable[ROTOR_SIZE] = '\0';

    char inbuf[9];
    inbuf[8] = '\0';
    strncpy (inbuf, username, 8);
    for (i = strlen (inbuf); i < 9; i++) inbuf[i] = '\0';
	
    char outbuf[9];
    for (i = 0; i < 8; i++) {
	seed += inbuf[i] + 29 - i*i;
	seed %= ROTOR_SIZE;
	outbuf[i] = cryptable[seed];
//	printf ("%d,%d,%c\n", i, seed, cryptable[seed]);
    }
    outbuf[8] = '\0';
    return Strdup (outbuf);
}    


static int validate_key (char *username, char *key)
{
    int i;
    return 1;
    for (i = 0; i < ROTOR_SIZE; i++)
    {
	char *test_key = make_key (username, i);
	if (!strcmp (test_key, key))
	{
	    free (test_key);
	    return 1;
	}
	free (test_key);
    }
    return 0;
}


int validate_solar (Tcl_Interp *interp)
{
    char filename[1024];

// See if already validated

    // **** modified by Sung 11/24/2015 - skip the key validation ****
    //if (key_validated) return 1;
	return 1;
    // ***************************************************************
// Set up key buffer

    char keyline[1024];  // from file, beware of whitespace
    char keybuf[9];     // after trimming and/or deputy conversion
    keybuf[0] = '\0';
    keybuf[8] = '\0';
    bool key_from_command = false;

// If key provided to command line, use that

    if (strlen(Key))
    {
	strncpy (keyline, Key, 8);
	key_from_command = true;
    }
    else
    {

// Key not provided to command line, so obtain from file (if not Condor or NIS)

	if (Condor)
	{
	    fprintf (stderr, "    *** Condor key must be specified with -key argument\n",
		     filename);
	    fprintf (stderr, "    *** solar -condor -key key\n\n",
		     filename);
	    return 0;
	}
	else if (NisKey)
	{
	    fprintf (stderr, "    *** NIS key must be specified with -niskey argument\n",
		     filename);
	    fprintf (stderr, "    *** solar -niskey key\n\n",
		     filename);
	    return 0;
	}

// Try to read key from ~/.solar_reg file

	int i;
	sprintf (filename, "%s/.solar_reg", getenv ("HOME"));
	FILE *keyfile = fopen (filename, "r");
	if (!keyfile)
	{
	    fprintf (stderr, "    *** Missing key file %s\n\n",filename);

// No keyfile
// Print registration instructions

	    if (TCL_OK != Solar_Eval (interp, "please_register"))
	    {

// Print backup registration instructions

		fprintf (stderr,
"Solar has not been registered, or has not been correctly registered.\n");
		fprintf (stderr,
"To register, send a mail message to solar@txbiomedgenetics.org specifying\n");
		fprintf (stderr,
"for each user: (1) username, (2) email address, (3) NIH grant number, if\n");
		fprintf (stderr,
"applicable.\n");
	    }

// Return unregistered status

	    return 0;
	}
	else
	{

// File opened, now read from it and close

	    if (!fgets (keyline, 256, keyfile))
	    {
		fprintf (stderr, "\n    *** Empty key file %s\n\n",filename);
		fclose (keyfile);
		return 0;
	    }
	    fclose (keyfile);
	}
    }

// Get username

    char* username;

// Condor username

    if (Condor)
    {
	int condorstatus = Solar_Eval (interp,"getcondoruid");
	if (TCL_OK != condorstatus)
	{
	    RESULT_LIT ("");
	    printf ("Error getting Condor userid from .job.ad file\n\n");
	    return 0;
	}
	username = Strdup (interp->result);
	RESULT_LIT ("");
    }
    else if (NisKey)
    {
	int niskeystatus = Solar_Eval (interp,"exec nisdomainname");
	if (TCL_OK != niskeystatus)
	{
	    RESULT_LIT ("");
	    printf ("Error getting nisdomainname\n");
	    return 0;
	}
	char* domainname = Strdup (interp->result);
//	printf ("nisdomainname is %s\n",domainname);
	username = domainname;
	char* lastperiod = strrchr (domainname, '.');
	if (lastperiod)
	{
	    *lastperiod = '\0';
	    char* nextperiod = strrchr (domainname, '.');
	    if (nextperiod)
	    {
		username = nextperiod+1;
	    }
	}
//	printf ("checking %s\n",username);
	username = Strdup (username);
	free (domainname);
	RESULT_LIT ("");
    }
    else
    {

// Regular username

#ifndef LINUX64
    username = ccuserid(NULL);
#else
    int userstatus = Solar_Eval(interp,
"exec whoami");
    if (TCL_OK != userstatus) {
	RESULT_LIT ("");
	printf ("Error accessing username with who command!\n\
Unable to validate SOLAR key!\n");
	return 0;
    }
    username = Strdup (interp->result);
    RESULT_LIT ("");
#endif

    }

// Apply deputy key conversion if this is a at key

    if (keyline[0] == '+')
    {
	int i;
	for (i = 1; i < 256; i++)
	{
	    if (!isgraph(keyline[i]))
	    {
		keyline[i] = '\0';
		break;
	    }
	}
	char* newkey = validate_key_from_deputy (username, keyline);
	if (newkey)
	{
	    strncpy (keyline, newkey, 256);
	    free (newkey);
	}
	else
	{
	    return 0;
	}
    }

    strncpy (keybuf, keyline, 8);

// Check for match
    if (strlen (keybuf))
    {
	char userbuf[9];
	strncpy (userbuf, username, 8);
	userbuf[8] = '\0';
	if (validate_key (userbuf, keybuf))
	{
	    key_validated = 1;
	    return 1;
	}
	fprintf (stderr,
        "\n    *** Username: %s   Key: %s\n",username,keybuf);
	if (key_from_command)
	{
	    if (Condor)
	    {
		fprintf (stderr, 
	 "    *** Key from -key argument doesn't match Condor username\n\n",filename);
	    }
	    else
	    {
	    fprintf (stderr, 
	 "    *** Key from -key argument doesn't match username\n\n",filename);
	    }
	}
	else
	{
	    fprintf (stderr, 
	 "    *** Key in file %s doesn't match username\n\n",filename);
	}	    
    }
    return 0;
}


// New functions to implement "deputy" registration scheme as documented
// for deputy command.

		
static int bigrandom[32][32] = {
    {12,23,19,0,16,15,1,6,8,2,26,28,24,27,30,14,17,11,5,18,10,4,7,22,31,3,13,25,29,9,20,21},
    {0,17,20,15,14,24,26,12,11,4,5,22,19,8,30,7,31,13,1,10,2,3,29,28,23,9,27,6,25,18,16,21},
    {5,31,12,27,9,3,22,0,8,13,18,20,26,21,6,11,24,28,1,19,30,17,10,4,15,2,7,29,25,14,16,23},
    {21,30,17,24,18,4,20,28,23,26,27,15,16,13,8,14,7,12,6,22,29,5,25,19,1,2,11,9,31,3,10,0},
    {15,23,8,31,21,13,22,6,14,11,10,3,26,27,1,19,0,18,25,4,2,9,28,30,5,29,7,20,16,24,17,12},
    {5,10,30,12,2,3,19,23,24,8,27,1,26,18,22,9,31,15,16,14,17,29,0,20,25,28,6,11,13,4,7,21},
    {5,28,13,4,23,2,31,27,20,7,8,0,16,26,21,22,10,25,9,30,12,18,29,1,11,19,14,24,3,6,17,15},
    {16,12,29,14,27,9,26,7,3,25,1,17,30,24,19,15,11,6,8,10,23,0,28,31,20,21,5,22,13,4,2,18},
    {30,24,21,29,20,15,10,22,1,4,3,7,26,16,18,6,19,28,8,17,0,11,27,2,14,13,9,12,25,31,5,23},
    {6,5,24,12,16,19,20,23,13,15,7,18,2,11,31,28,3,1,29,30,8,27,4,17,10,14,9,25,26,21,22,0},
    {17,30,19,11,13,7,26,27,12,22,2,25,10,20,24,16,15,0,3,28,21,1,18,29,31,9,5,4,23,6,14,8},
    {11,20,7,17,13,15,29,9,6,23,8,19,18,31,4,0,2,22,14,16,27,30,25,21,28,5,12,26,10,24,1,3},
    {13,29,25,26,8,10,17,15,6,4,2,11,0,12,27,16,28,22,7,18,31,5,19,14,20,1,9,3,24,30,23,21},
    {1,17,12,13,21,0,4,25,31,26,16,28,22,29,15,6,3,10,9,11,8,7,27,18,19,24,2,14,20,23,30,5},
    {17,14,23,27,16,6,7,12,22,24,18,21,30,28,20,3,13,29,5,11,19,25,10,2,26,8,1,4,15,9,31,0},
    {9,22,3,13,19,31,18,2,8,5,28,17,14,10,20,16,4,25,0,15,27,11,6,26,24,7,30,21,23,29,12,1},
    {0,12,16,9,14,4,17,3,8,29,10,7,13,6,31,5,18,20,26,28,1,11,15,24,2,25,22,21,27,30,19,23},
    {23,8,5,12,1,24,31,10,4,27,16,19,21,28,22,30,17,18,26,13,29,15,20,7,6,14,0,11,2,9,3,25},
    {15,3,8,16,12,27,21,20,31,10,26,22,2,30,5,13,7,14,25,11,6,29,9,19,18,0,4,23,28,17,1,24},
    {5,16,31,30,18,8,26,29,13,25,24,7,4,21,14,2,22,20,10,17,1,6,0,3,19,11,28,12,23,27,15,9},
    {30,27,12,28,15,26,2,22,9,21,1,4,24,25,31,14,13,5,3,23,29,7,18,20,17,10,16,8,19,0,11,6},
    {2,30,1,12,3,13,5,21,20,14,6,22,23,11,0,16,9,27,29,17,28,31,26,25,18,4,8,10,19,15,24,7},
    {6,25,9,7,22,24,1,0,26,16,5,20,10,12,30,14,27,8,19,28,31,13,17,29,18,11,21,2,23,4,15,3},
    {0,12,26,21,24,29,25,11,27,31,15,14,22,23,30,8,28,18,13,17,6,5,10,9,3,2,16,19,7,4,20,1},
    {19,29,5,6,11,14,21,26,8,25,20,0,4,18,9,2,27,3,10,30,24,17,13,23,1,28,31,12,22,7,16,15},
    {30,21,0,1,7,11,17,6,10,15,29,25,20,28,2,18,3,31,13,16,8,22,4,14,12,19,9,26,24,23,5,27},
    {25,11,1,31,16,0,24,22,14,30,9,2,7,15,21,29,6,5,26,13,8,12,4,20,27,28,3,18,19,23,17,10},
    {30,20,31,22,24,2,26,16,23,21,13,17,27,29,14,10,5,3,4,25,12,15,0,28,1,11,7,8,6,9,18,19},
    {26,30,9,8,15,31,0,17,22,24,13,18,3,21,12,11,14,20,6,2,16,5,29,28,4,27,19,25,10,7,23,1},
    {15,10,23,3,31,26,5,14,16,29,7,4,1,19,13,0,9,20,17,12,30,2,6,22,28,21,27,8,25,11,24,18},
    {5,7,4,20,13,31,15,29,16,0,12,21,11,9,26,14,3,25,10,23,1,30,6,8,24,19,2,28,27,17,18,22},
    {15,29,6,2,0,28,22,19,16,3,17,12,27,24,30,8,18,4,10,7,23,25,26,14,11,21,20,9,13,31,1,5}
};


// Function to find character in ezchar and return index
int dechar (char testchar)
{
    for (int i = 0; i < ezchar_len; i++)
    {
	if (testchar == ezchar[i])
	{
	    return i;
	}
    }
    return -1;
}

// Function to decode by finding number in random array and returning index
int decode (int* codez, int testnum)
{
    for (int i = 0; i < ezchar_len; i++)
    {
	if (testnum == codez[i])
	{
	    return i;
	}
    }
    return -1;
}

// Function to get ezchar from ezchar array while
// protecting against out-of-range

char ezchars(int index)
{
    if (index >= 0 && index < ezchar_len)
    {
	return ezchar[index];
    }
    else
    {
	return ezchar[ezchar_len]; // returns fake character
    }
}


// procedure to make a deputy-key

static char* make_deputy_key (char* username, unsigned seed)
{
    char *newkey = make_key (username, seed);
//  printf ("unencrypted key is %s\n", newkey);

// Modify into deputy key

    char deputykey[KEYLEN+1];  // Version of key to be validated
    for (int i = 0; i < KEYLEN; i++)
    {
	int keynum = dechar (newkey[i]);
	if (keynum < 0)
	{
	    free (newkey);
	    fprintf (stderr, "\n    *** Deputy key generation error\n");
	    return 0;
	}
	int newrandom = bigrandom[i+11][keynum];
//	printf ("Keynum is %d, random is %d\n", keynum, newrandom);
	deputykey[i] = ezchars(newrandom);
    }
    deputykey[KEYLEN] = 0;
    free (newkey);
    return Strdup (deputykey);
}


// procedure which makes limited keys for deputy
static char* deputy_make_key (char* username, unsigned seed)
{
    int i;
    char* deputyid = ccuserid(NULL);
    int didlen = strlen(deputyid);
    char deputykey[KEYLEN+1];
    char normal_userkey[KEYLEN+1];
    char interkey[KEYLEN+1];
    char outputkey[MAXID];

// retrieve deputy key from ~/.solar_deputy
    char filename[1024];
    sprintf (filename, "%s/.solar_deputy", getenv ("HOME"));
    FILE *keyfile = fopen (filename, "r");
    if (!keyfile)
    {
	fprintf (stderr, "\n    *** Missing key file ~/.solar_deputy\n\n");
	return 0;
    }

// File opened, now read from it and close
    char keyline[256];
    if (!fgets (keyline, 256, keyfile))
    {
	fprintf (stderr, "\n    *** Empty key file ~/.solar_deputy\n\n");
	fclose (keyfile);
	return 0;
    }
    strncpy (deputykey, keyline, 8);
    deputykey[KEYLEN] = 0;
    fclose (keyfile);

// Modify into standard user key for standard validation
    char deputyvalkey[KEYLEN+1];  // Version of key to be validated
    for (i = 0; i < KEYLEN; i++)
    {
	int keynum = dechar (deputykey[i]);
	if (keynum < 0)
	{
	    fprintf (stderr, "\n    *** Invalid key file ~/.solar_deputy\n\n");
	    return 0;
	}
	int decoded = decode (bigrandom[i+11], keynum);
//	printf ("Keynum is %d, random is %d\n", decoded, keynum);
	deputyvalkey[i] = ezchars(decoded);
    }
    deputyvalkey[KEYLEN] = 0;

//   printf ("unencrypted key is %s\n", deputyvalkey);

// Check for match with current username

    char userbuf[9];
    strncpy (userbuf, deputyid, 8);
    userbuf[8] = '\0';
    if (!validate_key (userbuf, deputyvalkey))
    {
	fprintf (stderr, 
  "\n    *** Key in key file ~/.solar_deputy doesn't match username\n\n");
	return 0;
    }

// make normal user key and free temp storage
    {
	char* userkey = make_key (username, seed);
	strncpy (normal_userkey, userkey, 8);
	free (userkey);
    }

    normal_userkey[KEYLEN] = 0;
//  printf ("Normal key is %s\n", normal_userkey);

// reshuffle seed for use in creating limited key

    seed = bigrandom[seed][dechar(normal_userkey[2])];
    int starting_seed = seed;

//* Now, compute visible interkey from seed, deputykey, and normal_userkey
    for (i = 0; i < KEYLEN; i++)
    {
	int start = dechar(deputykey[i]);
	seed = (start + seed) % ROTOR_SIZE;
	int step = dechar(normal_userkey[i]);
	int encoded = bigrandom[seed][step];
	interkey[i] = ezchars(encoded);
    }

    interkey[KEYLEN] = 0;
//  printf ("Interkey is %s\n", interkey);

//* Put it all together into composite limited key
    int oindex = 0;

//** First character is plus
    outputkey[oindex++] = '+';

//** First character is encoded seed
    outputkey[oindex++] = ezchars(starting_seed);
    seed = starting_seed;

    int dindex = 0;
    int iindex = 0;
    while ((dindex < didlen || iindex < 8) && oindex+2 < MAXID)
    {
//** Next (and alternating) is encoded deputy ID
	if (dindex < didlen)
	{
	    char didchar = deputyid[dindex++];
	    int ezindex;
	    if (-1 != (ezindex = dechar (didchar)))
	    {
		outputkey[oindex++] = ezchars(bigrandom[seed][ezindex]);
		seed++;
		seed %= 32;
	    }
	    else
	    {
//*** If ID character isn't in ezchar, it has to be put into key verbatim
		outputkey[oindex++] = didchar;
	    }
	}

//** Finally, (and alternating) is interkey
	if (iindex < KEYLEN)
	{
	    outputkey[oindex++] = interkey[iindex++];
	}
    }
    outputkey[oindex] = '\0';
    return Strdup (outputkey);
}


// Procedure which extracts normal key from key from deputy
static char* validate_key_from_deputy (char* username, char* key)
{
    int kindex = 1;  // skip over plus

// Get random seed from first character
    int seed = dechar (key[kindex++]);
    if (seed < 0)
    {
	return 0;
    }
    int starting_seed = seed;

    int didlen = strlen(key) - (KEYLEN + 2);
    char deputyid[MAXID];
    char interkey[KEYLEN+1];

    int dindex = 0;
    int iindex = 0;
    while ( kindex < strlen(key))
    {

// Extract deputyid

	if (dindex < didlen)
	{

//** if in character set, decode
	    int ezindex;
	    if (-1 != (ezindex = dechar(key[kindex++])))
	    {
		deputyid[dindex++] = ezchars(decode(bigrandom[seed], ezindex));
		seed++;
		seed %= 32;
	    }
	    else
	    {

//** if not in character set, it passes through
		deputyid[dindex++] = key[kindex-1];
	    }
	}

// Extract interkey
	if (kindex < strlen(key) && iindex < KEYLEN)
	{
	    interkey[iindex++] = key[kindex++];
	}
    }	    


// Read deputy file as documented for deputy command

    char filename[1024];
    char* deputy_home = getenv ("SOLAR_DEPUTY_HOME");
    if (deputy_home)
    {
	strncpy (filename, deputy_home, 1024);
    }
    else
    {
	deputyid[dindex] = '\0';
//	printf ("deputyid is %s\n", deputyid);
	strncpy (filename, getenv ("HOME"), 1024);
	char* lasts = strrchr (filename, '/');
	if (!lasts)
	{
	    error ("Invalid HOME shell variable for validation");
	}
	*(lasts+1) = '\0';
	strncat (filename, deputyid, 1024);
    }
    strncat (filename, "/.solar_deputy", 1024);
    FILE *dfile = fopen (filename, "r");
    if (!dfile)
    {
	fprintf (stderr, "\n    *** Registration key in deputy plus format\n");
	fprintf (stderr,   "    *** Expected filename: %s\n", filename);
	fprintf (stderr,   "    *** Deputy authorization file not found!\n\n");
	return 0;
    }
    char deputykey[256];
    if (!fgets (deputykey, 256, dfile))
    {
	fprintf (stderr, "\n    *** Empty deputy authorization file!");
	return 0;
    }

// Now decode user key using deputy key and interkey
    seed = starting_seed;
    char userkey[KEYLEN+1];
    int i;
    for (i = 0; i < KEYLEN; i++)
    {
	int start = dechar(deputykey[i]);
	seed = (start + seed) % ROTOR_SIZE;
	int target = dechar(interkey[i]);
	int decoded = decode (bigrandom[seed], target);
	userkey[i] = ezchars(decoded);
    }

// Return userkey to continue regular validation

    userkey[KEYLEN] = 0;
    return Strdup (userkey);
}





