/*
 * ccsearch.cc forms a C++ wrapper around the fisher 'search' subroutine
 * (and related subroutines) which form the core of fisher
 *
 * Written by Charles Peterson beginning on November 10, 1997
 * Copyright (c) 1997 Southwest Foundation for Biomedical Research
 */

#include <signal.h>
#include "safelib.h"
#include "solar.h"

// #define OLD86 1    // Define for old Solaris 4.x compiler used on x86
                   // Fortran dynamic arrays need to be allocated here

int *Vbinary = 0;  // Used here and by Covariate::savemeans


int howclose (const double d1, const double d2);
const int min_sig_digits = 5;

/*
 * FISHER constants
 */
#define STRLEN 40
#define NAMELEN 40
#define LINELEN 80
#define BIGNAMELEN 128
const int SUPERNAMELEN = 1024;
#define FORMATLEN 100
#define VNAMELEN 96 /* Using first three PEDSYS mnemonics */

#define MAX_VARIABLES MAX_PHENOTYPES
int TseIndex[MAX_VARIABLES] = {0};  /* Index to Standard Error for Trait */

/* Fortran interface to get index to trait std err */
extern "C" int getitse_ (int *ntrait) 
{
    return TseIndex[*ntrait-1];
}

static int vtraits = 0;
static int if_any_discrete = 0;

extern "C" int getvirt_ ()
{
    return vtraits;
}


void cfisher (void);

extern "C" void preped_ (...);
extern "C" void copyped_ (...);
extern "C" void optima_ (...);
extern "C" void search_ (...);
extern "C" void smpoutput_ (...);
extern "C" void dcopyped_ (...);
extern "C" double ppnd_ (double *mu, int *ierr);



// For machines (Alpha) w/o IEEE numbers (e.g. NaN) we must trap FP exceptions
#ifdef NOIEEE
static bool sigfpe_raised = FALSE;
void hsigfpe (int i)
{
    sigfpe_raised = TRUE;
}
#endif


double ccsearch (const char* outfilename, int get_who, bool sampledata)
{
/*
 * "Constants:"
 * Fisher says "Don't mess with these in ordinary circumstances."
 */
    if_any_discrete = 0;
    const int conout = 6;
    int mtrait = MAX_TRAITS;
    int mxtwin = 32000;
    double absent = MISSING_PHENOTYPE;
    double conv = Option::get_double ("Conv");
    int nconv = Option::get_int ("NConv");
    double tol = Option::get_double ("Tol");
    int mxstep = Option::get_int ("MaxStep");
    double dp = 1.0E-7;
    int unit1 = 1;
    int unit2 = 2;
    int unit3 = 3;
    int unit4 = 4;
    int unit7 = 7;
    const char *model_filename = "solar";
/*
 * Variables, variables, and more variables
 */    
    const int error_message_len = 1024;
    char error_message[error_message_len] = {0};
    int error_message_lnblnk = 0;

    size_t m0, m1, m2, m3, m4, m5, m6, m7, m8, m9;
    size_t m10, m11, m12, m13, m14, m15, m16, m17, m18, m19;
    size_t m20,m21,m22,m23,m24,mqd,m25, m_end;
#ifdef OLD86
    size_t m26;
#endif
    size_t par_start, par_se;
    int p2, p3, p4, p5;
    int n2, n3, n4, n5, n6, n7, n8, n9;
    int maxvar, mwork;
    int stand = Option::get_int ("StandLogLike");
    int travellen_ = 8;
    char travel[8+1];
    int npoint = Option::get_int ("GridPoints");
    int npar = Parameter::count ();
    int maxpar = (npar > 1) ? npar : 1;
    double *parest = (double *) Calloc (maxpar, sizeof (double));
    int ncnstr = Constraint::count();
    int mcnstr = (ncnstr > 1) ? ncnstr : 1;
    int maxtab = ncnstr + npar + 1;
    int diag = Option::get_int ("CMDiagonal");
    int asycv = Option::get_int ("StandErr");
    int irobst = Option::get_int ("RobustSE");
    int mxiter = Option::get_int ("MaxIter");
    int outlie = Option::get_int ("Outlier");
    double percut = Option::get_double ("CutPeople");
    double pedcut = Option::get_double ("CutPed");
    int normal = (Option::get_int ("TDist")) ? 0 : 1;
    int echo = 0;
    int verbosity = Verbosity::get();
    int formatlen_ = FORMATLEN;
    char format1[FORMATLEN+1];
    char format2[FORMATLEN+1];
    int namelen_ = NAMELEN;
    int bignamelen_ = BIGNAMELEN;    /* Should be migrating to this */
    char title[NAMELEN+1];
    char obs_pedfil[NAMELEN+1];
    char modfil[NAMELEN+1];
    char updfil[BIGNAMELEN+1];
    double qdform = 0.0L;
    Loglike::quadratic (0.0);
    double tarray[2] = {0.0, 0.0};
    int nsamp = 0;     // Sample size (w/o probands)
    int nprob = 0;     // Number of probands
    int nsampip = 0;   // Sample size _with_ probands
    int ibig = 0;
    int iter = 0;
    int last = 0;
    int inasyc = 0;
    double loglik = 0.0L;
    int maxpeo = FisherPedigree::Highest_Pedcnt ();
    int NumIncPED = 0;
    int momega = 0;
    int iprob = 1;
    int status = 0;
    double *vmin = NULL;
    double *vmax = NULL;
    double *vmean = NULL;
    double RawSD = 0;
    double RawMU = 0;
    int scorindx = 0;
    int TPeopIP = 0; // Total people in included pedigrees

    int nextra = Option::get_int ("LenExtra");
    int shortstringlen_ = STRLEN;
    double Big = -1e20;
    double *extra = (double *) Calloc (sizeof (double), nextra);
    int DiscreteTraits = 0;

#ifdef NOIEEE
    sigfpe_raised = FALSE;
    signal (SIGFPE,hsigfpe);
#endif
/*
 * We save the first person in each pedigree for later identification
 */
    int num_all_peds = FisherPedigree::Number_of_Pedigrees ();
    char *FirstIDs = FisherPedigree::Allocate_FirstIDs ();

/*
 * Setup variables: traits and covariates
 */
    bool premax = verbose ("PREMAX");
    if (premax) fprintf (stderr, "    **  Setting up traits and covariates\n");


// Get number of traits and vectors of trait names and indexes

    Phenotypes::start_setup ();

    int lastvar = 0;
    int ntrait = 0;
    int nvar = 0;
    char *vnames = (char *) Calloc (sizeof (char[VNAMELEN]), MAX_VARIABLES);
    char *sampledataheader = Strdup("id,famid");

    int nindep = 0;

    int i;
    for (i = 0; i < Trait::Number_Of(); i++)
    {
	const char* t = Trait::Name (i);
	strncpy (&vnames[VNAMELEN*nvar], t, VNAMELEN);
	fortspad (&vnames[VNAMELEN*nvar], VNAMELEN);
	Expression *exp = Trait::expression (i);
	if (!exp)
	{
	    Phenotypes::setup (t);
	}
	else
	{
	    Phenotypes::setup_virtual_position ("");
	}
	if (sampledata) {
	    char tempname[60];
	    sprintf (tempname, ",trait%d", i+1);
	    StringAppend (&sampledataheader, tempname);
	}
	nvar++;
	ntrait++;
    }

// If there is more than one trait, see if it is virtual or real
// Virtual traits are handled by creating individual-traits
// After reading in data, we roll back ntrait to 1

    vtraits = 0;
    if (-1 != Option::get_int ("UnbalancedTraits") && ntrait > 1)
    {
	vtraits = ntrait;
    }

// Setup IBDID pseudo-variable; filled in by PINPUT.F
// index: ntrait+1, as seen by pinput

    Phenotypes::setup_virtual_position ("");

    strncpy (&vnames[VNAMELEN*nvar], "IBDID", VNAMELEN);
    fortspad (&vnames[VNAMELEN*nvar], VNAMELEN);
    if (sampledata)
    {
	StringAppend (&sampledataheader, ",ibdid");
    }
    nindep++;
    nvar++;
    
// Setup PROBND field, if any (pseudo-variable created anyway if none)
// index: ntrait+2, as seen by pinput

    bool found_proband = false;
    if (Phenotypes::available (Field::Map ("PROBND")))
    {
	Phenotypes::setup_proband (Field::Map ("PROBND"));
	found_proband = true;
    }
    else if (!Strcmp ("PROBND", Field::Map ("PROBND")))
    {

//  If not field-mapped (maybe to -none) try also PROBAND and PRBAND

	if (Phenotypes::available ("PROBAND"))
	{
	    Phenotypes::setup_proband ("PROBAND");
	    found_proband = true;
	}
	else if (Phenotypes::available ("PRBAND"))
	{
	    Phenotypes::setup_proband ("PRBAND");
	    found_proband = true;
	}
    }

// Now we adjust for absence of proband in get_indexed_phenotype
// So we just setup a name placeholder with setup_virtual_position
    if (!found_proband)
    {
	Phenotypes::setup_virtual_position ("");
    }

    strncpy (&vnames[VNAMELEN*nvar], "PROBND", VNAMELEN);
    fortspad (&vnames[VNAMELEN*nvar], VNAMELEN);
    if (sampledata)
    {
	StringAppend (&sampledataheader, ",probnd");
    }
    nindep++;
    nvar++;

// Get intra-individual standard deviations, if any, and add to vectors

    for (i = 0; i < Trait::Number_Of(); i++)
    {
	const char* sv = Trait::Name_Stdev (i);
	if (sv)
	{
	    strncpy (&vnames[VNAMELEN*nvar], sv, VNAMELEN);
	    fortspad (&vnames[VNAMELEN*nvar], VNAMELEN);
	    if (sampledata)
	    {
		StringAppend (&sampledataheader, ",");
		StringAppend (&sampledataheader, sv);
	    }
	    Phenotypes::setup (sv);
	    nindep++;
	    nvar++;
	}
    }

// Add Covariates

    Covariate *c;
    Phenotypes::Initialize_Expressions();
    for (i = 0; c = Covariate::index(i); i++)
    {
	CovariateTerm *term;
	for (term = c->terms(); term; term = term->next)
	{
	    if (term->sex) continue;
	    if (term->expression)
	    {
		term->packed_number = nvar;
		Phenotypes::setup_virtual_position (term->name);
		Phenotypes::Setup_Expression (term->expression, nvar);
		strncpy (&vnames[VNAMELEN*nvar], term->name, VNAMELEN);
		fortspad (&vnames[VNAMELEN*nvar], VNAMELEN);
		if (sampledata)
		{
		    StringAppend (&sampledataheader, ",");
		    StringAppend (&sampledataheader, term->name);
		}
		nindep++;
		nvar++;
	    }
	    else
	    {
		int old_pos;
		if (-1 == (old_pos = Phenotypes::find_virtual_position (term->name)))
		{
		    term->packed_number = nvar;
		    Phenotypes::setup (term->name);
		    strncpy (&vnames[VNAMELEN*nvar], term->name, VNAMELEN);
		    fortspad (&vnames[VNAMELEN*nvar], VNAMELEN);
		    if (sampledata)
		    {
			StringAppend (&sampledataheader, ",");
			StringAppend (&sampledataheader, term->name);
		    }
		    nindep++;
		    nvar++;
		} 
		else 
		{
		    term->packed_number = old_pos;
		}
	    }

// If covar is qualified (but not null qualified)
// make note that it applies to only one trait
	    int j;
	    if (term->trait && strlen(term->trait))
	    {
// Find the trait that it is qualified to
		for (j = 0; j < Trait::Number_Of(); j++)
		{
		    if (!Strcmp (Trait::Name(j), term->trait))
		    {
			Phenotypes::Covar_Of_Trait (term->packed_number, j);
		    }
		}
	    }

// If covariate is unqualified, or null qualified, applies to all traits

	    else
	    {
		for (j = 0; j < Trait::Number_Of(); j++)
		{
		    Phenotypes::Covar_Of_Trait (term->packed_number, j);
		}
	    }
	}
    }

// Add variables required by Mu  and Omega equations and expression traits
//   unless they are already included
// Also setup packed indexes so they can be retrieved by Mu and Omega

    const char *vname;
    for (i = 0; vname = EqVar::index (i); i++)
    {
	int term_pos;
	if (-1 != (term_pos = Phenotypes::find_virtual_position (vname)))
	{
	    EqVar::set_packed_index (i, Phenotypes::find_virtual_position 
				     (vname));
	}
	else
	{
	    Phenotypes::setup (vname);
	    EqVar::set_packed_index (i, nvar);
	    strncpy (&vnames[VNAMELEN*nvar], vname, VNAMELEN);
	    fortspad (&vnames[VNAMELEN*nvar], VNAMELEN);
	    if (sampledata)
	    {
		StringAppend (&sampledataheader, ",");
		StringAppend (&sampledataheader, vname);
	    }
	    term_pos = nvar;
	    nindep++;
	    nvar++;
	}

// Check if this is needed by a particular trait.
// If so, flag only that partuclar trait

	if (-1 != EqVar::for_trait(i))
	{
	    Phenotypes::Covar_Of_Trait (term_pos, EqVar::for_trait(i));

	} else {	    

// Mu terms NOT already declared as covariates are assumed to be needed
// by ALL traits.  (Otherwise declare covariate specific to those traits
// that need it, and build up mu from scratch if necessary, which makes
// covariate statement moot except for phenotype requirement issue.)

	    for (int j = 0; j < Trait::Number_Of(); j++)
	    {
		Phenotypes::Covar_Of_Trait (term_pos, j);
	    }
	}
    }

    int NV2READ = 5 + lastvar;    // Last variable we need in each ped rec
/*
 * If only scoring one index,
 *   Determine FORTRAN index into actual parameter
 */
    int score_cindex =  Option::get_int ("ScoreOnlyIndex");
    if (score_cindex > -1)
    {
	char h2name[64];
	if (score_cindex == 0)
	{
	    sprintf (h2name, "h2r");
	}
	else
	{
	    sprintf (h2name, "h2q%d", score_cindex);
	}
	Parameter *p = Parameter::find (h2name);
	scorindx = p->number() + 1;
	mxiter = 1;
    }
/*
 * Now that we know how many variables we need to read, and the
 * number of parameters, we can allocate carray and nsampa (*new*)
 */
    if (premax) fprintf (stderr, "    **  Allocating CARRAY\n");
    int carray_size = 1 + shortstringlen_ * ((NV2READ > Parameter::count()) ? 
                                             NV2READ : Parameter::count());
    char *carray = (char*) Calloc (1, carray_size);
    int *nsampa = (int*) Calloc (sizeof (int), nvar);
    int *ttrasamp = (int*) Calloc (sizeof (int), ntrait);
    int *ntrasamp = (int*) Calloc (sizeof (int), ntrait);
    Vbinary = (int*) Calloc (sizeof (int), nvar);

/* #define CARRAY_DEBUG */
#ifdef CARRAY_DEBUG
    carray[carray_size-1] = -47;
#endif

// Determine size needed for iarray for PREPED
// (All rounded up to 4byte boundary, plus extra for test)
    int evenpeo = maxpeo + (maxpeo % 2);
    int LENI = (5*evenpeo)+(2*mxtwin)+nvar;
// Space for virtual trait arrays...father, mother, group
    LENI += vtraits*evenpeo*3;
    LENI += (LENI % 2) + 2;

    if (premax) fprintf (stderr, "    **  Allocating IARRAY\n");
    int *iarray = (int *) Calloc (sizeof (int), LENI);
    const int TEST_NUMBER = -1029478311;
    iarray[LENI-1] = TEST_NUMBER;   // Used for overwrite test
    if (verbosity & 512)
    {
	fprintf (stderr, "IARRAY allocation: %d\n",LENI);
    }

// "DEPVAR" put into iarray
// (In our case, DEPVAR is pretty primitive, since compaction is done.)
    for (i = 0; i < ntrait; i++)
    {
	iarray[i] = i+1;
    }

    strncpy (title, "Solar Model", NAMELEN);
    title[NAMELEN] = '\0';
    fortspad (title, NAMELEN);

// updfil, obs_pedfil, and modfil are now obsolete, but
// certain interfaces still require them

    strncpy (updfil, "solar.upd", BIGNAMELEN);
    updfil[BIGNAMELEN] = '\0';
    fortspad (updfil, BIGNAMELEN);

    strncpy (obs_pedfil, " ", NAMELEN);
    obs_pedfil[NAMELEN] = '\0';
    fortspad (obs_pedfil, NAMELEN);
    
    strncpy (modfil, model_filename, NAMELEN);
    modfil[NAMELEN] = '\0';
    fortspad (modfil, NAMELEN);
/*
 * TRAVEL
 */
    if (0 == Option::get_int ("Grid"))
    {
	strncpy (travel, "SEARCH", travellen_ + 1);
    }
    else
    {
	strncpy (travel, "GRID", travellen_ + 1);
    }
    fortspad (travel, travellen_);


/**************************************************************************\
 * Read the Pedigree File (some of this code originated from INPUT.F)
\**************************************************************************/

    p2 = nvar;        // Note: in FORTRAN, p2 was nvar+1
    p3 = nvar + p2;
    p4 = nvar + p3;
    p5 = nvar + p4;

// pm2 is effective m2 inside preped (where rarray starts at rarray[p5])
    int mvar = maxpeo*nvar;
    int pm2 = mxtwin + p5;
    int pm3 = 2*maxpeo + pm2;
    int pm4 = maxpeo + pm3;
    int pm5 = mvar + pm4;
    int pm6 = nvar + pm5;
    int pm7 = nvar + pm6;
    int pm8 = nvar + pm7;
    int pm9 = nvar + pm8;
    size_t Lenr = pm9 + maxpeo*nvar*vtraits + 2;
// Lenr for really big numbers
	
    if (premax) fprintf (stderr, "    **  Allocating RARRAY\n");
    double *rarray = (double *) Calloc (sizeof (double), Lenr);
    rarray[Lenr-1] = TEST_NUMBER;
    if (verbosity & 512)
    {
	fprintf (stderr, "RARRAY allocation: %lu\n", Lenr);
    }

    const char *pedfilename = FisherPedigree::pedfilename();
    const char *phenfilename = Phenotypes::filenames();

    int who = 0;
    if (get_who)
    {
	who = 1;
	FILE *whofile = fopen ("who.out", "w");
	if (!whofile)
	{
	    free (sampledataheader);
	    throw Safe_Error_Return ("Unable to open who.out");
	}
	fclose (whofile);
    }

    if (verbose ("ITERATE"))
    {
	cfprintf (stdout, "Pedigree:  %s\n", (void*) pedfilename);
	cfprintf (stdout, "Phenotypes:  %s\n", (void*) phenfilename);
	printf ("\n");
	cfprintf (stdout, "Merging pedigree and phenotype data...\n\n",
		  (void*) EmptyString);
    }

    const int titlelen_ = 255;
    char pedfile[titlelen_+1];
    strncpy (pedfile, pedfilename, titlelen_+1);
    fortspad (pedfile, titlelen_);

    char phenfile[SUPERNAMELEN+1];
    strncpy (phenfile, phenfilename, SUPERNAMELEN+1);
    fortspad (phenfile, SUPERNAMELEN);
    
    int writedat = 0;
    if (sampledata)
    {
	writedat = 1;

	FILE *datafile = fopen ("sampledata.out", "w");
	if (!datafile)
	{
	    free (sampledataheader);
	    throw Safe_Error_Return ("Unable to open sampledata.out");
	}
	fprintf (datafile, "%s\n", sampledataheader);
	fclose (datafile);
    }
    free (sampledataheader);

    int iquit = 0;
    int matrix_errors = 0;

// Write header for EVD output file if this is evd phase 1
    if (1==Option::get_int ("EVDPhase"))
    {
	FILE* evdoutfile = fopen ("evddata.out","w");
	fprintf (evdoutfile,"id,fa,mo,sex,sex_evd,tmean,lambda");
	for (int ivar = 0; ivar < nvar; ivar++)
	{
	    bool namedone = false;
	    char namebuf[VNAMELEN+1];
	    strncpy (namebuf,&vnames[ivar*VNAMELEN],VNAMELEN);
	    namebuf[VNAMELEN] = 0;
	    for (int ichar=0; ichar < VNAMELEN; ichar++)
	    {
		if (namebuf[ichar] == ' ')
		{
		    namebuf[ichar] = '\0';
		    if (Strcmp (namebuf, "IBDID") &&
			Strcmp (namebuf,"PROBND"))
		    {
			fprintf (evdoutfile,",%s_evd",namebuf);
		    }
		    namedone = true;
		    break;
		}
	    }
	    if (!namedone)
	    {
		fprintf (evdoutfile,",%s",namebuf);
	    }
	}
	fprintf (evdoutfile,"\n");
	fclose (evdoutfile);
    }

// Do phenotype input in Fortran

    try {
    preped_
      (&rarray[p5], rarray, &rarray[p2], &rarray[p3], &rarray[p4],
       iarray, carray, &absent, &LENI, &Lenr, &maxpeo, &mxtwin,
       &NumIncPED, &nsamp, &nvar, &unit1, &unit7, &unit3, format1, format2,
       &echo, &unit4, title, obs_pedfil, modfil,
       &conout, &NV2READ, &ntrait, &nindep, vnames, nsampa,
       &iquit, Vbinary, pedfile, phenfile,
       &nprob, &nsampip, FirstIDs, &num_all_peds,
       error_message, &error_message_lnblnk, &TPeopIP,
       &RawSD, &RawMU, &who, &vtraits, ntrasamp, ttrasamp, &writedat,
       shortstringlen_,   /* Each 'element' in carray is 8 characters */
       formatlen_,        /* format1 */
       formatlen_,        /* format2 */
       namelen_,          /* title */
       namelen_,          /* pedfil */
       namelen_,          /* modfil */
       VNAMELEN,          /* vnames */
       titlelen_,         /* pedfile */
       SUPERNAMELEN,      /* phenfile */
       PERMID_LEN,        /* FirstIDs */
       error_message_len
	  );
       if (!Zscore::If_Ready())
       {
	   throw Safe_Error_Return ("Trap Zscore Data Needed");
       }
    }
    catch (Safe_Error_Return& ser)
    {
	strncpy (error_message, ser.message(), error_message_len);
	error_message_lnblnk = (int) strlen (ser.message ());
	iquit = 1;
    }
    if ((matrix_errors = Matrix::Check_Matrices()))
    {
	iquit = 1;
    }

    if (!iquit && get_who!=1  && !sampledata) {
/*
 *  ** NEW **
 *  We set up two new arrays, explicitly named vmin and vmax, to
 *  save the minimum and maximum values of each variable.  These
 *  can be used later to set min/max bounds (they are passed to SEARCH)
 *
 *  vmean is also set up to contain variable mean values
 *  (previous VMEAN in rarray  is only set up after copyped/dcopyped)
 */
    vmin = (double *) Calloc (sizeof (double), nvar);
    vmax = (double *) Calloc (sizeof (double), nvar);
    vmean = (double *) Calloc (sizeof (double), nvar);

    try {

    {
	int i;
	for (i = 0; i < nvar; i++)
	{
	    vmin[i] = rarray[p3+i];
	    vmax[i] = rarray[i];
	    vmean[i] = rarray[p2+i];
	}
    }
    Mu::vmins = vmin;
    Mu::vmaxs = vmax;

/*
 * If there are virtual traits, we inflate the number of individuals
 * and deflate the number of traits
 */
    if (vtraits)
    {
	maxpeo = maxpeo * vtraits;
	nsampip = nsampip * vtraits;
	ntrait = 1;
    }

/*
 * Set up things we can set up now that we've read the pedigree file
 */
    momega = (maxpeo * ntrait) + 1;
    maxvar = (maxpeo * nvar);
    mwork = max3 (maxpar*(momega-1), momega, maxtab);
/*
 * Setup for discrete traits version and quantitative traits version differs
 *   quantitative trait version uses Ken's "search"
 *   discrete traits version uses John's "optima" (a "search" rewrite)
 */
    if (Strcmp("Default",Option::get_string("ModelType")) ||
	(Trait::Maximization_Type() == discrete && 
	 Option::get_int ("EnableDiscrete")))
    {
/*
 * This is the discrete traits version
 */
	DiscreteTraits = 1;
/*
 * Check for proper coding of discrete and non-single value (allowed by option)
 * now done in phenotypes.cc (was once done here).
 */
/*
 * Fill up data array with all variable data
 * (discrete traits version uses arrays instead of a scratch file)
 */
	int VLEN = nsampip * nvar;
	double* VARDATA = (double*) Calloc (VLEN, sizeof (double));
	int* NASCER = (int*) Calloc (NumIncPED, sizeof (int));
	int* NIND = (int*) Calloc (NumIncPED, sizeof (int));
	int* NCUMIND = (int*) Calloc (NumIncPED+1, sizeof (int));
	int* Male = (int*) Calloc (nsampip, sizeof (int));
#ifdef OLD86
	char *Disc = (char*) Calloc (maxpeo,1);
#endif
	if (LENI < maxpeo+1)
	{
	    if (iarray[LENI-1] != TEST_NUMBER)
	    {
		throw Safe_Error_Return (
     "IARRAY overrun occurred in preped; contact solar@txbiomedgenetics.org");
	    }
	    LENI = maxpeo+1;
            if (verbosity & 512)
            {
                fprintf (stderr, "\nResizing IARRAY before dcopyped to %d.\n", 
                         LENI);
            }
	    iarray = (int*) Realloc (iarray, LENI * sizeof (int));
	    iarray[LENI-1] = TEST_NUMBER;
	}

// Check that rarray is big enough for all the data for the largest pedigree

	if (maxvar+1 > Lenr)
	{
	    if (rarray[Lenr-1] != TEST_NUMBER)
	    {
		throw Safe_Error_Return (
     "RARRAY overrun occurred in preped; contact solar@txbiomedgenetics.org");
	    }
	    Lenr = maxvar+1;
            if (verbosity & 512)
            {
                fprintf (stderr,"\nResizing RARRAY before dcopyped to %lu.\n", 
                         Lenr);
            }
	    rarray = (double*) Realloc (rarray, Lenr * sizeof (double));
	    rarray[Lenr-1] = TEST_NUMBER;
	}


	char* Trttype = (char*) Malloc (nsampip);

	int n_traits = Trait::Number_Of();
	int* idiscrete = (int*) Calloc (n_traits, sizeof(int));
	for (int itrait = 0; itrait < n_traits; itrait++)
	{
	    idiscrete[itrait] = 0;
	    if (Option::get_int ("EnableDiscrete") && 
		discrete == Trait::Type(itrait))
	    {
		idiscrete[itrait] = 1;
		if_any_discrete = 1;
	    }
	}

// copy data into arrays for discrete likelihood search

	dcopyped_ (&unit1, &NumIncPED, &nvar, NASCER, NIND, VARDATA, 
		   &nsampip, rarray, &Lenr, NCUMIND, iarray, Male, &maxpeo,
		   vmean, vmin, &n_traits, idiscrete, Trttype, 1);

	int NOBS = 1;
	int NCASE = 1;
	double F = 0.0;
	double SMALL = 1.0E-10;
	int IOUNIT = unit3;           // output file
	int NPASS = 1;
	int STDERR = 0;
	if (asycv) STDERR = 2;
	int DIFFER[2];
	DIFFER[0] = 0;
	DIFFER[1] = 0;

// VMEAN and VVAR are now at the front of rarray (as for SEARCH)
// (Note that FORTRAN has 1-based arrays, C has 0-based arrays, so these
//  numbers are 1 less than they would need to be in FORTRAN.)

	m0=nvar;
	m9=m0+nvar;
	m7=maxpar+m9;
	m1=maxpar*maxpar+m7;
	m2=mcnstr*maxpar+m1;
	m3=mcnstr+m2;
	m4=maxpar+m3;
	m5=maxpar+m4;
	m6=maxpar+m5;
	m8=npoint*maxpar+m6;
	m10=NOBS*NCASE+m8;
	m11=maxpar+m10;
	m12=maxpar+m11;
	m13=maxpar+m12;
	m14=maxtab*maxtab+m13;
	m15=maxpar+m14;
	m16=maxtab+m15;
	m17=(maxpeo*maxpeo)+m16;
	m18=maxpeo+m17;
#ifndef OLD86
	m_end=maxpeo+m18;
#else
	m19=maxpeo+m18;
  	m20=NumIncPED+m19;
  	m21=NumIncPED+m20;
  	m22=maxpeo+m21;
  	m23=maxpeo+m22;
  	m24=maxpeo+m23;
  	m25=maxpeo+m24;
  	m26=maxpeo+m25;
  	m_end=maxpeo+m26;
#endif
	par_start = m9;
	par_se = m15;


        if (m_end+1 > Lenr)
        {
	    if (rarray[Lenr-1] != TEST_NUMBER)
	    {
		throw Safe_Error_Return (
     "RARRAY overrun occurred in preped; contact solar@txbiomedgenetics.org");
	    }
            Lenr = m_end+1;
            if (verbosity & 512)
            {
                fprintf (stderr, "\nResizing RARRAY before optima to %lu.\n", 
                         Lenr);
            }
            rarray = (double *) Realloc (rarray, Lenr * sizeof (double));
	    rarray[Lenr-1] = TEST_NUMBER;
        }
	Mu::vmeans = rarray;  // save pointer to VMEAN in mu module
	                      // (this is "VMEAN" in search)


	if (LENI < maxpar+1)
	{
	    if (iarray[LENI-1] != TEST_NUMBER)
	    {
		throw Safe_Error_Return (
     "IARRAY overrun occurred in preped; contact solar@txbiomedgenetics.org");
	    }
	    LENI = maxpar+1;
            if (verbosity & 512)
            {
                fprintf (stderr, "\nResizing IARRAY before optima to %d.\n", 
                         LENI);
            }
	    iarray = (int*) Realloc (iarray, LENI * sizeof (int));
	    iarray[LENI-1] = TEST_NUMBER;
	}
/*
 * Set starting value for Mu based on normal deviate (for discrete traits)
 *   (Only if start, upper, lower not already set)
 */
	
	for (i = 0; i < ntrait; i++)
	{
	    if (Trait::Type(i) == discrete)
	    {
		int imean = i * 2;
		Parameter *p = Parameter::index(imean);
		if (p->start == 0.0 && p->upper == 0.0 && p->lower == 0.0)
		{
		    int ierr = 0;

// -1.0 changed to 1.0 to conform with standard convention as per Jeff Williams

		    double disc_mu = rarray[i] - vmin[i];
		    double start =  0.0 - ppnd_ (&disc_mu, &ierr);
		    if (Verbosity::max())
		    {
			fprintf (stderr, 
				 "MU: %g   PPND: %g\n", disc_mu, start);
		    }
		    if (ierr != 0)
		    {
			throw Safe_Error_Return (
			    "Mu starting value search failed");
		    }
		    p->start = start;
		    p->upper = 3.0;
		    p->lower = (-3.0);
		    if (start > 0.0)
		    {
			p->upper = start + 3.0;
		    }
		    else
		    {
			p->lower = start - 3.0;
		    }
		}
	    }
	}

// Use separate option for discrete conv

	if (if_any_discrete)
	{
	    conv = Option::get_double ("Conv(Discrete)");
	}

// This is optima, the heart and soul of discrete traits analysis

	optima_(&rarray[m1],&rarray[m2],&rarray[m3],&rarray[m4],&rarray[m5],
		&rarray[m6],&rarray[m7],&rarray[m8],&rarray[m9],&rarray[m10],
		&rarray[m11],&rarray[m12],&rarray[m13],&rarray[m14],
		&rarray[m15],iarray,carray,&conv,&dp,&F,&SMALL,&tol,&IOUNIT,
		&maxpar,&maxtab,&mcnstr,&mxiter,&mxstep,&NCASE,&ncnstr,&nconv,
		&NOBS,&npar,&NPASS,&npoint,&STDERR,travel,DIFFER,
		rarray,&rarray[m0],vmin,vmax,&nvar,&ntrait,&unit3,&conout,
		parest,
		VARDATA,&nsampip,NIND,NASCER,&NumIncPED,&status,NCUMIND,
		&rarray[m16],&rarray[m17],&rarray[m18],&maxpeo,Male,&ibig,
		&loglik, &inasyc,&iter,&Big,
#ifdef OLD86
 		&rarray[m19],&rarray[m20],&rarray[m21],&rarray[m22],Disc,
 		&rarray[m23],&rarray[m24],&rarray[m25],&rarray[m26],
#endif
		Trttype, &vtraits,
		8,    // PNAME string length
		8     // travel string length
#ifdef OLD86
                ,1    // Disc string length (*1)
#endif
		,1    // Trttype string length
                );

	if (status == 0)
	{
	    if (loglik <= -1.0e+20)
	    {
		loglik = 0.0;
	    } 
	    else if (min_sig_digits > howclose (loglik, Big))
	    {
		char buf[512];
		sprintf (buf, 
	     "Last likelihood significantly smaller than highest likelihood:\n\
  Highest loglike: %-.11g   Last loglike: %-.11g",
			 loglik, Big);
		throw Safe_Error_Return (buf);
	    }
	    qdform = 1.0;
	    smpoutput_ 
		(modfil, &ncnstr, &npar, &rarray[m1], carray, &rarray[m9],
		 &ibig, &iter, &loglik, &qdform, &nsamp, &ntrait, &iter, 
		 &unit4, &rarray[m15], &maxtab, &maxpar, &inasyc, tarray,
		 &NumIncPED,
		 &rarray[m10], &rarray[m11], parest, &unit3, &conout, &irobst,
		 &DiscreteTraits,&Big,
		 namelen_,        // modfil
		 shortstringlen_  // Each 'element' in carray is 8 characters
		    );
	}
	free (Trttype);
#ifdef OLD86
	free (Disc);
#endif
	free (Male);
	free (VARDATA);
	free (NASCER);
	free (NIND);
	free (NCUMIND);
	double quadratic = 1.0;
	Loglike::quadratic (quadratic);

    } else {	
/*
 * This is the single quantitative trait version
 */
    int idiscq = 0;
    if (Trait::Maximization_Type() == discrete && 
	!Option::get_int ("EnableDiscrete"))
    {
	idiscq = 1;
    }

// COPYPED is more stuff from (Ken's original) INPUT.F
//   (data is packed into scratch files for re-reading, one pedigree at
//    a time.  Though slighly inefficient (unmeasurably so), this does
//    have the advantage that it is not necessary to have all data in
//    in memory at once...only one extended pedigree's worth)

	copyped_
	    (&unit1, &unit2, iarray, rarray, &LENI, &Lenr, &NumIncPED, &nvar);
/*
 * Now that we've done this, we need to reload info into the arrays
 * because it was clobbered during the preped process.
 * 
 * This duplicates what was done by FISHER.F, like it or not.
 * FISHER uses these arrays for disparate purposes at different
 * times (ouch).
 */
	int i;
	for (i = 0; i < mtrait; i++)
	{
	    iarray[i] = 0;
	}
/* 
 * Put DEPVAR back into iarray
 */
	for (i = 0; i < ntrait; i++)
	{
	    iarray[i] = i+1;
	}


	m0 = nvar;       // Note: In FORTRAN, m0 was nvar+1
	m1 = nvar + m0;
	m2 = mcnstr*maxpar + m1;
	m3 = mcnstr + m2;
	m4 = maxpar + m3;
	m5 = maxpar + m4;
	m6 = maxpar + m5;
	m7 = maxpar*npoint + m6;
	m8 = maxpar*maxpar + m7;
	m9 = maxpar*maxpar + m8;
	m10 = maxpeo + m9;
	m11 = momega*momega + m10;
	m12 = maxpar + m11;
	m13 = maxpar + m12;
	m14 = maxpar + m13;
	m15 = maxpar + m14;
	mqd = maxpar + m15;
	m16 = maxpar + mqd;
	m17 = maxpar + m16;
	m18 = maxtab*maxtab + m17;
	m19 = NumIncPED + m18;
	m20 = maxvar + m19;
	m21 = (NumIncPED + 1)*maxpar + m20;
	m22 = maxpar*maxpar + m21;
	m23 = maxpar*maxpar + m22;
	m24 = maxpar + m23;
	m25 = mwork + m24;
	par_start = m11;
	par_se = m24;
	size_t cnstr_start = m1;
	size_t cvalue_start = m2;
/*
 * Readjust Lenr if necessary (using long Lenr)
 */

	if (rarray[Lenr-1] != TEST_NUMBER)
	{
	    throw Safe_Error_Return (
     "RARRAY overrun occurred in preped; contact solar@txbiomedgenetics.org");
	}

	size_t m26a = m25 + maxpar + 1;
	if (Lenr < m26a)
	{
	    size_t LenrBytes = m26a * sizeof(double);
	    if (verbosity & 512)
	    {
		fprintf (stderr, "\nResizing RARRAY before search to %lu bytes\n", 
			  LenrBytes);
	    }
	    rarray = (double *) realloc (rarray, LenrBytes);
	    if (!rarray)
	    {
		printf ("Unable to allocate %lu bytes!\n",LenrBytes);
		size_t momegasize = (m11 - m10) * sizeof (double);
		printf ("momega was %lu bytes\n",momegasize);
		throw Safe_Error_Return (
		    "Inadequate memory to complete maximization");
	    }
	    Lenr =  m26a;
	    rarray[Lenr-1] = TEST_NUMBER;
	}
	Mu::vmeans = rarray;  // save pointer to VMEAN in mu module
	                      // (this is "VMEAN" in search)

/*
 * Set indexes into iarray
 */
	n2 = ntrait;       // Note: Was ntrait+1 in FORTRAN
	n3 = NumIncPED + n2;
	n4 = maxpeo + n3;
	n5 = maxpeo + n4;
	n6 = maxpeo + n5;
	n7 = maxpar + n6;
	n8 = NumIncPED + n7;
	n9 = maxpar + n8;
/*
 * Readjust LENI as necessary
 */
	if (iarray[LENI-1] != TEST_NUMBER)
	{
	    throw Safe_Error_Return (
     "IARRAY overrun occurred in preped; contact solar@txbiomedgenetics.org");
	}
	if (LENI < n9 + maxpeo + 1)
	{
	    LENI = n9 + maxpeo + 1;
	    if (verbosity & 512)
	    {
		cfprintf (stderr, "\nResizing IARRAY to %d.\n", (void *) LENI);
	    }
	    iarray = (int *) Realloc (iarray, LENI * sizeof (int));
	    iarray[LENI-1] = TEST_NUMBER;
	}
/*
 *    INSERT INTO rarray THE CONSTANT conv.  THIS IS NECESSARY
 *    BECAUSE SOME PC COMPILERS CHOKE WHEN A SUBROUTINE HAS MORE THAN
 *    63 ARGUMENTS.  IN SUBROUTINE SEARCH conv WILL THE FIRST ENTRY
 *    OF QDERIV.
 */

	rarray[mqd] = conv;
	
	status = 0;

// SEARCH is the heart and soul of fisher (quantitative trait)

	search_ (&rarray[m1],&rarray[m2],&rarray[m3],&rarray[m4],&rarray[m5],
    extra, &rarray[m6], &rarray[m7], &rarray[m8], &rarray[m9],
    &rarray[m10], &rarray[m11], &rarray[m12], &rarray[m13], &rarray[m14],
    &rarray[m15], &rarray[mqd], &rarray[m16], &rarray[m17], &rarray[m18],
    &rarray[m19], rarray, &rarray[m0], &rarray[m23], &rarray[m24],
    &iarray[n2], iarray, &iarray[n3], &iarray[n4], &iarray[n5], &iarray[n6],
    &iarray[n9], carray, &dp, &pedcut, &percut, &maxpar, &maxpeo, &maxtab,
    &maxvar, &mcnstr, &momega, &mwork, &mxiter, &mxstep, &ncnstr, &nconv,
    &nextra, &npar, &NumIncPED, &npoint, &nsamp, &ntrait, &nvar, &iprob,&unit2,
    &unit3, &travel, &asycv, &diag, &normal, &outlie, &stand, modfil,
    &unit4, &ibig, &iter, &last, &inasyc, &loglik, &qdform, &rarray[m20],
    &iarray[n7], &rarray[m21], &rarray[m22], &iarray[n8], &conout,
    vmin, vmax, &verbosity, &status, parest, updfil, &irobst, &tol,
    &rarray[m25], &scorindx, &RawSD, &RawMU, &Big, &idiscq,
		 shortstringlen_,   // Each 'element' in carray is 8 characters
		 travellen_,        // length of travel string
		 namelen_,          // length of modfil string
		 bignamelen_        // length of updfil string
	    );

	if (status != 0 && status != 4)
	{
	    if (verbose("MAXIMUM"))
	    {
		fprintf (stderr,"\nSEARCH Status is %d\n", status);
	    }
	    loglik = 0.0;  // If an error occurred, loglike is invalid
	} else {
//	    fprintf (stderr, "Digits: %d\n", howclose (loglik, Big));
	    if (loglik <= -1.0e+20)
	    {
		loglik = 0.0;
	    } 
	    else if (min_sig_digits > howclose (loglik, Big))
	    {
		char buf[512];
		sprintf (buf, 
	     "Last likelihood significantly smaller than highest likelihood:\n\
  Highest loglike: %-.11g   Last loglike: %-.11g",
			 loglik, Big);
		throw Safe_Error_Return (buf);
	    }
	    for (int cnum = 0; cnum < ncnstr; cnum++)
	    {
		double accum = rarray[cvalue_start+cnum];
		for (int tpar = 0; tpar < npar; tpar++)
		{
		    accum -= rarray[cnstr_start+(tpar*ncnstr)+cnum]
			*rarray[par_start+tpar];
		}
		if (0 <= Option::get_int ("EnforceConstraints") && 
		    fabs(accum) >= 1e-4)
		{
		    throw Safe_Error_Return ("Convergence failure (Numerical constraint failure)");
		}
	    }
	    smpoutput_ 
		(modfil, &ncnstr, &npar, &rarray[m1], carray, &rarray[m11],
		 &ibig, &iter, &loglik, &qdform, &nsamp, &ntrait, &last, 
		 &unit4, &rarray[m24], &maxtab, &maxpar, &inasyc, tarray,
		 &NumIncPED,
		 &rarray[m13], &rarray[m14], parest, &unit3, &conout, &irobst,
		 &DiscreteTraits, &Big,
		 namelen_,        // modfil
		 shortstringlen_  // Each 'element' in carray is 8 characters
		    );
	}
	if (status == 6)
	{
	    throw Safe_Error_Return ("Convergence failure (Loglikelihood NaN)");
	}
	double quadratic;
	if (vtraits && status != 7) 
	{
	    quadratic = qdform/(nsamp*vtraits);
	} else {
	    quadratic = qdform/(nsamp*ntrait);
	}
	Loglike::quadratic (quadratic);
    }
    if (iarray[LENI-1] != TEST_NUMBER)
    {
	throw Safe_Error_Return (
     "IARRAY overrun occurred in search; contact solar@txbiomedgenetics.org");
    }
    if (rarray[Lenr-1] != TEST_NUMBER)
    {
	throw Safe_Error_Return (
        "RARRAY overrun occurred in search; contact solar@txbiomedgenetics.org");
    }
    } // End try
    catch (Safe_Error_Return& ser)
    {
	if (Strcmp (ser.message(), "InitPar only"))
	{
	    strncpy (error_message, ser.message(), error_message_len);
	    error_message_lnblnk = (int) strlen (ser.message ());
	    iquit = 1;
	}
	else
	{
	    loglik = 0.0; // maximize knows this is initpar, but must skip ll
	}
    }
    if (!iquit)
    {
/*
 * Update parameters for both normal and initpar cases
 */
    for (i = 0; i < Parameter::count(); i++)
    {
	double result = rarray[par_start+i];
	if (fabs(result) < MIN_VISIBLE)
	{
	    result = 0.0L;
	}
	double se = rarray[par_se+i];
	if (inasyc == 0 || fabs(se) < MIN_VISIBLE)
	{
	    se = 0;
	}
	if (DiscreteTraits)
	{
            double score = 0.0;
	    Parameter::set (i, result, rarray[m11+i], rarray[m10+i], se,
                            score);
	}
	else
	{
	    Parameter::set (i, result, rarray[m14+i], rarray[m13+i], se,
			    rarray[m20+NumIncPED+i*(NumIncPED+1)]);
	}
    }
    } // End not quitting (i.e. useable parameters)

/*
 * Free huge arrays
 */
    free (vmin);
    free (vmax);
    free (vmean);

    }  // End if no pedigree error

    free (parest);
    free (vnames);
    free (Vbinary);
    free (ntrasamp);
    free (ttrasamp);
    free (nsampa);

    free (rarray);
    free (iarray);
    free (carray);
    free (extra);

    if (error_message_lnblnk)
    {
	error_message[(error_message_lnblnk<255)?error_message_lnblnk:255] = 
	    '\0';
	throw Safe_Error_Return (error_message);
    }
    if (status == 7)
    {
	throw Safe_Error_Return (
	    "Convergence failure (Restartable) (Loglikelihood NaN)");
    }
    if (iquit)
    {
	if (matrix_errors)
	{
	    char buf[256];
	    sprintf (buf,"Matrix file is missing person IBDID=%d (more details above)",matrix_errors);
	    throw Safe_Error_Return (buf);
	}
	throw Safe_Error_Return (
         "Error encountered in data; see above (or output files) for details");
    }
	

#ifdef NOIEEE
    if (sigfpe_raised)
    {
	loglik = 0.0;
    }
    signal (SIGFPE,SIG_DFL);
#endif

    return loglik;
}

/*
 * READ CAREFULLY
 *
 * fortspad is a replacement for Spadf in "safelib.c"
 *   Spadf had a serious design bug (actually, if you
 *   followed the comment in Spadf, it was OK, but it is
 *   very natural NOT to follow the comment...
 *
 *   Spadf required the input string to be null terminated.
 *   That would NOT work if input string was the exact length of
 *   the buffer.  Furthermore, if the input buffer was ultimately
 *   destined to be an array of FORTRAN strings (very common) there
 *   is no nice way to have an extra null-terminating byte.
 *
 *   So, this has been re-written NOT to require NULL termination
 *   if the string is the exact length specified.  Otherwise, null
 *   termination is STILL REQUIRED.
 *
 *   Here is the description of fortspad (modifed slightly from Spadf):
 *
 * fortspad pads out a c string with blanks to the specified length for use 
 * by FORTRAN.  FORTRAN strings are not terminated...they are blank padded.
 * Note: the length had better be correct or a segment violation will occur!
 * The string should be null terminated (as for a C string) on input, but
 * will probably not be null terminated on output.  If the string is
 * the exact length specified, null termination is not necessary.
 */

void fortspad (char *string, int length)
{
    int slen;
    for (slen = 0; slen < length; slen++)
    {
	if (string[slen] == '\0') break;
    }
    while (slen < length)
    {
	string[slen++] = ' ';
    }
}

/*
 * Callbacks for FORTRAN
 *
 */

extern "C" void zsderror_ ()
{
    throw Safe_Error_Return ("Zscore option activated with zero SD");
}


extern "C" int isnan_f_ (double* fdouble)
{
    return isnan (*fdouble);
}

extern "C" void ifanydisc_ (int* ifany)
{
    *ifany = if_any_discrete;
}
