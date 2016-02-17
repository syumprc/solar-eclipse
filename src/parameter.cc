#include <string.h>
#include <stdlib.h>
#include "safelib.h"
// tcl.h from solar.h
#include "solar.h"

extern void trim_blank_sides (char* buffer); // from tablefile.cc

extern "C" void cdfnor_ (int* which, double* p, double* q, double* x, double* mean, double* sd, int* status, double* bound);

int Parameter::_count = 0;
Parameter *Parameter::Parameters[MAX_PARAMETERS] = {0};

Parameter *Parameter::find (const char *name)
{
    for (int i = 0; i < _count; i++)
    {
	if (!StringCmp (name, Parameter::Parameters[i]->name(), case_ins))
	{
	    return Parameter::Parameters[i];
	}
    }
    return NULL;
}

int Parameter::Find0 (const char *name)
{
    Parameter* p = find (name);
    if (!p) return -1;
    return p->number();
}


Parameter *Parameter::add (const char *name)  // public interface
{
    Parameter *p = new Parameter (name);
    p->_number = _count++;
    Parameters[p->_number] = p;
    return p;
}

// Public interface used to insert trait parameters at beginning
Parameter *Parameter::insert (const char* name)
{
    int i;

// If parameter already exists, just move it to beginning
    Parameter* existing = find (name);
    if (existing)
    {
	for (i = existing->_number; i > 0; i--)
	{
	    Parameters[i] = Parameters[i-1];
	    Parameters[i]->_number = i;
	}
	Parameters[0] = existing;
	existing->_number = 0;
	return existing;
    }

// Make room for new parameter
    for (i = _count; i > 0; i--)
    {
	Parameters[i] = Parameters[i-1];
	Parameters[i]->_number = i;
    }

// Add new parameter at beginning
    Parameter *p = new Parameter (name);
    _count++;
    p->_number = 0;
    Parameters[p->_number] = p;
    return p;
}
    
int Parameter::longest_name ()
{
    int longest = 0;
    Parameter *p;
    for (int i = 0; p = Parameter::index (i); i++)
    {
	if (strlen(p->name()) > longest) longest = strlen (p->name());
    }
    return longest;
}

void Parameter::rename (const char *newname)
{
    free (_name);
    _name = Strdup (newname);
}


Parameter::Parameter (const char *nm)
{
    _name = Strdup (nm);
    start = 0.;
    last = 0.;
    lower = 0.;
    upper = 0.;
    se = 0.;
    score = 0.;
    covariate = 0;
    mod_start = mod_upper = mod_lower = false;
    fix_upper = false;
    fix_lower = false;
}

Parameter::~Parameter ()
{
    free (_name);
}

void Parameter::reset()
{
    while (count())
    {
	Parameter *p = index (0);
	p->remove ();
    }
}


void Parameter::remove ()
{
    if (covariate)
    {
	throw Covariate_Parameter ();
    }
    Constraint::remove_if_needs_par (_name);
    for (int i = _number+1; i < _count; i++)
    {
	Parameter *p = Parameters[i];
	Parameters[i-1] = p;
	p->_number--;
    }
    _count--;
    delete this;
}

char *Parameter::to_string (char *buf)
{
    char sbuf[128];
    char lbuf[128];
    char ubuf[128];
    char pformat[64];
    sprintf (pformat,"%%-.%ig",Option::get_int("ParameterFormat"));

//    printf ("parameter format is %s\n",pformat);
//    sprintf (sbuf, "%-.16g", start);

    sprintf (sbuf, pformat, start);
    sprintf (lbuf, pformat, lower);
    sprintf (ubuf, pformat, upper);
    const char* fixls = "";
    const char* fixus = "";
    if (fix_lower) fixls = "fix";
    if (fix_upper) fixus = "fix";
    sprintf (buf, "%8s = %-13s %sLower %-10s  %sUpper %-10s",
	     name(), sbuf, fixls, lbuf, fixus, ubuf);
    return buf;
}

void Parameter::return_names(Tcl_Interp* interp)
{
    char* nbuf = 0;
    for (int i = 0; i < count(); i++)
    {
	if (i > 0) {
	    StringAppend (&nbuf, " ");
	}
	StringAppend (&nbuf, Parameters[i]->name());
    }
    RESULT_BUF (nbuf);
    if (nbuf) free (nbuf);
}


void Parameter::display_all()
{
    for (int i = 0; i < count(); i++)
    {
	char buf[1024];
	printf ("%s\n", Parameters[i]->to_string (buf));
	if (Parameters[i]->se != 0.L)
	{
	    printf ("%8s                    se %-.10g  score %-.10g\n",
		    Parameters[i]->name(), Parameters[i]->se,
		    Parameters[i]->score);
	}
    }
}

void Parameter::write_commands (FILE *file)
{
    for (int i = 0; i < count(); i++)
    {
	char buf[256];
	Parameters[i]->to_string (buf);
	fprintf (file, "parameter %s\n", buf);
	if (Parameters[i]->se != 0.L)
	{
   fprintf (file, "parameter %8s                    se %-.10g  score %-.10g\n",
	    Parameters[i]->name(), Parameters[i]->se, Parameters[i]->score);
	}
    }
}    

int Parameter::return_all (Tcl_Interp* interp)
{
    char* buf = Strdup ("");
    char tbuf[1024];
    
    for (int i = 0; i < count(); i++)
    {
	string__append (&buf,"{");
	Parameters[i]->to_string (tbuf);
	trim_blank_sides (tbuf);
	string__append (&buf,tbuf);
	string__append (&buf,"}\n");
	
	if (Parameters[i]->se != 0.L)
	{
	    string__append (&buf,"{");
	    sprintf (tbuf,"%s  se %-.10g  score %-.10g}\n",
		    Parameters[i]->name(), Parameters[i]->se,
		    Parameters[i]->score);
	    string__append (&buf,tbuf);
	}
    }
    RESULT_BUF(buf);
    free (buf);
    return TCL_OK;
}


extern "C" int ParameterCmd (ClientData clientData, Tcl_Interp *interp,
		  int argc, char *argv[])
{
    int i;
    bool insert = false;

    if (argc == 1)
    {
	Parameter::display_all ();
	return TCL_OK;
    }

    if (argc == 2 && !Strcmp (argv[1], "-return"))
    {
	return Parameter::return_all(interp);
    }

    if (argc == 2 && !Strcmp (argv[1], "-names"))
    {
	Parameter::return_names (interp);
	return TCL_OK;
    }


    if (argc >= 3 && !Strcmp (argv[1], "-insert"))
    {
	for (i = 2; i < argc; i++)
	{
	    argv[i-1] = argv[i];
	}
	argc--;
	insert = true;
    }


// Now we have 2 or more arguments

// Check for invalid assignment operation such as parameter foobar=1, must
// be foobar = 1.

// The way this now works is following the fairly simple rule "a parameter
// name must not include the = unless it is enclosed in ()".  This rule allows
// expression trait variance components such as h2r(=log(q4)), which
// obviously must have = enclosed in parentheses.  In versions of SOLAR prior
// to 4.0, you could not have = in parameters at all.

// The test is then, find the last =.  If it is not followed by ), then it
// is invalid.

// Standalone = is OK, the one except to above rule.    

    int badeq = false;

    for (i = 1; i < argc; i++)
    {
	if (strchr(argv[i],'=') && argv[i][1] != '\0')
	{
	    int length = strlen (argv[i]);
	    bool seenparen = false;
	    for (int index = length-1; index >= 0; index--)
	    {
		if (argv[i][index] == '=')
		{
		    if (!seenparen)
		    {
			badeq = true;
			break;
		    }
		}
		else if (argv[i][index] == ')')
		{
		    seenparen = true;
		}
	    }
	}
    }
		    
    if (badeq)
    {
	RESULT_LIT (\
"\nparameter: Invalid syntax, see help parameter.\n\
If assigning parameter start value, = must be surrounded by spaces.\n\
Inside a parameter name, = only allowed inside parentheses.\n");
	return TCL_ERROR;
    }

    if (argc == 2 && !StringCmp (argv[1], "help", case_ins))
    {
	return Solar_Eval (interp, "help parameter");
    }

    if (argc == 2 && (!StringCmp (argv[1], "delete_all", case_ins)))
    {
	RESULT_LIT (\
        "Parameter command does not support delete_all; consider 'model new'");
	return TCL_ERROR;
    }


    if (argc == 3 && (!Strcmp (argv[1], "delete") ||
                      !Strcmp (argv[1], "delete_fully") ||
	              !Strcmp (argv[2], "delete") ||
	              !Strcmp (argv[2], "delete_fully")))
    {
	char *pname = argv[2];
	if (!Strcmp (argv[2], "delete") ||
	    !Strcmp (argv[2], "delete_fully"))
	{
	    pname = argv[1];
	}

	Parameter *p = Parameter::find (pname);
	if (!p)
	{
	    char buf[256];
	    sprintf (buf, "There is no parameter '%s' to remove", pname);
	    RESULT_BUF (buf);
	    return TCL_ERROR;
	}
	try
	{
	    p->remove ();
	}
	catch (Covariate_Parameter)
	{
	    RESULT_LIT ("Must remove covariate instead");
	    return TCL_ERROR;
	}
	return TCL_OK;
    }
		     
    Parameter *p;
    p = Parameter::find (argv[1]);
    if (!p)
    {
	if (argc == 3)
	{
	    char buf[256];
	    sprintf (buf,
		     "No parameter named %s has been created", argv[1]);
	    RESULT_BUF (buf);
	    return TCL_ERROR;
	}
	if (insert)
	{
	    p = Parameter::insert (argv[1]);
	}
	else
	{
	    p = Parameter::add (argv[1]);
	}
    }
    else
    {
	if (argc == 2)
	{
	    char buf[256];
	    RESULT_BUF (p->to_string (buf));
	    return TCL_OK;
	}
    }
    if (argc == 3)  // value reading command
    {
	char *selector = argv[2];
	char buf[256];
	char pformat[64];
	sprintf (pformat,"%%-.%ig",Option::get_int("ParameterFormat"));

	if (!StringCmp (selector, "start", case_ins) ||
	    !StringCmp (selector, "=", case_ins))
	{
	    sprintf (buf, pformat, p->start);
	}
	else if (!StringCmp (selector, "lower", case_ins))
	{
	    sprintf (buf, pformat, p->lower);
	}
	else if (!StringCmp (selector, "upper", case_ins))
	{
	    sprintf (buf, pformat, p->upper);
	}
	else if (!StringCmp (selector, "last", case_ins))
	{
	    sprintf (buf, pformat, p->last);
	}
	else if (!StringCmp (selector, "fixupper", case_ins))
	{
	    if (p->fix_upper)
	    {
		sprintf (buf, pformat, p->upper);
	    } else {
		sprintf (buf,"");
	    }
	}
	else if (!StringCmp (selector, "fixlower", case_ins))
	{
	    if (p->fix_lower)
	    {
		sprintf (buf, pformat, p->lower);
	    } else {
		sprintf (buf,"");
	    }
	}
	else if (!StringCmp (selector, "se", case_ins))
	{
	    sprintf (buf, pformat, p->se);
	}
	else if (!StringCmp (selector, "score", case_ins))
	{
	    sprintf (buf, pformat, p->score);
	}
	else
	{
	    sprintf (buf, "Invalid argument %s", selector);
	    RESULT_BUF (buf);
	    return TCL_ERROR;
	}
	RESULT_BUF (buf);
	return TCL_OK;
    }

// Now, we expect an arbitrary number of "field value" pairs

    for (i = 2; i < argc;) // increment in loop
    {
	char *selector = argv[i++];
	if (i >= argc)
	{
	    char buf[256];
	    sprintf (buf, "Missing value field for %s", selector);
	    RESULT_BUF (buf);
	    return TCL_ERROR;
	}
	double value;
	char junk[128];
	int count = sscanf (argv[i++], "%lg %s", &value, &junk);
	if (count != 1)
	{
	    char buf[256];
	    sprintf (buf, "Invalid value specification '%s'", argv[i-1]);
	    RESULT_BUF (buf);
	    return TCL_ERROR;
	}
	if (!StringCmp (selector, "start", case_ins) ||
	    !StringCmp (selector, "=", case_ins))
	{
	    p->start = value;
	    p->mod_start = true;
	}
	else if (!StringCmp (selector, "lower", case_ins))
	{
	    p->lower = value;
	    p->mod_lower = true;
	    p->fix_lower = false;
	}
	else if (!StringCmp (selector, "fixlower", case_ins))
	{
	    p->lower = value;
	    p->mod_lower = true;
	    p->fix_lower = true;
	}
	else if (!StringCmp (selector, "upper", case_ins))
	{
	    p->upper = value;
	    p->mod_upper = true;
	    p->fix_upper = false;
	}
	else if (!StringCmp (selector, "fixupper", case_ins))
	{
	    p->upper = value;
	    p->mod_upper = true;
	    p->fix_upper = true;
	}
	else if (!StringCmp (selector, "last", case_ins))
	{
	    p->last = value;
	}
	else if (!StringCmp (selector, "se", case_ins))
	{
	    p->se = value;
	}
	else if (!StringCmp (selector, "score", case_ins))
	{
	    p->score = value;
	}
	else
	{
	    char buf[256];
	    sprintf (buf, "Invalid parameter field selector %s", selector);
	    RESULT_BUF (buf);
	    return TCL_ERROR;
	}
    }
    return TCL_OK;
}

int Parameter::bind (Tcl_Interp *interp)
{
    Parameter *p;
    for (int i = 0; p = Parameter::index (i); i++)
    {
	p->last = p->start;
    }
    return TCL_OK;
}

/*
 * FORTRAN interface (see fun.f and newlik.f)
 */
extern "C" int tdistp_ ()
{
    Parameter *p = Parameter::find ("t_param");
    if (!p)
    {
	throw Safe_Error_Return 
	    ("Can't find t_param; use tdist command to set up ");
    }
    return p->number()+1;
}


/*
 * FORTRAN interface (see inital.f)
 * It is assumed that caller knows the correct number of parameters;
 *   this is not checked.
 */
extern "C" void getparam_ (int *i, char *pname, double *par, double *min, 
			   double *max, int pnamelen_)
{
    int i0 = *i - 1;
    Parameter *p = Parameter::index (i0);
    *par = p->start;
    *min = p->lower;
    *max = p->upper;
    strncpy (pname, p->name(), pnamelen_);
    fortspad (pname, pnamelen_);
}

/*
 * initparam FORTRAN interface (see inital.f)
 *   set some parameters to initial values if defaulted
 *     this cannot be done at "bind" time since means and sd's are not known
 *   this applies to mean and sd parameters, and covariate betas
 */

extern "C" void initparam_ (double* vmean, double* vvar, double* vmin,
			   double* vmax)
{
    int t;
    bool initialized[MAX_PARAMETERS+1] = {false};   // Use 1-based indexing

    for (t = 0; t < Trait::Number_Of(); t++)
    {
	int pi;
	Parameter* p;
	if (-1 != (pi = Trait::Find_Parameter (t, "mean")))
	{
	    initialized[pi] = true;
	    p = Parameter::index(pi);

// Special handling for discrete trait mean parameter

	    if (p->start == 0 && Trait::Type(t) == discrete &&
		     Option::get_int ("EnableDiscrete"))
	    {
		int which = 2;
		double q = vmean[t] - vmin[t];
		if (q == 0.5)
		{
		    p->start = 0.0;
		}
		else
		{
		    double pct = 1 - q;
		    double mean = 0;
		    double sd = 1;
		    int status = 0;
		    double bound = 0;
		    double x = 0;
		    cdfnor_ (&which, &pct, &q, &x, &mean, &sd, &status, &bound);
		    if (status == 0)
		    {
			p->start = x;
		    }
		    else
		    {
			fprintf (stderr, "\n    *** Error computing inverse normal for %g\n", q);
			fprintf (stderr, "    *** Defaulting mean starting value to 0\n\n");
		    }
		}
		if (!p->mod_lower && !p->fix_lower) p->lower = -8.0;
		if (!p->mod_upper && !p->fix_upper) p->upper = 8.0;
		p->mod_start = p->mod_upper = p->mod_lower = true;
	    }
	    else
	    {

// Now for quantitative trait mean parameter

		if (p->start == 0.0 && p->lower == 0.0 && p->upper == 0.0)
		{
		    p->start = vmean[t];
		    p->lower = vmin[t];
		    p->upper = vmax[t];
		}
		if (!p->mod_start) p->start = vmean[t];
		if (!p->mod_lower && !p->fix_lower) p->lower = vmin[t];
		if (!p->mod_upper && !p->fix_upper) p->upper = vmax[t];
		p->mod_start = p->mod_upper = p->mod_lower = true;
	    }
	    Constraint::check_parameter (pi, &p->start, &p->lower, &p->upper);
	}
	if (-1 != (pi = Trait::Find_Parameter (t, "SD")))
	{
	    initialized[pi] = true;
	    p = Parameter::index (pi);
	    double sdstart = p->start;
	    if (p->start == 0.0 && p->lower == 0.0 && p->upper == 0.0)
	    {
		p->start = sqrt(vvar[t]);
		p->lower = 0.0;
		p->upper = p->start * 5.0;
	    }
	    if (!p->mod_start) p->start = sqrt (vvar[t]);
	    if (!p->mod_lower && !p->fix_lower) p->lower = 0.0;
	    if (!p->fix_upper && (!p->mod_upper || 
		(p->upper == 0.0 && p->upper <= sdstart) ||
				  (p->upper < p->start))) {

		p->upper = p->start * 5.0;
	    }
	    p->mod_start = p->mod_upper = p->mod_lower = true;
	    Constraint::check_parameter (pi, &p->start, &p->lower, &p->upper);
	}
    }

// Now check the boundaries for all other parameters
// If not already initialized,
//   Initialize all beta parameter (which begin with 'b') boundaries

    int i;
    for (i = 0; i < Parameter::count(); i++)
    {
	Parameter* p = Parameter::index(i);
	if (!initialized[i])
	{
	    if (p->lower == 0.0 && p->upper == 0.0)
	    {
		const char* pname = p->name();

// Boundaries for covariates handled by covariate module

		if (pname[0] == 'b' || pname[0] == 'B')
		{
		    int pif = i+1;
		    Covariate::set_boundaries (&pif, &p->start, &p->lower, 
					   &p->upper, vvar, vmin, vmax);

// Assumed range for rho parameters is -1 to 1

		} else if ((pname[0] == 'r' || pname[0] == 'R') &&
			  (pname[1] == 'h' || pname[1] == 'H') &&
			  (pname[2] == 'o' || pname[2] == 'O'))
		{
		    p->lower = -0.9;
		    p->upper = 0.9;

// Assumed range (when not specified) for all other values is 0 to 1

		} else {

// What falls through here uninitialized will probably generate
// a maximization error.  However, that's better than "guessing" what kind
// of parameter this is.

		}
	    }
	    initialized[i] = true;
	}

// Now, one last check to be sure parameter is within boundaries.  If not,
// move boundary if not fixed.  Otherwise, this will force error.

	if (!p->fix_lower && (p->lower > p->start))
	{
	    p->lower = p->start;
	}
	if (!p->fix_upper && (p->upper < p->start))
	{
	    p->upper = p->start;
	}
    }
}


/*
 * FORTRAN interface (see inital.f)
 */
extern "C" int iscov_ (int *i)
{
    Parameter *p = Parameter::index(*i-1);
    if (p)
    {
	if (p->covariate)
	{
	    return 1;
	}
    }
    return 0;
}


void format_double (char *string, double number, int width, int prec, 
		    double min)
{
/* Format a double number into a fixed width string as good as possible,
 *   Using exponential notation only if necessary (would overflow field
 *   otherwise, or be less than min--causing precision loss).
 *   Never suppress decimal point.
 *
 * Essentially, a 'FORTRAN' compatible but intelligent output.  What
 *   %g should have been, in my opinion.
 *
 * Unlike %g, constant width is maintained between the decimal and
 * exponential notation forms.  (Here, it's assumed to be more important
 * to keep columns straight than to maintain constant precision.)
 */
    char formatstr[64];
    char buffer[1024];  /* %f format can be >308 digits long! */
    int count;

    sprintf (formatstr, "%s%d.%d%s", "%", width, prec, "f");
    count = sprintf (buffer, formatstr, number);
    if (count <= 0 || count > width || (fabs(number) < min && number != 0.L))
    {
	if (number < 0.L)
	{
	    sprintf (formatstr, "%s.%d%s", "%", width-7, "E");
	    sprintf (buffer, formatstr, number);
	}
	else
	{
	    sprintf (formatstr, "%s.%d%s", "%", width-6, "E");
	    sprintf (buffer, formatstr, number);
	}
    }
    strcpy (string, buffer);

#if 0
    fprintf (stderr, "number=%g\n", number);
    fprintf (stderr, "formatstr=%s\n", formatstr);
    fprintf (stderr, "string=%s\n", string);
#endif
}

static int param_rcol;  // Right column for "Parameter"
extern "C" void writephead_ (int *unitno, int *type)
{
    char buf[1000];
    int longest_name = Parameter::longest_name ();
    int namelen = longest_name;
    if (namelen < 9) namelen = 9;   // allow for word "Parameter"

    int nspaces = (80 - (namelen + 39))/2;  // spaces before longest name
    param_rcol = nspaces + namelen;         // right column
    int spaces_pcol = param_rcol - strlen ("Parameter");
    const char *spaces = 
"                                                                            ";
    switch (*type)
    {
    case 1:
	strncpy (buf, spaces, spaces_pcol);  // Space over to P column
	buf[spaces_pcol] = '\0';

	strcat (buf, "Parameter      Start        Lower        Upper");
	writestr_ (unitno, buf, strlen(buf));
	break;
    case 2:
	strncpy (buf, spaces, spaces_pcol);  // Space over to P column
	buf[spaces_pcol] = '\0';

	strcat (buf, "Parameter    Final Val     Std Err       Const");
	writestr_ (unitno, buf, strlen(buf));

	strncpy (buf, spaces, spaces_pcol);  // Space over to P column
	buf[spaces_pcol] = '\0';

	strcat (buf, "               Start        Lower        Upper");
	writestr_ (unitno, buf, strlen(buf));
	break;
    }
}

/*
 * For either writeparam or write_results, 
 * writephead must be called first to set
 * param_rcol
 */

extern "C" void writeparam_ (int *unitno, int *pnum, double *startp, 
			     double *lowerp, double *upperp)
{
    char buf[1000];
    int pnumber = *pnum - 1;

    char lower[64];
    char upper[64];
    char start[64];

    format_double (start, *startp, 12, 7, 0.01L);
    format_double (lower, *lowerp, 12, 7, 0.01L);
    format_double (upper, *upperp, 12, 7, 0.01L);

    char format[128];
    char per = '%';
    sprintf (format, "%c%ds%c13s%c13s%c13s", per, param_rcol,
	     per, per, per);
    
    sprintf (buf, format,
	     Parameter::index (pnumber)->name(),
	     start, lower, upper);
    writestr_ (unitno, buf, strlen(buf));
}


// write names for correlation matrix
extern "C" void wcnames_ (int* unitno)
{
    const int formatlen = 12; // 1x,a11
    const int maxcount = 6;   // max 6 fields per line

    const int bufsiz = 1024;
    char line[bufsiz];
    char parameter[bufsiz];
    int total = 0;
    bool done = false;
    int index = 0;

    strcpy (line, "");
    while (!done) {
	for (int perline = 0; perline < maxcount; perline++)
	{
	    Parameter* p = Parameter::index(index++);
	    if (!p)
	    {
		done = true;
		break;
	    }
	    
	    int used = 1 + strlen(p->name());
	    
	    total += used;
	    if (total+formatlen >= bufsiz)  // protect against hacking
	    {
		done = true;
		break;
	    }
	    
	    strcat (line, " ");
	    strcat (line, p->name());
	    

	    for (; used < formatlen; used++)
	    {
		strcat (line, " ");
		total++;
	    }	
	}
	if (total > 0)
	{
	    writestr_ (unitno, line, strlen(line));
	    strcpy (line, "");
	}
	total = 0;
    }
}


extern "C" void writeresults_ (int *unitno, int *pnum, 
				double *finalp,	double *sterrp, char *cnstring,
				double *startp, double *lowerp, double *upperp,
				int cnstring_len_)
{
    char buf[1000];
    int pnumber = *pnum - 1;

    char final[64];
    char sterr[64];
    char cnstr[64];
    char start[64];
    char lower[64];
    char upper[64];

    int evdopts = Option::get_int ("EvdOptimizations");

    format_double (final, *finalp, 12, 7, 0.01L);
    format_double (start, *startp, 12, 7, 0.01L);
    format_double (lower, *lowerp, 12, 7, 0.01L);
    format_double (upper, *upperp, 12, 7, 0.01L);

    if (*sterrp == -1.0)
    {
	strcpy (sterr, " ");
    }
    else
    {
	format_double (sterr, *sterrp, 12, 7, 0.01L);
    }

    int clen = cnstring_len_;
    if (clen < 1)
    {
	strcpy (cnstr, " ");
    }
    else if (clen > 11)
    {
	strcpy (cnstr, " ");
	strncat (cnstr, cnstring, 8);
	cnstr[9] = '\0';
	strcat (cnstr,"...");
    }
    else
    {
	strcpy (cnstr, " ");
	strncat (cnstr, cnstring, clen);
	cnstr[clen+1] = '\0';
    }

    char format[128];
    char per = '%';
    
// First line has pname, final value, std err, and constraints

    char pname[1000];
    strncpy (pname,Parameter::index (pnumber)->name(),1000);

    if (evdopts)
    {
// Skip mean parameter, it's bogus for EVD2

	if (!Strcmp(pname,"mean") || !strncmp(pname,"mean(",5)) {
	    return;
	}
// Rename btmean to mean for EVD2
	if (!Strcmp(pname,"btmean") || !strncmp(pname,"btmean(",7))
	{
	    int ib = 0;
	    while (pname[ib] != '\0')
	    {
		pname[ib] = pname[ib+1];
		ib++;
	    }
	    ib = 0;
	    while (pname[ib] != '\0')
	    {
		pname[ib] = pname[ib+1];
		ib++;
	    }
	    char* epos = strstr(pname,"_evd");
	    if (epos)
	    {
		while (*(epos-1) != '\0')
		{
		    *epos++ = *(epos+4);
		}
	    }
	}
// Remove evd from beta parameter names
	if (pname[0] == 'b')
	{
	    char* epos = strstr(pname,"_evd");
	    if (epos)
	    {
		while (*(epos-1) != '\0')
		{
		    *epos++ = *(epos+4);
		}
	    }
	    epos = strstr(pname,"_evd");
	    if (epos)
	    {
		while (*(epos-1) != '\0')
		{
		    *epos++ = *(epos+4);
		}
	    }
	}
    }


    sprintf (format, "%c%ds %c12s %c12s %c-12s", per, param_rcol,
	     per, per, per);
    sprintf (buf, format,
	     pname,
	     final, sterr, cnstr);
    writestr_ (unitno, buf, strlen(buf));

// Second line has starting value, lower bound, and upper bound

    sprintf (format, "%c%ds %c12s %c12s %c12s", per, param_rcol,
	     per, per, per);
    sprintf (buf, format,
	     " ",
	     start, lower, upper);
    writestr_ (unitno, buf, strlen(buf));

// Third line is blank

    writestr_ (unitno, " ", 1);
}

// This section is for viewing/debugging tableau

#if 0

static int Noshow = 1;

extern "C" void starttab_ ()
{
    Noshow = 0;
}

extern "C" void stoptab_ ()
{
    Noshow = 1;
}


extern "C" void showtab_ (double table[64], int* ntab, int* compare)
{
    static double lasttable[10000];
    char buf[1024];

    if (Noshow) return;
    printf ("\f");
    for (int i = 0; i < *ntab; i++)
    {
	for (int j = 0; j < *ntab; j++)
	{
	    if (j < i)
	    {
		sprintf (buf, ".");
	    }
	    else
	    {
		if (*compare &&
		    (lasttable[i+j*(*ntab)] == table[i+j*(*ntab)]))
		{
		    sprintf (buf, "*");
		}
		else
		{
		    sprintf (buf, "%9.2e", table[i+j*(*ntab)]);
		    lasttable[i+j*(*ntab)] = table[i+j*(*ntab)];
		}
	    }
	    printf ("%9s ",buf);
	}
	printf ("\n");
    }
    printf ("\n\n");
}

#endif
