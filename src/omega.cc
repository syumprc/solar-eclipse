/*
 * omega.cc implements the Omega class and command
 * Written by Charles Peterson beginning on October 23, 1997
 * Copyright (c) 1997 Southwest Foundation for Biomedical Research
 */

#include <string.h>
#include <stdlib.h>
#include "solar.h"
#include "expression.h"

#define UNDEFINED_OMEGA "Use_polygenic_to_set_standard_model_parameterization"

int VirtualMultivar = 0;  // used elsewhere also

static int DumpCount = 0;
static int DumpOmega = 0;
static int Ntraits = 0;

Expression *Omega::first_expression = new Expression (UNDEFINED_OMEGA);

Expression *Omega::_expression = first_expression;


extern "C" int OmegaCmd (ClientData clientData, Tcl_Interp *interp,
		  int argc, char *argv[])
{
    if (argc == 2 && !StringCmp ("help", argv[1], case_ins))
    {
	return Solar_Eval (interp, "help omega");
    }
    if (argc == 1)
    {
//	printf ("omega = %s\n", Omega::expression->string ());
	Solar_AppendResult2 (interp, "omega = ", Omega::expression()->string (), 
			  0);
	return TCL_OK;
    }

    if (argc == 2 && !StringCmp ("reset", argv[1], case_ins))
    {
	Omega::reset ();
	return TCL_OK;
    }

    if (argc == 2 && !StringCmp ("command", argv[1], case_ins))
    {
	char buf[1000];
	sprintf (buf, "omega = %s", Omega::expression()->string ());
	RESULT_BUF (buf);
	return TCL_OK;
    }

    if (argc == 2 && !StringCmp ("debug", argv[1], case_ins))
    {
	char buf[10000];  // this might get pretty big

	Omega::expression()->debug (buf);

	RESULT_BUF (buf);
	return TCL_OK;
    }

    if (argc > 1 && argv[1][0] == '=')
    {
    // Set a new Omega equation
	char* eq_string = Tcl_Concat (argc-1, &argv[1]);
	char* expr_string;
	if (eq_string[1] == ' ')
	{
	    expr_string = &eq_string[2];
	} 
	else
	{
	    expr_string = &eq_string[1];
	}
	try
	{
	    Omega::new_expression (expr_string);
	}
	catch (Syntax_Error)
	{
	    RESULT_LIT ("Syntax Error");
	    Tcl_Free (eq_string);
	    return TCL_ERROR;
	}
	catch (Undefined_Function)
	{
	    RESULT_LIT ("Undefined Function");
	    Tcl_Free (eq_string);
	    return TCL_ERROR;
	}
	Tcl_Free (eq_string);
	return TCL_OK;
    }

    RESULT_LIT ("Invalid omega command");
    return TCL_ERROR;
}

Context *Omega::context = new Context;

// COVAR function calls omegac function to evaluate omega function

static int I = 0;
static int J = 0;
static double Phi2 = 0.0;
static double Delta7 = 0.0;

static double Itsd = 0.0;
static double* Vari = 0;
static double* Varj = 0;
static double* Par = 0;
static int Ibdid_Position = 0;
static int TraitI = 0;
static int TraitJ = 0;
static int H2r[MAX_TRAITS];
static int E2[MAX_TRAITS];
static int SD[MAX_TRAITS];
static int H2q[MAX_TRAITS][MAX_H2QINDEX];
static int H2qe[MAX_TRAITS][MAX_H2QINDEX];
static int C2[MAX_TRAITS];

// Allow for arbitrary bivariate parameterization

static int GENERIC[MAX_TRAITS][MAX_PARAMETERS];
static char* Gname[MAX_PARAMETERS];

// Allow for standard trivariate parameterization

static int Rhoe[MAX_TRAITS][MAX_TRAITS];
static int Rhog[MAX_TRAITS][MAX_TRAITS];

static int Rhoa[MAX_TRAITS][MAX_TRAITS];
static int Rhob[MAX_TRAITS][MAX_TRAITS];
static int Rhoc[MAX_TRAITS][MAX_TRAITS];
static int Rhod[MAX_TRAITS][MAX_TRAITS];
static int Rhof[MAX_TRAITS][MAX_TRAITS];

// static int Rhoq[MAX_TRAITS][MAX_TRAITS][MAX_H2QINDEX];
static int* Rhoq = 0;

static int* Maleip = 0;         // 99% of the time these aren't used
static int* Malejp = 0;         // So only de-reference when needed
static int* Itsd_index = 0;


// Fortran interface to obtain SD offset to SD parameters
// Input: 1-based trait index
// Returns: 1-based parameter index
 
extern "C" int sdpar_ (int *trait)
{
    return (1 +SD[*trait - 1]);
}


extern "C" double omegac_ (int *i, int *j, double *vari, double *varj,
			    double *par, double *phi2, double *delta7, 
			    int *itsd_index, int *ntrait, int *malei,
			    int *malej, int *traiti, int *traitj)
{
    I = *i;
    J = *j;
    Par = par;
    Phi2 = *phi2;
    Delta7 = *delta7;
    Vari = vari;
    Varj = varj;
    Itsd_index = itsd_index;
    Ibdid_Position = *ntrait;
    Maleip = malei;
    Malejp = malej;
    if (Trait::Number_Of() == 1)
    {	
	TraitI = 0;
	TraitJ = 0;
    }
    else if (VirtualMultivar)
    {
	Ibdid_Position = VirtualMultivar;
	TraitI = ((int) vari[0]) - 1;
	TraitJ = ((int) varj[0]) - 1;
    }
    else
    {
	TraitI = (int) (*traiti - 1);
	TraitJ = (int) (*traitj - 1);
    } 
    double thisomega = Omega::eval ();
    if (DumpOmega < 0 || DumpOmega > DumpCount)
    {
	if (DumpOmega > 0) DumpCount++;
	int ibdid = (int) Vari[Ibdid_Position];
	int jbdid = (int) Varj[Ibdid_Position];
	printf ("Omega[%d,%d]<%d,%d>{%9.6f,%9.6f} = %15.12f\n",
		ibdid,jbdid,TraitI+1,TraitJ+1,Phi2,Delta7,thisomega);
    }
//    int ibdid = (int) Vari[Ibdid_Position];
//    printf ("ibdid is %d\n",ibdid);

    return thisomega;
}


#ifdef MATRIX_DEBUG
static int matrix_count = 1;
FILE *MDfile = 0;
#endif

double user_matrix (int key)  // See also user_second_matrix
{
    Matrix *m = Matrix::index (key);
    int i = (int) Vari[Ibdid_Position];
    int j = (int) Varj[Ibdid_Position];
    float test = m->get (i, j);

#ifdef MATRIX_DEBUG
    if (!MDfile) {
	MDfile = fopen ("MatrixDebugFile","w");
    }
    fprintf (MDfile, "%d (%d,%d) = %f\n", matrix_count, i, j, (double) test);
#endif

    if (test != -1.0 || !m->defaultable) return test;
    Matrix *m2 = m->default_matrix;
    if (m2)
    {
//	fprintf (stderr, "Defaulting Matrix %s\n", m->name());
	return m2->get (i, j);
    }
//   fprintf (stderr, "Default scalar for Matrix %s\n", m->name());
    return *(m->default_scalar);
}

double user_second_matrix (int key)
{
    Matrix *m1 = Matrix::index (key);
    Matrix *m = m1->second_matrix;
    int i = (int) Vari[Ibdid_Position];
    int j = (int) Varj[Ibdid_Position];
    float test = m->get (i, j);



#ifdef MATRIX_DEBUG
    if (!MDfile) {
	MDfile = fopen ("MatrixDebugFile","w");
    }
    fprintf (MDfile, "%d (%d,%d) = %f\n", matrix_count, i, j, (double) test);
#endif

    if (test != -1.0 || !m->defaultable) return test;
    Matrix *m2 = m->default_matrix;
    if (m2)
    {
//	fprintf (stderr, "Defaulting 2nd Matrix %s\n", m->name());
	return m2->get (i, j);
    }
//    fprintf (stderr, "Default scalar for 2nd Matrix %s\n", m->name());
    return *(m->default_scalar);
}

static double pvar (int key)
{
    return Par[SD[0]] * Par[SD[0]];
}

static double identity_matrix (int key)
{
    int i = (int) Vari[Ibdid_Position];
    int j = (int) Varj[Ibdid_Position];
    if (i != j) return 0.0;
    return 1.0;
}

static double default_kinship_matrix (int key)
{
    return Phi2;
}

static double default_delta7 (int key)
{
    return Delta7;
}

static double default_itsd (int key)
{
    return ((*Itsd_index) ? Vari[*Itsd_index] : 0.0);
}

static double get_parameter_value (int pindex_0)
{
    return Par[pindex_0];
}

static double get_variable_mean (int vindex)
{
    int pi = EqVar::get_packed_index (vindex);
    return Mu::vmeans[pi];
}

static double get_variable_min (int vindex)
{
    int pi = EqVar::get_packed_index (vindex);
    return Mu::vmins[pi];
}

static double get_variable_max (int vindex)
{
    int pi = EqVar::get_packed_index (vindex);
    return Mu::vmaxs[pi];
}

static double get_variable_i (int vindex)
{
    int pi = EqVar::get_packed_index (vindex);
    return Vari[pi];
}

static double get_trait_i (int vindex)
{
    return (double) (TraitI + 1);
}

static double get_trait_j (int vindex)
{
    return (double) (TraitJ + 1);
}

static double get_variable_j (int vindex)
{
    int pi = EqVar::get_packed_index (vindex);
    return Varj[pi];
}

static double get_male_i (int vindex)
{
    return *Maleip;
}

static double get_male_j (int vindex)
{
    return *Malejp;
}

static double get_female_i (int vindex)
{
    return (1 - (*Maleip));
}

static double get_female_j (int vindex)
{
    return (1 - (*Malejp));
}

static double get_h2qe (int vindex)
{
    return fabs (Par[H2qe[TraitI][vindex]]);
}

static double get_c2 (int vindex)
{
    return fabs (Par[C2[TraitI]]);
}

static double get_h2q_ti (int vindex)
{
    return fabs(Par[H2q[TraitI][vindex]]);
}

static double get_h2q_tj (int vindex)
{
    return fabs(Par[H2q[TraitJ][vindex]]);
}

static double get_h2r_ti (int vindex)
{
    return fabs(Par[H2r[TraitI]]);
}

static double get_h2r_tj (int vindex)
{
//    fprintf (stderr, "j = %d, h2r(tj) = %.15f\n", TraitJ, Par[H2r[TraitJ]]);
    return fabs(Par[H2r[TraitJ]]);
}

static double get_e2_ti (int vindex)
{
//    fprintf (stderr, "i = %d, e2(ti) = %.15f\n", TraitI, Par[E2[TraitI]]);
    return fabs(Par[E2[TraitI]]);
}

static double get_e2_tj (int vindex)
{
    return fabs(Par[E2[TraitJ]]);
}

static double get_sd_ti (int vindex)
{
    return Par[SD[TraitI]];
}

static double get_sd_tj (int vindex)
{
    return Par[SD[TraitJ]];
}

static double get_generic_ti (int vindex)
{
    return Par[GENERIC[TraitI][vindex]];
}

static double get_generic_tj (int vindex)
{
    return Par[GENERIC[TraitJ][vindex]];
}

static double get_rhoe_ij (int vindex)
{
    return Par[Rhoe[TraitI][TraitJ]];
}

static double get_rhog_ij (int vindex)
{
    return Par[Rhog[TraitI][TraitJ]];
}

static double get_rhoq_ij (int vindex)
{
    return Par[Rhoq[TraitI + (TraitJ*Ntraits) + (vindex*Ntraits*Ntraits)]];
}

static double get_rhoc_ij (int vindex)
{
    return Par[Rhoc[TraitI][TraitJ]];
}

static double get_rhoa_ij (int vindex)
{
    return Par[Rhoa[TraitI][TraitJ]];
}

static double get_rhob_ij (int vindex)
{
    return Par[Rhob[TraitI][TraitJ]];
}

static double get_rhod_ij (int vindex)
{
    return Par[Rhod[TraitI][TraitJ]];
}

static double get_rhof_ij (int vindex)
{
    return Par[Rhof[TraitI][TraitJ]];
}


// teq returns 1 if same trait, 0 otherwise
static double get_teq (int vindex)
{
    if (TraitI == TraitJ)
    {
	return 1.0;
    }
    else
    {
	return 0.0;
    }
}

// tne returns 1 if different trait, 0 otherwise
static double get_tne (int vindex)
{
    if (TraitI != TraitJ)
    {
	return 1.0;
    }
    else
    {
	return 0.0;
    }
}


// this used to be part of bind, but now
// precedes bind to allow check before checking SD constraint

int Omega::check_if_defined (Tcl_Interp *interp)
{
    char* fstring = free_string();
    int comp = strcmp(fstring, UNDEFINED_OMEGA);
    free (fstring);

    if (!comp)
    {
	RESULT_LIT ("Undefined omega. Use polygenic or polymod for standard parameters.");
        return TCL_ERROR;
    }
    return TCL_OK;
}

    

int Omega::bind (Tcl_Interp *interp)
{
    int i = 0;
    char buf[1024];

    if (Rhoq) free (Rhoq);
    Rhoq = 0;

    DumpOmega = 0;
    DumpCount = 0;
    const char *dostring;
    Ntraits = Trait::Number_Of();

    char* ds = Strdup ("SOLAR_DumpOmega");
    if (0 != (dostring = Tcl_GetVar (interp, ds, TCL_GLOBAL_ONLY)))
    {
	DumpOmega = atoi (dostring);
    }
    free (ds);

// See if VirtualMultivar
    VirtualMultivar = 0;
    if (-1 != Option::get_int ("UnbalancedTraits") && Ntraits > 1)
    {
	VirtualMultivar = Ntraits;
    }

// First, build the Omega context
    if (context) delete context;
    context = new Context;

// The built-in stuff
    
    context->add ("I", identity_matrix);
    context->add ("Phi2", default_kinship_matrix);
    context->add ("Delta7", default_delta7);
    context->add ("Itsd", default_itsd);
    context->add ("Male_i", get_male_i, 0);
    context->add ("Male_j", get_male_j, 0);
    context->add ("Female_i", get_female_i, 0);
    context->add ("Female_j", get_female_j, 0);

// The parameters (may be overriden by subsequent additions to context
// because what is added last to context supercedes what came earlier.

    Parameter *p;
    for (i = 0; p = Parameter::index (i); i++)
    {
	context->add (p->name(), get_parameter_value, i);
    }

// If Univariate, add functions to retrive standard variance components
//  (New for 2.0.0; prevents them from going negative and unconverging)
// Previously, these would simply be handled by get_parameter_value
// For non-standard variance components (different parameterization) it
// is up to user to apply abs() function where required.
// This feature can be disabled by setting AbsVarianceParms option to zero.

    int qindex = 1;
    bool found_h2q = true;
    int pindex = -1;

// smpoutput uses pointer to SD parameters, if positive, so set negative
// in case they don't exist, then set them

    for (i = 0; i < MAX_TRAITS; i++)
    {
	SD[i] = -1;
	if (i < Ntraits)
	{
	    if (-1 != (pindex = Trait::Find_Parameter (i,"SD")))
	    {
		SD[i] = pindex;
	    }
	}
    }
	

    if (Ntraits == 1)
    {

// Force setting of SD pointer if SD exists so that pvar will work

      if (-1 != (pindex = Trait::Find_Parameter (0,"SD")))
      {
	  context->add ("pvar", pvar);
      }

      if (Option::get_int ("AbsVarianceParms") < 0)
      {
	if (-1 != (pindex = Trait::Find_Parameter (0,"h2r")))
	{
//	    fprintf (stderr, "Setting up abs function for h2r\n");
	    context->add ("h2r", get_h2r_ti, 0);
	    H2r[0] = pindex;
	}
	if (-1 != (pindex = Trait::Find_Parameter (0,"e2")))
	{
//	    fprintf (stderr, "Setting up abs function for e2\n");
	    context->add ("e2", get_e2_ti, 0);
	    E2[0] = pindex;
	}
	if (-1 != (pindex = Trait::Find_Parameter (0,"SD")))
	{
//	    fprintf (stderr, "Setting up abs function for SD\n");
	    context->add ("SD", get_sd_ti, 0);
	    SD[0] = pindex;
	}
		
	for (int qindex = 1; found_h2q; qindex++) //H2Q index
	{
	    char h2name[1024];
	    found_h2q = false;
	    sprintf (h2name, "h2q%d", qindex);
	    pindex = Parameter::Find0 (h2name);
	    if (pindex != -1)
	    {
		H2q[0][qindex-1] = pindex;
		context->add (h2name, get_h2q_ti, qindex-1);
		found_h2q = true;
	    }

// Check for epistatic components

	    sprintf (h2name, "h2qe%d", qindex);
	    pindex = Parameter::Find0 (h2name);
	    if (pindex != -1)
	    {
		H2qe[0][qindex-1] = pindex;
		context->add (h2name, get_h2qe, qindex-1);
		found_h2q = true;
	    }
	}

// Check for C2

	pindex = Parameter::Find0 ("c2");
	if (pindex != -1)
	{
	    C2[0] = pindex;
	    context->add ("c2", get_c2, 0);
	}
      }  // End If AbsVarianceParms option

/****************************************************************************/
/*                      Multivariate Omegas                                 */
/****************************************************************************/

    } else {

      int highest_h2q = 0;
      context->add ("teq", get_teq, 0);    // traits the same
      context->add ("tne", get_tne, 0);    // traits different

      if (-1==Parameter::Find0("si"))
      {
	  context->add ("si", get_trait_i, 0);
      }

      if (-1==Parameter::Find0("sj"))
      {	  
	  context->add ("sj", get_trait_j, 0);
      }

// Now add ti and tj pseudoparameters for all parameters having parens
// User will need to add "abs()" if negative protection is required
// If AbsVarianceParms option negative, THIS is how the parameters are handled
// otherwise it will be overriden for the standard variance component
// parameters in a following section

	Parameter *p;
	char gname[1024];
	char tname[1024];
	int gindex = 0;
	int j;
	for (i = 0; p = Parameter::index (i); i++)
	{
	    const char* cptr;
	    if (0 != strstr (p->name(),"(") )
	    {

// Get parameter name and trait number
// If not a trait-indexed parameter, continue

		strcpy (gname, p->name());
		char* ptr = strstr (gname, "(");
		*ptr = '\0';
		strcpy (tname, ptr+1);
		char* ptr2 = strstr (tname, ")");
		if (!ptr2) break;
		*ptr2 = '\0';
		int traitnumber = -1;
		for (j = 0; j < Ntraits; j++)
		{
		    if (!Strcmp (Trait::Name(j), tname))
		    {
			traitnumber = j;
			break;  // 
		    }
		}
		if (traitnumber < 0) continue;

// Does prefix match any previously seen?

		bool done = false;
		for (j = 0; j < gindex; j++)
		{
		    if (!Strcmp (gname, Gname[j]))
		    {

// yes, append to sublist for that prefix and continue

			GENERIC[traitnumber][j] = i;
			done = true;
			break;
		    }
		}
		if (done) continue;

// no, add to list of names seen

		Gname[gindex] = Strdup (gname);

// clear out parameter indexes so that missing parameter can be found

//		fprintf (stderr, "Clearing out indexes for %s\n", gname);

		for (j = 0; j < Ntraits; j++)
		{
		    GENERIC[j][gindex] = -1;
		}
		GENERIC[traitnumber][gindex] = i;
		gindex++;
	    }
	}

// See if all trait variants have been found for each putative generic name

	for (i = 0; i < gindex; i++)
	{
	    bool missing = false;
	    for (j = 0; j < Ntraits; j++)
	    {
		if (GENERIC[j][i] < 0) {
		    missing = true;
		}
	    }

// If all traits have been found, now we can add ti and tj to context

	    if (!missing)
	    {
		strcpy (gname, Gname[i]);
		strcat (gname, "(ti)");
		context->add (gname, get_generic_ti, i);
//		fprintf (stderr, "Adding %s for trait %d\n",
//			 gname, j);

		strcpy (gname, Gname[i]);
		strcat (gname, "(tj)");
		context->add (gname, get_generic_tj, i);
//		fprintf (stderr, "Adding %s for trait %d\n",
//			 gname, j);
	    }
	}

// free generic names

	for (i = 0; i < gindex; i++)
	{
	    free (Gname[i]);
	}

// If AbsVarianceParms nonzero (1 is default) add special functions to
// retreive trait dependent standard parameters (if present) for i or
// j...also ensure the standard parameters remain positive

// NOTE: unlike in previous section, we do not test that variance parameters
// are available for all traits if available for one...we can avoid that messy
// test *because* it has already been done.  The above section is mandatory
// while the following section is optional, even though in normal usage
// the following section *overrides* the preceding one.

      if (Option::get_int ("AbsVarianceParms"))
      {

	int found_sd = 0;
	int found_e2 = 0;
	int found_h2r = 0;

	for (i = 0; i < Ntraits; i++)
	{
	    if (-1 != (SD[i] = Trait::Find_Parameter (i,"SD"))) found_sd++;
	    if (-1 != (E2[i] = Trait::Find_Parameter (i,"e2"))) found_e2++;
	    if (-1 != (H2r[i] = Trait::Find_Parameter (i,"h2r"))) found_h2r++;
	}

	if (found_h2r == Ntraits)
	{
	    context->add ("h2r(ti)", get_h2r_ti, 0);
	    context->add ("h2r(tj)", get_h2r_tj, 0);
	}
	if (found_e2 == Ntraits)
	{
	    context->add ("e2(ti)", get_e2_ti, 0);
	    context->add ("e2(tj)", get_e2_tj, 0);
	}
	if (found_sd == Ntraits)
	{
	    context->add ("sd(ti)", get_sd_ti, 0);
	    context->add ("sd(tj)", get_sd_tj, 0);
	}

// H2q's

	found_h2q = true;
	for (qindex = 1; found_h2q; qindex++)  // H2Q index
	{
	    found_h2q = false;
	    for (i = 0; i < Ntraits; i++)  // Trait number
	    {
		char h2name[1024];
		sprintf (h2name, "h2q%d(%s)", qindex, Trait::Name(i));
		pindex = Parameter::Find0 (h2name);
		if (pindex == -1)
		{
		    continue;
		}
		H2q[i][qindex-1] = pindex;
		if (!found_h2q)
		{
		    sprintf (h2name, "h2q%d(ti)", qindex);
		    context->add (h2name, get_h2q_ti, qindex-1);
		    sprintf (h2name, "h2q%d(tj)", qindex);
		    context->add (h2name, get_h2q_tj, qindex-1);
		}
		found_h2q = true;
		highest_h2q = qindex;
		if (highest_h2q > MAX_H2QINDEX)
		{
		    RESULT_LIT ("Maximum H2Q Index of 100 EXCEEDED");
		    return TCL_ERROR;
		}
	    }
	}
      }

// Trivariate+ Rhos
//  If all the required rho parameters are found, the corresponding
//    rho _ij pseuedovariable is created.  If not, the
//    pseudovariable is not created, and the user might get
//    an undefined omega term.  This is intended to allow
//    for custom parameterization using fewer or no rhos.

    if (Ntraits > 2)
    {
	int i, j, k;
	char pname[256];
	char errbuf[1024];

	int rhos_expected = 0;
	int rhoes_found = 0;
	int rhogs_found = 0;
	int rhocs_found = 0;
	int rhoqs_found[MAX_H2QINDEX];

	int rhoas_found = 0;
	int rhobs_found = 0;
	int rhods_found = 0;
	int rhofs_found = 0;

	for (i = 0; i < MAX_H2QINDEX; i++)
	{
	    rhoqs_found[i] = 0;
	}

// Zero all possible rhos (no NaN's please)
//   Should I NaN them to catch errors?

	for (i = 0; i < Ntraits; i++)
	{
	    for (j = 0; j < Ntraits; j++)
	    {
		Rhoe[i][j] = 0;
		Rhog[i][j] = 0;
		Rhoc[i][j] = 0;
// Extra Rhos
		Rhoa[i][j] = 0;
		Rhob[i][j] = 0;
		Rhod[i][j] = 0;
		Rhof[i][j] = 0;
	    }
	}

	for (i = 0; i < Ntraits - 1; i++)
	{
	    for (j = i+1; j < Ntraits; j++)
	    {
		rhos_expected++;

		sprintf (pname, "rhoe_%d%d", i+1, j+1);
	        if (-1 != (pindex = Parameter::Find0 (pname)))
	        {
		    rhoes_found++;
	    	    Rhoe[i][j] = pindex;
	            Rhoe[j][i] = pindex;
	        }

		pname[0] = '\0';
		sprintf (pname, "rhog_%d%d", i+1, j+1);

	        if (-1 != (pindex = Parameter::Find0 (pname)))
	        {
		    rhogs_found++;
		    Rhog[i][j] = pindex;
		    Rhog[j][i] = pindex;
	        }

		pname[0] = '\0';
		sprintf (pname, "rhoc_%d%d", i+1, j+1);
	        if (-1 != (pindex = Parameter::Find0 (pname)))
	        {
		    rhocs_found++;
	    	    Rhoc[i][j] = pindex;
	            Rhoc[j][i] = pindex;
	        }
	
		pname[0] = '\0';
		sprintf (pname, "rhoa_%d%d", i+1, j+1);
	        if (-1 != (pindex = Parameter::Find0 (pname)))
	        {
		    rhoas_found++;
	    	    Rhoa[i][j] = pindex;
	            Rhoa[j][i] = pindex;
	        }
	
		pname[0] = '\0';
		sprintf (pname, "rhob_%d%d", i+1, j+1);
	        if (-1 != (pindex = Parameter::Find0 (pname)))
	        {
		    rhobs_found++;
	    	    Rhob[i][j] = pindex;
	            Rhob[j][i] = pindex;
	        }
	
		pname[0] = '\0';
		sprintf (pname, "rhod_%d%d", i+1, j+1);
	        if (-1 != (pindex = Parameter::Find0 (pname)))
	        {
		    rhods_found++;
	    	    Rhod[i][j] = pindex;
	            Rhod[j][i] = pindex;
	        }
	
		pname[0] = '\0';
		sprintf (pname, "rhof_%d%d", i+1, j+1);
	        if (-1 != (pindex = Parameter::Find0 (pname)))
	        {
		    rhofs_found++;
	    	    Rhof[i][j] = pindex;
	            Rhof[j][i] = pindex;
	        }
	
		for (qindex = 1; qindex <= highest_h2q; qindex++)
		{
		    pname[0] = '\0';
		    sprintf (pname, "rhoq%d_%d%d", qindex, i+1, j+1);
		    if (-1 != (pindex = Parameter::Find0 (pname)))
		    {
			if (!Rhoq)
			{
			    Rhoq = (int*) Calloc (Ntraits*Ntraits*highest_h2q,
					   sizeof(int));
			}
			rhoqs_found[qindex-1]++;
			Rhoq[i+(j*Ntraits)+((qindex-1)*Ntraits*Ntraits)] =
			    pindex;
			Rhoq[j+(i*Ntraits)+((qindex-1)*Ntraits*Ntraits)] =
			    pindex;
		    }
		}
	    }
	}
	if (rhoes_found == rhos_expected)
	{
	    context->add ("rhoe_ij", get_rhoe_ij, 0);
	}
	if (rhogs_found == rhos_expected)
	{
	    context->add ("rhog_ij", get_rhog_ij, 0);
	}
	if (rhocs_found == rhos_expected)
	{
	    context->add ("rhoc_ij", get_rhoc_ij, 0);
	}
	if (rhoas_found == rhos_expected)
	{
	    context->add ("rhoa_ij", get_rhoa_ij, 0);
	}
	if (rhobs_found == rhos_expected)
	{
	    context->add ("rhob_ij", get_rhob_ij, 0);
	}
	if (rhods_found == rhos_expected)
	{
	    context->add ("rhod_ij", get_rhod_ij, 0);
	}
	if (rhofs_found == rhos_expected)
	{
	    context->add ("rhof_ij", get_rhof_ij, 0);
	}
	for (qindex = 1; qindex <= highest_h2q; qindex++)
	{
	    if (rhoqs_found[qindex-1] == rhos_expected)
	    {
		sprintf (pname, "rhoq%d_ij", qindex);
		context->add (pname, get_rhoq_ij, qindex-1);
	    }
	}
    }
    }

// The user-provided matrices

    Matrix *m;
    Matrix *phi2_matrix = 0;
    Matrix *delta7_matrix = 0;
    for (i = 0; m = Matrix::index (i); i++)
    {
	context->add (m->name(), user_matrix, i);
	if (!Strcmp ("phi2", m->name())) phi2_matrix = m;
	if (m->second_matrix)
	{
//	    fprintf (stderr, "Adding 2nd matrix\n");
	    Matrix *m2 = m->second_matrix;
	    context->add (m2->name(), user_second_matrix, i);
	    if (!Strcmp ("delta7", m2->name())) delta7_matrix = m2;
	}
    }

// Now, find the defaults for defaultable matrices

    for (i = 0; m = Matrix::index (i); i++)
    {
	if (m->ibd())
	{
	    m->defaultable = true;
	    m->default_matrix = phi2_matrix;
	    m->default_scalar = &Phi2;
	}
	if (m->second_matrix)
	{
	    Matrix *m2 = m->second_matrix;
	    if (m2->d7())
	    {
		m2->defaultable = true;
		m2->default_matrix = delta7_matrix;
		m2->default_scalar = &Delta7;
	    }
	}
    }

// Second, bind the Omega expression to the Omega context
//   Add variables as required

    bool ok = false;
    while (!ok)
    {
	try
	{
	    expression()->bind (context);
	    ok = true;
	}
	catch (Undefined_Name& un)
	{
	    if (!Strcmp (un.name, "sex_i") ||
		!Strcmp (un.name, "sex_j"))
	    {
		RESULT_LIT (\
     "In place of Sex_i, use Male_i or Female_i in Omega equation");
		return TCL_ERROR;
	    }
	    if (!Strcmp (un.name, "H2r"))
	    {
		RESULT_LIT (\
     "H2r not defined; Use polygenic, polymod, spormod, or automodel command");
		return TCL_ERROR;
	    }
	    char test[5];
	    strncpy (test, un.name, 2);
	    test[2] = '\0';
	    if (!Strcmp (test, "x_"))
	    {
		char *vname = &un.name[2];
		if (!Phenotypes::available (vname))
		{
		    char buf[1024];
		    sprintf (buf, 
		      "Omega includes %s, but there is no %s phenotype",
			     un.name, vname);
		    RESULT_BUF (buf);
		    return TCL_ERROR;
		}
		int vindex = EqVar::add (vname);
		context->add (un.name, get_variable_mean, vindex);
		continue;
	    }
	    int unnamelen = strlen (un.name);
	    if (unnamelen > 2)
	    {
		strcpy (test,&un.name[strlen(un.name) - 2]);
		if (!Strcmp (test, "_i") || !Strcmp (test, "_j"))
		{
		    char *vname = Strdup (un.name);
		    vname[unnamelen - 2] = '\0';
		    if (!Phenotypes::available (vname))
		    {
			char buf[1024];
			sprintf (buf, 
			"Omega includes %s, but there is no %s phenotype",
				 un.name, vname);
			free (vname);
			RESULT_BUF (buf);
			return TCL_ERROR;
		    }
		    int vindex = EqVar::add (vname);
		    free (vname);
		    if (!Strcmp (test, "_i"))
		    {
			context->add (un.name, get_variable_i, vindex);
		    } else {
			context->add (un.name, get_variable_j, vindex);
		    }
		    continue;
		}
	    }
	    if (unnamelen > 4)
	    {
		strncpy (test, un.name, 4);
		test[4] = '\0';
		if (!Strcmp (test, "min_") || !Strcmp (test, "max_"))
		{
		    char *vname = &un.name[4];
		    if (!Phenotypes::available (vname))
		    {
			char buf[1024];
			sprintf (buf, 
			  "Omega includes %s, but there is no %s phenotype",
				 un.name, vname);
			RESULT_BUF (buf);
			return TCL_ERROR;
		    }
		    int vindex = EqVar::add (vname);
		    if (!Strcmp (test, "min_"))
		    {
			context->add (un.name, get_variable_min, vindex);
		    } else {
			context->add (un.name, get_variable_max, vindex);
		    }
		    continue;
		}
	    }
	    if (strstr (un.name, "(ti)"))
	    {
		char pname[256];
		char tname[512];
		bool found_missing = false;

		strncpy (pname, un.name, 256);
		char *ptr = strstr (pname, "(ti)");
		if (ptr)
		{
		    *ptr = '\0';
		    int i;
		    for (i = 0; i < Ntraits; i++)
		    {
			sprintf (tname, "%s(%s)", pname, Trait::Name(i));
			if (-1 == Parameter::Find0 (tname))
			{
			    sprintf (buf, "Missing parameter %s\n", tname);
			    found_missing = true;
			}
		    }
		}
		if (!found_missing)  // if name >250 chars !!!
		{
		    sprintf (buf, 
		     "Missing parameter %s for at least one trait\n", pname);
		}
	    }
	    else if (strstr (un.name, "_ij"))
	    {
		char pname[256];
		char tname[512];
		bool found_missing = false;

		strncpy (pname, un.name, 256);
		char *ptr = strstr (pname, "_ij");
		if (ptr)
		{
		    *ptr = '\0';
		    int i;
		    int j;
		    for (i = 1; i < Ntraits; i++)
		    {
			for (j = i+1; j <= Ntraits; j++)
			{
			    sprintf (tname, "%s_%d%d", pname, i, j);
			    if (-1 == Parameter::Find0 (tname))
			    {
				sprintf (buf, "Missing parameter %s\n", tname);
				found_missing = true;
			    }
			}
		    }
		}
		if (!found_missing)   // if name >250 chars !!!
		{
		    sprintf (buf,
                     "Missing parameter %s for at least one pair of traits\n",
			     pname);
		}
	    } else {
		sprintf (buf,
			 "Undefined variable `%s' in Omega equation", un.name);
	    }
	    RESULT_BUF (buf);
	    return TCL_ERROR;
	}
    }
    return TCL_OK;
}

void Omega::write_commands (FILE *file)
{
    fprintf (file, "omega = %s\n", Omega::expression()->string ());
}

void Omega::new_expression (const char *expr_string)
{
    Expression *last_expression = expression();
    _expression = new Expression (expr_string);  // Might generate exception
    if (last_expression != first_expression) delete last_expression;
}

void Omega::reset()
{
    if (expression() != first_expression) delete expression();
    _expression = first_expression;
}

    
