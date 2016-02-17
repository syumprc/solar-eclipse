/*
 * evdlikc.cc forms a C++ wrapper around the evdlik subroutine when called
 * from ddfun.f, allowing it to use dynamically allocated and/or
 * reloaded arrays.  It also forms a wrapper around tred2 and tql2, the
 * fortran routines which fill in the required matrices, and provides a
 * way to free the matrices.
 *
 * Written by Charles Peterson beginning on May 26, 2010
 * Copyright (c) 2010 Southwest Foundation for Biomedical Research
 * 
 */

#include <string.h>
#include <stdlib.h>
#include "solar.h"

static int Max_allocated_index = -1;
static double** Eval = 0;
static double** Evec = 0;
static int** ID = 0;
static int* N = 0;
static int Max_valid_index = -1;
static bool PM_valid = false;
static bool PM2_valid = false;

// Callthrough interfaces to Fortran

extern "C" void evdlik_ (int* n, int* max, double* mu, double* cov, 
			 double* loglike, double* h2, double* eval, 
			 double* evec);

extern "C" void tred2_ (int *max, int* n, double* tphi2, double* eval,
			double* evali, double* evec);

extern "C" void tql2_ (int* max, int* n, double* eval, double* evali,
		       double* evec, int* ierr);


// Internal interface
extern "C" void evdout (int* iped, double* vardata, int* male, int* vtraits,
			 int* nvar,
			 int* ntot,int* nind,int* nascer,int* nped,
			 int* ncumind, double* cov);

// Tcl interface to EVD::Flush

extern "C" int EVDCmd (ClientData clientData, Tcl_Interp *interp,
		  int argc, char *argv[])
{
    if (argc == 2 && !Strcmp (argv[1],"flush"))
    {
	EVD::Flush();
	return TCL_OK;
    }
    RESULT_LIT ("Invalid EVD command");
    return TCL_ERROR;
}


// evdlikc is the EVD routine still called by ddfun
// to calculate a likelihood for a pedigree using EVD.  If the
// EVD matrices have not yet been created for this pedigree,
// they are created now.  If this is EVDPhase 1, write output file
// instead of computing likelihoods

extern "C" void evdlikc_ (int* iped, int* maxpeo, int* n,
			  double* mu, double* cov, double* loglike, 
			  double* h2, double* vardata, int* nvar, 
			  int* evdphase, int* male, int* vtraits,
			  int* ntot, int* nind, int* nascer, int* nped,
			  int* ncumind, int* ierr)
{
    int zindex = *iped - 1;
    if (zindex > Max_valid_index)
    {
	EVD::Make_Matrix (zindex,maxpeo,n,mu,vardata,nvar,ierr);
	if (*ierr)
	{
	    fprintf (stderr, "tql2 returned ierr=%d for pedigree %d",
		     *ierr,*iped);
	    throw Safe_Error_Return ("EVD Failed in tql2");
	}
    }
    double *eval = Eval[zindex];
    double *evec =  Evec[zindex];

    if (*evdphase==0)
    {
	evdlik_ (n, maxpeo, mu, cov, loglike, h2, eval, evec);
    }
    else
    {
	evdout (iped,vardata,male,vtraits,nvar,ntot,nind,nascer,nped,ncumind,
		cov);
	EVD::Flush();  // EVD storage no longer needed
    }
}


extern "C" void evdout (int* iped, double* vardata, int* male, int* vtraits,
			 int* nvar,
			 int* ntot,int* nind,int* nascer,int* nped,
			 int* ncumind, double* cov)
{
//    for (int i = 0; i < *nvar * *ntot; i++)
//	printf ("vardata[%d] = %g\n",i,vardata[i]);
//    printf ("nvar is %d=n\n",*nvar);

    bool show_traits = false;

// Note: vardata is offset by ddfun to the first person in current pedigree
// male array is not offset

// This is called for each pedigree during EVDPhase 1.

//  printf ("Entering EVDOUT\n");

    int npheno = *nvar;
    int ni;
    int ifirstper, ivar, iind, jind;
    char numbuf[80];;
    double z;
    int ntraits = Trait::Number_Of();
    FILE* outfile;
    int ped = *iped - 1;

    if (ped!=0)
    {
	outfile = fopen ("evddata.out","a");
    }
    else
    {
	outfile = fopen ("evddata.out","w");

// Write header

	char* hstring = Strdup ("");
	string__append (&hstring, "ID,fa,mo,sex,sex_evd,tmean,lambda");
	int i;
	
	int nvariables = *nvar - 2;
	for (i = 0; i < *nvar; i++)
	{
// ID and group are always at ntraits and ntraits+1
// for both uni and multi variate within phenotype names list
// vardata reflects this, but actual multivariate trait values are within
// successive trait-individuals, with the remaining allocated slots for
// ntraits=3 and above simply unused.  But in this loop, all we care about
// is the names
	    if (i==ntraits || i==(ntraits+1)) continue;
	    const char* phenname = Phenotypes::get_var_name(i);
//	    printf ("phenname %d is %s\n",i,phenname);
	    if (show_traits)
	    {
		sprintf (numbuf, ",%s",phenname);
		string__append (&hstring,numbuf);
	    }
	    sprintf (numbuf, ",%s_evd",phenname);
	    string__append (&hstring,numbuf);
	}
	fprintf (outfile,"%s\n",hstring);
	free (hstring);
    }
    
// ncumind is spec'd by fortran as 0 based, but that means 0 is -1
    ifirstper = ncumind[ped]+1;
    ni = ncumind[ped+1]-ncumind[ped];
    int lastper = ncumind[ped+1];
    double* eigenvec = Evec[ped];
    double* eigenval = Eval[ped];

    int idindex = ntraits;
    int covindex = ntraits + 2;
//
// Now this gets complicated because we need to handle both univariate
// and multivariate traits.  For univariate traits, the trait values are
// seen as ivar==0, and we write complete data record when ivar==npheno+1.
// For multivariate traits, the trait values are ivar==1, and
// are sequenced over N individuals.
//
// For univariate, we build up the output record within a single iind and
// output at the last ivar (npheno+1) which is fake.
//
// For multivariate, we "build up" the output record over several
// iind's, and we don't write out the record until we see THE NEXT IBDID.
// So the fake ivar isn't needed, but a fake iind is needed instead (at
// the end of pedigree) to ensure all individuals get written.
//
    char* outstring=0;
    char* covstring=0;
    int last_ibdid = 0;
    int this_ibdid = 0;
    bool empty_record = true;
    bool need_tmean = true;
    bool ibdid_changed = false;

    int lastiind = ni;
    int lastivar = npheno;
    if (ntraits>1)
    {
	lastiind++;  // this is one-after iind to force output
	lastivar--;  // not using one-after ivar in this case
    }

    for (iind = 1; iind <= lastiind; iind++) // doesn't include probands
    {
//	printf ("iind is %d\n",iind);
// Loop through each variable, at end, write out

	double inpheno;
	int vindex;
	int rawsex;

	last_ibdid = this_ibdid;
	if (iind <= ni)
	{
	    int this_index = (*nvar)*(iind-1)+idindex;
	    this_ibdid = (int) vardata[this_index];
	    ibdid_changed = (last_ibdid != this_ibdid && last_ibdid!=0);
	} else {
	    ibdid_changed = true;
	}

	for (ivar = 0; ivar <= lastivar; ivar++)
	{
//	    if (outstring) printf ("     outstring is %s\n",outstring);
//	    if (covstring) printf ("     covstring is %s\n",covstring);


	    z = 0;

	    if ((ivar==npheno && ntraits==1) || \
		((ntraits>1) && ibdid_changed))
	    {
//		printf ("Writing to file\n");
// All traits for previous individual have been seen, now write last record
// for this pedigree

		if (covstring){
		    string__append (&outstring, covstring);
		    free (covstring);
		    covstring = 0;
		}
		fprintf (outfile, "%s\n", outstring);
		free (outstring);
		outstring = 0;
		empty_record = true;
		need_tmean = true;
		ibdid_changed = false;
		if (ivar==npheno || iind>ni)
		{
		    break;  // break from ivar loop
		}
	    }

	    if (empty_record)
	    {
// start new record
		outstring = Strdup ("");
		covstring = Strdup ("");
		rawsex = male[iind+ifirstper-2];
		int sex = 2 - rawsex;
		sprintf (numbuf, "%d,0,0,%d",this_ibdid,sex);
		string__append (&outstring, numbuf); // seq ID
		empty_record = false;
	    }

// We only need to process covariates once, so if this is
// the second record for same person (but different trait) just
// skip re-doing the covariates into covstring (re-use covstring)

// Pull out the actual phenotype from vardata

	    vindex = (*nvar)*(iind-1)+ ivar;
	    inpheno = vardata[vindex];
//	    printf ("  ivar is %d, inpheno is %g\n",ivar,inpheno);

// Now skip non-phenotypic records

	    if (ntraits == 1)
	    {
		if (ivar == 1 || ivar == 2)
		{
		    continue; // ID or group, skip
		}
	    } 
	    else
	    {
//
// For multivariate traits, the covariates phenotypes for subsequent
// person-traits (after the first) are identical,
// so after the first has been fully read, the subsequent ones can be
// skipped (and must be somehow).
//
// This needs to be hinged upon a variable that encodes that
// this is a subsequent person-trait of the same person.
//
// The problem with using a last_ibdid comparison is that it is also used
// for writing out entire record, then gets reset to avoid writing out the
// entire record again on the next phenotype.
//
// 
		if (ivar>=covindex && this_ibdid==last_ibdid && last_ibdid!=0)
		{
//		    printf ("   Skipping covariates for subsequent trait\n");
//		    printf ("     this_ibdid is %d,  last_ibdid is %d\n",
//			    this_ibdid, last_ibdid);
		    break;
		}
		if (ivar == 0)
		{
		    continue; 
		}
		if (1 < ivar && ivar < covindex)
		{
		    continue;
		}
	    }
//
// Sum eigenvector products for this trait
//
// The eigenvectors for multivariate are re-sized to the actual number of
//   persons, rather than person-traits.

//	    printf ("    Summing eigenvalue products\n");
	    z = 0;
	    double tmean = 0;
	    double lambda = 0;
	    for (jind = 0; jind < ni; jind=jind+ntraits)
	    {
		    
// in univariate or continuous balanced, first variable(s) are
// trait(s) followed by ID and group record...  but in unbalanced
// (what we actually use for bivariate) the first variable is the
// trait number instead followed by actual trait value, followed by
// ibdid and group record then followed by other variable data.
//		int eindex = jind*ni + iind-1; WRONG

		int realpeople = ni / ntraits;
		int real_i = (iind-1) / ntraits;
		int real_j = jind / ntraits;

//		int eindex = (iind-1)*ni + jind;

		int eindex = real_i*realpeople + real_j;

		double ine = eigenvec[eindex];
		lambda = eigenval[real_i];
//		printf ("i,j = %d,%d\n",iind,jind);
//		printf ("  ri,rj = %d,%d\n",real_i,real_j);
//		printf ("  eindex=%d   eigenvec=%g\n",eindex,ine);
		z += inpheno*ine;
		if (need_tmean)
		{
		    tmean += ine;
		}
	    }

// Write actual trait or covariate value
	    if (need_tmean)
	    {
		double sex_evd;
		if (rawsex)
		{
		    sex_evd = rawsex * tmean;
		}
		else
		{
		    sex_evd = 0;
		}
		sprintf (numbuf, ",%018.12e",sex_evd);
		string__append (&outstring, numbuf);
		sprintf (numbuf, ",%018.12e", tmean);
		string__append (&outstring, numbuf);
		sprintf (numbuf,",%018.12e", lambda);
		string__append (&outstring, numbuf);
		need_tmean = false;
	    }
	    if (show_traits)
	    {
		if (ivar<covindex)
		{
		    sprintf (numbuf,",%018.12e", inpheno);
		    string__append (&outstring, numbuf);
		} else {
		    sprintf (numbuf,",%018.12e", inpheno);
		    string__append (&covstring, numbuf);
		}
	    }

// Write EVD sums

	    if (ivar<covindex)
	    {
		sprintf (numbuf,",%018.12e", z);
		string__append (&outstring, numbuf);
	    } else {
		sprintf (numbuf,",%018.12e", z);
		string__append (&covstring, numbuf);
	    }
	}
    }
    fclose (outfile);
}

extern "C" void evdtrap_ (int *maxibdid)
{
    char buf[126];
    sprintf (buf,"Trap EVD Phase 2  maxibdid is %d",*maxibdid);
    throw Safe_Error_Return (buf);
}



// Reset does the simple reset at the beginning of each maximize
// This forces a check during maximization if each pedigree is still valid
// It DOES NOT flush EVD storage.

void EVD::Reset ()
{
    bool premax = verbose ("PREMAX");
    if (premax) printf ("**** EVD reset\n");
    Max_valid_index = -1;
    PM_valid = false;
    PM2_valid = false;
}

void EVD::Flush()
{
    for (int i = 0; i <= Max_allocated_index; i++)
    {
	if (Eval && Eval[i]) free (Eval[i]);
	if (Evec && Evec[i]) free (Evec[i]);
	if (ID && ID[i]) free (ID[i]);
    }
    if (Eval) {free (Eval); Eval=0;}
    if (Evec) {free (Evec); Evec=0;}
    if (ID) {free (ID); ID=0;}
    if (N) {free (N); N=0;}
    Max_allocated_index = -1;
    EVD::Reset();
}

void EVD::Make_Matrix (int zindex, int* maxpeo, int* n,
                       double* mu, double* vardata,
		       int* nvar, int* ierr)
{
    bool premax  = verbose ("PREMAX");  // 0x80000
    int maxpeon = *n;

// get offset to ibdid position within each vardata record
// univariate: trait,ibdid,proband/famid,
// multivariate: traitno,trait,ibdid,proband/famid

    int ntraits = Trait::Number_Of();
    int idpos = ntraits;
    if (ntraits>1)
    {
	maxpeon = maxpeon / ntraits;
    }

    bool samplesame = (1==Option::get_int ("SampleSameTrustMe"));
    bool needid = true;

// SampleSameTrustMe eliminates all sample testing...not recommended anymore
// Sample testing is now based on ID test, this is N and not N*N comparison as
//   previous phi2 test was.
// However, now there is need to potentially raise error in case sample changes
// So there is new option for that, dontallowsamplechange

    bool dontallows = (1==Option::get_int ("DontAllowSampleChange"));

    if (premax) printf ("**** Max valid index is %d\n",Max_valid_index);

    Matrix *pm = Matrix::find("phi2");
    if (!pm) {
	error ("Didn't get phi2 matrix");
    }

// Make root matrices if not present
    if (!Eval)
    {
	Eval = (double**) Calloc (1, sizeof(double*));
	Evec = (double**) Calloc (1, sizeof(double*));
	if (needid) ID = (int**) Calloc (1, sizeof(int*));
	N = (int*) Calloc(1,sizeof(int));
    }

// Reallocate root matrices to required size
    bool new_matrix_needed = false;
    if (zindex > Max_allocated_index)
   {
       if (premax) printf ("**** Reallocating EVD matrices to zindex %d\n");
	new_matrix_needed = true;
	Eval = (double**) 
	    realloc (Eval, (zindex+1)*sizeof(double*));
	Evec = (double**) 
	    realloc (Evec, (zindex+1)*sizeof(double*));
	if (needid)
	    ID = (int**) 
		realloc (ID, (zindex+1)*sizeof(int*));
	N = (int*) realloc (N, (zindex+1)*sizeof(int));
	for (int i = Max_allocated_index+1; i <= zindex; i++)
	{
	    Eval[i] = 0;
	    Evec[i] = 0;
	    if (needid) ID[i] = 0;
	    N[i] = 0;
	}
    }
    else
    {

// We have already allocated matrices for this pedigree index
// First check if matrix if valid (created since the beginning of current
//   maximization)
// Next check if existing matrices are still valid or need to be change
// First check size, if size is different, no need to check further
// If sizes match, check all ID's

	if (Max_valid_index < zindex)
	{
	    if (premax) printf ("Testing pedigree %d\n",zindex);
	    if (N[zindex] != maxpeon)
	    {
		if (premax) printf ("Matrix size wrong for %d\n",zindex);
		new_matrix_needed = true;
		if (samplesame)
		{
		    throw Safe_Error_Return 
		       ("Matrix size changed so option SampleSame is invalid");
		}
	    }
	    else
	    {

// Matrix comparison should only be done up to nxn most likely
// but how does Fortran store them ???

		if (!samplesame)
		{
		    if (premax) 
			printf ("Matrix size correct, now testing IBDIDs\n", 
				zindex);
		    for (int i = 0; i < maxpeon; i++)
		    {
			int idindex = idpos + ((*nvar)*i);
			int ibdid = (int) vardata[idindex];
//		    printf ("IBDID[%d] = %d\n",idindex,ibdid);
			int lastid = ID[zindex][i];
			if (ibdid != lastid)
			{
			    if (premax)
			    {
				printf ("**** ID DIFFERS in ped %d\n",zindex+1);
				printf ("**** %d, %d\n",ibdid,lastid);	
			    }
			    new_matrix_needed = true;
			    break;
			}
		    }
		}
	    }
	    if (new_matrix_needed)
	    {
		if (dontallows)
		{
		    throw Safe_Error_Return 
			("Sample changed: ID's differ");
		}
// Free previous matrix, data no longer valid
		if (Eval[zindex]) free (Eval[zindex]);
		if (Evec[zindex]) free (Evec[zindex]);
		if (ID[zindex]) free (ID[zindex]);
	    }
	}
    }
    
    if (new_matrix_needed)
    {
// Allocate branch matrices for this pedigree
	if (premax) printf ("EVD: Allocating new matrix for ped %d\n",zindex+1);
	if (premax) printf ("***** size is %d\n",maxpeon);
	Eval[zindex] = (double*) Calloc (maxpeon,sizeof(double));
	Evec[zindex] = (double*) Calloc ((maxpeon)*(maxpeon),sizeof(double));
//	printf ("Evec allocated to %d\n",maxpeon*maxpeon);
	if (needid) ID[zindex] = (int*) 
			  Calloc ((maxpeon),sizeof(int));
	double* eval = Eval[zindex];
	double* evec = Evec[zindex];
	double* tphi2 = (double*) Calloc ((maxpeon)*(maxpeon),sizeof(double));

// Save the Phi2, N, and Max for later comparisons
// The original design saved phi2 values for comparison next time around.
// That design is adequate for earlier EVD approaches in which phenotypes
// were used separately.  Otherwise, it has problem that sample could change,
// say if one singleton pedigree is replaced by another.  If matrices include
// results derived from sample data, they would need to change, even if the
// genetic relationships are unchanged.
//
// So it has changed to now save IBDID's (ID).
//
// We still need to gather phi2 data for use in EVD calculations.  But
// we do not need to save it.  Instead, the old variable "needphi2" is
// now "needids" and a simple linear list of ID's (ibdid's) in
// each pedigree.  To protect against new ibdid's, a load pedigree op
// throws out saved id's.
// 
// Simultaneously, we now have greater need to detect sample changes.
// In fact, we may go down this path partly to validate the sample.


	N[zindex] = maxpeon;
	int reported = 0;
	if (premax) printf (" Creating matrix for one pedigree\n");

// For bi/multivariate traits:
//   subsequent real persons are nt apart
//   fully balanced assumption for evd
//   maxpeon is the size of arrays, or real people, not the "person-traits"


// i and j are real people, adjust indexes to match

	for (int i = 0; i < maxpeon; i++)
	{
	    if (needid)
	    {
		int idindex = idpos + ((*nvar)*i*ntraits);
//		printf ("idindex is %d\n",idindex);
		ID[zindex][i] = (int) vardata[idindex];
	    }
	    for (int j = i; j < maxpeon; j++)
	    {
		int vindex = (int) vardata[idpos+(*nvar)*i*ntraits];
		int vvindex = (int) vardata[idpos+(*nvar)*j*ntraits];

		double value = pm->get(vindex,vvindex);

		tphi2[j+(maxpeon)*i] = value;
		tphi2[i+(maxpeon)*j] = value;
//		printf ("phi2[%d,%d] = %g\n",vindex,vvindex,value);
	    }
	}

// Do the matrix decompositions

	double* evali = (double*) Calloc (maxpeon, sizeof(double));

	tred2_ (&maxpeon,&maxpeon,tphi2,eval,evali,evec);
	tql2_ (&maxpeon,&maxpeon,eval,evali,evec,ierr);
	if (*ierr)
	{
	    error ("error in tql2\n");
	}
	free (evali);
	free (tphi2);
    } else {
	if (premax) printf ("New matrix not needed.\n");
    }
    if (Max_valid_index < zindex)
    {
	Max_valid_index = zindex;
    }
    if (Max_allocated_index < zindex)
    {
	Max_allocated_index = zindex;
    }
}

// Fortran interface to phi2 for special purpose
// Get self phi2 for individual in current vardata (pedigree)

extern "C" void getphi2i_ (int* iind, double* vardata, int* nvar,
			   double* phi2)
{
    static Matrix* pm;
    if (!PM_valid)
    {
	pm = Matrix::find("phi2");
	if (!pm) {
	    error ("Didn't get phi2 matrix");
	}
	PM_valid = true;
    }
    int ntraits = Trait::Number_Of();
    int idoff = ntraits;
    int index = idoff +(*iind-1)*(*nvar);
    int vindex = (int) vardata[index];
//    for (int i = 0; i < *nvar; i++) {
//	printf ("vardata[%d] = %g\n",i,vardata[i]);
//    }
    *phi2 = pm->get(vindex,vindex);
//  printf ("index = %d; phi2 = %g\n",vindex,*phi2);
}

extern "C" void ibdid2phi_ (int* ibdid1, int* ibdid2, double *phi2)
{
    static Matrix* pm2;
    if (!PM2_valid)
    {
	pm2 = Matrix::find("phi2");
	if (!pm2) {
	    error ("Didn't get phi2 matrix");
	}
	PM2_valid = true;
    }
    *phi2 = pm2->get(*ibdid1, *ibdid2);
}

extern "C" void unitclof_ (int* unit);  // fortran close statement

// option=0 means close if open
// option=1 means flag this now open
// CURRENTLY ONLY SUPPORTS UNIT 26 !!!!!

extern "C" void unitactive_ (int* option, int* unit)
{
    static bool unitopen = false;

    if (*option == 1)
    {
	unitopen = true;
    }
    else
    {
	if (unitopen)
	{
	    unitclof_ (unit);
	    unitopen = false;
	}
    }
}

extern "C" void funitactive (int option, int unit)
{
    unitactive_ (&option, &unit);
}

