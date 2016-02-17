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


class EVD
{
public:
    static double* Eval;
    static double* Evec;
    static int Max_ped_index;
    static Free();
    static Make_Matrix (int zindex, int* n, double *tphi2,
                       double* mu, double* cov, int* ierr);
};

// Fortran interfaces

extern "C" void evdlik_ (int* n, int* max, double* mu, double* cov, 
			 double* loglike, double* h2, double* eval, 
			 double* evec);

extern "C" void tred2_ (int *max, int* n, double* tphi2, double* mu,
			double* cov, double* eval, double* evali,
                        double* evec);

extern "C" void tql2_ (int* max, int* n, double* eval, double* evali,
		       double* evec, int* ierr);

extern "C" void evdlikc_ (int* iped, int* n, double* tphi2, 
			  double* mu, double* cov, double* loglike, 
			  double* h2, int* ierr)
{
    zindex = *iped - 1;
    if (zindex > EVD::Max_ped_index || !EVD:Evec[zindex])
    {
	EVD::make_matrix (zindex,n,tphi2,mu,cov,ierr);
    }
    if (*ierr)
    {
	fprintf (stderr, "tql2 returned ierr=%d for pedigree %d",*ierr,*iped);
	throw Safe_Error_Return ("EVD Failed in tql2");
    }
    double *eval = EVD::Eval[zindex];
    double *evec =  EVD::Evec[zindex];


// Fortran needs separate n and max variables because arrays are allocated
// to max for all peds in main subroutine, then used for n elements in
// each ped.  However, we allocate each array to precisely required size,
// so use n for both

    evdlik_ (n, n, mu, cov, loglike, h2, eval, evec);
}

void EVD::free ()
{
    for (int i = 0; i < Max_ped_index; i++)
    {
	if (Eval[i]) free (Eval[i]);
	if (Evec[i]) free (Evec[i]);
    }
    Max_ped_index = 0;
    free (Eval);
    free (Evec);
    Eval = 0;
    Evec = 0;
}

void EVD::make_matrix (int zindex, int* n, double *tphi2,
                       double* mu, double* cov, int* ierr);
{
    if (!EVD::Eval)
    {
	EVD::Eval = Calloc (1, sizeof(double*));
	EVD::Evec = Calloc (1, sizeof(double*));
    }
    if (zindex > EVD::Max_ped_index)
    {
	EVD::Eval = realloc (EVD::Eval, (zindex+1)*sizeof(double*));
	EVD::Evec = realloc (EVD::Evec, (zindex+1)*sizeof(double*));
    }
    EVD::Eval[zindex] = Calloc (*n,*n * sizeof(double));
    EVD::Evec[zindex] = Calloc ((*n)*(*n),sizeof(double));

    double* eval = EVD::Eval[zindex];
    double* evec = EVD::Evec[zindex];
    double* evali[n];

    tred2_ (n,n,tphi2,mu,cov,eval,evali,evec);
    tql2_ (n,n,eval,evali,evec,ierr);
}
    

    
