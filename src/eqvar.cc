/*
 * eqvar.cc implements the eqvar class
 * Written by Charles Peterson beginning on August 3, 2000
 * Copyright (c) 2000 Southwest Foundation for Biomedical Research
 */

// EqVar is used to ensure variables used in Mu and Omega equations
// are included in analysis, and to get their indexes during evaluation


#include "solar.h"

char *EqVar::name[2*MAX_PHENOTYPES] = {0};
int EqVar::packed_index[2*MAX_PHENOTYPES] = {0};
int EqVar::count = 0;
int EqVar::_For_Trait[2*MAX_PHENOTYPES];

int EqVar::bind (Tcl_Interp *interp)
{

// Clear out previous values (avoid memory leak)

    int i;

    if (count > 0)
    {
	for (i = 0; i < count; i++)
	{
	    if (name[i] != 0)
		free (name[i]);
	    name[i] = 0;
	}
	count = 0;
    }

// Clear out trait indexes; -1 indicates unused

    for (i = 0; i < 2*MAX_PHENOTYPES; i++)
    {
	_For_Trait[i] = -1;
    }

    return TCL_OK;
}

int EqVar::add_for_trait (const char* vname, int trait)
{
    add (vname);
    _For_Trait[count-1] = trait;
    return count-1;
}


int EqVar::add (const char* vname)
{
    if (count >= 2*MAX_PHENOTYPES)
	throw Safe_Error_Return 
	    ("Unexpected internal SOLAR error in EqVar::add; please report");
    
    name[count] = Strdup (vname);
    return count++;
}

const char* EqVar::index (int i)
{
    if (i >= count) return 0;

    return name[i];
}
int EqVar::for_trait (int i)
{
    if (i >= count) return -1;

    return _For_Trait[i];
}

void EqVar::set_packed_index (int i, int packed_ind)
{
    if (i >= count)
	;

    packed_index[i] = packed_ind;
}



