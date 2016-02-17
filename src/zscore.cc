/*
 * zscore.cc implements the zscore function operator, not the zscore
 * command (implemented in Tcl) which is obsolescent now.
 *
 * Written by Charles Peterson beginning on October 25, 2010
 * Copyright (c) 2010 Southwest Foundation for Biomedical Research
 */

// Note that much of the implementation of zscore is in the binding
// section of the trait and covariate commands.


#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "solar.h"
#include "safelib.h"



Zscore* Zscore::Zscores[] = {0};
int Zscore::ZscoreCnt = 0;

Zscore* Zscore::Get (int index)
{
    if (index < ZscoreCnt)
    {
	Safe_Error_Return ("Zscore get failed");
    }
    return Zscores[index];
}

int Zscore::Setup (char* parname, int _vindex)
{
    Zscore* zs = new Zscore;
    zs->pname = Strdup (parname);
    zs->vindex = _vindex;
    return zs->enable(false);
}

int Zscore::enable (bool validity)
{
    Zscores[ZscoreCnt] = this;
    this->meanvalid = validity;
    return ZscoreCnt++;
}

int Zscore::SetupMean(char* parname, double _mean, double _sd, int _vindex)
{
    Zscore* zs = new Zscore;
    zs->pname = Strdup (parname);
    zs->mean = _mean;
    zs->sd = _sd;
    zs->vindex = _vindex;
    return zs->enable(true);
}

void Zscore::Reset ()
{
    for (int i = 0; i < ZscoreCnt; i++)
    {
//	printf ("deleting zscore #%d\n",i);
	free (Zscores[i]->pname);
	free (Zscores[i]);
	Zscores[i] = 0;
    }
    ZscoreCnt = 0;
    
}

// If_Ready determines if all zscore info ready

bool Zscore::If_Ready ()
{
    for (int i = 0; i < ZscoreCnt; i++)
    {
	if (!Zscores[i]->meanvalid)
	{
	    return false;
	}
    }
    return true;
}


// zcheck determines if this phenotype needs zscore info

bool Zscore::Zneeded (char* pname)
{
// compare this name to all the names in the Zscore object

    for (int i = 0; i < ZscoreCnt; i++)
    {
	if (!Strcmp (pname, Zscores[i]->pname))
	{
	    if (!Zscores[i]->meanvalid)
	    {
		return true;
	    }
	    else
	    {
		break;
	    }
	}
    }
    return false;
}

// savemsd saves only the required mean,sd values in Tcl variables
// maximize then traps and restarts with this information available


char* Phenotypes::maxphen_index (int index)
{
    return Setup_Names->index(index);
}


extern "C" void savemsd_ (char* pname, double* mean, double* sd, int _pnsize)
{
    char cname[VNAMELEN+1];
    strncpy (cname, pname, VNAMELEN);
    cname[VNAMELEN] = '\0';
    int i;
    for (i = VNAMELEN-1; i >= 0; i--)
    {
	if (!isspace(cname[i])) break;
	cname[i] = '\0';
    }

// Check if this is one of the parameters that needs saving

    if (Zscore::Zneeded (cname))
    {
	char cbuf[1024];
	sprintf (cbuf, "zscorexp set %s %18.12e %12.12e", cname, *mean,
	    *sd);
	Tcl_Interp* interp = get_maximize_interp();
	int errstatus = Solar_Eval (interp, cbuf);
	if (errstatus)
	{
	    throw Safe_Error_Return ("Unable to save zscore data");
	}
    }
}

