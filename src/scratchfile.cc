/*
 * scratchfile.cc implementes the ScratchFile class with FORTRAN
 * interfaces RSCRAT and ISCRAT to replace the old temporary-file-based
 * storage subsystem used by Fisher.
 *
 * Written by Charles Peterson
 * Copyright (c) 2001 Southwest Foundation for Biomedical Research
 */

#include "solar.h"

ScratchFile* ScratchFile::ScratchFiles[] = {0};
ScratchFile* ScratchFile::Ptrs[] = {0};

/****************************************************************************\
 *                  FORTRAN interfaces                                      *
\****************************************************************************/

extern "C" void rscrat_ (double* stuff, int* len, int* unit, int* rite,
			 int ignoresize)
{
    if (*rite)
    {
	ScratchFile::Write (*unit, *len, Real, 0, stuff);
    }
    else
    {
	ScratchFile::Read (*unit, *len, Real, 0, stuff);
    }
}

extern "C" void iscrat_ (int* stuff, int* len, int* unit, int* rite)
{
    if (*rite)
    {
	ScratchFile::Write (*unit, *len, Int, stuff, 0);
    }
    else
    {
	ScratchFile::Read (*unit, *len, Int, stuff, 0);
    }
}

extern "C" void rescratch_ (int* unit)
{
    ScratchFile::Rewind (*unit);
}

extern "C" void resetallsc_ ()
{
    ScratchFile::Reset();
}

// ***************************************************************************
//                            Implementation
// ***************************************************************************

void ScratchFile::Rewind (int unit)
{
    CheckActive (unit);
    ScratchFile *current = ScratchFiles[unit];
    Ptrs[unit] = current;
    current->index = 0;
}

void ScratchFile::Reset ()
{
    for (int i = 0; i < MAX_SCRATCH_UNIT; i++)
    {
	if (ScratchFiles[i] != 0)
	{
	    delete ScratchFiles[i];
	    ScratchFiles[i] = 0;
	}
    }
}


// WriteReal opens "file" if it doesn't exist
// If it does exist, it begins writing immediately after last write or read
// However, previous data beyond current pointer is deleted

void ScratchFile::Write (int unit, int len, ScratchType stype,
				    int* idata, double* rdata)
{
#if 0
    fprintf (stderr, "ScratchFile write unit %d  type %d  len %d\n", unit,
	     stype, len);
#endif

    CheckUnit (unit);
    ScratchFile* current;
    if (!ScratchFiles[unit])
    {
	Ptrs[unit] = ScratchFiles[unit] = new ScratchFile;
	current = ScratchFiles[unit];
    }
    else
    {
	Ptrs[unit]->size = Ptrs[unit]->index;  // Terminate last block
	if (Ptrs[unit]->next)  // Delete old content after this point
	{
	    delete Ptrs[unit]->next;
	}
	current = Ptrs[unit] = Ptrs[unit]->next = new ScratchFile;
    }
    if (stype == Int)
    {
	current->iarray = (int *) malloc (len*sizeof(int));
    } else {
	current->rarray = (double *) malloc (len*sizeof(double));
    }
    int i;
    for (i = 0; i < len; i++)
    {
	if (stype == Int)
	{
	    current->iarray[i] = idata[i];
	} else {
	    current->rarray[i] = rdata[i];
	}
    }
    current->size = current->index = i;
}

// Note: You cannot read what you have not written, and it must be same type

void ScratchFile::Read (int unit, int len, ScratchType stype,
				   int* idata, double* rdata)
{
#if 0
    fprintf (stderr, "ScratchFile read unit %d  type %d  len %d\n", unit,
	     stype, len);
#endif

    CheckUnit (unit);
    CheckActive (unit);
    ScratchFile *current = Ptrs[unit];
    int i;
    for (i = 0; i < len; i++)
    {
	if (current->index >= current->size)
	{
	    if (0 == (current = Ptrs[unit] = current->next))
	    {
		fprintf (stderr, "Scratchfile unit %d read past data\n", unit);
		exit (10);
	    }
	    current->index = 0;
	}
	if ((stype == Int && current->iarray == 0) ||
	    (stype == Real && current->rarray == 0))
	{
	    fprintf (stderr, "Scratchfile unit %d data type error\n", unit);
	    exit (10);
	}

	if (stype == Int)
	{
	    idata[i] = current->iarray[current->index++];
	} else {
	    rdata[i] = current->rarray[current->index++];
	}
    }
}

    



