/*
 * voxel implements voxel and mask commands
 * using RicVolumeSet
 * Written by Charles Peterson beginning on January 8, 2013
 * Copyright (c) 2013 Texas Biomedical Research Institute
 */

// Implements commands: voxel, mask, ricvolumeset
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "solar.h"
// tcl.h from solar.h
#include "safelib.h"
#include "RicVolumeSet.h"


bool Voxel::_Valid = false;
int Voxel::X;
int Voxel::Y;
int Voxel::Z;

RicVolumeSet* GlobalBin = 0;
char* GlobalBinFilename = 0;
class Imout {
public:
    static int Show(Tcl_Interp* interp);
    static int Close();
    static RicVolumeSet* File;
    static char* Filename;
    static bool MaxValid;
    static int MaxX;
    static int MaxY;
    static int MaxZ;
    static bool DeltaValid;
    static float DeltaX;
    static float DeltaY;
    static float DeltaZ;
    static bool NvolValid;
    static int Nvol;
    static bool Valid;
    static bool Rewrite;
    static char* Orientation;  // this is used as Valid indicator
};

RicVolumeSet* Imout::File = 0;
char* Imout::Filename = 0;
bool Imout::MaxValid = false;
bool Imout::DeltaValid = false;
bool Imout::NvolValid = false;
bool Imout::Rewrite = false;
char* Imout::Orientation = 0;

int Imout::MaxX;
int Imout::MaxY;
int Imout::MaxZ;
float Imout::DeltaX;
float Imout::DeltaY;
float Imout::DeltaZ;
int Imout::Nvol;

class Mask {
public:
    static int Move (Tcl_Interp* interp, int user_offset);
    static bool FileValid;
    static bool _Valid;
    static int Intensity;
    static int Index;
    static int Gt;
    static int Ge;
    static int Lt;
    static int Le;
    static RicVolumeSet* File;
    static char* Filename;
    static int MaxX;
    static int MaxY;
    static int MaxZ;
    static float DeltaX;
    static float DeltaY;
    static float DeltaZ;
};

bool Mask::FileValid = false;
bool Mask::_Valid = false;

RicVolumeSet* Mask::File = 0;
char* Mask::Filename = 0;
int Mask::MaxX;
int Mask::MaxY;
int Mask::MaxZ;
int Mask::Intensity;
int Mask::Gt;
int Mask::Ge;
int Mask::Lt;
int Mask::Le;
int Mask::Index;
float Mask::DeltaX;
float Mask::DeltaY;
float Mask::DeltaZ;

int volumes_required (int ntrait, int ncovar)
{
    if (ntrait > 1)
    {
	printf ("Volumes not assigned for bivariate models");
    }
    return 20 + (4 * ncovar);
}
    

int Imout::Show(Tcl_Interp* interp)
{
    char buf[10000];
    buf[0] = '\0';
    if (Imout::Rewrite) {
	sprintf (buf,"imout -r %s",Imout::Filename);
    } else {
	sprintf (buf,"imout %s",Imout::Filename);
    }	    
    if (Imout::MaxValid) {
	sprintf (&buf[strlen(buf)]," -nxyz %d:%d:%d", Imout::MaxX,
		 Imout::MaxY, Imout::MaxZ);
    }
    if (Imout::DeltaValid) {
	sprintf (&buf[strlen(buf)], " -dxyz %g:%g:%g", Imout::DeltaX,
		 Imout::DeltaY, Imout::DeltaZ);
    }
    if (Imout::NvolValid) {
	sprintf (&buf[strlen(buf)], " -nvol %d", Imout::Nvol);
    }
    if (Imout::Orientation) {
	sprintf (&buf[strlen(buf)], " -orient \"%s\"",Imout::Orientation);
    }
    RESULT_BUF (buf);
    return TCL_OK;
}    

int Imout::Close()
{
    if (Imout::Orientation) {
	free (Imout::Orientation);
	Imout::Orientation = 0;
    }
    if (Imout::File) {
	delete Imout::File;
	Imout::File = 0;
    }
    if (Imout::Filename) {
	free (Imout::Filename);
	Imout::Filename = 0;
    }
    Imout::MaxValid = false;
    Imout::DeltaValid = false;
    Imout::NvolValid = false;
    Imout::Rewrite = false;
}

extern "C" int ImoutCmd (ClientData clientData, Tcl_Interp *interp,
			 int argc, char* argv[])
{
    bool session_rewrite = false;
    if (argc == 1)
    {
	if (Imout::Filename)
	{
	    return Imout::Show(interp);
	}
	else
	{
	    RESULT_LIT ("");
	    return TCL_OK;
	}
    }
    int thisindex = 0;

// see if there is a Mask, and none of the arguments are -ignoremask,
// use the Mask to default all dimensions, but later

    bool at_eol = false;
    bool initialized_from_mask = false;
    bool use_mask = false;

    if (Mask::_Valid)
    {
	use_mask = true;
	int i = thisindex;
	while (++i < argc) {
	    char* thisarg = argv[i];
	    if (!strcmp (thisarg, "-ignoremask"))
	    {
		use_mask = false;
		break;
	    }
	}
    }

// Now process arguments one at a time

    while (++thisindex < argc || initialized_from_mask)
    {
	char* thisarg = 0;
	if (thisindex >= argc) {
	    at_eol = true;
	    initialized_from_mask = false;
	} else {
	    thisarg = argv[thisindex];
	}

// -orientation should be last argument, but is tested for
// first because this could be a mask-initialized imout lacking an actual
// -orientation argument.

	if (at_eol || !strcmp (thisarg, "-orient") ||
	         !strcmp (thisarg, "-orientation"))
	{
	    if (!at_eol)
	    {
		if (Imout::Orientation) {
		    free (Imout::Orientation);
		    Imout::Orientation = 0;
		}
		char *spec = argv[thisindex+1];
		if (!spec) {
		    RESULT_LIT ("missing -orient argument");
		}
	    }

// If rewrite object, just change the field
	    
	    if (Imout::Rewrite) {
		if (!Imout::File) {
		    RESULT_LIT ("Must open imout file first");
		    return TCL_ERROR;
		}
	    } else {
		
// if new object, Check that all required information has been set

		bool errorfound = false;
		if (!Imout::Filename) {
		    RESULT_LIT ("Missing imout filename\n");
		    errorfound = true;
		}
		if (!Imout::MaxValid) {
		    RESULT_LIT ("Missing imout -nxyz specification");
		    errorfound = true;
		}
		if (!Imout::DeltaValid) {
		    RESULT_LIT ("Missing imout -dxyz specification");
		    errorfound = true;
		    
		}
		if (!Imout::NvolValid) {
		    RESULT_LIT ("Missing imout -nvol specification");
		    errorfound = true;
		}
		if (errorfound)
		{
		    if (at_eol)
		    {
			Imout::Close();
		    }
		    return TCL_ERROR;
		}

// Now, we actually create binary image object, or open it
// Only if it is a new object

		Imout::File = new RicVolumeSet (Imout::MaxX,Imout::MaxY,
						Imout::MaxZ, Imout::Nvol);
	    }
	    if (at_eol)
	    {
		Imout::File->orientation = Imout::Orientation;
	    } else {
		Imout::File->orientation =Imout::Orientation = Strdup (argv[thisindex+1]);
	    }
	    Imout::Show(interp);
	    thisindex++;
	}
	else if (thisarg[0] != '-')
	{

// First (and only) non-hyphenated argument assumed to be filename
// First clear out last state
//   except for object itself, if filename has stayed the same

	    Imout::MaxValid = false;
	    Imout::DeltaValid = false;
	    Imout::NvolValid = false;
	    if (!session_rewrite) {
		Imout::Rewrite = false;
	    }
	    if (Imout::Orientation) {
		free (Imout::Orientation);
		Imout::Orientation = 0;
	    }
	    if (Imout::Filename)
	    {
		if (!strcmp (Imout::Filename,thisarg))
		{
		    thisindex++;
		    continue;
		}
		free (Imout::Filename);
		Imout::Filename = 0;
	    }
	    if (Imout::File)
	    {
		delete Imout::File;
		Imout::File = 0;
	    }
	    Imout::Filename = Strdup (thisarg);
	    if (session_rewrite) {
		Imout::File =  new RicVolumeSet (Imout::Filename);
		printf ("\nGot RicVolumeSet %s: %ld\n",Imout::Filename,
			(long) Imout::File);
		Imout::MaxX = Imout::File->nx;
		Imout::MaxY = Imout::File->ny;
		Imout::MaxZ = Imout::File->nz;
		Imout::DeltaX = Imout::File->dx;
		Imout::DeltaY = Imout::File->dy;
		Imout::DeltaZ = Imout::File->dz;
		Imout::Nvol = Imout::File->nvol;
		Imout::MaxValid = true;
		Imout::DeltaValid = true;
		Imout::NvolValid = true;
		const char* orient = Imout::File->orientation.c_str();
		Imout::Orientation = Strdup (orient);
		Imout::Show (interp);
	    } else if (use_mask) {
		printf ("setting up values from mask\n");
		Imout::MaxX = Mask::MaxX;
		Imout::MaxY = Mask::MaxY;
		Imout::MaxZ = Mask::MaxZ;
		Imout::DeltaX = Mask::DeltaX;
		Imout::DeltaY = Mask::DeltaY;
		Imout::DeltaZ = Mask::DeltaZ;
		Imout::Orientation = Strdup (Mask::File->orientation.c_str());

		Imout::DeltaValid = true;
		Imout::MaxValid = true;
		Imout::NvolValid = false;
		initialized_from_mask = true;
	    }
	}

// Regular hyphenated arguments	    

	else if (!strcmp (thisarg, "-r"))
	{
// set re-write flag
	    Imout::Rewrite = true;
	    session_rewrite = true;
	}
	else if (!strcmp (thisarg, "-valid"))
	{
	    if (Imout::Filename && Imout::File && Imout::MaxValid
		&& Imout::DeltaValid && Imout::NvolValid &&
		Imout::Orientation)
	    {
		RESULT_LIT ("1");
		return TCL_OK;
	    }
	    else
	    {
		RESULT_LIT ("0");
		return TCL_OK;
	    }
	}
	else if (!strcmp (thisarg, "-nxyz"))
	{
	    Imout::MaxValid = false;
	    char *spec = argv[thisindex+1];
	    if (!spec) {
		if (use_mask) Imout::Close();
		RESULT_LIT ("missing -nxyz argument");
		return TCL_ERROR;
	    }
	    int succeeded = sscanf (spec, "%d:%d:%d", &Imout::MaxX, 
				    &Imout::MaxY, &Imout::MaxZ);
	    if (succeeded != 3)
	    {
		if (use_mask) Imout::Close();
		RESULT_LIT ("Invalid imout -nxyz specification");
		return TCL_ERROR;
	    }
	    Imout::MaxValid = true;
	    thisindex++;
	} 
	else if (!strcmp (thisarg, "-dxyz"))
	{
	    Imout::DeltaValid = false;
	    char *spec = argv[thisindex+1];
	    if (!spec) {
		if (use_mask) Imout::Close();
		RESULT_LIT ("missing -dxyz argument");
		return TCL_ERROR;
	    }
	    int succeeded = sscanf (spec, "%f:%f:%f", &Imout::DeltaX,
				    &Imout::DeltaY,&Imout::DeltaZ);
	    if (succeeded != 3)
	    {
		if (use_mask) Imout::Close();
		RESULT_LIT ("Invalid imout -dxyz specification");
		return TCL_ERROR;
	    }
	    Imout::DeltaValid = true;
	    thisindex++;
	} 
	else if (!strcmp (thisarg, "-nvol"))
	{
	    Imout::NvolValid = false;
	    char *spec = argv[thisindex+1];
	    if (!spec) {
		if (use_mask) Imout::Close();
		RESULT_LIT ("missing -nvol argument");
		return TCL_ERROR;
	    }
	    int succeeded = sscanf (spec,"%d", &Imout::Nvol);
	    if (succeeded != 1)
	    {
		if (use_mask) Imout::Close();
		RESULT_LIT ("Invalid -nvol argument");
		return TCL_ERROR;
	    }
	    Imout::NvolValid = true;
	    thisindex++;
	}
	else if (!strcmp (thisarg, "-ncovar"))
	{
	    Imout::NvolValid = false;
	    char *spec = argv[thisindex+1];
	    if (!spec) {
		if (use_mask) Imout::Close();
		RESULT_LIT ("missing -ncovar argument");
		return TCL_ERROR;
	    }
	    int ncovar;
	    int succeeded = sscanf (spec,"%d", &ncovar);
	    if (succeeded != 1)
	    {
		if (use_mask) Imout::Close();
		RESULT_LIT ("Invalid -ncovar argument");
		return TCL_ERROR;
	    }
	    Imout::Nvol = volumes_required (1,ncovar);
	    Imout::NvolValid = true;
	    thisindex++;
	}
	else if (!strcmp (thisarg, "-ntrait"))
	{
	    Imout::NvolValid = false;
	    char *spec = argv[thisindex+1];
	    if (!spec) {
		if (use_mask) Imout::Close();
		RESULT_LIT ("missing -ntrait argument");
		return TCL_ERROR;
	    }
	    int ntrait;
	    int succeeded = sscanf (spec,"%d", &ntrait);
	    if (succeeded != 1)
	    {
		if (use_mask) Imout::Close();
		RESULT_LIT ("Invalid -ntrait argument");
		return TCL_ERROR;
	    }
	    if (ntrait != 1)
	    {
		if (use_mask) Imout::Close();
		RESULT_LIT ("Only one trait currently supported by imout");
		return TCL_ERROR;
	    }
	    Imout::NvolValid = true;
	    thisindex++;
	}

	else if (!strcmp (thisarg, "-puts") || !strcmp (thisarg, "-put"))
	{
	    if (!Imout::Orientation || !Imout::File) {
		RESULT_LIT ("Imout was not completely specified");
		return TCL_ERROR;
	    }
	    if (argc-thisindex < 3) {
		RESULT_LIT ("Missing arguments for imout -puts");
		return TCL_ERROR;
	    }
	    char *spec = argv[thisindex+1];
	    char *volkey = argv[thisindex+2];
	    char *vol = argv[thisindex+3];
	    if (strcmp(volkey,"-vol")) {
		RESULT_LIT ("Misplaced -vol in imout -puts command");
		return TCL_ERROR;
	    }
	    float value;
	    int succeeded = sscanf (spec,"%g",&value);
	    if (succeeded != 1) {
		RESULT_LIT ("Invalid value in imout -puts");
		return TCL_ERROR;
	    }
	    int volume;
	    succeeded = sscanf (vol, "%d",&volume);
	    if (succeeded != 1 || volume >= Imout::Nvol) {
		char buf[2000];
		sprintf (buf,"Invalid volume value %s in imout -puts",vol);
		RESULT_BUF (buf);
		return TCL_ERROR;
	    }
	    printf ("vol[%d] = %g \n",volume, value);
	    if (!Voxel::Valid()) {
		RESULT_LIT ("Voxel not defined");
		return TCL_ERROR;
	    }
	    int x = Voxel::X;
	    int y = Voxel::Y;
	    int z = Voxel::Z;
	    if (x < 0 || x > Imout::MaxX || y < 0 || y > Imout::MaxY ||
		z < 0 || z > Imout::MaxZ)
	    {
		RESULT_LIT ("Voxel out of range for Imout");
		return TCL_ERROR;
	    }

	    Imout::File->VolSet[volume].vox[x][y][z] = value;
	    thisindex++;
	    thisindex++;
	    thisindex++;
	}
	else if (!strcmp (thisarg, "-write"))
	{
	    if (!Imout::Orientation) {
		RESULT_LIT ("Incomplete imout specification");
		return TCL_ERROR;
	    }
	    int status =
		Imout::File->Write_NIFTI (Imout::Filename);
	    if (status != 1)
	    {
		RESULT_LIT ("Imout write failed");
		return TCL_ERROR;
	    }
	}
	else if (!strcmp (thisarg, "-close"))
	{
	    Imout::Close();
	}
	else if (strcmp (thisarg, "-ignoremask")) // ignore -ignoremask here
	{
	    if (use_mask) Imout::Close();
	    RESULT_LIT ("Invalid imout command");
	    return TCL_ERROR;
	}
    }
    return TCL_OK;
}


extern "C" int MaskCmd (ClientData clientData, Tcl_Interp *interp,
		  int argc, char *argv[])
{

// no arguments returns current state

    if (argc == 1)
    {
	if (Mask::Filename && Mask::_Valid)
	{
	    char buf[10000];
	    sprintf (buf,"mask %s -intensity=%d",Mask::Filename,Mask::Intensity);
	    RESULT_BUF(buf);
	    return TCL_OK;
	}
	else
	{
	    RESULT_LIT ("mask has not been defined\n");
	    return TCL_ERROR;
	}
    }

// one or more arguments sets or deletes mask

    int thisindex = 1;
    Mask::_Valid = false;
    bool invalid_found = false;
    int voxel_offset = 0;

    while (thisindex < argc)
    {
	char* thisarg = argv[thisindex];
	if (thisarg[0] != '-')
	{

// Non-hyphenated argument assumed to be filename

	    if (Mask::Filename && Mask::FileValid)
	    {
		if (!strcmp (Mask::Filename,thisarg))
		{
//		    printf ("Using already loaded mask\n");
		    continue;
		}
	    }
	    Mask::FileValid = false;
	    if (Mask::Filename)
	    {
		free (Mask::Filename);
		Mask::Filename = 0;
	    }
	    if (Mask::File)
	    {
		delete Mask::File;
		Mask::File = 0;
	    }

	    printf ("Opening mask file %s...\n",thisarg);
	    Mask::File = new RicVolumeSet (thisarg);
	    printf ("\nMask file opened %ld\n", (long) Mask::File);

	    if (!Mask::File)
	    {
		RESULT_LIT ("Unable to open RicVolumeSet");
		return TCL_ERROR;
	    }
	    Mask::MaxX = Mask::File->nx;
	    Mask::MaxY = Mask::File->ny;
	    Mask::MaxZ = Mask::File->nz;
	    Mask::DeltaX = Mask::File->dx;
	    Mask::DeltaY = Mask::File->dy;
	    Mask::DeltaZ = Mask::File->dz;
	    Mask::Filename = Strdup (thisarg);
	    Mask::FileValid = true;
	    Mask::Intensity = 0;
	    Mask::Index = 0;
	}
	else if (!Mask::FileValid)
	{
	    RESULT_LIT ("Mask filename must be specified first");
	    return TCL_ERROR;
	}
	else
	{
//	    printf ("checking for arguments\n");
	    if (!strcmp (thisarg, "-intensity"))
	    {
		if (++thisindex > argc)
		{
		    RESULT_LIT ("Missing mask intensity value");
		    return TCL_ERROR;
		}
		char* endptr;
		errno = 0;
		long larg = strtol (argv[thisindex], &endptr, 0);
		if (errno || larg > INT_MAX || larg < INT_MIN)
		{
		    RESULT_LIT ("Invalid mask intensity value");
		    return TCL_ERROR;
		}
		Mask::Intensity = (int) larg;
	    }
	    else if (!strcmp (thisarg, "-gt"))
	    {
		if (++thisindex > argc)
		{
		    RESULT_LIT ("Missing mask intensity value");
		    return TCL_ERROR;
		}
		char* endptr;
		errno = 0;
		long larg = strtol (argv[thisindex], &endptr, 0);
		if (errno || larg > INT_MAX || larg < INT_MIN)
		{
		    RESULT_LIT ("Invalid mask intensity value");
		    return TCL_ERROR;
		}
		Mask::Intensity = 0;
		Mask::Gt = 1;
		Mask::Ge = 0;
	    }
	    else if (!strcmp (thisarg, "-index"))
	    {
		if (++thisindex > argc)
		{
		    RESULT_LIT ("Missing mask index value");
		    return TCL_ERROR;
		}
		char* endptr;
		errno = 0;
		long larg = strtol (argv[thisindex], &endptr, 0);
		if (errno || larg > INT_MAX || larg < INT_MIN)
		{
		    RESULT_LIT ("Invalid mask index value");
		    return TCL_ERROR;
		}
		voxel_offset = (int) larg;
	    }
	    else if (!strcmp (thisarg, "-next"))
	    {
		voxel_offset=-1;
	    }
	    else if (!strcmp (thisarg, "-delete"))
	    {
		Mask::FileValid = false;
		Mask::_Valid = false;
		delete Mask::File;
		Mask::File = 0;
		Mask::Index = -1;
		Mask::Intensity = -1;
		return TCL_OK;
	    }
	    else
	    {
		RESULT_LIT ("Invalid mask command\n");
		Mask::_Valid = false;
		return TCL_ERROR;
	    }
	}
	Mask::_Valid = true;
	thisindex++;
    }
    if (Mask::FileValid)
    {
	Mask::_Valid = true;
	return Mask::Move (interp, voxel_offset);
    }
    return TCL_OK;
}


// If offset is negative, move to next offset from current position
// If offset is zero, move to first position
// If offset is positive, move to first position (zero) plus offset

int Mask::Move (Tcl_Interp* interp, int user_offset)
{
    printf ("Entering Mask::Move");
    int offset = 0;
    int passed = 0;
    int x=0;
    int y=0;
    int z=0;

    if (user_offset<0) // currently all negatives handled as -next
    {
	x=Voxel::X+1;
	y=Voxel::Y;
	z=Voxel::Z;
	offset = 0;
    }
    else
    {
	offset = user_offset;
    }
    
    printf ("offset is %d\n",offset);
    printf ("intensity is %d\n",Mask::Intensity);

    printf ("MaxX: %d  MaxY: %d  MaxZ: %d\n", MaxX, MaxY, MaxZ);

    for (; z < MaxZ; z++)
    {
	for (; y < MaxY; y++)
	{
	    for (; x < MaxX; x++)
	    {
		float val = Mask::File->VolSet[0].vox[x][y][z];
//		printf ("val[%d,%d,%d] = %d\n", x,y,z,(int) val);
		if ((!Mask::Intensity && val != 0) ||
		    (Mask::Intensity && Mask::Intensity == val))
		{
//		    printf ("checking passed %d offset %d\n",passed+1, offset);
		    if (passed++ >= offset)
		    {
//			printf ("setting voxel\n");
			
			char vbuf[128];
			return Voxel::Set (interp,x,y,z);
		    }
		}
	    }
	    x = 0;
	}
	y = 0;
    }
    if (user_offset < 0)
    {
	RESULT_LIT ("No more voxels in mask");
	return TCL_ERROR;
    }
    RESULT_LIT ("No matching voxels in mask");
    return TCL_ERROR;
}

	
extern "C" int VoxelCmd (ClientData clientData, Tcl_Interp *interp,
		  int argc, char *argv[])
{
    if (argc == 1)
    {
	return Voxel::Show (interp);
    }
    if (argc == 2)
    {
	return Voxel::Set (interp, argv[1]);
    }
    return 0;
}

int Voxel::Show (Tcl_Interp *interp)
{
    if (!Voxel::Valid())
    {
	RESULT_LIT ("Voxel not defined; mask command is usually used");
	return TCL_ERROR;
    }
    char buf[256];
    sprintf (buf,"%d:%d:%d",X,Y,Z);
    RESULT_BUF(buf);
    return TCL_OK;
}

void Voxel::write_commands (FILE *file)
{
    if (Voxel::Valid())
    {
	fprintf (file,"voxel %d:%d:%d",X,Y,Z);
    }
}


int Voxel::Set (Tcl_Interp* interp, char* voxel_spec)
{
    _Valid = false;
    int succeeded = sscanf (voxel_spec, "%d:%d:%d", &X, &Y, &Z);
    if (succeeded != 3)
    {
	RESULT_LIT ("Invalid Voxel specification");
	return TCL_ERROR;
    }
    _Valid = true;
    return Show (interp);
}

int Voxel::Set (Tcl_Interp* interp, int x, int y, int z)
{
    X=x;
    Y=y;
    Z=z;
    _Valid = true;
    return Show (interp);
}
