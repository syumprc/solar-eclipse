// ----------------------- RicVolumeSet.cpp -------------------------
/*! \file
Source file for the RicVolumeSet class. The class can open a variety
of different volume file formats.

Copyright (C) 2007 by Bill Rogers - Research Imaging Center - UTHSCSA
rogers@uthscsa.edu
 */
using namespace std;

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include "RicVolumeSet.h"
#include <nemardr.h>
#include <nifti1_io.h>

// functions from GameDev.net for swapping bytes
float FloatSwap( float f )
{
  union
  {
    float f;
    unsigned char b[4];
  } dat1, dat2;

  dat1.f = f;
  dat2.b[0] = dat1.b[3];
  dat2.b[1] = dat1.b[2];
  dat2.b[2] = dat1.b[1];
  dat2.b[3] = dat1.b[0];
  return dat2.f;
};

int LongSwap (int i)
{
  unsigned char b1, b2, b3, b4;

  b1 = i & 255;
  b2 = ( i >> 8 ) & 255;
  b3 = ( i>>16 ) & 255;
  b4 = ( i>>24 ) & 255;

  return ((int)b1 << 24) + ((int)b2 << 16) + ((int)b3 << 8) + b4;
};

short ShortSwap( short s )
{
  unsigned char b1, b2;

  b1 = s & 255;
  b2 = (s >> 8) & 255;

  return (b1 << 8) + b2;
};



/*!
Constructor for empty volume
 */
RicVolumeSet::RicVolumeSet()
{
	VolSet = NULL;
	nvol = 0;
	nz = 0;
	ny = 0;
	nx = 0;
	totvox = 0;
	scale = 1;

	orientation = "XYZ---";
	s_units = RIC_UNITS_MM;
	t_units = RIC_UNITS_SEC;
	dtype = RIC_INTEGER;
	dx = 1;
	dy = 1;
	dz = 1;
	dt = 1;
	xoffset = 0;
	yoffset = 0;
	zoffset = 0;
	TE1 = 0;
	TE2 = 0;
	TR = 0;

	filetype = RIC_NO_FILE;
}

/*!
Constructor to make a volume set with all voxel values set to zero

\param numx
number of x values (columns)

\param numy
number of y values (rows)

\param numz
number of z values (slices)

\param nv
number of volumes
 */
RicVolumeSet::RicVolumeSet(int numx, int numy, int numz, int nv)
{
	nvol = nv;
	nx = numx;
	ny = numy;
	nz = numz;
	totvox = nvol*nz*ny*nx;
	scale = 1;

	orientation = "XYZ---";
	s_units = RIC_UNITS_MM;
	t_units = RIC_UNITS_SEC;
	dtype = RIC_INTEGER;
	dx = 1;
	dy = 1;
	dz = 1;
	dt = 1;
	xoffset = 0;
	yoffset = 0;
	zoffset = 0;

	// file info
	TE1 = 0;
	TE2 = 0;
	TR = 0;
	scan_date = "";
	pname = "";
	pnum = "";

	filetype = RIC_NO_FILE;

	// allocate new volumes
	VolSet = new RicVolume[nv];

	for ( int i=0 ; i<nv ; ++i )
		VolSet[i].Init(nx,ny,nz);
}

/*!
Constructor to read in a volume from a file.
The filename is parsed to determine what kind of file to read in.

\param fname
Name of file containing volume set
 */
RicVolumeSet::RicVolumeSet(string fname)
{
	VolSet = NULL;
	// parse fname to figure out what kind of file we have

	// check to see if NEMA file
	if ( fname.find(".des") != string::npos
			|| fname.find(".DES") != string::npos )
	{
		this->Read_NEMA(fname);
	}
	// check for nifti or analyze
	else if ( fname.find(".hdr") != string::npos || fname.find(".HDR") != string::npos ||
				 fname.find(".nii") != string::npos || fname.find(".NII") != string::npos )
	{
		this->Read_NIFTI(fname);
	}
	else if ( fname.find(".dim") != string::npos
			|| fname.find(".DIM") != string::npos )
	{
		this->Read_GIS(fname);
	}
	else if ( fname.find(".ima") != string::npos
			|| fname.find(".IMA") != string::npos )
	{
		this->Read_GIS(fname);
	}
	else // cannot read file so zero everything out
	{
		cerr << "RicVolumeSet - Unknown file type " << fname << endl;
		VolSet = NULL;
		nvol = 0;
		nz = 0;
		ny = 0;
		nx = 0;
		totvox = 0;
		scale = 1;
		filetype = RIC_NO_FILE;
	}
}

/*!
Destructor frees memory from volume allocations
 */
RicVolumeSet::~RicVolumeSet()
{
	delete [] VolSet;
}

/*!
Initializer allocates memory for a set of volumes.

\param numx
number of x values (columns)

\param numy
number of y values (rows)

\param numz
number of z values (slices)

\param nv
number of volumes

\returns
Returns 1 on sucess, 0 on failure
*/
int RicVolumeSet::Init(int numx, int numy, int numz, int nv)
{
	// free memory if already allocated
	if ( VolSet ) delete [] VolSet;

	nvol = nv;
	nx = numx;
	ny = numy;
	nz = numz;
	totvox = nvol*nz*ny*nx;
	scale = 1;

	orientation = "XYZ---";
	s_units = RIC_UNITS_MM;
	t_units = RIC_UNITS_SEC;
	dtype = RIC_INTEGER;
	dx = 1;
	dy = 1;
	dz = 1;
	dt = 1;
	xoffset = 0;
	yoffset = 0;
	zoffset = 0;

	// file info
	TE1 = 0;
	TE2 = 0;
	TR = 0;
	scan_date = "";
	pname = "";
	pnum = "";

	filetype = RIC_NO_FILE;

	// allocate new volumes
	VolSet = new RicVolume[nv];

	for ( int i=0 ; i<nv ; ++i )
		VolSet[i].Init(nx,ny,nz);

	return 1;
}

/*!
Function to write a volume set to a file. The file extension
is used to determine what format to write. If there is no
extension then NEMA will be used.

\param fname
name of file to write volume set to

\returns
1 on success, 0 on failure
 */
int RicVolumeSet::Write(string fname)
{
	// parse fname to figure out what kind of file we have

	// check to see if NEMA file
	if ( fname.find(".des") != string::npos
			|| fname.find(".DES") != string::npos )
	{
		return this->Write_NEMA(fname);
	}
	// check for nifti or analyze
	else if ( fname.find(".hdr") != string::npos || fname.find(".HDR") != string::npos )
	{
		this->filetype = RIC_ANALYZE_FILE;
		return this->Write_NIFTI(fname);
	}
	else if ( fname.find(".nii") != string::npos || fname.find(".NII") != string::npos )
	{
		this->filetype = RIC_NIFTI_FILE;
		return this->Write_NIFTI(fname);
	}
	else if ( fname.find(".dim") != string::npos || fname.find(".DIM") != string::npos )
	{
		this->filetype = RIC_GIS_FILE;
		return this->Write_GIS(fname);
	}
	else // just do NEMA
	{
		return this->Write_NEMA(fname);
	}
}

/*!
Function to read in a volume from a NEMA (.des/.dat) file pair.
If there is more than one volume/time series in the file, only
the first will be read.

\param name
Name of file containing volume set

\returns
1 on success, 0 on failure
 */
int RicVolumeSet::Read_NEMA(string name)
{
	char fname[FILENAME_SIZE];
	strcpy(fname,name.c_str());
	NEMAFILE *nfile=NULL;

	// read in the NEMA header file
	this->nvol = NEMA_READ(&nfile, fname);

	if ( this->nvol == 0 )
		return 0;

	// note the file type
	this->filetype = RIC_NEMA_FILE;

	// copy meaningful stuff to member variables
	// assume all volumes have the same stuff
	this->nz = nfile[0].tot_scans();
	this->ny = nfile[0].num_rows();
	this->nx = nfile[0].num_cols();
	this->totvox = nvol*nz*ny*nx;
	this->orientation = nfile[0].get_orientation();
	string units;
	units = nfile[0].spatial_type(); // don't really need this as spec says units always mm
	this->s_units = RIC_UNITS_MM;
	this->t_units = RIC_UNITS_SEC;
	this->dx = nfile[0].get_col_mm();
	this->dy = nfile[0].get_row_mm();
	this->dz = nfile[0].get_slicevec();
	XYZPOINT offset;
	nfile[0].get_XYZoffset(offset);
	this->xoffset = offset[0];
	this->yoffset = offset[1];
	this->zoffset = offset[2];
	this->scale = 1.0f;	// data will be scaled on input

	// file info
	this->TE1 = nfile[0].get_echo1_time();
	this->TE2 = nfile[0].get_echo2_time();
	this->TR = nfile[0].get_rep_time1();
	this->scan_date = nfile[0].scan_date;
	this->pnum = nfile[0].patient_number;
	this->pname = nfile[0].name;

	//char date[32];
	//nfile[0].get_scan_date(date);
	//this->scan_date = date;
	//this->pname = nfile[0].get_patient_name();
	//this->pnum = nfile[0].get_patient_number();

	// allocate memory for a set of volumes
	this->VolSet = new RicVolume[nvol];

	// Low hanging fruit warning!!!!
	// at this time we are only reading data with XYZ orientation
	if ( this->orientation.compare(0,3,"XYZ") != 0 )
	{
		cerr << "We only read NEMA files with XYZ orientation" << endl;
		return 0;
	}

	// parse the orientation string
	int xdir=1,ydir=1,zdir=1;
	if ( this->orientation.compare(3,1,"-") != 0 ) xdir = -1;
	if ( this->orientation.compare(4,1,"-") != 0 ) ydir = -1;
	if ( this->orientation.compare(5,1,"-") != 0 ) zdir = -1;

	// force the orientation to be XYZ---
	this->orientation = "XYZ---";

	// get file type
	int ntype = nfile[0].ptype();
	int bswap=0;
	int highbit = nfile[0].hbp();
	int nbits = nfile[0].bpp();
	if (ntype == NEMA_SIGNED)
	{
		dtype = RIC_INTEGER;
		if ( highbit == 15 ) bswap = 1;
	}
	else if (ntype == NEMA_UNSIGNED)
	{
		dtype = RIC_INTEGER;
		if ( highbit == 15 ) bswap = 1;
	}
	else if (ntype == NEMA_IEEE_FLOAT)
	{
		dtype = RIC_FLOAT;
		if ( highbit == 31 ) bswap = 1;
	}
	else if (ntype == NEMA_ASCII)
	{
		dtype = RIC_INTEGER; // this is as small as we go
	}
	else
	{
		cerr << "Error unknown NEMA file type" << endl;
		return 0;
	}

	// copy data from NEMA data file to local matrix
	// NOTE - assuming a single data file
	FILE *fptr;
	fptr = fopen(nfile[0].fname(0),"rb");

	// repeat for every volume in file
	for (int n = 0; n < nvol; ++n)
	{
		// initialize the volume
		this->VolSet[n].Init(nx,ny,nz);
		this->VolSet[n].dx = dx;
		this->VolSet[n].dy = dy;
		this->VolSet[n].dz = dz;
		this->VolSet[n].xoffset = xoffset;
		this->VolSet[n].yoffset = yoffset;
		this->VolSet[n].zoffset = zoffset;
		this->filetype = RIC_NEMA_FILE; // init kills this

		// allocate temporary data for reading and organizing file
		float*** data;
		matrix3D(&data, nz, ny, nx);
		int i, j, k;

		// keep track of original pointers as we may swap them around
		// if they are in the wrong location then problems arise in freeing the matrix
		float *pntr0 = data[0][0];
		float **pntr1 = data[0];

		// read in the data in large chunks
		int nread;
		int ntoread = nx * ny;

		if (ntype == NEMA_IEEE_FLOAT) // float - assuming no scaling
		{
			for (i = 0; i < nz; ++i)
			{
				nread = fread(&data[i][0][0], sizeof(float), ntoread,fptr);
				if (ferror(fptr))
				{
					cerr << "Error reading file" << nfile[0].fname(0) << endl;
					return 0;
				}

				if (bswap)
				{
					float *tptr;
					tptr = &data[i][0][0];
					for (int l = 0; l < ntoread; ++l)
						tptr[l] = FloatSwap(tptr[l]);
				}
			}
		}
		else if ((ntype == NEMA_UNSIGNED && nbits == 8) || ntype == NEMA_ASCII)
		{
			// 8 bit character
			unsigned char **ctmp;
			matrix(&ctmp, ny, nx);
			for (i = 0; i < nz; ++i)
			{
				// read a slice of character data
				nread = fread(&ctmp[0][0], sizeof(unsigned char), ntoread,
								fptr);
				if (ferror(fptr))
				{
					cerr << "Error reading file" << nfile[0].fname(0) << endl;
					return 0;
				}

				// copy to float data array
				for (j = 0; j < ny; ++j)
					for (k = 0; k < nx; ++k)
						data[i][j][k] = nfile[n].data_scale[i]*ctmp[j][k];
			}
			free_matrix(ctmp);
		}
		else if (ntype == NEMA_UNSIGNED) // unsigned short
		{
			unsigned short **ctmp;
			matrix(&ctmp, ny, nx);
			for (i = 0; i < nz; ++i)
			{
				nread = fread(&ctmp[0][0], sizeof(unsigned short), ntoread,
						fptr);
				if (ferror(fptr))
				{
					cerr << "Error reading file" << nfile[0].fname(0) << endl;
					return 0;
				}

				if (bswap)
				{
					unsigned short *tptr;
					tptr = &ctmp[0][0];
					for (int l = 0; l < ntoread; ++l)
						tptr[l] = ShortSwap(tptr[l]);
				}

				// copy to float data array
				for (j = 0; j < ny; ++j)
					for (k = 0; k < nx; ++k)
						data[i][j][k] = nfile[n].data_scale[i]*ctmp[j][k];
			}
			free_matrix(ctmp);
		}
		else if (ntype == NEMA_SIGNED) // signed short
		{
			short **ctmp;
			matrix(&ctmp, ny, nx);
			for (i = 0; i < nz; ++i)
			{
				nread = fread(&ctmp[0][0], sizeof(short), ntoread,
						fptr);
				if (ferror(fptr))
				{
					cerr << "Error reading file" << nfile[0].fname(0) << endl;
					return 0;
				}

				if (bswap)
				{
					short *tptr;
					tptr = &ctmp[0][0];
					for (int l = 0; l < ntoread; ++l)
						tptr[l] = ShortSwap(tptr[l]);
				}

				// copy to float data array
				for (j = 0; j < ny; ++j)
					for (k = 0; k < nx; ++k)
						data[i][j][k] = nfile[n].data_scale[i]*ctmp[j][k];
			}
			free_matrix(ctmp);
		}
		else // unknown number of bits
		{
			cerr << "Error unknown number of bits in file" << endl;
			return 0;
		}

		// swap indices in arrays to compensate for orientation

		// see if we have to swap z indices
		if (zdir != 1)
		{
			float ***zptr;
			zptr = new float**[nz];
			for (i = 0; i < nz; ++i)
				zptr[i] = data[nz - i - 1];

			for (i = 0; i < nz; ++i)
				data[i] = zptr[i];
			delete[] zptr;
		}

		// see if we have to swap y indices
		if (ydir != 1)
		{
			float **yptr;
			yptr = new float*[ny];
			for (i = 0; i < nz; ++i)
			{
				for (j = 0; j < ny; ++j)
					yptr[j] = data[i][ny - j - 1];

				for (j = 0; j < ny; ++j)
					data[i][j] = yptr[j];
			}
			delete[] yptr;
		}

		// see if we have to swap x indices
		if (xdir != 1)
		{
			float *xptr;
			xptr = new float[nx];
			for (i = 0; i < nz; ++i)
			{
				for (j = 0; j < ny; ++j)
				{
					for (k = 0; k < nx; ++k)
						xptr[k] = data[i][j][nx - k - 1];

					for (k = 0; k < nx; ++k)
						data[i][j][k] = xptr[k];
				}
			}
			delete[] xptr;
		}

		// copy to volume array
		for (i = 0; i < nz; ++i)
		{
			for (j = 0; j < ny; ++j)
				for (k = 0; k < nx; ++k)
					VolSet[n].vox[k][j][i] = data[i][j][k];
		}


		// pointer to data array
		VolSet[n].varray = &VolSet[n].vox[0][0][0];

		VolSet[n].CalcMinMaxAvg();

		// remove temporary array - avoiding pointer swapping errors
		data[0] = pntr1;
		data[0][0] = pntr0;

		free_matrix3D(data);
	}
	fclose(fptr);

	// delete NEMA header structure
	delete [] nfile;

	return 1;
}

/*!
Function to write volumes to a NEMA (.des/.dat) file pair. Data can be
output as signed integer or float. The output file will always be little
endian (on Intel).

\param filename
name of file to write volume set to

\returns
1 on success, 0 on failure
 */
int RicVolumeSet::Write_NEMA(string fname)
{
	// just use the current orientation
	return this->Write_NEMA(fname,this->orientation);
}

/*!
Function to write volumes to a NEMA (.des/.dat) file pair. This varient will
write the data out in the orientation specified by the passed NEMA orientation
string. (Note must be in XYZ order) Data can be output as signed integer
or float. The output file will always be little endian (on Intel).

\param filename
name of file to write volume set to

\param orient_str
orientation string

\returns
1 on success, 0 on failure
 */
int RicVolumeSet::Write_NEMA(string fname, string orient_str)
{
	// remove extension from filename if necessary
	size_t pos,l;

	if ( (pos=fname.find(".des")) != string::npos )
	{
		l = fname.length();
		fname.erase(pos,l-pos);
	}
	if ( (pos=fname.find(".DES")) != string::npos )
	{
		l = fname.length();
		fname.erase(pos,l-pos);
	}

	// make sure there is something to write
	if ( VolSet[0].nz == 0 ) return 0;

	// make sure we have an orientation string we can work with
	if ( orient_str.length() != 6 ) // must be 6 characters
	{
		cerr << "Write_NEMA - invalid orientation string length"
			<< orient_str << endl;
		return 0;
	}
	string xyz = orient_str.substr(0,3); // must be of type XYZ
	if ( xyz.find("XYZ") == string::npos )
	{
		cerr << "Write_NEMA - invalid orientation string prefix"
			<< orient_str << endl;
		return 0;
	}

	// add des to the filename
	string desname = fname+".des";

	NEMA_Start_Header((char*) desname.c_str(), nvol);

	// output file name minus the path
	string dname;
	unsigned idx = fname.rfind("/");
	if (idx != string::npos)
		dname = fname.substr(idx + 1, string::npos);
	else
		dname = fname;

	string dataname = dname + ".dat";

	// open binary output file
	FILE *fptr;
	fptr = fopen(dataname.c_str(), "wb");

	// repeat for every volume in set
	for (int n = 0; n < nvol; ++n)
	{
		NEMAFILE *nfile = new NEMAFILE();

		// Allocate space for the filenames, offsets, and positions, slice_min, slice_max and data_scale
		nfile->init(nz);

		// copy local stuff to NEMA header
		nfile->set_desfname((char*) desname.c_str());
		nfile->set_image_min(VolSet[n].min);
		nfile->set_image_max(VolSet[n].max);
		nfile->set_dims(nx, ny);
		nfile->set_sizes(dx, dy, dz);
		if (dtype == RIC_FLOAT) // assume float
			nfile->set_type(NEMA_IEEE_FLOAT, 32, 32, 7); // force high bit to 7
		else
			// assume 16 bits
			nfile->set_type(NEMA_SIGNED, 16, 16, 7); // force high bit to 7
		nfile->set_orientation((char*) orient_str.c_str()); // used passed orientation
		nfile->set_spatial_type((char*)"MM"); // units always mm by spec
		XYZPOINT offset;
		offset[0] = xoffset;
		offset[1] = yoffset;
		offset[2] = zoffset;
		nfile->set_XYZoffset(offset);
		nfile->set_xyorigin(0, 0); // Set the origin (world coordinates) of the image

		// file info
		nfile->set_echo1_time(TE1);
		nfile->set_echo2_time(TE2);
		nfile->set_rep_time1(TR);
		strcpy(nfile->scan_date,scan_date.c_str());
		strcpy(nfile->patient_number,pnum.c_str());
		strcpy(nfile->name,pname.c_str());

		// get size of a slice
		long imgsize;
		if (dtype == RIC_FLOAT)
			imgsize = ny * nx * 4;
		else
			imgsize = ny * nx * 2;

		int i, j, k;
		for (i = 0; i < nz; i++)
		{
			nfile->set_pos(i * dz, i);
			nfile->set_fname((char*) dataname.c_str(), i);
			nfile->set_offset(i * imgsize, i);
			nfile->set_slice_min(0.0L, i);
			nfile->set_slice_max(0.0L, i);
			nfile->set_slice_data_scale(scale, i);
		}

		// append to header file
		nfile->write_vol(n+1);

		// copy data to temporary array
		// allocate data
		float*** data;
		matrix3D(&data, nz, ny, nx);

		// keep track of original pointers as we may swap them around
		// if they are in the wrong location then problems arise in freeing the matrix
		float *pntr0 = data[0][0];
		float **pntr1 = data[0];


		// copy from volume to temporary array
		for (i = 0; i < nz; ++i)
		{
			for (j = 0; j < ny; ++j)
			{
				for (k = 0; k < nx; ++k)
				{
					data[i][j][k] = VolSet[n].vox[k][j][i];
				}
			}
		}

		// now we swap the data indices around to match the orientation string

		// parse the orientation string
		int xdir = 1, ydir = 1, zdir = 1;
		if (orient_str.compare(3, 1, "-") == 0)
			xdir = -1;
		if (orient_str.compare(4, 1, "-") == 0)
			ydir = -1;
		if (orient_str.compare(5, 1, "-") == 0)
			zdir = -1;

		// see if we have to swap z indices
		if (zdir != -1)
		{
			float ***zptr;
			zptr = new float**[nz];
			for (i = 0; i < nz; ++i)
				zptr[i] = data[nz - i - 1];

			for (i = 0; i < nz; ++i)
				data[i] = zptr[i];
			delete[] zptr;
		}

		// see if we have to swap y indices
		if (ydir != -1)
		{
			float **yptr;
			yptr = new float*[ny];
			for (i = 0; i < nz; ++i)
			{
				for (j = 0; j < ny; ++j)
					yptr[j] = data[i][ny - j - 1];

				for (j = 0; j < ny; ++j)
					data[i][j] = yptr[j];
			}
			delete[] yptr;
		}

		// see if we have to swap x indices
		if (xdir != -1)
		{
			float *xptr;
			xptr = new float[nx];
			for (i = 0; i < nz; ++i)
			{
				for (j = 0; j < ny; ++j)
				{
					for (k = 0; k < nx; ++k)
						xptr[k] = data[i][j][nx - k - 1];

					for (k = 0; k < nx; ++k)
						data[i][j][k] = xptr[k];
				}
			}
			delete[] xptr;
		}

		// write the data a row at a time.
		if (dtype == RIC_INTEGER)
		{
			short *rowarray;
			rowarray = new short[nx]; // write a row at a time
			for (i = 0; i < nz; ++i)
			{
				for (j = 0; j < ny; ++j)
				{
					for (k = 0; k < nx; ++k)
					{
						rowarray[k] = (short) (data[i][j][k] + 0.5);
					}
					fwrite(rowarray, sizeof(short), nx, fptr);
				}
			}
			delete rowarray;
		}
		else // must be 32 bit float
		{
			float *rowarray;
			rowarray = new float[nx]; // write a row at a time
			for (i = 0; i < nz; ++i)
			{
				for (j = 0; j < ny; ++j)
				{
					for (k = 0; k < nx; ++k)
					{
						rowarray[k] = data[i][j][k];
					}
					fwrite(rowarray, sizeof(float), nx, fptr);
				}
			}
			delete rowarray;
		}


		// remove temporary array - avoiding pointer swapping errors
		data[0] = pntr1;
		data[0][0] = pntr0;

		free_matrix3D(data);

		// delete NEMA header
		delete nfile;
	}

	// close binary data file
	fclose(fptr);

	return 1;
}

/*!
Function to read in a volume from a NIFTI (.nii) file or ANALYSE (.hdr/.img) file pair.
This function will read all volumes in a time series.

\param filename
Name of file containing volume set

\returns
1 on success, 0 on failure
 */
int RicVolumeSet::Read_NIFTI(string fname)
{
	// see if the nifti file is real
	int nfiletype = is_nifti_file((char*)fname.c_str());
	if ( nfiletype < 0 ) return 0;

	// note the file type
	if ( nfiletype == NIFTI_FTYPE_ANALYZE )
		this->filetype = RIC_ANALYZE_FILE;
	else if ( nfiletype == NIFTI_FTYPE_NIFTI1_1 || nfiletype == NIFTI_FTYPE_NIFTI1_2 )
		this->filetype = RIC_NIFTI_FILE;
	else // we cannot handle it
		return 0;

	// read in header and image info
	nifti_image      * nim;
	nim = nifti_image_read((char*)fname.c_str(),1);

	// set up some nifti variables to read
	nifti_brick_list   brick_list;

	// load all bricks (volumes) at one time
	int nbricks = nifti_image_load_bricks( nim,nim->nt,NULL, &brick_list);

	// copy over a bunch of header info
	this->totvox = nbricks*nz*ny*nx;
	this->filename = fname;
	this->nvol = nbricks;
	this->dx = nim->dx;
	this->dy = nim->dy;
	this->dz = nim->dz;
	this->dt = nim->dt;

//	if ( nim->byteorder == 2 )
//		this->highbit = 15;
//	else
//		this->highbit = 0;
//	this->nbits =  nim->nbyper*8;

	if ( nim->datatype == DT_SIGNED_SHORT )
		this->dtype = RIC_INTEGER;
	else if (nim->datatype == DT_UINT16 )
		this->dtype = RIC_INTEGER;
	else if ( nim->datatype == DT_UNSIGNED_CHAR )
		this->dtype = RIC_INTEGER;
	else if (nim->datatype == DT_FLOAT)
		this->dtype = RIC_FLOAT;
	else
	{
		delete nim;
		return 0; // can't do any other format
	}

	// allocate new volumes
	this->Init(nim->nx, nim->ny, nim->nz, nbricks);

	// set scaling factor to one - data will be properly scaled on read
	this->scale = 1;

	// the the XYZ orientation
	this->orientation = convertNiftiSFormToNEMA(nim->qto_xyz);

	/// read in each volume of data
	int i,j,k,n;
	for (n = 0; n < this->nvol; ++n)
	{
		// allocate data
		float*** data;
		matrix3D(&data,nx,ny,nz);

		// keep track of original pointers as we may swap them around
		// if they are in the wrong location then problems arise in freeing the matrix
		float *pntr0 = data[0][0];
		float **pntr1 = data[0];

		// assign dx dy dz to each volume
		this->VolSet[n].dx = dx;
		this->VolSet[n].dy = dy;
		this->VolSet[n].dz = dz;

		/// Read in data to generic matrix

		if (nim->datatype == DT_SIGNED_SHORT)
		{
			short *dptr;
			// assign dx dy dz to each volume
			this->VolSet[n].dx = dx;
			this->VolSet[n].dy = dy;
			this->VolSet[n].dz = dz;

			dptr = (short*) brick_list.bricks[n];

			if (nfiletype == NIFTI_FTYPE_ANALYZE) // swap slice and row order
			{
				for (i = nim->nz - 1; i >= 0; --i)
				{
					for (j = nim->ny - 1; j >= 0; --j)
					{
						for (k = 0; k < nim->nx; ++k)
						{
							data[k][j][i] = *dptr;
							dptr++;
						}
					}
				}
			}
			else // just read the data - NIFTI
			{
				for (i = 0; i < nim->nz; ++i)
				{
					for (j = 0; j < nim->ny; ++j)
					{
						for (k = 0; k < nim->nx; ++k)
						{
							data[k][j][i] = *dptr
									* nim->scl_slope + nim->scl_inter;
							dptr++;
						}
					}
				}

			}

		}
		else if (nim->datatype == DT_UINT16)
		{
			unsigned short *uptr;

			uptr = (unsigned short*) brick_list.bricks[n];

			if (nfiletype == NIFTI_FTYPE_ANALYZE) // swap slice and row order
			{
				for (i = nim->nz - 1; i >= 0; --i)
				{
					for (j = nim->ny - 1; j >= 0; --j)
					{
						for (k = 0; k < nim->nx; ++k)
						{
							data[k][j][i] = *uptr;
							uptr++;
						}
					}
				}
			}
			else // just read the data - NIFTI
			{
				for (i = 0; i < nim->nz; ++i)
				{
					for (j = 0; j < nim->ny; ++j)
					{
						for (k = 0; k < nim->nx; ++k)
						{
							data[k][j][i] = *uptr
									* nim->scl_slope + nim->scl_inter;
							uptr++;
						}
					}
				}

			}

		}
		else if (nim->datatype == DT_UNSIGNED_CHAR)
		{
			unsigned char *uptr;

			uptr = (unsigned char*) brick_list.bricks[n];

			if (nfiletype == NIFTI_FTYPE_ANALYZE) // swap slice and row order
			{
				for (i = nim->nz - 1; i >= 0; --i)
				{
					for (j = nim->ny - 1; j >= 0; --j)
					{
						for (k = 0; k < nim->nx; ++k)
						{
							data[k][j][i] = *uptr;
							uptr++;
						}
					}
				}
			}
			else // just read the data - NIFTI
			{
				for (i = 0; i < nim->nz; ++i)
				{
					for (j = 0; j < nim->ny; ++j)
					{
						for (k = 0; k < nim->nx; ++k)
						{
							data[k][j][i] = *uptr
									* nim->scl_slope + nim->scl_inter;
							uptr++;
						}
					}
				}

			}

		}
		else if (nim->datatype == DT_FLOAT)
		{
			float *uptr;

			uptr = (float*) brick_list.bricks[n];

			if (nfiletype == NIFTI_FTYPE_ANALYZE) // swap slice and row order
			{
				for (i = nim->nz - 1; i >= 0; --i)
				{
					for (j = nim->ny - 1; j >= 0; --j)
					{
						for (k = 0; k < nim->nx; ++k)
						{
							data[k][j][i] = *uptr;
							uptr++;
						}
					}
				}
			}
			else // just read the data - NIFTI
			{
				for (i = 0; i < nim->nz; ++i)
				{
					for (j = 0; j < nim->ny; ++j)
					{
						for (k = 0; k < nim->nx; ++k)
						{
							data[k][j][i] = (*uptr)
									* nim->scl_slope + nim->scl_inter;
							uptr++;
						}
					}
				}

			}

		}

		/// now set orientation to internal format

		// now sort out the orientation to XYZ---
		// parse the orientation string
		int xdir=1,ydir=1,zdir=1;
		if ( this->orientation.compare(3,1,"-") != 0 ) xdir = -1;
		if ( this->orientation.compare(4,1,"-") != 0 ) ydir = -1;
		if ( this->orientation.compare(5,1,"-") != 0 ) zdir = -1;

		// force the orientation to be XYZ---
		this->orientation = "XYZ---";

		// see if we have to swap x indices
		if ( xdir != 1 )
		{

			float ***xptr;
			xptr = new float**[nx];
			for ( i=0 ; i<nx ; ++i )
				xptr[i] = data[nx-i-1];

			for ( i=0 ; i<nx ; ++i )
				data[i] = xptr[i];
			delete [] xptr;
		}

		// see if we have to swap y indices
		if ( ydir != 1 )
		{

			float **yptr;
			yptr = new float*[ny];
			for ( i=0 ; i<nx ; ++i )
			{
				for ( j=0 ; j<ny ; ++j )
					yptr[j] = data[i][ny-j-1];

				for ( j=0 ; j<ny ; ++j )
					data[i][j] = yptr[j];
			}
			delete [] yptr;
		}

		// see if we have to swap z indices
		if ( zdir != 1 )
		{

			float *zptr;
			zptr = new float[nz];
			for ( i=0 ; i<nx ; ++i )
			{
				for ( j=0 ; j<ny ; ++j )
				{
					for ( k=0 ; k<nz ; ++k )
						zptr[k] = data[i][j][nz-k-1];

					for ( k=0 ; k<nz ; ++k )
						data[i][j][k] = zptr[k];
				}
			}
			delete [] zptr;
		}

		// copy to volume array
		for (i = 0; i < nz; ++i)
			for (j = 0; j < ny; ++j)
				for (k = 0; k < nx; ++k)
					VolSet[n].vox[k][j][i] = data[k][j][i];

		// pointer to data array
		VolSet[n].varray = &VolSet[n].vox[0][0][0];

		VolSet[n].CalcMinMaxAvg();

		// remove temporary array - avoiding pointer swapping errors
		data[0] = pntr1;
		data[0][0] = pntr0;

		free_matrix3D(data);
	}

	// note the file type again as it was probably tromped on by Init
	if ( nfiletype == NIFTI_FTYPE_ANALYZE )
		this->filetype = RIC_ANALYZE_FILE;
	else if ( nfiletype == NIFTI_FTYPE_NIFTI1_1 || nfiletype == NIFTI_FTYPE_NIFTI1_2 )
		this->filetype = RIC_NIFTI_FILE;

	// clean up memory
	nifti_free_NBL(&brick_list);
	nifti_image_free(nim);

	return 1;
}

/*!
Function to write volumes to a NIFTI (.nii) file or ANALYSE (.hdr/.img) file pair.
This function will write all volumes in a time series.

\param fname
name of file to write volume set to

\returns
1 on success, 0 on failure
 */
int RicVolumeSet::Write_NIFTI(string fname)
{
	// remove extension from fname if necessary
	size_t pos,l;

	if ( (pos=fname.find(".nii")) != string::npos )
	{
		l = fname.length();
		fname.erase(pos,l-pos);
	}
	if ( (pos=fname.find(".NII")) != string::npos )
	{
		l = fname.length();
		fname.erase(pos,l-pos);
	}
	if ( (pos=fname.find(".hdr")) != string::npos )
	{
		l = fname.length();
		fname.erase(pos,l-pos);
	}
	else if ( (pos=fname.find(".HDR")) != string::npos )
	{
		l = fname.length();
		fname.erase(pos,l-pos);
	}

	nifti_image      * nim;
	int xdir=1,ydir=1,zdir=1;

       if(!this->filename.empty()){
       	cerr<<"Reading header from: "<< (char*)this->filename.c_str() <<endl;
		nim = nifti_image_read((char*)this->filename.c_str(),1);

		// the the XYZ orientation
		this->orientation = convertNiftiSFormToNEMA(nim->qto_xyz);

		if ( this->orientation.compare(3,1,"-") != 0 ) xdir = -1;
		if ( this->orientation.compare(4,1,"-") != 0 ) ydir = -1;
		if ( this->orientation.compare(5,1,"-") != 0 ) zdir = -1;
	}else{cerr<<"creating new nim"<<endl;
		nim = new nifti_image;

		nim->scl_slope = 1;          // scaling parameter - slope
		nim->scl_inter = 0;          // scaling parameter - intercept
		nim->cal_min = this->VolSet[0].min; // calibration parameter, minimum
		nim->cal_max = this->VolSet[0].max; // calibration parameter, maximum
		nim->qform_code = 1;           // codes for (x,y,z) space meaning
		nim->sform_code = 0;           // codes for (x,y,z) space meaning
		nim->freq_dim = 0 ;               // indexes (1,2,3, or 0) for MRI
		nim->xyz_units = this->s_units ;   // dx,dy,dz units
		nim->time_units = this->t_units;   // dt units
		nim->qfac = 1;
		nim->intent_p1 =0;             // intent parameters
		nim->intent_p2 =0;             // intent parameters
		nim->intent_p3 =0;             // intent parameters
		nim->quatern_b = 0;
		nim->quatern_c = 0;
		nim->quatern_d = 1;
		nim->num_ext = 0;              // number of extensions in ext_list
		nim->swapsize = 0;             // swap unit in image data (might be 0)

		if ( this->filetype == RIC_ANALYZE_FILE ){
			xdir=1;
			ydir=-1;
			zdir=-1;
		}else{
			xdir=1;
			ydir=1;
			zdir=-1;
		}

		// initialize some transform matrices
		for ( int l=0 ; l<4 ; ++l )
		{
			for ( int m=0 ; m<4 ; ++m )
			{
				nim->qto_xyz.m[l][m] = 0;
				nim->qto_ijk.m[l][m] = 0;
				nim->qto_xyz.m[l][m] = 0;
				nim->qto_ijk.m[l][m] = 0;
			}
		}
		nim->qto_xyz.m[0][0] = nim->qto_ijk.m[0][0] = 1;
		nim->qto_xyz.m[1][1] = nim->qto_ijk.m[1][1] = -1;
		nim->qto_xyz.m[2][2] = nim->qto_ijk.m[2][2] = -1;
		nim->qto_xyz.m[3][3] = nim->qto_ijk.m[3][3] = 1;

              //-----

		nim->ndim = 4; 				// last dimension greater than 1 (1..7)
		nim->nx = this->nx;       // dimensions of grid array
		nim->ny = this->ny;       // dimensions of grid array
		nim->nz = this->nz;     // dimensions of grid array
		nim->nt = this->nvol;       // dimensions of grid array
		nim->nu = 0;                // dimensions of grid array
		nim->nv = 0;                // dimensions of grid array
		nim->nw = 0;                // dimensions of grid array
		nim->dim[0] = 4;            // dim[0]=ndim, dim[1]=nx, etc.
		nim->dim[1] = this->nx;
		nim->dim[2] = this->ny;
		nim->dim[3] = this->nz;
		nim->dim[4] = this->nvol;
		nim->nvox = this->totvox;   // number of voxels = nx*ny*nz*...*nw

		// figure out data type (not doing 64 bits now)
		if ( this->dtype == RIC_FLOAT ) 	// assume float
		{
			nim->datatype = DT_FLOAT;
			nim->nbyper = 4; // bytes per voxel, matches datatype
		}
		else if ( this->dtype == RIC_INTEGER )	// assume short
		{
			nim->datatype = DT_SIGNED_SHORT;
			nim->nbyper = 2; // bytes per voxel, matches datatype
		}
		else	// last chance for type - default to unsigned char
		{
			nim->datatype = DT_UNSIGNED_CHAR;
			nim->nbyper = 1; // bytes per voxel, matches datatype
		}
		nim->byteorder = 1;    // byte order on disk (2=MSB_FIRST or 1=LSB_FIRST)

		nim->dx = this->dx;          // grid spacings
		nim->dy = this->dy;          // grid spacings
		nim->dz = this->dz;          // grid spacings
		nim->dt = this->dt;          // grid spacings
		nim->du = 0;                 // grid spacings
		nim->dv = 0;                 // grid spacings
		nim->dw = 0;                 // grid spacings
		nim->pixdim[0] =4;           // pixdim[1]=dx, etc.
		nim->pixdim[1] = this->dx;
		nim->pixdim[2] = this->dy;
		nim->pixdim[3] = this->dz;
		nim->pixdim[4] = this->dt;

		/* More nifti stuff to look at in the future
		nim->phase_dim ;               // directions in dim[]/pixdim[]
		nim->slice_dim ;               // directions in dim[]/pixdim[]
		nim->slice_code  ;           // code for slice timing pattern
		nim->slice_start ;           // index for start of slices
		nim->slice_end   ;           // index for end of slices
		nim->slice_duration ;        // time between individual slices
		nim->toffset ;               // time coordinate offset
		char  intent_name[16] ;       // optional description of intent data
		char descrip[80]  ;           // optional text to describe dataset
		char aux_file[24] ;           // auxiliary filename
		*/

	}

	


	// determine file type and set file names accordingly
	char nname[1024], iname[1024];
	if ( this->filetype == RIC_ANALYZE_FILE )
	{
		nim->nifti_type = NIFTI_FTYPE_ANALYZE;
		strcpy(nname,fname.c_str());
		strcat(nname,".hdr");
		nim->fname = nname;
		strcpy(iname,fname.c_str());
		strcat(iname,".img");
		nim->iname = iname;
		nim->iname_offset = 0; // offset into iname where data starts
	}
	else
	{
		nim->nifti_type = NIFTI_FTYPE_NIFTI1_1;
		strcpy(nname,fname.c_str());
		strcat(nname,".nii");
		nim->fname = nname;
		strcpy(iname,fname.c_str());
		strcat(iname,".nii");
		nim->iname = iname;
		nim->iname_offset = 352; // offset into iname where data starts
	}

	// copy data to array of unsigned short
	if ( nim->datatype == DT_FLOAT )
	{
		float *dbuff = new float[this->totvox];
		float *bptr;	// buffer pointer to increment

		// copy to buffer
		bptr = dbuff;
		int n,i,j,k;
		for ( n=0 ; n<this->nvol ; ++n )
		{
			if ( nim->nifti_type == NIFTI_FTYPE_ANALYZE ) // swap slice order
			{		
				if ( zdir != 1 ) // swap
				{
					for ( i=this->nz-1 ; i>=0 ; --i )
					{

						if ( ydir != 1 ) // swap
						{
							for ( j=this->ny-1 ; j>=0 ; --j )
							{
								if ( xdir != 1 ) // swap
								{
									for ( k=this->nx-1 ; k>=0 ; --k )
									{
										*bptr = (float)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}else{
									for ( k=0 ; k<this->nx ; ++k )
									{
										*bptr = (float)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}
							}
						}else{
							for ( j=0 ; j<this->ny ; ++j )
							{
								if ( xdir != 1 ) // swap
								{
									for ( k=this->nx-1 ; k>=0 ; --k )
									{
										*bptr = (float)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}else{
									for ( k=0 ; k<this->nx ; ++k )
									{
										*bptr = (float)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}
							}
						}
                            	}
				}else{
					for ( i=0 ; i<this->nz ; ++i )
					{
						if ( ydir != 1 ) // swap
						{
							for ( j=this->ny-1 ; j>=0 ; --j )
							{
								if ( xdir != 1 ) // swap
								{
									for ( k=this->nx-1 ; k>=0 ; --k )
									{
										*bptr = (float)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}else{
									for ( k=0 ; k<this->nx ; ++k )
									{
										*bptr = (float)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}
							}
						}else{
							for ( j=0 ; j<this->ny ; ++j )
							{
								if ( xdir != 1 ) // swap
								{
									for ( k=this->nx-1 ; k>=0 ; --k )
									{
										*bptr = (float)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}else{
									for ( k=0 ; k<this->nx ; ++k )
									{
										*bptr = (float)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}
							}
						}

                            	}

                            }
			}
			else // just use normal slice order
			{

				if ( zdir != 1 ) // swap
				{
					for ( i=this->nz-1 ; i>=0 ; --i )
					{

						if ( ydir != 1 ) // swap
						{
							for ( j=this->ny-1 ; j>=0 ; --j )
							{
								if ( xdir != 1 ) // swap
								{
									for ( k=this->nx-1 ; k>=0 ; --k )
									{
										*bptr = (float)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}else{
									for ( k=0 ; k<this->nx ; ++k )
									{
										*bptr = (float)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}
							}
						}else{
							for ( j=0 ; j<this->ny ; ++j )
							{
								if ( xdir != 1 ) // swap
								{
									for ( k=this->nx-1 ; k>=0 ; --k )
									{
										*bptr = (float)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}else{
									for ( k=0 ; k<this->nx ; ++k )
									{
										*bptr = (float)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}
							}
						}
                            	}
				}else{
					for ( i=0 ; i<this->nz ; ++i )
					{
						if ( ydir != 1 ) // swap
						{
							for ( j=this->ny-1 ; j>=0 ; --j )
							{
								if ( xdir != 1 ) // swap
								{
									for ( k=this->nx-1 ; k>=0 ; --k )
									{
										*bptr = (float)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}else{
									for ( k=0 ; k<this->nx ; ++k )
									{
										*bptr = (float)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}
							}
						}else{
							for ( j=0 ; j<this->ny ; ++j )
							{
								if ( xdir != 1 ) // swap
								{
									for ( k=this->nx-1 ; k>=0 ; --k )
									{
										*bptr = (float)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}else{
									for ( k=0 ; k<this->nx ; ++k )
									{
										*bptr = (float)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}
							}
						}

                            	}

                            }			}
		}
		nim->data = (void*)dbuff;      // pointer to data: nbyper*nvox bytes

		// now write to file
		nifti_image_write(nim);

		// clear up memory
		delete [] dbuff;
	}
	else if ( nim->datatype == DT_SIGNED_SHORT )
	{
		short *dbuff = new short[this->totvox];
		short *bptr;	// buffer pointer to increment

		// copy to buffer
		bptr = dbuff;
		int n,i,j,k;
		for ( n=0 ; n<this->nvol ; ++n )
		{
			if ( nim->nifti_type == NIFTI_FTYPE_ANALYZE ) // swap slice order
			{
				if ( zdir != 1 ) // swap
				{
					for ( i=this->nz-1 ; i>=0 ; --i )
					{

						if ( ydir != 1 ) // swap
						{
							for ( j=this->ny-1 ; j>=0 ; --j )
							{
								if ( xdir != 1 ) // swap
								{
									for ( k=this->nx-1 ; k>=0 ; --k )
									{
										*bptr = (short)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}else{
									for ( k=0 ; k<this->nx ; ++k )
									{
										*bptr = (short)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}
							}
						}else{
							for ( j=0 ; j<this->ny ; ++j )
							{
								if ( xdir != 1 ) // swap
								{
									for ( k=this->nx-1 ; k>=0 ; --k )
									{
										*bptr = (short)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}else{
									for ( k=0 ; k<this->nx ; ++k )
									{
										*bptr = (short)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}
							}
						}
                            	}
				}else{
					for ( i=0 ; i<this->nz ; ++i )
					{
						if ( ydir != 1 ) // swap
						{
							for ( j=this->ny-1 ; j>=0 ; --j )
							{
								if ( xdir != 1 ) // swap
								{
									for ( k=this->nx-1 ; k>=0 ; --k )
									{
										*bptr = (short)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}else{
									for ( k=0 ; k<this->nx ; ++k )
									{
										*bptr = (short)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}
							}
						}else{
							for ( j=0 ; j<this->ny ; ++j )
							{
								if ( xdir != 1 ) // swap
								{
									for ( k=this->nx-1 ; k>=0 ; --k )
									{
										*bptr = (short)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}else{
									for ( k=0 ; k<this->nx ; ++k )
									{
										*bptr = (short)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}
							}
						}

                            	}

                            }			
			}
			else // just use normal slice order
			{
				if ( zdir != 1 ) // swap
				{
					for ( i=this->nz-1 ; i>=0 ; --i )
					{

						if ( ydir != 1 ) // swap
						{
							for ( j=this->ny-1 ; j>=0 ; --j )
							{
								if ( xdir != 1 ) // swap
								{
									for ( k=this->nx-1 ; k>=0 ; --k )
									{
										*bptr = (short)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}else{
									for ( k=0 ; k<this->nx ; ++k )
									{
										*bptr = (short)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}
							}
						}else{
							for ( j=0 ; j<this->ny ; ++j )
							{
								if ( xdir != 1 ) // swap
								{
									for ( k=this->nx-1 ; k>=0 ; --k )
									{
										*bptr = (short)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}else{
									for ( k=0 ; k<this->nx ; ++k )
									{
										*bptr = (short)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}
							}
						}
                            	}
				}else{
					for ( i=0 ; i<this->nz ; ++i )
					{
						if ( ydir != 1 ) // swap
						{
							for ( j=this->ny-1 ; j>=0 ; --j )
							{
								if ( xdir != 1 ) // swap
								{
									for ( k=this->nx-1 ; k>=0 ; --k )
									{
										*bptr = (short)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}else{
									for ( k=0 ; k<this->nx ; ++k )
									{
										*bptr = (short)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}
							}
						}else{
							for ( j=0 ; j<this->ny ; ++j )
							{
								if ( xdir != 1 ) // swap
								{
									for ( k=this->nx-1 ; k>=0 ; --k )
									{
										*bptr = (short)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}else{
									for ( k=0 ; k<this->nx ; ++k )
									{
										*bptr = (short)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}
							}
						}

                            	}

                            }			
			}
		}
		nim->data = dbuff;      // pointer to data: nbyper*nvox bytes

		// now write to file
		nifti_image_write(nim);

		// clear up memory
		delete [] dbuff;
	}
	else //***** note, probably need scale factors for unsigned char
	{
		unsigned char *dbuff = new unsigned char[this->totvox];
		unsigned char *bptr;	// buffer pointer to increment

		// copy to buffer
		bptr = dbuff;
		int n,i,j,k;
		for ( n=0 ; n<this->nvol ; ++n )
		{
			if ( nim->nifti_type == NIFTI_FTYPE_ANALYZE ) // swap slice order
			{
				if ( zdir != 1 ) // swap
				{
					for ( i=this->nz-1 ; i>=0 ; --i )
					{

						if ( ydir != 1 ) // swap
						{
							for ( j=this->ny-1 ; j>=0 ; --j )
							{
								if ( xdir != 1 ) // swap
								{
									for ( k=this->nx-1 ; k>=0 ; --k )
									{
										*bptr = (unsigned char)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}else{
									for ( k=0 ; k<this->nx ; ++k )
									{
										*bptr = (unsigned char)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}
							}
						}else{
							for ( j=0 ; j<this->ny ; ++j )
							{
								if ( xdir != 1 ) // swap
								{
									for ( k=this->nx-1 ; k>=0 ; --k )
									{
										*bptr = (unsigned char)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}else{
									for ( k=0 ; k<this->nx ; ++k )
									{
										*bptr = (unsigned char)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}
							}
						}
                            	}
				}else{
					for ( i=0 ; i<this->nz ; ++i )
					{
						if ( ydir != 1 ) // swap
						{
							for ( j=this->ny-1 ; j>=0 ; --j )
							{
								if ( xdir != 1 ) // swap
								{
									for ( k=this->nx-1 ; k>=0 ; --k )
									{
										*bptr = (unsigned char)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}else{
									for ( k=0 ; k<this->nx ; ++k )
									{
										*bptr = (unsigned char)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}
							}
						}else{
							for ( j=0 ; j<this->ny ; ++j )
							{
								if ( xdir != 1 ) // swap
								{
									for ( k=this->nx-1 ; k>=0 ; --k )
									{
										*bptr = (unsigned char)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}else{
									for ( k=0 ; k<this->nx ; ++k )
									{
										*bptr = (unsigned char)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}
							}
						}

                            	}

                            }
			}
			else // just use normal slice order
			{
				if ( zdir != 1 ) // swap
				{
					for ( i=this->nz-1 ; i>=0 ; --i )
					{

						if ( ydir != 1 ) // swap
						{
							for ( j=this->ny-1 ; j>=0 ; --j )
							{
								if ( xdir != 1 ) // swap
								{
									for ( k=this->nx-1 ; k>=0 ; --k )
									{
										*bptr = (unsigned char)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}else{
									for ( k=0 ; k<this->nx ; ++k )
									{
										*bptr = (unsigned char)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}
							}
						}else{
							for ( j=0 ; j<this->ny ; ++j )
							{
								if ( xdir != 1 ) // swap
								{
									for ( k=this->nx-1 ; k>=0 ; --k )
									{
										*bptr = (unsigned char)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}else{
									for ( k=0 ; k<this->nx ; ++k )
									{
										*bptr = (unsigned char)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}
							}
						}
                            	}
				}else{
					for ( i=0 ; i<this->nz ; ++i )
					{
						if ( ydir != 1 ) // swap
						{
							for ( j=this->ny-1 ; j>=0 ; --j )
							{
								if ( xdir != 1 ) // swap
								{
									for ( k=this->nx-1 ; k>=0 ; --k )
									{
										*bptr = (unsigned char)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}else{
									for ( k=0 ; k<this->nx ; ++k )
									{
										*bptr = (unsigned char)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}
							}
						}else{
							for ( j=0 ; j<this->ny ; ++j )
							{
								if ( xdir != 1 ) // swap
								{
									for ( k=this->nx-1 ; k>=0 ; --k )
									{
										*bptr = (unsigned char)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}else{
									for ( k=0 ; k<this->nx ; ++k )
									{
										*bptr = (unsigned char)(this->VolSet[n].vox[k][j][i]);
										bptr++;
									}
								}
							}
						}

                            	}

                            }

			}
		}
		nim->data = dbuff;      // pointer to data: nbyper*nvox bytes

		// now write to file
		nifti_image_write(nim);

		// clear up memory
		delete [] dbuff;
	}

	delete nim;

	return 1;

}

/*!
Function to read a volume from a GIS (.dim/.ima) file pair.

\param filename
Name of file containing volume set

\returns
1 on success, 0 on failure
 */
int RicVolumeSet::Read_GIS(string fname)
{
	// remove extension from filename if necessary
	size_t pos,l;

	if ( (pos=fname.find(".dim")) != string::npos )
	{
		l = fname.length();
		fname.erase(pos,l-pos);
	}
	else if ( (pos=fname.find(".ima")) != string::npos )
	{
		l = fname.length();
		fname.erase(pos,l-pos);
	}

	// note the file type
	this->filetype = RIC_GIS_FILE;

	// add dim to the header filename
	string dimname = fname+".dim";

	// parse the header
	FILE *fptr;
	fptr = fopen(dimname.c_str(),"r");
	if ( fptr == NULL )
	{
		cerr << "RicVolumeSet::Read_GIS - cannot open file: " << dimname << endl;
		return 0;
	}

	char line[128], dum[128], bo[128], type[128], imode[128];

	// first get the size of the data
	fscanf(fptr,"%d %d %d %d",&nx,&ny,&nz,&nvol);

	// simple minded error check on array size
	if ( nvol < 1 || nx < 1 || ny < 1 || nz < 1 ) return 0;

	// now get the options
	int bswap=0;
	int binary = 1;
	float x=1,y=1,z=1,t=0;	// increments for x,y,z
	while ( fgets(line,100,fptr) != NULL )
	{
		if ( strncmp(line,"-bo",3 ) == 0 ) // byte order
		{
			sscanf(line,"%s %s",dum,bo);
			if ( strncmp(bo,"ABCD",4) == 0 ) // big endian
			{
				bswap = 1;
			}
			else if ( strncmp(bo,"DCBA",4) == 0 ) // big endian
			{
				bswap = 0;
			}
			else
			{
				cerr << "RicVolumeSet::Read_GIS - unknown byte order " << bo << endl;
				return 0;
			}
		}
		else if ( strncmp(line,"-type",5) == 0 )	// type
		{
			sscanf(line,"%s %s",dum,type);
			if ( strncmp(type,"S16",3) == 0  || strncmp(type,"s16",3) == 0 ) // unsigned 16
			{
				dtype = RIC_INTEGER;
			}
			else
			{
				cerr << "RicVolumeSet::Read_GIS - unknown data type " << type << endl;
				return 0;
			}
		}
		else if ( strncmp(line,"-om",3) == 0 ) // binary or ascii mode
		{
			sscanf(line,"%s %s",dum,imode);
			if ( strncmp(imode,"binar",5) == 0 )
			{
				binary = 1;
			}
			else
			{
				cerr << "RicVolumeSet::Read_GIS - unknown image mode " << imode << endl;
				return 0;
			}
		}

		// look for voxel dimensions - can be on same line
		char *sptr;
		if ( (sptr = strstr(line,"-dx") ) )
		{
			sscanf(sptr,"%s %f",dum,&x);
		}
		if ( (sptr = strstr(line,"-dy") ) )
		{
			sscanf(sptr,"%s %f",dum,&y);
		}
		if ( (sptr = strstr(line,"-dz") ) )
		{
			sscanf(sptr,"%s %f",dum,&z);
		}
		if ( (sptr = strstr(line,"-dt") ) )
		{
			sscanf(sptr,"%s %f",dum,&t);
		}
	}

	fclose(fptr);

	// output file name
	string dataname = fname+".ima";

	// now write the data file
	fptr = fopen(dataname.c_str(),"rb");

	if ( fptr == NULL )
	{
		cerr << "RicVolumeSet::Read_GIS - cannot open file: " << dataname << endl;
		return 0;
	}

	// allocate new volumes
	this->Init(nx, ny, nz, nvol);
	this->dx = x;
	this->dy = y;
	this->dz = z;
	this->dt = t;
	this->filetype = RIC_GIS_FILE; // init kills this

	int i,j,k,n;
	if ( dtype == RIC_INTEGER )
	{
		short val;
		unsigned short *rowarray;
		rowarray = new unsigned short[nx]; // array to read a row at a time
		for ( n=0 ; n< nvol ; ++n )
		{
			// set spacing for each volume
			this->VolSet[n].dx = x;
			this->VolSet[n].dy = y;
			this->VolSet[n].dz = z;

			for ( i=0 ; i<nz ; ++i )
			{
				for ( j=0 ; j< ny ; ++j )
				{
					fread(rowarray, sizeof(short), nx, fptr);
					for ( k=0 ; k<nx ; ++k )
					{

						if (bswap) val = (rowarray[k] >> 8) | (rowarray[k] << 8);
						else val = (short)rowarray[k];
						this->VolSet[n].vox[k][j][i] = val;
					}
				}
			}
		}
		delete rowarray;
	}

	fclose(fptr);

	return 1;

}

/*!
Function to write a volume to a GIS (.dim/.ima) file pair. It outputs
little endian data and assumes Intel.

\param filename
name of file to write volume set to

\returns
1 on success, 0 on failure
 */
int RicVolumeSet::Write_GIS(string fname)
{
	// remove extension from filename if necessary
	size_t pos,l;

	if ( (pos=fname.find(".dim")) != string::npos )
	{
		l = fname.length();
		fname.erase(pos,l-pos);
	}
	else if ( (pos=filename.find(".ima")) != string::npos )
	{
		l = fname.length();
		fname.erase(pos,l-pos);
	}

	// make sure there is something to write
	if ( VolSet[0].nz == 0 ) return 0;

	// add dim to the header filename
	string dimname = fname+".dim";

	// write the header
	FILE *fptr;
	fptr = fopen(dimname.c_str(),"w");

	fprintf(fptr,"%d %d %d %d\n",nx, ny, nz, nvol);
	fprintf(fptr,"-type S16\n");
	fprintf(fptr,"-dx %f -dy %f -dz %f -dt %f\n",dx,dy,dz,dt);
	fprintf(fptr,"-bo DCBA\n"); // assuming intel and Little_endian
	fprintf(fptr,"-om binar\n");

	fclose(fptr);

	// output file name
	string dataname = fname+".ima";

	// now write the data file
	fptr = fopen(dataname.c_str(),"wb");

	int i,j,k,n;
	short val;
	short *rowarray;
//	int bswap=0;	// if 1 then swap byte order
	rowarray = new short[nx]; // write a row at a time
	for ( n=0 ; n<nvol ; ++n )
	{
		for ( k=0 ; k<nz ; ++k )
		{
			for ( j=0 ; j< ny ; ++j )
			{
				for ( i=0 ; i<nx ; ++i )
				{
					val = (short)(VolSet[n].vox[i][j][k]+0.5); // round to nearest integer
					//if (bswap) rowarray[k] = (val >> 8) | (val << 8);
					rowarray[i] = val;
				}
				fwrite(rowarray, sizeof(short), nx, fptr);
			}
		}
	}

	delete rowarray;
	fclose(fptr);

	return 1;
}
