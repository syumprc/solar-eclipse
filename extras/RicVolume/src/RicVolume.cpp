// ----------------------- RicVolume.cpp -------------------------
/*! \mainpage
The library libRicVolume.a contains two classes for manupulating volume data.

The basic volume class, RicVolume, has basic maniuplations of 3D volumes. The
voxel data is stored in a three dimensional array where the origin, index [0][0][0],
is located at the left, anterior, superior. In other words, the array start at the
top, right, front of the head. This corresponds to a NEMA orientation string of 
XYZ---. Each voxel consists of a floating point intensity value and an integer flag 
that can be used for any purpose. The are overloads of common arithmatic operators 
so whole volumes can be operated on.

The voxel values can be accessed in two different manners. The functions get_at_index
and set_at_index allow accessing a voxel at a time using the x,y and z index values. 
In addition, the voxel array can be directly address by vox[ix][iy][iz].vv which is
more efficient.

The RicVolumeSet class is used to encapsulate sets of volumes such as a time series
of volumes. This class is responsible for reading and writing data files and keeping
track of input and output file orientation. When a file is read, it is converted to 
the internal XYZ--- orientation.

*/

/*! \file
Source file for the RicVolume class. The member functions allow
manipulation of voxel data.

Copyright (C) 2007 by Bill Rogers - UTHSCSA
rogers@uthscsa.edu
*/
using namespace std;

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "RicVolume.h"
#include "RicVolumeSet.h"

/*!
Constructor for empty volume
 */
RicVolume::RicVolume()
{
	vox = NULL;
	varray = NULL;
	nz = 0;
	ny = 0;
	nx = 0;
	nvox = 0;
	dx = dy  = dz = 1;
	xoffset = yoffset = zoffset = 0;
	min = 0;
	max = 1000;
	avg = 500;
};

/*!
Constructor to make a volume with all voxel values set to zero

\param nx
Number of x values in volume (columns in a row)

\param ny
Number of y values in volume (rows in a slice)

\param nz
Number of z values in volume (slices)
 */
RicVolume::RicVolume(int nx, int ny, int nz)
{
	vox = NULL;
	varray = NULL;
	dx = dy  = dz = 1;
	xoffset = yoffset = zoffset = 0;
	min = 0;
	max = 1000;
	avg = 500;
	this->Init(nx,ny,nz,0);
};

/*!
Constructor to make a volume with all voxel values set to the passed value

\param nx
Number of x values in volume (columns in a row)

\param ny
Number of y values in volume (rows in a slice)

\param nz
Number of z values in volume (slices)

\param val
Value to set all voxel intensity values to
 */
RicVolume::RicVolume(int nx, int ny, int nz, float val)
{
	vox = NULL;
	varray = NULL;
	dx = dy  = dz = 1;
	xoffset = yoffset = zoffset = 0;
	min = 0;
	max = 1000;
	avg = 500;
	this->Init(nx,ny,nz,val);
};

/*!
Constructor for volume to be read in from file. If more than one volume is in
the file then all be the first will be discarded.

\param filename
Name of file containing a volume
 */
RicVolume::RicVolume(string filename)
{
	vox = NULL;
	varray = NULL;
	nz = 0;
	ny = 0;
	nx = 0;
	nvox = 0;
	dx = dy  = dz = 1;
	xoffset = yoffset = zoffset = 0;
	
	RicVolumeSet Vset(filename);
	
	// return if no volumes read
	if ( Vset.VolSet == NULL) return;
	
	// since we must have at least one volume, copy it to a local structures
	this->Copy(Vset.VolSet[0]);
	
	// get the volume info
	this->CalcMinMaxAvg();
	
	
};

/*!
 * Destructor for RicVolume
 */
RicVolume::~RicVolume()
{
	if ( vox ) free_matrix3D(vox);
}

/*!
Initialization function that allocates memory for a volume and fills
it with zeros

\param nx
Number of x values in volume (columns in a row)

\param ny
Number of y values in volume (rows in a slice)

\param nz
Number of z values in volume (slices)

\returns
1 on success, 0 on failure
 */
int RicVolume::Init(int nx, int ny, int nz)
{
	return this->Init(nx,ny,nz,0);
}

/*!
Initialization function that allocates memory for a volume and fills
it with the passed value
 
\param numx
Number of x values in volume (columns in a row)

\param numy
Number of y values in volume (rows in a slice)

\param numz
Number of z values in volume (slices)

\param val
Value to assign to all voxel intensity values

\returns
1 on success, 0 on failure
*/
int RicVolume::Init(int numx, int numy, int numz, float val)
{
	// free memory if allocated
	if ( vox ) free_matrix3D(vox);
	
	nx = numx;
	ny = numy;
	nz = numz;
	nvox = nz*ny*nx;
	min = max = avg = val;
	
	// allocate a new voxel array
	int status=matrix3D(&vox,nx,ny,nz);
	if ( !status ) return 0;
	
	// assign point to first element of array
	varray = &vox[0][0][0];
	
	// fill volume with zeros
	for ( int i=0 ; i<nvox ; ++i ) varray[i] = val;
	
	return 1;
}

////////////////////////////// Overloaded Operators ////////////////////////////

/*!
Operator > overload - produces a new volume with voxels set to
1 where the original volume is greater than the passed threshold
value and 0 elsewhere

\param thresh - threshold value
\returns - pointer to new volume
 */
RicVolume *RicVolume::operator > (const float thresh)
{
	// create a new volume of the same size
	RicVolume *newVol = new RicVolume(this->nx,this->ny,this->nz);
	 
	// fill the volume with 1's where the voxel values are greater
	// than thresh and with 0's elsewhere
	for (int i=0; i<nvox; i++)
				newVol->varray[i] = this->varray[i] > thresh;   
	return (newVol);
}

/*!
Operator += overload - add a constant value to every voxel in a volume

\param con - value to be added to every voxel
\returns - volume
 */
RicVolume& RicVolume::operator+=(const float con)
{
	for (int i=0; i<nvox; i++) this->varray[i] += con;   
	return (*this);
}

/*!
Operator -= overload - subtract a constant value from every voxel in a volume

\param con - value to be subtracted from every voxel
\returns - volume
 */
RicVolume& RicVolume::operator-=(const float con)

{
	for (int i=0; i<nvox; i++) this->varray[i] -= con;   
	return (*this);
}

/*!
Operater *= overload - multiply every voxel in a volume by a constant value

\param con - constant value to multiply every voxel by
\returns - volume
 */
RicVolume& RicVolume::operator*=(const float con)
{
	for (int i=0; i<nvox; i++) this->varray[i] *= con;   
	return (*this);
}

/*!
Operater /= overload - divide every voxel in a volume by a constant value

\param con - value to divide every voxel by
\returns - volume
 */
RicVolume& RicVolume::operator/=(const float con)
{
	// just return if con == 0 - cannot divide by zero
	if ( con == 0 ) return (*this);
	
	for (int i=0; i<nvox; i++) this->varray[i] /= con;   
	return (*this);
}

/*!
Operator + overload - add each voxel from one volume to another

\param rhs - volume to be added to this volume
\returns - volume
 */
RicVolume &RicVolume :: operator + (RicVolume &rhs)
{
	for (int i=0; i<nvox; i++)
		this->varray[i] += rhs.varray[i];   
	return (*this);
}

/*!
Operator - overload - subtract voxels from the passed volume from this one

\param rhs - volume to be subtracted from this volume
\returns - volume
 */
RicVolume &RicVolume :: operator - (RicVolume &rhs)
{
	for (int i=0; i<nvox; i++)
		this->varray[i] -= rhs.varray[i];   
	return (*this);
}

/*!
Operator * overload - multiply voxels from the passed volume and this one

\param rhs - volume to multiply by this volume
\returns - volume
 */
RicVolume &RicVolume :: operator * (RicVolume &rhs)
{
	for (int i=0; i<nvox; i++)
		this->varray[i] *= rhs.varray[i];   
	return (*this);
}

/*!
Operator / overload - divide voxels from the passed volume from this one
There is no checking for divide by zero so use cautiously

\param - rhs - volume to divide this one by
\returns - volume
 */
RicVolume &RicVolume :: operator / (RicVolume &rhs)
{
	for (int i=0; i<nvox; i++)
		this->varray[i] /= rhs.varray[i];   
	return (*this);
}

/////////////////////////// Member Functions ///////////////////////////////

/*!
This function calculates the min, max, and average voxel values for a volume
 */
void RicVolume::CalcMinMaxAvg(void)
{
	min = 999999;
	max = -999999;
	avg = 0;
	
	// check every voxel
	for ( int i=0 ; i<nvox ; ++i )
	{
		if ( varray[i] < min ) min = varray[i];
		if ( varray[i] > max ) max = varray[i];
		avg += varray[i];
	}
	
	avg /= nvox;
}

/*!
This function fills a volume with the passed value

\param fillval
Intensity value to store in all voxel locations
 */
void RicVolume::Fill(float fillval)
{
	for (int i=0; i<nvox; i++)	this->varray[i] = fillval;
}

/*!
This function sets all voxels in a volume that have a value that is
above the threshold value to the passed fill value

\param threshold
Threshold value for voxels

\param fillval
Intensity value to store in all voxel locations
 */
void RicVolume::ThresholdAbove(float threshold, float fillval)
{
	for (int i=0; i<nvox; i++)
		if ( this->varray[i] > threshold ) this->varray[i] = fillval;   
}

/*!
This function sets all voxels in a volume that have a value that is
below the threshold value to the passed fill value

\param threshold
Threshold value for voxels

\param fillval
Intensity value to store in all voxel locations
 */
void RicVolume::ThresholdBelow(float threshold, float fillval)
{
	for (int i=0; i<nvox; i++)
				if ( this->varray[i] < threshold ) this->varray[i] = fillval;   
}

/*!
Copy from the source volume duplicating all setting and values.

\param source
Source volume to copy from
 */
void RicVolume::Copy(RicVolume &source)
{
	int i;
 
	if ( vox ) free_matrix3D(vox);
	nx = source.nx;
	ny = source.ny;
	nz = source.nz;
	nvox = source.nvox;
	dx = source.dx;
	dy = source.dy;
	dz = source.dz;
	xoffset = source.xoffset;
	yoffset = source.yoffset;
	zoffset = source.zoffset;
	min = source.min;
	max = source.max;
	avg = source.avg;

	matrix3D(&vox,nx,ny,nz);
	this->varray = &vox[0][0][0];
	for ( i=0; i<nvox; i++)	this->varray[i] = source.varray[i];   
}

/*!
Fills f6_buff with the values of the 6 non-diagonal adjacent pixel.  They
will be ordered as follows: 
Positive X neighbor, Negative X neighbor, +Y, -Y, +Z, -Z.
If the input indices are such that all values will not be in bounds then

\param ix
x index

\param iy
y index

\param iz
z index

\param buff6
buffer to store neighboring voxel values in

\returns
zero is returned, else 1 is returned on success.
 */
int RicVolume::get_6_neighbors(int ix, int iy, int iz, float *buff6)
{
	// check to see that indices are in bounds
	if ( ix < 1 || iy < 1 || iz < 1 ) return 0;
	if ( ix >= nx-1 || iy >= ny-1 || iz >= nz-1 ) return 0;
	
	// assign the neighbors
	buff6[0] = vox[iz][iy][ix+1];
	buff6[1] = vox[iz][iy][ix-1];
	buff6[2] = vox[iz][iy+1][ix];
	buff6[3] = vox[iz][iy-1][ix];
	buff6[4] = vox[iz+1][iy][ix];
	buff6[5] = vox[iz-1][iy][ix];
	return 1;
}

/*!
Fills d18_buff ordered as follows: slice iz-1 point (ix-1,iy) first then 
clockwise around but skipping the corners and lastly the center point of iz-1. 
Then slice iz - will all values except (ix,iy,iz).  Then iz+1 same order as iz-1. 
If the input indices are such that all values will not be in bounds then
zero is returned, else 1 is returned on success.

\param ix
x index

\param iy
y index

\param iz
z index

\param buff18
buffer to store neighboring voxel values in

\returns
zero is returned, else 1 is returned on success.
 */
int RicVolume::get_18_neighbors(int ix, int iy, int iz, float *buff18)
{
	// check to see that indices are in bounds
	if ( ix < 1 || iy < 1 || iz < 1 ) return 0;
	if ( ix >= nx-1 || iy >= ny-1 || iz >= nz-1 ) return 0;

	// bottom slice
	buff18[0] = vox[ix-1][iy][iz-1];
	buff18[1] = vox[ix][iy+1][iz-1];
	buff18[2] = vox[ix+1][iy][iz-1];
	buff18[3] = vox[ix][iy-1][iz-1];
	buff18[4] = vox[ix][iy][iz-1];

	// middle slice
	buff18[5] = vox[ix-1][iy-1][iz];
	buff18[6] = vox[ix-1][iy][iz];
	buff18[7] = vox[ix-1][iy+1][iz];
	buff18[8] = vox[ix][iy+1][iz];
	buff18[9] = vox[ix+1][iy+1][iz];
	buff18[10] = vox[ix+1][iy][iz];
	buff18[11] = vox[ix+1][iy-1][iz];
	buff18[12] = vox[ix][iy-1][iz];

	// top slice
	buff18[13] = vox[ix-1][iy][iz+1];
	buff18[14] = vox[ix][iy+1][iz+1];
	buff18[15] = vox[ix+1][iy][iz+1];
	buff18[16] = vox[ix][iy-1][iz+1];
	buff18[17] = vox[ix][iy][iz+1];

	return 1;
}

/*!
Fills d26_buff ordered as follows: slice iz-1 point (x-1,y-1) first then 
clockwise around and lastly the center point of iz-1. Then slice iz - same 
order but not (ix,iy,iz).  Then iz+1 same order as iz-1. 
If the input indicies are such that all values will not be in bounds then
zero is returned, else 1 is returned on success.

\param ix
x index

\param iy
y index

\param iz
z index

\param buff26
buffer to store neighboring voxel values in

\returns
zero is returned, else 1 is returned on success.
*/
int RicVolume::get_26_neighbors(int ix, int iy, int iz, float *buff26)
{
	// check to see that indices are in bounds
	if ( ix < 1 || iy < 1 || iz < 1 ) return 0;
	if ( ix >= nx-1 || iy >= ny-1 || iz >= nz-1 ) return 0;

	// bottom slice
	buff26[0] = vox[ix-1][iy-1][iz-1];
	buff26[1] = vox[ix-1][iy][iz-1];
	buff26[2] = vox[ix-1][iy+1][iz-1];
	buff26[3] = vox[ix][iy+1][iz-1];
	buff26[4] = vox[ix+1][iy+1][iz-1];
	buff26[5] = vox[ix+1][iy][iz-1];
	buff26[6] = vox[ix+1][iy-1][iz-1];
	buff26[7] = vox[ix][iy-1][iz-1];
	buff26[8] = vox[ix][iy][iz-1];

	// middle slice
	buff26[9] = vox[ix-1][iy-1][iz];
	buff26[10] = vox[ix-1][iy][iz];
	buff26[11] = vox[ix-1][iy+1][iz];
	buff26[12] = vox[ix][iy+1][iz];
	buff26[13] = vox[ix+1][iy+1][iz];
	buff26[14] = vox[ix+1][iy][iz];
	buff26[15] = vox[ix+1][iy-1][iz];
	buff26[16] = vox[ix][iy-1][iz];

	// top slice
	buff26[17] = vox[ix-1][iy-1][iz+1];
	buff26[18] = vox[ix-1][iy][iz+1];
	buff26[10] = vox[ix-1][iy+1][iz+1];
	buff26[20] = vox[ix][iy+1][iz+1];
	buff26[21] = vox[ix+1][iy+1][iz+1];
	buff26[22] = vox[ix+1][iy][iz+1];
	buff26[23] = vox[ix+1][iy-1][iz+1];
	buff26[24] = vox[ix][iy-1][iz+1];
	buff26[25] = vox[ix][iy][iz+1];

	return 1;
}

/*!
Use trilinear interpolation to calculate the intensity level at 
the passed x,y,z volume location. Values are ranged checked to make
sure the x,y,z location is in the volume. The interpolated 
intensity level is returned on success. Zero is return if the location 
is out of range.
Note - The passed x,y,z location is in real space and is converted to
voxel space by dividing by the x,y,z voxel size.

\param rx - x volume location
\param ry - y volume location
\param rz - z volume location
\returns - interpolated voxel value
 */
float RicVolume::interp3D(float rx, float ry, float rz)
{
	// scale
	float x = rx/dx;
	float y = ry/dy;
	float z = rz/dz;
	
	// range check input
	if (x<0 || x>=nx || y<0 || y>=ny || z<0 || z>=nz ) return 0;

	int x1,x2,y1,y2,z1,z2;
	x1 = (int)x;
	x2 = x1 + 1;
	y1 = (int)y;
	y2 = y1 + 1;
	z1 = (int)z;
	z2 = z1 + 1;

	// range check high index values
	if ( x2>=nx || y2>=ny || z2>=nz ) return 0;
	
	float answer, xfrac, xfrac2, yfrac, yfrac2, zfrac, zfrac2;
	xfrac = x - (float)x1;
	yfrac = y - (float)y1;
	zfrac = z - (float)z1;
	xfrac2 = 1.0f - xfrac;
	yfrac2 = 1.0f - yfrac;
	zfrac2 = 1.0f - zfrac;
	
	float values[8];
	values[0] = this->vox[x1][y1][z1];
	values[1] = this->vox[x1][y1][z2];
	values[2] = this->vox[x1][y2][z1];
	values[3] = this->vox[x1][y2][z2];
	values[4] = this->vox[x2][y1][z1];
	values[5] = this->vox[x2][y1][z2];
	values[6] = this->vox[x2][y2][z1];
	values[7] = this->vox[x2][y2][z2];

	float v01, v23, v45, v67, v0123, v4567;
	v01 = values[0] * xfrac2 + values[1] * xfrac;
	v23 = values[2] * xfrac2 + values[3] * xfrac;
	v45 = values[4] * xfrac2 + values[5] * xfrac;
	v67 = values[6] * xfrac2 + values[7] * xfrac;
	v0123 = v01 * yfrac2 + v23 * yfrac;
	v4567 = v45 * yfrac2 + v67 * yfrac;
	answer=v0123 * zfrac2 + v4567 * zfrac;

	return (answer);

}

/*!
Where ix, iy, and iz are column, row, and slice indices. The voxel value
at those indices is returned.

\param ix - x index
\param iy - y index
\param iz - z index
\returns - voxel value at location of indices
 */
float RicVolume::get_at_index(int ix, int iy, int iz, int flag)
{
	return this->vox[ix][iy][iz];
}

/*!
Where ix, iy, and iz are column, row, and slice indices. The voxel value
at those indices is returned. This version checks to see that the indices are
in the range of the volume.

\param ix - x index
\param iy - y index
\param iz - z index
\returns - voxel value at location of indices or 0 if out of range
 */
float RicVolume::range_checked_get_at_index(int ix, int iy, int iz)
{
	if ( ix < 0 || ix >= nx || iy < 0 || iy >= ny || iz < 0 || iz >= nz )
		return 0;
	
	return this->vox[ix][iy][iz];
}

/*!
Where ix, iy, and iz are column, row, and slice indices. The voxel value
at those indices is set to the passed value. If the indices are out of 
range then zero is returned. 1 is returned on success.

\param ix - x index
\param iy - y index
\param iz - z index
\returns - 1 on success or 0 if out of range
*/
int RicVolume::set_at_index(int ix, int iy, int iz, float value)
{
	if ( ix < 0 || ix >= nx || iy < 0 || iy >= ny || iz < 0 || iz >= nz )
		return 0;
	
	this->vox[ix][iy][iz] = value;
	return 1;
}

/*!
Where rx, ry, and rz make up a vertex in real space. The voxel value
nearest that location is set to the passed value. If the indices are out of
range then zero is returned. 1 is returned on success.

\param rx - x location
\param ry - y location
\param rz - z location
\returns - 1 on success or 0 if out of range
 */
int RicVolume::set_at_nearest(float rx, float ry, float rz, float value)
{
	// scale
	float x = rx/dx;
	float y = ry/dy;
	float z = rz/dz;
	
	// round to nearest integer
	int ix = (int)round(x);
	int iy = (int)round(y);
	int iz = (int)round(z);
	
	// check to see if in range
	if ( ix < 0 || ix >= nx || iy < 0 || iy >= ny || iz < 0 || iz >= nz )
		return 0;
	
	// set the proper voxel
	this->vox[ix][iy][iz] = value;
	return 1;
}

/*!
 * Function to read a volume from a file. If there is more than one volume in the
 * file then all be the first will be discarded. The file extension will determine
 * what kind of file it is.
 * 
 *\param filename
 * Name of file to read volume from
 * 
 * \returns
 * 1 on success, 0 on failure
 */
int RicVolume::Read(string filename)
{
	// free voxel memory if allocated
	if ( vox )
	{
		free_matrix(vox);
	}
	vox = NULL;
	nz = 0;
	ny = 0;
	nx = 0;
	nvox = 0;
	dx = dy  = dz = 1;
	xoffset = yoffset = zoffset = 0;
	
	RicVolumeSet Vset(filename);
	
	// return if no volumes read
	if ( Vset.VolSet == NULL) return 0;
	
	// since we must have at least one volume, copy it to a local structures
	this->Copy(Vset.VolSet[0]);
	
	// get the volume info
	this->CalcMinMaxAvg();

	return 1;
}

/*!
 * Function to write a volume to a file. The file extension will determine what kind
 * of file is written (i.e. dim, des, nii). If no extension is given then the NEMA
 * format will be used.
 * 
 * \param filename
 * Name of file to write to
 * 
 * \returns
 * 1 on success, 0 on failure
 */
int RicVolume::Write(string filename)
{
	// create a one volume volume set
	RicVolumeSet Vset(nx,ny,nz,1);
	
	// copy values to volume set
	Vset.dx = Vset.VolSet[0].dx = dx;
	Vset.dy = Vset.VolSet[0].dy = dy;
	Vset.dz = Vset.VolSet[0].dz = dz;
	Vset.xoffset = Vset.VolSet[0].xoffset = xoffset;
	Vset.yoffset = Vset.VolSet[0].yoffset = yoffset;
	Vset.zoffset = Vset.VolSet[0].zoffset = zoffset;
	Vset.VolSet[0].min = min;
	Vset.VolSet[0].max = max;
	
	
	// copy to the first volume
	int i,j,k;
	for ( i=0; i<nx; i++)
		for ( j=0; j<ny; j++)   
			for ( k=0; k<nz; k++)
				Vset.VolSet[0].vox[i][j][k] = vox[i][j][k];   
	
	return Vset.Write(filename);
	
}

/*!
 * Function to write a volume to a NEMA file. .
 * 
 * \param filename
 * Name of file to write to
 * 
 * \returns
 * 1 on success, 0 on failure
 */
int RicVolume::Write_NEMA(string filename)
{
	// create a one volume volume set
	RicVolumeSet Vset(nx,ny,nz,1);
	
	// copy values to volume set
	Vset.dx = Vset.VolSet[0].dx = dx;
	Vset.dy = Vset.VolSet[0].dy = dy;
	Vset.dz = Vset.VolSet[0].dz = dz;
	Vset.xoffset = Vset.VolSet[0].xoffset = xoffset;
	Vset.yoffset = Vset.VolSet[0].yoffset = yoffset;
	Vset.zoffset = Vset.VolSet[0].zoffset = zoffset;
	Vset.VolSet[0].min = min;
	Vset.VolSet[0].max = max;
	
	
	// copy to the first volume
	int i,j,k;
	for ( i=0; i<nx; i++)
		for ( j=0; j<ny; j++)   
			for ( k=0; k<nz; k++)
				Vset.VolSet[0].vox[i][j][k] = vox[i][j][k];   
	
	return Vset.Write_NEMA(filename);
	
}

/*!
 * Function to write a volume to a NEMA file. .
 * 
 * \param filename
 * Name of file to write to
 * 
 * \param orient_str - NEMA orientation string
 * \returns
 * 1 on success, 0 on failure
 */
int RicVolume::Write_NEMA(string filename, string orient_str)
{
	// create a one volume volume set
	RicVolumeSet Vset(nx,ny,nz,1);
	
	// copy values to volume set
	Vset.dx = Vset.VolSet[0].dx = dx;
	Vset.dy = Vset.VolSet[0].dy = dy;
	Vset.dz = Vset.VolSet[0].dz = dz;
	Vset.xoffset = Vset.VolSet[0].xoffset = xoffset;
	Vset.yoffset = Vset.VolSet[0].yoffset = yoffset;
	Vset.zoffset = Vset.VolSet[0].zoffset = zoffset;
	Vset.VolSet[0].min = min;
	Vset.VolSet[0].max = max;
	
	
	// copy to the first volume
	int i,j,k;
	for ( i=0; i<nx; i++)
		for ( j=0; j<ny; j++)   
			for ( k=0; k<nz; k++)
				Vset.VolSet[0].vox[i][j][k] = vox[i][j][k];   
	
	return Vset.Write_NEMA(filename,orient_str);
	
}
