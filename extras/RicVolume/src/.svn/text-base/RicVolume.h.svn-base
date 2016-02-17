// ------------------------ RicVolume.h --------------------------------
/*! \file
Header file for the RicVolume class. A set of member functions allows
manipulation of voxel data.
*/


#ifndef _RICVOLUME_H
#define _RICVOLUME_H

using namespace std;

#include <cstdlib>
#include <iostream>
#include <string>
#include "RicUtil.h"

class RicVolumeSet;

/*!
The voxel data is stored in a three dimensional array where the origin, index [0][0][0],
is located at the left, anterior, superior. In other words, the array start at the
top, left, front of the head. This corresponds to a NEMA orientation string of 
XYZ---.

Each voxel consists of a floating point intensity value. The are overloads of
common arithmatic operators so whole volumes can be operated on.

Copyright (C) 2007 by Bill Rogers - Research Imaging Center - UTHSCSA
rogers@uthscsa.edu     
*/
class RicVolume
{
	friend class RicVolumeSet;
	
	public:
	float 	***vox;		///< voxel matrix
	float	*varray;	///< pointer to first element of voxel matrix
	float	min;		///< minimum voxel value
	float	max;		///< maximum voxel value
	float	avg;		///< average voxel value
	float	dx;			///< spacing along x axis
	float	dy;			///< spacing along y axis
	float	dz;			///< spacing along z axis
	float	xoffset;	///< Talairach x offset
	float	yoffset;	///< Talairach y offset
	float	zoffset;	///< Talairach z offset
	int		nvox;		///< total number of voxels
	
	int		nx;			///< number of z values in volume (columns per row)
	int		ny;			///< number of y values in volume (rows per slice)
	int		nz;			///< number of x values in volume (slices)

	public:
	// constructors
	RicVolume();
	RicVolume(string filename);
	RicVolume(int nx, int ny, int nz);
	RicVolume(int nx, int ny, int nz, float val);
	~RicVolume();
	
	// read and write routines
	int Read(string filename);
	int Write(string filename);
	int Write_NEMA(string filename);
	int Write_NEMA(string filename,string orient_str);
	
	// operator overloads functions
	RicVolume &operator + (RicVolume &rhs);
	RicVolume &operator - (RicVolume &rhs);
	RicVolume &operator * (RicVolume &rhs);
	RicVolume & operator / (RicVolume &rhs);

	RicVolume& operator+=(const float con);
	RicVolume& operator-=(const float con);
	RicVolume& operator*=(const float con);
	RicVolume& operator/=(const float con);
	
	RicVolume *operator > (const float con);
	
	// member functions
	int Init(int nx, int ny, int nz, float val);
	int Init(int nx, int ny, int nz);
	
	/// returns number of x values
	int get_numx(){return nx;};
	
	/// returns number of y values
	int get_numy(){return ny;};
	
	///< returns number of z values
	int get_numz(){return nz;};	
	
	// for PeterK
	/// returns number of x values
	int num_cols(){return nx;};	
	
	/// returns number of y values
	int num_rows(){return ny;};	
	
	/// returns number of z values
	int num_slices(){return nz;};
	
	void CalcMinMaxAvg(void);
	void Fill(float fillval);
	void ThresholdAbove(float threshold, float fillval);
	void ThresholdBelow(float threshold, float fillval);
	void Copy(RicVolume &source);

	int get_6_neighbors(int i, int j, int k, float *buff6);
	int get_18_neighbors(int ix, int iy, int iz, float *buff18);
	int get_26_neighbors(int ix, int iy, int iz, float *buff26);

	float interp3D(float x, float y, float z);
	int set_at_nearest(float rx, float ry, float rz, float value);

	 //Where i, j, and k are column, row, and slice indices. 
	// obsolete functions put here for Peter's benefit
	float get_at_index(int i, int j, int k, int flag = 0);
	float range_checked_get_at_index(int i, int j, int k);
	int set_at_index(int i, int j, int k, float value); 
	
	private:
};

#endif // _RIC_VOLUME_H
