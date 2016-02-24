// 3D Matrix class from 3dKindoms
// Modified slightly for RIC use
// All angles in radians

#ifndef CMatrix_h
#define CMatrix_h

#include "RicPoint.h"
#include <math.h>

class CMatrix 
{
public:
// Data
	float mf[ 16 ];

// Functions
	CMatrix( const int bIdentity = true )
	{
		if ( bIdentity ) Identity();
	}

	void Identity( )
	{
		mf[ 0] = 1.0f;    mf[ 1] = 0.0f;      mf[ 2] = 0.0f;    mf[ 3] = 0.0f;  
		mf[ 4] = 0.0f;    mf[ 5] = 1.0f;      mf[ 6] = 0.0f;    mf[ 7] = 0.0f;  
		mf[ 8] = 0.0f;    mf[ 9] = 0.0f;      mf[10] = 1.0f;    mf[11] = 0.0f;  
		mf[12] = 0.0f;    mf[13] = 0.0f;      mf[14] = 0.0f;    mf[15] = 1.0f;
	}

	// Concatenate 2 matrices with the * operator
	inline CMatrix operator* (const CMatrix &InM) const
	{
		CMatrix Result( 0 );
		for (int i=0;i<16;i+=4)
		{
			for (int j=0;j<4;j++)
			{
			Result.mf[i + j] = mf[ i + 0] * InM.mf[ 0 + j] + mf[ i + 1] * InM.mf[ 4 + j]
				+ mf[ i + 2] * InM.mf[ 8 + j] + mf[ i + 3] * InM.mf[ 12 + j];
			}
		}
		return Result;
	}

	// Use a matrix to transform a 3D point with the * operator
	inline Point operator* (const Point &Pnt ) const
	{
		float x = Pnt.x*mf[0] + Pnt.y*mf[4] + Pnt.z*mf[8]  + mf[12];
		float y = Pnt.x*mf[1] + Pnt.y*mf[5] + Pnt.z*mf[9]  + mf[13];
		float z = Pnt.x*mf[2] + Pnt.y*mf[6] + Pnt.z*mf[10] + mf[14]; 
		return Point(x,y,z);
	}

	// Rotate the *this matrix fRad counter-clockwise around a single axis( either x, y, or z )
	void Rotate( float fRad, int x, int y, int z )
	{
		CMatrix Temp;
		if (x == 1) Temp.RotX( -fRad );
		if (y == 1) Temp.RotY( -fRad );
		if (z == 1) Temp.RotZ( -fRad );
		*this = Temp * (*this);
	}

	void Scale( float sx, float sy, float sz )
	{
		int x;
		for (x = 0; x <  4; x++) mf[x]*=sx;
		for (x = 4; x <  8; x++) mf[x]*=sy;
		for (x = 8; x < 12; x++) mf[x]*=sz;
	}

	void Translate( const Point &Test )
	{
		for (int j=0;j<4;j++)
		{
			mf[12+j] += Test.x * mf[j] + Test.y * mf[4+j] + Test.z * mf[8+j]; 
		}	 
	}
	
	Point GetTranslate( )
	{
		return Point( mf[12], mf[13], mf[14] );
	}
	
	// Zero out the translation part of the matrix
	CMatrix RotationOnly( )
	{
		CMatrix Temp = *this;
		Temp.mf[12] = 0;
		Temp.mf[13] = 0;
		Temp.mf[14] = 0;
		return Temp;
	}
	
	// transpose the rows and columns
	CMatrix Transpose()
	{
			CMatrix Temp;
			Temp.mf[0] = mf[0];
			Temp.mf[1] = mf[4];
			Temp.mf[2] = mf[8];
			Temp.mf[3] = mf[12];
			Temp.mf[4] = mf[1];
			Temp.mf[5] = mf[5];
			Temp.mf[6] = mf[9];
			Temp.mf[7] = mf[13];
			Temp.mf[8] = mf[2];
			Temp.mf[9] = mf[6];
			Temp.mf[10] = mf[10];
			Temp.mf[11] = mf[14];
			Temp.mf[12] = mf[3];
			Temp.mf[13] = mf[7];
			Temp.mf[14] = mf[11];
			Temp.mf[15] = mf[15];
			return Temp;
	}
	
	// Create a rotation matrix for a counter-clockwise rotation of fRad around an arbitrary axis(x, y, z)
	void RotateMatrix( float fRad, float x, float y, float z)
	{
		Identity();
		float cosA = cosf(fRad);
		float sinA = sinf(fRad);
		float m = 1.0f - cosA;
		mf[0] = cosA + x*x*m;
		mf[5] = cosA + y*y*m;
		mf[10]= cosA + z*z*m;
		
		float tmp1 = x*y*m;
		
		float tmp2 = z*sinA;
		mf[4] = tmp1 + tmp2;
		mf[1] = tmp1 - tmp2;
		
		tmp1 = x*z*m;
		tmp2 = y*sinA;
		mf[8] = tmp1 - tmp2;
		mf[2] = tmp1 + tmp2;
		
		tmp1 = y*z*m;
		tmp2 = x*sinA;
		mf[9] = tmp1 + tmp2;
		mf[6] = tmp1 - tmp2;
	}

	// Simple but not robust matrix inversion. (Doesn't work properly if there is a scaling or skewing transformation.)
	inline CMatrix InvertSimple()
	{
		CMatrix R(0);
		R.mf[0]  = mf[0]; 		R.mf[1]  = mf[4];		R.mf[2]  = mf[8];	R.mf[3]  = 0.0f;
		R.mf[4]  = mf[1];		R.mf[5]  = mf[5];		R.mf[6]  = mf[9];	R.mf[7]  = 0.0f;
		R.mf[8]  = mf[2];		R.mf[9]  = mf[6];		R.mf[10] = mf[10];	R.mf[11] = 0.0f;
		R.mf[12] = -(mf[12]*mf[0]) - (mf[13]*mf[1]) - (mf[14]*mf[2]);
		R.mf[13] = -(mf[12]*mf[4]) - (mf[13]*mf[5]) - (mf[14]*mf[6]);
		R.mf[14] = -(mf[12]*mf[8]) - (mf[13]*mf[9]) - (mf[14]*mf[10]);
		R.mf[15] = 1.0f;
		return R;
	}
	
	// Invert for only a rotation, any translation is zeroed out
	CMatrix InvertRot( )
	{
		CMatrix R( 0 );
		R.mf[0]  = mf[0]; 		R.mf[1]  = mf[4];		R.mf[2]  = mf[8];	R.mf[3]  = 0.0f;
		R.mf[4]  = mf[1];		R.mf[5]  = mf[5];		R.mf[6]  = mf[9];	R.mf[7]  = 0.0f;
		R.mf[8]  = mf[2];		R.mf[9]  = mf[6];		R.mf[10] = mf[10];	R.mf[11] = 0.0f;
		R.mf[12] = 0;			R.mf[13] = 0;			R.mf[14] = 0;		R.mf[15] = 1.0f;
		return R;
	}	


private:
	// helpers for Rotate
	void RotX(float angle)
	{  
		mf[5]  = cosf(angle);
		mf[6]  = sinf(angle);
		mf[9]  = -sinf(angle);
		mf[10] = cosf(angle);
	}
	void RotY(float angle)
	{
		mf[0]  =  cosf(angle);
		mf[2]  =  -sinf(angle);
		mf[8]  =  sinf(angle);
		mf[10] =  cosf(angle);    
	}
	void RotZ(float angle)
	{
		mf[0] =  cosf(angle);
		mf[1] =  sinf(angle);
		mf[4] =  -sinf(angle);
		mf[5] =  cosf(angle);
	}
};

#endif

