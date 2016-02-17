// --------------------------- NiftiOrientConvert.cpp -------------------------
/*! \file

Spike created a java NEMA orientation converter from the original nifti c++
source code. This version is derived from Spike's java code.

Copyright (C) 2007  UTHSCSA - Research Imaging Center

rogers@uthscsa.edu
*/


#include <cstdlib>
#include <math.h>
#include <nifti1_io.h>
#include <string>
using namespace std;

/*!
compute the (closest) orientation from a 4x4 ijk->xyz tranformation matrix

Input:  4x4 matrix that transforms (i,j,k) indexes to (x,y,z) coordinates,
	   where +x=Right, +y=Anterior, +z=Superior.
	   (Only the upper-left 3x3 corner of R is used herein.)
Output: 3 orientation codes that correspond to the closest "standard"
	   anatomical orientation of the (i,j,k) axes.
Method: Find which permutation of (x,y,z) has the smallest angle to the
	   (i,j,k) axes directions, which are the columns of the R matrix.
Errors: The codes returned will be zero.

For example, an axial volume might get return values of
 *icod = NIFTI_R2L   (i axis is mostly Right to Left)
 *jcod = NIFTI_P2A   (j axis is mostly Posterior to Anterior)
 *kcod = NIFTI_I2S   (k axis is mostly Inferior to Superior)


see "QUATERNION REPRESENTATION OF ROTATION MATRIX" in nifti1.h

see nifti_quatern_to_mat44, nifti_mat44_to_quatern,	nifti_make_orthog_mat44

Taken from nifti1_io.c function nifti_mat44_to_orientation()
converted to return a NEMA-style orientation string

@param R - nifti transformation matrix
@return - NEMA-style orientation string
*/
string convertNiftiSFormToNEMA(mat44 R)
{
	double xi, xj, xk, yi, yj, yk, zi, zj, zk, val, detQ, detP;
	mat33 P, Q, M;
	int i, j, k, p, q, r, ibest, jbest, kbest, pbest, qbest, rbest;
	k = 0;
	double vbest;
	string OString="bogus";

//	Q = new double[3][3];
//	P = new double[3][3];

	//if( icod == NULL || jcod == NULL || kcod == NULL ) return ; /* bad */

	//*icod = *jcod = *kcod = 0 ; /* error returns, if sh*t happens */

	/* load column vectors for each (i,j,k) direction from matrix */

	/*-- i axis --*/ /*-- j axis --*/ /*-- k axis --*/

	xi = R.m[0][0]; xj = R.m[0][1]; xk = R.m[0][2];
	yi = R.m[1][0]; yj = R.m[1][1]; yk = R.m[1][2];
	zi = R.m[2][0]; zj = R.m[2][1]; zk = R.m[2][2];

	/* normalize column vectors to get unit vectors along each ijk-axis */

	/* normalize i axis */
	val = sqrt( xi*xi + yi*yi + zi*zi ) ;
	if( val == 0.0 ) return NULL;                 /* stupid input */
	xi /= val; yi /= val; zi /= val;

	/* normalize j axis */
	val = sqrt( xj*xj + yj*yj + zj*zj ) ;
	if( val == 0.0 ) return NULL;                 /* stupid input */
	xj /= val; yj /= val; zj /= val;

	/* orthogonalize j axis to i axis, if needed */
	val = xi*xj + yi*yj + zi*zj ;    /* dot product between i and j */
	if( fabs(val) > 1.E-4 ) 
	{
		
		xj -= val*xi ; yj -= val*yi ; zj -= val*zi ;
		val = sqrt( xj*xj + yj*yj + zj*zj ) ;  /* must renormalize */
		if( val == 0.0 ) return OString;              /* j was parallel to i? */
		xj /= val; yj /= val; zj /= val;
	}

	/* normalize k axis; if it is zero, make it the cross product i x j */
	val = sqrt( xk*xk + yk*yk + zk*zk ) ;
	if( val == 0.0 ) 
	{ 
		xk = yi*zj-zi*yj; 
		yk = zi*xj-zj*xi; 
		zk=xi*yj-yi*xj; 
	}
	else
	{ 
		xk /= val; 
		yk /= val; 
		zk /= val; 
	}

	/* orthogonalize k to i */
	val = xi*xk + yi*yk + zi*zk;    /* dot product between i and k */
	if( fabs(val) > 1.E-4 )
	{
		xk -= val*xi; yk -= val*yi; zk -= val*zi;
		val = sqrt( xk*xk + yk*yk + zk*zk );
		if( val == 0.0 ) return OString;      /* bad */
		xk /= val; yk /= val; zk /= val;
	}

	/* orthogonalize k to j */
	val = xj*xk + yj*yk + zj*zk;    /* dot product between j and k */
	if( fabs(val) > 1.e-4 )
	{
		
		xk -= val*xj ; yk -= val*yj ; zk -= val*zj ;
		val = sqrt( xk*xk + yk*yk + zk*zk ) ;
		if( val == 0.0 ) return OString;      /* bad */
		xk /= val ; yk /= val ; zk /= val ;
	}

	Q.m[0][0] = xi; Q.m[0][1] = xj; Q.m[0][2] = xk;
	Q.m[1][0] = yi; Q.m[1][1] = yj; Q.m[1][2] = yk;
	Q.m[2][0] = zi; Q.m[2][1] = zj; Q.m[2][2] = zk;

	/* at this point, Q is the rotation matrix from the (i,j,k) to (x,y,z) axes */

	detQ = nifti_mat33_determ(Q) ;
	if( detQ == 0.0 ) return OString; /* shouldn't happen unless user is a DUFIS */

	/* Build and test all possible +1/-1 coordinate permutation matrices P;
		then find the P such that the rotation matrix M=PQ is closest to the
		identity, in the sense of M having the smallest total rotation angle. */

	/* Despite the formidable looking 6 nested loops, there are
		only 3*3*3*2*2*2 = 216 passes, which will run very quickly. */

	vbest = -666.0 ; ibest=pbest=qbest=rbest=1 ; jbest=2 ; kbest=3 ;
	for( i=1 ; i <= 3 ; i++ )
	{     /* i = column number to use for row #1 */
		for( j=1 ; j <= 3 ; j++ )
		{    /* j = column number to use for row #2 */
			if( i == j ) continue ;
			for( k=1 ; k <= 3 ; k++ )
			{  /* k = column number to use for row #3 */
				if( i == k || j == k ) continue ;
				P.m[0][0] = P.m[0][1] = P.m[0][2] = \
				P.m[1][0] = P.m[1][1] = P.m[1][2] = \
				P.m[2][0] = P.m[2][1] = P.m[2][2] = 0.0 ;
				for( p=-1 ; p <= 1 ; p+=2 )
				{    /* p,q,r are -1 or +1      */
					for( q=-1 ; q <= 1 ; q+=2 )
					{   /* and go into rows #1,2,3 */
						for( r=-1 ; r <= 1 ; r+=2 )
						{
							P.m[0][i-1] = p ; P.m[1][j-1] = q ; P.m[2][k-1] = r ;
							detP = nifti_mat33_determ(P) ;           /* sign of permutation */
							if( (detP * detQ) <= 0.0 ) continue ;  /* doesn't match sign of Q */
							M = nifti_mat33_mul(P,Q) ;
				
							/* angle of M rotation = 2.0*acos(0.5*sqrt(1.0+trace(M)))       */
							/* we want largest trace(M) == smallest angle == M nearest to I */
				
							val = M.m[0][0] + M.m[1][1] + M.m[2][2] ; /* trace */
							if( val > vbest ) 
							{
								vbest = val ;
								ibest = i ; jbest = j ; kbest = k ;
								pbest = p ; qbest = q ; rbest = r ;
							}
						}
					}
				}
			}
		}
	}

	/* At this point ibest is 1 or 2 or 3; pbest is -1 or +1; etc.

		The matrix P that corresponds is the best permutation approximation
		to Q-inverse; that is, P (approximately) takes (x,y,z) coordinates
		to the (i,j,k) axes.

		For example, the first row of P (which contains pbest in column ibest)
		determines the way the i axis points relative to the anatomical
		(x,y,z) axes.  If ibest is 2, then the i axis is along the y axis,
		which is direction P2A (if pbest > 0) or A2P (if pbest < 0).

		So, using ibest and pbest, we can assign the output code for
		the i axis.  Mutatis mutandis for the j and k axes, of course. */

	string iChar, jChar, kChar, iSense, jSense, kSense;
//	iChar = jChar = kChar = iSense = jSense = kSense = 0;

	switch( ibest*pbest ) 
	{
	
			case  1: /*i = NIFTI_L2R*/ iChar = "X"; iSense = "+"; break;
			case -1: /*i = NIFTI_R2L*/ iChar = "X"; iSense = "-"; break;
			case  2: /*i = NIFTI_P2A*/ iChar = "Y"; iSense = "+"; break;
			case -2: /*i = NIFTI_A2P*/ iChar = "Y"; iSense = "-"; break;
			case  3: /*i = NIFTI_I2S*/ iChar = "Z"; iSense = "+"; break;
			case -3: /*i = NIFTI_S2I*/ iChar = "Z"; iSense = "-"; break;
	}

	switch( jbest*qbest ) 
	{
	
			case  1: /*j = NIFTI_L2R*/ jChar = "X"; jSense = "+"; break;
			case -1: /*j = NIFTI_R2L*/ jChar = "X"; jSense = "-"; break;
			case  2: /*j = NIFTI_P2A*/ jChar = "Y"; jSense = "+"; break;
			case -2: /*j = NIFTI_A2P*/ jChar = "Y"; jSense = "-"; break;
			case  3: /*j = NIFTI_I2S*/ jChar = "Z"; jSense = "+"; break;
			case -3: /*j = NIFTI_S2I*/ jChar = "Z"; jSense = "-"; break;
	}

	switch( kbest*rbest ) 
	{
	
			case  1: /*k = NIFTI_L2R*/ kChar = "X"; kSense = "+"; break;
			case -1: /*k = NIFTI_R2L*/ kChar = "X"; kSense = "-"; break;
			case  2: /*k = NIFTI_P2A*/ kChar = "Y"; kSense = "+"; break;
			case -2: /*k = NIFTI_A2P*/ kChar = "Y"; kSense = "-"; break;
			case  3: /*k = NIFTI_I2S*/ kChar = "Z"; kSense = "+"; break;
			case -3: /*k = NIFTI_S2I*/ kChar = "Z"; kSense = "-"; break;
	}
	
	OString = iChar +jChar + kChar + iSense + jSense + kSense;
	return OString;
}

/*!
@brief
Given the quaternion parameters (etc.), compute a transformation matrix.

Given the quaternion parameters (etc.), compute a transformation matrix.

See comments in nifti1.h for details.
 - qb,qc,qd = quaternion parameters
 - qx,qy,qz = offset parameters
 - dx,dy,dz = grid stepsizes (non-negative inputs are set to 1.0)
 - qfac     = sign of dz step (< 0 is negative; >= 0 is positive)


If qx=qy=qz=0, dx=dy=dz=1, then the output is a rotation matrix.
For qfac >= 0, the rotation is proper.
For qfac <  0, the rotation is improper.

see "QUATERNION REPRESENTATION OF ROTATION MATRIX" in nifti1.h
see nifti_mat44_to_quatern, nifti_make_orthog_mat44, nifti_mat44_to_orientation

Taken from nifti1_io.c function nifti_quatern_to_mat44 (converted to Java)

@param qb - quaternion parameter
@param qc - quaternion parameter
@param qd - quaternion parameter
@param qx - offset parameter
@param qy - offset parameter
@param qz - offset parameter
@param dx - grid stepsize (non-negative inputs are set to 1.0)
@param dy - grid stepsize (non-negative inputs are set to 1.0)
@param dz - grid stepsize (non-negative inputs are set to 1.0)
@param qfac - sign of dz step (< 0 is negative; >= 0 is positive)
@return - transformation matrix
*/
mat44 convertNiftiQFormToNiftiSForm(double qb, double qc, double qd, double qx, double qy, double qz, double dx, 
		double dy, double dz, double qfac) 
{

	mat44 R;
	double a, b, c, d, xd, yd, zd;
	b = qb;
	c = qc;
	d = qd;

	// last row is always [ 0 0 0 1 ]
	R.m[3][0]=R.m[3][1]=R.m[3][2] = 0.0; R.m[3][3]= 1.0;

	// compute a parameter from b,c,d
	a = 1.0 - (b*b + c*c + d*d) ;
	if( a < 1.E-7 ) 
	{                   /* special case */
		
		a = 1.0 / sqrt(b*b+c*c+d*d) ;
		b *= a ; c *= a ; d *= a ;        /* normalize (b,c,d) vector */
		a = 0.0 ;                        /* a = 0 ==> 180 degree rotation */
	} 
	else 
	{
		
		a = sqrt(a) ;                     /* angle = 2*arccos(a) */
	}

	// load rotation matrix, including scaling factors for voxel sizes
	xd = (dx > 0.0) ? dx : 1.0 ;       /* make sure are positive */
	yd = (dy > 0.0) ? dy : 1.0 ;
	zd = (dz > 0.0) ? dz : 1.0 ;

	if( qfac < 0.0 ) zd = -zd ;         /* left handedness? */

	R.m[0][0] =       (a*a+b*b-c*c-d*d) * xd ;
	R.m[0][1] = 2.0 * (b*c-a*d        ) * yd ;
	R.m[0][2] = 2.0 * (b*d+a*c        ) * zd ;
	R.m[1][0] = 2.0 * (b*c+a*d        ) * xd ;
	R.m[1][1] =       (a*a+c*c-b*b-d*d) * yd ;
	R.m[1][2] = 2.0 * (c*d-a*b        ) * zd ;
	R.m[2][0] = 2.0 * (b*d-a*c        ) * xd ;
	R.m[2][1] = 2.0 * (c*d+a*b        ) * yd ;
	R.m[2][2] =       (a*a+d*d-c*c-b*b) * zd ;

	// load offsets
	R.m[0][3] = qx ; R.m[1][3] = qy ; R.m[2][3] = qz ;

	return R ;
}


	// multiply 2 3x3 matrices
	// Taken from nifti1_io.c (converted to Java)
/*mat33 nifti_mat33_mul( mat33 A , mat33 B )  // multiply 2 3x3 matrices
{
   mat33 C ; int i,j ;
   for( i=0 ; i < 3 ; i++ )
    for( j=0 ; j < 3 ; j++ )
      C.m[i][j] =  A.m[i][0] * B.m[0][j]
                 + A.m[i][1] * B.m[1][j]
                 + A.m[i][2] * B.m[2][j] ;
   return C ;
}*/
