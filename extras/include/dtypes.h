/////////////////////////////////////////////////////////////////////
/// University of Texas at San Antonio - Research Imaging Center  ///
/////////////////////////////////////////////////////////////////////
//  Module Name: dtypes.H
//  Birth Date: 
//  Programmer: Dan Nickerson / Hunter Downs
/////////////////////////////////////////////////////////////////////
/////////////////////////  Description  /////////////////////////////
// Some convenience typedefs used by VolClass and various
// other convex hull apps.
/////////////////////////////////////////////////////////////////////
//////////////////////////  Revisions  //////////////////////////////
//$Log: dtypes.H,v $
//Revision 1.3  1998/07/20 20:22:03  nickersd
//Moved rcsid into the #ifdef block.
//
//Revision 1.2  1998/07/20 19:05:27  nickersd
//Changed rcsid to static.
//
//Revision 1.1  1998/07/15 21:21:37  nickersd
//Initial revision
//
/////////////////////////////////////////////////////////////////////

#ifndef DTYPES
#define DTYPES
//static const char *dtypes_H_rcsid = "$Id: dtypes.H,v 1.3 1998/07/20 20:22:03 nickersd Exp $";

typedef unsigned short UNS ;
typedef unsigned long  UNL ;
typedef unsigned char  UNC ;
typedef double  XYPOINT[2];    /*  2D point (on a contour) */
typedef double  XYZPOINT[3];    /*  3D point in rectangular coordinates  */
typedef double RPTPOINT[3] ;   /*  3D point in spherical coordinates */
typedef double HOMOMATRIX[4][4];
typedef double HOMOPOINT[4];
typedef double MATRIX3x3[3][3];


#endif
