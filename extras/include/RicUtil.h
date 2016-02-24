/* ---------------------------- RicUtil.h ---------------------------- */
/*! \file
This is the header file for the math utilities. There are matrix allocation
routines, geometry routines, and coordinate conversion routines.

This as been extended for use with RIC projects
Bill Rogers
 */


#if !defined (_RICUTIL_H)
#define _RICUTIL_H 1

#include "RicPoint.h"
#include "RicMatrix.h"
#include <stdio.h>

#define MAX_STR			256 ///< max string length
//#define MSIZE 32
#define PI              3.14159265358979323846f
#define DTOR			PI/180.0f		///< degrees to radians
#define RTOD			180.0f/PI		///< radians to degrees
#define MM_PER_IN 		25.4f			///< 25.4 mm per inch
#define IN_PER_MM 		1.0f/MM_PER_IN	///< inches per mm
#define IN_PER_CM 		10.0f*IN_PER_MM	///< inches per cm
#define ERRVAL (float)9999 ///< table value for erroneous value

void error_msg(char *msg);


// typedefs

/// Point, IPoint, and DPoint are now classes in RicPoint.h
/// 3d point - float
//typedef struct { float x; float y; float z; }  Point;
/// 3d point - integer
///typedef struct { int x; int y; int z; }  IPnt;
/// 3d point - double
// typedef struct { double x; double y; double z;} DPoint;

/// 3d cylindrical point - float
typedef struct { float r; float theta; float z; } Cylind;
/// 3d cylindrical point - integer
typedef struct { int r; int theta; int z; } ICyl;
/// 3d cylindrical point - double
typedef struct { double r; double theta; double z; } DCylind;

/// 3d point in spherical coordinates - float
typedef struct 	{	float r; float theta; float phi; } Sphere;

/// 3d point in spherical coordinates - double
typedef struct 	{	double r; double theta; double phi; } DSphere;

/// 3d vector
typedef Point Vector;

/// complex number
typedef struct { float re, im;} Complex;

/// vertex on surface
typedef struct { int slice,pro; } vertx;
//typedef struct { int slice,pro; } Vertex;

/// triangle structure
typedef struct { Point p[3]; } Triangle;

/// color structure
typedef struct { float r; float g; float b; } Color;

/// class for rotating points about an axis
class Rotate
{
private:

	float angle;	///< rotation angle
	// point to rotate about
	float px,py,pz;	///< point to rotate about
	// vector for rotation axis
	float vx,vy,vz;	///< vector for rotation axis
	// variables for transformation matrix
	float r11, r12, r13, r21, r22, r23,r31, r32, r33, r41, r42, r43;

public:
	void InitRotate(float px, float py, float pz, float vx, float vy,
		float vz, float angle);
	void RotateIt(float *x, float *y, float *z);
	Point RotateIt(Point pnt);

};

/// Template version of 2D matrix memory allocation routines
template <class T>
int matrix(T*** mat, int nrow, int ncol)
{
	// allocate pointers for the start of each row
	T** m;
	m = new T *[nrow];
	if (!m)
	{
		char mes[MAX_STR];
		sprintf (mes, "couldn't allocate %d elements in matrix(1)",
		nrow * (int)sizeof (float *));
		error_msg (mes);
		return 0;
	}

	// allocate linear array for all points
	T* tmp;
	tmp = new T[nrow*ncol];
	if (!tmp) {
		char mes[MAX_STR];
		sprintf (mes, "couldn't allocate %d elements in matrix(2)",
		nrow * ncol * (int)sizeof (float));
		error_msg (mes);
		return 0;
	}
	for( int i=0 ; i<nrow ; i++) m[i]= &tmp[i*ncol];

	*mat = m;

	return 1;
}

/// Template version of 2D matrix memory free routines
template <class T>
void free_matrix(T **m)
{
	// first free linear array
	delete m[0];

	// free the array of pointers
	delete m;
}

/// Template version of 3D matrix memory allocation routines
template <class T>
		int matrix3D(T**** mat, int nslice, int nrow, int ncol)
{
	int i,j;
	char mes[MAX_STR];

	// allocate pointers for the start of each slice
	T*** m;
	m = new T **[nslice];
	if (!m)
	{
		sprintf (mes, "couldn't allocate %d elements in matrix3D(0)",
				 nslice);
		error_msg (mes);
		return 0;
	}

	// allocate pointers for the start of each row
	T** rp;
	rp = new T *[nslice*nrow];
	if (!rp)
	{
		sprintf (mes, "couldn't allocate %d elements in matrix3D(1)",
				 nrow);
		error_msg (mes);
		return 0;
	}
	for ( i=0 ; i<nslice ; ++i ) m[i] = &rp[i*nrow];

	// allocate linear array for all points
	T* rc;
	rc = new T[nslice*nrow*ncol];
	if (!rc) {
		sprintf (mes, "couldn't allocate %d elements in matrix3D(2)",
				 ncol);
		error_msg (mes);
		return 0;
	}
	for( i=0 ; i<nslice ; i++)
		for ( j=0 ; j<nrow ; ++j )
			m[i][j] = &rc[i*nrow*ncol + j*ncol];

	*mat = m;

	return 1;
}

/// Template version of 2D matrix memory free routines
template <class T>
void free_matrix3D(T **m)
{
	// first free linear array
	delete m[0][0];

	// free the array of pointers
	delete m[0];
	delete m;
}

// functions implemented as defines
#define dist_squ(a,b)	(Dist_squ(&(a),&(b)))
#define point_theta(p)	((float)atan2((p).y, (p).x))

// function defines
Cylind cart_to_cyl(Point cart);
Point cyl_to_cart(Cylind cyl);
Cylind cart_to_cyl_deg(Point cart);
Point cyl_to_cart_deg(Cylind cyl);
Sphere  cart_to_sphere(Point cart);
DSphere  cart_to_sphere(DPoint cart);
Point sphere_to_cart(Sphere sphere);
DPoint sphere_to_cart(DSphere sphere);
float dist(Point p1,Point p2);
double ddist(DPoint p1,DPoint p2);
float distsqu(Point p1,Point p2);
float cyldist(Cylind p1,Cylind p2);
float cyldist2(Cylind p1,Cylind p2,float *dval);
float cyldistsqu(Cylind p1,Cylind p2);
void int_sort(int a[],int n);
void float_sort(float a[],int n);
void theta_sort(Cylind a[],int n);
void y_sort(Point a[],int n);
void z_sort(Cylind a[],int n);
void z_sort(Point a[],int n);
Vector	cross_product(Point pnt[]);
DPoint	cross_product(DPoint pnt[]);
float dot_product(Vector *v1, Vector *v2);
Vector normal_pnt (Point p1, Point p2,Point p3);
DPoint normal_pnt (DPoint p1, DPoint p2,DPoint p3);
int file_exist(char *filename);
Point* catrom_curve (Point *pin,int nctrl, int step_size, float tension,
	int *noutpnts);
Point line_thru_plane (float a, float b, float c, float d, Point p0, Point p1);
int line_thru_triangle(Point t0,Point t1,Point t2,Point p0,Point p1,Point *pout);
int line_intersect_triangle(Point t0,Point t1,Point t2,Point p0,Point p1,
		Point *pout);
int line_intersect_plane(float a,float b,float c,float d,Point p0,
	Point p1,Point *pout);
void equ_plane(Point p1,Point p2,Point p3,float *a,float *b,float *c,float *d);
void equ_plane2(Point p1,Vector n1,float *a,float *b,float *c,float *d);
float Dist_squ (Point *p1, Point *p2);
void initialize_in_polygon (Point *v, int npnts);
int in_polygon (Point *v, int npnts, Point p);
int in_polygon2 (Point *v, int n, Point p, int add_pi);
Point interpolate_along_z(Point *p, float z);
float circle_from_3pnts(Point p1, Point p2, Point p3, Point *cp);
int scale_array(float newMin, float newMax, float *array, int l, int invert);

void splint(float xa[],float ya[],float y2a[],int n,float x,float *y);
void spline(float x[],float y[],int n,float yp1,float ypn,float y2[]);
void spline_d(double *, double *, int, double, double, double *);
void splint_d(double *, double *, double *, int, double, double *);
float angle_between_planes(float a1,float b1,float c1,float a2, float b2,
						   float c2);
float dist_point_plane (float a,float b,float c,float d,Point p0, Point *pout);
int dist_point_line(Point Pnt, Point LineStart, Point LineEnd, float *Distance, Point *Pout);
double distpntlinespace(Point,double,double,double,double,double,double);
float area_triangle (float a, float b, float c);
float area_triangle (Point p1, Point p2, Point p3);
int inside_triangle (Point t1, Point t2, Point t3, Point p);
int int_circle_line (Point c, float r, Point p1, Point p2);
int inbetween(Point lp1, Point lp2, Point p);
Point line_intersect (Point p11, Point p12, Point p21, Point p22);
int sphere_line_intersection (Point center, float r, Point l1, Point l2,
							  Point *p1, Point *p2);

// NR in C routines
void splint(float xa[],float ya[],float y2a[],int n,float x,float *y);
void spline(float x[],float y[],int n,float yp1,float ypn,float y2[]);
void splie2(float x1a[],float x2a[],float **ya,int m,int n,float **y2a);
void splin2(float x1a[],float x2a[],float **ya,float **y2a,int m,int n,
	float x1,float x2,float *y);


// ---------------------- function error_msg -----------------------
// Darn near everybody calls this one
#ifdef WINDOWS

/// Error message to generic Windows popup
inline void error_msg(char *msg)
{
	AfxMessageBox(msg,MB_OK,NULL);
}
#else

/// Error message to stderr
inline void error_msg(char *msg)
{
	fprintf(stderr,"%s",msg);
}

#endif

#endif /* _RICUTIL_H */

