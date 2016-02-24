// 3D Vector class from 3dKindoms - thanks
// Just slightly modified for RIC use
// Point, IPoint and DPoint classes are defined

#ifndef _RicPoint_h
#define _RicPoint_h

#include <math.h>

/// XYZ floating point class
class Point
{
public:
	// Data
	float x, y, z;

	// Ctors
	Point( float InX, float InY, float InZ ) : x( InX ), y( InY ), z( InZ )
	{
	}
	Point( ) : x(0), y(0), z(0)
	{
	}

	// Operator Overloads
	inline bool operator== (const Point& V2) const 
	{
		return (x == V2.x && y == V2.y && z == V2.z);
	}

	inline Point operator+ (const Point& V2) const 
	{
		return Point( x + V2.x,  y + V2.y,  z + V2.z);
	}
	
	inline Point operator- (const Point& V2) const
	{
		return Point( x - V2.x,  y - V2.y,  z - V2.z);
	}
	
	inline Point operator- ( ) const
	{
		return Point(-x, -y, -z);
	}

	inline Point operator/ (float S ) const
	{
		float fInv = 1.0f / S;
		return Point (x * fInv , y * fInv, z * fInv);
	}
	
	inline Point operator/ (const Point& V2) const
	{
		return Point (x / V2.x,  y / V2.y,  z / V2.z);
	}
	
	inline Point operator* (const Point& V2) const
	{
		return Point (x * V2.x,  y * V2.y,  z * V2.z);
	}
	
	inline Point operator* (float S) const
	{
		return Point (x * S,  y * S,  z * S);
	}

	inline void operator+= ( const Point& V2 )
	{
		x += V2.x;
		y += V2.y;
		z += V2.z;
	}
	
	inline void operator-= ( const Point& V2 )
	{
		x -= V2.x;
		y -= V2.y;
		z -= V2.z;
	}

	inline float operator[] ( int i )
	{
		if ( i == 0 ) return x;
		else if ( i == 1 ) return y;
		else return z;
	}

	// Functions
	inline float Dot( const Point &V1 ) const
	{
		return V1.x*x + V1.y*y + V1.z*z;
	}

	inline Point CrossProduct( const Point &V2 ) const
	{
		return Point(
			y * V2.z  -  z * V2.y,
			z * V2.x  -  x * V2.z,
			x * V2.y  -  y * V2.x 	);
	}

	// Return vector rotated by the 3x3 portion of matrix m
	// (provided because it's used by bbox.cpp in article 21)
	Point RotByMatrix( const float m[16] ) const
	{
		return Point( 
			x*m[0] + y*m[4] + z*m[8],
			x*m[1] + y*m[5] + z*m[9],
			x*m[2] + y*m[6] + z*m[10] );
 	}

	// These require math.h for the sqrtf function
	float Magnitude( ) const
	{
		return sqrtf( x*x + y*y + z*z );
	}

	float Distance( const Point &V1 ) const
	{
		return ( *this - V1 ).Magnitude();	
	}

	inline void Normalize()
	{
		float fMag = ( x*x + y*y + z*z );
		if (fMag == 0) {return;}

		float fMult = 1.0f/sqrtf(fMag);            
		x *= fMult;
		y *= fMult;
		z *= fMult;
		return;
	}
};

/// XYZ double class
class DPoint
{
public:
	// Data
	double x, y, z;

	// Ctors
	DPoint( double InX, double InY, double InZ ) : x( InX ), y( InY ), z( InZ )
	{
	}
	DPoint( ) : x(0), y(0), z(0)
	{
	}

	// Operator Overloads
	inline bool operator== (const DPoint& V2) const 
	{
		return (x == V2.x && y == V2.y && z == V2.z);
	}

	inline DPoint operator+ (const DPoint& V2) const 
	{
		return DPoint( x + V2.x,  y + V2.y,  z + V2.z);
	}
	
	inline DPoint operator- (const DPoint& V2) const
	{
		return DPoint( x - V2.x,  y - V2.y,  z - V2.z);
	}
	
	inline DPoint operator- ( ) const
	{
		return DPoint(-x, -y, -z);
	}

	inline DPoint operator/ (double S ) const
	{
		double fInv = 1.0f / S;
		return DPoint (x * fInv , y * fInv, z * fInv);
	}
	
	inline DPoint operator/ (const DPoint& V2) const
	{
		return DPoint (x / V2.x,  y / V2.y,  z / V2.z);
	}
	
	inline DPoint operator* (const DPoint& V2) const
	{
		return DPoint (x * V2.x,  y * V2.y,  z * V2.z);
	}
	
	inline DPoint operator* (double S) const
	{
		return DPoint (x * S,  y * S,  z * S);
	}

	inline void operator+= ( const DPoint& V2 )
	{
		x += V2.x;
		y += V2.y;
		z += V2.z;
	}
	
	inline void operator-= ( const DPoint& V2 )
	{
		x -= V2.x;
		y -= V2.y;
		z -= V2.z;
	}

	inline double operator[] ( int i )
	{
		if ( i == 0 ) return x;
		else if ( i == 1 ) return y;
		else return z;
	}

	// Functions
	inline double Dot( const DPoint &V1 ) const
	{
		return V1.x*x + V1.y*y + V1.z*z;
	}

	inline DPoint CrossProduct( const DPoint &V2 ) const
	{
		return DPoint(
			y * V2.z  -  z * V2.y,
			z * V2.x  -  x * V2.z,
			x * V2.y  -  y * V2.x 	);
	}

	// Return vector rotated by the 3x3 portion of matrix m
	// (provided because it's used by bbox.cpp in article 21)
	DPoint RotByMatrix( const double m[16] ) const
	{
		return DPoint( 
			x*m[0] + y*m[4] + z*m[8],
			x*m[1] + y*m[5] + z*m[9],
			x*m[2] + y*m[6] + z*m[10] );
 	}

	// These require math.h for the sqrtf function
	double Magnitude( ) const
	{
		return sqrt( x*x + y*y + z*z );
	}

	double Distance( const DPoint &V1 ) const
	{
		return ( *this - V1 ).Magnitude();	
	}

	inline void Normalize()
	{
		double fMag = ( x*x + y*y + z*z );
		if (fMag == 0) {return;}

		double fMult = 1.0f/sqrt(fMag);            
		x *= fMult;
		y *= fMult;
		z *= fMult;
		return;
	}
};

/// XYZ integer class
class IPoint
{
public:
	// Data
	int x, y, z;

	// Ctors
	IPoint( int InX, int InY, int InZ ) : x( InX ), y( InY ), z( InZ )
	{
	}
	IPoint( ) : x(0), y(0), z(0)
	{
	}

	// Operator Overloads
	inline bool operator== (const IPoint& V2) const 
	{
		return (x == V2.x && y == V2.y && z == V2.z);
	}

	inline IPoint operator+ (const IPoint& V2) const 
	{
		return IPoint( x + V2.x,  y + V2.y,  z + V2.z);
	}
	
	inline IPoint operator- (const IPoint& V2) const
	{
		return IPoint( x - V2.x,  y - V2.y,  z - V2.z);
	}
	
	inline IPoint operator- ( ) const
	{
		return IPoint(-x, -y, -z);
	}

	inline IPoint operator/ (int S ) const
	{
		return IPoint (x/S , y /S, z/S);
	}
	
	inline IPoint operator/ (const IPoint& V2) const
	{
		return IPoint (x / V2.x,  y / V2.y,  z / V2.z);
	}
	
	inline IPoint operator* (const IPoint& V2) const
	{
		return IPoint (x * V2.x,  y * V2.y,  z * V2.z);
	}
	
	inline IPoint operator* (int S) const
	{
		return IPoint (x * S,  y * S,  z * S);
	}

	inline void operator+= ( const IPoint& V2 )
	{
		x += V2.x;
		y += V2.y;
		z += V2.z;
	}
	
	inline void operator-= ( const IPoint& V2 )
	{
		x -= V2.x;
		y -= V2.y;
		z -= V2.z;
	}

	inline int operator[] ( int i )
	{
		if ( i == 0 ) return x;
		else if ( i == 1 ) return y;
		else return z;
	}

	// These require math.h for the sqrtf function
	float Magnitude( ) const
	{
		return (int)(0.5+sqrtf( x*x + y*y + z*z ));
	}

	float Distance( const IPoint &V1 ) const
	{
		return ( *this - V1 ).Magnitude();	
	}

};

inline Point DP2P(DPoint dp)
{
	Point p((float)dp.x,(float)dp.y,(float)dp.z);
	return p;
}

inline DPoint P2DP(Point p)
{
	DPoint dp((double)p.x,(double)p.y,(double)p.z);
	return dp;
}

#endif
