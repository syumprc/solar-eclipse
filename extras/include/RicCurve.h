// ------------------------------- RicCurve.h --------------------------
/*! \file
Header file for the RicCurve Class
 */
using namespace std;

#include <string>
#include <vector>
#include "RicUtil.h"
#include "sisl/sisl.h"


// RicCurve class members 
/*!
\brief
This class is a front end for the Sintef Spline Library which allows creation
of spline curves from control points.

This class is a front end for the Sintef Spline Library.
(http://www.sintef.no/sisl)
The class allows creation of spline curves from control points. It is also
possible to interpolate points along a spline curve. Curve points or
control points can be written to a BrainVisa mesh file.

 */
class RicCurve
{
	public:
		
	/// variables
	DPoint	*cpnts;		///< control vertices in curve
	DPoint	*pnts;		///< interpolated vertices in curve
	int		ncntl;		///< number of control vertices
	int		npnts;		///< number of interpolated vertices
	int		ival;		///< flag indicating valid interpolated vertices
	SISLCurve* scurve;	///< sisl curve structure
	double	*knots;		///< knot vector
	int		nknots;		///< should be order + ncntl
	int		order;		///< generally 4?
	int		kind;		///< (1-4) 1-bspline 2-nurbs 3-bezier 4-rat bezier
	double	*sislcpnts;	///< double pointer to cpnts - to fool sisl
	
	// constructors
	RicCurve(void);
	RicCurve(int nc);
	RicCurve(int nc, int ord, int kind);
	RicCurve(int nc, Point *cp);
	RicCurve(int nc, DPoint *cp);
	~RicCurve();
	
	// member functions
	int InitCurve();
	int Interpolate(int n);
	void WriteCurveGo(string filename);
	void WriteControlAsMesh(string filename);
	void WriteCurveAsMesh(string filename);
	int SortAlongLength();
	double CurveLength();
	
};

// sisl Go header-related info.
#define HEADER_SIZE 4
#define CURVE_INSTANCE_TYPE 100
#define SURFACE_INSTANCE_TYPE 200
#define POINTCLOUD_INSTANCE_TYPE 400
#define MAJOR_VERSION 1
#define MINOR_VERSION 0

// Sisl related functions for dealing with GO files

/*!
 * @ brief
 * does something with a GO file header
 * @param is - stream to read from
 * @return
 */
inline int determine_go_instance_type(istream& is) 
{
	int result;
	is >> result;
	for (int dummy, i = 1; i < HEADER_SIZE; ++i)
		is >> dummy;
	return result;
}

/*!
 * @brief
 * reads knots from a GO file
 *
 * @param is - input stream
 * @param n - number of vertices
 * @param k - order of curve
 * @param knots - array of knots
 */
inline void read_go_basis(istream& is, int& n, int& k, vector<double>& knots) 
{
	is >> n >> k;
	knots.resize(n + k);
	for (int i = 0; i < n + k; ++i) {
		is >> knots[i];
	}
}

/*!
 * @brief
 * writes knots to a GO file
 * @param os - output stream
 * @param n - number of vertices
 * @param k - order of curve
 * @param knots - array of knots
 */
inline void write_go_basis(ostream& os, const int&n, const int& k, const double* knots)
{
	os << n << ' ' << k << '\n';
	for (int i = 0; i < n + k; ++i) {
		os << knots[i] << ' ';
	}
	os << '\n';
}
