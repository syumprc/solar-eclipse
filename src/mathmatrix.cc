/*
 * matrixmath implements matrix math functions using Eigen
 *
 * Written by Charles Peterson beginning on January 16, 2015
 * Copyright (c) 2015 Texas Biomedical Research Institute
 */


// On older compilers, TR1 is required for latest C++ random
// #define TR1 1

#ifdef TR1
#define STDPRE std::tr1
#else
#define STDPRE std
#endif

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "solar.h"
// tcl.h from solar.h
#include "safelib.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#ifdef TR1
#include <tr1/random>
#else
#include <random>
#endif

using std::ofstream;
using namespace Eigen;

// Prefix for all math matrix handles, followed by index number
#define MathMatrix_PREFIX ".mm."
#define MathMatrix_PREFIX_LENGTH 4

int MathMatrixDebug = 0;
std::vector<Eigen::MatrixXd*> MathMatrixV(0);
int LastIndex = 0;
int EvaluesIndex = -1;
int EvectorsIndex = -1;
std::ptrdiff_t Min_i = INT_MIN;
std::ptrdiff_t Min_j = INT_MIN;
std::ptrdiff_t Max_i = INT_MIN;
std::ptrdiff_t Max_j = INT_MIN;
Eigen::MatrixXd* FortMatrix = 0;

class NumberSequence
{
public:
    NumberSequence () {position=0;}
    std::vector<double> numbers;
    int position;

    double next() {double rvalue = numbers[position];
	if(++position>=numbers.size())position=0;
	return rvalue;}

    void push_back (double value) {numbers.push_back(value);}
};

class SelectedColumn
{
public:
    int number;
    char* name;
    NumberSequence* sequence;
    SelectedColumn(const char* _name) {number=0; name=strdup(_name); 
	sequence=0;}
    SelectedColumn(int _number) {number=_number; name=0; sequence=0;}
    SelectedColumn(NumberSequence* seq) {sequence=seq; number=0; name=0;}
    ~SelectedColumn() {if (name) free (name); if (sequence) delete sequence;}
};

std::vector<SelectedColumn*> SelectedColumns;

// AdvancedRandom methods abstract the random layer
//   Now using Mersenne twister as best quality for fast generator


STDPRE::mt19937 mt;
int RandomSeeded = 0;

void AdvancedRandomSeed (int option)
{
    unsigned int reseeding = Option::get_int("ShuffleReseeding");
    if (reseeding == 0) {

// option 0 is "default all the way"
// mt19937 seeds itself at the start of program
// gives consistent result first time, then stochastic
// this is least useful option

	RandomSeeded = 0;
	return;
    }

    if (reseeding == 1)
    {
	mt.seed(5489u);  // this is the default seed for mt19937
                         // option 1 applys it to every shuffle
                         // for consistent results every time
                         // also best for comparing models same sample
                         // RECOMMENDED !
	RandomSeeded = 1;
    }
    else if (reseeding == -1)
    {
	mt.seed(time(NULL));  // reseeded stochastically each time
	RandomSeeded = -1;
    }
    else if (reseeding == -2)
    {
	if (!RandomSeeded != -2)
	{
	    mt.seed(time(NULL));  // stochastic first time, then best distr
	}
	RandomSeeded = -2;
    }
    else
    {
	mt.seed(reseeding);
	RandomSeeded = reseeding;
    }
}

int AdvancedRandomDraw (int toprange)  // 0 <= i <= toprange
{
    unsigned int mtval = mt();

// Sadly, C++ uniform_int_distribution is unavailable on gcc 4.4.1 or lower
// An alternative distribution-corrected approach follows
// Stack Overflow: generating-a-uniform-distribution-of-integers-in-c

//    return mtval % (toprange+1);   // fast but inaccurate cheat

    unsigned int fullrange = toprange + 1;
    unsigned int copies = mt.max() / fullrange;
    unsigned int limit = fullrange * copies;
    while ( mtval >= limit)
    {
	mtval = mt();
    }
    int result = mtval / copies;
    return result;
}

// Utility functions for MathMatrix
// Given MathMatrix id, return pointer to Eigen matrix

MatrixXd* MathMatrixGet (char* mmid, const char* fname, Tcl_Interp* interp)
{
    std::string Mmid(mmid);

    if (strncmp (mmid, MathMatrix_PREFIX, MathMatrix_PREFIX_LENGTH))
    {
	std::string errmess = "Invalid MathMatrixID or number \'" + Mmid + "\' or invalid matrix command";
	RESULT_BUF (errmess.c_str());
	return 0; // Caller checks for zero and returns TCL_ERROR;
    }
    errno = 0;
    char* eon;
    char* nptr = mmid + strlen(MathMatrix_PREFIX);
    long lindex = strtol (nptr, &eon, 10);
    if (errno || *eon!=0 || lindex < 1)
    {
	std::string errmess = "Invalid MathMatrixID: " + Mmid;
	RESULT_BUF (errmess.c_str());
	return 0;
    }
    if (lindex > MathMatrixV.size())
    {
	std::string errmess = "MathMatrixID " + Mmid + " not allocated";
	RESULT_BUF (errmess.c_str());
	return 0;
    }
    LastIndex = (int) lindex - 1;  // back to 0 based to access c++ vector
    if (MathMatrixV[LastIndex] == 0)
    {
	std::string errmess = "Empty (Deleted) MathMatrix: " + Mmid;
	RESULT_BUF (errmess.c_str());
	return 0;
    }
    return MathMatrixV[LastIndex];
}

// Power command for MathMatrix
//   power command in solar.tcl redirects here if first arg is mathmatrix

extern "C" int MatrixPowerCmd (ClientData clientData, Tcl_Interp *interp,
			 int argc, char* argv[])
{
    if (argc < 3)
    {
	RESULT_LIT ("power <MathMatrix> <integer>");
	return TCL_ERROR;
    }

// Get Matrix

    MatrixXd* m1 = MathMatrixGet (argv[1], "Power", interp);
    if (!m1)
    {
	RESULT_LIT ("matrixpower: Invalid matrix specified");
	return TCL_ERROR;
    }
    int rows = m1->rows();
    int cols = m1->cols();

// Get Integer

    errno=0;
    char* eptr=0;
    double power = strtod (argv[2], &eptr);
    if (errno || (*eptr!=0))
    {
	RESULT_LIT ("matrixpower: invalid power integer specified");
    }

// Allocate new matrix

    MatrixXd* newmat = new MatrixXd (rows, cols);

// transform while copying

    for (int i = 0; i < rows; i++)
    {
	for (int j = 0; j < cols; j++)
	{
	    {
		(*newmat)(i,j) = pow ( (*m1)(i,j), power);
	    }
	}
    }
    MathMatrixV.push_back (newmat);
    char namebuf[64];
    sprintf (namebuf, "%s%d", MathMatrix_PREFIX, MathMatrixV.size());
    RESULT_BUF (namebuf);
    return TCL_OK;
}

// Overloaded Matrix and Vector (1D Matrix) and Scalar multiplication
// Times has options -e (elementwise)
extern "C" int TimesCmd (ClientData clientData, Tcl_Interp *interp,
			 int argc, char* argv[])
{
    int option = 0;
    int first_index = 1;
    bool scalar_1 = false;
    bool scalar_2 = false;
    double s1, s2;

    if (argc < 3)
    {
	RESULT_LIT (
	    "Usage: Times <option> <matrix-or-scalar> <matrix-or-scalar>");
	return TCL_ERROR;
    }
    if (!Strcmp(argv[1], "-e"))
    {
	option = 1;
	first_index++;
	if (argc < 4)
	{
	    RESULT_LIT (
		"Usage: Times <option> <matrix-or-scalar> <matrix-or-scalar>");
	    return TCL_ERROR;
	}
    }

    char* eptr;

    MatrixXd* m1 = MathMatrixGet (argv[first_index], "Times", interp);
    if (!m1)
    {
	errno=0;
	eptr=0;
	s1 = strtod (argv[first_index], &eptr);
	if (errno || (*eptr!=0))
	{
// result set by MathMatrixGet
	    return TCL_ERROR;
	}
	RESULT_LIT ("");
	scalar_1 = true;
    }

    MatrixXd* m2 = MathMatrixGet (argv[first_index+1], "Times", interp);
    if (!m2)
    {
	errno=0;
	eptr=0;
	s2 = strtod (argv[first_index+1], &eptr);
	if (errno || (*eptr!=0))
	{
// result set by MathMatrixGet
	    return TCL_ERROR;
	}
	RESULT_LIT ("");
	scalar_2 = true;
    }

    int mindex;
    if (!scalar_1 && !scalar_2)
    {
	if (option==0)  // standard matrix multiply
	{
	    if (m1->cols() != m2->rows())
	    {
		RESULT_LIT (
		    "Matrix dimensions not suited to multiply");
		return TCL_ERROR;
	    }
	    MatrixXd* mr = new MatrixXd(m1->cols(),m2->rows());
	    mindex = MathMatrixV.size() - 1;
	    if (MathMatrixDebug) printf ("Beginning multiplication\n");
	    *mr = *m1 * *m2;
	    if (MathMatrixDebug) printf ("Finished multiplication\n");
	    MathMatrixV.push_back(mr);
	}
	else if (option==1)  // elementwise matrix multiply
	{
	    if (m1->rows() != m2->rows() || m1->cols() != m2->cols())
	    {
		RESULT_LIT (
		    "Matrix dimensions not suited to ewise multiply");
		return TCL_ERROR;
	    }
	    MathMatrixV.push_back(new MatrixXd(m1->rows(),m2->cols()));
	    mindex = MathMatrixV.size() - 1;
	    *MathMatrixV[mindex] = m1->cwiseProduct(*m2);
	}
    }
    else  // scalar 1 and/or 2
    {
	int rows, cols;
	if (scalar_1 && scalar_2)
	{
	    double result = s1 * s2;
	    char numbuf[64];
	    char showformat[64];
	    sprintf (showformat,"%%-.%ig",Option::get_int("MatrixNumberFormat"));
	    sprintf (numbuf, showformat,result);
	    RESULT_BUF (numbuf);
	    return TCL_OK;  // scalar result
	}

	if (scalar_1)
	{
	    rows = m2->rows();
	    cols = m2->cols();
	}
	else
	{
	    rows = m1->rows();
	    cols = m1->cols();
	}
	MathMatrixV.push_back(new MatrixXd(rows, cols));
	int mindex = MathMatrixV.size() - 1;

	if (scalar_1)
	{
	    *MathMatrixV[mindex] = *m2 * s1;
	}
	else
	{
	    *MathMatrixV[mindex] = *m1 * s2;
	}
    }

    char namebuf[64];
    sprintf (namebuf,"%s%d",MathMatrix_PREFIX,MathMatrixV.size());
    RESULT_BUF (namebuf);
    return TCL_OK;
}

extern "C" int ConcatenateCmd (ClientData clientData, Tcl_Interp *interp,
			 int argc, char* argv[])
{
    if (argc < 2)
    {
	RESULT_LIT ("Usage: concatenate <vertical | horizontal> [<matrix>]+");
	return TCL_ERROR;
    }

    int direction;

    if (!strcmp(argv[1],"1"))
    {
	direction = 1;
    }
    else if (!strcmp(argv[1],"2"))
    {
	direction = 2;
    }
    else if (!strcmp(argv[1],"vertical"))
    {
	direction = 1;
    }
    else if (!strcmp(argv[1],"horizontal"))
    {
	direction = 2;
    }
    else
    {
	RESULT_LIT ("Concatenation direction must be 1 or 2");
	return TCL_ERROR;
    }

// Row or Column Append?
//   User specification overrides, error if can't be done
//   Otherwise, if first matrix is a vector, subsequent vectors are appended
//     to extend it in it's unitary dimension.
//   If the first matrix is a non-vector, the second vector extends it in
//     it's unitary dimension.

    int rows, total_rows;
    int cols, total_cols;
    MatrixXd* mi[argc];

    int mati = -1;
    int argi;
    for (argi = 2; argi < argc; argi++)
    {
	mati++;
	mi[mati] = MathMatrixGet (argv[argi], "concatenate", interp);
	if (!mi[mati]) return TCL_ERROR;
	if (argi == 2)
	{
	    rows = mi[mati]->rows();
	    cols = mi[mati]->cols();
	    total_rows = rows;
	    total_cols = cols;
	}
	else
	{
	    if (direction==1)
	    {
		if (cols != mi[mati]->cols())
		{
		    std::string errmess("Column count mismatch in matrix number ");
		    errmess += argi-1;
		    RESULT_BUF(errmess.c_str());
		    return TCL_ERROR;
		}
		total_rows += mi[mati]->rows();
	    }
	    if (direction==2)
	    {
		if (rows != mi[mati]->rows())
		{
		    std::string errmess("Row count mismatch in matrix number ");
		    errmess += argi-1;
		    RESULT_BUF(errmess.c_str());
		    return TCL_ERROR;
		}
		total_cols += mi[mati]->cols();
	    }
	}
    }

// Allocate new matrix

    MatrixXd* mout = new MatrixXd(total_rows, total_cols);
    int torow = 0;  // or column if direction is 2

// Concatenate

    for (mati = 0; mati < argc-2; mati++)
    {
	if (direction==1)
	{
	    int rowshere = mi[mati]->rows();
	    for (int fromrow = 0; fromrow < rowshere; fromrow++)
	    {
		mout->block(torow,0,1,total_cols) <<
		    mi[mati]->block(fromrow,0,1,total_cols);
		torow++;
	    }
	}
	else
	{
	    int rowshere = mi[mati]->cols();
	    for (int fromrow = 0; fromrow < rowshere; fromrow++)
	    {
		mout->block(0,torow,total_rows,1) <<
		    mi[mati]->block(0,fromrow,total_rows,1);
		torow++;
	    }
	}
    }
    MathMatrixV.push_back (mout);
    char namebuf[64];
    sprintf (namebuf,"%s%d",MathMatrix_PREFIX,MathMatrixV.size());
    RESULT_BUF (namebuf);
    return TCL_OK;
}

extern "C" int InsertCmd (ClientData clientData, Tcl_Interp *interp,
			 int argc, char* argv[])
{
    if (argc != 5)
    {
	RESULT_LIT ("Usage: insert <matrix> <row> <col> <newvalue>");
	return TCL_ERROR;
    }

    MatrixXd* m1 = MathMatrixGet (argv[1], "insert", interp);
    if (!m1) return TCL_ERROR;

    char* eptr;
    errno = 0;
    long lrow = strtol (argv[2], &eptr, 10);
    if (errno || *eptr != '\0' || lrow < 1 || lrow > m1->rows())
    {
	std::string errmess("Invalid matrix row: ");
	errmess += argv[2];
	RESULT_BUF(errmess.c_str());
	return TCL_ERROR;
    }
    int row = (int) lrow - 1;
    long lcol = strtol (argv[3], &eptr, 10);
    if (errno || *eptr != '\0' || lcol < 1 || lcol > m1->cols())
    {
	std::string errmess("Invalid matrix col: ");
	errmess += argv[3];
	RESULT_BUF(errmess.c_str());
	return TCL_ERROR;
    }
    int col = (int) lcol - 1;

    errno = 0;
    double value = strtod (argv[4], &eptr);
    if (errno || *eptr != '\0')
    {
	std::string errmess("Invalid matrix value: ");
	errmess += argv[4];
	RESULT_BUF(errmess.c_str());
	return TCL_ERROR;
    }

// Do it

    (*m1)(row,col) = value;

    RESULT_LIT ("");
    return TCL_OK;
}

// ols uses the normal equations which are fastest (solve uses more robust QR)
extern "C" int OlsCmd (ClientData clientData, Tcl_Interp *interp,
			 int argc, char* argv[])
{
    int first_index = 1;

    if (argc < 3)
    {
	RESULT_LIT ("Usage: ols <y-Trait-Matrix> <X-Design-Matrix>");
	return TCL_ERROR;
    }

    MatrixXd* my = MathMatrixGet (argv[first_index], "OLS", interp);
    if (!my) return TCL_ERROR;

    MatrixXd* mx = MathMatrixGet (argv[first_index+1], "OLS", interp);
    if (!mx) return TCL_ERROR;

    if (my->rows() != mx->rows())
    {
	RESULT_LIT (
	    "Rows in y and X must be equal");
	return TCL_ERROR;
    }
    MathMatrixV.push_back(new MatrixXd(1,mx->cols()));
    int mindex = MathMatrixV.size() - 1;

// ?don't know where I got this
// not sure why I can't use the ldlt solve directly
    *MathMatrixV[mindex] = 
	(mx->transpose()*(*mx)).ldlt().solve(mx->transpose()*(*my));

// this documented version doesn't work
//    *MathMatrixV[mindex] =
//	mx->ldlt().solve(*my);

    char namebuf[64];
    sprintf (namebuf,"%s%d",MathMatrix_PREFIX,MathMatrixV.size());
    RESULT_BUF (namebuf);
    return TCL_OK;
}

extern "C" int SolveCmd (ClientData clientData, Tcl_Interp *interp,
			 int argc, char* argv[])
{
    int first_index = 1;

    if (argc < 3)
    {
	RESULT_LIT ("Usage: solve <Y-Trait-Matrix> <X-Design-Matrix>");
	return TCL_ERROR;
    }

    MatrixXd* my = MathMatrixGet (argv[first_index], "Solve", interp);
    if (!my) return TCL_ERROR;

    MatrixXd* mx = MathMatrixGet (argv[first_index+1], "Solve", interp);
    if (!mx) return TCL_ERROR;

    if (my->rows() != mx->rows())
    {
	RESULT_LIT (
	    "Rows in y and X must be equal");
	return TCL_ERROR;
    }
    MathMatrixV.push_back(new MatrixXd(1,mx->cols()));
    int mindex = MathMatrixV.size() - 1;

    *MathMatrixV[mindex] = 
	mx->colPivHouseholderQr().solve(*my);

    char namebuf[64];
    sprintf (namebuf,"%s%d",MathMatrix_PREFIX,MathMatrixV.size());
    RESULT_BUF (namebuf);
    return TCL_OK;
}

extern "C" int EigenCmd (ClientData clientData, Tcl_Interp *interp,
			 int argc, char* argv[])
{
    int first_index = 1;

    if (argc == 2)
    {

	MatrixXd* m1 = MathMatrixGet (argv[first_index], "EValues", interp);
	if (!m1) return TCL_ERROR;

	if (m1->rows() != m1->cols())
	{
	    RESULT_LIT ("Eigen decomposition for square matrix only");
	    return TCL_ERROR;
	}


	EigenSolver<MatrixXd> es(*m1);

	if (MathMatrixDebug) printf ("Matrix is decomposed\n");

// Push evalues and evectors matrixes into main vector, then
// return requested matrix    

	MathMatrixV.push_back (new MatrixXd (m1->rows(), m1->cols()));
	EvaluesIndex = MathMatrixV.size() - 1;
	*MathMatrixV[EvaluesIndex] = es.pseudoEigenvalueMatrix();

	MathMatrixV.push_back (new MatrixXd (m1->rows(), m1->cols()));
	EvectorsIndex = MathMatrixV.size() - 1;
	*MathMatrixV[EvectorsIndex] = es.pseudoEigenvectors();
    }

    int requested;
    if (!Strcmp (argv[0],"evalues"))
    {
	requested = EvaluesIndex;
    }
    else
    {
	requested = EvectorsIndex;
    }

    if (requested == -1)
    {
	RESULT_LIT ("Must use evalues or evectors with matrix argument first");
	return TCL_ERROR;
    }
    char namebuf[64];
    sprintf (namebuf,"%s%d",MathMatrix_PREFIX,requested+1);
    RESULT_BUF (namebuf);
    return TCL_OK;
}

// Variant of max command returns new matrix (or original, if unchanged)
//   with each element the max of current value or specified scalar
//   most often used for to truncate negative values to zero

extern "C" int MaxMatrix (Tcl_Interp* interp, MatrixXd* m1,
			      int argc, char* argv[])
{
    if (argc > 3)
    {
	RESULT_LIT ("Usage: max <MathMatrix> [<scalar>]");
	return TCL_ERROR;
    }

// Get scalar

    errno=0;
    char* eptr=0;
    double use_min = strtod (argv[2], &eptr);
    if (errno || *eptr != '\0')
    {
	RESULT_LIT ("Invalid integer for max");
	return TCL_ERROR;
    }

// Get matrix min value and see if new matrix needed

    int mrow;
    int mcol;
    double actual_min = m1->minCoeff(&mrow,&mcol);
    if (actual_min >= use_min)
    {
	char* mname = Strdup (argv[1]);
	RESULT_BUF (mname);
	free (mname);
	return TCL_OK;
    }

// Make new matrix and copy values while checking

    MatrixXd* m2 = new MatrixXd(m1->rows(),m1->cols());
    for (int i=0; i < m1->rows(); i++)
    {
	for (int j=0; j < m1->cols(); j++)
	{
	    double value = (*m1)(i,j);
	    if (value >= use_min)
	    {
		(*m2)(i,j) = (*m1)(i,j);
	    }
	    else
	    {
		(*m2)(i,j) = use_min;
	    }
	}
    }

    MathMatrixV.push_back (m2);
    char namebuf[64];
    sprintf (namebuf,"%s%d",MathMatrix_PREFIX,MathMatrixV.size());
    RESULT_BUF (namebuf);
    return TCL_OK;
}

//
// Calculate mean, sum, and similar operations on a vector
//   (note: for sum this is invoked through PlusCmd)

extern "C" int MeanSumCmd (ClientData clientData, Tcl_Interp *interp,
			 int argc, char* argv[])
{
    bool input_vector = false;
    bool do_where = false;

    if (argc < 2 && 0!=Strcmp(argv[0],"max"))
    {
	RESULT_LIT ("Usage: <unary-operation> [<MathMatrix>|-where]");
	return TCL_ERROR;
    }
    if (!Strcmp(argv[1],"-where"))
    {
	if (!Strcmp (argv[0], "min"))
	{
	    if (Min_i == INT_MIN)
	    {
		RESULT_LIT ("Must do min on matrix first without -where");
		return TCL_ERROR;
	    }
	    char buf[128];
	    sprintf (buf, "{%d %d}",Min_i+1, Min_j+1);
	    RESULT_BUF (buf);
	    return TCL_OK;
	}
	else if (!Strcmp (argv[0], "max"))
	{
	    if (Max_i == INT_MIN)
	    {
		RESULT_LIT ("Must do max on matrix first without -where");
		return TCL_ERROR;
	    }
	    char buf[128];
	    sprintf (buf, "{%d %d}",Max_i+1, Max_j+1);
	    RESULT_BUF (buf);
	    return TCL_OK;
	}
    }

    MatrixXd* m1 = MathMatrixGet (argv[1], "mean", interp);
    if (!m1) return TCL_ERROR;

    if (m1->cols() == 1)
    {
	input_vector = true;
    }
    else if (m1->rows() == 1)
    {
	input_vector = true;
    }

    double result;
    if (!Strcmp (argv[0], "mean"))
    {
	result = m1->mean();
    }
    else if (!Strcmp (argv[0], "min"))
    {
	result = m1->minCoeff(&Min_i, &Min_j);
    }
    else if (!Strcmp (argv[0], "max"))
    {
	if  (argc>2) return MaxMatrix (interp,m1,argc,argv);

	result = m1->maxCoeff(&Max_i, &Max_j);
    }
    else if (!Strcmp (argv[0], "plus"))
    {
	if (input_vector)
	{
	    result = m1->sum();
	}
	else
	{
	    RESULT_LIT ("unary plus not yet supported on non-vector");
	    return TCL_ERROR;
	}
    }
    else
    {
	RESULT_LIT ("Invalid matrix command");
	return TCL_ERROR;
    }
    char showformat[64];
    sprintf (showformat,"%%-.%ig",Option::get_int("MatrixNumberFormat"));
    char numbuf[64];
    sprintf (numbuf, showformat, result);
    RESULT_BUF (numbuf);
    return TCL_OK;
}




// Matrix plus is always elementwise
extern "C" int PlusCmd (ClientData clientData, Tcl_Interp *interp,
			 int argc, char* argv[])
{
    int first_index = 1;
    bool scalar_1 = false;
    bool scalar_2 = false;
    double s1, s2;
    char* eptr;

    if (argc < 3)
    {
	return MeanSumCmd (clientData, interp, argc, argv);
    }

    MatrixXd* m1 = MathMatrixGet (argv[first_index], "plus", interp);
    if (!m1)
    {
	errno=0;
	eptr=0;
	s1 = strtod (argv[first_index], &eptr);
	if (errno || (*eptr!=0))
	{
// result set by MathMatrixGet
	    return TCL_ERROR;
	}
	RESULT_LIT ("");
	scalar_1 = true;
    }

    MatrixXd* m2 = MathMatrixGet (argv[first_index+1], "plus", interp);
    if (!m2)
    {
	errno=0;
	eptr=0;
	s2 = strtod (argv[first_index+1], &eptr);
	if (errno || (*eptr!=0))
	{
// result set by MathMatrixGet
	    return TCL_ERROR;
	}
	RESULT_LIT ("");
	scalar_2 = true;
    }

    int mindex;
    
    if (!scalar_1 && !scalar_2)
    {
	if (m1->rows() != m2->rows() || m1->cols() != m2->cols())
	{
	    RESULT_LIT (
		"Matrix dimensions not equal as required for plus");
	    return TCL_ERROR;
	}
	MathMatrixV.push_back(new MatrixXd(m1->rows(),m1->cols()));
	mindex = MathMatrixV.size() - 1;

	*MathMatrixV[mindex] = *m1 + *m2;
    }
    else if (scalar_1 && !scalar_2)
    {
	MathMatrixV.push_back(new MatrixXd(m2->rows(),m2->cols()));
	mindex = MathMatrixV.size() - 1;

// Maddeningly, Eigen doesn't have scalar add and subtract for matrices but
//   would require conversion into and out of array.  Looping here seems
//     at least equally good.

	int i,j;
	MatrixXd* mo = MathMatrixV[mindex];

	for (i=0; i < m2->rows(); i++)
	{
	    for (j=0; j < m2->cols(); j++)
	    {
		(*mo)(i,j) = s1 + (*m2)(i,j);
	    }
	}
    }
    else if (!scalar_1 && scalar_2)
    {
	MathMatrixV.push_back(new MatrixXd(m1->rows(),m1->cols()));
	mindex = MathMatrixV.size() - 1;

	int i,j;
	MatrixXd* mo = MathMatrixV[mindex];

	for (i=0; i < m1->rows(); i++)
	{
	    for (j=0; j < m1->cols(); j++)
	    {
		(*mo)(i,j) = (*m1)(i,j) + s2;
	    }
	}
    }
    else // two scalars
    {
	char numbuf[64];
	char showformat[65];

	double result = s1 + s2;
	sprintf (showformat, "%%-.%ig",Option::get_int("MatrixNumberFormat"));
	sprintf (numbuf, showformat, result);
	RESULT_BUF (numbuf);
	return TCL_OK;
    }

    char namebuf[64];
    sprintf (namebuf,"%s%d",MathMatrix_PREFIX,MathMatrixV.size());
    RESULT_BUF (namebuf);
    return TCL_OK;
}

// Matrix minus is always elementwise
extern "C" int MinusCmd (ClientData clientData, Tcl_Interp *interp,
			 int argc, char* argv[])
{
    int first_index = 1;
    bool scalar_1 = false;
    bool scalar_2 = false;
    double s1, s2;
    char* eptr;

    if (argc < 3)
    {
	RESULT_LIT ("Usage: minus <MathMatrix1> <MathMatrix2>");
	return TCL_ERROR;
    }

    MatrixXd* m1 = MathMatrixGet (argv[first_index], "Minus", interp);
    if (!m1)
    {
	errno=0;
	eptr=0;
	s1 = strtod (argv[first_index], &eptr);
	if (errno || (*eptr!=0))
	{
// result set by MathMatrixGet
	    return TCL_ERROR;
	}
	RESULT_LIT ("");
	scalar_1 = true;
    }

    MatrixXd* m2 = MathMatrixGet (argv[first_index+1], "Minus", interp);
    if (!m2)
    {
	errno=0;
	eptr=0;
	s2 = strtod (argv[first_index+1], &eptr);
	if (errno || (*eptr!=0))
	{
// result set by MathMatrixGet
	    return TCL_ERROR;
	}
	RESULT_LIT ("");
	scalar_2 = true;
    }

    int mindex;

    if (!scalar_1 && !scalar_2)
    {
	if (m1->rows() != m2->rows() || m1->cols() != m2->cols())
	{
	    RESULT_LIT (
		"Matrix dimensions not equal as required for minus");
	    return TCL_ERROR;
	}
	MathMatrixV.push_back(new MatrixXd(m1->rows(),m1->cols()));
	mindex = MathMatrixV.size() - 1;
	*MathMatrixV[mindex] = *m1 - *m2;
    }
    else if (scalar_1 && !scalar_2)
    {
	MathMatrixV.push_back(new MatrixXd(m2->rows(),m2->cols()));
	mindex = MathMatrixV.size() - 1;

// Maddeningly, Eigen doesn't have scalar add and subtract for matrices but
//   would require conversion into and out of array.  Looping here seems
//     at least equally good.

	int i,j;
	MatrixXd* mo = MathMatrixV[mindex];

	for (i=0; i < m2->rows(); i++)
	{
	    for (j=0; j < m2->cols(); j++)
	    {
		(*mo)(i,j) = s1 - (*m2)(i,j);
	    }
	}
    }
    else if (!scalar_1 && scalar_2)
    {
	MathMatrixV.push_back(new MatrixXd(m1->rows(),m1->cols()));
	mindex = MathMatrixV.size() - 1;

	int i,j;
	MatrixXd* mo = MathMatrixV[mindex];

	for (i=0; i < m1->rows(); i++)
	{
	    for (j=0; j < m1->cols(); j++)
	    {
		(*mo)(i,j) = (*m1)(i,j) - s2;
	    }
	}
    }
    else  // two scalars
    {
	char numbuf[64];
	char showformat[65];

	double result = s1 - s2;
	sprintf (showformat, "%%-.%ig",Option::get_int("MatrixNumberFormat"));
	sprintf (numbuf, showformat, result);
	RESULT_BUF (numbuf);
	return TCL_OK;
    }

    char namebuf[64];
    sprintf (namebuf,"%s%d",MathMatrix_PREFIX,MathMatrixV.size());
    RESULT_BUF (namebuf);
    return TCL_OK;
}


extern "C" int ShuffleCmd (ClientData clientData, Tcl_Interp *interp,
			 int argc, char* argv[])
{
    static bool seeded = false;

    if (argc < 2 || argc > 3)
    {
	RESULT_LIT ("Usage: shuffle <vector>");
	return TCL_ERROR;
    }
    MatrixXd* m1 = MathMatrixGet (argv[1], "shuffle", interp);
    if (!m1) return TCL_ERROR;

    int width = 0;
    int basevector = 0;
    if (argc == 3)
    {
	errno = 0;
	char* eptr;
	long lwidth = strtol (argv[2], &eptr, 10);
	if (errno || *eptr != '\0' || lwidth < INT_MIN || lwidth > INT_MAX)
	{
	    RESULT_LIT ("shuffle: invalid dimension specified");
	    return TCL_ERROR;
	}
	if (lwidth > 0) basevector = 1;
	width = abs(lwidth);
    }

// if width is 0, we make vector from vector
// if width is 1+, we make matrix of original vector plus <width> shuffled cols
// if width is <0, we make matrix of -<width> shuffled cols but no original
// width is now abs(width) and basevector toggles the original col



    int rows = m1->rows();
    int cols = m1->cols();

    bool row_vector = (rows == 1);
    int nele;
    if (!row_vector)
    {
	if (cols != 1)
	{
	    RESULT_LIT ("shuffle: matrix must be row or column vector");
	    return TCL_ERROR;
	}
	nele = rows;
    }
    else
    {
	nele = cols;
    }

// If initial vector is row_vector, a shuffled row vector, or matrix of
//    shuffled row vectors is produced.

    int bwidth=1;
    if (width!=0)
    {
	bwidth = width + basevector;
	if (row_vector)
	{
	    rows = bwidth;
	}
	else
	{
	    cols = bwidth;
	}
    }

    MatrixXd* mout = new MatrixXd(rows,cols);

// Seed

    extern bool random_number_generator_seeded;
    if (!seeded && !random_number_generator_seeded)
    {
	srand48(time(NULL));
    }
    seeded = true;

    for (int iwidth=0; iwidth < bwidth; iwidth++)
    {
	if (basevector && iwidth == 0)
	{

// copy basevector into first row or column

	    for (unsigned int i = 0; i < nele; i++)
	    {
		if (row_vector)
		{
		    (*mout)(iwidth,i) = (*m1)(0,i);
		}
		else
		{
		    (*mout)(i,iwidth) = (*m1)(i,0);
		}
	    }
	}
	else
	{

// shuffle basevector into this row (if row vector) or column

// Knuth aka Fisher-Yates Shuffle
//   Inside-Out version from Wikipedia 6/2015

	    unsigned int j;
	    for (unsigned int i = 0; i < nele; i++)
	    {
		if (i == 0)
		{
		    j = 0;
		}
		else
		{
		    double drand = drand48();
		    while (drand >= 1) drand = drand48();  // we need < 1.0

// Generate a random number 0 <= j <= i

		    j = (unsigned int) floor (drand*(i+1));
		}

		if (j!=i)
		{
		    if (row_vector)
		    {
			(*mout)(iwidth,i) = (*mout)(iwidth,j);
		    }
		    else
		    {
			(*mout)(i,iwidth) = (*mout)(j,iwidth);
		    }
		}
		if (row_vector)
		{
		    (*mout)(iwidth,j) = (*m1)(0,i);
		}
		else
		{
		    (*mout)(j,iwidth) = (*m1)(i,0);
		}
	    }
	}
    }
    MathMatrixV.push_back (mout);
    char namebuf[64];
    sprintf (namebuf,"%s%d",MathMatrix_PREFIX,MathMatrixV.size());
    RESULT_BUF (namebuf);
    return TCL_OK;
}

// permutef produces the permutation Fp from unshuffled F, XB, and Hx
//   to save memory/cache and data movement, no intermediate matrices are used.
//   instead columns of Fp are created one at a time by shuffling columns of F,
//   adding XB, multiplying by rows of Hx

//  [permutef F XB Hx nP]

extern "C" int PermuteFCmd (ClientData clientData, Tcl_Interp *interp,
			 int argc, char* argv[])
{
    static bool seeded = false;

    if (argc != 5)
    {
	RESULT_LIT ("Usage: permutef <F> <XB> <Hx> <nP>");
	return TCL_ERROR;
    }
    MatrixXd* mF = MathMatrixGet (argv[1], "permuteF", interp);
    if (!mF) return TCL_ERROR;
    MatrixXd* mXB = MathMatrixGet (argv[2], "permuteF", interp);
    if (!mXB) return TCL_ERROR;
    MatrixXd* mHx = MathMatrixGet (argv[3], "permuteF", interp);
    if (!mHx) return TCL_ERROR;

// Get nP

    errno = 0;
    char* eptr;
    long np = strtol (argv[4], &eptr, 10);
    if (errno || *eptr != '\0' || np < INT_MIN || np > INT_MAX)
    {
	RESULT_LIT ("permutef: invalid nP");
	return TCL_ERROR;
    }
    int nP = (int) np;

    int rows = mHx->rows();
    int incols = mHx->cols();
    int outcols = nP + 1;

    MatrixXd* mFp = new MatrixXd(rows,outcols);


// First column of output matrix mFp is Hx * F+XB

    int i,j;

    {

// Cache f and xb

	double f[rows];
	double xb[rows];
	for (i = 0; i < rows; i++)
	{
	    f[i]  = (*mF)(i,0);
	    xb[i] = (*mXB)(i,0);
	}
	
// Compute first column

	for (i = 0; i < rows; i++)
	{
	    double accum = 0;
	    for (j = 0; j < incols; j++)
	    {
		accum += (*mHx)(i,j) * (f[i] + xb[i]);
	    }
	    (*mFp)(i,0) = accum;
	}

// Remaining columns of output matrix mFp have shuffled F+XB
// We first generate shuffle index, then do multiply as above

// Seed

	AdvancedRandomSeed (1); // fixed seed each time

	for (int jout = 1; jout < outcols; jout++)
	{

// create shuffle index vector using
// Knuth aka Fisher-Yates Shuffle
//   Inside-Out version from Wikipedia 6/2015

//	    if (jout % 100 == 0) printf ("Computing column %d\n",jout);

	    int shuffle[incols];
	    int r;

	    for (i = 0; i < incols; i++)
	    {
		if (i == 0)
		{
		    r = 0;
		}
		else
		{
		    r = AdvancedRandomDraw(i);  // 0 <= r <= i
		}

		if (r!=i)
		{
		    shuffle[i] = r;
		}
		shuffle[r] = i;
	    }

// Do multiply as before, but now with randomized F
	    for (i = 0; i < rows; i++)
	    {
		double accum = 0;
		for (j = 0; j < incols; j++)
		{
// this inner loop useful for test purposes
//		    accum += 1.0;
		    accum += (*mHx)(i,j) * (f[shuffle[i]] + xb[i]);
		}
		(*mFp)(i,jout) = accum;
	    }
	}
    }

    MathMatrixV.push_back (mFp);
    char namebuf[64];
    sprintf (namebuf,"%s%d",MathMatrix_PREFIX,MathMatrixV.size());
    RESULT_BUF (namebuf);
    return TCL_OK;
}


// PermuteY produces the shuffled matrix Y (aka Y*) for FPHI
// Inputs: srF XB nP
//      note: srF is unsquared F residuals: (Y - XB)
//            regular F is (Y-XB)^2
//      this allows us to shuffle the shuffle the srF part but leave XB alone

extern "C" int PermuteYCmd (ClientData clientData, Tcl_Interp *interp,
			 int argc, char* argv[])
{
    if (argc != 4)
    {
	RESULT_LIT ("Usage: permutey <XB> <srF> <nP>");
	return TCL_ERROR;
    }
    MatrixXd* mXB = MathMatrixGet (argv[1], "permuteY", interp);
    if (!mXB) return TCL_ERROR;
    MatrixXd* mSRF = MathMatrixGet (argv[2], "permuteY", interp);
    if (!mSRF) return TCL_ERROR;

// Get nP

    errno = 0;
    char* eptr;
    long np = strtol (argv[3], &eptr, 10);
    if (errno || *eptr != '\0' || np < INT_MIN || np > INT_MAX)
    {
	RESULT_LIT ("permuteY: invalid nP");
	return TCL_ERROR;
    }
    int nP = (int) np;

// Check that F and XB are column vectors
    if (mSRF->cols() != 1)
    {
	RESULT_LIT ("permuteY: srF is not a column vector");
	return TCL_ERROR;
    }
    if (mXB->cols() != 1)
    {
	RESULT_LIT ("permuteY: XB is not a column vector");
    }
    
    int rows = mXB->rows();
    int cols = nP + 1;

    MatrixXd* mY = new MatrixXd(rows,cols);

// First column of output matrix mY is srF + XB

    int i,j;

    mY->col(0) = mSRF->col(0) + mXB->col(0);

// Remaining columns of output matrix mSRFp have shuffled F+XB
// We first generate shuffle index, then do multiply as above

// Seed

    AdvancedRandomSeed (1); // use same fixed seed each time

    int shuffle[rows];
    for (int jout = 1; jout < cols; jout++)
    {

// create shuffle index using
// Knuth aka Fisher-Yates Shuffle
//   Inside-Out version from Wikipedia 6/2015

	int r;
	for (i = 0; i < rows; i++)
	{
	    if (i == 0)
	    {
		r = 0;
	    }
	    else
	    {
		r = AdvancedRandomDraw(i);  // 0 <= r <= i
	    }
	    
	    if (r!=i)
	    {
		shuffle[i] = r;
	    }
	    shuffle[r] = i;
	}

// Create this column following shuffle index for mSRF part

	for (i = 0; i < rows; i++)
	{
	    (*mY)(i,jout) = (*mXB)(i,0) + (*mSRF)(shuffle[i],0);
	}
    }

    MathMatrixV.push_back (mY);
    char namebuf[64];
    sprintf (namebuf,"%s%d",MathMatrix_PREFIX,MathMatrixV.size());
    RESULT_BUF (namebuf);
    return TCL_OK;
}

extern "C" int DinverseCmd (ClientData clientData, Tcl_Interp *interp,
			 int argc, char* argv[])
{
    if (argc != 2)
    {
	RESULT_LIT ("Usage: Dinverse <Diagonal-MathMatrix>");
	return TCL_ERROR;
    }
    MatrixXd* m1 = MathMatrixGet (argv[1], "Inverse", interp);
    if (!m1) return TCL_ERROR;

    int rows = m1->rows();
    int cols = m1->cols();
    if (rows != cols)
    {
	RESULT_LIT ("Matrix inversion requires square matrix");
	return TCL_ERROR;
    }
    MathMatrixV.push_back (new MatrixXd(rows,cols));
    int mindex = MathMatrixV.size() - 1;
    MatrixXd* mout = MathMatrixV[mindex];
    *mout  << MatrixXd::Zero(rows,cols);

    for (int i = 0; i < rows; i++)
    {
	double coeff = (*m1)(i,i);
	if (coeff == 0)
	{
	    RESULT_LIT ("Cannot inverse diagonal matrix with zero on diagonal");
	    return TCL_ERROR;
	}
	double recip = 1.0 / coeff;
	(*mout)(i,i) = recip;
    }
    char namebuf[64];
    sprintf (namebuf,"%s%d",MathMatrix_PREFIX,MathMatrixV.size());
    RESULT_BUF (namebuf);
    return TCL_OK;
}    


extern "C" int InverseCmd (ClientData clientData, Tcl_Interp *interp,
			 int argc, char* argv[])
{
    if (argc < 2)
    {
	RESULT_LIT ("Usage: Inverse <MathMatrix>");
	return TCL_ERROR;
    }
    MatrixXd* m1 = MathMatrixGet (argv[1], "Inverse", interp);
    if (!m1) return TCL_ERROR;

    if (m1->rows() != m1->cols())
    {
	RESULT_LIT ("Matrix inversion requires square matrix");
	return TCL_ERROR;
    }
    MathMatrixV.push_back (new MatrixXd (m1->rows(), m1->cols()));
    int mindex = MathMatrixV.size() - 1;
    *MathMatrixV[mindex] = m1->inverse();
    char namebuf[64];
    sprintf (namebuf,"%s%d",MathMatrix_PREFIX,MathMatrixV.size());
    RESULT_BUF (namebuf);
    return TCL_OK;
}

extern "C" int RowCmd (ClientData clientData, Tcl_Interp *interp,
			 int argc, char* argv[])
{
    int mpos = 1;
    if (argc != 3)
    {
	RESULT_LIT ("Usage: row <matrix> <row>");
	return TCL_ERROR;
    }
    MatrixXd* umatrix = MathMatrixGet (argv[mpos], "row", interp);
    if (!umatrix) return TCL_ERROR;
    char* eptr;
    errno = 0;
    long lrow = strtol (argv[3-mpos], &eptr, 10);
    if (errno || *eptr != '\0' || lrow < 1 || lrow > umatrix->rows())
    {
	std::string errmess("Invalid matrix row: ");
	errmess += argv[3-mpos];
	RESULT_BUF(errmess.c_str());
	return TCL_ERROR;
    }
    int row = (int) lrow - 1;
    MatrixXd* vec = new MatrixXd (1,umatrix->cols());
    *vec = umatrix->row(row);
    MathMatrixV.push_back (vec);
    char namebuf[64];
    sprintf (namebuf, "%s%d", MathMatrix_PREFIX, MathMatrixV.size());
    RESULT_BUF (namebuf);
    return TCL_OK;
}

extern "C" int ColCmd (ClientData clientData, Tcl_Interp *interp,
			 int argc, char* argv[])
{
    int mpos = 1;
    if (argc != 3)
    {
	RESULT_LIT ("Usage: col <matrix> <col>");
	return TCL_ERROR;
    }
    MatrixXd* umatrix = MathMatrixGet (argv[mpos], "col", interp);
    if (!umatrix) return TCL_ERROR;
    char* eptr;
    errno = 0;
    long lcol = strtol (argv[3-mpos], &eptr, 10);
    if (errno || *eptr != '\0' || lcol < 1 || lcol > umatrix->cols())
    {
	std::string errmess("Invalid matrix col: ");
	errmess += argv[3-mpos];
	RESULT_BUF(errmess.c_str());
	return TCL_ERROR;
    }
    int col = (int) lcol-1;
    MatrixXd* vec = new MatrixXd (1,umatrix->cols());
    *vec = umatrix->col(col);
    MathMatrixV.push_back (vec);
    char namebuf[64];
    sprintf (namebuf, "%s%d", MathMatrix_PREFIX, MathMatrixV.size());
    RESULT_BUF (namebuf);
    return TCL_OK;
}

extern "C" int IdentityCmd (ClientData clientData, Tcl_Interp *interp,
			 int argc, char* argv[])
{
    if (argc != 2)
    {
	RESULT_LIT ("Usage: identity <rows>");
	return TCL_ERROR;
    }
    char* eptr;
    errno = 0;
    long lrows = strtol (argv[1], &eptr, 10);
    if (errno || *eptr != '\0' || lrows < 0)
    {
	std::string errmess("Invalid identity matrix size: ");
	errmess += argv[2];
	RESULT_BUF(errmess.c_str());
	return TCL_ERROR;
    }
    int rows = (int) lrows;
    MatrixXd* imatrix = new MatrixXd (rows, rows);
    *imatrix = MatrixXd::Identity (rows, rows);
    MathMatrixV.push_back(imatrix);
    char namebuf[64];
    sprintf (namebuf, "%s%d", MathMatrix_PREFIX, MathMatrixV.size());
    RESULT_BUF (namebuf);
    return TCL_OK;
}


extern "C" int DiagonalCmd (ClientData clientData, Tcl_Interp *interp,
			 int argc, char* argv[])
{
    int mpos = 1;
    int diagonal = 0;
    if (argc < 2 || argc > 3)
    {
	RESULT_LIT ("Usage: diagonal <matrix> [<super-sub-diagonal-number>]");
	return TCL_ERROR;
    }
    MatrixXd* umatrix = MathMatrixGet (argv[mpos], "diagonal", interp);
    if (!umatrix) return TCL_ERROR;
    int rows = umatrix->rows();
    int cols = umatrix->cols();
    if (rows == 1 || cols == 1)
    {
	if (argc > 2)
	{
	   RESULT_LIT ("sub/super not allowed if creating matrix from vector");
	   return TCL_ERROR;
	}
	int size = (rows > cols) ? rows : cols;
	int cdim = (rows == 1) ? 1 : 2;
	MatrixXd *dmatrix = new MatrixXd (size,size);
	MathMatrixV.push_back(dmatrix);
	for (int i = 0; i < size; i++)
	{
	    for (int j = 0; j < size; j++)
	    {
		if (i==j)
		{
		    if (cdim==1)
		    {
			(*dmatrix)(i,i) = (*umatrix)(0,i);
		    }
		    else
		    {
			(*dmatrix)(i,i) = (*umatrix)(i,0);
		    }
		}
		else
		{
		    (*dmatrix)(i,j) = 0;
		}
	    }
	}
    }
    else
    {
	char* eptr;
	errno = 0;
	if (argc == 3)
	{
	   long ldiagonal = strtol (argv[3-mpos], &eptr, 10);
	   if (errno || *eptr != '\0' || ldiagonal < - ((umatrix->rows())-1) ||
		diagonal > (umatrix->cols()-1))
	   {
	       std::string errmess("Invalid matrix super or sub diagonal: ");
	       errmess += argv[3-mpos];
	       RESULT_BUF(errmess.c_str());
	       return TCL_ERROR;
	   }
	   diagonal = (int) ldiagonal;
	}
	MatrixXd* vec = new MatrixXd (1,umatrix->diagonal(diagonal).size());
	*vec = umatrix->diagonal(diagonal);
	MathMatrixV.push_back (vec);
    }

    char namebuf[64];
    sprintf (namebuf, "%s%d", MathMatrix_PREFIX, MathMatrixV.size());
    RESULT_BUF (namebuf);
    return TCL_OK;
}

extern "C" int MathMatrixOutputCmd (ClientData clientData, Tcl_Interp *interp,
			 int argc, char* argv[])
{
    if (argc != 3)
    {
	RESULT_LIT ("Usage: output <matrix> <filename>");
	return TCL_ERROR;
    }
    MatrixXd* omatrix = MathMatrixGet (argv[1], "output", interp);
    if (!omatrix) return TCL_ERROR;
    char* filename = argv[2];

    ofstream ofile;
    ofile.open (filename);
    if (!ofile.is_open())
    {
	std::string errmess("Unable to open output file: ");
	errmess += filename;
	RESULT_BUF (errmess.c_str());
	return TCL_ERROR;
    }
    int rows = omatrix->rows();
    int cols = omatrix->cols();

    for (int i = 0; i < rows; i++)
    {
	for (int j = 0; j < cols;)
	{
	    ofile << (*omatrix)(i,j);
	    if (++j < cols) ofile << ",";
	}
	ofile << "\n";
    }
    ofile.close();
    std::stringstream rmess;
    rmess << rows << " rows, " << cols << " cols written";
    RESULT_BUF (rmess.str().c_str());
    return TCL_OK;
}

extern "C" int MathMatrixShowCmd (ClientData clientData, Tcl_Interp *interp,
			 int argc, char* argv[])
{
    char* mmid;
    char* pstring = 0;
    char* eptr;
    int row = INT_MIN;  // INT_MIN if left blank by user
    int col = INT_MIN;

    if (MathMatrixDebug) printf ("Entering MathMatrixShowCmd\n");

    char showformat[64];
    sprintf (showformat,"%%-.%ig",Option::get_int("MatrixNumberFormat"));
    if (MathMatrixDebug)printf ("showformat is %s\n",showformat);

// Identify matrix to be printed, or non-matrix

    if (argc != 2 &&  argc != 4)
    {
	RESULT_LIT ("Usage: show <matrix> [<row> <col>]");
	return TCL_ERROR;
    }
    mmid = argv[1];
    MatrixXd* showmatrix = MathMatrixGet (mmid, "show", interp);
    if (!showmatrix) return TCL_ERROR;

    if (argc == 4)
    {
	errno = 0;
	long lrow = strtol (argv[2], &eptr, 10);
	if (errno || *eptr != 0 ||
	    lrow > showmatrix->rows() || lrow < 1)
	{
	    RESULT_LIT ("Row invalid");
	    return TCL_ERROR;
	}
	row = lrow;

	errno = 0;
	long lcol = strtol (argv[3], &eptr, 10);
	if (errno || *eptr != 0 ||
	    lcol > showmatrix->cols() || lcol < 1)
	{
	    RESULT_LIT ("Col invalid");
	    return TCL_ERROR;
	}
	col = lcol;
    }
//
// First, row and col returns scalar
//
    if (argc==4)
    {
	char buf[64];
	sprintf (buf, showformat, (*showmatrix)(row-1,col-1));
	RESULT_BUF(buf);
	return TCL_OK;
    }
//
// show matrix
// pretty print entire matrix
//
// build format string 

    int rows = showmatrix->rows();
    int cols = showmatrix->cols();

    pstring = Strdup ("{"); // Begin matrix

    for (int i = 0; i < rows; i++)
    {
	if (i > 0) string__append (&pstring, "\n ");
	
	string__append (&pstring,"{");  // Begin row
	
	for (int j = 0; j < cols; j++)
	{
	    char dbuf[64];
	    sprintf (dbuf, showformat, (*showmatrix)(i,j));
	    string__append (&pstring, dbuf);
	    if (j < cols-1)
	    {
		string__append (&pstring, " ");
	    }
	}
	string__append (&pstring, "}"); // End row
    }
    string__append (&pstring, "}");
    RESULT_BUF (pstring);
    free (pstring);
    return TCL_OK;
}


int MathMatrixLoad (int argc, char* argv[], Tcl_Interp* interp)
{
    if (MathMatrixDebug) printf ("entering MathMatrixLoad\n");
    bool skipheader = true;
    int firstindex = 2;
    int urows = 0;
    bool no_file_needed = false;
//
// Load matrix may select columns by name, number, or just be a NumberSequence
//   or any combination.  Names which look like numbers must be enclosed in
//     apostrophes.
//
    SelectedColumns.clear();

    while (argc-firstindex>0 && argv[firstindex][0] == '-')
    {
	if (!Strcmp (argv[firstindex], "-noheader"))
	{
	    skipheader = false;
	    firstindex++;
	}
	else if (!Strcmp (argv[firstindex], "-rows"))
	{
	    firstindex++;
	    char* eptr;
	    errno = 0;
	    size_t lurows = strtol (argv[firstindex], &eptr, 10);
	    if (errno || *eptr != '\0')
	    {
		std::string errmess("Invalid -rows specification: ");
		errmess += argv[firstindex];
		RESULT_BUF (errmess.c_str());
		return TCL_ERROR;
	    }
	    urows = lurows;
	    firstindex++;
	}
	else if (!Strcmp (argv[firstindex], "-cols"))
	{
	    no_file_needed = true;
	    firstindex++;
	    if (argc-firstindex<=0)
	    {
		RESULT_LIT ("Usage: matrix load -cols {name1 name2} <filename>");
		return TCL_ERROR;
	    }
	    std::stringstream argstream(argv[firstindex]);
	    firstindex++;
	    std::string cell;

	    while (std::getline (argstream, cell, ' '))
	    {
		if (MathMatrixDebug) printf ("Parsing -col cell %s\n", cell.c_str());
		errno = 0;
		char* eptr;
		const char* cellstr = cell.c_str();
		size_t ssize = 0;

		if (cellstr[0] == '{') // this is a sequence
		{

// Handle cells as sequence until terminated

		    NumberSequence* sc = new NumberSequence();
		    if (cellstr[1] != 0)
		    {
			errno=0;
			eptr=0;
			double value = strtod (&cellstr[1], &eptr);
			if (errno || (*eptr != 0 && *eptr!='}'))
			{
			    std::string errmess = "Invalid column specifier: "
				+ cell;
			    return TCL_ERROR;
			}
			sc->push_back (value);
		    }
		    while (*eptr != '}' && 
			   std::getline (argstream, cell, ' '))
		    {
			errno=0;
			eptr=0;
			double value = strtod (cell.c_str(), &eptr);
			if (errno || (*eptr != 0 && *eptr != '}'))
			{
			    std::string errmess = "Invalid column specifier: "
				+ cell;
			    return TCL_ERROR;	
			}
			sc->push_back (value);
		    }
		    if (*eptr=='}' && eptr[1] == '\0')
		    {
			SelectedColumn* scol = new SelectedColumn(sc);
			SelectedColumns.push_back (scol);
		    }
		    else
		    {
			delete sc;
			RESULT_LIT ("Invalid column specification");
			return TCL_ERROR;
		    }
		    continue;
		}

		no_file_needed = false;  // more than just sequences

// See if entire cell parses as an integer, then use this as number

		bool use_value = true;
		long lvalue = strtol (cellstr, &eptr, 10);
		if (errno || *eptr!=0 || lvalue > INT_MAX ||
		    lvalue < INT_MIN)
		{
		    use_value = false;
		}
		int ivalue = (int) lvalue;
		if (use_value)
		{
		    SelectedColumn* sc = new SelectedColumn(ivalue);
		    SelectedColumns.push_back (sc);
		    if (MathMatrixDebug) 
		    {
			printf ("adding %d to column number list\n",ivalue);
			printf ("sc[0]->value is %d\n",
			    SelectedColumns[0]->number);
		    }
		}
		else
		{

// Since didn't parse as sequence or number, use as name

		    if (cell[0] == '\'')
		    {

// This is name in apostrophes

			std::string cell2;
			int i;
			for (i = 1; cell[i] && cell[i]!='\''; i++)
			{
			    cell2.push_back(cell[i]);
			}
			if (MathMatrixDebug)
			{
			    printf ("adding %s to column name list\n",
				    cell2.c_str());
			}
			SelectedColumn* sc = new SelectedColumn 
			    (cell2.c_str());
			SelectedColumns.push_back (sc);
		    }
		    else // name not enclosed in apostrophes
		    {
			if (MathMatrixDebug)
			{
			    printf ("adding %s to column name list\n",
				    cell.c_str());
			}
			SelectedColumn* sc = new SelectedColumn (cell.c_str());
			SelectedColumns.push_back (sc);
		    }
		}
	    }
	}
	else
	{
	    std::string errmess ("Matrix: Invalid argument: ");
	    errmess.append (argv[firstindex]);
	    RESULT_BUF (errmess.c_str());
	    return TCL_ERROR;
	}
    }

    if (argc <= firstindex)
    {
	RESULT_LIT (
	    "Usage: matrix load [<option>]+ <filename> for mathmatrix");
	return TCL_ERROR;
    }

// Determine required matrix rows (resize for each line is very expensive)
//   by finding number of lines in file

    if (MathMatrixDebug)printf("argv[%d] is %s\n",firstindex,argv[firstindex]);

    std::ifstream matrix_file (argv[firstindex]); // open csv file

    if (MathMatrixDebug) printf ("file opened\n");

    if (matrix_file.fail())
    {
	std::string NoSuch("Matrix: No such file: ");
	NoSuch.append(argv[firstindex]);
	RESULT_BUF (NoSuch.c_str());
	return TCL_ERROR;
    }
    if (skipheader)
    {
	char* suffix;
	if ( (suffix = strstr (argv[firstindex],".mat.csv")) &&
	     suffix[8] == '\0')
	{
	    if (MathMatrixDebug) 
	      {printf (".mat.csv defaults to -noheader option\n");}
	    skipheader = false;
	}
	else
	{
	    if (MathMatrixDebug)
	      {printf ("First file line will NOT be included in matrix\n");}
	}
    }
    if (!skipheader && MathMatrixDebug)
    {
	printf ("First file line will be included in matrix\n");
    }

    int rows = std::count(std::istreambuf_iterator<char>(matrix_file),
			     std::istreambuf_iterator<char>(), '\n');

// now get number of columns in first line

    matrix_file.clear(); // clear EOF
    matrix_file.seekg (0); // back to beginning
    std::string line;
    getline (matrix_file, line);

// we need to reject trailing null column so can't do this
//    int file_cols = std::count(line.begin(),line.end(),',');
//    file_cols++; // one more field than commas
//
// Instead, count the non-null fields by parsing
//   (which we might have to do anyway)

    std::vector<std::string> all_names;
    std::stringstream cstream(line);
    std::string ccell;
    while (std::getline (cstream, ccell, ',') )
    {
	all_names.push_back (ccell);
    }
    int file_cols = all_names.size();

// Validate selected columns

    for (int i = 0; i < SelectedColumns.size(); i++)
    {
	if (MathMatrixDebug) printf ("Validating selected column[%d]\n",i);

	if (SelectedColumns[i]->number != 0)
	{
	    if (MathMatrixDebug) printf ("  is number\n",i);
	    if (SelectedColumns[i]->number > file_cols)
	    {
		char buf[128];
		sprintf (buf, "Selected column %d not in file\n",
			 SelectedColumns[i]->number );
		RESULT_BUF (buf);
		return TCL_ERROR;
	    }
	}
	else if (SelectedColumns[i]->name)
	{
	    if (MathMatrixDebug) printf ("  is name\n");
	    char* name = SelectedColumns[i]->name;
	    if (MathMatrixDebug) printf ("    name: %s\n",name);
	    bool found = false;
	    for (int j = 0; j < all_names.size(); j++)
	    {
		if (!Strcmp(name, all_names[j].c_str()))
		{
		    if (MathMatrixDebug) printf ("  name pos %d\n", j);
		    
		    SelectedColumns[i]->number = j+1;  // 1 based here
		    found = true;
		    break;
		}
	    }
	    if (!found)
	    {
		std::string Name(name);
		std::string errmess = "No such field in csv file: " + Name;
		RESULT_BUF (errmess.c_str());
		return TCL_ERROR;
	    }
	}
    }

// Skip first line if not part of matrix

    if (skipheader) {
	getline (matrix_file, line);
	if (MathMatrixDebug) printf ("Now Skipping to line past Header\n");
    } else {
	if (MathMatrixDebug) printf ("  Continuing with first line\n");
    }
	

// Allocate matrix

    int cols = file_cols;
    if (SelectedColumns.size()) cols = SelectedColumns.size();
    if (!skipheader) rows++;
    if (urows) rows = urows;

    if (MathMatrixDebug) printf ("making new matrix (%d,%d)\n",rows,cols);

    MatrixXd* newmatrix = new MatrixXd (rows, cols);

// Fill in one matrix row from each line in file
//   (first line already read)

    int row = 0;
    while (1) // reading lines from from file
    {
	std::vector<std::string> data;
	data.clear();
	std::string cell;

	std::stringstream lineStream(line); // makes line into lineStream

// Parse each comma delimited cell into vector data

	while (std::getline (lineStream, cell, ',') )
	{
	    data.push_back (cell);
	}
	
	if (data.size() < cols)
	{
	    std::stringstream errmess;;
	    errmess << "Line " << row+1 << " has too few columns";
	    RESULT_BUF(errmess.str().c_str());
	    return TCL_ERROR;
	}

// For each column in matrix, load column, using selections if made

	for (int col = 0; col < cols; col++)
	{
	    int getcol = col; // getcol is zero based, or -1 meaning sequence
	    if (SelectedColumns.size())
	    {
		if (MathMatrixDebug)
		{
		    printf ("getting selected column at position 0[%d]\n",col);
		}
		if (SelectedColumns[col]->number)
		{
		    getcol = SelectedColumns[col]->number - 1;
		    if (MathMatrixDebug)
		    {
			printf ("getting column offset %d\n", getcol);
		    }
		}
		else if (SelectedColumns[col]->sequence)
		{
		    if (MathMatrixDebug)
		    {
			printf ("getting column from sequence\n");
		    }
		    getcol = -1;
		}
	    }
	    if (MathMatrixDebug)
	    {
		printf ("selected column %d\n",getcol);
	    }

// Would have done this nice new shiny version
//   in which the code is perfectly clear and the catch catches all errors
// but our 4.4.7 compiler doesn't support c++11 for std::stod

#if 0
	    try
	    {
		(*newmatrix)(row,col) =  std::stod (data[getcol]);
	    }
	    catch (...)
	    {
	    }
#endif

// So instead I have to use ugly old strtod and check for size and error
//   and check for Fortran exponent if first attempt fails

	    double value;
	    if (getcol == -1)
	    {
		value = SelectedColumns[col]->sequence->next();
	    }
	    else
	    {
		bool number_ok = true;
		errno = 0;
		char* endptr;
		const char* beginptr = data[getcol].c_str();

		value = strtod (beginptr, &endptr);

		size_t dlen = endptr - beginptr;
		bool tooshort = dlen != data[getcol].size();
		if (errno || tooshort)
		{
		    number_ok = false;
		    if (tooshort)  // Check and convert Fortran exponent
		    {
			char* number_str = Strdup (data[getcol].c_str());
			char* nptr = number_str;
			while (*nptr != 0)
			{
			    if (*nptr == 'D' || *nptr == 'd')
			    {
				*nptr = 'e';
				errno = 0;
				value = strtod (number_str, &endptr);
				dlen = endptr - number_str;
				tooshort = dlen != data[getcol].size();
				if (errno || tooshort)
				{
				    free (number_str);
				    number_str = 0;
				    break;
				}
				number_ok = true;
				break;
			    }
			    nptr++;
			}
			if (number_str) free (number_str);
		    }
		}
		if (!number_ok)
		{
		    std::stringstream errmess;
		    errmess << "Invalid number " <<data[getcol] <<
			" in csv matrix file line " << row+1;
		    RESULT_BUF (errmess.str().c_str());
		    return TCL_ERROR;
		}
	    }
	    if (MathMatrixDebug) 
	      {printf ("setting value(%d,%d) = %g\n", row, col, value);}

	    (*newmatrix)(row,col) = value;

	}
	row++;
	if (MathMatrixDebug) printf ("row incremented to %d\n",row);
	if (row>=rows) break;
	if (std::getline (matrix_file, line)) continue;
	break;
    }

    matrix_file.close();

// resize back to actual number of rows if needed

    if (row <= rows && !urows)
    {
	rows--;
	if (MathMatrixDebug) printf ("resizing matrix back to %d\n",rows);
	(*newmatrix).conservativeResize(rows,cols);
    }

// Now add this to matrix vector and return handle
    if (MathMatrixDebug) {
	printf ("Vsize before new is %d\n",MathMatrixV.size());
    }

    MathMatrixV.push_back(newmatrix);

    if (MathMatrixDebug) {
	printf ("Vsize before new is %d\n",MathMatrixV.size());
    }

    char namebuf[64];
    sprintf (namebuf, "%s%d", MathMatrix_PREFIX, MathMatrixV.size());
    RESULT_BUF (namebuf);
    return TCL_OK;
}



extern "C" int TransposeCmd (ClientData clientData, Tcl_Interp *interp,
			 int argc, char* argv[])
{
    if (argc != 2 || strncmp (MathMatrix_PREFIX,argv[1],
			      strlen(MathMatrix_PREFIX)))
    {
	RESULT_LIT ("transpose <MathMatrix>");
	return TCL_ERROR;
    }

// Identify matrix to be transposed

    MatrixXd* oldmatrix = MathMatrixGet (argv[1],"Transpose",interp);
    if (!oldmatrix) return TCL_ERROR;

// Transpose is done by making new matrix of the correct size
//   and assigning to it.

    int rows = oldmatrix->rows();
    int cols = oldmatrix->cols();
    if (MathMatrixDebug) printf ("making new matrix(%d,%d)\n",
				 cols,rows);
    MatrixXd* newmatrix = new MatrixXd (cols,rows);
    *newmatrix = oldmatrix->transpose();

    MathMatrixV.push_back(newmatrix);

    char namebuf[64];
    sprintf (namebuf,"%s%d",MathMatrix_PREFIX,MathMatrixV.size());
    RESULT_BUF (namebuf);
    return TCL_OK;
}

int MathMatrixNew (int argc, char* argv[], Tcl_Interp* interp)
{
    if (argc < 3)
    {
	RESULT_LIT ("matrix new {{1 2} {3 4}}");
	return TCL_ERROR;
    }

    char* inbuf = argv[2];
    char* inptr = inbuf;
    char* endptr;
    int rows = 0;
    int cols = 0;

// Scan list-style matrix to determine size
// Scan first row to get column size
	
    while (*inptr != '\0' && isspace(*inptr)) inptr++;
    if (*inptr != '{')
    {
	RESULT_LIT ("Math matrix missing leading {");
	return TCL_ERROR;
    }
    inptr++;
    while (1)
    {

// skip over leading space	    

	while (*inptr != '\0' && isspace(*inptr)) inptr++;
	
// count value token(s) if non-space exists within braces

	char* beginptr = inptr;
	while (*inptr != '\0' && *inptr != '}' && !isspace(*inptr)) inptr++;
	if (inptr > beginptr) cols++;
	if (*inptr == '\0')
	{
	    RESULT_LIT ("Math matrix list missing row terminator }");
	    return TCL_ERROR;
	}
	if (*inptr == '}') break;
    }

// exit loop with pointer at terminating brace of first row

    rows++;
    inptr++;

// Scan second and later rows

    while (1)
    {

// skip over leading space

	while (*inptr != '\0' && isspace(*inptr)) inptr++;
	if (*inptr == '\0') break;

// first non-space must be open brace

	if (*inptr != '{')
	{
	    RESULT_LIT ("Math Matrix list missing row initiator {");
	    return TCL_ERROR;
	}
	inptr++;
	char* beginrow = inptr;
	    
	while (*inptr != '\0' && *inptr != '}') inptr++;
	if (inptr > beginrow) rows++;
	if (*inptr == '\0') break;
	inptr++;  // skip over terminating brace for next loop
    }

// Allocate matrix

    if (rows == 0)
    {
	RESULT_LIT ("Math matrix list has zero rows");
	return TCL_ERROR;
    }
    if (cols == 0)
    {
	RESULT_LIT ("Math matrix list has zero columns");
	return TCL_ERROR;
    }
    if (MathMatrixDebug)
	printf ("Allocating %d rows and %d columns\n", rows, cols);
    
    MathMatrixV.push_back(new MatrixXd (rows,cols));
    
    int mindex = MathMatrixV.size()-1;
    if (MathMatrixDebug)
	printf ("new matrix has %d rows and %d columns\n",
		MathMatrixV[mindex]->rows(),
		MathMatrixV[mindex]->cols());
    
// Re-read lists to store values in matrix
    
    inptr = inbuf;
    int row = 0;
    int col = 0;
    
    while (1) // All Rows
    {
	
// Skip over leading space, quitting at EOL
	
	while (*inptr != '\0' && isspace(*inptr)) inptr++;
	if (*inptr == '\0') break;
	
// First character must be open brace
	
	if (*inptr != '{')
	{
	    RESULT_LIT ("Math matrix list missing leading {");
	    return TCL_ERROR;
	}
	inptr++;
	
// Parse value(s)
	
	while (1) // All Columns in Row
	{
	    errno = 0;
	    
// Parse double value (strtod skips leading space too)
	    
	    double value = strtod (inptr, &endptr);
	    if (errno || inptr == endptr)
	    {
		RESULT_LIT ("Invalid value in math matrix");
		return TCL_ERROR;
	    }

// Assign value

	    if (MathMatrixDebug) printf 
			     ("setting matrix(%d,%d) = %g\n",row, col, value);
	    (*MathMatrixV[mindex])(row,col) = value;
	    col++;
	    inptr = endptr;

// Skip over trailing space

	    while (*inptr != '\0' && isspace (*inptr)) inptr++;
	    
// If close brace, reset column and loop		
	    
	    if (*inptr == '}')
	    {
		inptr++;
		col = 0;
		row++;
		break;
	    }
	    
// If no close brace, bad format
	    
	    if (*inptr == '\0')
	    {
		RESULT_LIT ("Math matrix missing row terminator }");
		return TCL_ERROR;
	    }
	}
    }
    char namebuf[64];
    sprintf (namebuf,"%s%d",MathMatrix_PREFIX,MathMatrixV.size());
    RESULT_BUF (namebuf);
    return TCL_OK;
}

// MathMatrixCmd handles some short MathMatrix commands and dipatches others

extern "C" int MathMatrixCmd (ClientData clientData, Tcl_Interp *interp,
			 int argc, char* argv[])
{
    if (argc>1 && !Strcmp (argv[1], "debug"))
    {
	if (argc < 3)
	{
	    RESULT_LIT ("mathmatrix debug [<on>|<off>");
	    return TCL_ERROR;
	}
	if (!Strcmp (argv[2],"on"))
	{
	    MathMatrixDebug = 1;
	} else if (!Strcmp (argv[2],"off"))
	{
	    MathMatrixDebug = 0;
	}
	else
	{
	    RESULT_LIT ("mathmatrix debug [<on>|<off>");
	    return TCL_ERROR;
	}
	return TCL_OK;
    }

    if (MathMatrixDebug) printf ("MathMatrixCmd argc is %d,argv[0],argv[1] is %s,%s\n",argc, argv[0], argv[1]);

    if (argc>1 && !Strcmp (argv[1], "delete"))
    {
	for (int i = 2; i < argc; i++)
	{
	    MatrixXd* m1 = MathMatrixGet (argv[i], "Delete", interp);
	    if (!m1) return TCL_ERROR;
	    delete m1;
	    MathMatrixV[LastIndex] = 0;
	}
	return TCL_OK;
    }

    if (argc>2 && !Strcmp (argv[2], "delete"))
    {
	MatrixXd* m1 = MathMatrixGet (argv[1], "Delete", interp);
	if (!m1) return TCL_ERROR;
	delete m1;
	MathMatrixV[LastIndex] = 0;
	return TCL_OK;
    }

    if (argc==2 && !Strcmp (argv[1], "reset"))
    {
	for (int i = 0; i < MathMatrixV.size(); i++)
	{
	    if (MathMatrixV[i])
	    {
		if (MathMatrixDebug) printf ("Deleting MathMatrixID[%d]\n",i);
		delete MathMatrixV[i];
	    }
	}
	MathMatrixV.clear();
	EvaluesIndex = -1;
	EvectorsIndex = -1;
	Min_i = INT_MIN;
	Min_j = INT_MIN;
	Max_i = INT_MIN;
	Max_j = INT_MIN;
	FortMatrix = 0;

	return TCL_OK;
    }

// "rows" is a command that dispatches here, so in argv[0]
    if (argc==2 && !Strcmp (argv[0], "rows"))
    {
	MatrixXd* m1 = MathMatrixGet (argv[1], "Matrix", interp);
	if (!m1) return TCL_ERROR;
	char buf[64];
	sprintf (buf, "%d", m1->rows() );
	RESULT_BUF (buf);
	return TCL_OK;
    }
    if (argc==2 && !Strcmp (argv[0], "cols"))
    {
	MatrixXd* m1 = MathMatrixGet (argv[1], "Matrix", interp);
	if (!m1) return TCL_ERROR;
	char buf[64];
	sprintf (buf, "%d", m1->cols() );
	RESULT_BUF (buf);
	return TCL_OK;
    }
    if (argc>1 && !Strcmp (argv[1],"load"))
    {
	return MathMatrixLoad (argc, argv, interp);
    }

    if (argc>1 && !Strcmp (argv[1],"new"))
    {
	return MathMatrixNew (argc, argv, interp);
    }

    if (argc==2 && !Strcmp (argv[1],"lastid"))
    {
	int idnum = MathMatrixV.size();
	if (idnum == 0)
	{
	    RESULT_LIT ("No matrix has been created");
	    return TCL_ERROR;
	}
	char namebuf[64];
	sprintf (namebuf,"%s%d",MathMatrix_PREFIX,idnum);
	RESULT_BUF (namebuf);
	return TCL_OK;
    }
	
    RESULT_LIT ("Invalid MathMatrix command");
    return TCL_ERROR;
}




// Fortran interfaces

extern "C" void fmatrix_start_ (int* rowsp, int* colsp, int* idreturn)
{
//    printf ("allocating matrix at %d %d\n",*rowsp,*colsp);
    MatrixXd* fortmatrix = new MatrixXd (*rowsp, *colsp);
    MathMatrixV.push_back(fortmatrix);
    *idreturn = MathMatrixV.size();
}


extern "C" void fmatrix_morerows_ (int* id, int* rowsp)
{
    int mindex = (*id)-1;
    if (mindex < 0 || mindex >= MathMatrixV.size() || MathMatrixV[mindex] == 0)
    {
	printf ("Internal Fortran Matrix not allocated\n");
    }
    else
    {
	printf ("Resizing matrix to %d rows",*rowsp);
	MathMatrixV[mindex]->conservativeResize(*rowsp,NoChange);
    }
}


extern "C" void fmatrix_insert_ (int* id, int* rowsp, int* colsp, double* value)
{
    int mindex = (*id)-1;
    if (mindex < 0 || mindex >= MathMatrixV.size() || MathMatrixV[mindex] == 0)
    {
	printf ("Internal Fortran Matrix not allocated\n");
    }
    else
    {
	if (*value>0.1 && *id==1)
	    printf ("storing at %d %d = %g\n",(*rowsp)-1,(*colsp)-1,*value);
	(*MathMatrixV[mindex])((*rowsp)-1,(*colsp)-1) = *value;
    }
}


    



	
	


