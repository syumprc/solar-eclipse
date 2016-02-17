/* =============================================================================
file:		utils.h
description:	file with general utility functions
author:		Harald Göring
============================================================================= */



#include "homo.h"



/* =============================================================================
function:	ln()
version:	1.0
date:		20000308
author:		Harald Göring
============================================================================= */
/* -----------------------------------------------------------------------------
computes ln(x) if x > 0. 
otherwise, prints error message to stderr & teminates program.
*/
double		/*ln(x)*/
ln
(
double x	/*number of which ln is to be taken*/
)
{
	if (x <= 0.0) {
		fprintf(stderr, 
			"\a\n"
			"ERROR in function ln():\n"
			"Attempting to take natural logarithm of non-positive number: %f\n" 
			"Terminating.\n\n",
			x);
		exit(ERRORFATAL);
	};

	return (log(x));
}



/* =============================================================================
function:	fopen2()
version:	1.1
date:		20000223
author:		Harald Göring
============================================================================= */
/* -----------------------------------------------------------------------------
opens a file.
otherwise, prints error message to stderr & teminates program.
*/
FILE*			/*opened file*/
fopen2
(
const char	*fn,	/*file name*/
const char	*mode	/*mode specifying how to open file*/
)
{
	FILE	*fp;

	fp = fopen(fn, mode);
	if (fp == NULL) {
		fprintf(stderr, 
			"\a\n"
			"ERROR in function fopen2():\n"
			"Failure opening file: %s\n"
			"Terminating.\n\n", 
			fn);
		exit(ERRORFATAL);
	};

	return fp;
}



/* =============================================================================
function:	malloc2()
version	:	1.0
date:		20000223
author:		Harald Göring
============================================================================= */
/* -----------------------------------------------------------------------------
allocates space.
if unsuccessful, calls outOfMemory().
calls:		outOfMemory()
*/
void*		/*allocated memory*/
malloc2
(
int	size	/*size of memory to be allocated*/
)
{
	void	*ptr;

	if ((ptr = (void*)malloc(size)) == NULL) {
		outOfMemory();
	};

	return ptr;
}



/* =============================================================================
function:	mallocMatrixDouble()
version:	1.0
date:		20000223
author:		Harald Göring
============================================================================= */
/* -----------------------------------------------------------------------------
allocates space for a matrix of type double **.
calls:		malloc2()
*/
double**	/*allocated matrix*/
mallocMatrixDouble
(
int	nrows,	/*# of rows*/
int	ncols	/*# of columns*/
)
{
	double	**matrix;	/*matrix*/
	int	irow;		/*row index*/
	
	/*allocate space for rows*/
	matrix = (double**)malloc2(nrows*sizeof(double*));

	/*allocate space for columns*/
	for (irow=0; irow < nrows; irow++) {
		matrix[irow] = (double*)malloc2(ncols*sizeof(double));
	};

	return matrix;
}



/* =============================================================================
function:	outOfMemory()
version:	1.0
date:		20000224
author:		Harald Göring
============================================================================= */
/* -----------------------------------------------------------------------------
prints error message to stderr & terminates program.
called if sufficient memory cannot be allocated.
*/
void
outOfMemory
(
void
)
{
	fprintf(stderr, 
		"\a\n"
		"ERROR:\n"
		"Unable to allocate sufficient memory.\n"
		"Terminating.\n\n");
	exit(ERRORFATAL);

	return;
}



/* =============================================================================
function:       outEgo()
version:        1.0
date:           20001014
author:         Harald Göring
============================================================================= */
/* -----------------------------------------------------------------------------
prints ego stuff to open file.
*/
void
outEgo
(
FILE	*fp,		/*output file*/
char	*border,
char	*program,
char	*author,
char	*readme,
char	*email
)
{
	time_t	currTime;	/*current calendar time*/
	char	*currLocTime;	/*current local time*/

	time(&currTime);
	currLocTime = ctime(&currTime);

        fprintf(fp, "%s\n", border);
        fprintf(fp, "%s\n", program);
        fprintf(fp, "%s\n", author);
        fprintf(fp, "%s\n", readme);
        fprintf(fp, "%s\n", email);
        fprintf(fp, "%s\n", border);
        fprintf(fp, "%s\n", currLocTime);
        fprintf(fp, "\n");

        return;
}



/* =============================================================================
function:	countSpaceObjectsOnLine()
version:	1.0
date:		19970726
author:         Harald Göring
============================================================================= */
/* -----------------------------------------------------------------------------
counts no. of objects on line in a file 
(object: sequence of non-white chars; objects separated by white space).
assumes file is opened for reading & pointer is at beginning of line.
*/
int		/*no. of objects on line*/
countSpaceDelimitedObjectsOnLine
(
FILE	*fp	/*file pointer, positioned at beginning of line to be examined*/
)
{
	int	nobj;		/*no. of objects on line*/
	int	ch;		/*char*/
	int	iicurrChWhite, 	/*indicator: 1 if current  char is white char*/
		iiprevChWhite;	/*indicator: 1 if previous char is white char*/

	for (nobj=0, iiprevChWhite=1, ch=fgetc(fp); ch != '\n' && ch != EOF; ch=fgetc(fp)) {

		iicurrChWhite = isspace(ch);
		if (iiprevChWhite && !iicurrChWhite) {
			nobj++;
		};
		iiprevChWhite = iicurrChWhite;
	};

	return nobj;
}



/* =============================================================================
function:	countSingleCharacterDelimitedObjectsOnLine()
version:	1.0
date:		20020913
author:         Harald Göring
============================================================================= */
/* -----------------------------------------------------------------------------
counts no. of objects on line in a file delimited by a single character.
assumes file is opened for reading & pointer is at beginning of line.
*/
int			/*no. of objects on line*/
countSingleCharacterDelimitedObjectsOnLine
(
FILE	*fp,		/*file pointer, positioned at beginning of line to be examined*/
char	delimitor	/*character delimiting line*/
)
{
	int	nobj;		/*no. of objects on line*/
	int	ch;		/*char*/

	for (nobj=1, ch=fgetc(fp); ch != '\n' && ch != EOF; ch=fgetc(fp)) {
		if (ch == delimitor) {
			nobj++;
		};
	};

	return nobj;
}



/* =============================================================================
function:	countLinesInFile()
version:	1.0
date:		20011124
author:         Harald Göring
============================================================================= */
/* -----------------------------------------------------------------------------
counts # of lines in a file from current position of file pointer.
*/
int		/*# of lines in file*/
countLinesInFile
(
FILE	*fp	/*file pointer*/
)
{
	int	ch;		/*char. in file*/
	int	nlines;		/*# of lines in file*/

	nlines = 0;
	while ((ch = fgetc(fp)) != EOF) {
		if (ch == '\n') {
			nlines++;
		};
	};

	return nlines;
}



/* =============================================================================
function:	normalCdf()
version	:	1.0
date:		20000916
author:		Harald Göring
============================================================================= */
/* -----------------------------------------------------------------------------
computes cdf of N(mean,var) distribution for a given normal deviate.
after 26.2.16 in Abramowitz M, Stegun I.A. (eds.)
Handbook of mathematical functions. 9th ed.. Dover Publications, New York. 1964.
abs(error) < ? (look it up)
issues error message to stderr & terminates program after error.
*/
double		/*cdf(deviate)*/
normalCdf
(
double	m,	/*mean*/
double	v,	/*var*/
double	x	/*deviate*/
)
{
	double	abs_x;	/*absolute value of x*/
	double	prob;	/*tail prob. (lower or upper depending on sign of x)*/
	double	f;	/*f(y): density of N(0,1) at y*/
	double	t;	/*for computation only*/

	if (v <= 0.0) {
		fprintf(stderr, 
			"\a\n"
			"ERROR in function normalCdf():\n"
			"Non-positive variance: %f\n"
			"Terminating.\n\n", 
			v);
		exit(ERRORFATAL);
	};

	/*transform into N(0,1) deviate*/
	x = (x - m) / sqrt(v);

	/*take absolute value*/
	if (x < 0.0) {
		abs_x = -x;
	} else {
		abs_x = x;
	};

	if (abs_x > 12.0) {
		prob = 0.0;
	} else {
		f = 0.3989423 * exp(-0.5 * abs_x * abs_x);
		t = 1.0 / (1.0 + 0.33267 * abs_x);
		prob = 1.0 - f * t * (0.4361836 + t * (-0.1201676 + t * 0.9372980));
	};

	if (x >= 0.0) {
		return prob;
	} else {
		return (1.0 - prob);
	};
}



/* =============================================================================
function:	chiCdf()
version	:	1.0
date:		200020920
author:		Jurg Ott?
		(taken from his Linkage Utility Programs
		[Ott J (1985) Analysis of human genetic linkage. Johns Hopkins
		University Press])
		modified by Harald Göring
============================================================================= */
/* -----------------------------------------------------------------------------
computes cdf of chi-squared (ndf) distribution for a given chi-squared deviate.
after 26.4.4 and 26.2.21 in Abramowitz M, Stegun I.A. (eds.)
Handbook of mathematical functions. 9th ed.. Dover Publications, New York. 1964.
abs(error) < ? (look it up)
issues error message to stderr & terminates program after error.
calls:		UNIFORM()
*/
double		/*cdf(deviate)*/
chiCdf
(
double	x2,	/*deviate*/
int	ndf	/*no. of degrees of freedom*/
)
{
   double	z, p, x = 0.0, sum, re, ch, chp;
   long int	i, n1, n2;

	if (x2 < 0.0) {
		fprintf(stderr, 
			"\a\n"
			"ERROR in function chiCdf():\n"
			"Negative deviate: %f\n"
			"Terminating.\n\n", 
			x);
		exit(ERRORFATAL);
	};
	if (ndf <= 0) {
		fprintf(stderr, 
			"\a\n"
			"ERROR in chiCdf():\n"
			"Non-positive number of degrees of freedom: %d\n"
			"Terminating.\n\n", 
			ndf);
		exit(ERRORFATAL);
	};

	if (ndf == 1) {					/*ndf = 1*/
		x = sqrt(x2);
		return (1.0 - 2.0 * (1.0 - normalCdf(0.0, 1.0, x)));
	};

	if (x2 > 169.0) {  /*formula 26.4.14, p.941*/	/*x2 very large*/
		ch = 2.0 / (9.0 * ndf);
		x = (exp(ln(x2 / ndf) / 3.0) - 1.0 + ch) / sqrt(ch);
		return (normalCdf(0.0, 1.0, x));
	};

	if (ndf == 2) {					/*ndf = 2*/
		return (1.0 - exp(-0.5 * x2));
	};

	n1 = (ndf - 1) / 2;
	n2 = (ndf - 2) / 2;

   	if (n1 == n2) {					/*ndf = 4, 6, 8, ...*/
		sum = 1.0;
		re  = 0.5 * x2;
		ch  = 1.0;
		for (i = 1; i <= n2; i++) {
			ch = ch * re / i;
			sum += ch;
		};
		return (1.0 - exp(-re) * sum);
	};

	ch = sqrt(x2);
	z  = 0.39894228 * exp(-0.5 * x2);
	p  = 1.0 - normalCdf(0.0, 1.0, ch);

	if (ndf == 3) {					/*ndf = 3*/
		return (1.0 - 2.0 * (p + z * ch));
	};
	
							/*ndf = 5, 7, 9, ...*/
	chp = ch;
	re = 1.0;
	for (i = 2; i <= n1; i++) {
		re += 2.0;
		chp = chp * x2 / re;
		ch += chp;
	};
	return (1.0 - 2.0 * (p + z * ch));
}



/* =============================================================================
function:	chisq2lod()
version:	1.0
date:		20010108
author:		Harald Göring
============================================================================= */
/* -----------------------------------------------------------------------------
converts chi-squared statistic (0.5(0) + 0.5(1) df.) to lod score.
*/
double		/*lod score*/
chisq2lod
(
double	chisq	/*chi-squared statistic*/
)
{
	return (chisq / (2.0*ln(10.0)));
}
