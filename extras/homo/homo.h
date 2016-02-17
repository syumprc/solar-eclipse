/* =============================================================================
file:		homo.h
description:	header file for program homo
author:		Harald Göring
============================================================================= */



#ifndef homo_h
#define homo_h



#include "utils.h"



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
global preprocessor constants
*/
/*program message*/
#define PROGRAM		"program homo, test version 0.2, 22 October 2002"
#define AUTHOR		"by Harald Göring"
#define README		"See README file for documentation."
#define EMAIL		"For bug reports, comments or questions, send email to hgoring@darwin.sfbr.org."

/*output stuff*/
#define BORDERSINGLE	"--------------------------------------------------------------------------------"
#define BORDERDOUBLE	"================================================================================"

/*program return values*/
#define	NORMAL		0			/*program return value after   normal termination*/
#define	ABNORMAL	-1			/*program return value after abnormal termination*/

/*file names*/
#define	FNIN		"homo.in"		/*name of input  file*/
#define FNOUT		"homo.out"		/*name of output file*/

/*memory allocation*/
#define MAXCHLINE	2000			/*max. + 1 no. of characters in a line*/
#define MAXCHWORD	100			/*max. + 1 no. of characters in a word*/

/*misc.*/
#define	ALPHAINCREMENT	0.001			/*stepsize in homogeneity parameter over which likelihood is maximized*/



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
global derived data structures
*/
typedef struct {	/*info. about a pedigree*/
	/*as in input*/
	char		*id;		/*name*/
	double		lnL0,		/*     ln L under H0 (no linkage)*/
			lnL1,		/*max. ln L under H1 (linkage, homogeneity)*/
			lnL2;		/*max. ln L under H2 (linkage, heterogeneity)*/
	double		*lnL1grid;	/*ln L's under H1 for different values of linkage parameter*/
	double		pLinked;	/*probability that pedigree is of linked type*/
} PedS;



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
global variable declarations;
corresponding definitions are in file containing main()
*/
extern long int		seed;		/*seed for random no. generator*/



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
global function declarations
*/
/*in file containing main():*/
int	main				(int, char*[]);
void	inData				(char*, PedS*, int, double*, int, char);
double	computeLnL0			(PedS*, int);
double	computeLnL1			(PedS*, int, double*, int, double*);
double	computeLnL2			(PedS*, int, double*, int, double*, double*, double);
void	outData				(FILE*, PedS*, int, double, double, double, double, double, double, double);



#endif
