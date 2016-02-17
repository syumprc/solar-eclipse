/* =============================================================================
file:		utils.h
description:	header file for utils.h which contains general utility functions.
author:		Harald Göring
============================================================================= */



#ifndef _utils_h
#define _utils_h



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <time.h>
#include <string.h>
#include <float.h>
#include <ctype.h>



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
global preprocessor constants
*/
#define	ERRORFATAL	-99	/*return code of program
				  if execution is terminated by a function
				  upon fatal error*/



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
declarations of external variables.
corresponding definitions must be in other files for use.
*/



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
declarations of functions in this file
*/
/*mathematical functions*/
double	ln			(double);

/*file functions*/
FILE*	fopen2			(const char*, const char*);

/*memory functions*/
void*	malloc2			(int);
double**mallocMatrixDouble	(int, int);
void	outOfMemory		(void);

/*misc. functions*/
void	outEgo			(FILE*, char*, char*, char*, char*, char*);
int	countSpaceDelimitedObjectsOnLine		(FILE*);
int	countSingleCharacterDelimitedObjectsOnLine	(FILE*, char);
int	countLinesInFile	(FILE*);

/*statistical functions*/
double	normalCdf		(double, double, double);
double	chiCdf			(double, int);

/*genetics functions*/
double	chisq2lod		(double);



#endif
