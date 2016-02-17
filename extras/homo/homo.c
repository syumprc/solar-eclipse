/* =============================================================================
file:		homo.c
description:	main() + computational functions of program homo
author:		Harald Göring
============================================================================= */



#include "homo.h"



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
global variable definitions
*/



/* -----------------------------------------------------------------------------
controls overall program flow. 
*/
int	/*NORMAL after normal termination, ABNORMAL otherwise*/
main
(
int	argc,
char	*argv[]
)
{
	char	delimChar;	/*delimiting character in input file*/
	FILE	*fpIn,		/*input  file*/
		*fpOut;		/*input  file*/
	PedS	*peds;		/*pedigrees*/
	int	npeds;		/*no. of pedigrees*/
	int	iped;		/*index for pedigrees*/
	double	*thetas;	/*       different values of linkage parameter given in input file*/
	int	nthetas;	/*no. of different values of linkage parameter given in input file*/
	double	thetaMLE_H1,	/*ML estimate of linkage parameter under H1 (linkage, homogeneity)*/	
		thetaMLE_H2;	/*ML estimate of linkage parameter under H2 (linkage, heterogeneity)*/	
	double	alphaMLE_H2;	/*ML estimate of homogeneity parameter under H2 (linkage, heterogeneity)
				  alpha=1: all pedigrees are of linked type*/
	double	lnL0,		/*ln L under H0 (no linkage)*/
		lnL1,		/*ln L under H1 (linkage, homogeneity)*/
		lnL2;		/*ln L under H0 (linkage, heterogeneity)*/

	/*deal with command line option*/
	if (!(argc == 1 || (argc == 2 && strcmp(argv[1], "-c") == 0))) {
		fprintf(stderr, "\a\n"
			"program usage:\n"
			"  homo    : space-delimited input file\n"
			"  homo -c : comma-delimited input file\n\n");
		exit(ABNORMAL);
	} else {
		if (argc == 1) {
			/*default: space-delimited*/
			delimChar = ' ';
		} else {
				 /*comma-delimited*/
			delimChar = ',';
		};
	};

	fpOut = fopen2(FNOUT, "w");

	outEgo(fpOut, BORDERDOUBLE, PROGRAM, AUTHOR, README, EMAIL);

	/*how many pedigrees?*/
	fpIn = fopen2(FNIN, "r");
	npeds = countLinesInFile(fpIn) - 1;
	fclose(fpIn);

	/*how many different values of linkage parameter in input file?*/
	fpIn = fopen2(FNIN, "r");
	if (delimChar == ',') {
		nthetas = countSingleCharacterDelimitedObjectsOnLine(fpIn, ',') - 2;
	} else {
		nthetas = countSpaceDelimitedObjectsOnLine(fpIn) - 2;
	};
	fclose(fpIn);

	/*allocate some space*/
	peds = (PedS*)malloc2(npeds * sizeof(PedS));
	for (iped = 0; iped < npeds; iped++) {
		peds[iped].id             = (char  *)malloc2(MAXCHWORD * sizeof(char  ));
		peds[iped].lnL1grid       = (double*)malloc2(nthetas   * sizeof(double));
	};
	thetas = (double*)malloc2(nthetas * sizeof(double));

	/*in input*/
	inData(FNIN, peds, npeds, thetas, nthetas, delimChar);

	/*perform calculations*/
	lnL0 = computeLnL0(peds, npeds);
	lnL1 = computeLnL1(peds, npeds, thetas, nthetas, &thetaMLE_H1);
	lnL2 = computeLnL2(peds, npeds, thetas, nthetas, &thetaMLE_H2, &alphaMLE_H2, ALPHAINCREMENT);

	/*out output*/
	outData(fpOut, peds, npeds, lnL0, lnL1, lnL2, thetas[0], thetaMLE_H1, thetaMLE_H2, alphaMLE_H2);

	fclose(fpOut);

	/*clean up*/
	for (iped = 0; iped < npeds; iped++) {
		free(peds[iped].id);
		free(peds[iped].lnL1grid);
	};
	free(peds);
	free(thetas);

	return NORMAL;
}



/* -----------------------------------------------------------------------------
reads in data from input file, which is assumed to be error-free.
see README file for file format.
*/
void
inData
(
char	*fn,
PedS	*peds,	/*(modified by reading in from input file)*/
int	npeds,
double	*thetas,/*(modified by reading in from input file)*/
int	nthetas,
char	delimChar
)
{
	FILE	*fp;
	char	*line;
	char	*word;
	int	iped;
	int	itheta;

	line = (char*)malloc2(MAXCHLINE * sizeof(char));
	/*word = (char*)malloc2(MAXCHWORD * sizeof(char));*/	
	/*TODO: why is this not necessary? does strtok do this already?*/

	fp = fopen2(fn, "r");

	if (delimChar == ',') {	/*comma-delimited file*/

		/*header line*/
		line = fgets(line, MAXCHLINE, fp);
		word = strtok(line, ",");
		word = strtok(NULL, ",");
		for (itheta=0; itheta < nthetas; itheta++) {
			word = strtok(NULL, ",");
			sscanf(word, "%lf", &thetas[itheta]);
		};

		/*pedigree data*/
		for (iped = 0; iped < npeds; iped++) {
			line = fgets(line, MAXCHLINE, fp);
			word = strtok(line, ",");
			sscanf(word, "%s", peds[iped].id);
			word = strtok(NULL, ",");
			sscanf(word, "%lf", &peds[iped].lnL0);
			for (itheta=0; itheta < nthetas; itheta++) {
				word = strtok(NULL, ",");
				sscanf(word, "%lf", &peds[iped].lnL1grid[itheta]);
			};
		};

	} else {	/*space-delimited file*/

		/*header line*/
		fscanf(fp, "%*s %*s");
		for (itheta=0; itheta < nthetas; itheta++) {
			fscanf(fp, "%lf", &thetas[itheta]);
		};
		line = fgets(line, MAXCHLINE, fp);

		/*pedigree data*/
		for (iped = 0; iped < npeds; iped++) {
			fscanf(fp, "%s %lf", peds[iped].id, &peds[iped].lnL0);
			for (itheta=0; itheta < nthetas; itheta++) {
				fscanf(fp, "%lf", &peds[iped].lnL1grid[itheta]);
			};
			line = fgets(line, MAXCHLINE, fp);
		};
	};

	fclose(fp);

	free(line);
	/*free(word);*/	/*TODO: why does this not work?*/

	return;
}



/* -----------------------------------------------------------------------------
computes ln likelihood under H0 (no linkage).
*/
double	/*ln L0*/
computeLnL0
(
PedS	*peds,	/*(modified: pedigree-specific ln L0's)*/
int	npeds
)
{
	double	lnL;
	int	iped;

	lnL = 0.0;
	for (iped = 0; iped < npeds; iped++) {
		lnL += peds[iped].lnL0;
	};

	return lnL;
}



/* -----------------------------------------------------------------------------
computes ln likelihood under H1 (linkage, homogeneity).
*/
double	/*ln L1*/
computeLnL1
(
PedS	*peds,		/*(modified: pedigree-specific ln L1's)*/
int	npeds,
double	*thetas,
int	nthetas,
double	*thetaMLE	/*(modified by initialization)*/
)
{
	double	lnL = 0.0,	/*max. ln L under H1 (linkage, homogeneity)*/
		lnLcurrent;
	int	positionOfThetaMLE = 0;
	int	itheta;
	int	iped;

	/*likelihood maximization*/
	for (itheta=0; itheta < nthetas; itheta++) {
		lnLcurrent = 0.0;
		for (iped = 0; iped < npeds; iped++) {
			lnLcurrent += peds[iped].lnL1grid[itheta];
		};
		if (itheta == 0 || lnL < lnLcurrent) {
			lnL = lnLcurrent;
			positionOfThetaMLE = itheta;
			*thetaMLE = thetas[positionOfThetaMLE];
		};
	};

	/*store pedigree-specific ln likelihoods at parameter MLEs*/
	for (iped = 0; iped < npeds; iped++) {
		peds[iped].lnL1 = peds[iped].lnL1grid[positionOfThetaMLE];
	};

	return lnL;
}



/* -----------------------------------------------------------------------------
computes ln likelihood under H2 (linkage, homogeneity).
*/
double	/*ln L2*/
computeLnL2
(
PedS	*peds,		/*(modified: pedigree-specific ln L2's)*/
int	npeds,
double	*thetas,
int	nthetas,
double	*thetaMLE,	/*(modified by initialization)*/
double	*alphaMLE,	/*(modified by initialization)*/
double	alphaIncrement
)
{
	double	lnL = 0.0,	/*max. ln L under H2 (linkage, heterogeneity)*/
		lnLcurrent;
	double	**scaledPedL1s;	/*scaled likelihood under H1 (linkage, homogeneity)
				  for each pedigree (1st dim.) and each considered value of linkage parameter (2nd dim.);
				  scaling is done to avoid under- or overflow*/
	int	positionOfThetaMLE = 0;
	int	itheta;
	double	alphaCurrent,
		oneMinusAlphaCurrent;
	int	iped;

	scaledPedL1s = mallocMatrixDouble(npeds, nthetas);

	/*scale pedigree-specific ln likelihoods to avoid under- or overflow;
	  take exp. to get pedigree-specific likelihoods*/
	for (iped = 0; iped < npeds; iped++) {
		for (itheta=0; itheta < nthetas; itheta++) {
			scaledPedL1s[iped][itheta] = exp(peds[iped].lnL1grid[itheta] - peds[iped].lnL1grid[0]);
		};
	};

	/*likelihood maximization*/
	for (alphaCurrent = 0.0; alphaCurrent <= 1.0; alphaCurrent+=alphaIncrement) {
		oneMinusAlphaCurrent = 1.0 - alphaCurrent;
		for (itheta=0; itheta < nthetas; itheta++) {
			lnLcurrent = 0.0;
			for (iped = 0; iped < npeds; iped++) {
				lnLcurrent += ln(alphaCurrent         * scaledPedL1s[iped][itheta] + 
						 oneMinusAlphaCurrent * scaledPedL1s[iped][0]);
			};
			if ((alphaCurrent <= 0.0 && itheta == 0) || lnL < lnLcurrent) {
				lnL = lnLcurrent;
				*alphaMLE = alphaCurrent;
				positionOfThetaMLE = itheta;
				*thetaMLE = thetas[positionOfThetaMLE];
			};
		};
	};

	/*un-scaled ln L*/
	for (iped = 0; iped < npeds; iped++) {
		lnL += peds[iped].lnL1grid[0];
	};

	/*store pedigree-specific ln likelihoods (un-scaled) at parameter MLEs*/
	for (iped = 0; iped < npeds; iped++) {
		peds[iped].lnL2 = ln(       *alphaMLE  * scaledPedL1s[iped][positionOfThetaMLE] + 
				     (1.0 - *alphaMLE) * scaledPedL1s[iped][0])
				  + peds[iped].lnL1grid[0];	/*un-scale*/
	};

	/*compute probabilities of linkage for all pedigrees*/
	for (iped = 0; iped < npeds; iped++) {
		peds[iped].pLinked =         *alphaMLE  * scaledPedL1s[iped][positionOfThetaMLE] /
				     (       *alphaMLE  * scaledPedL1s[iped][positionOfThetaMLE] + 
				      (1.0 - *alphaMLE) * scaledPedL1s[iped][0]);
	};

	/*clean up*/
	for (iped = 0; iped < npeds; iped++) {
		free(scaledPedL1s[iped]);
	};
	free(scaledPedL1s);

	return lnL;
}



/* -----------------------------------------------------------------------------
outputs output to output file.
*/
void
outData
(
FILE	*fp,
PedS	*peds,
int	npeds,
double	lnL0,		/*ln L under H0 (no linkage)*/
double	lnL1,		/*ln L under H1 (linkage, homogeneity)*/
double	lnL2,		/*ln L under H0 (linkage, heterogeneity)*/
double	theta_H0,	/*value of linkage parameter under H0 (no linkage)*/
double	thetaMLE_H1,	/*ML estimate of linkage parameter under H1 (linkage, homogeneity)*/	
double	thetaMLE_H2,	/*ML estimate of linkage parameter under H2 (linkage, heterogeneity)*/	
double	alphaMLE_H2	/*ML estimate of homogeneity parameter under H2 (linkage, heterogeneity)
			  alpha=1: all pedigrees are of linked type*/
)
{
	int	iped;
	char	*word;
	double	chi;
	double	p;

	word = (char*)malloc2(MAXCHWORD * sizeof(char));

	fprintf(fp, "description of hypotheses:\n");
	fprintf(fp, "\n");
	fprintf(fp, "%s %s\n", "abbr.", "       description");
	fprintf(fp, "%s %s\n", "-----", "-------------------------");
	fprintf(fp, "%s %s\n", "  H0 ", "no linkage");
	fprintf(fp, "%s %s\n", "  H1 ", "   linkage, homogeneity");
	fprintf(fp, "%s %s\n", "  H2 ", "   linkage, heterogeneity");
	fprintf(fp, "\n");
	fprintf(fp, "\n");

	fprintf(fp, "description of tests:\n");
	fprintf(fp, "\n");
	fprintf(fp, "%s %s %s %s %s\n", "  abbr.  ", "        description        ", "   lod score   ", " chi^2 statistic ", " theoretical asymptotic distribution");
	fprintf(fp, "%s %s %s %s %s\n", "---------", "---------------------------", "---------------", "-----------------", "------------------------------------");
	fprintf(fp, "%s %s %s %s %s\n", "H0 vs. H1", "  linkage under homogeneity", "lg[L(H1)/L(H0)]", "-2ln[L(H0)/L(H1)]", " .5 (0) + .5 chi^2(1)");
	fprintf(fp, "%s %s %s %s %s\n", "H1 vs. H2", "heterogeneity given linkage", "lg[L(H2)/L(H1)]", "-2ln[L(H1)/L(H2)]", " .5 (0) + .5 chi^2(1)");
	fprintf(fp, "%s %s %s %s %s\n", "H0 vs. H2", "  linkage and heterogeneity", "lg[L(H2)/L(H0)]", "-2ln[L(H0)/L(H2)]", ".25 (0) + .5 chi^2(1) + .25 chi^2(2) (?)");
	fprintf(fp, "\n");
	fprintf(fp, "\n");

	fprintf(fp, "test results:\n");
	fprintf(fp, "\n");
	fprintf(fp, "%s %s %s %s\n", "   test  ", "lod score", "chi^2 statistic", "   p-value");
	fprintf(fp, "%s %s %s %s\n", "---------", "---------", "---------------", "------------");

	chi = 2.0 * (lnL1 - lnL0);
	if (chi < 0.0) {	/*due to round-off errors or due to presence of nuisance parameters*/
		chi = 0.0;
	};
	p = 0.5 * (1.0 - chiCdf(chi, 1));
	fprintf(fp, "%s %9.6f %15.6f %12.6f\n", "H0 vs. H1", chisq2lod(chi), chi, p);

	chi = 2.0 * (lnL2 - lnL1);
	if (chi < 0.0) {	/*due to round-off errors or due to presence of nuisance parameters*/
		chi = 0.0;
	};
	p = 0.5 * (1.0 - chiCdf(chi, 1));
	fprintf(fp, "%s %9.6f %15.6f %12.6f", "H1 vs. H2", chisq2lod(chi), chi, p);
	if ((lnL1 - lnL0) / ln(10.0) >= 3.0) {
		fprintf(fp, "\n");
	} else {
		fprintf(fp, " (TEST MAKES NO SENSE AS H0 VS. H1 IS NOT SIGNIFICANT!)\n");
	};

	chi = 2.0 * (lnL2 - lnL0);
	if (chi < 0.0) {	/*due to round-off errors or due to presence of nuisance parameters*/
		chi = 0.0;
	};
	p   = 0.5 * (1.0 - chiCdf(chi, 1)) + 0.25 * (1.0 - chiCdf(chi, 2));
	fprintf(fp, "%s %9.6f %15.6f %12.6f\n",    "H0 vs. H2", chisq2lod(chi), chi, p);
	fprintf(fp, "\n");
	fprintf(fp, "\n");

	fprintf(fp, "maximum likelihoods and maximum likelihood estimates of parameters:\n");
	fprintf(fp, "\n");
	fprintf(fp, "%10s %13s %12s %12s\n", "hypothesis", "ln likelihood", "linkage par.", " homog. par.");
	fprintf(fp, "%10s %13s %12s %12s\n", "----------", "-------------", "------------", "------------");
	sprintf(word, "(%.6f)", theta_H0);
	fprintf(fp, "%10s %13.6f %12s %12s\n", "    H0    ", lnL0, word, "(undefined)");
	fprintf(fp, "%10s %13.6f  %10.6f   (%8.6f)\n", "    H1    ", lnL1, thetaMLE_H1, 1.0);
	fprintf(fp, "%10s %13.6f  %10.6f  %10.6f\n", "    H2    ", lnL2, thetaMLE_H2, alphaMLE_H2);
	fprintf(fp, "\n");
	fprintf(fp, "\n");

	fprintf(fp, "pedigree-specific results:\n");
	fprintf(fp, "\n");
	fprintf(fp, "%8s %38s %29s %10s\n", 
		"        ", "             ln likelihood            ", "          lod score          ",             "prob. that");
	fprintf(fp, "%8s %38s %29s %10s\n", 
		"        ", "--------------------------------------", "-----------------------------",             "  pedigree");
	fprintf(fp, "%8s %12s %12s %12s %9s %9s %9s %10s\n",
		"pedigree", "     H0     ", "     H1     ", "     H2     ", "H0 vs. H1", "H1 vs. H2", "H0 vs. H2", " is linked");
	fprintf(fp, "%8s %11s %11s %11s %9s %9s %9s %10s\n",
		"--------", "------------", "------------", "------------", "---------", "---------", "---------", "----------"); 
	for (iped = 0; iped < npeds; iped++) {
		fprintf(fp, "%8s %12.6f %12.6f %12.6f %9.6f %9.6f %9.6f %8.3f\n",
			peds[iped].id, 
			peds[iped].lnL0, 
			peds[iped].lnL1, 
			peds[iped].lnL2, 
			(peds[iped].lnL1 - peds[iped].lnL0) / ln(10.0), 
			(peds[iped].lnL2 - peds[iped].lnL1) / ln(10.0), 
			(peds[iped].lnL2 - peds[iped].lnL0) / ln(10.0), 
			peds[iped].pLinked);
	};
	fprintf(fp, "%8s %12.6f %12.6f %12.6f %9.6f %9.6f %9.6f\n",
		"total",
		lnL0,
		lnL1,
		lnL2,
		(lnL1 - lnL0) / ln(10.0), 
		(lnL2 - lnL1) / ln(10.0), 
		(lnL2 - lnL0) / ln(10.0));
	
	free(word);

	return;
}
