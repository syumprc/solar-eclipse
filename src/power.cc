/*
 * power.cc implements the Power command
 * Written by Thomas Dyer beginning on September 5, 2003
 * Copyright (c) 2003 Southwest Foundation for Biomedical Research
 */

#include "safelib.h"
// tcl.h from solar.h
#include "solar.h"

extern "C" void fit2dp_ (int*, double*, double*, double*, int*, double*,
                         int*, double*);

extern "C" int PowerCmd (ClientData clientData, Tcl_Interp *interp,
	    int argc, char *argv[])
{
    if (argc == 2 && !StringCmp ("help", argv[1], case_ins))
    {
	return Solar_Eval (interp, "help power");
    }

    if (argc == 3 && !StringCmp ("smooth_elods", argv[1], case_ins))
    {
        int nh2q;
        if (sscanf(argv[2], "%d", &nh2q) != 1) {
            RESULT_LIT ("Number of grid points must be an integer");
            return TCL_ERROR;
        }

        FILE *infp = fopen("power.out.tmp", "r");
        if (!infp) {
            RESULT_LIT ("Cannot open file power.out.tmp");
            return TCL_ERROR;
        }

        double *h2q = (double *) malloc(3*nh2q*sizeof(double));
        double *elod = (double *) malloc(nh2q*sizeof(double));

        int i;
        double th2q;
        for (i = 0; i < nh2q; i++) {
            if (fscanf(infp, "%lf %lf", &th2q, &elod[i]) != 2) {
                fclose(infp);
                free(elod);
                free(h2q);
                RESULT_LIT ("Error reading power.out.tmp");
                return TCL_ERROR;
            }
            h2q[i] = 1;
            h2q[i+nh2q] = th2q;
            h2q[i+2*nh2q] = pow(th2q, 2.);
        }
        fclose(infp);

        double coeff[3], work[9];
        double *work2 = (double *) malloc(3*nh2q*sizeof(double));
        int info, *ipvt = (int *) malloc(nh2q*sizeof(int));

        fit2dp_(&nh2q, h2q, elod, work, ipvt, work2, &info, coeff);

        free(ipvt);
        free(work2);
        free(elod);
        free(h2q);

        if (info) {
            RESULT_LIT ("Matrix inversion failed");
            return TCL_ERROR;
        }

        FILE *outfp = fopen("power.out.tmp", "w");
        if (!outfp) {
            RESULT_LIT ("Cannot write to file power.out.tmp");
            return TCL_ERROR;
        }

        for (th2q = 0.01; th2q <= 1; th2q += 0.01) {
            fprintf(outfp, "%g %g\n", th2q,
                    coeff[0] + coeff[1]*th2q + coeff[2]*pow(th2q, 2.));
        }
        fclose(outfp);

	return TCL_OK;
    }

    RESULT_LIT ( "Invalid power command");
    return TCL_ERROR;
}
