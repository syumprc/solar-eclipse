/*
 * drand.cc provides a Tcl interface to the C library random number generator
 * Written by Thomas Dyer October 1999
 * Copyright (c) 1999 Southwest Foundation for Biomedical Research
 */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include "solar.h"
#include "safelib.h"

bool random_number_generator_seeded = false;

extern "C" int DrandCmd (ClientData clientData, Tcl_Interp *interp,
		         int argc, char *argv[])
{
    if (argc == 2 && !StringCmp ("help", argv[1], case_ins)) {
        return Solar_Eval (interp, "help drand");
    }

    else if (argc == 2) {
        long seed;
        if (sscanf(argv[1], "%ld", &seed) != 1) {
            RESULT_LIT ("drand: seed must be an integer");
            return TCL_ERROR;
        }

        if (!seed)
            srand48(time(NULL));
        else
            srand48(seed);

        random_number_generator_seeded = true;
        return TCL_OK;
    }

    else if (argc == 1) {
        if (!random_number_generator_seeded)
            Solar_Eval(interp, "drand 0");

        char buf[128];
        sprintf(buf, "%f", drand48());
        RESULT_BUF (buf);
        return TCL_OK;
    }

    RESULT_LIT ("Usage: drand [ <seed> ]");
    return TCL_ERROR;
}
