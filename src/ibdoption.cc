/*
 * ibdoption.cc implements the ibdoption command
 * Written by Thomas Dyer January 1998
 * Copyright (c) 1998 Southwest Foundation for Biomedical Research
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <sys/stat.h>
#include "solar.h"
#include "tcl.h"
#include "safelib.h"

bool  NoMLE = false;
bool  XLinked = false;
bool  MCarlo = false;
bool  MaxRisk = true;
int   NumImp = 200;
float MibdWin = 1000.;
bool  MMSibs = false;

extern "C" int IbdOptCmd (ClientData clientData, Tcl_Interp *interp, int argc,
                          char *argv[])
{
    if (argc == 2 && !StringCmp ("help", argv[1], case_ins)) {
        return Solar_Eval (interp, "help ibdoption");
    }

    try {
        loadedPed();
    }
    catch (Safe_Error_Return& ser) {
        RESULT_BUF (ser.message());
        Solar_AppendResult1(interp,
                    "\nReloading the pedigree data is advised.", NULL);
        return TCL_ERROR;
    }

    if (argc == 1) {
    // display current settings of all IBD options
        char buf[1024];
        buf[0] = '\0';
        if (Solar_Eval(interp, "ibdoption xlinked ?") == TCL_OK)
        {
            strcat(buf, Tcl_GetStringResult (interp));
            strcat(buf, "\n");
        }
        if (Solar_Eval(interp, "ibdoption nomle ?") == TCL_OK)
        {
            strcat(buf, Tcl_GetStringResult (interp));
            strcat(buf, "\n");
        }
        if (Solar_Eval(interp, "ibdoption mcarlo ?") == TCL_OK)
        {
            strcat(buf, Tcl_GetStringResult (interp));
            strcat(buf, "\n");
        }
        if (Solar_Eval(interp, "ibdoption mcarlo #") == TCL_OK)
        {
            strcat(buf, Tcl_GetStringResult (interp));
            strcat(buf, "\n");
        }
        if (Solar_Eval(interp, "ibdoption mcarlo max ?") == TCL_OK)
        {
            strcat(buf, Tcl_GetStringResult (interp));
            strcat(buf, "\n");
        }
        if (Solar_Eval(interp, "ibdoption mibdwin") == TCL_OK)
//        {
            strcat(buf, Tcl_GetStringResult (interp));
//            strcat(buf, "\n");
//        }
//        if (Solar_Eval(interp, "ibdoption mmsibs ?") == TCL_OK)
//            strcat(buf, Tcl_GetStringResult (interp));
        RESULT_BUF (buf);
        return TCL_OK;
    }

    else if (argc >= 2 && !StringCmp ("xlinked", argv[1], case_ins)) {
        if (argc == 2) {
        // toggle XLinked option, display new setting
            if (XLinked) {
                XLinked = false;
                RESULT_LIT ("XLinked = no");
            }
            else {
                XLinked = true;
                RESULT_LIT ("XLinked = yes");
            }
            return TCL_OK;
        }

        else if (argc == 3 && !strcmp(argv[2], "y")) {
        // set XLinked option to yes
            XLinked = true;
            RESULT_LIT ("XLinked = yes");
            return TCL_OK;
        }

        else if (argc == 3 && !strcmp(argv[2], "n")) {
        // set XLinked option to no
            XLinked = false;
            RESULT_LIT ("XLinked = no");
            return TCL_OK;
        }

        else if (argc == 3 && !strcmp(argv[2], "?")) {
        // display current setting of XLinked option
            if (XLinked)
                RESULT_LIT ("XLinked = yes");
            else
                RESULT_LIT ("XLinked = no");
            return TCL_OK;
        }

        RESULT_LIT ("Usage: ibdoption xlinked [y | n | ?]");
        return TCL_ERROR;
    }

    else if (argc >= 2 && !StringCmp ("nomle", argv[1], case_ins)) {
        if (argc == 2) {
        // toggle NoMLE option, display new setting
            if (NoMLE) {
                NoMLE = false;
                RESULT_LIT ("NoMLE = no");
            }
            else {
                NoMLE = true;
                RESULT_LIT ("NoMLE = yes");
            }
            return TCL_OK;
        }

        else if (argc == 3 && !strcmp(argv[2], "y")) {
        // set NoMLE option to yes
            NoMLE = true;
            RESULT_LIT ("NoMLE = yes");
            return TCL_OK;
        }

        else if (argc == 3 && !strcmp(argv[2], "n")) {
        // set NoMLE option to no
            NoMLE = false;
            RESULT_LIT ("NoMLE = no");
            return TCL_OK;
        }

        else if (argc == 3 && !strcmp(argv[2], "?")) {
        // display current setting of NoMLE option
            if (NoMLE)
                RESULT_LIT ("NoMLE = yes");
            else
                RESULT_LIT ("NoMLE = no");
            return TCL_OK;
        }

        RESULT_LIT ("Usage: ibdoption nomle [y | n | ?]");
        return TCL_ERROR;
    }

    else if (argc >= 2 && !StringCmp ("mcarlo", argv[1], case_ins)) {
        if (argc == 2) {
        // toggle MCarlo option, display new setting
            if (MCarlo) {
                MCarlo = false;
                RESULT_LIT ("MCarlo = no");
            }
            else {
                MCarlo = true;
                RESULT_LIT ("MCarlo = yes");
            }
            return TCL_OK;
        }

        else if (argc == 3 && !strcmp(argv[2], "y")) {
        // set MCarlo option to yes
            MCarlo = true;
            RESULT_LIT ("MCarlo = yes");
            return TCL_OK;
        }

        else if (argc == 3 && !strcmp(argv[2], "n")) {
        // set MCarlo option to no
            MCarlo = false;
            RESULT_LIT ("MCarlo = no");
            return TCL_OK;
        }

        else if (argc == 3 && !strcmp(argv[2], "?")) {
        // display current setting of MCarlo option
            if (MCarlo)
                RESULT_LIT ("MCarlo = yes");
            else
                RESULT_LIT ("MCarlo = no");
            return TCL_OK;
        }

        else if (argc == 4 && !strcmp(argv[2], "#")) {
        // set number of imputations
            int n;
            if (sscanf(argv[3], "%d", &n) != 1 || n <= 0) {
                RESULT_LIT ("Invalid number of imputations");
                return TCL_ERROR;
            }
            NumImp = n;
            char buf[1024];
            sprintf(buf, "NumImp = %d", NumImp);
            RESULT_BUF (buf);
            return TCL_OK;
        }

        else if (argc == 3 && !strcmp(argv[2], "#")) {
        // display number of imputations
            char buf[1024];
            sprintf(buf, "NumImp = %d", NumImp);
            RESULT_BUF (buf);
            return TCL_OK;
        }

        else if (argc == 4 && !strcmp(argv[2], "max") &&
                 !strcmp(argv[3], "y")) {
        // set max risk option to yes
            MaxRisk = true;
            RESULT_LIT ("MaxRisk = yes");
            return TCL_OK;
        }

        else if (argc == 4 && !strcmp(argv[2], "max") &&
                 !strcmp(argv[3], "n")) {
        // set max risk option to no
            MaxRisk = false;
            RESULT_LIT ("MaxRisk = no");
            return TCL_OK;
        }

        else if (argc == 4 && !strcmp(argv[2], "max") &&
                 !strcmp(argv[3], "?")) {
        // display current setting of max risk option
            if (MaxRisk)
                RESULT_LIT ("MaxRisk = yes");
            else
                RESULT_LIT ("MaxRisk = no");
            return TCL_OK;
        }

        RESULT_LIT (
        "Usage: ibdoption mcarlo [y | n | ? | # [<#imp>] | max <y|n|?>]");
        return TCL_ERROR;
    }

    else if (argc >= 2 && !StringCmp ("mibdwin", argv[1], case_ins)) {
        if (argc == 3) {
        // set multipoint IBD window
            float w;
            if (sscanf(argv[2], "%f", &w) != 1 || w <= 0) {
                RESULT_LIT ("Invalid window size");
                return TCL_ERROR;
            }
            MibdWin = w;
            char buf[1024];
            sprintf(buf, "MibdWin = %g", MibdWin);
            RESULT_BUF (buf);
            return TCL_OK;
        }

        else if (argc == 2) {
        // display multipoint IBD window
            char buf[1024];
            sprintf(buf, "MibdWin = %g", MibdWin);
            RESULT_BUF (buf);
            return TCL_OK;
        }

        RESULT_LIT ("Usage: ibdoption mibdwin [<size>]");
        return TCL_ERROR;
    }

    else if (argc >= 2 && !StringCmp ("mmsibs", argv[1], case_ins)) {

// MAPMAKER/SIBS processing not yet fully implemented
        RESULT_LIT ("Sorry, MAPMAKER/SIBS interface not yet implemented.");
        return TCL_ERROR;
/*
        if (argc == 2) {
        // toggle MMSibs option, display new setting
            if (MMSibs) {
                MMSibs = false;
                RESULT_LIT ("MMSibs = no");
            }
            else {
                MMSibs = true;
                RESULT_LIT ("MMSibs = yes");
            }
            return TCL_OK;
        }

        else if (argc == 3 && !strcmp(argv[2], "y")) {
        // set MMSibs option to yes
            MMSibs = true;
            RESULT_LIT ("MMSibs = yes");
            return TCL_OK;
        }

        else if (argc == 3 && !strcmp(argv[2], "n")) {
        // set MMSibs option to no
            MMSibs = false;
            RESULT_LIT ("MMSibs = no");
            return TCL_OK;
        }

        else if (argc == 3 && !strcmp(argv[2], "?")) {
        // display current setting of MMSibs option
            if (MMSibs)
                RESULT_LIT ("MMSibs = yes");
            else
                RESULT_LIT ("MMSibs = no");
            return TCL_OK;
        }

        RESULT_LIT ("Usage: ibdoption mmsibs [y | n | ?]");
        return TCL_ERROR;
*/
    }

    RESULT_LIT ("Invalid ibdoption command");
    return TCL_ERROR;
}
