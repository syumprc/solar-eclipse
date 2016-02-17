/*
 * mibd.cc implements the mibd command
 * Written by Thomas Dyer January 1998
 * Copyright (c) 1998 Southwest Foundation for Biomedical Research
 */

#define MXNREL	10	/* default maximum number of times a pair of
                           individuals can share the same relationship  */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <sys/stat.h>
#include "solar.h"
#include "tcl.h"
#include "safelib.h"

extern bool  XLinked;
extern float MibdWin;
extern bool  MMSibs;

static int run_relate (Tcl_Interp*, int);
static int run_merge (Tcl_Interp*);
static int run_means (Tcl_Interp*, bool);

// DEC doesn't have this header file
//#include <sunmath.h>
#include "math.h"

extern "C" int MibdCmd (ClientData clientData, Tcl_Interp *interp, int argc,
                        char *argv[])
{
    if (argc == 2 && !StringCmp ("help", argv[1], case_ins)) {
        return Solar_Eval (interp, "help mibd");
    }
 
    else if (argc == 2 && !StringCmp ("merge", argv[1], case_ins)) {
    // merge marker IBDs with relative-class data
        if (MMSibs) {
            RESULT_LIT (\
"This operation is not defined when the MMSibs option has been chosen.");
            return TCL_ERROR;
        }

        if (XLinked) {
            RESULT_LIT (\
"SOLAR does not yet support multipoint IBD calculation for X-linked loci\n\
in extended pedigrees.");
            return TCL_ERROR;
        }

        try {
            if (!loadedPed()) {
                RESULT_LIT ("Pedigree data have not been loaded.");
                return TCL_ERROR;
            }
        }
        catch (Safe_Error_Return& ser) {
            RESULT_BUF (ser.message());
            Solar_AppendResult1(interp,
                        "\nReloading the pedigree data is advised.", NULL);
            return TCL_ERROR;
        }
 
        if (!currentMap) {
            RESULT_LIT ("Map data have not been loaded.");
            return TCL_ERROR;
        }

        if (run_merge(interp) == TCL_ERROR) {
            return TCL_ERROR;
        }
 
        return TCL_OK;
    }

    else if (argc == 2 && !StringCmp ("means", argv[1], case_ins) ||
             argc == 3 && !StringCmp ("means", argv[1], case_ins) &&
             (!strcmp(argv[2], "-all") || !strcmp(argv[2], "-typed")))
    {
    // compute mean IBD by relative-class
        bool typed_only = false;
        if (argc == 3 && !strcmp(argv[2], "-typed"))
            typed_only = true;

        if (MMSibs) {
            RESULT_LIT (\
      "This operation is not defined when the MMSibs option has been chosen.");
            return TCL_ERROR;
        }

        if (XLinked) {
            RESULT_LIT (\
"SOLAR does not yet support multipoint IBD calculation for X-linked loci\n\
in extended pedigrees.");
            return TCL_ERROR;
        }

        try {
            if (!loadedPed()) {
                RESULT_LIT ("Pedigree data have not been loaded.");
                return TCL_ERROR;
            }
        }
        catch (Safe_Error_Return& ser) {
            RESULT_BUF (ser.message());
            Solar_AppendResult1(interp,
                        "\nReloading the pedigree data is advised.", NULL);
            return TCL_ERROR;
        }
 
        if (!currentMap) {
            RESULT_LIT ("Map data have not been loaded.");
            return TCL_ERROR;
        }

        if (run_means(interp, typed_only) == TCL_ERROR) {
            return TCL_ERROR;
        }
 
        return TCL_OK;
    }

    else if (argc == 2 || argc == 4) {
    // compute multipoint IBDs
        try {
            if (!loadedPed()) {
                RESULT_LIT ("Pedigree data have not been loaded.");
                return TCL_ERROR;
            }
        }
        catch (Safe_Error_Return& ser) {
            RESULT_BUF (ser.message());
            Solar_AppendResult1(interp,
                        "\nReloading the pedigree data is advised.", NULL);
            return TCL_ERROR;
        }
 
        if (!currentMap) {
            RESULT_LIT ("Map data have not been loaded.");
            return TCL_ERROR;
        }

        if (MMSibs) {
            try {
                if (!currentPed->marker()) {
                    RESULT_LIT ("Marker data have not been loaded.");
                    return TCL_ERROR;
                }
            }
            catch (Safe_Error_Return& ser) {
                RESULT_BUF (ser.message());
                Solar_AppendResult1(interp,
                            "\nReloading the marker data is advised.", NULL);
                return TCL_ERROR;
            }
        }

        int i;
        char mibddir[1024];

        if (!MMSibs) {
            if (Solar_Eval(interp, "mibddir -session") == TCL_ERROR)
                return TCL_ERROR;
            else
                strcpy(mibddir, Tcl_GetStringResult (interp));

            FILE *infp;
            infp = fopen("mibdrel.ped", "r");
            if (!infp) {
                printf("Creating relative-class file ...\n");
                fflush(stdout);
                if (run_relate(interp, MXNREL) == TCL_ERROR)
                    return TCL_ERROR;
            }
            else
                fclose(infp);

            char infile[1024];
            struct stat mrgsbuf;
            sprintf(infile, "%s/mibdchr%s.mrg.gz", mibddir,
                    currentMap->chrnum());
            infp = fopen(infile, "r");
            if (!infp)
            {
                printf("Merging marker IBDs ...\n");
                fflush(stdout);
                if (run_merge(interp) == TCL_ERROR)
                    return TCL_ERROR;
            }
            else
            {
                fclose(infp);
                stat(infile, &mrgsbuf);

                char ibddir[1024];
                if (Solar_Eval(interp, "ibddir") == TCL_ERROR)
                    return TCL_ERROR;
                else
                    strcpy(ibddir, Tcl_GetStringResult (interp));

                struct stat ibdsbuf;
                for (i = 0; i < currentMap->nloci(); i++) {
                    sprintf(infile, "%s/ibd.%s.gz", ibddir,
                            currentMap->mrkname(i));
                    infp = fopen(infile, "r");
                    if (infp)
                    {
                        fclose(infp);
                        stat(infile, &ibdsbuf);
                        if (mrgsbuf.st_mtime < ibdsbuf.st_mtime)
                        {
                            RESULT_LIT (\
"One or more IBD files are newer than the merged marker-IBD file.\n\
Previously computed multipoint IBD files may be out of date.\n\
You must enter the 'mibd merge' command before computing multipoint IBDs.");
                            return TCL_ERROR;
                        }
                    }
                }
            }

            struct stat meansbuf;
            sprintf(infile, "%s/mibdchr%s.mean", mibddir,
                    currentMap->chrnum());
            infp = fopen(infile, "r");
            if (!infp)
            {
                printf("Computing mean IBD by relative-class ...\n");
                fflush(stdout);
                if (run_means(interp, false) == TCL_ERROR)
                    return TCL_ERROR;
            }
            else
            {
                fclose(infp);
                stat(infile, &meansbuf);
                if (meansbuf.st_mtime < mrgsbuf.st_mtime)
                {
                    printf("Computing mean IBD by relative-class ...\n");
                    fflush(stdout);
                    if (run_means(interp, false) == TCL_ERROR)
                        return TCL_ERROR;
                }
            }
        }

        if (XLinked) {
            RESULT_LIT (\
"SOLAR does not yet support multipoint IBD calculation for X-linked loci\n\
in extended pedigrees.");
            return TCL_ERROR;
        }

        double last = ceil(currentMap->mrklocn(currentMap->nloci() - 1));

        double from, to, incr;
        char errmsg[1024];

        if (argc == 4) {
            if (sscanf(argv[1], "%lf", &from) != 1 || from < 0)
            {
                RESULT_LIT ("<from> must be a non-negative number");
                return TCL_ERROR;
            }
            if (sscanf(argv[2], "%lf", &to) != 1 || to < 0)
            {
                RESULT_LIT ("<to> must be a non-negative number");
                return TCL_ERROR;
            }
            if (sscanf(argv[3], "%lf", &incr) != 1 || incr <= 0)
            {
                RESULT_LIT ("<incr> must be a positive number");
                return TCL_ERROR;
            }
        }

        else {
            if (sscanf(argv[1], "%lf", &incr) != 1 || incr <= 0)
            {
                RESULT_LIT ("<incr> must be a positive number");
                return TCL_ERROR;
            }
            from = 0;
            to = last;
        }

        char mibd_cmd[1024];
        if (MMSibs) {
            sprintf(mibd_cmd, "exec dommsibs y %s %f %f %f", 
                    currentMap->chrnum(), from, to, incr);

            if (Solar_Eval(interp, mibd_cmd) == TCL_ERROR) {
                char errfile[1024];
                sprintf(errfile, "dommsibs.err");
                FILE *errfp = fopen(errfile, "r");
                if (errfp && fgets(errmsg, sizeof(errmsg), errfp)) {
                    fclose(errfp);
                    RESULT_BUF (strtok (errmsg, "\n"));
                }
                else
                    RESULT_LIT ("dommsibs failed");
                return TCL_ERROR;
            }

            return TCL_OK;
        }

        char fname[1024];
        sprintf(fname, "%s/mibdchr%s.loc", mibddir, currentMap->chrnum());

        FILE *fp = fopen(fname, "r");
        if (fp) {
            char rec[1024], errmsg[1024];
            if (!fgets(rec, sizeof(rec), fp)) {
                sprintf(errmsg, "Error reading %s.", fname);
                RESULT_BUF (errmsg);
                return TCL_ERROR;
            }
            fclose(fp);

            char map_func = 'k';
            if (strncmp(rec, "NLOCI", 5)) {
                if (!strcmp(rec, "Kosambi\n"))
                    map_func = 'k';
                else if (!strcmp(rec, "Haldane\n"))
                    map_func = 'h';
                else {
                    sprintf(errmsg,
                            "Unrecognized mapping function in %s.", fname);
                    RESULT_BUF (errmsg);
                    return TCL_ERROR;
                }
            }

            if (currentMap->mapfunc() != map_func) {
                if (map_func == 'k') {
                    RESULT_LIT (
"Existing multipoint IBDs were computed using the Kosambi mapping function.\n\
If you want to use the Kosambi mapping function, reload the map file.\n\
If not, you must delete the old multipoint IBDs or specify a new mibddir.");
                    return TCL_ERROR;
                }
                else if (map_func == 'h') {
                    RESULT_LIT (
"Existing multipoint IBDs were computed using the Haldane mapping function.\n\
If you want to use the Haldane mapping function, reload the map file.\n\
If not, you must delete the old multipoint IBDs or specify a new mibddir.");
                    return TCL_ERROR;
                }
            }
        }

        fp = fopen(fname, "w");
        if (!fp) {
            sprintf(errmsg, "Cannot open %s.", fname);
            RESULT_BUF (errmsg);
            return TCL_ERROR;
        }

        if (currentMap->mapfunc() == 'h')
            fprintf(fp, "Haldane\n");
        else
            fprintf(fp, "Kosambi\n");

        for (i = 0; i < currentMap->nloci(); i++)
            fprintf(fp, "%-11.11s  %6.1f\n",
                    currentMap->mrkname(i), currentMap->mrklocn(i));
        fclose(fp);

        double locn;
        for (locn = from; locn <= to + 0.000001; locn += incr)
        {
            for (i = 0; i < currentMap->nloci(); i++) {
                if (currentMap->mrklocn(i) >= locn - .5*MibdWin &&
                    currentMap->mrklocn(i) <= locn + .5*MibdWin)
                {
                    break;
                }
            }
            if (i == currentMap->nloci()) {
                sprintf(errmsg,
"There are no markers inside the MIBD window centered at location %g.\n\
Use the ibdoption command to increase MibdWin.",
                        locn);
                RESULT_BUF (errmsg);
                return TCL_ERROR;
            }
        }

        char show_status = 'n';
        FILE *ttyfp = fopen("/dev/tty", "w");
        if (ttyfp) {
            show_status = 'y';
            fprintf(ttyfp, "Computing multi-point IBDs: ");
            fclose(ttyfp);
        }

        for (locn = from; locn <= to + 0.000001; locn += incr)
        {
            char *p, clocn[100];
            sprintf(clocn, "%f", locn);
            for (p = clocn + strlen(clocn) - 1; p > clocn; p--) {
                if (*p == '0')
                    *p = '\0';
                else
                    break;
            }
            for ( ; p > clocn; p--) {
                if (*p == '.')
                    *p = '\0';
                else
                    break;
            }

            sprintf(mibd_cmd,
"exec multipnt %s/mibdchr%s.loc %s/mibdchr%s.mrg.gz %s/mibdchr%s.mean %f %s %f %c >& multipnt.out",
                    mibddir, currentMap->chrnum(),
                    mibddir, currentMap->chrnum(),
                    mibddir, currentMap->chrnum(),
                    locn, clocn, MibdWin, show_status);

            if (Solar_Eval(interp, mibd_cmd) == TCL_ERROR) {
                printf("\n"); fflush(stdout);
                if (Solar_Eval(interp, "exec cat multipnt.out") == TCL_ERROR)
                {
                    RESULT_LIT ("Cannot cat multipnt.out");
                    return TCL_ERROR;
                }
                if (!strlen(Tcl_GetStringResult (interp)))
                    RESULT_LIT ("Program multipnt did not run.");
                return TCL_ERROR;
            }

            unlink("multipnt.out");

            sprintf(fname, "mibd.%s.%s", currentMap->chrnum(), clocn);
            sprintf(mibd_cmd, "exec mv mibd.out %s/%s", mibddir, fname);
            if (Solar_Eval(interp, mibd_cmd) == TCL_ERROR) {
                RESULT_LIT ("Cannot rename mibd.out");
                return TCL_ERROR;
            }

            sprintf(mibd_cmd, "exec gzip -f %s/%s", mibddir, fname);
            if (Solar_Eval(interp, mibd_cmd) == TCL_ERROR) {
                RESULT_LIT ("gzip failed");
                return TCL_ERROR;
            }

	    sprintf (mibd_cmd, "matcrc %s/%s.gz", mibddir, fname);
	    if (Solar_Eval(interp, mibd_cmd) == TCL_ERROR) {
		return TCL_ERROR;
	    }

            FILE *ttyfp = fopen("/dev/tty", "w");
            if (ttyfp) {
                int i;
                for (i = 0; i < 10 + strlen(clocn); i++)
                    fputc('\b', ttyfp);
                for (i = 0; i < 10 + strlen(clocn); i++)
                    fputc(' ', ttyfp);
                for (i = 0; i < 10 + strlen(clocn); i++)
                    fputc('\b', ttyfp);
                fclose(ttyfp);
            }
        }

        printf("\n"); fflush(stdout);
        return TCL_OK;
    }

    RESULT_LIT ("Invalid mibd command");
    return TCL_ERROR;
}

int run_relate (Tcl_Interp *interp, int mxnrel)
{
    char run_cmd[1024];
    sprintf(run_cmd, "exec relate %d > relate.out", mxnrel);

    if (Solar_Eval(interp, run_cmd) == TCL_ERROR) {
        char buf[1024];
        FILE *fp = fopen("relate.out", "r");
        if (fp) {
            if (fgets(buf, sizeof(buf), fp)) {
                unlink("mibdrel.ped");
                unlink("relate.unk");
                RESULT_LIT (\
               "Program relate failed. Check the file relate.out for errors.");
                fclose(fp);
            }
            else {
                unlink("relate.out");
                if (fp = fopen("relate.unk", "r")) {
                    RESULT_LIT (\
"Relationships of an unknown type were found. SOLAR cannot compute multipoint\n\
IBDs for this pedigree structure. A detailed list of these relationships can\n\
be found in the file \"relate.unk\".");
                    fclose(fp);
                }
            }
        }
        else
            RESULT_LIT ("Program relate did not run.");

        return TCL_ERROR;
    }

    unlink("relate.out");

    return TCL_OK;
}

int run_merge (Tcl_Interp *interp)
{
    char ibddir[1024];
    if (Solar_Eval(interp, "ibddir") == TCL_ERROR)
        return TCL_ERROR;
    else
        strcpy(ibddir, Tcl_GetStringResult (interp));

    char mibddir[1024];
    if (Solar_Eval(interp, "mibddir -session") == TCL_ERROR)
        return TCL_ERROR;
    else
        strcpy(mibddir, Tcl_GetStringResult (interp));

    char merge_cmd[1024];
    sprintf(merge_cmd, "exec mrgibd %s %s %s > mrgibd.out",
            currentMap->filename(), ibddir, mibddir);

    if (Solar_Eval(interp, merge_cmd) == TCL_ERROR) {
        if (Solar_Eval(interp, "exec cat mrgibd.out") == TCL_ERROR) {
            RESULT_LIT ("Cannot cat mrgibd.out");
            return TCL_ERROR;
        }
        if (!strlen(Tcl_GetStringResult (interp)))
            RESULT_LIT ("Program mrgibd did not run.");
        return TCL_ERROR;
    }

    unlink("mrgibd.out");
    return TCL_OK;
}

int run_means (Tcl_Interp *interp, bool typed_only)
{
    char ibddir[1024];
    if (Solar_Eval(interp, "ibddir") == TCL_ERROR)
        return TCL_ERROR;
    else
        strcpy(ibddir, Tcl_GetStringResult (interp));

    char mibddir[1024];
    if (Solar_Eval(interp, "mibddir -session") == TCL_ERROR)
        return TCL_ERROR;
    else
        strcpy(mibddir, Tcl_GetStringResult (interp));

    char means_cmd[1024];
    sprintf(means_cmd,
    "exec getmeans %s/mibdchr%s.mrg.gz %s/mibdchr%s.mean %d %c > getmeans.out",
            mibddir, currentMap->chrnum(), mibddir, currentMap->chrnum(),
            currentMap->nloci(), typed_only ? 'y' : 'n');

    if (Solar_Eval(interp, means_cmd) == TCL_ERROR) {
        if (Solar_Eval(interp, "exec cat getmeans.out") == TCL_ERROR) {
            RESULT_LIT ("Cannot cat getmeans.out");
            return TCL_ERROR;
        }
        if (!strlen(Tcl_GetStringResult (interp)))
            RESULT_LIT ("Program getmeans did not run.");
        return TCL_ERROR;
    }

    unlink("getmeans.out");
    return TCL_OK;
}
