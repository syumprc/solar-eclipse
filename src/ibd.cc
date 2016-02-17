/*
 * ibd.cc implements the ibd command
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

extern bool  NoMLE;
extern bool  XLinked;
extern bool  MCarlo;
extern bool  MaxRisk;
extern int   NumImp;
extern bool  MMSibs;

#define MXMCALL	99	/* max alleles handled by Monte Carlo IBD method */

static int inf_ibd (const char*, const char*, Tcl_Interp*);
static int mito_ibd (const char*, Tcl_Interp*);
static int run_ibd (const char*, bool, const char *, Tcl_Interp*);
static int run_mc (const char*, bool, const char *, Tcl_Interp*);
static int (*ibd_func)(const char*, bool, const char*, Tcl_Interp*) = 0;

extern "C" int IbdCmd (ClientData clientData, Tcl_Interp *interp, int argc,
                       char *argv[])
{
    if (argc == 2 && !StringCmp ("help", argv[1], case_ins)) {
        return Solar_Eval (interp, "help ibd");
    }
 
    else if (argc == 2 && !StringCmp ("mito", argv[1], case_ins)) {
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
 
        char ibddir[1024];
        if (Solar_Eval(interp, "ibddir -session") == TCL_ERROR)
            return TCL_ERROR;
        else
            strcpy(ibddir, Tcl_GetStringResult (interp));

        return mito_ibd(ibddir, interp);
    }

    else if (argc <= 2 || argc == 3 && !strcmp("-nomle", argv[1])) {
        char mrkname[1024];
        bool doall = false;
        bool nomle = NoMLE;

        if (argc == 1)
            doall = true;
        else if (argc == 2 && !strcmp("-nomle", argv[1])) {
            doall = true;
            nomle = true;
        }
        else if (argc == 2)
            strcpy(mrkname, argv[1]);
        else {
            nomle = true;
            strcpy(mrkname, argv[2]);
        }

    // compute marker-specific IBDs
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
 
        Marker *marker;
        try {
            marker = currentPed->marker();
            if (!marker) {
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
 
        if (MMSibs && !currentMap) {
            RESULT_LIT ("Map data have not been loaded.");
            return TCL_ERROR;
        }

        char ibddir[1024];
        if (Solar_Eval(interp, "ibddir -session") == TCL_ERROR)
            return TCL_ERROR;
        else
            strcpy(ibddir, Tcl_GetStringResult (interp));
 
        if (MMSibs) {

            char ibd_cmd[1024];
            sprintf(ibd_cmd,
                    "exec dommsibs n %s", currentMap->filename());

            if (Solar_Eval(interp, ibd_cmd) == TCL_ERROR) {
                char errfile[1024], errmsg[1024];
                sprintf(errfile, "dommsibs.err");
                FILE *errfp = fopen(errfile, "r");
                if (errfp && fgets(errmsg, sizeof(errmsg), errfp)) {
                    RESULT_BUF (strtok(errmsg, "\n"));
                    fclose(errfp);
                }
                else
                    RESULT_LIT ("dommsibs failed");
                return TCL_ERROR;
            }

            return TCL_OK;
        }

        ibd_func = MCarlo ? run_mc : run_ibd;
        int i, retval = TCL_OK;

        if (doall) {
            for (i = 0; i < marker->nloci(); i++) {
                if (marker->ntyped(i) == currentPed->nind()
                    && currentFreq->xlinked(i) == 'n')
                {
                    if (run_mc(marker->mrkname(i),
                               true, ibddir, interp) == TCL_ERROR) {
                        printf("%s\n", Tcl_GetStringResult (interp)); 
                        fflush(stdout);
                        Tcl_ResetResult (interp);
                        retval = TCL_ERROR;
                    }
                }
                else if ((*ibd_func)(marker->mrkname(i),
                                     nomle, ibddir, interp) == TCL_ERROR) {
                    printf("%s\n", Tcl_GetStringResult (interp)); 
                    fflush(stdout);
                    Tcl_ResetResult (interp);
                    retval = TCL_ERROR;
                }
            }
        }
        else {
            if ((i = marker->get_marker(mrkname)) >= 0
                && marker->ntyped(i) == currentPed->nind()
                && currentFreq->xlinked(i) == 'n')
            {
                retval = run_mc(mrkname, true, ibddir, interp);
            }
            else {
                retval = (*ibd_func)(mrkname, nomle, ibddir, interp);
            }
        }

        return retval;
    }

    else if (argc == 4 && !strcmp("-inform", argv[1])) {
        char mrkfile[1024], ibdfile[1024];
        strcpy(mrkfile, argv[2]);
        strcpy(ibdfile, argv[3]);

    // compute IBDs for a fully-informative marker
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

        return inf_ibd(mrkfile, ibdfile, interp);
    }

    RESULT_LIT ("Usage: ibd [-nomle] [<marker>]");
    return TCL_ERROR;
}

int inf_ibd (const char *mrkfile, const char *ibdfile, Tcl_Interp *interp)
{
    char errmsg[1024];

    FILE *mrkfp = fopen(mrkfile, "r");
    if (!mrkfp) {
        sprintf(errmsg, "Cannot open %s", mrkfile);
        RESULT_BUF (errmsg);
        return TCL_ERROR;
    }

    FILE *ibdfp = fopen(ibdfile, "w");
    if (!ibdfp) {
        sprintf(errmsg, "Cannot open %s", ibdfile);
        RESULT_BUF (errmsg);
        return TCL_ERROR;
    }

    int *all1 = new int[currentPed->nind()];
    int *all2 = new int[currentPed->nind()];

    float ibd, delta7;
    char rec[1024];
    int id, a1, a2;

    int i, j;
    for (i = 0; i < currentPed->nind(); i++) {
        if (!fgets(rec, sizeof(rec), mrkfp) ||
            sscanf(rec, "%d %d/%d", &id, &a1, &a2) != 3 || id != i+1)
        {
            delete[] all1;
            delete[] all2;
            fclose(ibdfp);
            fclose(mrkfp);

            sprintf(errmsg, "Error reading %s", mrkfile);
            RESULT_BUF (errmsg);
            return TCL_ERROR;
        }

        all1[i] = a1;
        all2[i] = a2;
        for (j = 0; j <= i; j++) {
            ibd = 0;
            delta7 = 0;
            if (all1[i] == all1[j]) ibd += 0.5;
            if (all1[i] == all2[j]) ibd += 0.5;
            if (all2[i] == all1[j]) ibd += 0.5;
            if (all2[i] == all2[j]) ibd += 0.5;
            if (ibd >= 1) delta7 = 1;
            if (ibd)
                fprintf(ibdfp, "%5d %5d %10.7f %10.7f\n", i+1, j+1, ibd,
                        delta7);
        }
    }

    delete[] all1;
    delete[] all2;

    fclose(ibdfp);
    fclose(mrkfp);

    char cmd[1024];
    sprintf(cmd, "exec gzip -f %s", ibdfile);
    if (Solar_Eval(interp, cmd) == TCL_ERROR) {
        RESULT_LIT ("gzip failed");
        return TCL_ERROR;
    }

    sprintf (cmd, "matcrc %s.gz", ibdfile);
    if (Solar_Eval(interp, cmd) == TCL_ERROR) {
	return TCL_ERROR;
    }

    return TCL_OK;
}

int mito_ibd (const char *ibddir, Tcl_Interp *interp)
{
    printf("Computing mitochondrial IBDs ... ");
    fflush(stdout);

    FILE *fp = fopen("pedindex.out", "r");
    if (!fp) {
        RESULT_LIT ("\nCannot open pedindex.out");
        return TCL_ERROR;
    }

    float *mito = new float[currentPed->nind()];

    int id, fa, mo;
    char rec[1024];

    while (fgets(rec, sizeof(rec), fp)) {
        if (sscanf(rec, "%d %d %d", &id, &fa, &mo) != 3 ||
            id < 1 || id > currentPed->nind() ||
            mo < 0 || mo > currentPed->nind())
        {
            delete[] mito;
            RESULT_LIT ("\nError reading pedindex.out");
            return TCL_ERROR;
        }
        if (mo)
            mito[id-1] = mito[mo-1];
        else
            mito[id-1] = id;
    }

    fclose(fp);

    char fname[1024], errmsg[1024];
    sprintf(fname, "%s/ibd.mito", ibddir);

    fp = fopen(fname, "w");
    if (!fp ) {
        delete[] mito;
        sprintf(errmsg, "\nCannot open %s", fname);
        RESULT_BUF (errmsg);
        return TCL_ERROR;
    }

    for (int i = 0; i < currentPed->nind(); i++) {
        for (int j = 0; j <= i; j++) {
            if (mito[j] == mito[i])
                fprintf(fp, "%5d %5d  1.0000000  1.0000000\n", i+1, j+1);
        }
    }

    fclose(fp);
    delete[] mito;

    char cmd[1024];
    sprintf(cmd, "exec gzip -f %s", fname);
    if (Solar_Eval(interp, cmd) == TCL_ERROR) {
        RESULT_LIT ("\ngzip failed");
        return TCL_ERROR;
    }

    sprintf (cmd, "matcrc %s.gz", fname);
    if (Solar_Eval(interp, cmd) == TCL_ERROR) {
	return TCL_ERROR;
    }

    printf("\n");
    fflush(stdout);

    return TCL_OK;
}

int run_ibd (const char *mrkname, bool nomle, const char *ibddir, 
             Tcl_Interp *interp)
{
    char ibd_cmd[1024], errmsg[1024];

    int mrk = currentFreq->get_marker(mrkname);
    if (mrk < 0) {
        sprintf(errmsg, "%s: No such marker.", mrkname);
        RESULT_BUF (errmsg);
        return TCL_ERROR;
    }

    if (!nomle && currentFreq->mle_status(mrk) == 'n'
               && currentFreq->whence(mrk) == 'm')
    {
        sprintf(errmsg,
"The allele freqs for %s are not MLEs. Enter 'ibd -nomle' to use them.",
                currentFreq->mrkname(mrk));
        RESULT_BUF (errmsg);
        return TCL_ERROR;
    }

    if (!nomle && currentFreq->mle_status(mrk) == 'o'
               && currentFreq->whence(mrk) == 'f')
    {
        sprintf(errmsg,
"The allele freqs for %s are old MLEs. Enter 'ibd -nomle' to use them.",
                currentFreq->mrkname(mrk));
        RESULT_BUF (errmsg);
        return TCL_ERROR;
    }

    printf("Computing IBDs for %s ... ", currentFreq->mrkname(mrk));
    fflush(stdout);

    sprintf(ibd_cmd, "cd d_%s", currentFreq->mrkname(mrk));
    if (Solar_Eval(interp, ibd_cmd) == TCL_ERROR) {
        sprintf(errmsg, "\nCannot change to directory d_%s.",
                currentFreq->mrkname(mrk));
        RESULT_BUF (errmsg);
        return TCL_ERROR;
    }

    char fname[1024];
    if (ibddir[0] == '/')
        sprintf(fname, "%s/ibd.%s", ibddir, currentFreq->mrkname(mrk));
    else
        sprintf(fname, "../%s/ibd.%s", ibddir, currentFreq->mrkname(mrk));

    char show_status = 'n';
    FILE *fp = fopen("/dev/tty", "w");
    if (fp) {
        show_status = 'y';
        fclose(fp);
    }

    sprintf(ibd_cmd, "exec ibdmat %c %c %s >& ibdmat.out",
            currentFreq->xlinked(mrk), show_status, fname);
    if (Solar_Eval(interp, ibd_cmd) == TCL_ERROR) {
        sprintf(ibd_cmd, "cd ..");
        if (Solar_Eval(interp, ibd_cmd) == TCL_ERROR) {
            RESULT_LIT ("\nCannot return to current working directory.");
            return TCL_ERROR;
        }

        char errfile[1024], errbuf[1024];
        sprintf(errfile, "d_%s/dolink.err", currentFreq->mrkname(mrk));
        FILE *errfp = fopen(errfile, "r");
        if (errfp && fgets(errbuf, sizeof(errbuf), errfp)) {
            sprintf(errmsg, "\n%s", strtok(errbuf, "\n"));
            RESULT_BUF (errmsg);
            fclose(errfp);
        }
        else {
            sprintf(errfile, "d_%s/ibdmat.out", currentFreq->mrkname(mrk));
            errfp = fopen(errfile, "r");
            if (errfp && fgets(errbuf, sizeof(errbuf), errfp)) {
                sprintf(errmsg, "\n%s", strtok(errbuf, "\n"));
                RESULT_BUF (errmsg);
                fclose(errfp);
            }
            else
                RESULT_LIT ("\nProgram ibdmat did not run.");
        }
        return TCL_ERROR;
    }

    sprintf(ibd_cmd, "cd ..");
    if (Solar_Eval(interp, ibd_cmd) == TCL_ERROR) {
        RESULT_LIT ("\nCannot return to current working directory.");
        return TCL_ERROR;
    }

    sprintf (ibd_cmd, "matcrc %s.gz", fname);
    if (Solar_Eval(interp, ibd_cmd) == TCL_ERROR) {
	return TCL_ERROR;
    }

    printf("\n");
    fflush(stdout);

    return TCL_OK;
}

int run_mc (const char *mrkname, bool nomle, const char *ibddir, 
            Tcl_Interp *interp)
{
    char ibd_cmd[1024], errmsg[1024];

    int mrk = currentFreq->get_marker(mrkname);
    if (mrk < 0) { 
        sprintf(errmsg, "%s: No such marker.", mrkname);
        RESULT_BUF (errmsg);
        return TCL_ERROR;
    }

    if (currentFreq->nall(mrk) > MXMCALL) {
        sprintf(errmsg, "Cannot compute IBDs for %s, more than %d alleles.",
                currentFreq->mrkname(mrk), MXMCALL);
        RESULT_BUF (errmsg);
        return TCL_ERROR;
    }
 
    if (currentFreq->xlinked(mrk) == 'y') {
        RESULT_LIT (\
  "SOLAR does not yet support Monte Carlo IBD calculation for X-linked loci.");
        return TCL_ERROR;
    }

    if (!nomle && currentFreq->mle_status(mrk) == 'n'
               && currentFreq->whence(mrk) == 'm')
    {
        sprintf(errmsg,
"The allele freqs for %s are not MLEs. Enter 'ibd -nomle' to use them.",
                currentFreq->mrkname(mrk));
        RESULT_BUF (errmsg);
        return TCL_ERROR;
    }

    if (!nomle && currentFreq->mle_status(mrk) == 'o'
               && currentFreq->whence(mrk) == 'f')
    {
        sprintf(errmsg,
"The allele freqs for %s are old MLEs. Enter 'ibd -nomle' to use them.",
                currentFreq->mrkname(mrk));
        RESULT_BUF (errmsg);
        return TCL_ERROR;
    }

    printf("Computing IBDs for %s ... ", currentFreq->mrkname(mrk));
    fflush(stdout);

    sprintf(ibd_cmd, "cd d_%s", currentFreq->mrkname(mrk));
    if (Solar_Eval(interp, ibd_cmd) == TCL_ERROR) {
        sprintf(errmsg, "\nCannot change to directory d_%s.",
                currentFreq->mrkname(mrk));
        RESULT_BUF (errmsg);
        return TCL_ERROR;
    }

    char fname[1024];
    if (ibddir[0] == '/')
        sprintf(fname, "%s/ibd.%s", ibddir, currentFreq->mrkname(mrk));
    else
        sprintf(fname, "../%s/ibd.%s", ibddir, currentFreq->mrkname(mrk));

    char maxrisk = MaxRisk ? 'y' : 'n';
    char show_status = 'n';
    FILE *fp = fopen("/dev/tty", "w");
    if (fp) {
        show_status = 'y';
        fclose(fp);
    }

    sprintf(ibd_cmd, "exec domcibd %s %d %c %c", fname, NumImp, maxrisk,
            show_status);
    if (Solar_Eval(interp, ibd_cmd) == TCL_ERROR) {
        sprintf(ibd_cmd, "cd ..");
        if (Solar_Eval(interp, ibd_cmd) == TCL_ERROR) {
            RESULT_LIT ("\nCannot return to current working directory.");
            return TCL_ERROR;
        }

        char errfile[1024], errbuf[1024];
        sprintf(errfile, "d_%s/ibdmc.err", currentFreq->mrkname(mrk));
        FILE *errfp = fopen(errfile, "r");
        if (errfp && fgets(errbuf, sizeof(errbuf), errfp)) {
            sprintf(errmsg, "\n%s", strtok(errbuf, "\n"));
            RESULT_BUF (errmsg);
            fclose(errfp);
        }
        else
            RESULT_LIT ("\nProgram domcibd did not run.");
        return TCL_ERROR;
    }

    sprintf(ibd_cmd, "cd ..");
    if (Solar_Eval(interp, ibd_cmd) == TCL_ERROR) {
        RESULT_LIT ("\nCannot return to current working directory.");
        return TCL_ERROR;
    }

    sprintf (ibd_cmd, "matcrc %s.gz", fname);
    if (Solar_Eval(interp, ibd_cmd) == TCL_ERROR) {
	return TCL_ERROR;
    }

    printf("\n");
    fflush(stdout);

    return TCL_OK;
}
