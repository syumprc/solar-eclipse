/*
 * ibs.cc implements the ibs command
 * Written by Thomas Dyer July 1998
 * Copyright (c) 1998 Southwest Foundation for Biomedical Research
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "solar.h"
#include "tcl.h"
#include "safelib.h"
#include "pipeback.h"

static int run_ibs (const char*, const char*, Tcl_Interp*);

extern "C" int IbsCmd (ClientData clientData, Tcl_Interp *interp, int argc,
            char *argv[])
{
    if (argc == 2 && !StringCmp ("help", argv[1], case_ins)) {
        return Solar_Eval (interp, "help ibd");
    }
 
    else {
    // compute IBSs
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

	if (XLinked) {
            RESULT_LIT ("Not yet implemented for X-linked markers.");
            return TCL_ERROR;
        }

        char ibddir[1024];
        if (Solar_Eval(interp, "ibddir -session") == TCL_ERROR)
            return TCL_ERROR;
        else
            strcpy(ibddir, Tcl_GetStringResult (interp));
 
        int i, retval = TCL_OK;
        int nloci = currentFreq->nloci();

        if (argc == 1) {
            for (i = 0; i < nloci; i++) {
                if (run_ibs(currentFreq->mrkname(i), ibddir, interp) == TCL_ERROR)
                {
                    printf("%s\n", interp->result);
                    fflush(stdout);
                    retval = TCL_ERROR;
                }
            }
            if (retval == TCL_ERROR)
                RESULT_LIT ("An error has occurred for one or more markers.");
        }
        else {
            for (i = 0; i < nloci; i++) {
                if (!StringCmp(argv[1], currentFreq->mrkname(i), case_ins))
                    break;
            }
            if (i == nloci) {
                char errmsg[1000];
                sprintf(errmsg, "%s: No such marker.", argv[1]);
                RESULT_BUF (errmsg);
                retval = TCL_ERROR;
            }
            else
                retval = run_ibs(argv[1], ibddir, interp);
        }
 
        return retval;
    }

    RESULT_LIT ("Invalid ibs command");
    return TCL_ERROR;
}

int run_ibs (const char *mrkname, const char *ibddir, Tcl_Interp *interp)
{
    int mrk = currentFreq->get_marker(mrkname);
    if (mrk < 0) {
        RESULT_LIT ("No such marker.");
        return TCL_ERROR;
    }

    printf("Computing IBSs for %s ... ", currentFreq->mrkname(mrk));
    fflush(stdout);

    int nids = currentPed->nind();
    int npairs = nids * nids;

    unsigned char *related;
    related = (unsigned char *) malloc(npairs*sizeof(unsigned char));
    if (!related) {
        RESULT_LIT ("Out of memory.");
        return TCL_ERROR;
    }
    memset(related, 0, npairs);
   
    FILE *phi2fp = fopen("phi2.gz", "r");
    if (!phi2fp) {
        RESULT_LIT ("Cannot open phi2.gz");
        return TCL_ERROR;
    }
    fclose(phi2fp);

    const char *pbarg[4];
    pbarg[0] = "gunzip";
    pbarg[1] = "-c";
    pbarg[2] = "phi2";
    pbarg[3] = 0;
    phi2fp = pipeback_shell_open ("gunzip", pbarg);
    if (!phi2fp) {
        RESULT_LIT ("Unable to uncompress phi2");
        return TCL_ERROR;
    }

    int id1, id2;
    double phi2, delta7;
    char rec[1000];
    while (fgets(rec, sizeof(rec), phi2fp)) {
        if (sscanf(rec, "%d %d %lf %lf", &id1, &id2, &phi2, &delta7) != 4)
        {
            free(related);
            RESULT_LIT ("Error reading phi2");
            return TCL_ERROR;
        }
        if (phi2 != 0.)
            *(related + (id1-1)*nids + (id2-1)) = (unsigned char) 1;
        else
            *(related + (id1-1)*nids + (id2-1)) = (unsigned char) 0;
    }
    pipeback_shell_close(phi2fp);

    int *all1 = (int *) malloc(currentPed->nind()*sizeof(int));
    if (!all1) {
        RESULT_LIT ("Out of memory.");
        return TCL_ERROR;
    }

    int *all2 = (int *) malloc(currentPed->nind()*sizeof(int));
    if (!all2) {
        free(all1); free(related);
        RESULT_LIT ("Out of memory.");
        return TCL_ERROR;
    }

    char filename[100], errmsg[1000];
    sprintf(filename, "d_%s/translat.tab", currentFreq->mrkname(mrk));
    FILE *pedfp = fopen(filename, "r");
    if (!pedfp) {
        sprintf(errmsg, "Cannot open %s", filename);
        RESULT_BUF (errmsg);
        return TCL_ERROR;
    }

    if (!fgets(rec, sizeof(rec), pedfp) || !fgets(rec, sizeof(rec), pedfp))
    {
        free(all1); free(all2); free(related);
        sprintf(errmsg, "Error reading %s", filename);
        RESULT_BUF (errmsg);
        return TCL_ERROR;
    }

    int id, i = 0;
    while (fgets(rec, sizeof(rec), pedfp)) {
        if (sscanf(rec, "%d", &nids) != 1)
        {
            free(all1); free(all2); free(related);
            sprintf(errmsg, "Error reading %s", filename);
            RESULT_BUF (errmsg);
            return TCL_ERROR;
        }

        while (nids) {
            if (!fgets(rec, sizeof(rec), pedfp) || sscanf(rec, "%d", &id) != 1)
            {
                free(all1); free(all2); free(related);
                sprintf(errmsg, "Error reading %s", filename);
                RESULT_BUF (errmsg);
                return TCL_ERROR;
            }
            if (!strncmp(rec+19, "      ", 6))
                all1[i] = all2[i] = 0;
            else if (sscanf(rec+19, "%d %d", all1+i, all2+i) != 2)
            {
                free(all1); free(all2); free(related);
                sprintf(errmsg, "Error reading %s", filename);
                RESULT_BUF (errmsg);
                return TCL_ERROR;
            }
            i++;
            nids--;
        }
    }

    fclose(pedfp);

    if (i != currentPed->nind())
    {
        free(all1); free(all2); free(related);
        sprintf(errmsg, "Error reading %s", filename);
        RESULT_BUF (errmsg);
        return TCL_ERROR;
    }

    sprintf(filename, "%s/ibs.%s", ibddir, currentFreq->mrkname(mrk));
    FILE *ibsfp = fopen(filename, "w");
    if (!ibsfp) {
        free(all1); free(all2); free(related);
        sprintf(errmsg, "Cannot open %s", filename);
        RESULT_BUF (errmsg);
        return TCL_ERROR;
    }

    int j, is, js;
    double ibs;
    nids = currentPed->nind();
    for (i = 0; i < currentPed->nind(); i++) {
        if (all1[i] != 0 && all2[i] != 0) {
            for (j = 0; j < i; j++) {
                if (all1[j] == 0 && all2[j] == 0) continue;
                if (*(related + i*nids + j)) continue;
                is = 0; js = 0;
                if (all1[i] == all1[j] || all1[i] == all2[j]) is++;
                if (all2[i] == all1[j] || all2[i] == all2[j]) is++;
                if (all1[j] == all1[i] || all1[j] == all2[i]) js++;
                if (all2[j] == all1[i] || all2[j] == all2[i]) js++;
                if (is == 2 && js == 2) ibs = 1;
                else if (is == 0 && js == 0) ibs = 0;
                else ibs = .5;
                if (ibs)
                    fprintf(ibsfp, "%5d %5d %10.7f %10.7f\n", i+1, j+1,
                            ibs, 0.);
            }
        }
        fprintf(ibsfp, "%5d %5d %10.7f %10.7f\n", i+1, i+1, 0., 0.);
    }

    fclose(ibsfp);
    free(all1);
    free(all2);
    free(related);

    char cmd[1024];
    sprintf(cmd, "exec gzip -f %s/ibs.%s", ibddir, currentFreq->mrkname(mrk));
    if (Tcl_Eval(interp, cmd) == TCL_ERROR) {
        RESULT_LIT ("gzip failed");
        return TCL_ERROR;
    }

    sprintf (cmd, "matcrc %s/ibs.%s.gz", ibddir, currentFreq->mrkname(mrk));
    if (Solar_Eval(interp, cmd) == TCL_ERROR) {
        return TCL_ERROR;
    }

    printf("\n");
    fflush(stdout);

    return TCL_OK;
}
