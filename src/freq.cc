/*
 * freq.cc implements the freq command and class
 * Written by Thomas Dyer January 1998
 * Copyright (c) 1998 Southwest Foundation for Biomedical Research
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "solar.h"
// tcl.h from solar.h
#include "safelib.h"

Freq *currentFreq = 0;

extern "C" int FreqCmd (ClientData clientData, Tcl_Interp *interp,
                  int argc, char *argv[])
{
    if (argc == 2 && !StringCmp ("help", argv[1], case_ins)) {
        return Solar_Eval (interp, "help freq");
    }

    else if (argc >= 2 && !StringCmp ("load", argv[1], case_ins)) {
        if (argc != 3 && argc != 4) {
            RESULT_LIT ("Usage: freq load [-nosave] <filename>");
            return TCL_ERROR;
        }

        bool nosave;
        char fname[1024];

        if (argc == 4) { 
            if (strcmp("-nosave", argv[2])) {
                RESULT_LIT ("Usage: freq load [-nosave] <filename>");
                return TCL_ERROR;
            }
            else {
                nosave = true;
                strcpy(fname, argv[3]);  
            }
        }
        else {
            if (!strcmp("-nosave", argv[2])) {
                RESULT_LIT ("Usage: freq load [-nosave] <filename>");
                return TCL_ERROR;
            }
            else {
                nosave = false;
                strcpy(fname, argv[2]);  
            }
        }

    // load new frequency data
        char xlinked = XLinked ? 'y' : 'n';

        bool freq_loaded;
        try {
            freq_loaded = loadedFreq();
        }
        catch (Safe_Error_Return& ser) {
            RESULT_BUF (ser.message());
            return TCL_ERROR;
        }

        if (freq_loaded) {

//  unless nosave option is in effect, abort the load if any of
//  the markers have unsaved MLE allele freqs
 
            bool ped_loaded;
            try {
                ped_loaded = loadedPed();
            }
            catch (Safe_Error_Return& ser) {
                RESULT_BUF (ser.message());
                Solar_AppendResult1(interp,
                            "\nReloading the pedigree data is advised.", NULL);
                return TCL_ERROR;
            }

            Marker *marker;
            if (ped_loaded) {
                try {
                    marker = currentPed->marker();
                }
                catch (Safe_Error_Return& ser) {
                    RESULT_BUF (ser.message());
                    Solar_AppendResult1(interp,
                            "\nReloading the marker data is advised.", NULL);
                    return TCL_ERROR;
                }
            }

            if (ped_loaded && marker) {
                if (currentFreq->nloci())
                    xlinked = currentFreq->xlinked(0);

                if (!nosave) {
                    int i;
                    for (i = 0; i < currentFreq->nloci(); i++) {
                        if (marker->get_marker(currentFreq->mrkname(i)) >= 0
                                && currentFreq->mle_status(i) == 'c')
                        {
                                RESULT_LIT (\
"Use freq save to save MLE allele frequencies, then do freq load again.\n\
Or use freq load -nosave to load allele freqs without saving MLEs.");
                            return TCL_ERROR;
                        }
                    }
                }
            }

            if (currentFreq->unload(interp)) {
                delete currentFreq;
                currentFreq = new Freq (fname);
            }
            else
                currentFreq->change_filename(fname);
        }
        else
            currentFreq = new Freq (fname);

        if (currentFreq->load(xlinked, interp) == TCL_ERROR) {
            delete currentFreq;
            try {
                if (!loadedFreq()) {
                    return TCL_ERROR;
                }
            }
            catch (Safe_Error_Return& ser) {
                RESULT_BUF (ser.message());
                return TCL_ERROR;
            }
            if (currentFreq->unload(interp))
                delete currentFreq;
            return TCL_ERROR;
        }

        return TCL_OK;
    }

    else if (argc == 2 && !StringCmp ("unload", argv[1], case_ins)) {
    // unload frequency data (except for loci from a current marker load)
 
        try {
            bool ped_loaded = loadedPed();
        }
        catch (Safe_Error_Return& ser) {
            RESULT_BUF (ser.message());
            Solar_AppendResult1(interp,
                        "\nReloading the pedigree data is advised.", NULL);
            return TCL_ERROR;
        }

        try {
            if (!loadedFreq()) {
                RESULT_LIT ("Frequency data have not been loaded.");
                return TCL_ERROR;
            }
        }
        catch (Safe_Error_Return& ser) {
            RESULT_BUF (ser.message());
            return TCL_ERROR;
        }

        if (currentFreq->unload(interp))
            delete currentFreq;

        return TCL_OK;
    }

    else if (argc >= 2 && !StringCmp ("mle", argv[1], case_ins)) {
    // compute MLE allele frequencies

        try {
            if (!loadedPed()) {
                RESULT_LIT ("Pedigree and marker data have not been loaded.");
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

        try {
            if (!loadedFreq()) {
                RESULT_LIT ("Frequency data have not been loaded.");
                return TCL_ERROR;
            }
        }
        catch (Safe_Error_Return& ser) {
            RESULT_BUF (ser.message());
            return TCL_ERROR;
        }

        if (argc == 5 && !StringCmp ("-nose", argv[2], case_ins)
                      && !StringCmp ("-hwe", argv[3], case_ins))
        {
            if (currentFreq->mlefreq(argv[4], false, true, interp) == TCL_ERROR)
                return TCL_ERROR;
            if (currentFreq->write_info(interp) == TCL_ERROR)
                return TCL_ERROR;
        }
        else if (argc == 4 && !StringCmp ("-nose", argv[2], case_ins))
        {
            if (currentFreq->mlefreq(argv[3], false, false, interp) == TCL_ERROR)
                return TCL_ERROR;
            if (currentFreq->write_info(interp) == TCL_ERROR)
                return TCL_ERROR;
        }
        else if (argc == 4 && !StringCmp ("-hwe", argv[2], case_ins))
        {
            if (currentFreq->mlefreq(argv[3], true, true, interp) == TCL_ERROR)
                return TCL_ERROR;
            if (currentFreq->write_info(interp) == TCL_ERROR)
                return TCL_ERROR;
        }
        else if (argc == 3) {
            if (currentFreq->mlefreq(argv[2], true, false, interp) == TCL_ERROR)
                return TCL_ERROR;
            if (currentFreq->write_info(interp) == TCL_ERROR)
                return TCL_ERROR;
        }
        else {
            RESULT_LIT ("Usage: freq mle [-nose] [-hwe] <marker>");
            return TCL_ERROR;
        }

        return TCL_OK;
    }

    else if (argc >= 2 && !StringCmp ("save", argv[1], case_ins)) {
        if (argc != 3) {
            RESULT_LIT ("Usage: freq save <filename>");
            return TCL_ERROR;
        }

    // save frequency data to file
        try {
            if (!loadedFreq()) {
                RESULT_LIT ("Frequency data have not been loaded.");
                return TCL_ERROR;
            }
        }
        catch (Safe_Error_Return& ser) {
            RESULT_BUF (ser.message());
            return TCL_ERROR;
        }

        if (currentFreq->save(argv[2], interp) == TCL_ERROR)
            return TCL_ERROR;
        if (currentFreq->write_info(interp) == TCL_ERROR)
            return TCL_ERROR;

        return TCL_OK;
    }

    else if (argc >= 2 && !StringCmp ("show", argv[1], case_ins)) {
        if (argc > 3) {
            RESULT_LIT ("Usage: freq show [<marker>]");
            return TCL_ERROR;
        }

    // display frequency data
        try {
            if (!loadedFreq()) {
                RESULT_LIT ("Frequency data have not been loaded.");
                return TCL_ERROR;
            }
        }
        catch (Safe_Error_Return& ser) {
            RESULT_BUF (ser.message());
            return TCL_ERROR;
        }
 
        int i, bsize = 80;
        for (i = 0; i < currentFreq->nloci(); i++)
            bsize += 80 + (currentFreq->nall(i)/4 + 1)*80;
        char *buf = (char *) malloc(bsize);

        if (argc == 2)
            strcpy(buf, "");
        else
            strcpy(buf, argv[2]);
 
        currentFreq->show(buf);
        RESULT_BUF (buf);
        free(buf);
        return TCL_OK;
    }

    else if (argc == 3 && !StringCmp ("nall", argv[1], case_ins)) {

    // display number of alleles
        try {
            if (!loadedFreq()) {
                RESULT_LIT ("Frequency data have not been loaded.");
                return TCL_ERROR;
            }
        }
        catch (Safe_Error_Return& ser) {
            RESULT_BUF (ser.message());
            return TCL_ERROR;
        }

        int mrk = currentFreq->get_marker(argv[2]);
        if (mrk < 0) {
            RESULT_LIT ("No such marker.");
            return TCL_ERROR;
        }

        int nall = currentFreq->nall(mrk);
        char buf[100];
        sprintf(buf, "%d", nall);
        RESULT_BUF (buf);
        return TCL_OK;
    }

    else if (argc == 3 && !StringCmp ("chi2_hwe", argv[1], case_ins)) {

    // display chi-square statistic for test of HWE
        try {
            if (!loadedFreq()) {
                RESULT_LIT ("Frequency data have not been loaded.");
                return TCL_ERROR;
            }
        }
        catch (Safe_Error_Return& ser) {
            RESULT_BUF (ser.message());
            return TCL_ERROR;
        }

        int mrk = currentFreq->get_marker(argv[2]);
        if (mrk < 0) {
            RESULT_LIT ("No such marker.");
            return TCL_ERROR;
        }

        double chi2 = currentFreq->chi2_hwe(mrk);
        char buf[100];
        sprintf(buf, "%f", chi2);
        RESULT_BUF (buf);
        return TCL_OK;
    }

    else if (argc == 4 && !StringCmp ("name", argv[1], case_ins)) {

    // display allele name
        try {
            if (!loadedFreq()) {
                RESULT_LIT ("Frequency data have not been loaded.");
                return TCL_ERROR;
            }
        }
        catch (Safe_Error_Return& ser) {
            RESULT_BUF (ser.message());
            return TCL_ERROR;
        }

        int mrk = currentFreq->get_marker(argv[2]);
        if (mrk < 0) {
            RESULT_LIT ("No such marker.");
            return TCL_ERROR;
        }

        int all;
        if (sscanf(argv[3], "%d", &all) != 1 ||
                all < 0 || all >= currentFreq->nall(mrk))
        {
            RESULT_LIT ("Illegal allele index.");
            return TCL_ERROR;
        }

        char buf[100];
        strcpy(buf, currentFreq->name(mrk, all));
        RESULT_BUF (buf);
        return TCL_OK;
    }

    else if (argc == 4 && !StringCmp ("freq", argv[1], case_ins)) {

    // display allele frequency
        try {
            if (!loadedFreq()) {
                RESULT_LIT ("Frequency data have not been loaded.");
                return TCL_ERROR;
            }
        }
        catch (Safe_Error_Return& ser) {
            RESULT_BUF (ser.message());
            return TCL_ERROR;
        }

        int mrk = currentFreq->get_marker(argv[2]);
        if (mrk < 0) {
            RESULT_LIT ("No such marker.");
            return TCL_ERROR;
        }

        int all;
        if (sscanf(argv[3], "%d", &all) != 1 ||
                all < 0 || all >= currentFreq->nall(mrk))
        {
            RESULT_LIT ("Illegal allele index.");
            return TCL_ERROR;
        }

        double freq = currentFreq->freq(mrk, all);
        char buf[100];
        sprintf(buf, "%f", freq);
        RESULT_BUF (buf);
        return TCL_OK;
    }

    else if (argc == 4 && !StringCmp ("mle_se", argv[1], case_ins)) {

    // display standard error for allele frequency MLE
        try {
            if (!loadedFreq()) {
                RESULT_LIT ("Frequency data have not been loaded.");
                return TCL_ERROR;
            }
        }
        catch (Safe_Error_Return& ser) {
            RESULT_BUF (ser.message());
            return TCL_ERROR;
        }

        int mrk = currentFreq->get_marker(argv[2]);
        if (mrk < 0) {
            RESULT_LIT ("No such marker.");
            return TCL_ERROR;
        }

        int all;
        if (sscanf(argv[3], "%d", &all) != 1 ||
                all < 0 || all >= currentFreq->nall(mrk))
        {
            RESULT_LIT ("Illegal allele index.");
            return TCL_ERROR;
        }

        double se = currentFreq->mle_se(mrk, all);
        char buf[100];
        sprintf(buf, "%f", se);
        RESULT_BUF (buf);
        return TCL_OK;
    }

    RESULT_LIT ("Invalid freq command");
    return TCL_ERROR;
}

Freq::~Freq ()
{
    int i, j;

    for (i = 0; i < _nloci; i++) {
        for (j = 0; j < _locus[i]->nall; j++)
            free(_locus[i]->all[j]);
        free(_locus[i]);
    }

    currentFreq = 0;
}
 
bool loadedFreq (void)
{
    if (currentFreq)
        return true;
 
    FILE *fp = fopen("freq.info", "r");
    if (!fp)
        return false;

    char fname[1024];
    fscanf(fp, "%s", fname);
    fclose(fp);

    currentFreq = new Freq (fname);

    const char *errmsg = 0;
    if (!currentFreq->read_info(&errmsg)) {
        delete currentFreq;
        throw Safe_Error_Return(errmsg);
    }

    return true;
}

bool Freq::read_info (const char **errmsg)
{
    char rec[10000], *recp;
    char name[1024];

    FILE *fp = fopen("freq.info", "r");
    if (!fp) {
        *errmsg = "Cannot open freq.info";
        return false;
    }

    if (!fgets(rec, sizeof(rec), fp) || sscanf(rec, "%s", name) != 1) {
        *errmsg = "Read error on freq.info, line 1";
        return false;
    }

    _nloci = 0;
    while (fgets(rec, sizeof(rec), fp)) {
        if (!(recp = strtok(rec, " \t\n")) || sscanf(recp, "%s", name) != 1) {
            sprintf(error_message,
                    "Read error on freq.info, line %d", _nloci + 2);
            *errmsg = error_message;
            return false;
        }

        if (_nloci == MAXLOC) {
            sprintf(error_message,
                    "ERROR: More than the maximum of %d loci in freq.info",
                    MAXLOC);
            *errmsg = error_message;
            return false;
        }

        _locus[_nloci] = (struct LocusFreq *) malloc(sizeof(struct LocusFreq));
        strcpy(_locus[_nloci]->name, name);

        // increment now so all heap space is freed if read error occurs
        _nloci++;
        _locus[_nloci-1]->nall = 0;
        _locus[_nloci-1]->chi2_hwe = -1;

        if (!(recp = strtok(NULL, " \t\n")) ||
                sscanf(recp, "%c", &_locus[_nloci-1]->xlinked) != 1) {
            sprintf(error_message,
                    "Read error on freq.info, line %d", _nloci + 1);
            *errmsg = error_message;
            return false;
        }

        if (!(recp = strtok(NULL, " \t\n")) ||
                sscanf(recp, "%c", &_locus[_nloci-1]->mle_status) != 1) {
            sprintf(error_message,
                    "Read error on freq.info, line %d", _nloci + 1);
            *errmsg = error_message;
            return false;
        }

        if (!(recp = strtok(NULL, " \t\n")) ||
                sscanf(recp, "%c", &_locus[_nloci-1]->whence) != 1) {
            sprintf(error_message,
                    "Read error on freq.info, line %d", _nloci + 1);
            *errmsg = error_message;
            return false;
        }

        while (recp = strtok(NULL, " \t\n")) {
            if (sscanf(recp, "%s", name) != 1) {
                sprintf(error_message,
                        "Read error on freq.info, line %d", _nloci + 1);
                *errmsg = error_message;
                return false;
            }

            if (!strcmp(name, "se=")) {
                if (!(recp = strtok(NULL, " \t\n")) ||
                    sscanf(recp, "%lf",
                           &_locus[_nloci-1]->all[_locus[_nloci-1]->nall-1]->mle_se) != 1)
                {
                    sprintf(error_message,
                            "Read error on freq.info, line %d", _nloci + 1);
                    *errmsg = error_message;
                    return false;
                }
            }

            else if (!strcmp(name, "hwe=")) {
                if (!(recp = strtok(NULL, " \t\n")) ||
                    sscanf(recp, "%lf", &_locus[_nloci-1]->chi2_hwe) != 1)
                {
                    sprintf(error_message,
                            "Read error on freq.info, line %d", _nloci + 1);
                    *errmsg = error_message;
                    return false;
                }
            }

            else {
                _locus[_nloci-1]->all[_locus[_nloci-1]->nall] =
                            (struct Allele *) malloc(sizeof(struct Allele));
                strcpy(_locus[_nloci-1]->all[_locus[_nloci-1]->nall]->name, name);

                // increment now so all heap space is freed if read error occurs
                _locus[_nloci-1]->nall++;
                _locus[_nloci-1]->all[_locus[_nloci-1]->nall-1]->mle_se = -1;

                if (!(recp = strtok(NULL, " \t\n")) || sscanf(recp, "%lf",
                    &_locus[_nloci-1]->all[_locus[_nloci-1]->nall-1]->freq) != 1)
                {
                    sprintf(error_message,
                            "Read error on freq.info, line %d", _nloci + 1);
                    *errmsg = error_message;
                    return false;
                }
            }
        }
    }

    return true;
}

int Freq::load (char xlinked, Tcl_Interp *interp)
{
    return get_freqs(_filename, xlinked, interp);
}

int Freq::get_freqs (const char *fname, char xlinked, Tcl_Interp *interp)
{
    char rec[10000], *recp;
    char errmsg[1024];
    char cmd[1024];
    double sumfreq;
    bool ok, errors, numeric, check_names;
    struct Allele *allp;
    int old_loc;
    int i, j, loc, lineno;
    char mle_status, whence;

    FILE *locfp = fopen(fname, "r");
    if (!locfp) {
        Solar_AppendResult2(interp, "Unable to open freq file: ",
                         (char*) fname, NULL);
        return TCL_ERROR;
    }

    bool ped_loaded;
    try {
        ped_loaded = loadedPed();
    }
    catch (Safe_Error_Return& ser) {
        RESULT_BUF (ser.message());
        Solar_AppendResult1(interp,
                    "\nReloading the pedigree data is advised.", NULL);
        return TCL_ERROR;
    }

    Marker *marker = NULL;
    if (ped_loaded) {
        try {
            marker = currentPed->marker();
        }
        catch (Safe_Error_Return& ser) {
            RESULT_BUF (ser.message());
            Solar_AppendResult1(interp,
                        "\nReloading the marker data is advised.", NULL);
            return TCL_ERROR;
        }
    }

    errors = false;
    lineno = 0;
    while (fgets(rec, sizeof(rec), locfp)) {
        lineno++;
        if (!(recp = strtok(rec, " \t\n"))) {
            fclose(locfp);
            sprintf(errmsg,
                    "Load failed: Invalid record, line %d of freq file",
                    lineno);
            RESULT_BUF (errmsg);
            return TCL_ERROR;
        }

        if (strlen(recp) > MMRKNM) {
            sprintf(errmsg,
        "ERROR: Marker name [%s] longer than the maximum of %d characters\n",
                    recp, MMRKNM);
            Solar_AppendResult1(interp, errmsg, NULL);
            errors = true;
            continue;
        }

        check_names = false;
        old_loc = get_marker(recp);
        if (old_loc >= 0) {
            if (!strcmp(fname, "ibdprep.loc")) {
                mle_status = _locus[old_loc]->mle_status;
                whence = _locus[old_loc]->whence;
            }
            else {
                mle_status = 'n';
                whence = 'f';
                check_names = true;
            }
        }
        else {
            if (_nloci == MAXLOC) {
                sprintf(errmsg,
                        "Load failed: More than the maximum of %d loci",
                        MAXLOC);
                RESULT_BUF (errmsg);
                return TCL_ERROR;
            }
            mle_status = 'n';
            if (!strcmp(fname, "ibdprep.loc"))
                whence = 'm';
            else
                whence = 'f';
        }

        loc = _nloci;
        _locus[loc] = (struct LocusFreq *) malloc(sizeof(struct LocusFreq));
        strcpy(_locus[loc]->name, recp);

        // increment now so all heap space is freed if read error occurs
        _nloci++;
 
        _locus[loc]->nall = 0;
        _locus[loc]->chi2_hwe = -1;
        numeric = true;
        sumfreq = 0;

        ok = true;
        while (recp = strtok(NULL, " \t\n")) {
            if (_locus[loc]->nall >= MAXALL) {
                sprintf(errmsg,
                "ERROR: Marker %s has more than the maximum of %d alleles\n",
                        _locus[loc]->name, MAXALL);
                Solar_AppendResult1(interp, errmsg, NULL);
                errors = true;
                break;
            }
 
            if (strlen(recp) > MALLNM) {
                sprintf(errmsg,
        "ERROR: Marker %s allele name [%s] longer than maximum of %d characters\n",
                        _locus[loc]->name, recp, MALLNM);
                Solar_AppendResult1(interp, errmsg, NULL);
                ok = false;
                recp = strtok(NULL, " \t\n");	/* skip allele freq */
                continue;
            }

            if (!strcmp(recp, "se=")) {
                if (!(recp = strtok(NULL, " \t\n")) || sscanf(recp, "%lf",
                    &_locus[loc]->all[_locus[loc]->nall-1]->mle_se) != 1)
                {
                    fclose(locfp);
                    sprintf(errmsg,
                            "Load failed: Invalid record, line %d of freq file",
                            lineno);
                    RESULT_BUF (errmsg);
                    return TCL_ERROR;
                }
                continue;
            }

            if (!strcmp(recp, "hwe=")) {
                if (!(recp = strtok(NULL, " \t\n")) ||
                    sscanf(recp, "%lf", &_locus[loc]->chi2_hwe) != 1)
                {
                    fclose(locfp);
                    sprintf(errmsg,
                            "Load failed: Invalid record, line %d of freq file",
                            lineno);
                    RESULT_BUF (errmsg);
                    return TCL_ERROR;
                }
                continue;
            }

            _locus[loc]->all[_locus[loc]->nall] =
                        (struct Allele *) malloc(sizeof(struct Allele));
            strcpy(_locus[loc]->all[_locus[loc]->nall]->name, recp);

            // increment now so all heap space is freed if read error occurs
            _locus[loc]->nall++;
            _locus[loc]->all[_locus[loc]->nall-1]->mle_se = -1;

            for (char *p = recp; *p; p++)
                if (!isdigit(*p))
                    numeric = false;
 
            if (!(recp = strtok(NULL, " \t\n")) || sscanf(recp, "%lf",
                &_locus[loc]->all[_locus[loc]->nall-1]->freq) != 1)
            {
                fclose(locfp);
                sprintf(errmsg,
                        "Load failed: Invalid record, line %d of freq file",
                        lineno);
                RESULT_BUF (errmsg);
                return TCL_ERROR;
            }

            if (_locus[loc]->all[_locus[loc]->nall-1]->freq <= 0) {
                sprintf(errmsg,
                        "ERROR: Marker %s has an allele with frequency <= 0\n",
                        _locus[loc]->name);
                Solar_AppendResult1(interp, errmsg, NULL);
                ok = false;
                continue;
            }

            sumfreq += _locus[loc]->all[_locus[loc]->nall-1]->freq;
        }

        if (!ok) {
            errors = true;
            continue;
        }

        if (_locus[loc]->nall &&
                (sumfreq <= 0.9999 || sumfreq >= 1.0001)) {
            sprintf(errmsg,
        "ERROR: Allele frequencies don't sum to 1 for marker %s, sum = %g\n",
                    _locus[loc]->name, sumfreq);
            Solar_AppendResult1(interp, errmsg, NULL);
            errors = true;
            continue;
        }

        for (i = 0; i < _locus[loc]->nall; i++) {
            for (j = i + 1; j < _locus[loc]->nall; j++) {
                if (numeric) {
                    if (atoi(_locus[loc]->all[j]->name)
                        < atoi(_locus[loc]->all[i]->name))
                    {
                        allp = _locus[loc]->all[i];
                        _locus[loc]->all[i] = _locus[loc]->all[j];
                        _locus[loc]->all[j] = allp;
                    }
                }
                else {
                    if (strcmp(_locus[loc]->all[j]->name,
                               _locus[loc]->all[i]->name) < 0)
                    {
                        allp = _locus[loc]->all[i];
                        _locus[loc]->all[i] = _locus[loc]->all[j];
                        _locus[loc]->all[j] = allp;
                    }
                }
            }
        }

        if (check_names) {
            bool found;
            for (i = 0; i < _locus[old_loc]->nall; i++) {
                found = false;
                for (j = 0; j < _locus[loc]->nall; j++) {
                    if (!strcmp(_locus[loc]->all[j]->name,
                                _locus[old_loc]->all[i]->name))
                    {
                        found = true;
                        break;
                    }
                }

                if (!found) {
                    sprintf(errmsg,
        "ERROR: Marker %s alleles don't include all those found in the genotype data\n",
                            _locus[loc]->name);
                    Solar_AppendResult1(interp, errmsg, NULL);
                    errors = true;
                }

            }
        }

        _locus[loc]->xlinked = xlinked;
        _locus[loc]->mle_status = mle_status;
        _locus[loc]->whence = whence;

        if (old_loc >= 0)
            remove_marker(old_loc);
    }

    fclose(locfp);
    if (errors) {
        Solar_AppendResult1(interp,
                            "Errors encountered. Frequency data not loaded.",
                            NULL);
        return TCL_ERROR;
    }

    if (!lineno) {
        RESULT_LIT ("Can't load empty freq file!");
        return TCL_ERROR;
    }

    int mrk;
    for (loc = 0; loc < _nloci; loc++) {
        if (marker && (mrk = marker->get_marker(_locus[loc]->name)) >= 0) {
            strcpy(_locus[loc]->name, marker->mrkname(mrk));

            sprintf(cmd, "cd d_%s", _locus[loc]->name);
            if (Solar_Eval(interp, cmd) == TCL_ERROR) {
                sprintf(errmsg, "\nCannot change to directory d_%s.",
                        _locus[loc]->name);
                RESULT_BUF (errmsg);
                return TCL_ERROR;
            }

            write_locfiles(_locus[loc], true);

            sprintf(cmd, "cd ..");
            if (Solar_Eval(interp, cmd) == TCL_ERROR) {
                RESULT_LIT ("\nCannot return to current working directory.");
                return TCL_ERROR;
            }
        }
    }

    if (write_info(interp) == TCL_ERROR)
        return TCL_ERROR;

    return TCL_OK;
}

int Freq::write_info (Tcl_Interp *interp)
{
    FILE *fp = fopen("freq.info", "w");
    if (!fp) {
        RESULT_LIT ("Cannot open freq.info");
        return TCL_ERROR;
    }

    fprintf(fp, "%s\n", _filename);

    int i, j;
    for (i = 0; i < _nloci; i++) {
        fprintf(fp, "%s %c %c %c", _locus[i]->name, _locus[i]->xlinked,
                _locus[i]->mle_status, _locus[i]->whence);
        for (j = 0; j < _locus[i]->nall; j++) {
            fprintf(fp, " %s %f", _locus[i]->all[j]->name,
                    _locus[i]->all[j]->freq);
            if (_locus[i]->all[j]->mle_se >= 0)
                fprintf(fp, " se= %f", _locus[i]->all[j]->mle_se);
        }
        if (_locus[i]->chi2_hwe >= 0)
            fprintf(fp, " hwe= %f", _locus[i]->chi2_hwe);
        fprintf(fp, "\n");
    }

    fclose(fp);
    return TCL_OK;
}

int Freq::get_marker (const char *name)
{
    int i;

    for (i = 0; i < _nloci; i++)
        if (!StringCmp(_locus[i]->name, name, case_ins))
            return i;

    return -1;
}

bool Freq::unload (Tcl_Interp *interp)
{
    if (currentPed) {
        Marker *marker;
        try {
            marker = currentPed->marker();
        }
        catch (Safe_Error_Return& ser) {
            printf("%s\nReloading the marker data is advised.\n", ser.message());
            return false;
        }

        if (marker) {
            int i, j, n = 0;

            for (i = 0; i < _nloci; i++) {

//  if no genotype data currently loaded, remove freq data
                if (marker->get_marker(_locus[i]->name) < 0) {
                    for (j = 0; j < _locus[i]->nall; j++)
                        free(_locus[i]->all[j]);
                    free(_locus[i]);
                }

//  else, retain freq data but set whence to show freq file unloaded
                else {
                    _locus[n] = _locus[i];
                    if (_locus[n]->whence == 'f')
                        _locus[n]->whence = 'u';
                    n++;
                }
            }

            _nloci = n;
            if (_nloci) {

//  freq data for loaded markers has been retained, but no freq file now
                change_filename("[marker-data-file]");
                if (write_info(interp) == TCL_ERROR) {
                    printf("%s\n", Tcl_GetStringResult (interp));
                    return true;
                }
                return false;
            }
        }
    }

    unlink("freq.info");

    return true;
}

void Freq::remove_marker (int loc)
{
    int i;

    for (i = 0; i < _locus[loc]->nall; i++)
        free(_locus[loc]->all[i]);
    free(_locus[loc]);

    _nloci--;
    for (i = loc; i < _nloci; i++)
        _locus[i] = _locus[i+1];
}

void Freq::remove_marker (const char *name)
{
    int i, j, n = 0;

    for (i = 0; i < _nloci; i++) {
        if (!strcmp(_locus[i]->name, name)) {
            for (j = 0; j < _locus[i]->nall; j++)
                free(_locus[i]->all[j]);
            free(_locus[i]);
        }
        else {
            _locus[n] = _locus[i];
            n++;
        }
    }

    _nloci = n;
}

int Freq::mlefreq (const char *mrkname, bool get_stderr, bool test_hwe,
                   Tcl_Interp *interp)
{
    char errmsg[1024];
    int loc = get_marker(mrkname);
    if (loc < 0) {
        sprintf(errmsg, "%s: No such marker", mrkname);
        RESULT_BUF (errmsg);
        return TCL_ERROR;
    }

    int mrk = currentPed->marker()->get_marker(mrkname);
    if (mrk < 0) {
        sprintf(errmsg,
                "Genotype data have not been loaded for marker %s.", mrkname);
        RESULT_BUF (errmsg);
        return TCL_ERROR;
    }

    struct LocusFreq *locus = _locus[loc];
    if (locus->mle_status == 'c' || locus->mle_status == 's') {
        if (!test_hwe) {
            printf(
        "MLE allele frequencies have already been computed for marker %s.\n",
                   locus->name);
            fflush(stdout);
            return TCL_OK;
        }
        else if (currentFreq->chi2_hwe(currentFreq->get_marker(mrkname)) >= 0) {
            printf(
        "The HWE test statistic has already been computed for marker %s.\n",
                   locus->name);
            fflush(stdout);
            return TCL_OK;
        }
    }

    char mle_cmd[1024];
    FILE *infp, *outfp;

    sprintf(mle_cmd, "cd d_%s", locus->name);
    if (Solar_Eval(interp, mle_cmd) == TCL_ERROR) {
        sprintf(errmsg, "\nCannot change to directory d_%s.", locus->name);
        RESULT_BUF (errmsg);
        return TCL_ERROR;
    }

    if (test_hwe) {
        unlink("genfreq.out");
        unlink("genfreq.frq");
    }

    if (locus->mle_status != 'c' && locus->mle_status != 's') {
        printf("Running allfreq for marker %s ... ", locus->name);
        fflush(stdout);

        unlink("allfreq.out");
        unlink("allfreq.frq");
        if (get_stderr)
            unlink("allfreq.se");

        char show_status = 'n';
        FILE *fp = fopen("/dev/tty", "w");
        if (fp) {
            show_status = 'y';
            fclose(fp);
        }

        sprintf(mle_cmd, "exec allfreq %c %c", show_status, get_stderr?'y':'n');
        if (Solar_Eval(interp, mle_cmd) == TCL_ERROR) {
            sprintf(mle_cmd, "cd ..");
            if (Solar_Eval(interp, mle_cmd) == TCL_ERROR) {
                RESULT_LIT ("\nCannot return to current working directory.");
                return TCL_ERROR;
            }

            FILE *fp;
            char file[1024];
            sprintf(file, "d_%s/allfreq.out", locus->name);

            fp = fopen(file, "r");
            if (!fp)
                RESULT_LIT ("Program allfreq failed. Not enough memory?");
            else {
                sprintf(errmsg,
                    "\nProgram allfreq failed. Check the file %s for errors.",
                        file);

                char rec[1024], cnum[10];
                int ped, ind;
                while (fgets(rec, sizeof(rec), fp)) {
                    if (strncmp(rec, " *** ERROR ***", 14))
                        continue;
                    strncpy(cnum, &rec[46], 5);
                    cnum[5] = 0;
                    ped = atoi(cnum);
                    if (!fgets(rec, sizeof(rec), fp))
                        continue;
                    if (!strncmp(rec, " LOCUS ", 7)) {
                        strncpy(cnum, &rec[34], 4);
                        cnum[4] = 0;
                        ind = atoi(cnum);
                    }
                    else if (!strncmp(rec, " NEAR PERSON ", 13)) {
                        strncpy(cnum, &rec[19], 4);
                        if (sscanf(cnum, "%d", &ind) != 1)
                            continue;
                    }
                    else
                        continue;

                    for (int i = 0; i < ped - 1; i++)
                        ind += currentPed->nind(i);
                    sprintf(errmsg, "%d", ind);
                    break;
                }

                RESULT_BUF (errmsg);
                fclose(fp);
            }

            return TCL_ERROR;

        }
        printf("\n");
        fflush(stdout);
    }

    infp = fopen("allfreq.frq", "r");
    if (!infp) {
        sprintf(mle_cmd, "cd ..");
        if (Solar_Eval(interp, mle_cmd) == TCL_ERROR) {
            RESULT_LIT ("Cannot return to current working directory.");
            return TCL_ERROR;
        }
        sprintf(errmsg, "Cannot open d_%s/allfreq.frq.", locus->name);
        RESULT_BUF (errmsg);
        return TCL_ERROR;
    }

    char rec[1024];
    char aname[MALLNM+1], cfreq[10];
    double afreq, sumfreq = 0;
    int maxfreq = 0;

    for (int i = 0; i < locus->nall; i++) {
        if (!fgets(rec, sizeof(rec), infp) ||
                sscanf(rec, "%s %lf", aname, &afreq) != 2 ||
                strcmp(aname, locus->all[i]->name))
        {
            sprintf(mle_cmd, "cd ..");
            if (Solar_Eval(interp, mle_cmd) == TCL_ERROR) {
                RESULT_LIT ("Cannot return to current working directory.");
                return TCL_ERROR;
            }
            fclose(infp);
            sprintf(errmsg, "Error reading d_%s/allfreq.frq.", locus->name);
            RESULT_BUF (errmsg);
            return TCL_ERROR;
        }
        sprintf(cfreq, "%8.6f", afreq);
        sscanf(cfreq, "%lf", &locus->all[i]->freq);
        if (locus->all[i]->freq < 0)
            locus->all[i]->freq = 0;
        sumfreq += locus->all[i]->freq;
        if (locus->all[i]->freq > locus->all[maxfreq]->freq)
            maxfreq = i;
    }

    double loglike0 = 0;
    if (fgets(rec, sizeof(rec), infp) && sscanf(rec, "%lf", &afreq) == 1)
        loglike0 = afreq;

    fclose(infp);

    if (sumfreq > 1) {
        sprintf(cfreq, "%8.6f", locus->all[maxfreq]->freq - sumfreq + 1);
        sscanf(cfreq, "%lf", &locus->all[maxfreq]->freq);
    }

    if (get_stderr) {
        infp = fopen("allfreq.se", "r");
        if (!infp) {
            sprintf(mle_cmd, "cd ..");
            if (Solar_Eval(interp, mle_cmd) == TCL_ERROR) {
                RESULT_LIT ("Cannot return to current working directory.");
                return TCL_ERROR;
            }
            sprintf(errmsg, "Cannot open d_%s/allfreq.se.", locus->name);
            RESULT_BUF (errmsg);
            return TCL_ERROR;
        }

        for (int i = 0; i < locus->nall; i++) {
            if (!fgets(rec, sizeof(rec), infp) ||
                    sscanf(rec, "%lf", &afreq) != 1)
            {
                sprintf(mle_cmd, "cd ..");
                if (Solar_Eval(interp, mle_cmd) == TCL_ERROR) {
                    RESULT_LIT ("Cannot return to current working directory.");
                    return TCL_ERROR;
                }
                fclose(infp);
                sprintf(errmsg, "Error reading d_%s/allfreq.se.", locus->name);
                RESULT_BUF (errmsg);
                return TCL_ERROR;
            }
            locus->all[i]->mle_se = afreq;
        }

        fclose(infp);
    }

    if (test_hwe) {
        printf("Running genfreq for marker %s ... ", locus->name);
        fflush(stdout);

        char show_status = 'n';
        FILE *fp = fopen("/dev/tty", "w");
        if (fp) {
            show_status = 'y';
            fclose(fp);
        }

        sprintf(mle_cmd, "exec genfreq %c n", show_status);
        if (Solar_Eval(interp, mle_cmd) == TCL_ERROR) {
            sprintf(mle_cmd, "cd ..");
            if (Solar_Eval(interp, mle_cmd) == TCL_ERROR) {
                RESULT_LIT ("\nCannot return to current working directory.");
                return TCL_ERROR;
            }

            FILE *fp;
            char file[1024];
            sprintf(file, "d_%s/genfreq.out", locus->name);

            fp = fopen(file, "r");
            if (!fp)
                RESULT_LIT ("Program genfreq failed. Not enough memory?");
            else {
                sprintf(errmsg,
                        "\nProgram genfreq failed. Check the file %s for errors.",
                        file);
                RESULT_BUF (errmsg);
                fclose(fp);
            }
            return TCL_ERROR;
        }

        infp = fopen("genfreq.frq", "r");
        if (!infp) {
            sprintf(mle_cmd, "cd ..");
            if (Solar_Eval(interp, mle_cmd) == TCL_ERROR) {
                RESULT_LIT ("\nCannot return to current working directory.");
                return TCL_ERROR;
            }
            sprintf(errmsg, "\nCannot open d_%s/genfreq.frq.", locus->name);
            RESULT_BUF (errmsg);
            return TCL_ERROR;
        }

        for (int i = 0; i < locus->nall*(locus->nall+1)/2; i++) {
            if (!fgets(rec, sizeof(rec), infp)) {
                sprintf(mle_cmd, "cd ..");
                if (Solar_Eval(interp, mle_cmd) == TCL_ERROR) {
                    RESULT_LIT ("\nCannot return to current working directory.");
                    return TCL_ERROR;
                }
                fclose(infp);
                sprintf(errmsg, "\nError reading d_%s/genfreq.frq.", locus->name);
                RESULT_BUF (errmsg);
                return TCL_ERROR;
            }
        }

        if (!fgets(rec, sizeof(rec), infp) || sscanf(rec, "%lf", &afreq) != 1)
        {
            sprintf(mle_cmd, "cd ..");
            if (Solar_Eval(interp, mle_cmd) == TCL_ERROR) {
                RESULT_LIT ("\nCannot return to current working directory.");
                return TCL_ERROR;
            }
            fclose(infp);
            sprintf(errmsg, "\nError reading d_%s/genfreq.frq.", locus->name);
            RESULT_BUF (errmsg);
            return TCL_ERROR;
        }

        locus->chi2_hwe = 2*(afreq - loglike0);
        fclose(infp);

        printf("\n");
        fflush(stdout);
    }

    write_locfiles(locus, false);

    sprintf(mle_cmd, "cd ..");
    if (Solar_Eval(interp, mle_cmd) == TCL_ERROR) {
        RESULT_LIT ("Cannot return to current working directory.");
        return TCL_ERROR;
    }

    locus->mle_status = 'c';

    return TCL_OK;
}

void Freq::write_locfiles (struct LocusFreq *locus, bool allfreq)
{
    int i, j;
    FILE *outfp;

    outfp = fopen("datafile.dat", "r");
    if (!outfp)
        outfp = fopen("datain.dat", "r");

    if (outfp) {
        fclose(outfp);
        outfp = fopen("datafile.dat", "w");

        if (locus->xlinked == 'y')
            fprintf(outfp, "2 1 1 5\n");
        else
            fprintf(outfp, "2 1 0 5\n");

        fprintf(outfp, "0 0.00000000 0.00000000 0\n");
        fprintf(outfp, " 1 2\n");
        fprintf(outfp, "1 2\n");
        fprintf(outfp, " 0.99999999 0.00000001\n");
        fprintf(outfp, " 1\n");
        fprintf(outfp, " 0.00000000 0.00000000 1.00000000\n");
        if (locus->xlinked == 'y')
            fprintf(outfp, " 0.00000000 0.50000000\n");
        fprintf(outfp, "2\n\n");

        if (locus->nall == 1) {
            fprintf(outfp, "3 2\n");
            fprintf(outfp, " 0.90000000 0.10000000\n");
        }
        else {
            fprintf(outfp, "3 %d\n", locus->nall);
            if (locus->nall)
                fprintf(outfp, "%11.8f", locus->all[0]->freq);
            for (i = 1; i < locus->nall; i++)
                fprintf(outfp, "%11.8f", locus->all[i]->freq);
            fprintf(outfp, "\n");
        }

        fprintf(outfp, "\n0 0\n");
        fprintf(outfp, " 0.00000000\n");
        fprintf(outfp, "1 0.10000000 0.09000000\n");

        fclose(outfp);
    }

    outfp = fopen("ibd.loc", "r");
    if (outfp) {
        fclose(outfp);
        outfp = fopen("ibd.loc", "w");

        if (locus->xlinked == 'y')
            fprintf(outfp,"%-8.8sX-LINKED%2d%3d\n", locus->name, locus->nall,
                    locus->nall*(locus->nall+1)/2);
        else
            fprintf(outfp,"%-8.8sAUTOSOME%2d%3d\n", locus->name, locus->nall,
                    locus->nall*(locus->nall+1)/2);
 
        for (i = 0; i < locus->nall; i++)
            fprintf(outfp, "%2d      %8.7f\n", i + 1, locus->all[i]->freq);
 
        for (i = 0; i < locus->nall; i++) {
            for (j = i; j < locus->nall; j++) {
                fprintf(outfp, " %2d %2d   1\n", i + 1, j + 1);
                fprintf(outfp, "%2d/%2d\n", i + 1, j + 1);
            }
        }

        fclose(outfp);
    }

    if (!allfreq)
        return;

    outfp = fopen("allfreq.loc", "w");
    if (locus->xlinked == 'y')
        fprintf(outfp,"%-8.8sX-LINKED%2d\n", locus->name, locus->nall);
    else
        fprintf(outfp,"%-8.8sAUTOSOME%2d\n", locus->name, locus->nall);
    for (i = 0; i < locus->nall; i++)
        fprintf(outfp, "%5d   %8.7f\n", i + 1, locus->all[i]->freq);
    fclose(outfp);

    outfp = fopen("allfreq.bat", "w");
    fprintf(outfp,"9\n");
    fprintf(outfp,"%-8.8s\n", locus->name);
    fprintf(outfp,"17\n");
    fprintf(outfp,"%2d\n", locus->nall);
    fprintf(outfp,"21\n");
    fprintf(outfp,"n\n");
    fclose(outfp);

    outfp = fopen("allfreq.mod", "w");
    for (i = 0; i < locus->nall; i++)
        fprintf(outfp,
                "%2d%-5s      %8.6fD+00   0.100000D-05   0.100000D+01\n",
                i + 1, locus->all[i]->name, locus->all[i]->freq);
    fprintf(outfp, "CNS LINES=%2d\n", locus->nall);
    for (i = 0; i < locus->nall; i++)
        fprintf(outfp, "  1 %2d 0.1D+01\n", i + 1);
    fprintf(outfp, "CVALUES  = 1\n");
    fprintf(outfp, "     1 0.1D+01\n");
    fclose(outfp);
}

int Freq::save (const char *fname, Tcl_Interp *interp)
{
    FILE *locfp = fopen(fname, "w");
    if (!locfp) {
        Solar_AppendResult2(interp, "Unable to write freq file:",
                         fname, NULL);
        return TCL_ERROR;
    }

    int i, j;
    for (i = 0; i < _nloci; i++) {
        fprintf(locfp, "%s", _locus[i]->name);
        for (j = 0; j < _locus[i]->nall; j++) {
            fprintf(locfp, " %s %8.6f", _locus[i]->all[j]->name,
                    _locus[i]->all[j]->freq);
            if (_locus[i]->all[j]->mle_se >= 0)
                fprintf(locfp, " se= %8.6f", _locus[i]->all[j]->mle_se);
        }
        if (_locus[i]->chi2_hwe >= 0)
            fprintf(locfp, " hwe= %8.6f", _locus[i]->chi2_hwe);
        fprintf(locfp, "\n");
        if (_locus[i]->mle_status == 'c')
            _locus[i]->mle_status = 's';
    }

    fclose(locfp);
    return TCL_OK;
}

const char *Freq::show (char *buf)
{
    int i, j;
    char buf1[1024], status[1024];

    if (!strcmp(buf, "")) {
        sprintf(buf, "\nfrequency data file: %s\n", _filename);
        for (i = 0; i < _nloci; i++) {
            switch (_locus[i]->xlinked) {
            case 'n': status[0] = '\0';
                      break;
            case 'y': strcpy(status, ", X-linked");
                      break;
            default:  strcpy(status, ", **illegal X-linked code**");
            }

            switch (_locus[i]->mle_status) {
            case 'n': break;
            case 'c': strcat(status, ", MLEs computed");
                      break;
            case 's': strcat(status, ", MLEs saved");
                      break;
            case 'o': strcat(status, ", stale MLEs");
                      break;
            default:  strcat(status, ", **illegal MLE status**");
            }

            sprintf(buf1, "\n%s: %d alleles%s\n", _locus[i]->name,
                    _locus[i]->nall, status);
            strcat(buf, buf1);

            j = 0;
            while (j < _locus[i]->nall) {
                sprintf(buf1, "\t%s  %8.6f", _locus[i]->all[j]->name,
                        _locus[i]->all[j]->freq);
                strcat(buf, buf1);
                if (!(++j%4)) strcat(buf, "\n");
            }
            if (j%4) strcat(buf, "\n");
        }
    }

    else {
        if ((i = get_marker(buf)) < 0)
            strcat(buf, ": No such marker.");
        else {
            switch (_locus[i]->xlinked) {
            case 'n': status[0] = '\0';
                      break;
            case 'y': strcpy(status, ", X-linked");
                      break;
            default:  strcpy(status, ", **illegal X-linked code**");
            }

            switch (_locus[i]->mle_status) {
            case 'n': break;
            case 'c': strcat(status, ", MLEs computed");
                      break;
            case 's': strcat(status, ", MLEs saved");
                      break;
            case 'o': strcat(status, ", stale MLEs");
                      break;
            default:  strcat(status, ", **illegal MLE status**");
            }

            sprintf(buf, "%s: %d alleles%s\n", _locus[i]->name,
                    _locus[i]->nall, status);

            j = 0;
            while (j < _locus[i]->nall) {
                sprintf(buf1, "\t%s  %8.6f", _locus[i]->all[j]->name,
                        _locus[i]->all[j]->freq);
                strcat(buf, buf1);
                if (!(++j%4)) strcat(buf, "\n");
            }
            if (j%4) strcat(buf, "\n");
        }
    }

    return buf;
}
