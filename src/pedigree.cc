/*
 * pedigree.cc implements the pedigree command and class
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

Pedigree *currentPed = 0;

void delete_ped_state (void);

extern "C" int PedigreeCmd (ClientData clientData, Tcl_Interp *interp,
		  int argc, char *argv[])
{
    if (argc == 2 && !StringCmp ("help", argv[1], case_ins)) {
	return Solar_Eval (interp, "help pedigree");
    }

    else if (argc >= 2 && !StringCmp ("load", argv[1], case_ins)) {
        if (argc == 2) {
            RESULT_LIT ("Usage: pedigree load <filename>");
            return TCL_ERROR;
        }

        bool all_founders = false;
        if (argc == 4 && !StringCmp ("all_founders", argv[3], case_ins))
            all_founders = true;

    // load a new pedigree

	// Notify other classes
	Pedigree::Changing_Pedigree();

        // Flush all working EVD storage
	EVD::Flush();

        bool ped_loaded;
        try {
            ped_loaded = loadedPed();
        }
        catch (Safe_Error_Return) {
            ped_loaded = true;
        }

        if (ped_loaded) {
            printf("Unloading current pedigree data ...\n");
            fflush(stdout);
            if (currentPed) {
                try {
                    if (currentPed->marker())
                        currentPed->marker()->unload(true, interp);
                }
                catch (Safe_Error_Return) {
                }

                delete currentPed;
            }
            delete_ped_state();
            Phenotypes::reset();
        }

        currentPed = new Pedigree (argv[2]);
        printf("Loading pedigree data from the file %s ...\n", argv[2]);
        fflush(stdout);

        if (currentPed->load(all_founders, interp) == TCL_ERROR) {
            delete currentPed;
            delete_ped_state();
            return TCL_ERROR;
        }

        if (MMSibs) {
            int i;
            for (i = 0; i < currentPed->nped(); i++) {
                if (currentPed->nfam(i) > 1) {
                    Solar_AppendResult2(interp,
                        "MAPMAKER/SIBS can only be used for sibships! ",
                        "Set MMSibs option to no.", NULL);
                    return TCL_ERROR;
                }
            }
        }

        return TCL_OK;
    }

    else if (argc >= 2 && !StringCmp ("show", argv[1], case_ins)) {
    // display pedigree info
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

        char *buf = (char *) malloc(1200+40*(currentPed->nped()+1));
        if (argc == 2)
            strcpy(buf, "totals");
        else
            strcpy(buf, argv[2]);

        char tmpfile[10];
        strcpy(tmpfile, "tmpXXXXXX");
        int tmpfd = mkstemp(tmpfile);
        if (tmpfd == -1)
            return TCL_ERROR;

        FILE *tmpfp = fdopen(tmpfd, "w");
        if (tmpfp) {
            fprintf(tmpfp, "%s\n", currentPed->show(buf));
            fclose(tmpfp);
            char more_cmd[1024];
            sprintf(more_cmd, "exec >&@stdout <@stdin more %s", tmpfile);
            if (Solar_Eval(interp, more_cmd) == TCL_ERROR)
                return TCL_ERROR;
            remove(tmpfile);
        }
        else
            printf("%s\n", currentPed->show(buf));

        free(buf);
        return TCL_OK;
    }

    else if (argc == 2 && !StringCmp ("fname", argv[1], case_ins)) {
    // return pedigree file name
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

        char buf[1024];
        sprintf(buf, "%s", currentPed->filename());
        RESULT_BUF (buf);
        return TCL_OK;
    }

    RESULT_LIT ("Invalid pedigree command");
    return TCL_ERROR;
}

Pedigree::Pedigree (const char *fname)
{
    strcpy(_filename, fname);
    Tfile = 0;
    _widths = 0;
    _count = 0;
    _id_len = 0;
    _sex_len = 0;
    _mztwin_len = 0;
    _hhid_len = 0;
    _famid_len = 0;
    _nped = 0;
    _marker = 0;
}

Pedigree::~Pedigree ()
{
    int i;
    for (i = 0; i < _nped; i++)
        free(_ped[i]);

    delete_Tfile();
    delete_marker();

    currentPed = 0;
}

void delete_ped_state ()
{
    unlink("mibdrel.in");
    unlink("mibdrel.ped");
    unlink("mibdrel.cde");
    unlink("mibdrel.tab");
    unlink("pedindex.out");
    unlink("pedindex.cde");
    unlink("pedigree.info");
    unlink("phi2.gz");
    unlink("house.gz");
}

bool loadedPed ()
{
    if (currentPed)
        return true;

    FILE *fp = fopen("pedindex.out", "r");
    if (!fp)
        return false;
    fclose(fp);

    char fname[1024];
    fp = fopen("pedigree.info", "r");
    if (fp) {
        fscanf(fp, "%s", fname);
        fclose(fp);
    }
    else
        strcpy(fname, "unknown");

    currentPed = new Pedigree (fname);

    const char *errmsg = 0;
    if (!currentPed->get_stats(&errmsg)) {
        delete currentPed;
        throw Safe_Error_Return(errmsg);
    }

    return true;
}

int Pedigree::load (bool all_founders, Tcl_Interp *interp)
{
    char load_cmd[1024], buferr[1024];
    const char *errmsg = 0;
    FILE *fp;

    Tfile = SolarFile::open ("pedigree data", _filename, &errmsg);
    if (errmsg) {
        Solar_AppendResult2(interp, "Error opening pedigree data file:\n",
                            errmsg, NULL);
        return TCL_ERROR;
    }

    _widths = Tfile->widths (&_count, &errmsg);
    if (errmsg) {
        Solar_AppendResult2(interp, "Error loading pedigree data file:\n",
                            errmsg, NULL);
        return TCL_ERROR;
    }

// Set up fields and get their widths

    int i;
    int fidlen = 0;  // These aren't saved in object
    int midlen = 0;
    Tfile->start_setup(&errmsg);

    if (Tfile->test_name ("famid", &errmsg))
    {
	i = Tfile->setup ("famid", &errmsg);
	_famid_len = _widths[i];
    }

    i = Tfile->setup ("id", &errmsg);
    _id_len = _widths[i];

    if (!all_founders) {
        i = Tfile->setup ("fa", &errmsg);
        fidlen = _widths[i];

        i = Tfile->setup ("mo", &errmsg);
        midlen = _widths[i];
    }

    if (!all_founders) {
        i = Tfile->setup ("sex", &errmsg);
        _sex_len = _widths[i];
    } else if (Tfile->test_name ("sex", &errmsg)) {
        i = Tfile->setup ("sex", &errmsg);
        _sex_len = _widths[i];
    } else {
        _sex_len = 0;
    }

    if (Tfile->test_name ("mztwin", &errmsg))
    {
	i = Tfile->setup ("mztwin", &errmsg);
	_mztwin_len = _widths[i];
        if (!_mztwin_len) _mztwin_len = 1;
    }

    if (Tfile->test_name ("hhid", &errmsg))
    {
	i = Tfile->setup ("hhid", &errmsg);
	_hhid_len = _widths[i];
        if (!_hhid_len) _hhid_len = 1;
    }

// Check for error(s)  (missing fields, etc.)

    if (errmsg) {
	char mbuf[1024];
	sprintf (mbuf, "Pedigree data file error: %s", errmsg);
	RESULT_BUF (mbuf);
        return TCL_ERROR;
    }

    if ((!all_founders && _sex_len != 1) ||
        (all_founders && _sex_len > 1))
    {
        char mbuf[1024];
        sprintf (mbuf,
        "The sex field has a max width of %d, but should be a single character.",
                 _sex_len);
	RESULT_BUF (mbuf);
        return TCL_ERROR;
    }

// ID fields "should" all be same width, so make it so

    if (fidlen > _id_len) _id_len = fidlen;
    if (midlen > _id_len) _id_len = midlen;

    fp = fopen("ibdprep.ped", "w");
    if (!fp) {
        RESULT_LIT ("Cannot open ibdprep.ped");
        return TCL_ERROR;
    }

    Tfile->rewind(&errmsg);
    char **record;
    while (1) {
        record = Tfile->get(&errmsg);
        if (errmsg && !strcmp("EOF", errmsg)) break;
        if (errmsg) {
	    char mbuf[1024];
	    sprintf (mbuf, "Pedigree data file error: %s", errmsg);
	    RESULT_BUF (mbuf);
            return TCL_ERROR;
        }
        i = 0;
        if (_famid_len)
            fprintf(fp, "%*s", _famid_len, record[i++]);
        fprintf(fp, "%*s", _id_len, record[i++]);
        if (!all_founders) {
            fprintf(fp, "%*s", _id_len, record[i++]);
            fprintf(fp, "%*s", _id_len, record[i++]);
        } else {
            fprintf(fp, "%*s", _id_len, " ");
            fprintf(fp, "%*s", _id_len, " ");
        }
        if (_sex_len)
            fprintf(fp, "%*s", _sex_len, record[i++]);
        else
            fprintf(fp, "U");
        if (_mztwin_len)
            fprintf(fp, "%*s", _mztwin_len, record[i++]);
        if (_hhid_len)
            fprintf(fp, "%*s", _hhid_len, record[i++]);
        fprintf(fp, "\n");
    }
    fclose(fp);
    delete_Tfile();

    if (all_founders && _sex_len == 0)
        _sex_len = 1;

// Create initial pedigree.info, contains name of pedigree data file only

    fp = fopen("pedigree.info", "w");
    if (!fp) {
        RESULT_LIT ("Cannot open pedigree.info");
        return TCL_ERROR;
    }
    fprintf(fp, "%s\n", _filename);
    fclose(fp);

// Run program ibdprep with appropriate command line arguments

    if (_famid_len) {
        sprintf(load_cmd,
                "exec ibdprep y ibdprep.ped %d %d %d %d %d >& ibdprep.out",
                _id_len, _sex_len, _mztwin_len, _hhid_len, _famid_len);
    }
    else {
        sprintf(load_cmd,
                "exec ibdprep y ibdprep.ped %d %d %d %d >& ibdprep.out",
                _id_len, _sex_len, _mztwin_len, _hhid_len);
    }

// Check for error output from program ibdprep

    if (Solar_Eval(interp, load_cmd) == TCL_ERROR) {
        if (Solar_Eval(interp, "exec cat ibdprep.out") == TCL_ERROR) {
            RESULT_LIT ("Cannot cat ibdprep.out");
            return TCL_ERROR;
        }
        if (!strlen(Tcl_GetStringResult (interp)))
            RESULT_LIT ("Program ibdprep failed to run");
        return TCL_ERROR;
    } else {
        Solar_Eval(interp, "exec cat ibdprep.out");
    }

// Read pedigree stats from pedigree.info

    if (!get_stats(&errmsg)) {
        RESULT_LIT (errmsg);
        return TCL_ERROR;
    }

    unlink("ibdprep.ped");
    unlink("ibdprep.out");
    return TCL_OK;
}

bool Pedigree::get_stats (const char **errmsg)
{
    FILE *fp = fopen("pedigree.info", "r");
    if (!fp) {
        *errmsg = "Cannot open pedigree.info";
        return false;
    }

    int i;
    char rec[1024], fname[1024];

    if (!fgets(rec, sizeof(rec), fp) ||
            sscanf(rec, "%s", fname) != 1) {
        fclose(fp);
        *errmsg = "Read error on pedigree.info, line 1";
        return false;
    }

    if (!fgets(rec, sizeof(rec), fp) ||
            sscanf(rec, "%d %d %d %d %d", &_id_len, &_sex_len, &_mztwin_len,
                   &_hhid_len, &_famid_len) != 5) {
        fclose(fp);
        *errmsg = "Read error on pedigree.info, line 2";
        return false;
    }

    if (!fgets(rec, sizeof(rec), fp) ||
            sscanf(rec, "%d %d %d %d", &_nped, &_nfam, &_nind,
                   &_nfou) != 4) {
        fclose(fp);
        _nped = 0;     // keeps Pedigree destructor from trying to free
                       // unallocated _ped[] structs
        *errmsg = "Read error on pedigree.info, line 3";
        return false;
    }

    _ped = (struct Ped **) malloc(_nped * sizeof(struct Ped *));
    for (i = 0; i < _nped; i++)
        _ped[i] = 0;

    for (i = 0; i < _nped; i++) {
        _ped[i] = (struct Ped *) malloc(sizeof(struct Ped));
        if (!_ped[i]) { 
            fclose(fp);
            *errmsg = "Read error on pedigree.info: not enough memory";
            return false; 
        }

        if (!fgets(rec, sizeof(rec), fp) ||
                sscanf(rec, "%d %d %d %d %c", &_ped[i]->nfam,
                       &_ped[i]->nind, &_ped[i]->nfou, &_ped[i]->nlbrk,
                       &_ped[i]->inbred) != 5) {
            fclose(fp);
            sprintf(error_message,
                    "Read error on pedigree.info, line %d", i + 4);
            *errmsg = error_message;
            return false;
        }

        if (_ped[i]->nlbrk > 1 || _ped[i]->inbred == 'y')
            MCarlo = true;
    }

    fclose(fp);
    *errmsg = 0;
    return true;
}

char *Pedigree::show (char *buf)
{
    int i;
    char clbrk[10], line[1024];

    if (!strcmp(buf, "totals")) {
        sprintf(buf, "\npedigree data file: %s\n\n", _filename);
        sprintf(line, "%5d pedigrees\n", _nped);
        strcat(buf, line);
        sprintf(line, "%5d nuclear families\n", _nfam);
        strcat(buf, line);
        sprintf(line, "%5d individuals\n", _nind);
        strcat(buf, line);
        sprintf(line, "%5d founders\n", _nfou);
        strcat(buf, line);

        int mlbrk = 0;
        bool inbred = false;
        for (i = 0; i < _nped; i++) {
            if (_ped[i]->nlbrk > mlbrk)
                mlbrk = _ped[i]->nlbrk;
            if (_ped[i]->inbred == 'y')
                inbred = true;
        }

        if (mlbrk > 1 && inbred)
            strcat(buf,
        "\nMultiple loops, inbreeding present in 1 or more pedigrees.\n");
        else if (mlbrk > 1)
            strcat(buf,
                   "\nMultiple loops present in 1 or more pedigrees.\n");
        else if (inbred)
            strcat(buf,
                   "\nInbreeding present in 1 or more pedigrees.\n");
    }

    else if (!strcmp(buf, "all")) {
        sprintf(buf, "\npedigree data file: %s\n\n", _filename);
        strcat(buf, "ped#\t#nfam\t #ind\t #fou\t#bits\t#lbrk\tinbred?\n");
        strcat(buf, "----\t-----\t-----\t-----\t-----\t-----\t-------\n");

        for (i = 0; i < _nped; i++) {
            if (_ped[i]->nlbrk)
                sprintf(clbrk, "%5d", _ped[i]->nlbrk);
            else
                strcpy(clbrk, "     ");

            if (_ped[i]->nind > 1) {
                sprintf(line, "%4d\t%5d\t%5d\t%5d\t%5d\t%5s\t   %c\n", i + 1,
                        _ped[i]->nfam, _ped[i]->nind, _ped[i]->nfou,
                        2*_ped[i]->nind - 3*_ped[i]->nfou,
                        clbrk, _ped[i]->inbred == 'y' ? 'y' : ' ');
            }
            else {
                sprintf(line, "%4d\t%5d\t%5d\t%5d\t     \t%5s\t   %c\n", i + 1,
                        _ped[i]->nfam, _ped[i]->nind, _ped[i]->nfou,
                        clbrk, _ped[i]->inbred == 'y' ? 'y' : ' ');
            }
            strcat(buf, line);
        }

        strcat(buf, "\t-----\t-----\t-----\n");
        sprintf(line, "\t%5d\t%5d\t%5d\n", _nfam, _nind, _nfou);
        strcat(buf, line);
    }

    else if (sscanf(buf, "%d", &i) != 1 || i < 1 || i > _nped)
        strcpy(buf, "No such pedigree.");

    else {
        strcat(buf, "ped#\t#nfam\t #ind\t #fou\t#bits\t#lbrk\tinbred?\n");
        strcat(buf, "----\t-----\t-----\t-----\t-----\t-----\t-------\n");

        if (_ped[i-1]->nlbrk)
            sprintf(clbrk, "%5d", _ped[i-1]->nlbrk);
        else
            strcpy(clbrk, "     ");

        if (_ped[i-1]->nind > 1) {
            sprintf(line, "%4d\t%5d\t%5d\t%5d\t%5d\t%5s\t   %c\n", i,
                    _ped[i-1]->nfam, _ped[i-1]->nind, _ped[i-1]->nfou,
                    2*_ped[i-1]->nind - 3*_ped[i-1]->nfou,
                    clbrk, _ped[i-1]->inbred == 'y' ? 'y' : ' ');
        }
        else {
            sprintf(line, "%4d\t%5d\t%5d\t%5d\t     \t%5s\t   %c\n", i,
                    _ped[i-1]->nfam, _ped[i-1]->nind, _ped[i-1]->nfou,
                    clbrk, _ped[i-1]->inbred == 'y' ? 'y' : ' ');
        }
        strcat(buf, line);
    }

    return buf;
}

Marker *Pedigree::marker ()
{
    if (!_marker) {
        FILE *fp = fopen("marker.info", "r");
        if (fp) {
            char fname[1024];
            fscanf(fp, "%s", fname);
            fclose(fp);
            Marker *marker = new Marker (fname);
            add_marker(marker);

            const char *errmsg = 0;
            if (!marker->get_stats(&errmsg)) {
                delete_marker();
                throw Safe_Error_Return(errmsg);
            }

            try {
                if (!loadedFreq()) {
                    delete_marker();
                    throw Safe_Error_Return("The file freq.info is missing.");
                }
            }
            catch (Safe_Error_Return& ser) {
                delete_marker();
                throw Safe_Error_Return(ser.message());
            }
        }
    }

    return _marker;
}

void Pedigree::delete_Tfile ()
{
    if (Tfile) {
        delete Tfile;
        Tfile = 0;
    }
}

void Pedigree::delete_marker ()
{
    if (_marker) {
        delete _marker;
        _marker = 0;
    }
}
