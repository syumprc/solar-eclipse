/*
 * marker.cc mplements the marker command and class
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

extern "C" int MarkerCmd (ClientData clientData, Tcl_Interp *interp,
		  int argc, char *argv[])
{
    if (argc == 2 && !StringCmp ("help", argv[1], case_ins)) {
	return Solar_Eval (interp, "help marker");
    }

    else if (argc >= 2 && !StringCmp ("load", argv[1], case_ins)) {
        if (argc != 3 && argc != 4) {
            RESULT_LIT ("Usage: marker load [-xlinked] <filename>");
            return TCL_ERROR;
        }

        char fname[1024];
        char xlinked = XLinked ? 'y' : 'n';

        if (argc == 4) {
            if (strcmp("-xlinked", argv[2])) {
                RESULT_LIT ("Usage: marker load [-xlinked] <filename>");
                return TCL_ERROR;
            }
            else {
                xlinked = 'y';
                strcpy(fname, argv[3]);
            }
        }
        else {
            if (!strcmp("-xlinked", argv[2])) {
                RESULT_LIT ("Usage: marker load [-xlinked] <filename>");
                return TCL_ERROR;
            }
            else
                strcpy(fname, argv[2]);
        }

    // load new marker data
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

            if (!currentMap) {
                RESULT_LIT ("Map data have not been loaded.");
                return TCL_ERROR;
            }
        }

        try {
            if (currentPed->marker() &&
                    !currentPed->marker()->unload(false, interp)) {
                RESULT_LIT (\
"Use freq save to save MLE allele frequencies, then do marker load again.\n\
Or use marker unload -nosave to unload without saving MLEs, then marker load.");
                return TCL_ERROR;
            }
        }
        catch (Safe_Error_Return) {
        }

        currentPed->delete_marker();
        Marker *marker = new Marker (fname);
        currentPed->add_marker(marker);

        if (marker->load(xlinked, interp) == TCL_ERROR) {
            currentPed->delete_marker();
            return TCL_ERROR;
        }

        return TCL_OK;
    }

    else if (argc >= 2 && !StringCmp ("unload", argv[1], case_ins)) {
        if (argc > 3 || argc == 3 && strcmp("-nosave", argv[2])) {
            RESULT_LIT ("Usage: marker unload [-nosave]");
            return TCL_ERROR;
        }

    // unload marker data
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

        bool nosave = false;
        if (argc == 3 && !strcmp("-nosave", argv[2]))
            nosave = true;

        if (!marker->unload(nosave, interp)) {
	    RESULT_LIT (\
"Use freq save to save MLE allele frequencies, then do marker unload again.\n\
Or use marker unload -nosave to unload markers without saving MLEs.");
            return TCL_ERROR;
        }

        currentPed->delete_marker();

        return TCL_OK;
    }

    else if (argc >= 2 && !StringCmp ("discrep", argv[1], case_ins)) {
        if (argc > 3) {
            RESULT_LIT ("Usage: marker discrep [<marker>]");
            return TCL_ERROR;
        }

    // check for marker discrepancies
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
 
        if (argc == 2) {
            int i, retval = TCL_OK;
            for (i = 0; i < marker->nloci(); i++) {
                if (marker->discrep(marker->mrkname(i), interp)
                        == TCL_ERROR) {
                    printf("%s\n", Tcl_GetStringResult (interp)); 
		    fflush(stdout);
		    Tcl_ResetResult (interp);
                    retval = TCL_ERROR;
                }
            }
            return retval;
        }
        else {
            if (marker->discrep(argv[2], interp) == TCL_ERROR)
                return TCL_ERROR;
        }
 
        return TCL_OK;
    }

    else if (argc >= 2 && !StringCmp ("show", argv[1], case_ins)) {
        if (argc > 3) {
            RESULT_LIT ("Usage: marker show [<marker>]");
            return TCL_ERROR;
        }

    // display marker data
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

        char *buf = (char *) malloc(70*(marker->nloci()+3));
        if (argc == 2)
            strcpy(buf, "");
        else
            strcpy(buf, argv[2]);
 
        marker->show(buf, interp);
        RESULT_BUF (buf);
        free(buf);
        return TCL_OK;
    }

    else if (argc == 2 && !StringCmp ("names", argv[1], case_ins)) {
    // display marker names
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

        marker->show_names(interp);

        return TCL_OK;
    }

    else if (argc == 2 && !StringCmp ("fname", argv[1], case_ins)) {

    // return marker file name
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

        char buf[1024];
        sprintf(buf, "%s", marker->filename());
        RESULT_BUF (buf);
        return TCL_OK;
    }

    else if (argc == 3 && !StringCmp ("name", argv[1], case_ins)) {

    // return marker name
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

        char buf[1024];
        int mrk = currentFreq->get_marker(argv[2]);
        if (mrk < 0) {
            sprintf(buf, "%s: No such marker.", argv[2]);
            RESULT_BUF (buf);
            return TCL_ERROR;
        }

        sprintf(buf, "%s", marker->mrkname(mrk));
        RESULT_BUF (buf);
        return TCL_OK;
    }

    RESULT_LIT ("Invalid marker command");
    return TCL_ERROR;
}

Marker::Marker (const char *fname)
{
    strcpy(_filename, fname);
    Tfile = 0;
    _names = 0;
    _widths = 0;
    _count = 0;
    _id_len = 0;
    _famid_len = 0;
    _gtype_len = 0;
    _nloci = 0;
}

Marker::~Marker ()
{
    int i;
    for (i = 0; i < _nloci; i++)
        free(_mrk[i]);

    delete_Tfile();
}
 
int Marker::load (const char xlinked, Tcl_Interp *interp)
{
    char load_cmd[1024];
    const char *errmsg = 0;
    int i, j;
    int fidpos = -1;
    int idpos = -1;
    FILE *fp;

// Open marker data file and get _names and _widths

    Tfile = SolarFile::open ("marker data", _filename, &errmsg);
    if (errmsg) {
        Solar_AppendResult2(interp, "Error opening marker data file:\n", 
			 errmsg, NULL);
        return TCL_ERROR;
    }
    _names = Tfile->names (&_count, &errmsg);
    _widths = Tfile->widths (&_count, &errmsg);
    if (errmsg) {
        Solar_AppendResult2(interp, "Error loading marker data file:\n", 
			 errmsg, NULL);
        return TCL_ERROR;
    }

// Setup famid (if present here and in pedigree file) and id
// get _famid_len and _id_len from marker file

    int loc0 = 0;
    Tfile->start_setup(&errmsg);
    if (currentPed->famid_len() && Tfile->test_name ("famid", &errmsg)) {
	fidpos = Tfile->setup("famid", &errmsg);
	_famid_len = _widths[fidpos];
        loc0++;
    }
    idpos = Tfile->setup("id", &errmsg);
    _id_len = _widths[idpos];
    loc0++;

    int *use_fld = new int[_count];
    for (i = 0; i < _count; i++)
        use_fld[i] = 1;
    use_fld[idpos] = 0;
    if (fidpos != -1) use_fld[fidpos] = 0;
    if (Tfile->test_name("fa", &errmsg)) {
        use_fld[Tfile->setup("fa", &errmsg)] = 0;
        loc0++;
    }
    if (Tfile->test_name("mo", &errmsg)) {
        use_fld[Tfile->setup("mo", &errmsg)] = 0;
        loc0++;
    }
    if (Tfile->test_name("sex", &errmsg)) {
        use_fld[Tfile->setup("sex", &errmsg)] = 0;
        loc0++;
    }
    if (Tfile->test_name("mztwin", &errmsg)) {
        use_fld[Tfile->setup("mztwin", &errmsg)] = 0;
        loc0++;
    }
    if (Tfile->test_name("hhid", &errmsg)) {
        use_fld[Tfile->setup("hhid", &errmsg)] = 0;
        loc0++;
    }
    if (Tfile->test_name("age", &errmsg)) {
        use_fld[Tfile->setup("age", &errmsg)] = 0;
        loc0++;
    }
    if (Tfile->test_name("pedno", &errmsg)) {
        use_fld[Tfile->setup("pedno", &errmsg)] = 0;
        loc0++;
    }
    if (Tfile->test_name("gen", &errmsg)) {
        use_fld[Tfile->setup("gen", &errmsg)] = 0;
        loc0++;
    }
    if (Tfile->test_name("pedsiz", &errmsg)) {
        use_fld[Tfile->setup("pedsiz", &errmsg)] = 0;
        loc0++;
    }
        

// Find maximum genotype length

    for (i = 0; i < _count; i++) {
	if (use_fld[i]) {
	    if (_widths[i] > _gtype_len)
		_gtype_len = _widths[i];
	}
    }

// Setup markers

    _nloci = 0;
    for (i = 0; i < _count; i++) {
	if (use_fld[i]) {
            Tfile->setup(_names[i], &errmsg);
            if (strlen(_names[i]) > MMRKNM) {
		char buf[256];
                sprintf(buf,
            "Marker data file error:\nMarker name longer than %d characters",
                        MMRKNM);
		RESULT_BUF (buf);
                return TCL_ERROR;
            }
            for (j = 0; j < _nloci; j++) {
                if (!strcmp(_names[i], _mrk[j]->name)) {
		    char buf[256];
                    sprintf(buf,
            "Marker data file error:\nMarker %s appears more than once",
                        _names[i]);
		RESULT_BUF (buf);
                return TCL_ERROR;
                }
            }
            if (_nloci == MAXLOC) {
		char buf[256];
                sprintf(buf,
            "Marker data file contains too many loci: MAXLOC = %d", MAXLOC);
		RESULT_BUF (buf);
                return TCL_ERROR;
            }
            _mrk[_nloci] = (struct Mrk *) malloc(sizeof(struct Mrk));
            if (!_mrk[_nloci]) {
                RESULT_LIT ("Not enough memory");
                return TCL_ERROR;
            }
            strcpy(_mrk[_nloci]->name, _names[i]);
            _nloci++;
        }
    }

// Now that we've set everything up, check for errors

    if (errmsg) {
	Solar_AppendResult2 (interp, "Marker data file error:\n",
			  errmsg, NULL);
        return TCL_ERROR;
    }

    if (!_nloci) {
	char buf[256];
        sprintf(buf, "No fields containing marker data are present in %s",
                _filename);
        RESULT_BUF (buf);
        return TCL_ERROR;
    }

// Open ibdprep.mrk

    fp = fopen("ibdprep.mrk", "w");
    if (!fp) {
        RESULT_LIT ("Cannot open ibdprep.mrk");
        return TCL_ERROR;
    }

// Read marker data records and copy to ibdprep.mrk

    Tfile->rewind(&errmsg);
    char **record;
    while (1) {
        record = Tfile->get(&errmsg);
        if (errmsg && !strcmp("EOF", errmsg)) break;
        if (errmsg) {
	    Solar_AppendResult2(interp, "Error reading marker data file:\n",
			     errmsg, NULL);
            fclose(fp);
            return TCL_ERROR;
        }

// Marker file records with IDs/family IDs wider than those in the pedigree
// file can't be for individuals in the pedigree data so don't include them.
// This is a kludge for the sake of ibdprep's fixed-width fields!

        i = 0;
        if (_famid_len) {
            if (strlen(record[i]) > currentPed->famid_len()
                   || strlen(record[i+1]) > currentPed->id_len())
                continue;
            fprintf(fp, "%*s", currentPed->famid_len(), record[i++]);
        }

        if (strlen(record[i]) > currentPed->id_len())
            continue;

        fprintf(fp, "%*s", currentPed->id_len(), record[i++]);

        for (i = loc0 ; i < _count; i++)
            fprintf(fp, "%*s", _gtype_len, record[i]);
        fprintf(fp, "\n");
    }
    fclose(fp);

// Open ibdprep.loc

    fp = fopen("ibdprep.loc", "w");
    if (!fp) {
        RESULT_LIT ("Cannot open ibdprep.loc");
        return TCL_ERROR;
    }

// For each marker
//   if freq info has been loaded, write marker name and allele
//      names/freqs to ibdprep.loc
//   else
//      just write marker name to ibdprep.loc

    if (currentFreq) {
        printf("Using allele freqs loaded from the file %s ...\n",
               currentFreq->filename());
        fflush(stdout);
        bool stale = false;
        for (i = 0; i < _nloci; i++) {
            int j = currentFreq->get_marker(_mrk[i]->name);
            if (j < 0)
                fprintf(fp, "%s\n", _mrk[i]->name);
            else {
                fprintf(fp, "%s", _mrk[i]->name);
                int k;
                for (k = 0; k < currentFreq->nall(j); k++)
                    fprintf(fp, " %s %8.6f", currentFreq->all(j,k)->name,
                            currentFreq->all(j,k)->freq);
                fprintf(fp, "\n");
                if (currentFreq->mle_status(j) == 'o') stale = true;
            }
        }

        if (stale) {
            printf(
"NOTE: Some or all of these allele freqs are stale MLEs, i.e. the MLEs\n\
      were computed with respect to a previously loaded marker file.\n");
            fflush(stdout);
        }
    }
    else {
        printf(
        "Getting initial estimates of allele freqs from the marker data ...\n");
        fflush(stdout);
        for (i = 0; i < _nloci; i++)
            fprintf(fp, "%s\n", _mrk[i]->name);
    }
    fclose(fp);

// Create initial marker.info, contains name of marker data file only
 
    fp = fopen("marker.info", "w");
    if (!fp) {
        RESULT_LIT ("Cannot open marker.info");
        return TCL_ERROR;
    }
    fprintf(fp, "%s\n", _filename);
    fclose(fp);

// Run program ibdprep with appropriate command line arguments

    printf("Loading marker data from the file %s ...\n", _filename);
    fflush(stdout);

    char doMCarlo = MCarlo ? 'y' : 'n';
    char inMrkFile = _famid_len ? 'y' : 'n';

    if (MMSibs) {
        if (currentPed->famid_len()) {
            sprintf(load_cmd,
"exec ibdprep n ibdprep.mrk %d %d %c %d ibdprep.loc %c y %s %d %c >& ibdprep.out",
                    currentPed->id_len(), _gtype_len, xlinked, _nloci, doMCarlo,
                    currentMap->filename(), currentPed->famid_len(), inMrkFile);
        }
        else {
            sprintf(load_cmd,
"exec ibdprep n ibdprep.mrk %d %d %c %d ibdprep.loc %c y %s >& ibdprep.out",
                    currentPed->id_len(), _gtype_len, xlinked, _nloci, doMCarlo,
                    currentMap->filename());
        }
    }
    else {
        if (currentPed->famid_len()) {
            sprintf(load_cmd,
"exec ibdprep n ibdprep.mrk %d %d %c %d ibdprep.loc %c n %d %c >& ibdprep.out",
                    currentPed->id_len(), _gtype_len, xlinked, _nloci, doMCarlo,
                    currentPed->famid_len(), inMrkFile);
        }
        else {
            sprintf(load_cmd,
"exec ibdprep n ibdprep.mrk %d %d %c %d ibdprep.loc %c n >& ibdprep.out",
                    currentPed->id_len(), _gtype_len, xlinked, _nloci, doMCarlo);
        }
    }

// Check for error output from program ibdprep

    if (Solar_Eval(interp, load_cmd) == TCL_ERROR) {
        if (Solar_Eval(interp, "exec cat ibdprep.out") == TCL_ERROR) {
            RESULT_LIT ("Cannot cat ibdprep.out");
            return TCL_ERROR;
        }
        if (!strlen(Tcl_GetStringResult (interp)))
            RESULT_LIT ("Program ibdprep did not run");
        return TCL_ERROR;
    }

// Update frequency data from ibdprep.loc

    if (!currentFreq)
        currentFreq = new Freq ("[marker-data-file]");

    if (currentFreq->get_freqs("ibdprep.loc", xlinked, interp) == TCL_ERROR)
        return TCL_ERROR;

// Read marker stats from marker.info

    if (!get_stats(&errmsg)) {
        RESULT_LIT (errmsg);
        return TCL_ERROR;
    }

// If FASTLINK-based method, run makeped command for each marker

    if (!MCarlo && !MMSibs && run_makeped(interp) == TCL_ERROR)
        return TCL_ERROR;

    unlink("ibdprep.out");
    unlink("ibdprep.loc");
    unlink("ibdprep.mrk");

    return TCL_OK;
}

bool Marker::get_stats (const char **errmsg)
{
    int i;
    char rec[1024], fname[1024];

    FILE *fp = fopen("marker.info", "r");
    if (!fp) {
        *errmsg = "Cannot open marker.info";
        return false;
    }

    if (!fgets(rec, sizeof(rec), fp) ||
            sscanf(rec, "%s", fname) != 1) {
        fclose(fp);
        *errmsg = "Read error on marker.info, line 1";
        return false;
    }

    if (_nloci) {
        for (i = 0; i < _nloci; i++) {
            if (!fgets(rec, sizeof(rec), fp) ||
                    sscanf(rec, "%s %d %d", _mrk[i]->name,
                           &_mrk[i]->ntyped, &_mrk[i]->nfoutyp) != 3) {
                fclose(fp);
                sprintf(error_message,
                        "Read error on marker.info, line %d", i + 2);
                *errmsg = error_message;
                return false;
            }
        }
    }
    else {
        while (fgets(rec, sizeof(rec), fp)) {
            if (_nloci == MAXLOC) {
                sprintf(error_message,
                        "ERROR: More than the maximum of %d loci in marker.info",
                        MAXLOC);
                *errmsg = error_message;
                return false;
            }
            _mrk[_nloci] = (struct Mrk *) malloc(sizeof(struct Mrk));
            if (!_mrk[_nloci]) {
                fclose(fp);
                *errmsg = "Read error on marker.info: not enough memory";
                return false;
            }
            if (sscanf(rec, "%s %d %d", _mrk[_nloci]->name,
                       &_mrk[_nloci]->ntyped, &_mrk[_nloci]->nfoutyp) != 3) {
                fclose(fp);
                sprintf(error_message,
                        "Read error on marker.info, line %d", _nloci + 2);
                *errmsg = error_message;
                return false;
            }
            _nloci++;
        }
    }

    fclose(fp);
    *errmsg = 0;
    return true;
}

int Marker::run_makeped (Tcl_Interp *interp)
{
    int i;
    FILE *fp;
    char infile[1024];
    char load_cmd[1024];
    const char *makeped_cmd =
         "exec makeped <makeped.cmd >&/dev/null";

    for (i = 0; i < _nloci; i++) {
        sprintf(infile, "d_%s/makeped.cmd", _mrk[i]->name);
        fp = fopen(infile, "r");
        if (!fp)
            continue;
        fclose(fp);

        sprintf(load_cmd, "cd d_%s", _mrk[i]->name);
        if (Solar_Eval(interp, load_cmd) == TCL_ERROR) {
            RESULT_LIT ("run_makeped: cd command failed");
            return TCL_ERROR;
        }

        if (Solar_Eval(interp, makeped_cmd) == TCL_ERROR) {
            sprintf(load_cmd, "cd ..");
            if (Solar_Eval(interp, load_cmd) == TCL_ERROR) {
                RESULT_LIT ("run_makeped: cd command failed");
                return TCL_ERROR;
            }
            RESULT_LIT ("Program makeped failed");
            return TCL_ERROR;
        }

        sprintf(load_cmd, "cd ..");
        if (Solar_Eval(interp, load_cmd) == TCL_ERROR) {
            RESULT_LIT ("run_makeped: cd command failed");
            return TCL_ERROR;
        }
    }

    return TCL_OK;
}

int Marker::get_marker (const char *name)
{
    int i;

    for (i = 0; i < _nloci; i++)
        if (!StringCmp(_mrk[i]->name, name, case_ins))
            return i;

    return -1;
}

bool Marker::unload (bool nosave, Tcl_Interp *interp)
{
    char cmd[1024];
    int i, rc;

//  unless nosave option is in effect, abort the unload if any of
//  the markers have unsaved MLE allele freqs

    if (!nosave) {
        for (i = 0; i < _nloci; i++) {
            int m = currentFreq->get_marker(_mrk[i]->name);
            if (m >= 0 && currentFreq->mle_status(m) == 'c')
                return false;
        }
    }

    printf("Unloading current marker data ...\n");

    for (i = 0; i < _nloci; i++) {
        sprintf(cmd, "exec rm -r d_%s", _mrk[i]->name);
        rc = Solar_Eval(interp, cmd);

        int m = currentFreq->get_marker(_mrk[i]->name);
        if (m >= 0) {

//  if this marker's allele freqs were computed by load marker,
//  remove the freqs from freq.info

            if (currentFreq->whence(m) == 'm')
                currentFreq->remove_marker(_mrk[i]->name);

//  if this marker's allele freqs were loaded from a freq file
//  that was subsequently unloaded while genotype data was still
//  loaded, remove the freqs from freq.info now

            else if (currentFreq->whence(m) == 'u')
                currentFreq->remove_marker(_mrk[i]->name);

//  if this marker's allele freqs were loaded from a freq file
//  and subsequently maximized, mark the MLE freqs as "old"

            else if (currentFreq->mle_status(m) == 'c' ||
                     currentFreq->mle_status(m) == 's')
                currentFreq->set_mle_status(m, 'o');
        }
    }

    if (!currentFreq->nloci()) {
        currentFreq->unload(interp);
        delete currentFreq;
    }
    else
        rc = currentFreq->write_info(interp);	// update freq.info

    unlink("marker.info");

    return true;
}

int Marker::discrep (const char *mrkname, Tcl_Interp *interp)
{
    char discrep_cmd[1024], errmsg[1024];

    int mrk = get_marker(mrkname);
    if (mrk < 0) {
        sprintf(errmsg, "%s: No such marker.", mrkname);
        RESULT_BUF (errmsg);
        return TCL_ERROR;
    }

    printf("Running allfreq for marker %s ... ", _mrk[mrk]->name);
    fflush(stdout);

    sprintf(discrep_cmd, "cd d_%s", _mrk[mrk]->name);
    if (Solar_Eval(interp, discrep_cmd) == TCL_ERROR) {
        sprintf(errmsg, "\nCannot change to directory d_%s.", _mrk[mrk]->name);
        RESULT_BUF (errmsg);
        return TCL_ERROR;
    }

    unlink("allfreq.out");

    sprintf(discrep_cmd, "exec allfreq n 1");
    if (Solar_Eval(interp, discrep_cmd) == TCL_ERROR) {
        sprintf(discrep_cmd, "cd ..");
        if (Solar_Eval(interp, discrep_cmd) == TCL_ERROR) {
            RESULT_LIT ("\nCannot return to current working directory.");
            return TCL_ERROR;
        }

        FILE *fp;
        char file[1024];
        sprintf(file, "d_%s/allfreq.out", _mrk[mrk]->name);

        fp = fopen(file, "r");
        if (!fp)
            RESULT_LIT ("\nallfreq: Not enough memory.");
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

    unlink("allfreq.out");

    sprintf(discrep_cmd, "cd ..");
    if (Solar_Eval(interp, discrep_cmd) == TCL_ERROR) {
        RESULT_LIT ("\nCannot return to current working directory.");
        return TCL_ERROR;
    }

    printf("\n");
    fflush(stdout);

    return TCL_OK;
}

void Marker::show_names (Tcl_Interp *interp)
{
    int i;
    char *buf = (char *) malloc(_nloci*(MMRKNM+1));

    strcpy(buf, _mrk[0]->name);
    for (i = 1; i < _nloci; i++) {
        strcat(buf, " ");
        strcat(buf, _mrk[i]->name);
    }

    RESULT_BUF (buf);
    free(buf);
}

void Marker::show (char *buf, Tcl_Interp *interp)
{
    int i, mrk, nall;
    char buf1[1024], cmd[100];
    float nind = (float) currentPed->nind();
    float nfou = (float) currentPed->nfou();

    bool hwe_done = false;
    if (!strcmp(buf, "")) {
        for (i = 0; i < _nloci; i++) {
            if (currentFreq->chi2_hwe(currentFreq->get_marker(_mrk[i]->name)) >= 0)
            {
                hwe_done = true;
                break;
            }
        }
        sprintf(buf, "\nmarker data file: %s\n\n", _filename);
        strcat(buf, "marker           #typed  %typed  #foutyp  %foutyp");
        if (hwe_done) strcat(buf, "  HWE p-val");
        strcat(buf, "\n");
        strcat(buf, "---------------  ------  ------  -------  -------");
        if (hwe_done) strcat(buf, "  ---------");
        strcat(buf, "\n");
        for (i = 0; i < _nloci; i++) {
            sprintf(buf1, "%-15s  %6d  %6.1f  %7d  %7.1f", _mrk[i]->name,
                    _mrk[i]->ntyped, _mrk[i]->ntyped*100/nind,
                    _mrk[i]->nfoutyp, _mrk[i]->nfoutyp*100/nfou);
            strcat(buf, buf1);
            mrk = currentFreq->get_marker(_mrk[i]->name);
            if (currentFreq->chi2_hwe(mrk) >= 0) {
                nall = currentFreq->nall(mrk);
                if (nall > 1) {
                    sprintf(cmd,
                            "chi -number %g %d", currentFreq->chi2_hwe(mrk),
                            nall*(nall + 1)/2 - nall);
                    Solar_Eval(interp, cmd);
                    sprintf(buf1, "  %s", interp->result);
                    strcat(buf, buf1);
                }
            }
            strcat(buf, "\n");
        }
    }

    else {
        for (i = 0; i < _nloci; i++) {
            if (!StringCmp(buf, _mrk[i]->name, case_ins))
                break;
        }
        if (i == _nloci)
            strcat(buf, ": No such marker.\n");
        else {
            sprintf(buf1, "%-15s  %6d  %6.1f  %7d  %7.1f", _mrk[i]->name,
                    _mrk[i]->ntyped, _mrk[i]->ntyped*100/nind,
                    _mrk[i]->nfoutyp, _mrk[i]->nfoutyp*100/nfou);
            strcpy(buf, buf1);
            mrk = currentFreq->get_marker(_mrk[i]->name);
            if (currentFreq->chi2_hwe(mrk) >= 0) {
                nall = currentFreq->nall(mrk);
                if (nall > 1) {
                    sprintf(cmd,
                            "chi -number %g %d", currentFreq->chi2_hwe(mrk),
                            nall*(nall + 1)/2 - nall);
                    Solar_Eval(interp, cmd);
                    sprintf(buf1, "  %s", interp->result);
                    strcat(buf, buf1);
                }
            }
            strcat(buf, "\n");
        }
    }
}

void Marker::delete_Tfile ()
{
    if (Tfile) {
        delete Tfile;
        Tfile = 0;
    }
}
