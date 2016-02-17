/*
 * map.cc implements the map command and class
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

Map *currentMap = 0;

extern "C" int MapCmd (ClientData clientData, Tcl_Interp *interp,
		  int argc, char *argv[])
{
    if (argc == 2 && !StringCmp ("help", argv[1], case_ins)) {
        return Solar_Eval (interp, "help map");
    }

    else if (argc == 2 && !StringCmp ("names", argv[1], case_ins)) {

    // display marker names
        if (!currentMap) {
            RESULT_LIT ("Map data have not been loaded.");
            return TCL_ERROR;
        }

        currentMap->show_names(interp);

        return TCL_OK;
    }

    else if (argc == 2 && !StringCmp ("fname", argv[1], case_ins)) {

    // return map file name
        if (!currentMap) {
            RESULT_LIT ("Map data have not been loaded.");
            return TCL_ERROR;
        }

        char buf[1024];
        sprintf(buf, "%s", currentMap->filename());
        RESULT_BUF (buf);
        return TCL_OK;
    }

    else if (argc == 2 && !StringCmp ("chrnum", argv[1], case_ins)) {

    // return chromosome identifier
        if (!currentMap) {
            RESULT_LIT ("Map data have not been loaded.");
            return TCL_ERROR;
        }

        char buf[1024];
        sprintf(buf, "%s", currentMap->chrnum());
        RESULT_BUF (buf);
        return TCL_OK;
    }

    else if (argc == 2 && !StringCmp ("nloci", argv[1], case_ins)) {

    // return number of marker loci
        if (!currentMap) {
            RESULT_LIT ("Map data have not been loaded.");
            return TCL_ERROR;
        }

        char buf[1024];
        sprintf(buf, "%d", currentMap->nloci());
        RESULT_BUF (buf);
        return TCL_OK;
    }

    else if (argc == 2 && !StringCmp ("func", argv[1], case_ins)) {

    // return mapping function
        if (!currentMap) {
            RESULT_LIT ("Map data have not been loaded.");
            return TCL_ERROR;
        }

        char buf[1024];
        sprintf(buf, "%c", currentMap->mapfunc());
        RESULT_BUF (buf);
        return TCL_OK;
    }

    else if (argc >= 2 && !StringCmp ("unload", argv[1], case_ins)) {

    // unload map data
        if (currentMap)
            delete currentMap;
        return TCL_OK;
    }

    else if (argc >= 2 && !StringCmp ("load", argv[1], case_ins)) {
        char mapfunc;
        if (argc == 3)
            mapfunc = ' ';
        else if (argc == 4 && !StringCmp ("-kosambi", argv[2], case_ins))
            mapfunc = 'k';
        else if (argc == 4 && !StringCmp ("-haldane", argv[2], case_ins))
            mapfunc = 'h';
        else if (argc == 4 && !StringCmp ("-basepair", argv[2], case_ins))
            mapfunc = 'b';
        else {
            RESULT_LIT (
            "Usage: map load [-haldane | -kosambi | -basepair] <filename>");
            return TCL_ERROR;
        }

    // load new map data
        if (currentMap)
            delete currentMap;
        if (argc == 3)
            currentMap = new Map (argv[2]);
        else
            currentMap = new Map (argv[3]);

        if (currentMap->load(interp, mapfunc) == TCL_ERROR) {
            delete currentMap;
            return TCL_ERROR;
        }

        return TCL_OK;
    }

    else if (argc >= 2 && !StringCmp ("show", argv[1], case_ins)) {
        if (argc > 3) {
            RESULT_LIT ("Usage: map show [<marker>]");
            return TCL_ERROR;
        }

    // display map info
        if (!currentMap) {
            RESULT_LIT ("Map data have not been loaded.");
            return TCL_ERROR;
        }

        char *buf = (char *) malloc(30*(currentMap->nloci()+2)+150);
        if (argc == 2)
            strcpy(buf, "");
        else
            strcpy(buf, argv[2]);

        currentMap->show(buf);
        RESULT_BUF (buf);
        free(buf);
        return TCL_OK;
    }

    else if (argc == 3 && !StringCmp ("locn", argv[1], case_ins)) {

    // return marker location
        if (!currentMap) {
            RESULT_LIT ("Map data have not been loaded.");
            return TCL_ERROR;
        }

        int mrk = currentMap->get_marker(argv[2]);

        char buf[10000];
        if (mrk < 0) {
            sprintf(buf, "%s: No such marker.", argv[2]);
            RESULT_BUF (buf);
            return TCL_ERROR;
        }
        else {
            if (currentMap->mapfunc() == 'b')
                sprintf(buf, "%8d", (int) currentMap->mrklocn(mrk));
            else
                sprintf(buf, "%8.4f", currentMap->mrklocn(mrk));
            RESULT_BUF (buf);
            return TCL_OK;
        }
    }

    RESULT_LIT ("Invalid map command");
    return TCL_ERROR;
}

int Map::load (Tcl_Interp *interp, char mapfunc)
{
    _mapfunc = mapfunc;

    FILE *mapfp = fopen("map.file.tmp", "r");
    if (!mapfp) {
        RESULT_LIT ("Cannot open map.file.tmp");
        return TCL_ERROR;
    }

    char rec[1024], *recp;
    char errmsg[1024];
    double locn;
    int i;

    if (!fgets(rec, sizeof(rec), mapfp)) {
        fclose(mapfp);
        RESULT_LIT ("Map data file is empty.");
        return TCL_ERROR;
    }

    if (!(recp = strtok(rec, " \t\n")) || sscanf(recp, "%s", _chrnum) != 1)
    {
        fclose(mapfp);
        RESULT_LIT ("Invalid chromosome number.");
        return TCL_ERROR;
    }

    if (recp = strtok(NULL, " \t\n")) {
        if (!StringCmp("kosambi", recp, case_ins)) {
            if (_mapfunc == ' ') _mapfunc = 'k';
            if (_mapfunc == 'h') {
                printf(
"Map file specifies Kosambi mapping, but Haldane mapping will be assumed.\n");
            }
            if (_mapfunc == 'b') {
                printf(
"Map file specifies Kosambi mapping, but basepair locations will be assumed.\n");
            }
        }
        else if (!StringCmp("haldane", recp, case_ins)) {
            if (_mapfunc == ' ') _mapfunc = 'h';
            if (_mapfunc == 'k') {
                printf(
"Map file specifies Haldane mapping, but Kosambi mapping will be assumed.\n");
            }
            if (_mapfunc == 'b') {
                printf(
"Map file specifies Haldane mapping, but basepair locations will be assumed.\n");
            }
        }
        else if (!StringCmp("basepair", recp, case_ins)) {
            if (_mapfunc == ' ') _mapfunc = 'b';
            if (_mapfunc == 'k') {
                printf(
"Map file specifies basepair locations, but Kosambi mapping will be assumed.\n");
            }
            if (_mapfunc == 'h') {
                printf(
"Map file specifies basepair locations, but Haldane mapping will be assumed.\n");
            }
        }
        else {
            fclose(mapfp);
            RESULT_LIT ("Invalid mapping function.");
            return TCL_ERROR;
        }
    }

    if (_mapfunc == ' ')
        _mapfunc = 'k';

    _nloci = 0;
    while (fgets(rec, sizeof(rec), mapfp)) {
        if (_nloci == MAXLOC) {
            sprintf(errmsg,
                    "ERROR: More than the maximum of %d loci in map data file",
                    MAXLOC);
            RESULT_BUF (errmsg);
            return TCL_ERROR;
        }

        if (!(recp = strtok(rec, " \t\n"))) {
            fclose(mapfp);
            sprintf(errmsg, "Invalid record, line %d of map data file",
                    _nloci + 2);
            RESULT_BUF (errmsg);
            return TCL_ERROR;
        }

        for (i = 0; i < _nloci; i++) {
            if (!strcmp(recp, _locus[i]->name)) {
                fclose(mapfp);
                sprintf(errmsg,
                        "Marker %s appears more than once in map data file",
                        recp);
                RESULT_BUF (errmsg);
                return TCL_ERROR;
            }
        }

        _locus[_nloci] = (struct LocusMap *) malloc(sizeof(struct LocusMap));
        strcpy(_locus[_nloci]->name, recp);

        if (!(recp = strtok(NULL, " \t\n")) || sscanf(recp, "%lf", &locn) != 1) {
            fclose(mapfp);
            sprintf(errmsg, "Invalid record, line %d of map data file",
                    _nloci + 2);
            RESULT_BUF (errmsg);
            return TCL_ERROR;
        }

        if (_nloci && locn < mrklocn(_nloci-1)) {
            fclose(mapfp);
            RESULT_LIT ("Markers must be ordered by chromosome location.");
            return TCL_ERROR;
        }

        _locus[_nloci]->locn = locn;
        _nloci++;
    }
 
    if (!_nloci) {
        fclose(mapfp);
        RESULT_LIT ("Can't load empty map data file!");
        return TCL_ERROR;
    }

    fclose(mapfp);

    return TCL_OK;
}
 
Map::~Map ()
{
    int i;

    for (i = 0; i < _nloci; i++)
        free(_locus[i]);

    currentMap = 0;
}

char *Map::show (char *buf)
{
    int i;
    char buf1[1024];
 
    if (!strcmp(buf, "")) {
        sprintf(buf, "\nmap data file: %s\n\n", _filename);
        sprintf(buf1, "chromosome: %s", _chrnum);
        strcat(buf, buf1);
        if (_mapfunc == 'b')
            sprintf(buf1, "");
        else if (_mapfunc == 'k')
            sprintf(buf1, "\tmapping function: Kosambi");
        else if (_mapfunc == 'h')
            sprintf(buf1, "\tmapping function: Haldane");
        strcat(buf, buf1);
        sprintf(buf1, "\n\n");
        strcat(buf, buf1);
        if (_mapfunc == 'b')
            strcat(buf, "marker         \tlocn(bp)\n");
        else
            strcat(buf, "marker         \tlocn(cM)\n");
        strcat(buf, "--------------\t--------\t\n");
        for (i = 0; i < _nloci; i++) {
            if (_mapfunc == 'b')
                sprintf(buf1, "%-15s\t%8d\n", mrkname(i), (int) mrklocn(i));
            else
                sprintf(buf1, "%-15s\t%8.4f\n", mrkname(i), mrklocn(i));
            strcat(buf, buf1);
        }
    }    

    else {
        for (i = 0; i < _nloci; i++) {
            if (!StringCmp(buf, mrkname(i), case_ins))
                break;
        }
        if (i == _nloci)
            strcat(buf, ": No such marker.\n");
        else {
            if (_mapfunc == 'b')
                sprintf(buf, "%-15s\t%8d\n", mrkname(i), (int) mrklocn(i));
            else
                sprintf(buf, "%-15s\t%8.4f\n", mrkname(i), mrklocn(i));
        }
    }    

    return buf;
}

void Map::show_names (Tcl_Interp *interp)
{
    int i;
    char *buf = (char *) malloc(_nloci*(MMRKNM+1));

    strcpy(buf, _locus[0]->name);
    for (i = 1; i < _nloci; i++) {
        strcat(buf, " ");
        strcat(buf, _locus[i]->name);
    }

    RESULT_BUF (buf);
    free(buf);
}

int Map::get_marker (const char *name)
{
    int i;

    for (i = 0; i < _nloci; i++)
        if (!StringCmp(_locus[i]->name, name, case_ins))
            return i;

    return -1;
}
