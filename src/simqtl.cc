/*
 * Written by Thomas Dyer beginning on October 18, 1999
 * Copyright (c) 1999 Southwest Foundation for Biomedical Research
 */

#include <stdlib.h>
#include <math.h>
#include "solar.h"
#include "pipeback.h"

#define MXLOC  500
#define MXTRT  500
#define MXCOV  500
#define MIDLEN 20

void prtmat (double *x, int n)
{
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf(" %g", x[j*n+i]);
        }
        printf("\n");
    }
}

struct SimLocPars
{
    int nall;
    double *freq;
    SimLocPars()   {nall = 0; freq = 0;}
    ~SimLocPars()  {if (freq) delete[] freq;}
};

struct SimTrtPars {
    double *mean;
    double sdev;
    double h2r;
    double *cmean;
    int nbeta;
    double *beta[MXCOV+3];
    double *rg;
    double *re;
    SimTrtPars(int n);
    ~SimTrtPars();
};

struct SimPars {
    int nloc;
    SimLocPars *loc[MXLOC];
    int ngen;
    int ntrt;
    int ncov;
    SimTrtPars *trt[MXTRT];
    int nmrk;
    int nmall;
    double *mfreq;
    double theta;
    SimPars();
    ~SimPars();
    void get_pars(FILE*);
    void prt_pars(FILE*);
};

struct Ego {
    char famid[MIDLEN+1];
    char id[MIDLEN+1];
    int ibdid;
    int fibdid;
    int mibdid;
    int sex;
    int mztwin;
    double age;
    double *cov;
    bool *missing;
    int hap[2];
    Ego()   {cov = 0;}
    ~Ego()  {if (cov) delete[] cov; if (missing) delete[] missing;}
};

static int do_sim (Tcl_Interp*, bool, bool);
static int get_ego (Tcl_Interp*, char**, int, int, int, int, bool, SimPars*,
                    Ego*);
static int sim_ped (Tcl_Interp*, FILE*, bool, FILE*, FILE*, bool, int, int,
                    int, int, int*, SimPars*, int, Ego*, double*, double*,
                    FILE*);
static void factor (double*, int, int*, int);
static int kincoef (double*, int*, int*, int*, int, int, FILE*);
static void simva (double*, double*, double**, int, int);
static void simve (double*, double**, int, int);
static void dropgenel (int*, int*, int*, int, double*, double*, double,
                       double, int*, int*, int, int);
static void dropgene (int*, int*, int*, int, double*, int*, int*, int);
static double gasdev (void);
static char *fp2str (double);

extern "C" void dppfa_ (double*, int*, int*);
extern "C" void eigstruc_ (double*, int*, int*);

extern bool random_number_generator_seeded;

extern "C" int SimqtlCmd (ClientData clientData, Tcl_Interp *interp,
                          int argc, char *argv[])
{
    if (argc == 2 && !StringCmp ("help", argv[1], case_ins))
    {
        return Solar_Eval (interp, "help simqtl");
    }

    if (argc == 2 && !StringCmp ("-model", argv[1], case_ins))
    {
        SimPars pars;

        FILE *fp = fopen("simqtl.par", "r");
        if (!fp) {
            RESULT_LIT ("Cannot open parameter file simqtl.par");
            return TCL_ERROR;
        }

        try {
            pars.get_pars(fp);
            fclose(fp);
        }
        catch (Safe_Error_Return &ser) {
            fclose(fp);
            RESULT_BUF (ser.message());
            return TCL_ERROR;
        }

        char tmpfile[10];
        strcpy(tmpfile, "tmpXXXXXX");
        int tmpfd = mkstemp(tmpfile);
        if (tmpfd == -1)
            return TCL_ERROR;

        FILE *tmpfp = fdopen(tmpfd, "w");
        if (tmpfp) {
            pars.prt_pars(tmpfp);
            fclose(tmpfp);
            char more_cmd[1024];
            sprintf(more_cmd, "exec >&@stdout <@stdin more %s", tmpfile);
            if (Solar_Eval(interp, more_cmd) == TCL_ERROR)
                return TCL_ERROR;
            remove(tmpfile);
        }
        return TCL_OK;
    }

    if (argc == 1)
    {
        return do_sim(interp, false, false);
    }

    if (argc == 2 && !StringCmp ("-inform", argv[1], case_ins))
    {
        return do_sim(interp, true, false);
    }

    if (!StringCmp ("-seed", argv[1], case_ins) &&
        (argc == 3 || argc == 4 && !StringCmp ("-inform", argv[3], case_ins)))
    {
        long seed;
        if (sscanf(argv[2], "%ld", &seed) == 1)
        {
            char cmd[1024];
            sprintf(cmd, "drand %d\n", seed);
            Solar_Eval(interp, cmd);
            if (argc == 4)
                return do_sim(interp, true, false);
            else
                return do_sim(interp, false, false);
        }
        RESULT_LIT ("The seed must be an integer");
        return TCL_ERROR;
    }

    if (!StringCmp ("-gfile", argv[1], case_ins) &&
        (argc == 2 || argc == 3 && !StringCmp ("-inform", argv[2], case_ins)))
    {
        if (argc == 3)
            return do_sim(interp, true, true);
        else
            return do_sim(interp, false, true);
    }

    if (!StringCmp ("-gfile", argv[1], case_ins) &&
        !StringCmp ("-seed", argv[2], case_ins) &&
        (argc == 4 || argc == 5 && !StringCmp ("-inform", argv[4], case_ins)))
    {
        long seed;
        if (sscanf(argv[3], "%ld", &seed) == 1)
        {
            char cmd[1024];
            sprintf(cmd, "drand %d\n", seed);
            Solar_Eval(interp, cmd);
            if (argc == 5)
                return do_sim(interp, true, true);
            else
                return do_sim(interp, false, true);
        }
        RESULT_LIT ("The seed must be an integer");
        return TCL_ERROR;
    }

    RESULT_LIT ("Invalid simqtl command");
    return TCL_ERROR;
}

SimTrtPars::SimTrtPars (int n)
{
    int i;
    nbeta = n;
    for (i = 0; i < nbeta; i++)
        beta[i] = 0;
    mean = 0;
    sdev = 0;
    h2r = 0;
    cmean = 0;
    rg = 0;
    re = 0;
}

SimTrtPars::~SimTrtPars ()
{
    int i;
    if (mean) delete[] mean;
    if (cmean) delete[] cmean;
    if (beta) {
        for (i = 0; i < nbeta; i++)
            if (beta[i]) delete[] beta[i];
    }
    if (rg) delete[] rg;
    if (re) delete[] re;
}

SimPars::SimPars ()
{
    int i;
    for (i = 0; i < MXLOC; i++)
        loc[i] = 0;
    for (i = 0; i < MXTRT; i++)
        trt[i] = 0;
    nloc = 0;
    ngen = 0;
    ntrt = 0;
    ncov = 0;
    nmrk = 0;
    nmall = 0;
    mfreq = 0;
    theta = 0;
}

SimPars::~SimPars ()
{
    int i;
    for (i = 0; i < nloc; i++)
        if (loc[i]) delete loc[i];
    for (i = 0; i < ntrt; i++)
        if (trt[i]) delete trt[i];
    if (nmrk) delete[] mfreq;
}

int do_sim (Tcl_Interp *interp, bool inform, bool gfile)
{
    struct storage {
        double *gmat, *emat;
        Ego *ego;
        storage() {}
        ~storage() {
            if (ego) delete[] ego;
            if (emat) delete[] emat;
            if (gmat) delete[] gmat;
        }
    } s;

    int i, j;
    SimPars pars;
    FILE *fp;

    if (!random_number_generator_seeded)
        Solar_Eval(interp, "drand 0");

    fp = fopen("simqtl.par", "r");
    if (!fp) {
        RESULT_LIT ("Cannot open parameter file simqtl.par");
        return TCL_ERROR;
    }

    try {
        pars.get_pars(fp);
        fclose(fp);
    }
    catch (Safe_Error_Return &ser) {
        fclose(fp);
        RESULT_BUF (ser.message());
        return TCL_ERROR;
    }

    if (pars.nloc > 1 && !gfile) {
        RESULT_LIT ("A genotype file is required for a multilocus simulation");
        return TCL_ERROR;
    }

    try { s.gmat = new double[pars.ntrt*pars.ntrt]; }
    catch (...) {
        RESULT_LIT ("Out of memory");
        return TCL_ERROR;
    }

    try { s.emat = new double[pars.ntrt*pars.ntrt]; }
    catch (...) {
        RESULT_LIT ("Out of memory");
        return TCL_ERROR;
    }

    for (i = 0; i < pars.ntrt; i++) {
        s.gmat[i*pars.ntrt+i] = pars.trt[i]->h2r * pow(pars.trt[i]->sdev, 2.);
        s.emat[i*pars.ntrt+i] =
                        (1 - pars.trt[i]->h2r) * pow(pars.trt[i]->sdev, 2.);
        for (j = i + 1; j < pars.ntrt; j++) {
            s.gmat[i*pars.ntrt+j] = s.gmat[j*pars.ntrt+i] =
                pars.trt[j]->rg[i] *
                sqrt(pars.trt[i]->h2r * pow(pars.trt[i]->sdev, 2.)
                     * pars.trt[j]->h2r * pow(pars.trt[j]->sdev, 2.));
            s.emat[i*pars.ntrt+j] = s.emat[j*pars.ntrt+i] =
                pars.trt[j]->re[i] *
                sqrt((1 - pars.trt[i]->h2r) * pow(pars.trt[i]->sdev, 2.)
                     * (1 - pars.trt[j]->h2r) * pow(pars.trt[j]->sdev, 2.));
        }
    }

    int info;
    factor(s.gmat, pars.ntrt, &info, 0);
    if (info) {
        if (info == -1) {
            RESULT_LIT ("Out of memory");
            return TCL_ERROR;
        }
        char mbuf[1024];
        sprintf (mbuf,
            "Cannot factor genetic correlation matrix, info = %d", info);
        RESULT_BUF (mbuf);
        return TCL_ERROR;
    }

    factor(s.emat, pars.ntrt, &info, 0);
    if (info) {
        if (info == -1) {
            RESULT_LIT ("Out of memory");
            return TCL_ERROR;
        }
        char mbuf[1024];
        sprintf (mbuf,
            "Cannot factor environmental correlation matrix, info = %d", info);
        RESULT_BUF (mbuf);
        return TCL_ERROR;
    }

    TableFile *Tfile;
    const char *errmsg = 0;
    const char **names;
    int count, ncov, need_famid, age_avail;

    Tfile = TableFile::open ("simqtl.dat", &errmsg);
    if (errmsg) {
        if (!strcmp(errmsg, "File not found"))
            RESULT_LIT ("Cannot open data file simqtl.dat");
        else
            Solar_AppendResult2(interp, "do_sim: ", (char*) errmsg, NULL);
        delete Tfile;
        return TCL_ERROR;
    }

    Tfile->start_setup(&errmsg);
    names = Tfile->names(&count, &errmsg);

    need_famid = 0;
    if (Tfile->test_name ("FAMID", &errmsg)) {
        Tfile->setup ("FAMID", &errmsg);
        need_famid = 1;
    }
    Tfile->setup ("ID", &errmsg);

    Tfile->setup ("PEDNO", &errmsg);
    Tfile->setup ("IBDID", &errmsg);
    Tfile->setup ("FIBDID", &errmsg);
    Tfile->setup ("MIBDID", &errmsg);
    Tfile->setup ("SEX", &errmsg);
    Tfile->setup ("MZTWIN", &errmsg);

    age_avail = 0;
    if (Tfile->test_name ("AGE", &errmsg)) {
        Tfile->setup ("AGE", &errmsg);
        age_avail = 1;
    }

    if (gfile)
        ncov = count - (7 + age_avail + need_famid + pars.nloc);
    else
        ncov = count - (7 + age_avail + need_famid);

    if (ncov != pars.ncov - 2) {
        RESULT_LIT ("do_sim: simqtl.dat: wrong number of covariates");
        delete Tfile;
        return TCL_ERROR;
    }

    if (gfile) {
        for (i = 0; i < ncov + pars.nloc; i++)
            Tfile->setup (names[count-(pars.nloc+ncov)+i], &errmsg);
    }
    else {
        for (i = 0; i < ncov; i++)
            Tfile->setup (names[count-ncov+i], &errmsg);
    }

    if (errmsg) {
        char mbuf[1024];
        sprintf (mbuf, "do_sim: simqtl.dat: %s", errmsg);
        RESULT_BUF (mbuf);
        delete Tfile;
        return TCL_ERROR;
    }

    Tfile->rewind(&errmsg);
    char **record;
    int ped, nrec = 0, nped = 0;
    int nind = 0, mxind = 0;
    while (1) {
        record = Tfile->get(&errmsg);
        if (errmsg && !strcmp("EOF", errmsg)) {
            if (nind > mxind) mxind = nind;
            break;
        }
        if (errmsg) {
            char mbuf[1024];
            sprintf (mbuf, "do_sim: simqtl.dat: %s", errmsg);
            RESULT_BUF (mbuf);
            delete Tfile;
            return TCL_ERROR;
        }
        if (sscanf(record[1+need_famid], "%d", &ped) != 1) {
            RESULT_LIT ("do_sim: simqtl.dat: Pedigree numbers must be integer");
            delete Tfile;
            return TCL_ERROR;
        }
        if (ped != nped) {
            nped = ped;
            if (nind > mxind) mxind = nind;
            nind = 0;
        }
        nind++;
        nrec++;
    }

    try { s.ego = new Ego[mxind]; }
    catch (...) {
        RESULT_LIT ("Out of memory");
        delete Tfile;
        return TCL_ERROR;
    }

    if (ncov) {
        for (i = 0; i < mxind; i++) {
            try { s.ego[i].cov = new double[ncov]; }
            catch (...) {
                RESULT_LIT ("Out of memory");
                delete Tfile;
                return TCL_ERROR;
            }
        }
    }

    for (i = 0; i < mxind; i++) {
        try { s.ego[i].missing = new bool[ncov+1]; }
        catch (...) {
            RESULT_LIT ("Out of memory");
            delete Tfile;
            return TCL_ERROR;
        }
    }

    Tfile->rewind(&errmsg);
    unlink("simqtl.phn");
    unlink("simqtl.qtl");
    unlink("simqtl.mrk");

    FILE *phnfp = fopen("simqtl.phn", "w");
    if (!phnfp) {
        RESULT_LIT ("Cannot open output file simqtl.phn");
        delete Tfile;
        return TCL_ERROR;
    }

    if (need_famid)
        fprintf(phnfp, "%s,", Field::Map("FAMID"));
    fprintf(phnfp, "%s", Field::Map("ID"));
    if (age_avail)
        fprintf(phnfp, ",AGE");
    for (i = 0; i < ncov; i++)
        fprintf(phnfp, ",%s", names[count-ncov+i]);
    if (pars.ntrt == 1)
        fprintf(phnfp, ",SIMQT");
    else {
        for (i = 0; i < pars.ntrt; i++)
            fprintf(phnfp, ",SIMQT%d", i+1);
    }
    fprintf(phnfp, "\n");

    FILE *qtlfp = fopen("simqtl.qtl", "w");
    if (!qtlfp) {
        RESULT_LIT ("Cannot open output file simqtl.qtl");
        fclose(phnfp);
        delete Tfile;
        return TCL_ERROR;
    }

    if (need_famid)
        fprintf(qtlfp, "%s,", Field::Map("FAMID"));
    fprintf(qtlfp, "%s,QTL\n", Field::Map("ID"));

    FILE *mrkfp;
    if (pars.nmrk) {
        mrkfp = fopen("simqtl.mrk", "w");
        if (!fp) {
            RESULT_LIT ("Cannot open output file simqtl.mrk");
            fclose(qtlfp);
            fclose(phnfp);
            delete Tfile;
            return TCL_ERROR;
        }

        if (!inform) {
            if (need_famid)
                fprintf(mrkfp, "%s,", Field::Map("FAMID"));
            fprintf(mrkfp, "%s,SIMMRK\n", Field::Map("ID"));
        }
    }

    FILE *kinfp = fopen("phi2.gz", "r");
    if (!kinfp) {
        RESULT_LIT ("Cannot open phi2.gz");
        if (pars.nmrk) {
            fclose(mrkfp);
        }
        fclose(qtlfp);
        fclose(phnfp);
        delete Tfile;
        return TCL_ERROR;
    }
    fclose(kinfp);

    const char *pbarg[4];
    pbarg[0] = "gunzip";
    pbarg[1] = "-c";
    pbarg[2] = "phi2";
    pbarg[3] = 0;

    kinfp = pipeback_shell_open("gunzip", pbarg);
    if (!kinfp) {
        RESULT_LIT ("Cannot uncompress phi2");
        if (pars.nmrk) {
            fclose(mrkfp);
        }
        fclose(qtlfp);
        fclose(phnfp);
        delete Tfile;
        return TCL_ERROR;
    }

    nind = 0;
    nped = 1;
    int tnind = 0;
    while (1) {
        record = Tfile->get(&errmsg);
        if (errmsg && !strcmp("EOF", errmsg)) {
            if (sim_ped(interp, phnfp, gfile, qtlfp, mrkfp, inform, need_famid,
                        age_avail, nped, nind, &tnind, &pars, ncov, s.ego,
                        s.gmat, s.emat, kinfp) == TCL_ERROR)
            {
                if (pars.nmrk) {
                    fclose(mrkfp);
                }
                fclose(qtlfp);
                fclose(phnfp);
                delete Tfile;
                return TCL_ERROR;
            }
            break;
        }

        if (errmsg) {
            char mbuf[1024];
            sprintf (mbuf, "do_sim: simqtl.dat: %s", errmsg);
            RESULT_BUF (mbuf);
            delete Tfile;
            return TCL_ERROR;
        }

        if (sscanf(record[1+need_famid], "%d", &ped) != 1) {
            RESULT_LIT ("do_sim: simqtl.dat: Pedigree numbers must be integer");
            delete Tfile;
            return TCL_ERROR;
        }

        if (ped != nped) {
            if (sim_ped(interp, phnfp, gfile, qtlfp, mrkfp, inform, need_famid,
                        age_avail, nped, nind, &tnind, &pars, ncov, s.ego,
                        s.gmat, s.emat, kinfp) == TCL_ERROR)
            {
                if (pars.nmrk) {
                    fclose(mrkfp);
                }
                fclose(qtlfp);
                fclose(phnfp);
                delete Tfile;
                return TCL_ERROR;
            }
            nped = ped;
            nind = 0;
        }

        if (get_ego(interp, record, need_famid, age_avail, nind, ncov, gfile,
                    &pars, s.ego) == TCL_ERROR)
        {
            if (pars.nmrk) {
                fclose(mrkfp);
            }
            fclose(qtlfp);
            fclose(phnfp);
            delete Tfile;
            return TCL_ERROR;
        }

        nind++;
    }

    pipeback_shell_close(kinfp);

    if (pars.nmrk) {
        fclose(mrkfp);
    }
    fclose(qtlfp);
    fclose(phnfp);
    delete Tfile;

    return TCL_OK;
}

int get_ego(Tcl_Interp *interp, char **record, int need_famid, int age_avail,
            int nind, int ncov, bool gfile, SimPars *pars, Ego *ego)
{
    if (need_famid)
        strcpy(ego[nind].famid, record[0]);

    strcpy(ego[nind].id, record[0+need_famid]);

    ego[nind].fibdid = ego[nind].mibdid = 0;
    if (*record[2+need_famid] &&
          sscanf(record[2+need_famid], "%d", &ego[nind].ibdid) != 1 ||
        *record[3+need_famid] &&
          sscanf(record[3+need_famid], "%d", &ego[nind].fibdid) != 1 ||
        *record[4+need_famid] &&
          sscanf(record[4+need_famid], "%d", &ego[nind].mibdid) != 1)
    {
        RESULT_LIT ("do_sim: simqtl.dat: IDs must be integer");
        return TCL_ERROR;
    }

    if (sscanf(record[5+need_famid], "%d", &ego[nind].sex) != 1) {
        RESULT_LIT ("do_sim: simqtl.dat: Sex must be coded (1,2)");
        return TCL_ERROR;
    }

    ego[nind].mztwin = 0;
    if (*record[6+need_famid] &&
          sscanf(record[6+need_famid], "%d", &ego[nind].mztwin) != 1)
    {
        RESULT_LIT ("do_sim: simqtl.dat: Invalid MZTWIN ID");
        return TCL_ERROR;
    }

    ego[nind].age = 0;
    ego[nind].missing[0] = false;
    if (age_avail) {
        if (!*record[7+need_famid])
            ego[nind].missing[0] = true;
        else if (sscanf(record[7+need_famid], "%lf", &ego[nind].age) != 1)
        {
            RESULT_LIT ("do_sim: simqtl.dat: Invalid age");
            return TCL_ERROR;
        }
    }

    for (int i = 0; i < ncov; i++) {
        ego[nind].cov[i] = 0;
        ego[nind].missing[i+1] = false;
        if (!*record[7+age_avail+need_famid+i])
            ego[nind].missing[i+1] = true;
        else if (sscanf(record[7+age_avail+need_famid+i], "%lf",
                        &ego[nind].cov[i]) != 1)
        {
            RESULT_LIT ("do_sim: simqtl.dat: Invalid covariate");
            return TCL_ERROR;
        }
    }

    if (gfile) {
        int a1, a2;
        ego[nind].hap[0] = ego[nind].hap[1] = 0;
        for (int i = 0; i < pars->nloc; i++) {
            if (*record[7+age_avail+need_famid+ncov+i] &&
                  sscanf(record[7+age_avail+need_famid+ncov+i],
                         "%d/%d", &a1, &a2) != 2)
            {
                RESULT_LIT ("do_sim: simqtl.dat: Invalid QTL genotype");
                return TCL_ERROR;
            }
            if (i) {
                ego[nind].hap[0] *= pars->loc[i]->nall;
                ego[nind].hap[1] *= pars->loc[i]->nall;
            }
            ego[nind].hap[0] += a1 - 1;
            ego[nind].hap[1] += a2 - 1;
        }
    }

    return TCL_OK;
}

int sim_ped (Tcl_Interp *interp, FILE *phnfp, bool gfile, FILE *qtlfp,
             FILE *mrkfp, bool inform, int need_famid, int age_avail,
             int ped, int nind, int *tnind, SimPars *pars, int ncov,
             Ego *ego, double *gmat, double *emat, FILE *kinfp)
{
    struct storage {
        int *ifa, *imo, *itwin;
        int *twinid, *twin1;
        double *amat;
        int ntrt;
        double **gdev, **edev;
        int *patgene, *matgene;
        storage() {}
        ~storage() {
            if (patgene) delete[] patgene;
            if (matgene) delete[] matgene;
            for (int i = 0; i < ntrt; i++) {
                if (edev[i]) delete[] edev[i];
                if (gdev[i]) delete[] gdev[i];
            }
            if (edev) delete[] edev;
            if (gdev) delete[] gdev;
            if (amat) delete[] amat;
            if (twin1) delete[] twin1;
            if (twinid) delete[] twinid;
            if (itwin) delete[] itwin;
            if (imo) delete[] imo;
            if (ifa) delete[] ifa;
        }
    } s;

    int i, j;
    int ntwin, nfou;

    try { s.ifa = new int[nind]; }
    catch (...) {
        RESULT_LIT ("Out of memory");
        return TCL_ERROR;
    }

    try { s.imo = new int[nind]; }
    catch (...) {
        RESULT_LIT ("Out of memory");
        return TCL_ERROR;
    }

    try { s.itwin = new int[nind]; }
    catch (...) {
        RESULT_LIT ("Out of memory");
        return TCL_ERROR;
    }

    try { s.twinid = new int[nind/2]; }
    catch (...) {
        RESULT_LIT ("Out of memory");
        return TCL_ERROR;
    }

    try { s.twin1 = new int[nind/2]; }
    catch (...) {
        RESULT_LIT ("Out of memory");
        return TCL_ERROR;
    }

    nfou = 0;
    ntwin = 0;
    for (i = 0; i < nind; i++) {
        s.ifa[i] = 0;
        if (ego[i].fibdid) {
            for (j = 0; j < nind; j++) {
                if (ego[j].ibdid == ego[i].fibdid) {
                    s.ifa[i] = j + 1;
                    break;
                }
            }
            if (j == nind) {
                RESULT_LIT ("Internal error");
                return TCL_ERROR;
            }
        }
        else
            nfou++;

        s.imo[i] = 0;
        if (ego[i].mibdid) {
            for (j = 0; j < nind; j++) {
                if (ego[j].ibdid == ego[i].mibdid) {
                    s.imo[i] = j + 1;
                    break;
                }
            }
            if (j == nind) {
                RESULT_LIT ("Internal error");
                return TCL_ERROR;
            }
        }

        s.itwin[i] = i + 1;
        if (ego[i].mztwin) {
            for (j = 0; j < ntwin; j++)
                if (s.twinid[j] == ego[i].mztwin)
                    break;
            if (j == ntwin) {
                s.twinid[ntwin] = ego[i].mztwin;
                s.twin1[ntwin] = i + 1;
                ntwin++;
            }
            else
                s.itwin[i] = s.twin1[j];
        }
    }

    try { s.amat = new double[nind*nind]; }
    catch (...) {
        RESULT_LIT ("Out of memory");
        return TCL_ERROR;
    }

    if (!kincoef(s.amat, s.ifa, s.imo, s.itwin, nind, *tnind, kinfp)) {
        RESULT_LIT ("Error encountered reading phi2.");
        return TCL_ERROR;
    }

    int info;
    factor(s.amat, nind, &info, ntwin);
    if (info) {
        if (info == -1) {
            RESULT_LIT ("Out of memory");
            return TCL_ERROR;
        }
        char mbuf[1024];
        sprintf (mbuf,
            "Cannot factor kinship coefficient matrix, info = %d", info);
        RESULT_BUF (mbuf);
        return TCL_ERROR;
    }

    try {
        s.ntrt = pars->ntrt;
        s.gdev = new double*[pars->ntrt];
        s.edev = new double*[pars->ntrt];
        for (i = 0; i < pars->ntrt; i++) {
            s.gdev[i] = new double[nind];
            s.edev[i] = new double[nind];
        }
        s.patgene = new int[2*nind];
        s.matgene = new int[2*nind];
    }
    catch (...) {
        RESULT_LIT ("Out of memory");
        return TCL_ERROR;
    }

    try {
        simva(s.amat, gmat, s.gdev, nind, pars->ntrt);
        simve(emat, s.edev, nind, pars->ntrt);
        if (gfile) {
            for (i = 0; i < nind; i++) {
                s.patgene[2*i] = ego[i].hap[0];
                s.matgene[2*i] = ego[i].hap[1];
            }
        }
        else {
            if (pars->nmrk) {
                dropgenel(s.ifa, s.imo, s.itwin, nind, pars->loc[0]->freq,
                          pars->mfreq, pars->theta, pars->theta, s.patgene,
                          s.matgene, pars->loc[0]->nall, pars->nmall);
            } else {
                dropgene(s.ifa, s.imo, s.itwin, nind, pars->loc[0]->freq,
                         s.patgene, s.matgene, pars->loc[0]->nall);
            }
        }
    }
    catch (Safe_Error_Return &ser) {
        RESULT_BUF (ser.message());
        return TCL_ERROR;
    }

    for (i = 0; i < nind; i++) {
        int qtlgen;
        if (s.patgene[2*i] >= s.matgene[2*i])
            qtlgen = s.patgene[2*i]*(s.patgene[2*i] + 1)/2 + s.matgene[2*i];
        else
            qtlgen = s.matgene[2*i]*(s.matgene[2*i] + 1)/2 + s.patgene[2*i];

        if (need_famid) {
            fprintf(phnfp, "%s,", ego[i].famid);
            fprintf(qtlfp, "%s,", ego[i].famid);
        }
        fprintf(phnfp, "%s", ego[i].id);
        fprintf(qtlfp, "%s", ego[i].id);

        bool has_data = true;
        if (age_avail) has_data = !ego[i].missing[0];
        for (j = 0; j < ncov; j++)
            if (!ego[i].cov[j]) has_data = has_data && !ego[i].missing[j+1];

        double trait;
        if (age_avail) {
            if (!ego[i].missing[0])
                fprintf(phnfp, ",%s", fp2str(ego[i].age));
            else
                fprintf(phnfp, ",");
        }
        for (j = 0; j < ncov; j++) {
            if (!ego[i].missing[j+1])
                fprintf(phnfp, ",%s", fp2str(ego[i].cov[j]));
            else
                fprintf(phnfp, ",");
        }
        if (has_data) {
            for (j = 0; j < pars->ntrt; j++) {
                if (s.itwin[i] != i + 1)
                    trait = pars->trt[j]->mean[qtlgen] +
                            s.gdev[j][s.itwin[i]-1] + s.edev[j][i];
                else
                    trait = pars->trt[j]->mean[qtlgen] +
                            s.gdev[j][i] + s.edev[j][i];
                trait += pars->trt[j]->beta[0][qtlgen] * (ego[i].sex - 1);

                if (!age_avail &&
                    (pars->trt[j]->beta[1][qtlgen] != 0 ||
                     pars->trt[j]->beta[2][qtlgen] != 0 ||
                     pars->trt[j]->beta[3][qtlgen] != 0 ||
                     pars->trt[j]->beta[4][qtlgen] != 0))
                {
                     RESULT_LIT (
"The model contains a non-zero regression coefficient for an age term,\n\
but there is no AGE field in the phenotypes file.");
                     return TCL_ERROR;
                }

                double age = ego[i].age - pars->trt[j]->cmean[1];
                if (ego[i].sex == 1) {
                    trait += pars->trt[j]->beta[1][qtlgen] * age;
                    trait += pars->trt[j]->beta[3][qtlgen] * age * age;
                }
                else {
                    trait += pars->trt[j]->beta[2][qtlgen] * age;
                    trait += pars->trt[j]->beta[4][qtlgen] * age * age;
                }
                for (int k = 0; k < ncov; k++) {
                    trait += pars->trt[j]->beta[5+k][qtlgen]
                             * (ego[i].cov[k] - pars->trt[j]->cmean[5+k]);
                }
                fprintf(phnfp, ",%s", fp2str(trait));
            }
        }
        else {
            for (j = 0; j < pars->ntrt; j++)
                fprintf(phnfp, ",");
        }
        fprintf(phnfp, "\n");

        if (s.patgene[2*i] >= s.matgene[2*i])
            fprintf(qtlfp, ",%d/%d\n", s.matgene[2*i] + 1, s.patgene[2*i] + 1);
        else
            fprintf(qtlfp, ",%d/%d\n", s.patgene[2*i] + 1, s.matgene[2*i] + 1);

        if (pars->nmrk) {
            int mrkgen, patmrk, matmrk;
            patmrk = s.patgene[2*i+1] % pars->nmall;
            matmrk = s.matgene[2*i+1] % pars->nmall;
            if (patmrk >= matmrk)
                mrkgen = patmrk*(patmrk + 1)/2 + matmrk;
            else
                mrkgen = matmrk*(matmrk + 1)/2 + patmrk;

            if (inform) {
                fprintf(mrkfp, "%d ", ego[i].ibdid);
                if (s.patgene[2*i+1] >= s.matgene[2*i+1])
                    fprintf(mrkfp, "%d/%d\n",
                            2*(*tnind) + s.matgene[2*i+1]/pars->nmall + 1,
                            2*(*tnind) + s.patgene[2*i+1]/pars->nmall + 1);
                else
                    fprintf(mrkfp, "%d/%d\n",
                            2*(*tnind) + s.patgene[2*i+1]/pars->nmall + 1,
                            2*(*tnind) + s.matgene[2*i+1]/pars->nmall + 1);
            }
            else {
                if (need_famid)
                    fprintf(mrkfp, "%s,", ego[i].famid);
                fprintf(mrkfp, "%s,", ego[i].id);
                if (patmrk >= matmrk)
                    fprintf(mrkfp, "%d/%d\n", matmrk + 1, patmrk + 1);
                else
                    fprintf(mrkfp, "%d/%d\n", patmrk + 1, matmrk + 1);
            }
        }
    }

    *tnind += nind;

    return TCL_OK;
}

void factor (double *mat, int dim, int *info, int ntwin)
{
    int n = dim;

    if (n == 1) {
        mat[0] = pow(mat[0], 0.5);
        *info = 0;
        return;
    }

    if (ntwin)
        eigstruc_ (mat, &n, info);

    else {
        double *pmat;
        try { pmat = new double[n*(n+1)/2]; }
        catch (...) {
            *info = -1;
            return;
        }

        int i, j, k;
        k = 0;
        for (j = 0; j < n; j++) {
            for (i = 0; i <= j; i++) {
                pmat[k] = mat[i*n+j];
                k++;
            }
        }

        dppfa_ (pmat, &n, info);

        k = 0;
        for (j = 0; j < n; j++) {
            for (i = 0; i <= j; i++) {
                mat[j*n+i] = 0;
                mat[i*n+j] = pmat[k];
                k++;
            }
        }

        delete[] pmat;
    }
}

int kincoef (double *kin2, int *fa, int *mo, int *itwin, int nind, int tnind,
             FILE *kinfp)
{
    int i, j;
    static double phi2, delta7;
    static int id1 = 0, id2;
    char rec[100];

    for (i = 0; i < nind*nind; i++)
        kin2[i] = 0;

    if (tnind && id1 <= tnind+nind) {
        i = id1 - tnind;
        j = id2 - tnind;
        kin2[(i-1)*nind+j-1] = phi2;
        kin2[(j-1)*nind+i-1] = kin2[(i-1)*nind+j-1];
    }

    while (fgets(rec, sizeof(rec), kinfp)) {
        if (sscanf(rec, "%d %d %lf %lf", &id1, &id2, &phi2, &delta7) != 4)
            return 0;
        if (id1 <= tnind)
            continue;
        if (id1 > tnind+nind)
            return 1;
        i = id1 - tnind;
        j = id2 - tnind;
        kin2[(i-1)*nind+j-1] = phi2;
        kin2[(j-1)*nind+i-1] = kin2[(i-1)*nind+j-1];
    }

    return 1;
}

/*
 *
 *     Multi-trait simulation of additive genetic deviations:
 *
 *     devA  =  (chol(A)' .*. chol(G)') * T    where  T ~ N(0,1)
 *
 *
 *     amat is the upper-triangular Cholesky of two times the kinship
 *     coefficient matrix, i.e.
 *
 *          amat(i,j) = CholA(i,j)     i <= j
 *          amat(i,j) = 0              i >  j
 *
 *     gmat is the upper-triangular Cholesky of the additive genetic
 *     covariance matrix, i.e.
 *
 *          gmat(i,j) = CholG(i,j)     i <= j
 *          gmat(i,j) = 0              i >  j
 *
 */

void simva (double *amat, double *gmat, double **gdev, int nind, int ntrt)
{
    int i, j, it, jt;

    double *t;
    try { t = new double[ntrt]; }
    catch (...) { throw Safe_Error_Return("Out of memory"); }

    double **z;
    try { z = new double*[ntrt]; }
    catch (...) { throw Safe_Error_Return("Out of memory"); }
    for (i = 0; i < ntrt; i++) {
        try { z[i] = new double[nind]; }
        catch (...) { throw Safe_Error_Return("Out of memory"); }
    }

    for (i = 0; i < nind; i++) {
        for (it = 0; it < ntrt; it++)
            t[it] = gasdev();
        for (it = 0; it < ntrt; it++) {
            z[it][i] = 0;
            for (jt = 0; jt <= it; jt++)
                z[it][i] += gmat[jt*ntrt+it]*t[jt];
        }
    }

    for (i = 0; i < nind; i++) {
        for (it = 0; it < ntrt; it++) {
            gdev[it][i] = 0;
            for (j = 0; j <= i; j++)
                gdev[it][i] += amat[j*nind+i]*z[it][j];
        }
    }

    for (i = 0; i < ntrt; i++)
        delete[] z[i];
    delete[] z;
    delete[] t;
}

/*
 *
 *     Multi-trait simulation of environmental deviations:
 *
 *     devE  =  (I .*. chol(E)') * T    where  T ~ N(0,1)
 *
 *
 *     emat is the upper-triangular Cholesky of the environmental
 *     covariance matrix, i.e.
 *
 *          emat(i,j) = CholE(i,j)     i <= j
 *          emat(i,j) = 0              i >  j
 *
 */

void simve (double *emat, double **edev, int nind, int ntrt)
{
    int i, j, it, jt;

    double *t;
    try { t = new double[ntrt]; }
    catch (...) { throw Safe_Error_Return("Out of memory"); }

    for (i = 0; i < nind; i++) {
        for (it = 0; it < ntrt; it++)
            t[it] = gasdev();
        for (it = 0; it < ntrt; it++) {
            edev[it][i] = 0;
            for (jt = 0; jt <= it; jt++)
                edev[it][i] += emat[jt*ntrt+it]*t[jt];
        }
    }

    delete[] t;
}

void dropgenel (int *fa, int *mo, int *twin, int nind, double *gfreq,
                double *mfreq, double xtheta, double ytheta, int *patgene,
                int *matgene, int nall, int nmall)
{
    double *hfreq, *sumfreq;
    try {
        hfreq = new double[nall*nmall];
        sumfreq = new double[nall*nmall + 1];
    }
    catch (...) {
        throw Safe_Error_Return("Out of memory");
    }

    int i, j;
    for (i = 0; i < nall; i++)
        for (j = 0; j < nmall; j++)
            hfreq[i*nmall+j] = gfreq[i]*mfreq[j];

    sumfreq[0] = 0;
    for (i = 1; i < nall*nmall + 1; i++)
        sumfreq[i] = sumfreq[i-1] + hfreq[i-1];

    double sumytheta[3], sumxtheta[3];
    sumytheta[0] = .5*(1 - ytheta);
    sumytheta[1] = sumytheta[0] + .5*ytheta;
    sumytheta[2] = sumytheta[1] + .5*ytheta;
    sumxtheta[0] = .5*(1 - xtheta);
    sumxtheta[1] = sumxtheta[0] + .5*xtheta;
    sumxtheta[2] = sumxtheta[1] + .5*xtheta;

    int ifa, imo;
    double z;
    for (i = 0; i < nind; i++) {
        if (twin[i] != i + 1) {
            patgene[2*i] = patgene[2*(twin[i]-1)];
            patgene[2*i+1] = patgene[2*(twin[i]-1)+1];
            matgene[2*i] = matgene[2*(twin[i]-1)];
            matgene[2*i+1] = matgene[2*(twin[i]-1)+1];
            continue;
        }

        ifa = fa[i];
        imo = mo[i];

        if (!ifa) {
            int haplo = 0;
            z = drand48();
            for (j = 0; j < nall*nmall; j++) {
                if (z > sumfreq[j] && z <= sumfreq[j+1]) {
                    haplo = j;
                    break;
                }
            }
            patgene[2*i] = haplo/nmall;
            patgene[2*i+1] = i*2*nmall + haplo%nmall;

            z = drand48();
            for (j = 0; j < nall*nmall; j++) {
                if (z > sumfreq[j] && z <= sumfreq[j+1]) {
                    haplo = j;
                    break;
                }
            }
            matgene[2*i] = haplo/nmall;
            matgene[2*i+1] = (i*2 + 1)*nmall + haplo%nmall;
        }

        else {
            ifa--;
            imo--;

            int fhet = 0, mhet = 0;
            for (j = 0; j < 2; j++) {
                if (patgene[2*ifa+j] != matgene[2*ifa+j])
                    fhet++;
                if (patgene[2*imo+j] != matgene[2*imo+j])
                    mhet++;
            }

            switch (fhet) {
            case 0:
                patgene[2*i] = patgene[2*ifa];
                patgene[2*i+1] = patgene[2*ifa+1];
                break;
            case 1:
                z = drand48();
                if (z <= 0.5) {
                    patgene[2*i] = patgene[2*ifa];
                    patgene[2*i+1] = patgene[2*ifa+1];
                }
                else {
                    patgene[2*i] = matgene[2*ifa];
                    patgene[2*i+1] = matgene[2*ifa+1];
                }
                break;
            default:
                z = drand48();
                if (z <= sumytheta[0]) {
                    patgene[2*i] = patgene[2*ifa];
                    patgene[2*i+1] = patgene[2*ifa+1];
                }
                else if (z <= sumytheta[1]) {
                    patgene[2*i] = patgene[2*ifa];
                    patgene[2*i+1] = matgene[2*ifa+1];
                }
                else if (z <= sumytheta[2]) {
                    patgene[2*i] = matgene[2*ifa];
                    patgene[2*i+1] = patgene[2*ifa+1];
                }
                else {
                    patgene[2*i] = matgene[2*ifa];
                    patgene[2*i+1] = matgene[2*ifa+1];
                }
            }

            switch (mhet) {
            case 0:
                matgene[2*i] = patgene[2*imo];
                matgene[2*i+1] = patgene[2*imo+1];
                break;
            case 1:
                z = drand48();
                if (z <= 0.5) {
                    matgene[2*i] = patgene[2*imo];
                    matgene[2*i+1] = patgene[2*imo+1];
                }
                else {
                    matgene[2*i] = matgene[2*imo];
                    matgene[2*i+1] = matgene[2*imo+1];
                }
                break;
            default:
                z = drand48();
                if (z <= sumxtheta[0]) {
                    matgene[2*i] = patgene[2*imo];
                    matgene[2*i+1] = patgene[2*imo+1];
                }
                else if (z <= sumxtheta[1]) {
                    matgene[2*i] = patgene[2*imo];
                    matgene[2*i+1] = matgene[2*imo+1];
                }
                else if (z <= sumxtheta[2]) {
                    matgene[2*i] = matgene[2*imo];
                    matgene[2*i+1] = patgene[2*imo+1];
                }
                else {
                    matgene[2*i] = matgene[2*imo];
                    matgene[2*i+1] = matgene[2*imo+1];
                }
            }
        }
    }

    delete[] sumfreq;
    delete[] hfreq;
}

void dropgene (int *fa, int *mo, int *twin, int nind, double *freq,
               int *patgene, int *matgene, int nall)
{
    double *sumfreq;
    try {
        sumfreq = new double[nall+1];
    }
    catch (...) {
        throw Safe_Error_Return("Out of memory");
    }

    int i, j;
    sumfreq[0] = 0;
    for (i = 0; i < nall; i++)
        sumfreq[i+1] = sumfreq[i] + freq[i];

    double z;
    for (i = 0; i < nind; i++) {
        if (twin[i] != i + 1) {
            patgene[2*i] = patgene[2*(twin[i]-1)];
            matgene[2*i] = matgene[2*(twin[i]-1)];
            continue;
        }

        if (!fa[i]) {
            patgene[2*i] = 0;
            z = drand48();
            for (j = 0; j < nall; j++) {
                if (z > sumfreq[j] && z <= sumfreq[j+1]) {
                    patgene[2*i] = j;
                    break;
                }
            }

            matgene[2*i] = 0;
            z = drand48();
            for (j = 0; j < nall; j++) {
                if (z > sumfreq[j] && z <= sumfreq[j+1]) {
                    matgene[2*i] = j;
                    break;
                }
            }
        }

        else {
            z = drand48();
            if (z <= 0.5)
                patgene[2*i] = patgene[2*(fa[i]-1)];
            else
                patgene[2*i] = matgene[2*(fa[i]-1)];

            z = drand48();
            if (z <= 0.5)
                matgene[2*i] = patgene[2*(mo[i]-1)];
            else
                matgene[2*i] = matgene[2*(mo[i]-1)];
        }
    }

    delete[] sumfreq;
}

double gasdev (void)
{
    static double gset;		/* save to return at next call */
    static int iset = 0;	/* saved return value? */

    double v1, v2, r, fac, gret;

    if (!iset) {
        do {
            v1 = 2*drand48() - 1;
            v2 = 2*drand48() - 1;
            r = v1*v1 + v2*v2;
        } while (r >= 1);
        fac = sqrt(-2*log(r)/r);
        gset = v1*fac;
        gret = v2*fac;
        iset = 1;
    }

    else {
        gret = gset;
        iset = 0;
    }

    return gret;
}

void SimPars::get_pars (FILE *fp)
{
    char rec[1024];
    int i, j, k;
    double sum;

    if (!fgets(rec, sizeof(rec), fp))
        throw Safe_Error_Return(
                "Error reading parameter file, premature end-of-file");

    char *p = strtok(rec, " \t");
    if (sscanf(p, "%d", &nloc) != 1)
        throw Safe_Error_Return("Error reading parameter file, line 1");
    if (nloc <= 0)
        throw Safe_Error_Return("Invalid number of QTLs");
    if (ntrt > MXLOC) {
        char buf[256];
        sprintf(buf, "Maximum number of QTLs is %d", MXLOC);
        throw Safe_Error_Return(buf);
    }

    ngen = 1;
    for (i = 0; i < nloc; i++) {
        try { loc[i] = new SimLocPars; }
        catch (...) { throw Safe_Error_Return("Out of memory"); }

        p = strtok(NULL, " \t");
        if (!p || sscanf(p, "%d", &loc[i]->nall) != 1)
            throw Safe_Error_Return("Error reading parameter file, line 1");
        if (loc[i]->nall <= 0)
            throw Safe_Error_Return("Invalid number of QTL alleles");
        ngen *= loc[i]->nall;
    }
    ngen = ngen*(ngen + 1)/2;

    p = strtok(NULL, " \t");
    if (!p || sscanf(p, "%d", &nmrk) != 1)
        throw Safe_Error_Return("Error reading parameter file, line 1");
    if (nmrk < 0)
        throw Safe_Error_Return("Invalid number of markers");
    if (nmrk != 0) {
        if (nmrk != 1 || nloc != 1)
            throw Safe_Error_Return(
                    "Only 1 QTL and 1 linked marker can be simulated");
        p = strtok(NULL, " \t");
        if (!p || sscanf(p, "%d", &nmall) != 1)
            throw Safe_Error_Return("Error reading parameter file, line 1");
        if (nmall <= 0)
            throw Safe_Error_Return("Invalid number of marker alleles");
    }

    p = strtok(NULL, " \t");
    if (!p || sscanf(p, "%d", &ntrt) != 1)
        throw Safe_Error_Return("Error reading parameter file, line 1");
    if (ntrt <= 0)
        throw Safe_Error_Return("Invalid number of traits");
    if (ntrt > MXTRT) {
        char buf[256];
        sprintf(buf, "Maximum number of traits is %d", MXTRT);
        throw Safe_Error_Return(buf);
    }

    p = strtok(NULL, " \t");
    if (!p || sscanf(p, "%d", &ncov) != 1)
        throw Safe_Error_Return("Error reading parameter file, line 1");
    if (ncov < 0)
        throw Safe_Error_Return("Invalid number of covariates");
    if (ncov > MXCOV - 3) {
        char buf[256];
        sprintf(buf, "Maximum number of covariates is %d", MXCOV - 3);
        throw Safe_Error_Return(buf);
    }

    for (i = 0; i < ncov - 2; i++) {
        if (!fgets(rec, sizeof(rec), fp))
            throw Safe_Error_Return(
                    "Error reading parameter file, premature end-of-file");
    }

    if (nloc == 1) {
        try { loc[0]->freq = new double[loc[0]->nall]; }
        catch (...) { throw Safe_Error_Return("Out of memory"); }

        if (!fgets(rec, sizeof(rec), fp))
            throw Safe_Error_Return(
                    "Error reading parameter file, premature end-of-file");
        p = strtok(rec, " \t");
        if (sscanf(p, "%lf", &loc[0]->freq[0]) != 1) {
            char buf[256];
            sprintf(buf, "Error reading parameter file, line %d", ncov);
            throw Safe_Error_Return(buf);
        }
        if (loc[0]->freq[0] <= 0 || loc[0]->freq[0] > 1)
            throw Safe_Error_Return("Invalid QTL allele frequency");

        i = 0;
        sum = loc[0]->freq[0];
        while (++i < loc[0]->nall - 1) {
            p = strtok(NULL, " \t");
            if (!p || sscanf(p, "%lf", &loc[0]->freq[i]) != 1) {
                char buf[256];
                sprintf(buf, "Error reading parameter file, line %d", ncov);
                throw Safe_Error_Return(buf);
            }
            if (loc[0]->freq[i] <= 0 || loc[0]->freq[i] > 1)
                throw Safe_Error_Return("Invalid QTL allele frequency");
            sum += loc[0]->freq[i];
        }
        loc[0]->freq[loc[0]->nall-1] = 1 - sum;
    }

    if (nmrk) {
        try { mfreq = new double[nmall]; }
        catch (...) { throw Safe_Error_Return("Out of memory"); }

        if (!fgets(rec, sizeof(rec), fp))
            throw Safe_Error_Return(
                    "Error reading parameter file, premature end-of-file");

        p = strtok(rec, " \t");
        if (sscanf(p, "%lf", &mfreq[0]) != 1) {
            char buf[256];
            sprintf(buf, "Error reading parameter file, line %d", ncov + 1);
            throw Safe_Error_Return(buf);
        }
        if (mfreq[0] <= 0 || mfreq[0] > 1)
            throw Safe_Error_Return("Invalid marker allele frequency");

        i = 0;
        sum = mfreq[0];
        while (++i < nmall - 1) {
            p = strtok(NULL, " \t");
            if (!p || sscanf(p, "%lf", &mfreq[i]) != 1) {
                char buf[256];
                sprintf(buf, "Error reading parameter file, line %d", ncov + 1);
                throw Safe_Error_Return(buf);
            }
            if (mfreq[i] <= 0 || mfreq[i] > 1)
                throw Safe_Error_Return("Invalid marker allele frequency");
            sum += mfreq[i];
        }
        mfreq[nmall-1] = 1 - sum;

        if (!fgets(rec, sizeof(rec), fp) || sscanf(rec, "%lf", &theta) != 1)
        {
            char buf[256];
            sprintf(buf, "Error reading parameter file, line %d", ncov + 2);
            throw Safe_Error_Return(buf);
        }
        if (theta < 0 || theta > .5)
            throw Safe_Error_Return("Invalid recombination fraction");
    }

    for (i = 0; i < ntrt; i++) {
        try { trt[i] = new SimTrtPars(ncov+3); }
        catch (...) { throw Safe_Error_Return("Out of memory"); }

        try { trt[i]->mean = new double[ngen]; }
        catch (...) { throw Safe_Error_Return("Out of memory"); }

        if (fscanf(fp, "%lf", &trt[i]->mean[0]) != 1) {
            char buf[256];
            sprintf(buf, "Error reading genotypic mean 1 for trait %d", i + 1);
            throw Safe_Error_Return(buf);
        }
        j = 0;
        while (++j < ngen) {
            if (fscanf(fp, "%lf", &trt[i]->mean[j]) != 1) {
                char buf[256];
                sprintf(buf, "Error reading genotypic mean %d for trait %d",
                        j + 1, i + 1);
                throw Safe_Error_Return(buf);
            }
        }

        if (fscanf(fp, "%lf", &trt[i]->sdev) != 1) {
            char buf[256];
            sprintf(buf,
                    "Error reading phenotypic standard deviation for trait %d",
                    i + 1);
            throw Safe_Error_Return(buf);
        }
        if (trt[i]->sdev <= 0)
            throw Safe_Error_Return("Invalid phenotypic standard deviation");

        if (fscanf(fp, "%lf", &trt[i]->h2r) != 1) {
            char buf[256];
            sprintf(buf, "Error reading residual heritability for trait %d",
                    i + 1);
            throw Safe_Error_Return(buf);
        }
        if (trt[i]->h2r < 0 || trt[i]->h2r > 1)
            throw Safe_Error_Return("Invalid residual heritability");

        try { trt[i]->cmean = new double[ncov + 3]; }
        catch (...) { throw Safe_Error_Return("Out of memory"); }

        for (j = 0; j < ncov + 3; j++) {
            if (fscanf(fp, "%lf", &trt[i]->cmean[j]) != 1) {
                char buf[256];
                sprintf(buf,
                        "Error reading mean value for covariate %d, trait %d",
                        j + 1, i + 1);
                throw Safe_Error_Return(buf);
            }

            try { trt[i]->beta[j] = new double[ngen]; }
            catch (...) { throw Safe_Error_Return("Out of memory"); }

            for (k = 0; k < ngen; k++) {
                if (fscanf(fp, "%lf", &trt[i]->beta[j][k]) != 1) {
                    char buf[256];
                    sprintf(buf,
                            "Error reading beta %d for covariate %d, trait %d",
                            k + 1, j + 1, i + 1);
                    throw Safe_Error_Return(buf);
                }
            }
        }

        if (i) {
            try { trt[i]->rg = new double[i]; }
            catch (...) { throw Safe_Error_Return("Out of memory"); }

            try { trt[i]->re = new double[i]; }
            catch (...) { throw Safe_Error_Return("Out of memory"); }

            j = 0;
            while (j < i) {
                if (fscanf(fp, "%lf %lf",
                           &trt[i]->rg[j], &trt[i]->re[j]) != 2)
                {
                    char buf[256];
                    sprintf(buf,
                            "Error reading correlations for traits %d and %d",
                            i + 1, j + 1);
                    throw Safe_Error_Return(buf);
                }
                if (trt[i]->rg[j] < -1 || trt[i]->rg[j] > 1)
                    throw Safe_Error_Return("Invalid genetic correlation");
                if (trt[i]->re[j] < -1 || trt[i]->re[j] > 1)
                    throw Safe_Error_Return(
                            "Invalid environmental correlation");
                j++;
            }
        }
    }
}

void SimPars::prt_pars (FILE *fp)
{
    int i, j, k;

    fprintf(fp, "#qtl = %d  #all = %d", nloc, loc[0]->nall);
    for (i = 1; i < nloc; i++)
        fprintf(fp, ", %d", loc[i]->nall);
    fprintf(fp, "\n");

    if (nmrk)
        fprintf(fp, "#mrk = %d  #all = %d\n", nmrk, nmall);

    fprintf(fp, "#gen = %d\n", ngen);
    fprintf(fp, "#trt = %d\n", ntrt);
    fprintf(fp, "#cov = %d (including sex and age)\n", ncov);

    if (nloc == 1) {
        fprintf(fp, "qtl:\n");
        fprintf(fp, "  freq = %s", fp2str(loc[0]->freq[0]));
        for (j = 1; j < loc[0]->nall; j++)
            fprintf(fp, " %s", fp2str(loc[0]->freq[j]));
        fprintf(fp, "\n");
    }

    if (nmrk) {
        fprintf(fp, "marker:\n");
        fprintf(fp, "  freq = %s", fp2str(mfreq[0]));
        for (j = 1; j < nmall; j++)
            fprintf(fp, " %s", fp2str(mfreq[j]));
        fprintf(fp, "\n");
        fprintf(fp, "theta = %s\n", fp2str(theta));
    }

    for (i = 0; i < ntrt; i++) {
        fprintf(fp, "trait %d:\n", i+1);
        fprintf(fp, "  mean = %s", fp2str(trt[i]->mean[0]));
        for (j = 1; j < ngen; j++)
            fprintf(fp, " %s", fp2str(trt[i]->mean[j]));
        fprintf(fp, "\n");
        fprintf(fp, "  sdev = %s\n", fp2str(trt[i]->sdev));
        fprintf(fp, "  h2r  = %s\n", fp2str(trt[i]->h2r));

        for (j = 0; j < ncov + 3; j++) {
            switch (j) {
            case 0:
                fprintf(fp, "  covariate sex:\n");
                break;
            case 1:
                fprintf(fp, "  covariate male age:\n");
                break;
            case 2:
                fprintf(fp, "  covariate female age:\n");
                break;
            case 3:
                fprintf(fp, "  covariate male age^2:\n");
                break;
            case 4:
                fprintf(fp, "  covariate female age^2:\n");
                break;
            default:
                fprintf(fp, "  covariate %d:\n", j - 4);
            }
            if (j)
                fprintf(fp, "    mean = %s\n", fp2str(trt[i]->cmean[j]));
            fprintf(fp, "    beta = %s", fp2str(trt[i]->beta[j][0]));
            for (k = 1; k < ngen; k++)
                fprintf(fp, " %s", fp2str(trt[i]->beta[j][k]));
            fprintf(fp, "\n");
        }

        if (i) {
            fprintf(fp, "  rhog = %s", fp2str(trt[i]->rg[0]));
            for (j = 1; j < i; j++)
                fprintf(fp, " %s", fp2str(trt[i]->rg[j]));
            fprintf(fp, "\n");
            fprintf(fp, "  rhoe = %s", fp2str(trt[i]->re[0]));
            for (j = 1; j < i; j++)
                fprintf(fp, " %s", fp2str(trt[i]->re[j]));
            fprintf(fp, "\n");
        }
    }

    fflush(stdout);
}

char *fp2str (double num)
{
    static char str[128];
    char *p;
    sprintf(str, "%f", num);
    p = str + strlen(str);
    while (*--p == '0')
        *p = 0;
    return str;
}
