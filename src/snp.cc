/*
 * snp.cc implements the snp command
 * Written by Thomas Dyer beginning on September 5, 2003
 * Copyright (c) 2003 Southwest Foundation for Biomedical Research
 */

#include <math.h>
#include "solar.h"
// tcl.h from solar.h
#include "safelib.h"

#define MISSVAL -1

struct SnpHap {
    double freq;
    int nsnp;
    char **all;
    SnpHap() {}
    ~SnpHap();
};

SnpHap::~SnpHap ()
{
    if (all) {
        for (int i = 0; i < nsnp; i++) {
            if (all[i]) delete[] all[i];
        }
        delete[] all;
    }
}

struct HapTree {
    int hap;
    char *lall;
    struct HapTree *left;
    char *rall;
    struct HapTree *right;
    HapTree() {lall = 0; left = 0; rall = 0; right = 0;}
    ~HapTree();
};

HapTree::~HapTree ()
{
    if (lall) {
        delete[] lall;
        delete left;
    }
    if (rall) {
        delete[] rall;
        delete right;
    }
}

int make_genocovars (bool, Tcl_Interp *);
int make_haplocovars (Tcl_Interp *);
int compute_ld (int, Tcl_Interp *);
int mult_corr_liji (int, double, Tcl_Interp *);
void prt_haplo_tree (int, struct HapTree *);

extern "C" void symeig_ (int*, double*, double*, double*, double*, int*);

extern "C" int SnpCmd (ClientData clientData, Tcl_Interp *interp,
                       int argc, char *argv[])
{
    if (argc == 2 && !StringCmp ("help", argv[1], case_ins))
    {
        return Solar_Eval (interp, "help snp");
    }

    if (argc >= 2 && !StringCmp ("covar", argv[1], case_ins))
    {
        bool impute = false;

        if (argc == 3 && !StringCmp ("-impute", argv[2], case_ins))
            impute = true;

        return make_genocovars(impute, interp);
    }

    if (argc >= 2 && !StringCmp ("hapcov", argv[1], case_ins))
    {
        return make_haplocovars(interp);
    }

    if (argc >= 2 && !StringCmp ("ld", argv[1], case_ins))
    {
        int ld_window = 0;

        if (argc == 4 && !StringCmp ("-window", argv[2], case_ins))
        {
            ld_window = atoi(argv[3]);
        }

        return compute_ld(ld_window, interp);
    }

    if (argc >= 2 && !StringCmp ("liji", argv[1], case_ins))
    {
        if (argc == 4)
        {
            int nsnp = atoi(argv[2]);
            double alpha = atof(argv[3]);
            return mult_corr_liji(nsnp, alpha, interp);
        }
    }

    RESULT_LIT ("Invalid snp command");
    return TCL_ERROR;
}

int
make_genocovars (bool impute, Tcl_Interp *interp)
{
    int i, j, nhap;
    char **record;
    const char *errmsg = 0;
    SolarFile *Tfile;
    SnpHap *haplo;

    int nsnp;

    if (impute) {
        Tfile = SolarFile::open("SNP haplotype frequencies",
                                "snp.haplofreqs", &errmsg);
        if (errmsg) {
            Solar_AppendResult2(interp,
                           "Error opening SNP haplotype frequencies file:\n",
                                errmsg, NULL);
            return TCL_ERROR;
        }

        Tfile->start_setup(&errmsg);
        Tfile->setup("freq", &errmsg);
        nsnp = 0;
        for (i = 0; i < currentMap->nloci(); i++) {
            int mrk = currentFreq->get_marker(currentMap->mrkname(i));
            if (currentFreq->nall(mrk) < 2)
                continue;
            Tfile->setup(currentMap->mrkname(i), &errmsg);
            nsnp++;
        }

        if (errmsg) {
            Solar_AppendResult2(interp,
                                "SNP haplotype frequency file error:\n",
    	                        errmsg, NULL);
            return TCL_ERROR;
        }

        Tfile->rewind(&errmsg);
        nhap = 0;
        while (1) {
            record = Tfile->get(&errmsg);
            if (errmsg && !strcmp("EOF", errmsg)) break;
            if (errmsg) {
                Solar_AppendResult2(interp,
                            "Error reading SNP haplotype frequency file:\n",
                                    errmsg, NULL);
                return TCL_ERROR;
            }
            nhap++;
        }

        try { haplo = new SnpHap[nhap]; }
        catch (...) {
            RESULT_LIT ("Out of memory");
            return TCL_ERROR;
        }
        for (i = 0; i < nhap; i++) {
            haplo[i].nsnp = nsnp;
            try {
                haplo[i].all = new char*[nsnp];
                for (j = 0; j < nsnp; j++) {
                    haplo[i].all[j] = new char[MALLNM+1];
                }
            }
            catch (...) {
                delete[] haplo;
                RESULT_LIT ("Out of memory");
                return TCL_ERROR;
            }
        }

        Tfile->rewind(&errmsg);
        for (i = 0; i < nhap; i++) {
            record = Tfile->get(&errmsg);
            if (errmsg && !strcmp("EOF", errmsg)) break;
            if (errmsg) {
                Solar_AppendResult2(interp,
                            "Error reading SNP haplotype frequency file:\n",
                                    errmsg, NULL);
                return TCL_ERROR;
            }
            for (j = 0; j < nsnp; j++)
                strcpy(haplo[i].all[j], record[j+1]);
            sscanf(record[0], "%lf", &haplo[i].freq);
        }
        delete Tfile;
    }

    Tfile = SolarFile::open("SNP haplotypes", "snp.haplotypes", &errmsg);
    if (errmsg) {
        Solar_AppendResult2(interp, "Error opening SNP haplotypes file:\n", 
                            errmsg, NULL);
        return TCL_ERROR;
    }

    char *id, *famid;
    char *currid, *currfamid;
    int count, pos;
    bool need_famid = false;
    int *widths = Tfile->widths(&count, &errmsg);

    Tfile->start_setup(&errmsg);
    if (currentPed->famid_len() && Tfile->test_name("FAMID", &errmsg)) {
        need_famid = true;
        pos = Tfile->setup("FAMID", &errmsg);
        famid = new char[widths[pos] + 1];
        currfamid = new char[widths[pos] + 1];
    }
    pos = Tfile->setup("ID", &errmsg);
    id = new char[widths[pos] + 1];
    currid = new char[widths[pos] + 1];

    int maxwid = 0;
    nsnp = 0;
    for (i = 0; i < currentMap->nloci(); i++) {
        int mrk = currentFreq->get_marker(currentMap->mrkname(i));
        if (currentFreq->nall(mrk) < 2)
            continue;
        pos = Tfile->setup(currentMap->mrkname(i), &errmsg);
        if (widths[pos] > maxwid) maxwid = widths[pos];
        nsnp++;
    }

    if (errmsg) {
        Solar_AppendResult2(interp, "SNP haplotypes file error:\n",
                            errmsg, NULL);
        return TCL_ERROR;
    }

    FILE *fp = fopen("snp.genocov", "w");
    if (!fp) {
        RESULT_LIT ("Cannot open file snp.genocov");
        return TCL_ERROR;
    }

    FILE *fp2 = fopen("snp.geno-list", "w");
    if (!fp2) {
        RESULT_LIT ("Cannot open file snp.geno-list");
        return TCL_ERROR;
    }

    if (need_famid) fprintf(fp, "FAMID,");
    fprintf(fp, "ID");
    for (i = 0; i < currentMap->nloci(); i++) {
        int mrk = currentFreq->get_marker(currentMap->mrkname(i));
        if (currentFreq->nall(mrk) < 2)
            continue;
        fprintf(fp, ",snp_%s", currentMap->mrkname(i));
        fprintf(fp2, "snp_%s\n", currentMap->mrkname(i));
    }
    fprintf(fp, "\n");
    fclose(fp2);

    int ncons;
    bool *cons;
    if (impute) {
        cons = new bool[nhap];
    }

    int snp1, which;
    bool typed;
    double sumf, *prob[2];

    prob[0] = new double[nsnp];
    prob[1] = new double[nsnp];

    char **snp = new char*[nsnp];
    for (i = 0; i < nsnp; i++)
        snp[i] = new char[maxwid + 1];
    char **rare = new char*[nsnp];
    for (i = 0; i < nsnp; i++)
        rare[i] = new char[maxwid + 1];

    nsnp = 0;
    for (i = 0; i < currentMap->nloci(); i++) {
        int mrk = currentFreq->get_marker(currentMap->mrkname(i));
        if (currentFreq->nall(mrk) < 2)
            continue;
        if (currentFreq->all(mrk,0)->freq > currentFreq->all(mrk,1)->freq)
            strcpy(rare[nsnp], currentFreq->all(mrk,1)->name);
        else
            strcpy(rare[nsnp], currentFreq->all(mrk,0)->name);
        nsnp++;
    }

    if (need_famid) {
        currfamid[0] = '\0';
        snp1 = 2;
    } else {
        currid[0] = '\0';
        snp1 = 1;
    }

    while (1) {
        record = Tfile->get(&errmsg);
        if (errmsg && !strcmp("EOF", errmsg)) break;
        if (errmsg) {
            Solar_AppendResult2(interp,
                "Error reading SNP haplotypes file:\n", errmsg, NULL);
            return TCL_ERROR;
        }

        which = 1;
        if (need_famid) {
            strcpy(famid, record[0]);
            strcpy(id, record[1]);
            if (strcmp(currfamid, famid) || strcmp(currid, id)) {
                strcpy(currfamid, famid);
                strcpy(currid, id);
                which = 0;
                fprintf(fp, "%s,%s", famid, id);
            }
        }
        else {
            strcpy(id, record[0]);
            if (strcmp(currid, id)) {
                strcpy(currid, id);
                which = 0;
                fprintf(fp, "%s", id);
            }
        }

        for (i = 0; i < nsnp; i++)
            strcpy(snp[i], record[snp1+i]);

        if (impute) {
            ncons = 0;
            sumf = 0;
            for (i = 0; i < nhap; i++) {
                cons[i] = true;
                for (j = 0; j < nsnp; j++) {
                    if (!strcmp(snp[j], ""))
                        continue;
                    if (strcmp(snp[j], haplo[i].all[j])) {
                        cons[i] = false;
                        break;
                    }
                }
                if (cons[i]) {
                    ncons++;
                    sumf += haplo[i].freq;
                }
            }

            if (!ncons) sumf = 1;
        }

        for (i = 0; i < nsnp; i++) {
            typed = true;
            prob[which][i] = 0;
            if (!strcmp(snp[i], rare[i])) {
                prob[which][i] = 1;
            }
            else if (!strcmp(snp[i], "")) {
                if (!impute) {
                    typed = false;
                }
                else {
                    for (j = 0; j < nhap; j++) {
                        if (ncons && !cons[j])
                            continue;
                        if (!strcmp(haplo[j].all[i], rare[i]))
                            prob[which][i] += haplo[j].freq/sumf;
                    }
                }
            }
            if (which) {
                if (impute || typed)
                    fprintf(fp, ",%g", prob[0][i] + prob[1][i]);
                else
                    fprintf(fp, ",");
            }
        }
        if (which)
            fprintf(fp, "\n");
    }
    fclose(fp);
    delete Tfile;

    for (i = 0; i < nsnp; i++)
        delete[] rare[i];
    delete[] rare;

    for (i = 0; i < nsnp; i++)
        delete[] snp[i];
    delete[] snp;

    delete prob[1], prob[0];

    if (impute) {
        delete cons;
        delete[] haplo;
    }

    return TCL_OK;
}

int
make_haplocovars (Tcl_Interp *interp)
{
    int i, j;
    char **record;
    const char *errmsg = 0;
    SolarFile *Tfile;

    Tfile = SolarFile::open("SNP haplotype frequencies",
                            "snp.haplofreqs", &errmsg);
    if (errmsg) {
        Solar_AppendResult2(interp,
                            "Error opening SNP haplotype frequencies file:\n",
                            errmsg, NULL);
        return TCL_ERROR;
    }

    int count, pos;
    int *widths = Tfile->widths(&count, &errmsg);

    Tfile->start_setup(&errmsg);
    int nsnp = 0;
    int maxwid = 0;
    for (i = 0; i < currentMap->nloci(); i++) {
        int mrk = currentFreq->get_marker(currentMap->mrkname(i));
        if (currentFreq->nall(mrk) < 2)
            continue;
        pos = Tfile->setup(currentMap->mrkname(i), &errmsg);
        if (widths[pos] > maxwid) maxwid = widths[pos];
        nsnp++;
    }

    if (errmsg) {
        Solar_AppendResult2(interp, "SNP haplotype frequency file error:\n",
	                    errmsg, NULL);
        return TCL_ERROR;
    }

    HapTree *p, *haplos = new HapTree;
    int nhap = 0;
    while (1) {
        record = Tfile->get(&errmsg);
        if (errmsg && !strcmp("EOF", errmsg)) break;
        if (errmsg) {
            Solar_AppendResult2(interp,
                        "Error reading SNP haplotype frequency file:\n",
                                errmsg, NULL);
            return TCL_ERROR;
        }
        p = haplos;
        try {
            for (i = 0; i < nsnp; i++) {
                if (p->lall == 0) {
                    p->lall = new char[maxwid+1];
                    strcpy(p->lall, record[i]);
                    p->left = new HapTree;
                    p = p->left;
                    p->hap = nhap;
                }
                else if (!strcmp(p->lall, record[i])) {
                    p = p->left;
                }
                else if (p->rall == 0) {
                    p->rall = new char[maxwid+1];
                    strcpy(p->rall, record[i]);
                    p->right = new HapTree;
                    p = p->right;
                    p->hap = nhap;
                }
                else if (!strcmp(p->rall, record[i])) {
                    p = p->right;
                }
                else {
//                    prt_haplo_tree(haplos);
                    RESULT_LIT ("Haplotype tree corrupted");
                    return TCL_ERROR;
                }
            }
        }
        catch (...) {
            RESULT_LIT ("Out of memory");
            return TCL_ERROR;
        }
        nhap++;
    }
    delete Tfile;

    Tfile = SolarFile::open("SNP haplotypes",
                                       "snp.haplotypes", &errmsg);
    if (errmsg) {
        Solar_AppendResult2(interp, "Error opening SNP haplotypes file:\n", 
                            errmsg, NULL);
        return TCL_ERROR;
    }

    bool need_famid = false;
    char *id, *famid;
    char *currid, *currfamid;

    widths = Tfile->widths(&count, &errmsg);
    Tfile->start_setup(&errmsg);
    if (currentPed->famid_len() && Tfile->test_name("FAMID", &errmsg)) {
        need_famid = true;
        pos = Tfile->setup("FAMID", &errmsg);
        famid = new char[widths[pos] + 1];
        currfamid = new char[widths[pos] + 1];
    }
    pos = Tfile->setup("ID", &errmsg);
    id = new char[widths[pos] + 1];
    currid = new char[widths[pos] + 1];

    for (i = 0; i < currentMap->nloci(); i++) {
        int mrk = currentFreq->get_marker(currentMap->mrkname(i));
        if (currentFreq->nall(mrk) < 2)
            continue;
        Tfile->setup(currentMap->mrkname(i), &errmsg);
    }

    if (errmsg) {
        Solar_AppendResult2(interp, "SNP haplotypes file error:\n",
                            errmsg, NULL);
        return TCL_ERROR;
    }

    FILE *fp = fopen("snp.haplocov", "w");
    if (!fp) {
        RESULT_LIT ("Cannot open file snp.haplocov");
        return TCL_ERROR;
    }

    FILE *fp2 = fopen("snp.haplo-list", "w");
    if (!fp2) {
        RESULT_LIT ("Cannot open file snp.haplo-list");
        return TCL_ERROR;
    }

    if (need_famid) fprintf(fp, "FAMID,");
    fprintf(fp, "ID");
    for (i = 1; i <= nhap; i++) {
        fprintf(fp, ",hap_%d", i);
        fprintf(fp2, "hap_%d\n", i);
    }
    fprintf(fp, "\n");
    fclose(fp2);

    int snp1, which;
    int ncopies, hap[2];
    bool missing;

    if (need_famid) {
        currfamid[0] = '\0';
        snp1 = 2;
    } else {
        currid[0] = '\0';
        snp1 = 1;
    }

    while (1) {
        record = Tfile->get(&errmsg);
        if (errmsg && !strcmp("EOF", errmsg)) break;
        if (errmsg) {
            Solar_AppendResult2(interp,
                "Error reading SNP haplotypes file:\n", errmsg, NULL);
            return TCL_ERROR;
        }

        which = 1;
        if (need_famid) {
            strcpy(famid, record[0]);
            strcpy(id, record[1]);
            if (strcmp(currfamid, famid) || strcmp(currid, id)) {
                strcpy(currfamid, famid);
                strcpy(currid, id);
                which = 0;
            }
        }
        else {
            strcpy(id, record[0]);
            if (strcmp(currid, id)) {
                strcpy(currid, id);
                which = 0;
            }
        }

        if (!which)
            missing = false;

        p = haplos;
        for (i = 0; !missing && i < nsnp; i++) {
            if (!strcmp(record[snp1+i], "")) {
                missing = true;
                break;
            }
            else if (p->lall && !strcmp(record[snp1+i], p->lall)) {
                p = p->left;
                hap[which] = p->hap;
            }
            else if (p->rall && !strcmp(record[snp1+i], p->rall)) {
                p = p->right;
                hap[which] = p->hap;
            }
            else {
                missing = true;
                break;
            }
        }
        if (which && !missing) {
            if (need_famid)
                fprintf(fp, "%s,%s", famid, id);
            else
                fprintf(fp, "%s", id);
            for (i = 0; i < nhap; i++) {
                ncopies = 0;
                if (hap[0] == i) ncopies++;
                if (hap[1] == i) ncopies++;
                fprintf(fp, ",%d", ncopies);
            }
            fprintf(fp, "\n");
        }
    }
    fclose(fp);
    delete Tfile;

    delete haplos;

    return TCL_OK;
}

int
compute_ld (int ld_window, Tcl_Interp *interp)
{
    int i, j, k;
    int nsnp, nind;
    char snpname[MMRKNM+5];
    char *id, *famid;
    char **record;
    const char *errmsg = 0;

    SolarFile *Tfile = SolarFile::open("SNP genotypes covariates",
                                       "snp.genocov", &errmsg);
    if (errmsg) {
        Solar_AppendResult2(interp, "Error opening file snp.genocov:\n",
                            errmsg, NULL);
        return TCL_ERROR;
    }

    Tfile->start_setup(&errmsg);
    nsnp = 0;
    int count;
    const char **names;
    names = Tfile->names(&count, &errmsg);
    for (i = 0; i < count; i++) {
        if (strncmp("snp_", names[i], 4))
            continue;
        Tfile->setup(names[i], &errmsg);
        nsnp++;
    }

    nind = 0;
    while (1) {
        record = Tfile->get(&errmsg);
        if (errmsg && !strcmp("EOF", errmsg)) break;
        if (errmsg) {
            Solar_AppendResult2(interp, "Error reading file snp.genocov:\n",
                                errmsg, NULL);
            return TCL_ERROR;
        }
        nind++;
    }

    double **snp = new double*[nind];
    for (i = 0; i < nind; i++) {
        snp[i] = new double[nsnp];
    }

    Tfile->rewind(&errmsg);
    for (i = 0; i < nind; i++) {
        record = Tfile->get(&errmsg);
        if (errmsg && !strcmp("EOF", errmsg)) break;
        if (errmsg) {
            Solar_AppendResult2(interp, "Error reading file snp.genocov:\n",
                                errmsg, NULL);
            return TCL_ERROR;
        }

        for (j = 0; j < nsnp; j++) {
            if (!strcmp(record[j], "")) {
                snp[i][j] = MISSVAL;
                continue;
            }
            sscanf(record[j], "%lf", &snp[i][j]);
        }
    }

    int n;
    int* locn = (int*) Calloc (nsnp, sizeof (int));
    double ld, xy, sx, sy, sxx, syy;
    FILE *fp;

    nsnp = 0;
    for (i = 0; i < count; i++) {
        if (strncmp("snp_", names[i], 4))
            continue;
        sscanf(names[i], "snp_%s", snpname);
        if (currentMap) {
            int mrk = currentMap->get_marker(snpname);
            if (mrk < 0) {
                char buf[1024];
                sprintf(buf, "SNP %s not found in map.", snpname);
                RESULT_BUF (buf);
                return TCL_ERROR;
            }
            locn[nsnp] = int(currentMap->mrklocn(mrk));
        }
        else {
            locn[nsnp] = 0;
        }
        nsnp++;
    }

    fp = fopen("snp.ld.dat", "w");
    if (!fp) {
        RESULT_LIT ("Cannot open file snp.ld.dat");
        return TCL_ERROR;
    }

    fprintf(fp, "M1 M2 DISEQ\n");
    for (i = 0; i < nsnp; i++) {
        for (j = i; j < nsnp; j++) {
            if (ld_window && locn[j] - locn[i] > ld_window)
                continue;
            n = 0;
            xy = 0;
            sx = 0; sxx = 0;
            sy = 0; syy = 0;
            for (k = 0; k < nind; k++) {
                if (snp[k][i] != MISSVAL && snp[k][j] != MISSVAL) {
                    n++;
                    xy += snp[k][i]*snp[k][j];
                    sx += snp[k][i];
                    sy += snp[k][j];
                    sxx += snp[k][i]*snp[k][i];
                    syy += snp[k][j]*snp[k][j];
                }
            }
            if (n == 0) {
                ld = 0;
            }
            else if (sxx/n - pow(sx/n, 2.) == 0 || syy/n - pow(sy/n, 2.) == 0) {
                ld = 0;
            }
            else {
                ld = (xy/n  - sx/n*sy/n) / 
                     sqrt((sxx/n - pow(sx/n, 2.))*(syy/n - pow(sy/n, 2.)));
            }
            fprintf(fp, "%d %d %f\n", i+1, j+1, ld);
        }
    }
    fclose(fp);
    delete Tfile;

    return TCL_OK;
}

int
mult_corr_liji (int nsnp, double alpha, Tcl_Interp *interp)
{
    int i, j;
    double ld;
    char buf[1024];
    FILE *fp;

    double* r = (double*) Calloc (nsnp*nsnp, sizeof (double));

    fp = fopen("snp.ld.dat", "r");
    if (!fp) {
        RESULT_LIT ("Cannot open snp.ld.dat");
        return TCL_ERROR;
    }

    fgets(buf, sizeof(buf), fp);
    while (fgets(buf, sizeof(buf), fp)) {
        sscanf(buf, "%d %d %lf", &i, &j, &ld);
        i--; j--;
        r[i*nsnp+j] = r[j*nsnp+i] = ld;
    }
    fclose(fp);

    int info = 0;
    double* z = (double*) Calloc (nsnp*nsnp, sizeof (double));
    double* e =  (double*) Calloc (nsnp, sizeof (double));
    double* e2 =  (double*) Calloc (nsnp, sizeof (double));

    symeig_(&nsnp, r, e, e2, z, &info);
    if (info) {
        char buf[1024];
        sprintf(buf, "Eigenvalue computation failed (info = %d)", info);
        RESULT_BUF (buf);
        return TCL_ERROR;
    }

    double neff = 0;
    for (i = 0; i < nsnp; i++) {
        if (e[i] < 0) e[i] = 0;
        if (e[i] >= 1) neff += 1;
        neff += e[i] - floor(e[i]);
    }

    double newalpha = 1 - exp(log(1-alpha)/neff);

    printf("Experiment-Wide Alpha = %g\n", alpha);
    printf("Number of SNPs = %d\n", nsnp);
    printf("Effective Number of SNPs = %g\n", neff);
    printf("Effective Number/Total Number = %g\n", neff/nsnp);
    printf("Required Alpha = %g\n", newalpha);
    printf("Target -log(p) = %g\n", -log(newalpha)/log(10.));
    printf("Effective number of SNPs and required alpha written to the file snp.effnum.\n");

    fp = fopen("snp.effnum", "w");
    if (!fp) {
        RESULT_LIT ("Cannot open file snp.effnum");
        return TCL_ERROR;
    }
    fprintf(fp, "method = Li&Ji\n");
    fprintf(fp, "nsnp = %d\n", nsnp);
    fprintf(fp, "effnum = %g\n", neff);
    fprintf(fp, "alpha = %g\n", newalpha);
    fclose(fp);

    return TCL_OK;
}

void
prt_haplo_tree (int level, struct HapTree *tree)
{
    level++;
    if (!tree->lall && !tree->rall) {
        printf("%d: hap=%d\n", level, tree->hap);fflush(stdout);
        return;
    }
    if (tree->lall) {
        printf("%d: lall=%s\n", level, tree->lall);fflush(stdout);
        prt_haplo_tree(level, tree->left);
    }
    if (tree->rall) {
        printf("%d: rall=%s\n", level, tree->rall);fflush(stdout);
        prt_haplo_tree(level, tree->right);
    }
}
