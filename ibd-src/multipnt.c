#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <math.h>
#ifdef __sun
#ifndef __GNUC__
#include <sunmath.h>
#endif
#endif
#include "pipeback.h"

/* make sure that Int4 is a 4-byte integer! */
typedef int Int4;

#define NCOEFF 14
#define MCLASS 1000
#define MLOCI  5000

struct _class {
    int ikin;
    int cnstr;
    char relation[1000];
    double phi2;
    double evar;
    double c[NCOEFF];
} *classes[MCLASS];

void *allocMem (size_t);

#define ri(a,b)       (ri[(a)*nloci + (b)])
#define difmat(a,b)   (difmat[(a)*nloci + (b)])

main (int argc, char **argv)
{
    int i, j, line, id1, id2, class;
    Int4 nloci;
    double qtloc, mrklocn[MLOCI];
    double window;
    int inwindow[MLOCI];
    char map_func = 'k';
    double phi2, *ibd, *tmu, *tsd;
    double *mu[MCLASS], *sd[MCLASS];
    double theta, *beta, *cov, dis, *difmat;
    double betasum, ibdadj, corr();
    int *mrk;
    int done[MCLASS];
    double *s_beta[MCLASS];
    Int4 nloc, info, one = 1;
    Int4 *ipvt;
    double *ri, *work, det[2];
    char rec[100000];
    FILE *locfp, *ibdfp, *meanfp, *ttyfp, *mibdfp;
    const char *pbarg[4];

    if (argc != 8) {
        printf(
"usage: multipnt locus_file ibd_file mean_file qtloc cqtloc window showstat\n");
        exit(1);
    }

    locfp = fopen(argv[1], "r");
    if (!locfp) {
        printf("Locus file not found.\n");
        exit(1);
    }

    ibdfp = fopen(argv[2], "r");
    if (!ibdfp) {
        printf(
"Merged IBD file not found. Use the command \"mibd merge\" to create it.\n");
        exit(1);
    }
    fclose(ibdfp);

    pbarg[0] = "gunzip";
    pbarg[1] = "-c";
    pbarg[2] = argv[2];
    pbarg[3] = 0;

    ibdfp = pipeback_shell_open("gunzip", pbarg);
    if (!ibdfp) {
        printf("Cannot uncompress merged IBD file.\n");
        exit(1);
    }

    meanfp = fopen(argv[3], "r");
    if (!meanfp) {
        printf(
"Mean IBD file not found. Use the command \"mibd means\" to create it.\n");
        exit(1);
    }

    mibdfp = fopen("mibd.out", "w");
    if (!mibdfp) {
        printf("Cannot open mibd.out.\n");
        exit(1);
    }

    if (sscanf(argv[4], "%lf", &qtloc) != 1) {
        printf("Invalid QTL location [%s]\n", argv[4]);
        exit(1);
    }

    if (sscanf(argv[6], "%lf", &window) != 1) {
        printf("Invalid window size [%s]\n", argv[6]);
        exit(1);
    }

    if (argv[7][0] == 'y') {
        ttyfp = fopen("/dev/tty", "w");
        if (!ttyfp) {
            printf("Cannot open /dev/tty\n");
            exit(1);
        }
        fprintf(ttyfp, "location %s ", argv[5]);
        fclose(ttyfp);
    }

    if (!(nloci = getLocInfo(locfp, mrklocn, &map_func))) {
        exit(1);
    }

    for (i = 0; i < nloci; i++) {
        inwindow[i] = 0;
        if (fabs(mrklocn[i] - qtloc) <= .5*window)
            inwindow[i] = 1;
    }

    ibd = (double *) allocMem(nloci*sizeof(double));
    tmu = (double *) allocMem(nloci*sizeof(double));
    tsd = (double *) allocMem(nloci*sizeof(double));
    beta = (double *) allocMem(nloci*sizeof(double));
    cov = (double *) allocMem(nloci*sizeof(double));
    mrk = (int *) allocMem(nloci*sizeof(int));
    ipvt = (Int4 *) allocMem(nloci*sizeof(Int4));
    work = (double *) allocMem(nloci*sizeof(double));

    difmat = (double *) allocMem(nloci*nloci*sizeof(double));
    ri = (double *) allocMem(nloci*nloci*sizeof(double));

    if (!getClasses()) {
        exit(1);
    }

    for (i = 0; i < MCLASS; i++) {
        done[i] = 0;
        mu[i] = (double *) allocMem(nloci*sizeof(double));
        sd[i] = (double *) allocMem(nloci*sizeof(double));
        s_beta[i] = (double *) allocMem(nloci*sizeof(double));
        for (j = 0; j < nloci; j++) {
            mu[i][j] = 0;
            sd[i][j] = 0;
        }
    }

    line = 0;
    while (fgets(rec, sizeof(rec), meanfp)) {
        line++;
        if (!readMeanRec(rec, line, nloci, &class, tmu, tsd))
            exit(1);
        for (j = 0; j < nloci; j++) {
            mu[class][j] = tmu[j];
            sd[class][j] = tsd[j];
        }
    }
    fclose(meanfp);

    for (i = 0; i < nloci; i++) {
        for (j = i; j < nloci; j++) {
            difmat(i,j) = fabs(mrklocn[i] - mrklocn[j]) / 100.;
            difmat(j,i) = difmat(i,j);
        }
    }

    line = 0;
    while (fgets(rec, sizeof(rec), ibdfp)) {
        line++;
        if (!readIbdRec(rec, line, nloci, &id1, &id2, &class, &phi2, ibd))
            exit(1);

        nloc = 0;
        j = 0;
        for (i = 0; i < nloci; i++) {
            if (inwindow[i]) {
                if (ibd[i] != -1) {
                    mrk[j] = i;
                    nloc++;
                    j++;
                }
            }
        }

        if (class <= 2 || class == 300) {
            ibdadj = phi2;
        }

        else if (nloc == 0) {
            ibdadj = phi2;
        }

        else if (nloci == 1) {
            dis = fabs(qtloc - mrklocn[0]) / 100.;
            if (map_func == 'h')
                theta = .5 * (1 - exp(-2*dis));
            else
                theta = .5 * (exp(4*dis) - 1) / (exp(4*dis) + 1);
            ibdadj = phi2 + (ibd[0] - phi2)*corr(class, theta);
        }

        else if (nloc == nloci && done[class]) {
            for (i = 0; i < nloc; i++)
                beta[i] = s_beta[class][i];

            betasum = 0;
            ibdadj = 0;
            for (i = 0; i < nloc; i++) {
                betasum += beta[i] * mu[class][mrk[i]];
                ibdadj += beta[i] * ibd[mrk[i]];
            }

            ibdadj += phi2 - betasum;
            if (ibdadj < 0) ibdadj = 0;
            if (ibdadj > 1) ibdadj = 1;
            if (classes[class]->cnstr && ibdadj > .5) ibdadj = .5;
        }

        else {
            for (i = 0; i < nloc; i++) {
                for (j = i; j < nloc; j++) {
                    dis = difmat(mrk[i],mrk[j]);
                    if (dis == 0) dis = 0.0001;
                    if (map_func == 'h')
                        theta = .5 * (1 - exp(-2*dis));
                    else
                        theta = .5 * (exp(4*dis) - 1) / (exp(4*dis) + 1);
                    if (i == j)
                        ri(i,j) = sd[class][mrk[i]] * sd[class][mrk[i]];
                    else {
                        ri(i,j) = (1 / classes[class]->evar)
                                   * sd[class][mrk[i]] * sd[class][mrk[i]]
                                   * sd[class][mrk[j]] * sd[class][mrk[j]]
                                   * corr(class, theta);
                        ri(j,i) = ri(i,j);
                    }
                }
            }

            dgefa_(ri, &nloci, &nloc, ipvt, &info);
            if (info) {
                printf("Matrix inversion failed, info = %d\n", info);
                exit(1);
            }
            dgedi_(ri, &nloci, &nloc, ipvt, det, work, &one);

            for (i = 0; i < nloc; i++) {
                dis = fabs(qtloc - mrklocn[mrk[i]]) / 100.;
                if (map_func == 'h')
                    theta = .5 * (1 - exp(-2*dis));
                else
                    theta = .5 * (exp(4*dis) - 1) / (exp(4*dis) + 1);
                cov[i] = sd[class][mrk[i]] * sd[class][mrk[i]]
                         * corr(class, theta);
            }

            for (i = 0; i < nloc; i++) {
                beta[i] = 0;
                for (j = 0; j < nloc; j++)
                    beta[i] += ri(i,j) * cov[j];
            }

            if (nloc == nloci) {
                done[class] = 1;
                for (i = 0; i < nloc; i++)
                    s_beta[class][i] = beta[i];
            }

            betasum = 0;
            ibdadj = 0;
            for (i = 0; i < nloc; i++) {
                betasum += beta[i] * mu[class][mrk[i]];
                ibdadj += beta[i] * ibd[mrk[i]];
            }

            ibdadj += phi2 - betasum;
            if (ibdadj < 0) ibdadj = 0;
            if (ibdadj > 1) ibdadj = 1;
            if (classes[class]->cnstr && ibdadj > .5) ibdadj = .5;
        }

        fprintf(mibdfp, "%5d %5d %10.7f %10.7f\n", id1, id2, ibdadj, phi2);
    }

    pipeback_shell_close(ibdfp);
    fclose(mibdfp);

    exit(0);
}

int getLocInfo (FILE *locfp, double *mrklocn, char *map_func)
{
    int nloci;
    char rec[100];

    if (!fgets(rec, sizeof(rec), locfp)) {
        printf("Read error on locus file, line 1\n");
        return 0;
    }

    if (!strncmp(rec, "NLOCI", 5)) {
        if (!fgets(rec, sizeof(rec), locfp)) {
            printf("Read error on locus file, line 2\n");
            return 0;
        }
        *map_func = 'k';
    }
    else {
        if (!strcmp(rec, "Kosambi\n"))
            *map_func = 'k';
        else if (!strcmp(rec, "Haldane\n"))
            *map_func = 'h';
        else {
            printf("Unrecognized mapping function\n");
            return 0;
        }
    }

    nloci = 0;
    while (fgets(rec, sizeof(rec), locfp)) {
        if (nloci == MLOCI) {
            printf("Too many loci, MLOCI = %d\n", MLOCI);
            return 0;
        }
        if (sscanf(&rec[13], "%lf", &mrklocn[nloci]) != 1 ||
            mrklocn[nloci] < 0)
        {
            printf("Read error on locus file, line %d\n", nloci + 3);
            return 0;
        }
        nloci++;
    }

    return nloci;
}

int getClasses (void)
{
    FILE *fp;
    int i, icl, icoeff, line;
    char *dir, fname[1000];
    char *recp, rec[10000];

    for (i = 0; i < MCLASS; i++)
        classes[i] = 0;

    dir = getenv("SOLAR_LIB");
    if (!dir) {
        printf("Environment variable SOLAR_LIB is not defined.\n");
        return 0;
    }

    sprintf(fname, "%s/classes.tab", dir);
    fp = fopen(fname, "r");
    if (!fp) {
        printf("Cannot open %s.\n", fname);
        return 0;
    }

    line = 0;
    while (fgets(rec, sizeof(rec), fp)) {
        line++;
        if (!(recp = strtok(rec, "|")) || sscanf(recp, "%d", &icl) != 1 ||
            icl < 0 || icl >= MCLASS)
        {
            printf("Illegal class number on line %d of %s\n", line, fname);
            return 0;
        }
        classes[icl] = (struct _class *) allocMem(sizeof(struct _class));

        if (!(recp = strtok(NULL, "|")))
        {
            printf("Illegal class name on line %d of %s\n", line, fname);
            return 0;
        }
        strcpy(classes[icl]->relation, recp);

        if (!(recp = strtok(NULL, "|")) ||
            sscanf(recp, "%d", &classes[icl]->ikin) != 1 ||
            classes[icl]->ikin < 0 || classes[icl]->ikin >= MCLASS)
        {
            printf("Illegal kinship class on line %d of %s\n", line, fname);
            return 0;
        }

        if (!(recp = strtok(NULL, "|")) ||
            sscanf(recp, "%d", &classes[icl]->cnstr) != 1 ||
            classes[icl]->cnstr < -1 || classes[icl]->cnstr > 1)
        {
            printf("Illegal constraint indicator on line %d of %s\n", line,
                   fname);
            return 0;
        }

        if (!(recp = strtok(NULL, "|")) ||
            !get_fraction(recp, &classes[icl]->phi2))
        {
            printf("Illegal kinship coefficient on line %d of %s\n", line,
                   fname);
            return 0;
        }

        if (!(recp = strtok(NULL, "|")) ||
            !get_fraction(recp, &classes[icl]->evar))
        {
            printf("Illegal expected variance on line %d of %s\n", line,
                   fname);
            return 0;
        }

        for (i = 0; i < NCOEFF; i++)
            classes[icl]->c[i] = 0;

        icoeff = 0;
        while (recp = strtok(NULL, "|")) {
            if (icoeff == NCOEFF)
            {
                printf("More than %d coefficients on line %d of %s\n",
                       NCOEFF, line, fname);
                return 0;
            }
            if (!get_fraction(recp, &classes[icl]->c[icoeff]))
            {
                printf("Coefficient #%d is illegal on line %d of %s\n",
                       icoeff + 1, line, fname);
                return 0;
            }
            icoeff++;
        }
    }

    fclose(fp);
    return 1;
}

int get_fraction (char *str, double *num)
{
    char *q, *p, cnum[64], fract[1024];
    double numer, denom;

    q = str;
    p = fract;
    while (*q && *q != '/') *p++ = *q++;
    *p = '\0';

    if (sscanf(fract, "%lf", &numer) != 1)
        return 0;

    if (!*q) {
        sprintf(cnum, "%.7g", numer);
        sscanf(cnum, "%lf", num);
        return 1;
    }

    q++;
    if (!*q || sscanf(q, "%lf", &denom) != 1 || denom == 0)
        return 0;

    sprintf(cnum, "%.7g", numer / denom);
    sscanf(cnum, "%lf", num);
    return 1;
}

int readMeanRec (char *rec, int line, int nloci, int *class, double *mu,
                 double *sd)
{
    int i;
    double phi2;
    char *recp;

    if (!(recp = strtok(rec, " \n")) || sscanf(recp, "%lf", &phi2) != 1) {
        printf("Illegal kinship coefficient on line %d of mean IBD file\n",
               line);
        return 0;
    }

    if (!(recp = strtok(NULL, " \n")) || sscanf(recp, "%d", class) != 1) {
        printf("Illegal class number on line %d of mean IBD file\n", line);
        return 0;
    }

    for (i = 0; i < nloci; i++) {
        if (!(recp = strtok(NULL, " \n")) ||
                sscanf(recp, "%lf", &mu[i]) != 1) {
            printf("Illegal class mean on line %d of mean IBD file\n", line);
            return 0;
        }
        if (!(recp = strtok(NULL, " \n")) ||
                sscanf(recp, "%lf", &sd[i]) != 1) {
            printf("Illegal class std dev on line %d of mean IBD file\n",
                   line);
            return 0;
        }
    }

    if (*class < 0 || *class >= MCLASS) {
        printf(
"Class is outside the range 0 to 999, line %d of mean IBD file\n", line);
        return 0;
    }

    if (!classes[*class]) {
        printf(
"Undefined class [%d], line %d of mean IBD file\n", *class, line);
        return 0;
    }

    return 1;
}

int readIbdRec (char *rec, int line, int nloci, int *id1, int *id2,
                int *class, double *phi2, double *ibd)
{
    int i;
    short typed;
    char *recp;

    if (!(recp = strtok(rec, " \n")) || sscanf(recp, "%d", id1) != 1) {
        printf("Read error on merged IBD file, line %d\n", line);
        return 0;
    }

    if (!(recp = strtok(NULL, " \n")) || sscanf(recp, "%d", id2) != 1) {
        printf("Read error on merged IBD file, line %d\n", line);
        return 0;
    }

    if (!(recp = strtok(NULL, " \n")) || sscanf(recp, "%lf", phi2) != 1) {
        printf("Read error on merged IBD file, line %d\n", line);
        return 0;
    }

    if (!(recp = strtok(NULL, " \n")) || sscanf(recp, "%d", class) != 1) {
        printf("Read error on merged IBD file, line %d\n", line);
        return 0;
    }

    recp = strtok(NULL, " \n");
    for (i = 0; i < nloci; i++) {
        if (!recp || sscanf(recp, "%lf", &ibd[i]) != 1) {
            printf("Read error on merged IBD file, line %d\n", line);
            return 0;
        }

        recp = strtok(NULL, " \n");
        if (recp && strlen(recp) == 1) {
            if (sscanf(recp, "%hd", &typed) != 1 ||
                (typed != 0 && typed != 1))
            {
                printf("Read error on merged IBD file, line %d\n", line);
                return 0;
            }
            recp = strtok(NULL, " \n");
        }
    }

    if (*class < 0 || *class >= MCLASS) {
        printf(
"Class is outside the range 0 to 999, line %d of merged IBD file\n", line);
        return 0;
    }

    if (!classes[*class]) {
        printf(
"Undefined class [%d], line %d of merged IBD file\n", *class, line);
        return 0;
    }

    if (classes[*class]->cnstr == -1) {
        printf("Multipoint not supported for class %d, relation \"%s\"\n",
               *class, classes[*class]->relation);
        return 0;
    }

    if (*phi2 != classes[*class]->phi2)
    {
        printf(
"Pair's phi2 does not match class phi2, line %d of merged IBD file\n", line);
        return 0;
    }

    return 1;
}

double corr (int class, double theta)
{
    int i;
    double t, corr;

    corr = 0;
    t = 1;
    for (i = 0; i < NCOEFF; i++) {
        corr += classes[class]->c[i] * t;
        t *= theta;
    }

    return corr;
}

void *allocMem (size_t nbytes)
{
    void *ptr;
    ptr = (void *) malloc(nbytes);
    if (!ptr) {
        printf("Not enough memory.\n");
        exit(1);
    }
    return ptr;
}
