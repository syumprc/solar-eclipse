#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "pipeback.h"

#define NCOEFF 14
#define MCLASS 1000

struct _class {
    int ikin;
    int cnstr;
    char relation[1000];
    double phi2;
    double evar;
    double c[NCOEFF];
} *classes[MCLASS];

void *allocMem (size_t);

main (int argc, char **argv)
{
    int nloci, line, class;
    int i, j, k;
    double phi2, *ibd;
    short *typed;
    int typed_only;
    double *xn[MCLASS], *xsum[MCLASS], *x2sum[MCLASS];
    double *xmu1[MCLASS], *xsd1[MCLASS];
    char rec[100000];
    const char *pbarg[4];
    FILE *ibdfp, *meanfp;

    if (argc != 5) {
        printf("usage: getmeans ibd_file mean_file nloci typed_only(y/n)\n");
        exit(1);
    }

    ibdfp = fopen(argv[1], "r");
    if (!ibdfp) {
        printf(
"Merged IBD file not found. Use the command \"mibd merge\" to create it.\n");
        exit(1);
    }
    fclose(ibdfp);

    pbarg[0] = "gunzip";
    pbarg[1] = "-c";
    pbarg[2] = argv[1];
    pbarg[3] = 0;

    ibdfp = pipeback_shell_open("gunzip", pbarg);
    if (!ibdfp) {
        printf("Cannot uncompress merged IBD file.");
        exit(1);
    }

    if (sscanf(argv[3], "%d", &nloci) != 1) {
        printf("Invalid number of loci [%s]\n", argv[3]);
        exit(1);
    }

    typed_only = 0;
    if (argv[4][0] == 'y')
        typed_only = 1;

    if (!getClasses()) {
        exit(1);
    }

    ibd = (double *) allocMem(nloci*sizeof(double));
    typed = (short *) allocMem(nloci*sizeof(short));

    for (i = 0; i < MCLASS; i++) {
        xn[i] = (double *) allocMem(nloci*sizeof(double));
        xsum[i] = (double *) allocMem(nloci*sizeof(double));
        x2sum[i] = (double *) allocMem(nloci*sizeof(double));
        xmu1[i] = (double *) allocMem(nloci*sizeof(double));
        xsd1[i] = (double *) allocMem(nloci*sizeof(double));
        for (j = 0; j < nloci; j++) {
            xn[i][j] = 0;
            xsum[i][j] = 0;
            x2sum[i][j] = 0;
            xmu1[i][j] = 0;
            xsd1[i][j] = 0;
        }
    }

    line = 0;
    while (fgets(rec, sizeof(rec), ibdfp)) {
        line++;
        if (!readIbdRec(rec, line, nloci, &class, &phi2, ibd, typed))
            exit(1);
        k = classes[class]->ikin;
        for (i = 0; i < nloci; i++) {
            if (typed_only && typed[i] && ibd[i] != -1.) {
                xsum[k][i] += ibd[i];
                x2sum[k][i] += ibd[i]*ibd[i];
                xn[k][i]++;
            }
            else if (!typed_only && ibd[i] != -1.) {
                xsum[k][i] += ibd[i];
                x2sum[k][i] += ibd[i]*ibd[i];
                xn[k][i]++;
            }
        }
    }

    pipeback_shell_close(ibdfp);

    for (i = 0; i < MCLASS; i++) {
        if (!classes[i]) continue;
        k = classes[i]->ikin;
        for (j = 0; j < nloci; j++) {
            if (xn[k][j] < 100 || xsum[k][j] == 0) {
                xmu1[i][j] = classes[i]->phi2;
                xsd1[i][j] = sqrt(classes[i]->evar);
            }
            else {
                xmu1[i][j] = xsum[k][j] / xn[k][j];
                xsd1[i][j] = x2sum[k][j] / xn[k][j] - xmu1[i][j] * xmu1[i][j];
                if (xsd1[i][j] > 0)
                    xsd1[i][j] = sqrt(xsd1[i][j]);
                else
                    xsd1[i][j] = sqrt(classes[i]->evar);
            }
        }
    }

    meanfp = fopen(argv[2], "w");
    if (!meanfp) {
        printf("Cannot create mean IBD file [%s]\n", argv[2]);
        exit(1);
    }

    for (i = 0; i < MCLASS; i++) {
        if (!classes[i]) continue;

        if (i == 1) {
            for (j = 0; j < nloci; j++) {
                if (xmu1[i][j] != 1 || xsd1[i][j] != 0) {
                    printf(
"The mean IBDs for class = 1 (self) are not valid.\n\
Make sure you are using the correct marker-specific IBDs.\n");
                    fclose(meanfp);
                    unlink(argv[2]);
                    exit(1);
                }
            }
        }

        if (i == 2) {
            for (j = 0; j < nloci; j++) {
                if (xmu1[i][j] != 0.5 || xsd1[i][j] != 0) {
                    printf(
"The mean IBDs for parent-offspring pairs are not valid.\n\
Make sure you are using the correct marker-specific IBDs.\n");
                    fclose(meanfp);
                    unlink(argv[2]);
                    exit(1);
                }
            }
        }

        fprintf(meanfp, "%.7g %d", classes[i]->phi2, i);
        for (j = 0; j < nloci; j++)
            fprintf(meanfp, " %.7g %.7g", xmu1[i][j], xsd1[i][j]);
        fprintf(meanfp, "\n");
    }

    fclose(meanfp);
    exit(0);
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

int readIbdRec (char *rec, int line, int nloci, int *class, double *phi2,
                double *ibd, short *typed)
{
    int id, i;
    char *recp;

    if (!(recp = strtok(rec, " \n")) || sscanf(recp, "%d", &id) != 1) {
        printf("Read error on merged IBD file, line %d\n", line);
        return 0;
    }

    if (!(recp = strtok(NULL, " \n")) || sscanf(recp, "%d", &id) != 1) {
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
            if (sscanf(recp, "%hd", &typed[i]) != 1 ||
                (typed[i] != 0 && typed[i] != 1))
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
        printf("Undefined class [%d], line %d of merged IBD file\n",
               *class, line);
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
