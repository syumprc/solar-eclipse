#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pipeback.h"

#define MAXLOC 3000
#define MMRKNM 20

main (int argc, char **argv)
{
    FILE *infp, *outfp;
    int i, loc, nloci, npairs, nids, cdfile;
    char chrnum[1024], fmt[32];
    int *id1, *id2, *class;
    double *kin2, *ibd[MAXLOC];
    short *typed[MAXLOC];
    unsigned char *missing[MAXLOC];
    int tid1, tid2;
    double tibd, tdlt7;
    short ttype;
    char mrknam[MAXLOC][MMRKNM+1];
    char rec[1024], fname[1024], cmd[1024];
    const char *pbarg[4];

    if (argc != 4) {
        printf("usage: mrgibd mapfile ibddir mibddir\n");
        exit(1);
    }

    infp = fopen(argv[1], "r");
    if (!infp) {
        printf("Cannot open map-data file \"%s\".\n", argv[1]);
        exit(1);
    }

    if (!fgets(rec, sizeof(rec), infp) || sscanf(rec, "%s", chrnum) != 1) {
        printf("Error reading map-data file.\n");
        exit(1);
    }

    nloci = 0;
    while (fgets(rec, sizeof(rec), infp)) {
        if (nloci >= MAXLOC) {
            printf("Too many loci for mrgibd: MAXLOC = %d\n", MAXLOC);
            exit(1);
        }
        if (sscanf(rec, "%s ", mrknam[nloci]) != 1) {
            printf("Error reading map-data file.\n");
            exit(1);
        }
        nloci++;
    }
    fclose(infp);

    infp = fopen("mibdrel.ped", "r");
    if (!infp) {
        printf("Relative-class file not found. ");
        printf("Use the command \"mibd relate\" to create it.\n");
        exit(1);
    }

    fgets(rec, sizeof(rec), infp);
    if (!strncmp(rec, "IBDID1", 6)) {
        npairs = 0;
        cdfile = 1;
    }
    else {
        npairs = 1;
        cdfile = 0;
    }

    while (fgets(rec, sizeof(rec), infp))
        npairs++;
    if (!npairs) {
        printf("Relative-class file is empty. ");
        printf("Use the command \"mibd relate\" to create it.\n");
        exit(1);
    }

    id1 = (int *) malloc(npairs*sizeof(int));
    if (!id1) {
        printf("mrgibd: Not enough memory.\n");
        exit(1);
    }

    id2 = (int *) malloc(npairs*sizeof(int));
    if (!id2) {
        printf("mrgibd: Not enough memory.\n");
        exit(1);
    }

    kin2 = (double *) malloc(npairs*sizeof(double));
    if (!kin2) {
        printf("mrgibd: Not enough memory.\n");
        exit(1);
    }

    class = (int *) malloc(npairs*sizeof(int));
    if (!class) {
        printf("mrgibd: Not enough memory.\n");
        exit(1);
    }

    for (loc = 0; loc < nloci; loc++) {
        ibd[loc] = (double *) malloc(npairs*sizeof(double));
        if (!ibd[loc]) {
            printf("mrgibd: Not enough memory.\n");
            exit(1);
        }
        for (i = 0; i < npairs; i++)
            ibd[loc][i] = 0;
    }

    rewind(infp);
    strcpy(fmt, "%d %d %lf %d\n");
    if (cdfile) {
        fgets(rec, sizeof(rec), infp);
        strcpy(fmt, "%d,%d,%lf,%d\n");
    }

    nids = 0;
    for (i = 0; i < npairs; i++) {
        if (fscanf(infp, fmt, id1+i, id2+i, kin2+i, class+i) != 4) {
            printf("Error reading relative-class file.\n");
            exit(1);
        }
        if (id1[i] == id2[i]) nids++;
    }
    fclose(infp);

    for (loc = 0; loc < nloci; loc++) {
        missing[loc] = (unsigned char *) malloc(nids*sizeof(unsigned char));
        if (!missing[loc]) {
            printf("mrgibd: Not enough memory.\n");
            exit(1);
        }
        for (i = 0; i < nids; i++)
            missing[loc][i] = (unsigned char) 0;

        typed[loc] = (short *) malloc(npairs*sizeof(short));
        if (!typed[loc]) {
            printf("mrgibd: Not enough memory.\n");
            exit(1);
        }
        for (i = 0; i < npairs; i++)
            typed[loc][i] = (short) 1;
    }

    for (loc = 0; loc < nloci; loc++) {
        sprintf(fname, "%s/ibd.%s.gz", argv[2], mrknam[loc]);
        infp = fopen(fname, "r");
        if (!infp) {
            printf("Cannot find marker-specific IBDs for marker %s.\n",
                   mrknam[loc]);
            exit(1);
        }
        fclose(infp);

        pbarg[0] = "gunzip";
        pbarg[1] = "-c";
        pbarg[2] = fname;
        pbarg[3] = 0;
        infp = pipeback_shell_open("gunzip", pbarg);
        if (!infp) {
            printf("Cannot uncompress IBD file for marker %s.\n",
                   mrknam[loc]);
            exit(1);
        }

        i = 0;
        while (fgets(rec, sizeof(rec), infp)) {
            if (sscanf(rec, "%d %d %lf %lf %hd\n",
                       &tid1, &tid2, &tibd, &tdlt7, &ttype) != 5)
            {
                if (sscanf(rec, "%d %d %lf %lf\n",
                           &tid1, &tid2, &tibd, &tdlt7) != 4)
                {
                    printf("Error reading IBD file \"%s\".\n", fname);
                    exit(1);
                }
                ttype = (short) 1;
            }
            while (i < npairs && (id1[i] != tid1 || id2[i] != tid2))
                i++;
            if (i == npairs) {
                printf("Error reading IBD file \"%s\".\n", fname);
                exit(1);
            }
            ibd[loc][i] = tibd;
            typed[loc][i] = ttype;
            if (tibd == -1.)
                missing[loc][id1[i]-1] = (unsigned char) 1;
        }
        if (tid1 != tid2 || tid1 != id1[npairs-1]) {
            printf("Error reading IBD file \"%s\".\n", fname);
            exit(1);
        }
        pipeback_shell_close(infp);
    }

    sprintf(fname, "%s/mibdchr%s.mrg", argv[3], chrnum);
    outfp = fopen(fname, "w");
    if (!outfp) {
        printf("Cannot open merged IBD file \"%s\".\n", fname);
        exit(1);
    }

    for (i = 0; i < npairs; i++) {

        /* skip pairs with kinship = 0 */
        if (!kin2[i])
            continue;

/*        fprintf(outfp, "%5d %5d %12.10lf %3d",*/
        fprintf(outfp, "%d %d %.7g %d",
                       id1[i], id2[i], kin2[i], class[i]);
        for (loc = 0; loc < nloci; loc++) {
            if (missing[loc][id1[i]-1] || missing[loc][id2[i]-1])
/*                fprintf(outfp, " %10.7f %1d", -1., 0);*/
                fprintf(outfp, " %.7g %d", -1., 0);
            else
/*                fprintf(outfp, " %10.7f %1d", ibd[loc][i], typed[loc][i]);*/
                fprintf(outfp, " %.7g %d", ibd[loc][i], typed[loc][i]);
        }
        fprintf(outfp, "\n");
    }
    fclose(outfp);

    sprintf(cmd, "gzip -f %s", fname);
    if (system(cmd)) {
        printf("merge: gzip failed\n");
        exit(1);
    }

    exit(0);
}
