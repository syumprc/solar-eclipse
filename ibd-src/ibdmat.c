#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>

#define MAXPED 10000
#define MXTWIN 999
#define MAXBRK 1

int nind[MAXPED];
int twin[MXTWIN+1];
int brk1[MAXBRK], brk2[MAXBRK];
int *id;

int nids, ntot, nbrk, ntyped, currped;
int xlinked, showstat;

struct _ego {
    int fa;
    int mo;
    int kid1;
    int psib;
    int msib;
    int sex;
    int prbnd;
    int tid;
    int a1;
    int a2;
    int rfa;
    int rmo;
    int same;
} *ego;

char *solarbin;
char cwd[1024];
char ibdfile[1024];
FILE *ibdfp;

double *kin2;
float *ibd, *delta7;
short *typed;

#define kin2(a, b)      kin2[((a)*nids) + (b)]
#define ibd(a, b)       ibd[((a)*nids) + (b)]
#define delta7(a, b)    delta7[((a)*nids) + (b)]
#define typed(a, b)     typed[((a)*nids) + (b)]

#define max(a, b)       ( (a) > (b) ? (a) : (b) )
#define min(a, b)       ( (a) < (b) ? (a) : (b) )

void add_ego(int, int, int, int, int, int, int, int, int, int, int);
void processPed (void);
void calcKin2 (void);
int readPedInfo (void);
void create_datafile (FILE*);
void *allocMem (size_t);
void fatal_error (void);

main (int argc, char **argv)
{
    int i, line, lastped;
    int iped, iid, ifa, imo, ikid1, ipsib, imsib, isex, iprbnd;
    int ia1, ia2, itwin, junk;
    char buf[1024];
    FILE *infp, *ttyfp;

    if (argc != 4) {
        printf("Usage: ibdmat xLinked(y/n) showStatus(y/n) ibdFile\n");
        exit(1);
    }

    getcwd(cwd, sizeof(cwd));

    switch (argv[1][0]) {
    case 'y':
        xlinked = 1;
        break;
    case 'n':
        xlinked = 0;
        break;
    default:
        printf("Enter y or n for xLinked.\n");
        exit(1);
    }

    switch (argv[2][0]) {
    case 'y':
        showstat = 1;
        break;
    case 'n':
        showstat = 0;
        break;
    default:
        printf("Enter y or n for showStat.\n");
        exit(1);
    }

    strcpy(ibdfile, argv[3]);
    ibdfp = fopen(ibdfile, "w");
    if (!ibdfp) {
        printf("Cannot open IBD file: %s\n", ibdfile);
        exit(1);
    }

    solarbin =  getenv("SOLAR_BIN");
    if (!solarbin) {
        printf("Environment variable SOLAR_BIN not defined.\n");
        fatal_error();
    }

    if (!readPedInfo())
        fatal_error();

    unlink("checkpoint.MLINK");

    infp = fopen("datafile.dat", "r");
    if (!infp) {
        infp = fopen("datain.dat", "r");
        if (!infp) {
            printf("Cannot open LINKAGE input file \"datafile.dat\".\n");
            fatal_error();
        }
        create_datafile(infp);
    }
    fclose(infp);

    infp = fopen("pedin.dat", "r");
    if (!infp) {
        printf("Cannot open LINKAGE input file \"pedin.dat\".\n");
        fatal_error();
    }

    if (showstat) {
        ttyfp = fopen("/dev/tty", "w");
        if (!infp) {
            printf("Cannot open /dev/tty\n");
            fatal_error();
        }
    }

    currped = 0;
    lastped = 0;
    line = 0;
    ntyped = 0;
    ntot = 0;
    nbrk = 0;
    brk1[0] = 0;
    brk2[0] = 0;
    for (i = 0; i < MXTWIN+1; i++)
        twin[i] = 0;

    while (!lastped) {
        if (!fgets(buf, sizeof(buf), infp)) {
            lastped = 1;
            fclose(infp);
        }
        else if (sscanf(buf, "%d %d %d %d %d %d %d %d %d %d %d %d %d",
                        &iped, &iid, &ifa, &imo, &ikid1, &ipsib, &imsib,
                        &isex, &iprbnd, &junk, &itwin, &ia1, &ia2) != 13)
        {
            printf("Error reading pedin.dat, line %d\n", line + 1);
            fatal_error();
        }
        else {
            line++;
            if (!currped) {
                currped = iped;
                ego = (struct _ego *)
                      allocMem(nind[currped-1] * sizeof(struct _ego));
                id = (int *) allocMem(nind[currped-1] * sizeof(int));
            }

            if (!ia1 && ia2 || ia1 && !ia2) {
                printf(
                "Only one allele typed for individual %d, pedigree %d\n",
                       iid, iped);
                fatal_error();
            }

            if (itwin > MXTWIN) {
                printf("Too many sets of twins, MXTWIN = %d\n", MXTWIN);
                fatal_error();
            }

            if (iprbnd == 2 && nbrk == MAXBRK) {
                printf("More than %d loopbreakers in pedigree %d\n",
                       MAXBRK, iped);
                fatal_error();
            }
        }

        if (!lastped && iped == currped) {
            if (ntot == nind[currped-1]) {
                printf("More individuals than expected, pedigree %d\n",
                       currped);
                fatal_error();
            }

            --iid;
            add_ego(iid, ifa, imo, ikid1, ipsib, imsib, isex, iprbnd, ia1,
                    ia2, itwin);
        }
        else {
            if (showstat) {
                if (currped > 1)
                    fprintf(ttyfp, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
                fprintf(ttyfp, "pedigree%5d ", currped);
                fflush(ttyfp);
            }

            processPed();

            if (!lastped) {
                if (ntot != nind[currped-1]) {
                    printf("Fewer individuals than expected, pedigree %d\n",
                           currped);
                    fatal_error();
                }

                currped = iped;
                free(ego);
                ego = (struct _ego *)
                      allocMem(nind[currped-1] * sizeof(struct _ego));
                free(id);
                id = (int *) allocMem(nind[currped-1] * sizeof(int));

                ntyped = 0;
                ntot = 0;
                nbrk = 0;
                brk1[0] = 0;
                brk2[0] = 0;
                for (i = 0; i < MXTWIN+1; i++)
                    twin[i] = 0;

                --iid;
                add_ego(iid, ifa, imo, ikid1, ipsib, imsib, isex, iprbnd,
                        ia1, ia2, itwin);
            }
        }
    }

    if (showstat)
        fclose(ttyfp);
    fclose(ibdfp);

    sprintf(buf, "gzip -f %s", ibdfile);
    if (system(buf)) {
        printf("merge: gzip failed\n");
        exit(1);
    }

    exit(0);
}

void add_ego(int iid, int ifa, int imo, int ikid1, int ipsib, int imsib,
             int isex, int iprbnd, int ia1, int ia2, int itwin)
{
    ego[iid].fa = ifa;
    ego[iid].mo = imo;
    ego[iid].kid1 = ikid1;
    ego[iid].psib = ipsib;
    ego[iid].msib = imsib;
    ego[iid].sex = isex;
    ego[iid].a1 = ia1;
    ego[iid].a2 = ia2;

    if (iprbnd == 1)
        ego[iid].prbnd = 0;
    else
        ego[iid].prbnd = iprbnd;

    if (itwin) {
        if (twin[itwin-1])
            ego[iid].tid = twin[itwin-1];
        else {
            twin[itwin-1] = iid + 1;
            ego[iid].tid = 0;
        }
    }
    else
        ego[iid].tid = 0;

    if (iprbnd == 2) {
        if (brk1[nbrk]) {
            brk2[nbrk] = iid + 1;
            nbrk++;
        }
        else
            brk1[nbrk] = iid + 1;
    }

    if ((iprbnd != 2 || !brk2[nbrk-1]) && ia1 && ia2)
        ntyped++;

    id[ntot] = iid;
    ntot++;
}

void processPed (void)
{
    static id0 = 1;
    static char *fmt = "%d %d %d %d %d %d %d %d %d %d %d %d\n";

    int i, j, k, ii, kk, done;
    int ifa, imo;
    float p1, p2, ptmp;
    char buf[1024];
    FILE *fp;

    for (i = 0; i < ntot; i++) {
        ego[i].rfa = ego[i].fa;
        ego[i].rmo = ego[i].mo;
        for (j = 0; j < nbrk; j++) {
            if (ego[i].fa == brk2[j]) ego[i].rfa = brk1[j];
            if (ego[i].mo == brk2[j]) ego[i].rmo = brk1[j];
        }

        ego[i].same = 0;
        if (!xlinked && ego[i].rfa && ego[i].a1) {
            for (j = 0; j < i; j++) {
                if (ego[j].rfa == ego[i].rfa && ego[j].rmo == ego[i].rmo &&
                    (ego[j].a1 == ego[i].a1 && ego[j].a2 == ego[i].a2 ||
                     ego[j].a1 == ego[i].a2 && ego[j].a2 == ego[i].a1))
                {
                    ego[i].same = j;
                }
            }
        }
    }

    nids = ntot - nbrk;
    kin2 = (double *) allocMem(nids * nids * sizeof(double));
    calcKin2();

    ibd = (float *) allocMem(nids * nids * sizeof(float));
    delta7 = (float *) allocMem(nids * nids * sizeof(float));
    typed = (short *) allocMem(nids * nids * sizeof(short));

    for (i = 0; i < nids; i++) {
        if (ego[i].a1)
            typed(i,i) = (short) 1;
        else
            typed(i,i) = (short) 0;

        if (!ntyped) {
            ibd(i,i) = -1;
            delta7(i,i) = -1;
        }
        else {
            ibd(i,i) = 1;
            delta7(i,i) = 1;
        }

        if (ego[i].rfa) {
            ifa = ego[i].rfa - 1;
            imo = ego[i].rmo - 1;
        }

        for (j = 0; j < i; j++) {
            if (ego[i].a1 && ego[j].a1)
                typed(i,j) = (short) 1;
            else
                typed(i,j) = (short) 0;

            if (!ntyped) {
                if (ego[i].tid == j+1) {
                    ibd(i,j) = -1;
                    delta7(i,j) = -1;
                }
                else {
                    ibd(i,j) = 0;
                    delta7(i,j) = 0;
                }
            }

            else if (!xlinked && ego[i].rfa &&
                     (ego[i].rfa == j+1 || ego[i].rmo == j+1 ||
                      ego[ifa].tid == j+1 || ego[imo].tid == j+1 ||
                      (ego[ifa].tid && ego[ifa].tid == ego[j].tid) ||
                      (ego[imo].tid && ego[imo].tid == ego[j].tid)))
            {
                ibd(i,j) = 0.5;
                delta7(i,j) = 0;
            }

            else if (ego[i].tid && !ego[j].tid) {
                ibd(i,j) = ibd(max(ego[i].tid-1,j), min(ego[i].tid-1,j));
                delta7(i,j) = delta7(max(ego[i].tid-1,j), min(ego[i].tid-1,j));
            }

            else if (ego[i].tid && ego[j].tid) {
                ibd(i,j) = ibd(max(ego[i].tid-1,ego[j].tid-1),
                               min(ego[i].tid-1,ego[j].tid-1));
                delta7(i,j) = delta7(max(ego[i].tid-1,ego[j].tid-1),
                                     min(ego[i].tid-1,ego[j].tid-1));
            }

            else if (ego[j].tid) {
                ibd(i,j) = ibd(max(i,ego[j].tid-1), min(i,ego[j].tid-1));
                delta7(i,j) = delta7(max(i,ego[j].tid-1), min(i,ego[j].tid-1));
            }

            else if (ego[i].rfa &&
                     !ibd(max(ego[i].rfa-1,j), min(ego[i].rfa-1,j)) &&
                     !ibd(max(ego[i].rmo-1,j), min(ego[i].rmo-1,j)))
            {
                ibd(i,j) = 0;
                delta7(i,j) = 0;
            }

            else if (j < ego[i].same &&
                     (!ego[ego[i].same].tid ||
                         ego[ego[i].same].tid != j+1 &&
                         ego[ego[i].same].tid != ego[j].tid))
            {
                ibd(i,j) = ibd(ego[i].same,j);
                delta7(i,j) = delta7(ego[i].same,j);
            }

            else if (ego[j].same && ibd(j,ego[j].same) == 1) {
                ibd(i,j) = ibd(i,ego[j].same);
                delta7(i,j) = delta7(i,ego[j].same);
            }

            else if (kin2(i,j) &&
                     (!ego[i].a1 || !ego[j].a1 ||
                      ego[i].a1 == ego[j].a1 || ego[i].a1 == ego[j].a2 ||
                      ego[i].a2 == ego[j].a1 || ego[i].a2 == ego[j].a2))
            {
                fp = fopen("pedfile.dat", "w");
                if (!fp) {
                    printf("Cannot open pedfile.dat\n");
                    fatal_error();
                }

                for (k = 0; k < ntot; k++) {
                    kk = id[k];

                    if (kk == i) {
                        fprintf(fp, fmt, currped, kk+1, ego[kk].fa,
                                ego[kk].mo, ego[kk].kid1, ego[kk].psib,
                                ego[kk].msib, ego[kk].sex, 1, 0, ego[kk].a1,
                                ego[kk].a2);
                    }
                    else if (kk == j) {
                        fprintf(fp, fmt, currped, kk+1, ego[kk].fa,
                                ego[kk].mo, ego[kk].kid1, ego[kk].psib,
                                ego[kk].msib, ego[kk].sex, ego[kk].prbnd, 2,
                                ego[kk].a1, ego[kk].a2);
                    }
                    else {
                        done = 0;
                        for (ii = 0; !done && ii < nbrk; ii++) {
                            if (brk1[ii] == i+1 && brk2[ii] == kk+1) {
                                fprintf(fp, fmt, currped, kk+1, ego[kk].fa,
                                        ego[kk].mo, ego[kk].kid1, ego[kk].psib,
                                        ego[kk].msib, ego[kk].sex, 2, 0,
                                        ego[kk].a1, ego[kk].a2);
                                done = 1;
                            }
                        }

                        for (ii = 0; !done && ii < nbrk; ii++) {
                            if (brk1[ii] == j+1 && brk2[ii] == kk+1) {
                                fprintf(fp, fmt, currped, kk+1, ego[kk].fa,
                                        ego[kk].mo, ego[kk].kid1, ego[kk].psib,
                                        ego[kk].msib, ego[kk].sex, 2, 2,
                                        ego[kk].a1, ego[kk].a2);
                                done = 1;
                            }
                        }

                        if (!done) {
                            fprintf(fp, fmt, currped, kk+1, ego[kk].fa,
                                    ego[kk].mo, ego[kk].kid1, ego[kk].psib,
                                    ego[kk].msib, ego[kk].sex, ego[kk].prbnd,
                                    0, ego[kk].a1, ego[kk].a2);
                        }
                    }
                }

                fclose(fp);
                sprintf(buf, "%s/unknown > unknown.out 2>&1", solarbin);
                if (system(buf)) {
                printf("Error in program unknown: See %s/unknown.out.\n", cwd);
                    fatal_error();
                }

                sprintf(buf, "%s/mlink > mlink.out 2>&1", solarbin);
                if (system(buf)) {
                    printf("Error in program mlink: See %s/mlink.out.\n", cwd);
                    fatal_error();
                }

                fp = fopen("mlink.out", "r");
                if (!fp) {
                    printf("Cannot open mlink.out\n");
                    fatal_error();
                }

                ptmp = -1;
                while (fgets(buf, sizeof(buf), fp)) {
                    if (!strncmp(buf, "HOMOZYGOTE CARRIER", 18)
                        || !strncmp(buf, "HETEROZYGOTE CARRIER", 20)
                        || !strncmp(buf, "MALE CARRIER", 12))
                    {
                        sscanf(buf, "%*[A-Z: ]%f", &ptmp);
                        break;
                    }
                }

                if (ptmp == -1) {
                    printf("Error reading mlink.out\n");
                    fatal_error();
                }

                if (xlinked && ego[i].sex == 1) {
                    p1 = ptmp;
                    p2 = 0;
                }
                else {
                    if (!fgets(buf, sizeof(buf), fp)
                        || strncmp(buf, "HETEROZYGOTE CARRIER", 20))
                    {
                        printf("Error reading mlink.out\n");
                        fatal_error();
                    }
                    sscanf(buf, "%*[A-Z: ]%f", &p1);
                    p2 = ptmp;
                }

                fclose(fp);
                unlink("mlink.out");
                unlink("unknown.out");

                if (!xlinked || ego[i].sex == 2 && ego[j].sex == 2)
                    ibd(i,j) = .5 * p1 + p2;
                else if (ego[i].sex == 1 && ego[j].sex == 1)
                    ibd(i,j) = p1;
                else
                    ibd(i,j) = .70710678 * p1;

                delta7(i,j) = p2;
            }

            else {
                ibd(i,j) = 0;
                delta7(i,j) = 0;
            }

            if (ibd(i,j))
                fprintf(ibdfp, "%5d %5d %10.7f %10.7f %1d\n", id0 + i, id0 + j,
                        ibd(i,j), delta7(i,j), typed(i,j));
        }

        fprintf(ibdfp, "%5d %5d %10.7f %10.7f %1d\n", id0 + i, id0 + i,
                ibd(i,i), delta7(i,i), typed(i,i));
    }

    id0 += nids;

    unlink("ipedfile.dat");
    unlink("loopfile.dat");
    unlink("newspeedfile.dat");
    unlink("outfile.dat");
    unlink("pedfile.dat");
    unlink("speedfile.dat");
    unlink("stream.dat");
    unlink("FASTLINK.err");
    unlink("checkpoint.MLINK");

    free(typed);
    free(delta7);
    free(ibd);
    free(kin2);
}

void calcKin2 (void)
{
    int i, j, ifa, imo, count;

    count = 0;
    for (i = 0; i < nids; i++) {
        for (j = 0; j < i; j++)
            kin2(i,j) = 0;
        if (!ego[i].rfa) {
            count++;
            kin2(i,i) = 1;
        }
        else
            kin2(i,i) = 0;
    }

    do {
        for (i = 0; i < nids; i++) {
            if (kin2(i,i))
                continue;
            if (!ego[i].rfa)
                continue;

            ifa = ego[i].rfa - 1;
            imo = ego[i].rmo - 1;
            if (!kin2(ifa,ifa) || !kin2(imo,imo))
                continue;

            for (j = 0; j < nids; j++) {
                if (!kin2(j,j))
                    continue;
                kin2(max(i,j), min(i,j)) =
                               .5 * (kin2(max(ifa,j), min(ifa,j)) +
                                     kin2(max(imo,j), min(imo,j)));
            }
            count++;
            kin2(i,i) = 1 + .5 * kin2(max(ifa,imo), min(ifa,imo));
        }
    } while (count < nids);

    for (i = 0; i < nids; i++)
        for (j = 0; j < i; j++)
            kin2(j,i) = kin2(i,j);
}

int readPedInfo (void)
{
    int i, junk, nloop;
    char buf[1024];
    FILE *infp;

    infp = fopen("../pedigree.info", "r");
    if (!infp) {
        printf("Cannot open pedigree.info\n");
        return 0;
    }

    if (!fgets(buf, sizeof(buf), infp)) {
        printf("Error reading pedigree.info\n");
        return 0;
    }
    if (!fgets(buf, sizeof(buf), infp)) {
        printf("Error reading pedigree.info\n");
        return 0;
    }
    if (!fgets(buf, sizeof(buf), infp)) {
        printf("Error reading pedigree.info\n");
        return 0;
    }

    i = 0;
    while (fgets(buf, sizeof(buf), infp)) {
        if (i == MAXPED) {
            printf("Too many pedigrees, MAXPED = %d\n", MAXPED);
            return 0;
        }
        if (sscanf(buf, "%d %d %d %d", &junk, &nind[i], &junk, &nloop) != 4
            || nind[i] < 1 || nloop < 0)
        {
            printf("Error reading pedigree.info\n");
            return 0;
        }
        nind[i] += nloop;
        i++;
    }

    fclose(infp);
    return 1;
}

void create_datafile (FILE *infp)
{
    int i, junk, nall, xlinked;
    double *freq;
    char *p, buf[1024];

    FILE *outfp = fopen("datafile.dat", "w");
    if (!outfp) {
        printf("Cannot create LINKAGE input file \"datafile.dat\".\n");
        fatal_error();
    }

    if (!fgets(buf, sizeof(buf), infp) ||
        sscanf(buf, "%d %d %d %d", &junk, &junk, &xlinked, &junk) != 4)
    {
        printf("Error reading datain.dat\n");
        fatal_error();
    }

    for (i = 0; i < 7; i++) {
        if (!fgets(buf, sizeof(buf), infp)) {
            printf("Error reading datain.dat\n");
            fatal_error();
        }
    }

    if (!fgets(buf, sizeof(buf), infp) ||
        sscanf(buf, "%d %d", &junk, &nall) != 2)
    {
        printf("Error reading datain.dat\n");
        fatal_error();
    }

    if (!fgets(buf, sizeof(buf), infp)) {
        printf("Error reading datain.dat\n");
        fatal_error();
    }

    freq = (double *) malloc(nall*sizeof(double));
    if (!freq) {
        printf("Not enough memory.\n");
        fatal_error();
    }

    if (!(p = strtok(buf, " \t\n")) || sscanf(p, "%lf", &freq[0]) != 1) {
        printf("Error reading datain.dat\n");
        fatal_error();
    }

    for (i = 1; i < nall; i++) {
        if (!(p = strtok(NULL, " \t\n")) || sscanf(p, "%lf", &freq[i]) != 1) {
            printf("Error reading datain.dat\n");
            fatal_error();
        }
    }

    if (xlinked)
        fprintf(outfp, "2 1 1 5\n");
    else
        fprintf(outfp, "2 1 0 5\n");

    fprintf(outfp, "0 0.00000000 0.00000000 0\n");
    fprintf(outfp, " 1 2\n");
    fprintf(outfp, "1 2\n");
    fprintf(outfp, " 0.99999999 0.00000001\n");
    fprintf(outfp, " 1\n");
    fprintf(outfp, " 0.00000000 0.00000000 1.00000000\n");

    if (xlinked)
        fprintf(outfp, " 0.00000000 0.50000000\n");
    fprintf(outfp, "2\n\n");

    if (nall == 1) {
        fprintf(outfp, "3 2\n");
        fprintf(outfp, " 0.90000000 0.10000000\n");
    }
    else {
        fprintf(outfp, "3 %d\n", nall);
        if (nall)
            fprintf(outfp, "%11.8f", freq[0]);
        for (i = 1; i < nall; i++)
            fprintf(outfp, "%11.8f", freq[i]);
        fprintf(outfp, "\n");
    }

    fprintf(outfp, "\n0 0\n");
    fprintf(outfp, " 0.00000000\n");
    fprintf(outfp, "1 0.10000000 0.09000000\n");
    fclose(outfp);

    free(freq);
}
 
void *allocMem (size_t nbytes)
{
    void *ptr;
    ptr = (void *) malloc(nbytes);
    if (!ptr) {
        printf("Not enough memory.\n");
        fatal_error();
    }
    return ptr;
}

void fatal_error (void)
{
    if (ibdfp) {
        fclose(ibdfp);
        unlink(ibdfile);
    }

    exit(1);
}
