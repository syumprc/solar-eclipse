/*
 *   This program implements an iterative version of the recursive algorithm
 *   due to Alun Thomas (Annals of Human Biology, 1988, 15(3):229-235).
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>

int currped, *nind;

#define NCOEFF 14	/* number of regression coeffs per class	*/
int NClass;

struct _class {
    int used;
    int ikin;
    int cnstr;
    char relation[1000];
    double phi2;
    double evar;
    double c[NCOEFF];
} *classes;

int *class_cnt;

int NGen;
int NDeg;
int NTwin;
int *fa, *mo, *twin;
int *first;

int MxNRel;		/* max num of ways a pair are related	*/

struct Rel {
    int sum;
    char *deg;		/* degree of relationship (plus 2)	*/
    char *gen;		/* number of removals			*/
    int *cnt;
} **relate, *classDef;

int readPedInfo (void);
int getClasses (void);
int getClassDefs (void);
void set_relate (void);
void relcount (char, char, int, int);
void write_relate (int, FILE*, FILE*);
int get_fraction (char *str, double *num);

int sum (int, int);
void set_sum (int, int, int);
int full (char, char, int, int);
void set_full (char, char, int, int, int);
int half (char, char, int, int);
void set_half (char, char, int, int, int);

void *allocMem (size_t);

int main (int argc, char **argv)
{
    int id, ifa, imo, isex, itwin, iped;
    int igen;
    int i, j, n, id0, nindt, lastped;
    char rec[1000];
    FILE *infp, *outfp, *unkfp;

    if (argc != 2) {
        printf("usage: relate mxnrel\n");
        exit(1);
    }

    if (sscanf(argv[1], "%d", &MxNRel) != 1 || MxNRel < 1) {
        printf("invalid argument: mxnrel = %s\n", argv[1]);
        exit(1);
    }

    /* MxNRel must be large enough to accomodate relate.tab
    if (MxNRel < 3) MxNRel = 3;	*/

    if (!readPedInfo())
        exit(1);

    if (!getClasses())
        exit(1);

    if (!getClassDefs())
        exit(1);

    outfp = fopen("mibdrel.ped", "w");
    if (!outfp) {
        printf("Cannot open mibdrel.ped\n");
        exit(1);
    }

    unkfp = fopen("relate.unk", "w");
    if (!unkfp) {
        printf("Cannot open relate.unk\n");
        exit(1);
    }

    fprintf(unkfp,
"This file contains an entry for each pair of individuals whose relationship\n\
is not currently defined in SOLAR. Note that the IDs shown are SOLAR indexed\n\
IDs (IBDIDs) rather than permanent IDs from your pedigree file.\n\n");

    infp = fopen("pedindex.out", "r");
    if (!infp) {
        printf("Cannot open pedindex.out\n");
        exit(1);
    }

    currped = 0;
    lastped = 0;
    nindt = 0;
    id0 = 0;
    n = 0;

    fprintf(outfp, "IBDID1,IBDID2,PHI2,CLASS\n");

    while (!lastped) {
        if (!fgets(rec, sizeof(rec), infp)) {
            lastped = 1;
            fclose(infp);
        }
        else if (sscanf(rec, "%d %d %d %d %d %d %d",
                        &id, &ifa, &imo, &isex, &itwin, &iped, &igen) != 7)
        {
            printf("Error reading pedindex.out, line %d\n", nindt + 1);
            exit(1);
        }
        else {
            if (!currped) {
                currped = iped;
                fa = (int *) allocMem((nind[currped-1]+1) * sizeof(int));
                mo = (int *) allocMem((nind[currped-1]+1) * sizeof(int));
                twin = (int *) allocMem((nind[currped-1]+1) * sizeof(int));
                relate = (struct Rel **)
                         allocMem(nind[currped-1]*sizeof(struct Rel *));
                for (i = 0; i < nind[currped-1]; i++) {
                    relate[i] = (struct Rel *)
                                allocMem(nind[currped-1]*sizeof(struct Rel));
                    for (j = 0; j < nind[currped-1]; j++) {
                        relate[i][j].deg = (char *) allocMem(MxNRel);
                        relate[i][j].gen = (char *) allocMem(MxNRel);
                        relate[i][j].cnt = (int *) allocMem(MxNRel*sizeof(int));
                    }
                }
            }
            nindt++;
            if (id != nindt) {
                printf("IDs must be sequential starting at 1.\n");
                exit(1);
            }
            if (id < ifa || id < imo) {
                printf("Parents must appear before offspring.\n");
                exit(1);
            }
        }

        if (!lastped && iped == currped) {
            n++;
            if (n > nind[currped-1]) {
                printf("More individuals than expected, pedigree %d\n",
                       currped);
                exit(1);
            }
            fa[n] = ifa > id0 ? ifa - id0 : 0;
            mo[n] = imo > id0 ? imo - id0 : 0;
            twin[n] = itwin;
        }
        else {
            set_relate();
            write_relate(id0, outfp, unkfp);

            if (!lastped) {
                if (n != nind[currped-1]) {
                    printf("Fewer individuals than expected, pedigree %d\n",
                           currped);
                    exit(1);
                }

                for (i = 0; i < nind[currped-1]; i++) {
                    for (j = 0; j < nind[currped-1]; j++) {
                        free(relate[i][j].deg);
                        free(relate[i][j].gen);
                        free(relate[i][j].cnt);
                    }
                    free(relate[i]);
                }
                free(relate);

                currped = iped;

                free(fa);
                fa = (int *) allocMem((nind[currped-1]+1) * sizeof(int));
                free(mo);
                mo = (int *) allocMem((nind[currped-1]+1) * sizeof(int));
                free(twin);
                twin = (int *) allocMem((nind[currped-1]+1) * sizeof(int));

                relate = (struct Rel **)
                         allocMem(nind[currped-1]*sizeof(struct Rel *));
                for (i = 0; i < nind[currped-1]; i++) {
                    relate[i] = (struct Rel *)
                                allocMem(nind[currped-1]*sizeof(struct Rel));
                    for (j = 0; j < nind[currped-1]; j++) {
                        relate[i][j].deg = (char *) allocMem(MxNRel);
                        relate[i][j].gen = (char *) allocMem(MxNRel);
                        relate[i][j].cnt = (int *) allocMem(MxNRel*sizeof(int));
                    }
                }

                id0 += n;
                n = 1;
                fa[n] = ifa > id0 ? ifa - id0 : 0;
                mo[n] = imo > id0 ? imo - id0 : 0;
                twin[n] = itwin;
            }
        }
    }

    fclose(unkfp);
    fclose(outfp);

    if (class_cnt[99])
        exit(1);

    unlink("relate.unk");
    exit(0);
}

void set_relate (void)
{
    int i, j, k;
    int id1, id2;
    int n = nind[currped-1];
    char deg, gen;

    fa[0] = 0;
    mo[0] = 0;

    for (i = 0; i < NTwin + 1; i++)
        first[i] = 0;
    twin[0] = 0;
    for (i = 1; i <= n; i++) {
        if (twin[i]) {
            if (first[twin[i]]) {
                fa[i] = 0;
                mo[i] = 0;
            }
            else
                first[twin[i]] = i;
        }
        if (twin[fa[i]])
            fa[i] = first[twin[fa[i]]];
        if (twin[mo[i]])
            mo[i] = first[twin[mo[i]]];
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            relate[i][j].sum = 0;
            for (k = 0; k < MxNRel; k++) {
                relate[i][j].deg[k] = 0;
                relate[i][j].gen[k] = 0;
                relate[i][j].cnt[k] = 0;
            }
        }
    }

    for (id1 = 1; id1 <= n; id1++) {
        for (id2 = id1; id2 <= n; id2++) {
            for (deg = 1; deg <= NGen + 1; deg++) {
                for (gen = 0; gen <= NGen; gen++) {
                    relcount(deg, gen, id1, id2);
                    if (id1 != id2 && gen > 0)
                        relcount(deg, gen, id2, id1);
                }
            }
        }
    }
}

void relcount (char deg, char gen, int id1, int id2)
{
    int fcnt, hcnt;

    if (gen == 0) {

        if (deg == 1) {
            if (id1 == id2) {
                set_full(1, 0, id1, id1, 1);
                set_sum(id1, id1, sum(id1, id1) + 1);
            }
        }

        else if (deg == 2) {
            if (id1 != id2 && fa[id1] && fa[id2]) {
                if (fa[id1] == fa[id2]) {
                    if (mo[id1] == mo[id2]) {
                        set_full(2, 0, id1, id2, 1);
                        set_full(2, 0, id2, id1, 1);
                        set_sum(id1, id2, sum(id1, id2) + 1);
                    }
                    else {
                        set_half(2, 0, id1, id2, 1);
                        set_half(2, 0, id2, id1, 1);
                        set_sum(id1, id2, sum(id1, id2) + 1);
                    }
                }
                else if (mo[id1] == mo[id2]) {
                    set_half(2, 0, id1, id2, 1);
                    set_half(2, 0, id2, id1, 1);
                    set_sum(id1, id2, sum(id1, id2) + 1);
                }
            }
        }

        else if (fa[id1] && fa[id2]) {
            fcnt = full(deg-1, 0, fa[id1], fa[id2])
                       + full(deg-1, 0, fa[id1], mo[id2])
                       + full(deg-1, 0, mo[id1], fa[id2])
                       + full(deg-1, 0, mo[id1], mo[id2]);
            if (fcnt) {
                set_full(deg, 0, id1, id2, fcnt);
                if (id1 != id2)
                    set_full(deg, 0, id2, id1, fcnt);
                set_sum(id1, id2, sum(id1, id2) + fcnt);
            }

            hcnt = half(deg-1, 0, fa[id1], fa[id2])
                       + half(deg-1, 0, fa[id1], mo[id2])
                       + half(deg-1, 0, mo[id1], fa[id2])
                       + half(deg-1, 0, mo[id1], mo[id2]);
            if (hcnt) {
                set_half(deg, 0, id1, id2, hcnt);
                if (id1 != id2)
                    set_half(deg, 0, id2, id1, hcnt);
                set_sum(id1, id2, sum(id1, id2) + hcnt);
            }
        }
    }

    else if (fa[id2]) {
        fcnt = full(deg, gen-1, id1, fa[id2])
                   + full(deg, gen-1, id1, mo[id2]);
        if (fcnt) {
            set_full(deg, gen, id1, id2, fcnt);
            set_sum(id1, id2, sum(id1, id2) + fcnt);
        }

        hcnt = half(deg, gen-1, id1, fa[id2])
                   + half(deg, gen-1, id1, mo[id2]);
        if (hcnt) {
            set_half(deg, gen, id1, id2, hcnt);
            set_sum(id1, id2, sum(id1, id2) + hcnt);
        }
    }
}

void write_relate (int id0, FILE *outfp, FILE *unkfp)
{
    int i, j, k, id1, id2;
    int class, done, icl, irel;
    int ifa, imo, psum;
    char deg, gen;
    int fcnt, hcnt;
    int **ufcnt, **uhcnt;
    double kin;
    struct Rel r;

    ufcnt = (int **) allocMem((NDeg+1)*sizeof(int *));
    uhcnt = (int **) allocMem((NDeg+1)*sizeof(int *));
    for (i = 0; i < NDeg+1; i++) {
        ufcnt[i] = (int *) allocMem((NGen+1)*sizeof(int));
        uhcnt[i] = (int *) allocMem((NGen+1)*sizeof(int));
    }

    for (i = 1; i <= nind[currped-1]; i++) {
        if (twin[i] && first[twin[i]] != i)
            id2 = first[twin[i]];
        else
            id2 = i;

        for (j = 1; j <= i; j++) {
            if (twin[i] && i != j && twin[i] == twin[j]) {
                class = 300;
                kin = 1;
            }
            else {
                if (twin[j] && first[twin[j]] != j)
                    id1 = first[twin[j]];
                else
                    id1 = j;

                class = 99;
                kin = 0;

                if (!sum(id1, id2)) {
                    class = 0;
                    kin = 0;
                }

                else if (id1 == id2) {
                    ifa = fa[id1];
                    imo = mo[id1];
                    if (!ifa || sum(ifa, imo) == 0) {
                        class = 1;
                        kin = 1;
                    }
                    else if (sum(ifa, imo) == 1) {
                        if (full(2, 1, ifa, imo)
                            + full(2, 1, imo, ifa) == 1)
                        {
                            class = 200;
                            kin = 1.125;
                        }
                        else if (full(3, 0, ifa, imo) == 1) {
                            class = 201;
                            kin = 1.0625;
                        }
                        else if (full(4, 0, ifa, imo) == 1) {
                            class = 202;
                            kin = 1.015625;
                        }
                        else if (full(2, 0, ifa, imo) == 1) {
                            class = 203;
                            kin = 1.25;
                        }
                        else if (full(1, 1, ifa, imo)
                                 + full(1, 1, imo, ifa) == 1)
                        {
                            class = 204;
                            kin = 1.25;
                        }
                        else if (half(2, 0, ifa, imo) == 1) {
                            class = 205;
                            kin = 1.125;
                        }
                        else {
                            class = 210;
                            kin = 0;
                        }
                    }
                    else {
                        class = 210;
                        kin = 0;
                    }
                }

                else {
                    done = 0;
                    for (icl = 0; !done && icl < NClass; icl++) {
                        if (classes[icl].used &&
                                classDef[icl].sum == sum(id1, id2))
                        {
                            irel = 0;
                            deg = classDef[icl].deg[irel];
                            done = 1;
                            while (done && deg) {
                                gen = classDef[icl].gen[irel];
                                if (deg > 0) {
                                    fcnt = full(deg, gen, id1, id2);
                                    if (gen > 0)
                                        fcnt += full(deg, gen, id2, id1);
                                    if (classDef[icl].cnt[irel] != fcnt)
                                        done = 0;
                                }
                                if (deg < 0) {
                                    hcnt = half(-deg, gen, id1, id2);
                                    if (gen > 0)
                                        hcnt += half(-deg, gen, id2, id1);
                                    if (classDef[icl].cnt[irel] != hcnt)
                                        done = 0;
                                }
                                irel++;
                                if (irel == MxNRel)
                                    break;
                                deg = classDef[icl].deg[irel];
                            }
                            if (done) {
                                class = icl;
                                kin = classes[icl].phi2;
                            }
                        }
                    }
                }
            }

            if (class == 99) {
                for (deg = 1; deg <= NDeg; deg++) {
                    for (gen = 0; gen <= NGen; gen++) {
                        ufcnt[deg][gen] = 0;
                        uhcnt[deg][gen] = 0;
                    }
                }
                fprintf(unkfp, "ibdid1=%d ibdid2=%d", j+id0, i+id0);
                r = relate[id1-1][id2-1];
                fprintf(unkfp, " %d", r.sum);
                k = 0;
                while (k < MxNRel && r.deg[k]) {
                    if (r.deg[k] > 0)
                        ufcnt[r.deg[k]][r.gen[k]] = r.cnt[k];
                    if (r.deg[k] < 0)
                        uhcnt[-r.deg[k]][r.gen[k]] = r.cnt[k];
                    k++;
                }
                r = relate[id2-1][id1-1];
                k = 0;
                while (k < MxNRel && r.deg[k]) {
                    if (r.deg[k] > 0 && r.gen[k] != 0)
                        ufcnt[r.deg[k]][r.gen[k]] += r.cnt[k];
                    if (r.deg[k] < 0 && r.gen[k] != 0)
                        uhcnt[-r.deg[k]][r.gen[k]] += r.cnt[k];
                    k++;
                }
                for (deg = 1; deg <= NDeg; deg++) {
                    for (gen = 0; gen <= NGen; gen++) {
                        if (ufcnt[deg][gen])
                            fprintf(unkfp, " %d %d %d",
                                    deg, gen, ufcnt[deg][gen]);
                        if (uhcnt[deg][gen])
                            fprintf(unkfp, " %d %d %d",
                                    -deg, gen, uhcnt[deg][gen]);
                    }
                }
                fprintf(unkfp, "\n");
            }

            fprintf(outfp, "%d,%d,%.15g,%d\n", i+id0, j+id0, kin, class);
            class_cnt[class]++;
        }
    }
}

int readPedInfo (void)
{
    int i, pedno, gen, twinid, junk;
    char rec[1000];
    FILE *infp;

    infp = fopen("pedindex.out", "r");
    if (!infp) {
        printf("Cannot open pedindex.out\n");
        return 0;
    }

    while (fgets(rec, sizeof(rec), infp)) {
        if (sscanf(rec, "%d %d %d %d %d %d", &junk, &junk, &junk, &junk, &junk, &pedno) != 6) {
            printf("Error reading pedigree.info\n");
            return 0;
        }
    }

    nind = (int *) allocMem(pedno*sizeof(int));
    for (i = 0; i < pedno; i++) {
        nind[i] = 0;
    }

    rewind(infp);
    NTwin = 0;
    NGen = 0;
    while (fgets(rec, sizeof(rec), infp)) {
        if (sscanf(rec, "%d %d %d %d %d %d %d", &junk, &junk, &junk, &junk, &twinid, &pedno, &gen) != 7) {
            printf("Error reading pedigree.info\n");
            return 0;
        }
        nind[pedno-1]++;
        if (gen > NGen)
            NGen = gen;
        if (twinid > NTwin)
            NTwin = twinid;
    }
    fclose(infp);

    first = (int *) allocMem((NTwin+1)*sizeof(int));
    NDeg = NGen + 1;

    return 1;
}

int getClasses (void)
{
    FILE *fp;
    int i, icl, icoeff, line;
    char *dir, fname[1000];
    char *recp, rec[10000];

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

    NClass = 0;
    line = 0;
    while (fgets(rec, sizeof(rec), fp)) {
        line++;
        if (!(recp = strtok(rec, "|")) || sscanf(recp, "%d", &icl) != 1 || icl < 0)
        {
            printf("Illegal class number on line %d of %s\n", line, fname);
            return 0;
        }
        if (icl > NClass)
            NClass = icl;
    }
    NClass++;

    class_cnt = (int *) allocMem(NClass*sizeof(int));
    classes = (struct _class *) allocMem(NClass*sizeof(struct _class));
    for (i = 0; i < NClass; i++) {
        class_cnt[i] = 0;
        classes[i].used = 0;
    }

    rewind(fp);
    line = 0;
    while (fgets(rec, sizeof(rec), fp)) {
        line++;
        if (!(recp = strtok(rec, "|")) || sscanf(recp, "%d", &icl) != 1 ||
            icl < 0 || icl >= NClass)
        {
            printf("Illegal class number on line %d of %s\n", line, fname);
            return 0;
        }
        classes[icl].used = 1;

        if (!(recp = strtok(NULL, "|")))
        {
            printf("Illegal class name on line %d of %s\n", line, fname);
            return 0;
        }
        strcpy(classes[icl].relation, recp);

        if (!(recp = strtok(NULL, "|")) ||
            sscanf(recp, "%d", &classes[icl].ikin) != 1 ||
            classes[icl].ikin < 0 || classes[icl].ikin >= NClass)
        {
            printf("Illegal kinship class on line %d of %s\n", line, fname);
            return 0;
        }

        if (!(recp = strtok(NULL, "|")) ||
            sscanf(recp, "%d", &classes[icl].cnstr) != 1 ||
            classes[icl].cnstr < -1 || classes[icl].cnstr > 1)
        {
            printf("Illegal constraint indicator on line %d of %s\n", line,
                   fname);
            return 0;
        }

        if (!(recp = strtok(NULL, "|")) ||
            !get_fraction(recp, &classes[icl].phi2))
        {
            printf("Illegal kinship coefficient on line %d of %s\n", line,
                   fname);
            return 0;
        }

        if (!(recp = strtok(NULL, "|")) ||
            !get_fraction(recp, &classes[icl].evar))
        {
            printf("Illegal expected variance on line %d of %s\n", line,
                   fname);
            return 0;
        }

        for (i = 0; i < NCOEFF; i++)
            classes[icl].c[i] = 0;

        icoeff = 0;
        while (recp = strtok(NULL, "|")) {
            if (icoeff == NCOEFF)
            {
                printf("More than %d coefficients on line %d of %s\n",
                       NCOEFF, line, fname);
                return 0;
            }
            if (!get_fraction(recp, &classes[icl].c[icoeff]))
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
    char *q, *p, fract[1024];
    double numer, denom;

    q = str;
    p = fract;
    while (*q && *q != '/') *p++ = *q++;
    *p = '\0';

    if (sscanf(fract, "%lf", &numer) != 1)
        return 0;

    if (!*q) {
        *num = numer;
        return 1;
    }

    q++;
    if (!*q || sscanf(q, "%lf", &denom) != 1 || denom == 0)
        return 0;

    *num = numer / denom;
    return 1;
}

int getClassDefs (void)
{
    FILE *fp;
    int i, j, icl, irel, line, num;
    char *dir, fname[1000], rec[10000];
    char *recp;

    dir = getenv("SOLAR_LIB");
    if (!dir) {
        printf("Environment variable SOLAR_LIB is not defined.\n");
        return 0;
    }

    sprintf(fname, "%s/relate.tab", dir);
    fp = fopen(fname, "r");
    if (!fp) {
        printf("Cannot open %s.\n", fname);
        return 0;
    }

    classDef = (struct Rel *) allocMem(NClass*sizeof(struct Rel));
    for (i = 0; i < NClass; i++) {
        classDef[i].deg = (char *) allocMem(MxNRel);
        classDef[i].gen = (char *) allocMem(MxNRel);
        classDef[i].cnt = (int *) allocMem(MxNRel*sizeof(int));
        for (j = 0; j < MxNRel; j++) {
            classDef[i].deg[j] = 0;
            classDef[i].gen[j] = 0;
            classDef[i].cnt[j] = 0;
        }
    }

    line = 0;
    while (fgets(rec, sizeof(rec), fp)) {
        line++;
        if (!(recp = strtok(rec, " \t\n")) ||
                sscanf(rec, "%d", &icl) != 1 || icl < 0 || icl >= NClass)
        {
            printf("Illegal class number on line %d of %s\n", line, fname);
            return 0;
        }

        if (!classes[icl].used)
        {
            printf(
            "Relative class %d appears in %s but not in %s/classes.tab.\n",
                   icl, fname, dir);
            return 0;
        }

        if (!(recp = strtok(NULL, " \t\n")) ||
            sscanf(recp, "%d", &num) != 1)
        {
            printf("Illegal sum on line %d of %s\n", line, fname);
            return 0;
        }
        classDef[icl].sum = num;

        irel = 0;
        while (recp = strtok(NULL, " \t\n")) {
            if (irel == MxNRel) {
                classes[icl].used = 0;
                break;;
            }
            if (sscanf(recp, "%d", &num) != 1)
            {
                printf("Degree #%d is illegal on line %d of %s\n", irel + 1,
                       line, fname);
                return 0;
            }
            classDef[icl].deg[irel] = num;

            if (!(recp = strtok(NULL, " \t\n")) ||
                sscanf(recp, "%d", &num) != 1 || num < 0)
            {
                printf("Generation #%d is illegal on line %d of %s\n",
                       irel + 1, line, fname);
                return 0;
            }
            classDef[icl].gen[irel] = num;

            if (!(recp = strtok(NULL, " \t\n")) ||
                sscanf(recp, "%d", &num) != 1 || num < 0 ||
                classDef[icl].sum > 0 && num > classDef[icl].sum)
            {
                printf("Count #%d is illegal on line %d of %s\n", irel + 1,
                       line, fname);
                return 0;
            }
            classDef[icl].cnt[irel] = num;

            irel++;
        }
    }

    fclose(fp);
    return 1;
}

int sum (int id1, int id2)
{
    struct Rel r;

    if (id2 >= id1)
        r = relate[id1-1][id2-1];
    else
        r = relate[id2-1][id1-1];

    return r.sum;
}

void set_sum (int id1, int id2, int sum)
{
    struct Rel *r;

    if (id2 >= id1)
        r = &relate[id1-1][id2-1];
    else
        r = &relate[id2-1][id1-1];

    r->sum = sum;
}

int full (char deg, char gen, int id1, int id2)
{
    int i;
    struct Rel r = relate[id1-1][id2-1];

    for (i = 0; i < MxNRel; i++)
        if (r.deg[i] == deg && r.gen[i] == gen)
            return r.cnt[i];

    return 0;
}

void set_full (char deg, char gen, int id1, int id2, int cnt)
{
    int i;
    struct Rel *r = &relate[id1-1][id2-1];

    for (i = 0; i < MxNRel && r->deg[i]; i++) ;
    if (i == MxNRel) {
        printf(\
"The maximum number of ways two individuals are related exceeds %d.\n\
Try running 'mibd relate' with argument 'mxnrel' set to a larger value.\n", MxNRel);
/*
printf("id1=%d id2=%d sum=%d deg=%d gen=%d cnt=%d\n",id1,id2,r->sum,deg,gen,cnt);
for(int j=0;j<MxNRel;j++)printf("j=%d deg=%d gen=%d cnt=%d\n",j,r->deg[j],r->gen[j],r->cnt[j]);
*/
        exit(1);
    }

    r->deg[i] = deg;
    r->gen[i] = gen;
    r->cnt[i] = cnt;
}

int half (char deg, char gen, int id1, int id2)
{
    int i;
    struct Rel r = relate[id1-1][id2-1];

    for (i = 0; i < MxNRel; i++)
        if (r.deg[i] == -deg && r.gen[i] == gen)
            return r.cnt[i];

    return 0;
}

void set_half (char deg, char gen, int id1, int id2, int cnt)
{
    int i;
    struct Rel *r = &relate[id1-1][id2-1];

    for (i = 0; i < MxNRel && r->deg[i]; i++) ;
    if (i == MxNRel) {
        printf(\
"The maximum number of ways two individuals are related exceeds %d.\n\
Try running 'mibd relate' with argument 'mxnrel' set to a larger value.\n", MxNRel);
/*
printf("id1=%d id2=%d sum=%d deg=%d gen=%d cnt=%d\n",id1,id2,r->sum,deg,gen,cnt);
for(int j=0;j<MxNRel;j++)printf("j=%d deg=%d gen=%d cnt=%d\n",j,r->deg[j],r->gen[j],r->cnt[j]);
*/
        exit(1);
    }

    r->deg[i] = -deg;
    r->gen[i] = gen;
    r->cnt[i] = cnt;
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
