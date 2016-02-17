#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#define MAXID   30000
#define MAXALL  99
#define MAXGEN  MAXALL*(MAXALL+1)/2

struct _ego {
    short id;
    short fa;
    short mo;
    short fam;
    short gener;
    char sex[2];
    char twin[4];
    short twin1;
    char geno[7];
    char igeno[7];
    short lgen;
    short npgen;
    short tnpgen;
    unsigned char posgen[MAXGEN];
    unsigned char tposgen[MAXGEN];
    unsigned char savegen[MAXGEN];
} *ego;

short locid[MAXID+1];

struct _fam {
    short fa;
    short mo;
    short nkid;
    short kid1;
} *fam;

int maxrisk;
short *infam;
int stack;
short *istack, *jstack;
double *lstack;

int ntot, nfound, nfam;
int nall, ngen, nuntyped;
short *untyped, *imputed;
double freq[MAXALL], *gfreq;
short *iall, *jall;
char *genotyp[MAXGEN];
char *solarbin;

#define max(a,b)        ( (a) > (b) ? (a) : (b) )
#define min(a,b)        ( (a) < (b) ? (a) : (b) )

#define ign(a,b)        max(a,b)*(max(a,b) - 1)/2 + min(a,b) - 1

#define typed(a,b)      typed[a*(a + 1)/2 + b]
#define ibd(a,b)        ibd[a*(a + 1)/2 + b]
#define cumibd(a,b)     cumibd[a*(a + 1)/2 + b]

#define alpha(a,b,c,d)  alpha[(((a*ntot) + b)*2 + c)*ntot + d]

void set_posgen (void);
double calc_risk (int, int);
double sumlog (double, double);
void getibd (double*, double*);
void getalpha (double*, short**, short*, short*, short, short, short, short);
void wgt (double*, short, short, short, short, short, short);
void wgt2 (double*, short, short, short, short, short, short);

void *allocMem (size_t);

main (int argc, char **argv)
{
    int i, j, k, nimp, ndx;
    int imp, tries, ok, done;
    int nlast, lastgen;
    short *typed;
    double *alpha, *cumibd, *ibd;
    double *risk, *cumrskl, *cumrsku;

    if (argc != 3) {
        printf("usage: %s #imputations maxrisk(y/n)\n", argv[0]);
        exit(1);
    }

    if (sscanf(argv[1], "%d", &nimp) != 1 || nimp < 1) {
        printf("Number of imputations must be a positive integer.\n");
        exit(1);
    }

    switch (argv[2][0]) {
    case 'y':
        maxrisk = 1;
        break;
    case 'n':
        maxrisk = 0;
        break;
    default:
        printf("Enter y or n for maxrisk.\n");
        exit(1);
    }

    solarbin =  getenv("SOLAR_BIN");
    if (!solarbin) {
        printf("Environment variable SOLAR_BIN not defined.\n");
        exit(1);
    }

    srand48(time(NULL));

    if (!get_freqs())
        exit(1);

    if (!read_ped())
        exit(1);

    set_posgen();

    if (nuntyped == ntot) {
        FILE *fp = fopen("ibd.mat", "w");
        if (!fp) {
            printf("Cannot open ibd.mat\n");
            exit(1);
        }

        fprintf(fp, "%d %f\n", 1, 1.);
        for (i = 0; i < ntot; i++)
            fprintf(fp, "%5d %5d %10.7f %10.7f %1d\n", ego[i].id, ego[i].id,
                    -1., -1., 0);

        fclose(fp);
        exit(0);
    }

    infam = (short *) allocMem(ntot * sizeof(short));
    istack = (short *) allocMem(ntot * sizeof(short));
    jstack = (short *) allocMem(ntot * sizeof(short));
    lstack = (double *) allocMem(ntot * sizeof(double));

    typed = (short *) allocMem(ntot*(ntot + 1)/2 * sizeof(short));
	for (i = 0; i < ntot; i++) {
	    for (j = 0; j <= i; j++) {
            if (!strcmp(ego[i].geno, "      ")
                || !strcmp(ego[j].geno, "      "))
            {
                typed(i,j) = (short) 0;
            }
            else
                typed(i,j) = (short) 1;
        }
    }

    alpha = (double *) allocMem(4*ntot*ntot * sizeof(double));
    ibd = (double *) allocMem(ntot*(ntot + 1)/2 * sizeof(double));
    cumibd = (double *) allocMem(ntot*(ntot + 1)/2 * sizeof(double));
    for (i = 0; i < ntot*(ntot + 1)/2; i++)
        cumibd[i] = 0;

    risk = (double *) allocMem(ngen * sizeof(double));
    cumrskl = (double *) allocMem(ngen * sizeof(double));
    cumrsku = (double *) allocMem(ngen * sizeof(double));

    for (imp = 0; imp < nimp; imp++)
    {
        set_posgen();

        for (j = 0; j < nuntyped; j++)
            imputed[j] = 0;
        i = nuntyped;
        if (!check_genotypes()) {
            printf("impossible initial genotype vector, bailing out\n");
            exit(1);
        }

        tries = 0;
        nlast = 0;
        while (i > 0) {
            if (!nlast) {
                lastgen = 0;
                for (k = 0; k < nuntyped; k++) {
                    if (!imputed[k] && ego[untyped[k]].gener > lastgen)
                        lastgen = ego[untyped[k]].gener;
                }
                for (k = 0; k < nuntyped; k++) {
                    if (!imputed[k] && ego[untyped[k]].gener == lastgen)
                        nlast++;
                }
            }
            j = drand48() * nlast;
            if (j == nlast) j--;
            for (k = 0; k < nuntyped; k++) {
                if (!imputed[k] && ego[untyped[k]].gener == lastgen) {
                    if (!j) {
                        ndx = untyped[k];
                        ok = impute_genotype(ndx, risk, cumrskl, cumrsku, !imp);
                        imputed[k] = 1;
                        break;
                    }
                    j--;
                }
            }
            if (!ok || !(ok = check_genotypes())) {
                strcpy(ego[ndx].igeno, ego[ndx].geno);
                imputed[k] = 0;
                tries++;
                if (tries == 10) break;
            }
            else {
                tries = 0;
                nlast--;
                i--;
            }
        }

        if (nuntyped && !ok) {
            imp--;
            continue;
        }

        if (!make_ibds(alpha, cumibd, ibd))
            exit(1);

        if (!write_ibds(imp + 1, cumibd, typed))
            exit(1);

        if (!nuntyped)
            exit(0);
    }

    if (!write_ibds(nimp, cumibd, typed))
        exit(1);

    exit(0);
}

int write_ibds (int nimp, double *cumibd, short *typed)
{
    int i, j;
    FILE *fp = fopen("ibd.mat", "w");
    if (!fp) {
        printf("Cannot open ibd.mat\n");
        return 0;
    }

    fprintf(fp, "%d\n", nimp);
    for (i = 0; i < ntot; i++) {
        for (j = 0; j <= i; j++) {
            if (cumibd(i,j) != 0)
                fprintf(fp, "%5d %5d %10.7f %10.7f %1d\n", ego[i].id, ego[j].id,
                        cumibd(i,j) / nimp, 0., typed(i,j));
        }
    }

    fclose(fp);
    return 1;
}

int get_freqs (void)
{
    int i, j, k;
    double sum;
    char rec[1024];
    FILE *fp;

    fp = fopen("ibd.loc", "r");
    if (!fp) {
        printf("Cannot open ibd.loc\n");
        return 0;
    }

    if (!fgets(rec, sizeof(rec), fp)) {
        printf("Error reading ibd.loc, line 1\n");
        return 0;
    }

    if (sscanf(&rec[16], "%2d", &nall) != 1) {
        printf("Illegal #alleles, line 1 of ibd.loc\n");
        return 0;
    }

    if (nall > MAXALL) {
        printf("Too many alleles, MAXALL = %d\n", MAXALL);
        return 0;
    }

    sum = 0;
    for (i = 0; i < nall; i++) {
        if (!fgets(rec, sizeof(rec), fp)) {
            printf("Error reading ibd.loc, line %d\n", i+2);
            return 0;
        }
        if (sscanf(&rec[8], "%lf", &freq[i]) != 1) {
            printf("Illegal allele frequency, line %d of ibd.loc\n", i+2);
            return 0;
        }
        sum += freq[i];
    }

    fclose(fp);
    if (sum < 0.99999 || sum > 1.00001) {
        printf("Allele frequencies in ibd.loc do not sum to 1.\n");
        return 0;
    }

    ngen = nall*(nall + 1)/2;
    gfreq = allocMem(ngen * sizeof(double));
    iall = allocMem(ngen * sizeof(short));
    jall = allocMem(ngen * sizeof(short));
    for (i = 0; i < ngen; i++)
        genotyp[i] = (char *) allocMem(7);

    for (i = 1; i <= nall; i++) {
        for (j = 1; j < i; j++) {
            k = i*(i - 1)/2 + j - 1;
            iall[k] = i;
            jall[k] = j;
            sprintf(genotyp[k], "%3d%3d", j, i);
            genotyp[k][6] = '\0';
            gfreq[k] = 2*freq[i-1]*freq[j-1];
        }
        k = i*(i + 1)/2 - 1;
        iall[k] = i;
        jall[k] = i;
        sprintf(genotyp[k], "%3d%3d", i, i);
        genotyp[k][6] = '\0';
        gfreq[k] = freq[i-1]*freq[i-1];
    }

    return 1;
}

int read_ped (void)
{
    int i, j, ifa, imo, num;
    char rec[1024];
    char id[6], fa[6], mo[6], sex[2], twin[4], geno[7];
    FILE *fp;

    fp = fopen("ibd.in", "r");
    if (!fp) {
        printf("Cannot open ibd.in\n");
        return 0;
    }

    if (!fgets(rec, sizeof(rec), fp)) {
        printf("Error reading ibd.in, line 1\n");
        return 0;
    }

    if (!fgets(rec, sizeof(rec), fp)) {
        printf("Error reading ibd.in, line 2\n");
        return 0;
    }

    if (!fgets(rec, sizeof(rec), fp)) {
        printf("Error reading ibd.in, line 3\n");
        return 0;
    }

    if (sscanf(rec, "%d", &ntot) != 1 || ntot < 1) {
        printf("Illegal #individuals, line 3 of ibd.loc\n");
        return 0;
    }

    ego = (struct _ego *) allocMem(ntot * sizeof(struct _ego));
    locid[0] = -1;

    fam = (struct _fam *) allocMem(ntot * sizeof(struct _fam));
    nfam = 0;

    untyped = (short *) allocMem(ntot * sizeof(short));
    imputed = (short *) allocMem(ntot * sizeof(short));

    nfound = 0;
    for (i = 0; i < ntot; i++) {
        if (!fgets(rec, sizeof(rec), fp)) {
            printf("Error reading ibd.in, line %d\n", i+4);
            return 0;
        }

        if (sscanf(rec, "%5c%5c%5c%1c%3c%6c", id, fa, mo, sex, twin,
                   geno) != 6)
        {
            printf("Error reading ibd.in, line %d\n", i+4);
            return 0;
        }

        id[5] = '\0';
        num = atoi(id);
        if (num > MAXID) {
            printf("Ego's ID too large, line %d of ibd.in, MAXID = %d\n",
                   i+4, MAXID);
            return 0;
        }
        ego[i].id = num;

        fa[5] = '\0';
        if (strcmp(fa, "     ")) {
            num = atoi(fa);
            if (num > MAXID) {
                printf(
                "Father's ID too large, line %d of ibd.in, MAXID = %d\n",
                       i+4, MAXID);
                return 0;
            }
            ego[i].fa = num;
        }
        else
            ego[i].fa = 0;

        mo[5] = '\0';
        if (strcmp(mo, "     ")) {
            num = atoi(mo);
            if (num > MAXID) {
                printf(
                "Mother's ID too large, line %d of ibd.in, MAXID = %d\n",
                       i+4, MAXID);
                return 0;
            }
            ego[i].mo = num;
        }
        else
            ego[i].mo = 0;

        sex[1] = '\0';
        strcpy(ego[i].sex, sex);

        twin[3] = '\0';
        strcpy(ego[i].twin, twin);

        ego[i].twin1 = -1;
        if (strcmp(ego[i].twin, "   ")) {
            for (j = 0; j < i; j++) {
                if (!strcmp(ego[i].twin, ego[j].twin)) {
                    ego[i].twin1 = j;
                    break;
                }
            }
        }

        geno[6] = '\0';
        strcpy(ego[i].geno, geno);

        locid[ego[i].id] = i;

        if (!ego[i].fa) {
            ego[i].gener = 0;
            nfound++;
        }
        else {
            ifa = locid[ego[i].fa];
            imo = locid[ego[i].mo];
            ego[i].gener = max(ego[ifa].gener, ego[imo].gener) + 1;

            for (j = 0; j < nfam; j++) {
                if (fam[j].fa == ifa && fam[j].mo == imo)
                {
                    break;
                }
            }

            if (j == nfam) {
                fam[nfam].fa = ifa;
                fam[nfam].mo = imo;
                fam[nfam].nkid = 1;
                fam[nfam].kid1 = i;
                nfam++;
            }
            else {
                fam[j].nkid++;
            }
        }
    }

    fclose(fp);
    return 1;
}

void set_posgen (void)
{
    int i, j, ia, ja;

    nuntyped = 0;
    for (i = 0; i < ntot; i++) {
        if (strcmp(ego[i].geno, "      ")) {
            for (j = 0; j < ngen; j++)
                ego[i].posgen[j] = 0;
            sscanf(ego[i].geno, "%3d%3d", &ia, &ja);
            ego[i].posgen[ja*(ja-1)/2+ia-1] = 1;
            ego[i].npgen = 1;
        }
        else {
            untyped[nuntyped] = i;
            nuntyped++;
            for (j = 0; j < ngen; j++)
                ego[i].posgen[j] = 1;
            ego[i].npgen = ngen;
        }
        strcpy(ego[i].igeno, ego[i].geno);
        ego[i].tnpgen = ego[i].npgen;
        for (j = 0; j < ngen; j++)
            ego[i].tposgen[j] = ego[i].posgen[j];
    }
}

int check_genotypes (void)
{
    int i, j, k, ig;
    int ifa, imo, ifg, img, ikg;
    int ia, ja, iaf, jaf, iam, jam;
    int redo, ok;

    do
    {
        redo = 0;
        for (i = 0; i < nfam; i++)
        {
            ifa = fam[i].fa;
            imo = fam[i].mo;
            for (ig = 0; ig < ngen; ig++) {
                ego[ifa].savegen[ig] = 0;
                ego[imo].savegen[ig] = 0;
                for (j = 0; j < fam[i].nkid; j++) {
                    k = fam[i].kid1 + j;
                    ego[k].savegen[ig] = 0;
                }
            }

            for (ifg = 0; ifg < ngen; ifg++) {
                if (!ego[ifa].tposgen[ifg]) continue;
                for (img = 0; img < ngen; img++) {
                    if (!ego[imo].tposgen[img]) continue;
                    for (j = 0; j < fam[i].nkid; j++) {
                        ok = 0;
                        k = fam[i].kid1 + j;
                        for (ikg = 0; ikg < ngen; ikg++) {
                            if (!ego[k].tposgen[ikg]) continue;
                            if (ikg == ign(iall[ifg],iall[img]) ||
                                ikg == ign(iall[ifg],jall[img]) ||
                                ikg == ign(jall[ifg],iall[img]) ||
                                ikg == ign(jall[ifg],jall[img]))
                            {
                                ok = 1;
                            }
                        }
                        if (!ok)
                            break;
                    }

                    if (!ok)
                        continue;

                    ego[ifa].savegen[ifg] = 1;
                    ego[imo].savegen[img] = 1;
                    for (j = 0; j < fam[i].nkid; j++) {
                        k = fam[i].kid1 + j;
                        for (ikg = 0; ikg < ngen; ikg++) {
                            if (!ego[k].tposgen[ikg]) continue;
                            if (ikg == ign(iall[ifg],iall[img]) ||
                                ikg == ign(iall[ifg],jall[img]) ||
                                ikg == ign(jall[ifg],iall[img]) ||
                                ikg == ign(jall[ifg],jall[img]))
                            {
                                ego[k].savegen[ikg] = 1;
                            }
                        }
                    }
                }
            }

            for (ig = 0; ig < ngen; ig++) {
                if (!ego[ifa].savegen[ig]) {
                    if (ego[ifa].tposgen[ig]) {
                        ego[ifa].tposgen[ig] = 0;
                        ego[ifa].tnpgen--;
                        redo++;
                    }
                }

                if (!ego[imo].savegen[ig]) {
                    if (ego[imo].tposgen[ig]) {
                        ego[imo].tposgen[ig] = 0;
                        ego[imo].tnpgen--;
                        redo++;
                    }
                }

                for (j = 0; j < fam[i].nkid; j++) {
                    k = fam[i].kid1 + j;
                    if (!ego[k].savegen[ig]) {
                        if (ego[k].tposgen[ig]) {
                            ego[k].tposgen[ig] = 0;
                            ego[k].tnpgen--;
                            redo++;
                        }
                    }
                }
            }
        }

    } while (redo);

    for (i = 0; i < ntot; i++) {
        if (!ego[i].tnpgen)
            return 0;
    }

    for (i = 0; i < ntot; i++) {
        ego[i].npgen = ego[i].tnpgen;
        for (ig = 0; ig < ngen; ig++)
            ego[i].posgen[ig] = ego[i].tposgen[ig];
    }

    return 1;
}

int impute_genotype (int ndx, double *risk, double *cumrskl, double *cumrsku,
                     int first_imp)
{
    int i, j, k, imax, nmax;
    double zmxrisk, gprb, sumrisk;
    FILE *fp;

    for (i = 0; i < ntot; i++) {
        ego[i].tnpgen = ego[i].npgen;
        for (j = 0; j < ngen; j++)
            ego[i].tposgen[j] = ego[i].posgen[j];
    }

    if (strcmp(ego[ndx].twin, "   ")) {
        for (i = 0; i < ndx; i++) {
            if (!strcmp(ego[ndx].twin, ego[i].twin)) {
                for (j = 0; j < ngen; j++) {
                    ego[ndx].tposgen[j] = ego[i].posgen[j];
                    if (ego[i].posgen[j])
                        strcpy(ego[ndx].igeno, genotyp[j]);
                }
                ego[ndx].tnpgen = ego[i].npgen;
                return 1;
            }
        }
    }

    for (i = 0; i < ngen; i++)
    {
        risk[i] = -DBL_MAX;
        if (!gfreq[i])
            continue;

        if (!ego[ndx].posgen[i])
            continue;

        if (ego[ndx].npgen == 1) {
            risk[i] = 0;
            continue;
        }

        if (ego[ndx].npgen == ngen) {
            risk[i] = log(gfreq[i]);
            continue;
        }

        risk[i] = calc_risk(ndx, i);
    }

    sumrisk = -DBL_MAX;
    for (i = 0; i < ngen; i++) {
        if (risk[i] > -DBL_MAX) sumrisk = sumlog(sumrisk, risk[i]);
    }
    if (sumrisk == -DBL_MAX)
        return 0;

    for (i = 0; i < ngen; i++) {
        if (risk[i] > -DBL_MAX) risk[i] = exp(risk[i] - sumrisk);
        else risk[i] = 0;
    }

    cumrskl[0] = 0;
    cumrsku[0] = risk[0];
    zmxrisk = risk[0];

    imax = 0;
    for (j = 1; j < ngen; j++) {
        cumrskl[j] = cumrsku[j-1];
        cumrsku[j] = cumrsku[j-1] + risk[j];
        if (risk[j] > zmxrisk) {
            zmxrisk = risk[j];
            imax = j;
        }
    }

    nmax = 0;
    for (j = 0; j < ngen; j++) {
        if (risk[j] == zmxrisk) nmax++;
        ego[ndx].tposgen[j] = 0;
    }

    gprb = drand48();
    if ((maxrisk && first_imp) || zmxrisk >= 1) {
        if (nmax == 1) {
            strcpy(ego[ndx].igeno, genotyp[imax]);
            ego[ndx].tposgen[imax] = 1;
        }
        else {
            gprb *= nmax;
            k = 0;
            for (j = 0; j < ngen; j++) {
                if (risk[j] == zmxrisk) {
                    k++;
                    if (gprb >= k-1 && gprb < k) {
                        strcpy(ego[ndx].igeno, genotyp[j]);
                        ego[ndx].tposgen[j] = 1;
                    }
                }
            }
        }
    }
    else {
        for (j = 0; j < ngen; j++) {
            if (gprb >= cumrskl[j] && gprb < cumrsku[j]) {
                strcpy(ego[ndx].igeno, genotyp[j]);
                ego[ndx].tposgen[j] = 1;
            }
        }
    }

    ego[ndx].tnpgen = 1;

    return 1;
}

double calc_risk (int ndx, int gt)
{
    int i, j;
    int i0, j0;
    double lnlik, tlnlik;
    int ia, ja, iaf, jaf, iam, jam;
    double trnprb;

    int ftyped = 0, mtyped = 0;
    if (ego[ndx].fa && ego[locid[ego[ndx].fa]].npgen == 1)
        ftyped = 1;
    if (ego[ndx].mo && ego[locid[ego[ndx].mo]].npgen == 1)
        mtyped = 1;

    for (i = 0; i < ntot; i++)
        infam[i] = 0;

    for (i = 0; i < ntot; i++) {
        if (i == ndx)
        {
            infam[i] = 1;
        }
        if (ego[ndx].fa &&
            (ego[i].id == ego[ndx].fa || ego[i].id == ego[ndx].mo))
        {
            infam[i] = 1;
        }
        else if (ego[ndx].fa && (!ftyped || !mtyped) &&
                 ego[i].fa == ego[ndx].fa && ego[i].mo == ego[ndx].mo &&
                 ego[i].npgen == 1)
        {
            infam[i] = 1;
        }
        else if (ego[ndx].fa && ego[i].fa == ego[ndx].fa && !ftyped &&
                 ego[i].mo != ego[ndx].mo && ego[i].npgen == 1 &&
                 ego[locid[ego[i].mo]].npgen == 1)
        {
            infam[i] = 1;
            infam[locid[ego[i].mo]] = 1;
        }
        else if (ego[ndx].fa && ego[i].mo == ego[ndx].mo && !mtyped &&
                 ego[i].fa != ego[ndx].fa && ego[i].npgen == 1 &&
                 ego[locid[ego[i].fa]].npgen == 1)
        {
            infam[i] = 1;
            infam[locid[ego[i].fa]] = 1;
        }
        else if (ego[i].fa == ego[ndx].id && ego[i].npgen < ngen &&
                 ego[locid[ego[i].mo]].npgen == 1)
        {
            infam[i] = 1;
            infam[locid[ego[i].mo]] = 1;
        }
        else if (ego[i].mo == ego[ndx].id && ego[i].npgen < ngen &&
                 ego[locid[ego[i].fa]].npgen == 1)
        {
            infam[i] = 1;
            infam[locid[ego[i].fa]] = 1;
        }
    }

    tlnlik = -DBL_MAX;

    istack[0] = 0;
    jstack[0] = 0;
    lstack[0] = 0;
    stack = 1;

    do {
        stack--;
        i0 = istack[stack];
        j0 = jstack[stack];
        lnlik = lstack[stack];

        for (i = i0; i < ntot; i++) {
            if (!infam[i]) continue;
            if (i != ndx) {
                for (j = j0; j < ngen; j++) {
                    if (ego[i].posgen[j]) {
                        if (ego[i].npgen > 1) {
                            istack[stack] = i;
                            jstack[stack] = j + 1;
                            lstack[stack] = lnlik;
                            stack++;
                            if (stack==ntot) {
                                printf("stack overflow!!!\n");
                                exit(1);
                            }
                        }
                        break;
                    }
                }
                if (j == ngen) {
                    lnlik = -DBL_MAX;
                    break;
                }
            }
            else
                j = gt;

            ego[i].lgen = j;

            if (!ego[i].fa ||
                !infam[locid[ego[i].fa]] || !infam[locid[ego[i].mo]])
            {
                lnlik += log(gfreq[ego[i].lgen]);
            }
            else {
                ia = iall[ego[i].lgen];
                ja = jall[ego[i].lgen];
                iaf = iall[ego[locid[ego[i].fa]].lgen];
                jaf = jall[ego[locid[ego[i].fa]].lgen];
                iam = iall[ego[locid[ego[i].mo]].lgen];
                jam = jall[ego[locid[ego[i].mo]].lgen];
                trnprb = 0;
                if (ia == iaf && ja == iam || ia == iam && ja == iaf)
                    trnprb += .25;
                if (ia == iaf && ja == jam || ia == jam && ja == iaf)
                    trnprb += .25;
                if (ia == jaf && ja == iam || ia == iam && ja == jaf)
                    trnprb += .25;
                if (ia == jaf && ja == jam || ia == jam && ja == jaf)
                    trnprb += .25;
                if (trnprb)
                    lnlik += log(trnprb);
                else {
                    lnlik = -DBL_MAX;
                    break;
                }
            }

            j0 = 0;
        }

        tlnlik = sumlog(tlnlik, lnlik);

    } while (stack);

    return tlnlik;
}

int make_ibds (double *alpha, double *cumibd, double *ibd)
{
    int i, j, k, ii, jj;

    getibd(ibd, alpha);

    for (i = 0; i < ntot; i++) {

        /* if ind. i has an MZ twin, use the twin's IBDs */
        ii = i;
        if (ego[i].twin1 != -1) ii = ego[i].twin1;

        for (j = 0; j <= i; j++) {

            /* if ind. j has an MZ twin, use the twin's IBDs ... */
            jj = j;
            if (ego[j].twin1 != -1) jj = ego[j].twin1;

            /* ... unless ind. j is one of ind. i's parents, in which
               case these IBDs are correct and any previously computed
               IBDs involving ind. i and one of ind. j's MZ twins must
               be replaced */

            if (ego[ii].fa == ego[j].id || ego[ii].mo == ego[j].id) {
                if (ibd(max(ii,j),min(ii,j)))
                    cumibd(i,j) = cumibd(i,j) + ibd(max(ii,j),min(ii,j));
                if (jj != j) {
                    cumibd(i,jj) = cumibd(i,j);
                    ibd(i,jj) = ibd(i,j);
                    for (k = jj + 1; k < j; k++) {
                        if (ego[k].twin1 == jj) {
                            cumibd(i,k) = cumibd(i,j);
                            ibd(i,k) = ibd(i,j);
                        }
                    }
                }
            }
            else
                if (ibd(max(ii,jj),min(ii,jj)))
                    cumibd(i,j) = cumibd(i,j) + ibd(max(ii,jj),min(ii,jj));
        }
    }

    return 1;
}

double sumlog (double x, double y)
{
    double diff = x - y;
    if (diff < -32.236191302)
        return y;
    else if (diff < 32.236191302)
        return y + log(exp(diff) + 1);
    else
        return x;
}

void getibd (double *ibd, double *alpha)
{
    short i, j, ia, ja;
    short *gene[2], *fa, *mo;

    gene[0] = (short *) allocMem(ntot * sizeof(short));
    gene[1] = (short *) allocMem(ntot * sizeof(short));
    fa = (short *) allocMem(ntot * sizeof(short));
    mo = (short *) allocMem(ntot * sizeof(short));

    for (i = 0; i < ntot; i++) {
        sscanf(ego[i].igeno, "%3hd%3hd", &gene[0][i], &gene[1][i]);
        fa[i] = -1;
        mo[i] = -1;

        for (j = 0; j < ntot; j++) {
            alpha(0,i,0,j) = -1;
            alpha(0,i,1,j) = -1;
            alpha(1,i,0,j) = -1;
            alpha(1,i,1,j) = -1;
        }
    }

    for (i = nfound; i < ntot; i++) {
        fa[i] = locid[ego[i].fa];
        mo[i] = locid[ego[i].mo];
    }

    for (i = ntot - 1; i >= 0; i--) {
        for (j = i; j >= 0; j--) {
            for (ia = 0; ia < 2; ia++) {
                for (ja = 0; ja < 2; ja++) {
                    getalpha(alpha,gene,fa,mo,j,ja,i,ia);
                }
            }
        }
    }

    for (i = 0; i < ntot; i++) {
        for (j = 0; j <= i; j++) {
            ibd(i,j) = (alpha(0,i,0,j) + alpha(1,i,0,j) +
                        alpha(0,i,1,j) + alpha(1,i,1,j)) / 2;
        }
    }

    free(mo);
    free(fa);
    free(gene[1]);
    free(gene[0]);
}

void getalpha (double *alpha, short *gene[2], short *fa, short *mo, short i,
               short ia, short j, short ja)
{
    double w[4];

    if (i > j) {
        getalpha(alpha,gene,fa,mo,j,ja,i,ia);
        return;
    }

    if (alpha(ja,j,ia,i) != -1)
        return;

    if (fa[i] == -1 && mo[i] == -1 && fa[j] == -1 && mo[j] == -1) {
        if (i != j) {
            alpha(ja,j,ia,i) = alpha(ia,i,ja,j) = 0;
            return;
        }
        if (ia != ja) {
            alpha(ja,j,ia,i) = alpha(ia,i,ja,j) = 0;
            return;
        }
    }

    if (i == j) {
        if (ia == ja) {
            alpha(ja,j,ia,i) = alpha(ia,i,ja,j) = 1;
            return;
        }
        if (gene[ia][i] != gene[ja][j]) {
            alpha(ja,j,ia,i) = alpha(ia,i,ja,j) = 0;
            return;
        }

        getalpha(alpha,gene,fa,mo,mo[i],0,fa[i],0);
        getalpha(alpha,gene,fa,mo,mo[i],0,fa[i],1);
        getalpha(alpha,gene,fa,mo,mo[i],1,fa[i],0);
        getalpha(alpha,gene,fa,mo,mo[i],1,fa[i],1);

        wgt2(w, gene[0][i], gene[1][i], gene[0][fa[i]],
             gene[1][fa[i]], gene[0][mo[i]], gene[1][mo[i]]);

        alpha(1,i,0,i) =
            alpha(0,fa[i],0,mo[j]) * w[0] +
            alpha(1,fa[i],0,mo[j]) * w[1] +
            alpha(0,fa[i],1,mo[j]) * w[2] +
            alpha(1,fa[i],1,mo[j]) * w[3];
        alpha(0,i,1,i) = alpha(1,i,0,i);
        return;
    }

    if (gene[ia][i] != gene[ja][j]) {
        alpha(ja,j,ia,i) = alpha(ia,i,ja,j) = 0;
        return;
    }

    getalpha(alpha,gene,fa,mo,i,ia,fa[j],0);
    getalpha(alpha,gene,fa,mo,i,ia,fa[j],1);
    getalpha(alpha,gene,fa,mo,i,ia,mo[j],0);
    getalpha(alpha,gene,fa,mo,i,ia,mo[j],1);

    if (ja)
        wgt(w, gene[1][j], gene[0][j], gene[0][fa[j]],
            gene[1][fa[j]], gene[0][mo[j]], gene[1][mo[j]]);
    else
        wgt(w, gene[0][j], gene[1][j], gene[0][fa[j]],
            gene[1][fa[j]], gene[0][mo[j]], gene[1][mo[j]]);

    alpha(ja,j,ia,i) =
        alpha(0,fa[j],ia,i) * w[0] +
        alpha(1,fa[j],ia,i) * w[1] +
        alpha(0,mo[j],ia,i) * w[2] +
        alpha(1,mo[j],ia,i) * w[3];
    alpha(ia,i,ja,j) = alpha(ja,j,ia,i);
}

void wgt (double *w, short c1, short c2, short f1, short f2, short m1,
          short m2)
{
    int i;
    double sumw;

    for (i = 0; i < 4; i++)
        w[i] = 0;

    if (c1 == f1 && (c2 == m1 || c2 == m2)) w[0] = 1;
    if (c1 == f2 && (c2 == m1 || c2 == m2)) w[1] = 1;
    sumw = w[0] + w[1];
    if (sumw) {
        w[0] = w[0] / sumw;
        w[1] = w[1] / sumw;
    }

    if (c1 == m1 && (c2 == f1 || c2 == f2)) w[2] = 1;
    if (c1 == m2 && (c2 == f1 || c2 == f2)) w[3] = 1;
    sumw = w[2] + w[3];
    if (sumw) {
        w[2] = w[2] / sumw;
        w[3] = w[3] / sumw;
    }

    for (sumw = 0, i = 0; i < 4; i++)
        sumw += w[i];

    for (i = 0; i < 4; i++)
        w[i] = w[i] / sumw;
}

void wgt2 (double *w, short c1, short c2, short f1, short f2, short m1,
           short m2)
{
    int i;
    double sumw;

    for (i = 0; i < 4; i++)
        w[i] = 0;

    if (c1 == f1 && c2 == m1) w[0] = 1;
    if (c1 == f2 && c2 == m1) w[1] = 1;
    if (c1 == f1 && c2 == m2) w[2] = 1;
    if (c1 == f2 && c2 == m2) w[3] = 1;

    for (sumw = 0, i = 0; i < 4; i++)
        sumw += w[i];

    for (i = 0; i < 4; i++)
        w[i] = w[i] / sumw;
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
