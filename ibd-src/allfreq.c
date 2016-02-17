#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>

/* make sure that Int4 is a 4-byte integer! */
typedef int Int4;

void memory_error (void);

main (int argc, char **argv)
{
    Int4 LENC = 10000;
    Int4 LENI = 200000;
    Int4 LENL = 10;
    Int4 LENR = 200000;

    char showstat, getse;
    Int4 mxiter = 100;

    double *rarray;
    Int4 *iarray, *larray;
    char *carray;
    Int4 cneed, ineed, lneed, rneed;

    if (argc < 2 || argc > 3) {
        printf("usage: %s showStatus(y/n) [ computeSE(y/n) | maxIter ]\n",
               argv[0]);
        exit(1);
    }

    showstat = argv[1][0];
    switch (showstat) {
    case 'n':
    case 'y':
        break;
    default:
        printf("showStatus must be y or n.\n");
        exit(1);
    }

    if (argc == 3) {
        if (argv[2][0] == 'y')
            mxiter = -mxiter;
        else if (argv[2][0] != 'n' &&
                   (sscanf(argv[2], "%d", &mxiter) != 1 || mxiter < 1))
        {
            printf("maxIter must be a postive integer.\n");
            exit(1);
        }
    }

    rarray = (double *) malloc(LENR * sizeof(double));
    iarray = (Int4 *) malloc(LENI * sizeof(Int4));
    larray = (Int4 *) malloc(LENL * sizeof(Int4));
    carray = (char *) malloc(8 * LENC);

    if (!rarray || !iarray || !larray || !carray)
        memory_error();

    allfrq_(&showstat, &mxiter, carray, &LENC, iarray, &LENI, larray,
            &LENL, rarray, &LENR, &cneed, &ineed, &lneed, &rneed, 1);

    while (cneed || ineed || lneed || rneed) {
        if (cneed) {
            LENC = 1.0 * cneed;
            carray = (char *) realloc((void *) carray, 8 * LENC);
            if (!carray)
                memory_error();
        }

        if (ineed) {
            LENI = 1.0 * ineed;
            iarray = (Int4 *) realloc((void *) iarray, LENI * sizeof(Int4));
            if (!iarray)
                memory_error();
        }

        if (lneed) {
            LENL = 1.0 * lneed;
            larray = (Int4 *) realloc((void *) larray, LENL * sizeof(Int4));
            if (!larray)
                memory_error();
        }

        if (rneed) {
            LENR = 1.0 * rneed;
            rarray = (double *) realloc((void *) rarray, LENR * sizeof(double));
            if (!rarray)
                memory_error();
        }

        allfrq_(&showstat, &mxiter, carray, &LENC, iarray, &LENI, larray,
                &LENL, rarray, &LENR, &cneed, &ineed, &lneed, &rneed, 1);
    }

    exit(0);
}

void memory_error (void)
{
    unlink("allfreq.out");
    exit(1);
}
