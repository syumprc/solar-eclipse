#include <stdio.h>
#include "pipeback.h"

main ()
{
    unsigned int cksum;
    unsigned int noct;
    int status;

    status = getcrc ("any_list.c", &cksum, &noct);

    printf ("status is %d\n", status);
    printf ("checksum is %u\n", cksum);
    printf ("noct is %u\n", noct);
}

