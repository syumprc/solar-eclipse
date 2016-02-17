/* This file contains old nuclear family update routines */
/* used in the MLINK program */

#include "commondefs.h"
#include "mldefs.h"

Local Void getapprox(LINK)
struct LOC_seg *LINK;
{
  long first;
  double maxval;
  thisarray *WITH;

  maxval = (*LINK->p)->gen->genarray[0];
  WITH = (*LINK->p)->gen;
  for (first = 0; first < fgeno; first++) {
    if (WITH->genarray[first] > maxval)
      maxval = WITH->genarray[first];
  }
  WITH = (*LINK->p)->gen;
  for (first = 0; first < fgeno; first++) {
    approxarray[LINK->LINK->thisped - 1][first] =
      WITH->genarray[first] > (maxval *epsilon);
  }
  if (lasttime)
    return;
  WITH = (*LINK->p)->gen;
  for (first = 0; first < fgeno; first++) {
    if (!(approxarray[LINK->LINK->thisped - 1][first]))
      WITH->genarray[first] = 0.0;
  }
}  /*getapprox*/

#include "comnuclear.c"
