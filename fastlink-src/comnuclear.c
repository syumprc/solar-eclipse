/* This file contains some of the old nuclear family update routines */
/* shared by modified versions of the LODSCORE, ILINK, LINKMAP, and
MLINK programs*/

Local double msegsex(LINK)
struct LOC_seg *LINK;
{
  int g1, g2, g3, g4, g5, g6, g7, g8, ms, mf, ms1, ms2, mf1, mf2, j, k, f1,
       f2, s1, s2;
  double val, temp2;
  double temp[maxchild];
  int FORLIM;
  thetavalues *WITH1;
  thisarray *WITH2;

  for (k = 0; k < nchild; k++)
    temp[k] = 0.0;
  if ((*LINK->p)->male) {
    mf = muthap[LINK->fseg - 1];
    LINK->secondseg = LINK->sseg;
    WITH1 = LINK->secondsex;
    FORLIM = LINK->send;
    for (j = LINK->sstart - 1; j < FORLIM; j++) {
      if (WITH1->segprob[j] == 0.0)
	LINK->secondseg++;
      else {
	temp2 = WITH1->segprob[j];
	s1 = invgenenum1[LINK->secondseg - 1];
	s2 = invgenenum2[LINK->secondseg - 1];
	ms1 = muthap[s1 - 1];
	ms2 = muthap[s2 - 1];
	if (s1 != s2) {
	  g1 = genenumber[LINK->fseg - 1][s1 - 1];
	  g2 = genenumber[LINK->fseg - 1][s2 - 1];
	  g3 = genenumber[LINK->fseg - 1][ms1 - 1];
	  g4 = genenumber[LINK->fseg - 1][ms2 - 1];
	  g5 = genenumber[mf - 1][s1 - 1];
	  g6 = genenumber[mf - 1][s2 - 1];
	  g7 = genenumber[mf - 1][ms1 - 1];
	  g8 = genenumber[mf - 1][ms2 - 1];
	  for (k = 0; k < nchild; k++) {
	    WITH2 = thischild[k];
	    if (malechild[k])
	      temp[k] += temp2 * ((1 - LINK->ps) * (WITH2->genarray[s1 - 1] +
				    WITH2->genarray[s2 - 1]) + LINK->ps *
		      (WITH2->genarray[ms1 - 1] + WITH2->genarray[ms2 - 1]));
	    else
	      temp[k] += temp2 * ((1 - LINK->pf) * (1 - LINK->ps) *
		      (WITH2->genarray[g1 - 1] + WITH2->genarray[g2 - 1]) +
		    (1 - LINK->pf) * LINK->ps * (WITH2->genarray[g3 - 1] +
			WITH2->genarray[g4 - 1]) + LINK->pf * (1 - LINK->ps) *
		      (WITH2->genarray[g5 - 1] + WITH2->genarray[g6 - 1]) +
		    LINK->pf * LINK->ps *
		      (WITH2->genarray[g7 - 1] + WITH2->genarray[g8 - 1]));
	  }
	} else {
	  g1 = genenumber[LINK->fseg - 1][s1 - 1];
	  g3 = genenumber[LINK->fseg - 1][ms1 - 1];
	  g5 = genenumber[mf - 1][s1 - 1];
	  g7 = genenumber[mf - 1][ms1 - 1];
	  for (k = 0; k < nchild; k++) {
	    WITH2 = thischild[k];
	    if (malechild[k])
	      temp[k] += 2.0 * temp2 * ((1 - LINK->ps) * WITH2->genarray[s1 - 1] +
					LINK->ps * WITH2->genarray[ms1 - 1]);
	    else
	      temp[k] += 2.0 * temp2 *
		  ((1 - LINK->pf) * (1 - LINK->ps) * WITH2->genarray[g1 - 1] +
		   LINK->ps * (1 - LINK->pf) * WITH2->genarray[g3 - 1] +
		   LINK->pf * (1 - LINK->ps) * WITH2->genarray[g5 - 1] +
		   LINK->pf * LINK->ps * WITH2->genarray[g7 - 1]);
	  }
	}
	LINK->secondseg++;
      }
    }
  } else {
    LINK->firstseg = LINK->fseg;
    ms = muthap[LINK->sseg - 1];
    WITH1 = LINK->firstsex;
    FORLIM = LINK->fend;
    for (j = LINK->fstart - 1; j < FORLIM; j++) {
      if (WITH1->segprob[j] == 0.0)
	LINK->firstseg++;
      else {
	temp2 = WITH1->segprob[j];
	f1 = invgenenum1[LINK->firstseg - 1];
	f2 = invgenenum2[LINK->firstseg - 1];
	mf1 = muthap[f1 - 1];
	mf2 = muthap[f2 - 1];
	if (f1 != f2) {
	  g1 = genenumber[LINK->sseg - 1][f1 - 1];
	  g2 = genenumber[LINK->sseg - 1][f2 - 1];
	  g3 = genenumber[LINK->sseg - 1][mf1 - 1];
	  g4 = genenumber[LINK->sseg - 1][mf2 - 1];
	  g5 = genenumber[ms - 1][f1 - 1];
	  g6 = genenumber[ms - 1][f2 - 1];
	  g7 = genenumber[ms - 1][mf1 - 1];
	  g8 = genenumber[ms - 1][mf2 - 1];
	  for (k = 0; k < nchild; k++) {
	    WITH2 = thischild[k];
	    if (malechild[k])
	      temp[k] += temp2 * ((1 - LINK->pf) * (WITH2->genarray[f1 - 1] +
				    WITH2->genarray[f2 - 1]) + LINK->pf *
		      (WITH2->genarray[mf1 - 1] + WITH2->genarray[mf2 - 1]));
	    else
	      temp[k] += temp2 * ((1 - LINK->pf) * (1 - LINK->ps) *
		      (WITH2->genarray[g1 - 1] + WITH2->genarray[g2 - 1]) +
		    (1 - LINK->ps) * LINK->pf * (WITH2->genarray[g3 - 1] +
			WITH2->genarray[g4 - 1]) + LINK->ps * (1 - LINK->pf) *
		      (WITH2->genarray[g5 - 1] + WITH2->genarray[g6 - 1]) +
		    LINK->pf * LINK->ps *
		      (WITH2->genarray[g7 - 1] + WITH2->genarray[g8 - 1]));
	  }
	} else {
	  g1 = genenumber[LINK->sseg - 1][f1 - 1];
	  g3 = genenumber[LINK->sseg - 1][mf1 - 1];
	  g5 = genenumber[ms - 1][f1 - 1];
	  g7 = genenumber[ms - 1][mf1 - 1];
	  for (k = 0; k < nchild; k++) {
	    WITH2 = thischild[k];
	    if (malechild[k])
	      temp[k] += 2.0 * temp2 * ((1 - LINK->pf) * WITH2->genarray[f1 - 1] +
					LINK->pf * WITH2->genarray[mf1 - 1]);
	    else
	      temp[k] += 2.0 * temp2 *
		  ((1 - LINK->pf) * (1 - LINK->ps) * WITH2->genarray[g1 - 1] +
		   LINK->pf * (1 - LINK->ps) * WITH2->genarray[g3 - 1] +
		   LINK->ps * (1 - LINK->pf) * WITH2->genarray[g5 - 1] +
		   LINK->ps * LINK->pf * WITH2->genarray[g7 - 1]);
	  }
	}
	LINK->firstseg++;
      }
    }
  }
  val = 1.0;
  for (k = 0; k < nchild; k++)
    val *= temp[k];
  return val;
}  /*msegsex*/


Local double msegsexf(LINK)
struct LOC_seg *LINK;
{
  int g1, g2, g3, g4, g5, g6, g7, g8, mf, ms1, ms2, j, k, l, s1, s2, slength;
  double val, temp1;
  double temp[maxchild][maxseg];
  double temp2[maxseg];
  int FORLIM;
  thetavalues *WITH1;
  thisarray *WITH2;

  mf = muthap[LINK->fseg - 1];
  slength = LINK->send - LINK->sstart + 1;
  for (k = 0; k < nchild; k++) {
    for (l = 0; l < slength; l++)
      temp[k][l] = 0.0;
  }
  LINK->secondseg = LINK->sseg;
  WITH1 = LINK->secondsex;
  FORLIM = LINK->send;
  for (j = LINK->sstart - 1; j < FORLIM; j++) {
    for (l = 0; l < slength; l++)
      temp2[l] = WITH1->segprob[j + l * slength];
    s1 = invgenenum1[LINK->secondseg - 1];
    s2 = invgenenum2[LINK->secondseg - 1];
    ms1 = muthap[s1 - 1];
    ms2 = muthap[s2 - 1];
    if (s1 != s2) {
      g1 = genenumber[LINK->fseg - 1][s1 - 1];
      g2 = genenumber[LINK->fseg - 1][s2 - 1];
      g3 = genenumber[LINK->fseg - 1][ms1 - 1];
      g4 = genenumber[LINK->fseg - 1][ms2 - 1];
      g5 = genenumber[mf - 1][s1 - 1];
      g6 = genenumber[mf - 1][s2 - 1];
      g7 = genenumber[mf - 1][ms1 - 1];
      g8 = genenumber[mf - 1][ms2 - 1];
      for (k = 0; k < nchild; k++) {
	WITH2 = thischild[k];
	if (malechild[k])
	  val = (1 - LINK->ps) * (WITH2->genarray[s1 - 1] + WITH2->genarray[s2 - 1]) +
		LINK->ps * (WITH2->genarray[ms1 - 1] + WITH2->genarray[ms2 - 1]);
	else
	  val = (1 - LINK->pf) * (1 - LINK->ps) * (WITH2->genarray[g1 - 1] +
		  WITH2->genarray[g2 - 1]) + (1 - LINK->pf) * LINK->ps *
		(WITH2->genarray[g3 - 1] + WITH2->genarray[g4 - 1]) +
	      LINK->pf * (1 - LINK->ps) * (WITH2->genarray[g5 - 1] + WITH2->
		    genarray[g6 - 1]) + LINK->pf * LINK->ps *
		(WITH2->genarray[g7 - 1] + WITH2->genarray[g8 - 1]);
	for (l = 0; l < slength; l++)
	  temp[k][l] += temp2[l] * val;
      }
    } else {
      g1 = genenumber[LINK->fseg - 1][s1 - 1];
      g3 = genenumber[LINK->fseg - 1][ms1 - 1];
      g5 = genenumber[mf - 1][s1 - 1];
      g7 = genenumber[mf - 1][ms1 - 1];
      for (k = 0; k < nchild; k++) {
	WITH2 = thischild[k];
	if (malechild[k])
	  val = 2.0 * ((1 - LINK->ps) * WITH2->genarray[s1 - 1] +
		       LINK->ps * WITH2->genarray[ms1 - 1]);
	else
	  val = 2.0 * ((1 - LINK->pf) * (1 - LINK->ps) * WITH2->genarray[g1 - 1] +
		       LINK->ps * (1 - LINK->pf) * WITH2->genarray[g3 - 1] +
		       LINK->pf * (1 - LINK->ps) * WITH2->genarray[g5 - 1] +
		       LINK->ps * LINK->pf * WITH2->genarray[g7 - 1]);
	for (l = 0; l < slength; l++)
	  temp[k][l] += temp2[l] * val;
      }
    }
    LINK->secondseg++;
  }
  temp1 = 0.0;
  for (l = 0; l < slength; l++) {
    val = 1.0;
    for (k = 0; k < nchild; k++)
      val *= temp[k][l];
    temp1 += val;
  }
  return temp1;
}  /*msegsexf*/


Local double segsex(LINK)
struct LOC_seg *LINK;
{
  int g1, g2, j, k, f1, f2, s1, s2;
  double val, temp2;
  double temp[maxchild];
  int FORLIM;
  thetavalues *WITH1;
  thisarray *WITH2;

  for (k = 0; k < nchild; k++)
    temp[k] = 0.0;
  if ((*LINK->p)->male) {
    LINK->secondseg = LINK->sseg;
    WITH1 = LINK->secondsex;
    FORLIM = LINK->send;
    for (j = LINK->sstart - 1; j < FORLIM; j++) {
      if (WITH1->segprob[j] == 0.0)
	LINK->secondseg++;
      else {
	temp2 = WITH1->segprob[j];
	s1 = invgenenum1[LINK->secondseg - 1];
	s2 = invgenenum2[LINK->secondseg - 1];
	if (s1 != s2) {
	  g1 = genenumber[LINK->fseg - 1][s1 - 1];
	  g2 = genenumber[LINK->fseg - 1][s2 - 1];
	  for (k = 0; k < nchild; k++) {
	    WITH2 = thischild[k];
	    if (malechild[k])
	      temp[k] += temp2 *
			 (WITH2->genarray[s1 - 1] + WITH2->genarray[s2 - 1]);
	    else
	      temp[k] += temp2 *
			 (WITH2->genarray[g1 - 1] + WITH2->genarray[g2 - 1]);
	  }
	} else {
	  g1 = genenumber[LINK->fseg - 1][s1 - 1];
	  for (k = 0; k < nchild; k++) {
	    WITH2 = thischild[k];
	    if (malechild[k])
	      temp[k] += 2.0 * temp2 * WITH2->genarray[s1 - 1];
	    else
	      temp[k] += 2.0 * temp2 * WITH2->genarray[g1 - 1];
	  }
	}
	LINK->secondseg++;
      }
    }
  } else {
    LINK->firstseg = LINK->fseg;
    WITH1 = LINK->firstsex;
    FORLIM = LINK->fend;
    for (j = LINK->fstart - 1; j < FORLIM; j++) {
      if (WITH1->segprob[j] == 0.0)
	LINK->firstseg++;
      else {
	temp2 = WITH1->segprob[j];
	f1 = invgenenum1[LINK->firstseg - 1];
	f2 = invgenenum2[LINK->firstseg - 1];
	if (f1 != f2) {
	  g1 = genenumber[LINK->sseg - 1][f1 - 1];
	  g2 = genenumber[LINK->sseg - 1][f2 - 1];
	  for (k = 0; k < nchild; k++) {
	    WITH2 = thischild[k];
	    if (malechild[k])
	      temp[k] += temp2 *
			 (WITH2->genarray[f1 - 1] + WITH2->genarray[f2 - 1]);
	    else
	      temp[k] += temp2 *
			 (WITH2->genarray[g1 - 1] + WITH2->genarray[g2 - 1]);
	  }
	} else {
	  g1 = genenumber[LINK->sseg - 1][f1 - 1];
	  for (k = 0; k < nchild; k++) {
	    WITH2 = thischild[k];
	    if (malechild[k])
	      temp[k] += 2.0 * temp2 * WITH2->genarray[f1 - 1];
	    else
	      temp[k] += 2.0 * temp2 * WITH2->genarray[g1 - 1];
	  }
	}
	LINK->firstseg++;
      }
    }
  }
  val = 1.0;
  for (k = 0; k < nchild; k++)
    val *= temp[k];
  return val;
}  /*segsex*/

Local double segsexf(LINK)
struct LOC_seg *LINK;
{
  int g1, g2, j, k, l, s1, s2, slength;
  double val, temp1;
  double temp[maxchild][maxseg];
  double temp2[maxseg];
  int FORLIM;
  thetavalues *WITH1;
  thisarray *WITH2;

  slength = LINK->send - LINK->sstart + 1;
  for (k = 0; k < nchild; k++) {
    for (l = 0; l < slength; l++)
      temp[k][l] = 0.0;
  }
  LINK->secondseg = LINK->sseg;
  WITH1 = LINK->secondsex;
  FORLIM = LINK->send;
  for (j = LINK->sstart - 1; j < FORLIM; j++) {
    for (l = 0; l < slength; l++)
      temp2[l] = WITH1->segprob[j + l * slength];
    s1 = invgenenum1[LINK->secondseg - 1];
    s2 = invgenenum2[LINK->secondseg - 1];
    if (s1 != s2) {
      g1 = genenumber[LINK->fseg - 1][s1 - 1];
      g2 = genenumber[LINK->fseg - 1][s2 - 1];
      for (k = 0; k < nchild; k++) {
	WITH2 = thischild[k];
	if (malechild[k])
	  val = WITH2->genarray[s1 - 1] + WITH2->genarray[s2 - 1];
	else
	  val = WITH2->genarray[g1 - 1] + WITH2->genarray[g2 - 1];
	if (val != 0.0) {
	  for (l = 0; l < slength; l++)
	    temp[k][l] += temp2[l] * val;
	}
      }
    } else {
      g1 = genenumber[LINK->fseg - 1][s1 - 1];
      for (k = 0; k < nchild; k++) {
	WITH2 = thischild[k];
	if (malechild[k])
	  val = 2.0 * WITH2->genarray[s1 - 1];
	else
	  val = 2.0 * WITH2->genarray[g1 - 1];
	if (val != 0.0) {
	  for (l = 0; l < slength; l++)
	    temp[k][l] += temp2[l] * val;
	}
      }
    }
    LINK->secondseg++;
  }
  temp1 = 0.0;
  for (l = 0; l < slength; l++) {
    val = 1.0;
    for (k = 0; k < nchild; k++)
      val *= temp[k][l];
    temp1 += val;
  }
  return temp1;
}

Local double segfun(LINK)
struct LOC_seg *LINK;
{
  int g1, g2, g3, g4, i, j, k, f1, f2, s1, s2;
  double val, temp1, temp2;
  double temp[maxchild];
  int FORLIM, FORLIM1;
  thetavalues *WITH1, *WITH2;
  thisarray *WITH3;

  LINK->firstseg = LINK->fseg;
  for (k = 0; k < nchild; k++)
    temp[k] = 0.0;
  WITH1 = LINK->firstsex;
  FORLIM = LINK->fend;
  for (i = LINK->fstart - 1; i < FORLIM; i++) {
    if (WITH1->segprob[i] == 0.0)
      LINK->firstseg++;
    else {
      temp1 = WITH1->segprob[i];
      f1 = invgenenum1[LINK->firstseg - 1];
      f2 = invgenenum2[LINK->firstseg - 1];
      LINK->secondseg = LINK->sseg;
      WITH2 = LINK->secondsex;
      if (f1 != f2) {
	FORLIM1 = LINK->send;
	for (j = LINK->sstart - 1; j < FORLIM1; j++) {
	  if (WITH2->segprob[j] == 0.0)
	    LINK->secondseg++;
	  else {
	    temp2 = temp1 * WITH2->segprob[j];
	    s1 = invgenenum1[LINK->secondseg - 1];
	    s2 = invgenenum2[LINK->secondseg - 1];
	    if (s1 != s2) {
	      g1 = genenumber[f1 - 1][s1 - 1];
	      g2 = genenumber[f1 - 1][s2 - 1];
	      g3 = genenumber[f2 - 1][s1 - 1];
	      g4 = genenumber[f2 - 1][s2 - 1];
	      for (k = 0; k < nchild; k++) {
		WITH3 = thischild[k];
		temp[k] +=
		  temp2 * (WITH3->genarray[g1 - 1] + WITH3->genarray[g2 - 1] +
			   WITH3->genarray[g3 - 1] + WITH3->genarray[g4 - 1]);
	      }
	    } else {
	      g1 = genenumber[f1 - 1][s1 - 1];
	      g3 = genenumber[f2 - 1][s1 - 1];
	      for (k = 0; k < nchild; k++) {
		WITH3 = thischild[k];
		temp[k] += 2 * temp2 *
			   (WITH3->genarray[g1 - 1] + WITH3->genarray[g3 - 1]);
	      }
	    }
	    LINK->secondseg++;
	  }
	}
      } else {
	FORLIM1 = LINK->send;
	for (j = LINK->sstart - 1; j < FORLIM1; j++) {
	  if (WITH2->segprob[j] == 0.0)
	    LINK->secondseg++;
	  else {
	    temp2 = temp1 * WITH2->segprob[j];
	    s1 = invgenenum1[LINK->secondseg - 1];
	    s2 = invgenenum2[LINK->secondseg - 1];
	    if (s1 != s2) {
	      g1 = genenumber[f1 - 1][s1 - 1];
	      g2 = genenumber[f1 - 1][s2 - 1];
	      for (k = 0; k < nchild; k++) {
		WITH3 = thischild[k];
		temp[k] += 2 * temp2 *
			   (WITH3->genarray[g1 - 1] + WITH3->genarray[g2 - 1]);
	      }
	    } else {
	      g1 = genenumber[f1 - 1][s1 - 1];
	      for (k = 0; k < nchild; k++) {
		WITH3 = thischild[k];
		temp[k] += 4 * temp2 * WITH3->genarray[g1 - 1];
	      }
	    }
	    LINK->secondseg++;
	  }
	}
      }
      LINK->firstseg++;
    }
  }
  val = 1.0;
  for (k = 0; k < nchild; k++)
    val *= temp[k];
  return val;
}

Local double msegfast(LINK)
struct LOC_seg *LINK;
{
  int g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, g13, g14, g15, g16,
       i, j, k, l, f1, f2, s1, s2, ms1, ms2, mf1, mf2, slength;
  double val, temp1;
  double temp[maxchild][maxseg];
  double temp2[maxseg];
  int FORLIM, FORLIM1;
  thetavalues *WITH1, *WITH2;
  thisarray *WITH3;

  LINK->firstseg = LINK->fseg;
  slength = LINK->send - LINK->sstart + 1;
  for (k = 0; k < nchild; k++) {
    for (l = 0; l < slength; l++)
      temp[k][l] = 0.0;
  }
  WITH1 = LINK->firstsex;
  FORLIM = LINK->fend;
  for (i = LINK->fstart - 1; i < FORLIM; i++) {
    if (WITH1->segprob[i] == 0.0)
      LINK->firstseg++;
    else {
      temp1 = WITH1->segprob[i];
      f1 = invgenenum1[LINK->firstseg - 1];
      f2 = invgenenum2[LINK->firstseg - 1];
      mf1 = muthap[f1 - 1];
      mf2 = muthap[f2 - 1];
      LINK->secondseg = LINK->sseg;
      WITH2 = LINK->secondsex;
      if (f1 != f2) {
	FORLIM1 = LINK->send;
	for (j = LINK->sstart - 1; j < FORLIM1; j++) {
	  for (l = 0; l < slength; l++)
	    temp2[l] = temp1 * WITH2->segprob[j + l * slength];
	  s1 = invgenenum1[LINK->secondseg - 1];
	  s2 = invgenenum2[LINK->secondseg - 1];
	  ms1 = muthap[s1 - 1];
	  ms2 = muthap[s2 - 1];
	  if (s1 != s2) {
	    g1 = genenumber[f1 - 1][s1 - 1];
	    g2 = genenumber[f1 - 1][s2 - 1];
	    g3 = genenumber[f2 - 1][s1 - 1];
	    g4 = genenumber[f2 - 1][s2 - 1];
	    g5 = genenumber[f1 - 1][ms1 - 1];
	    g6 = genenumber[f1 - 1][ms2 - 1];
	    g7 = genenumber[f2 - 1][ms1 - 1];
	    g8 = genenumber[f2 - 1][ms2 - 1];
	    g9 = genenumber[mf1 - 1][s1 - 1];
	    g10 = genenumber[mf1 - 1][s2 - 1];
	    g11 = genenumber[mf2 - 1][s1 - 1];
	    g12 = genenumber[mf2 - 1][s2 - 1];
	    g13 = genenumber[mf1 - 1][ms1 - 1];
	    g14 = genenumber[mf1 - 1][ms2 - 1];
	    g15 = genenumber[mf2 - 1][ms1 - 1];
	    g16 = genenumber[mf2 - 1][ms2 - 1];
	    for (k = 0; k < nchild; k++) {
	      WITH3 = thischild[k];
	      val = (1 - LINK->ps) * (1 - LINK->pf) *
		    (WITH3->genarray[g1 - 1] + WITH3->genarray[g2 - 1] +
		      WITH3->genarray[g3 - 1] + WITH3->genarray[g4 - 1]) +
		  LINK->ps * (1 - LINK->pf) * (WITH3->genarray[g5 - 1] +
		      WITH3->genarray[g6 - 1] + WITH3->genarray[g7 - 1] +
		      WITH3->genarray[g8 - 1]) + LINK->pf * (1 - LINK->ps) *
		    (WITH3->genarray[g9 - 1] + WITH3->genarray[g10 - 1] +
		      WITH3->genarray[g11 - 1] + WITH3->genarray[g12 - 1]) +
		  LINK->pf * LINK->ps *
		    (WITH3->genarray[g13 - 1] + WITH3->genarray[g14 - 1] +
		      WITH3->genarray[g15 - 1] + WITH3->genarray[g16 - 1]);
	      if (val != 0.0) {
		for (l = 0; l < slength; l++)
		  temp[k][l] += temp2[l] * val;
	      }
	    }
	  } else {
	    g1 = genenumber[f1 - 1][s1 - 1];
	    g3 = genenumber[f2 - 1][s1 - 1];
	    g5 = genenumber[f1 - 1][ms1 - 1];
	    g7 = genenumber[f2 - 1][ms1 - 1];
	    g9 = genenumber[mf1 - 1][s1 - 1];
	    g11 = genenumber[mf2 - 1][s1 - 1];
	    g13 = genenumber[mf1 - 1][ms1 - 1];
	    g15 = genenumber[mf2 - 1][ms1 - 1];
	    for (k = 0; k < nchild; k++) {
	      WITH3 = thischild[k];
	      val = 2 * ((1 - LINK->ps) * (1 - LINK->pf) *
		      (WITH3->genarray[g1 - 1] + WITH3->genarray[g3 - 1]) +
		    LINK->ps * (1 - LINK->pf) * (WITH3->genarray[g5 - 1] +
			WITH3->genarray[g7 - 1]) + LINK->pf * (1 - LINK->ps) *
		      (WITH3->genarray[g9 - 1] + WITH3->genarray[g11 - 1]) +
		    LINK->pf * LINK->ps *
		      (WITH3->genarray[g13 - 1] + WITH3->genarray[g15 - 1]));
	      for (l = 0; l < slength; l++)
		temp[k][l] += temp2[l] * val;
	    }
	  }
	  LINK->secondseg++;
	}
      } else {
	FORLIM1 = LINK->send;
	for (j = LINK->sstart - 1; j < FORLIM1; j++) {
	  for (l = 0; l < slength; l++)
	    temp2[l] = temp1 * WITH2->segprob[j + l * slength];
	  s1 = invgenenum1[LINK->secondseg - 1];
	  s2 = invgenenum2[LINK->secondseg - 1];
	  ms1 = muthap[s1 - 1];
	  ms2 = muthap[s2 - 1];
	  if (s1 != s2) {
	    g1 = genenumber[f1 - 1][s1 - 1];
	    g2 = genenumber[f1 - 1][s2 - 1];
	    g5 = genenumber[f1 - 1][ms1 - 1];
	    g6 = genenumber[f1 - 1][ms2 - 1];
	    g9 = genenumber[mf1 - 1][s1 - 1];
	    g10 = genenumber[mf1 - 1][s2 - 1];
	    g13 = genenumber[mf1 - 1][ms1 - 1];
	    g14 = genenumber[mf1 - 1][ms2 - 1];
	    for (k = 0; k < nchild; k++) {
	      WITH3 = thischild[k];
	      val = 2 * ((1 - LINK->ps) * (1 - LINK->pf) *
		      (WITH3->genarray[g1 - 1] + WITH3->genarray[g2 - 1]) +
		    LINK->ps * (1 - LINK->pf) * (WITH3->genarray[g5 - 1] +
			WITH3->genarray[g6 - 1]) + LINK->pf * (1 - LINK->ps) *
		      (WITH3->genarray[g9 - 1] + WITH3->genarray[g10 - 1]) +
		    LINK->pf * LINK->ps *
		      (WITH3->genarray[g13 - 1] + WITH3->genarray[g14 - 1]));
	      if (val != 0.0) {
		for (l = 0; l < slength; l++)
		  temp[k][l] += temp2[l] * val;
	      }
	    }
	  } else {
	    g1 = genenumber[f1 - 1][s1 - 1];
	    g5 = genenumber[f1 - 1][ms1 - 1];
	    g9 = genenumber[mf1 - 1][s1 - 1];
	    g13 = genenumber[mf1 - 1][ms1 - 1];
	    for (k = 0; k < nchild; k++) {
	      WITH3 = thischild[k];
	      val = 4 * ((1 - LINK->ps) * (1 - LINK->pf) * WITH3->genarray[g1 - 1] +
			 LINK->ps * (1 - LINK->pf) * WITH3->genarray[g5 - 1] +
			 LINK->pf * (1 - LINK->ps) * WITH3->genarray[g9 - 1] +
			 LINK->pf * LINK->ps * WITH3->genarray[g13 - 1]);
	      if (val != 0.0) {
		for (l = 0; l < slength; l++)
		  temp[k][l] += temp2[l] * val;
	      }
	    }
	  }
	  LINK->secondseg++;
	}
      }
      LINK->firstseg++;
    }
  }
  temp1 = 0.0;
  for (l = 0; l < slength; l++) {
    val = 1.0;
    for (k = 0; k < nchild; k++)
      val *= temp[k][l];
    temp1 += val;
  }
  return temp1;
}  /*msegfast*/


Local double msegfun(LINK)
struct LOC_seg *LINK;
{
  int g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, g13, g14, g15, g16,
       i, j, k, f1, f2, s1, s2, ms1, ms2, mf1, mf2;
  double val, temp1, temp2;
  double temp[maxchild];
  int FORLIM, FORLIM1;
  thetavalues *WITH1, *WITH2;
  thisarray *WITH3;

  LINK->firstseg = LINK->fseg;
  for (k = 0; k < nchild; k++)
    temp[k] = 0.0;
  WITH1 = LINK->firstsex;
  FORLIM = LINK->fend;
  for (i = LINK->fstart - 1; i < FORLIM; i++) {
    if (WITH1->segprob[i] == 0.0)
      LINK->firstseg++;
    else {
      temp1 = WITH1->segprob[i];
      f1 = invgenenum1[LINK->firstseg - 1];
      f2 = invgenenum2[LINK->firstseg - 1];
      mf1 = muthap[f1 - 1];
      mf2 = muthap[f2 - 1];
      LINK->secondseg = LINK->sseg;
      WITH2 = LINK->secondsex;
      if (f1 != f2) {
	FORLIM1 = LINK->send;
	for (j = LINK->sstart - 1; j < FORLIM1; j++) {
	  if (WITH2->segprob[j] == 0.0)
	    LINK->secondseg++;
	  else {
	    temp2 = temp1 * WITH2->segprob[j];
	    s1 = invgenenum1[LINK->secondseg - 1];
	    s2 = invgenenum2[LINK->secondseg - 1];
	    ms1 = muthap[s1 - 1];
	    ms2 = muthap[s2 - 1];
	    if (s1 != s2) {
	      g1 = genenumber[f1 - 1][s1 - 1];
	      g2 = genenumber[f1 - 1][s2 - 1];
	      g3 = genenumber[f2 - 1][s1 - 1];
	      g4 = genenumber[f2 - 1][s2 - 1];
	      g5 = genenumber[f1 - 1][ms1 - 1];
	      g6 = genenumber[f1 - 1][ms2 - 1];
	      g7 = genenumber[f2 - 1][ms1 - 1];
	      g8 = genenumber[f2 - 1][ms2 - 1];
	      g9 = genenumber[mf1 - 1][s1 - 1];
	      g10 = genenumber[mf1 - 1][s2 - 1];
	      g11 = genenumber[mf2 - 1][s1 - 1];
	      g12 = genenumber[mf2 - 1][s2 - 1];
	      g13 = genenumber[mf1 - 1][ms1 - 1];
	      g14 = genenumber[mf1 - 1][ms2 - 1];
	      g15 = genenumber[mf2 - 1][ms1 - 1];
	      g16 = genenumber[mf2 - 1][ms2 - 1];
	      for (k = 0; k < nchild; k++) {
		WITH3 = thischild[k];
		temp[k] += temp2 * ((1 - LINK->ps) * (1 - LINK->pf) *
			(WITH3->genarray[g1 - 1] + WITH3->genarray[g2 - 1] +
			  WITH3->genarray[g3 - 1] + WITH3->genarray[g4 - 1]) +
		      LINK->ps * (1 - LINK->pf) * (WITH3->genarray[g5 - 1] +
			  WITH3->genarray[g6 - 1] + WITH3->genarray[g7 - 1] +
			  WITH3->genarray[g8 - 1]) + LINK->pf *
			(1 - LINK->ps) * (WITH3->genarray[g9 - 1] + WITH3->
			    genarray[g10 - 1] + WITH3->genarray[g11 - 1] +
			  WITH3->genarray[g12 - 1]) + LINK->pf * LINK->ps *
			(WITH3->genarray[g13 - 1] + WITH3->genarray[g14 - 1] +
			  WITH3->genarray[g15 - 1] +
			  WITH3->genarray[g16 - 1]));
	      }
	    } else {
	      g1 = genenumber[f1 - 1][s1 - 1];
	      g3 = genenumber[f2 - 1][s1 - 1];
	      g5 = genenumber[f1 - 1][ms1 - 1];
	      g7 = genenumber[f2 - 1][ms1 - 1];
	      g9 = genenumber[mf1 - 1][s1 - 1];
	      g11 = genenumber[mf2 - 1][s1 - 1];
	      g13 = genenumber[mf1 - 1][ms1 - 1];
	      g15 = genenumber[mf2 - 1][ms1 - 1];
	      for (k = 0; k < nchild; k++) {
		WITH3 = thischild[k];
		temp[k] += 2 * temp2 * ((1 - LINK->ps) * (1 - LINK->pf) *
			(WITH3->genarray[g1 - 1] + WITH3->genarray[g3 - 1]) +
		      LINK->ps * (1 - LINK->pf) *
			(WITH3->genarray[g5 - 1] + WITH3->genarray[g7 - 1]) +
		      LINK->pf * (1 - LINK->ps) *
			(WITH3->genarray[g9 - 1] + WITH3->genarray[g11 - 1]) +
		      LINK->pf * LINK->ps * (WITH3->genarray[g13 - 1] +
			  WITH3->genarray[g15 - 1]));
	      }
	    }
	    LINK->secondseg++;
	  }
	}
      } else {
	FORLIM1 = LINK->send;
	for (j = LINK->sstart - 1; j < FORLIM1; j++) {
	  if (WITH2->segprob[j] == 0.0)
	    LINK->secondseg++;
	  else {
	    temp2 = temp1 * WITH2->segprob[j];
	    s1 = invgenenum1[LINK->secondseg - 1];
	    s2 = invgenenum2[LINK->secondseg - 1];
	    ms1 = muthap[s1 - 1];
	    ms2 = muthap[s2 - 1];
	    if (s1 != s2) {
	      g1 = genenumber[f1 - 1][s1 - 1];
	      g2 = genenumber[f1 - 1][s2 - 1];
	      g5 = genenumber[f1 - 1][ms1 - 1];
	      g6 = genenumber[f1 - 1][ms2 - 1];
	      g9 = genenumber[mf1 - 1][s1 - 1];
	      g10 = genenumber[mf1 - 1][s2 - 1];
	      g13 = genenumber[mf1 - 1][ms1 - 1];
	      g14 = genenumber[mf1 - 1][ms2 - 1];
	      for (k = 0; k < nchild; k++) {
		WITH3 = thischild[k];
		temp[k] += 2 * temp2 * ((1 - LINK->ps) * (1 - LINK->pf) *
			(WITH3->genarray[g1 - 1] + WITH3->genarray[g2 - 1]) +
		      LINK->ps * (1 - LINK->pf) *
			(WITH3->genarray[g5 - 1] + WITH3->genarray[g6 - 1]) +
		      LINK->pf * (1 - LINK->ps) *
			(WITH3->genarray[g9 - 1] + WITH3->genarray[g10 - 1]) +
		      LINK->pf * LINK->ps * (WITH3->genarray[g13 - 1] +
			  WITH3->genarray[g14 - 1]));
	      }
	    } else {
	      g1 = genenumber[f1 - 1][s1 - 1];
	      g5 = genenumber[f1 - 1][ms1 - 1];
	      g9 = genenumber[mf1 - 1][s1 - 1];
	      g13 = genenumber[mf1 - 1][ms1 - 1];
	      for (k = 0; k < nchild; k++) {
		WITH3 = thischild[k];
		temp[k] += 4 * temp2 * ((1 - LINK->ps) * (1 - LINK->pf) *
			WITH3->genarray[g1 - 1] +
		      LINK->ps * (1 - LINK->pf) * WITH3->genarray[g5 - 1] +
		      LINK->pf * (1 - LINK->ps) * WITH3->genarray[g9 - 1] +
		      LINK->pf * LINK->ps * WITH3->genarray[g13 - 1]);
	      }
	    }
	    LINK->secondseg++;
	  }
	}
      }
      LINK->firstseg++;
    }
  }
  val = 1.0;
  for (k = 0; k < nchild; k++)
    val *= temp[k];
  return val;
}  /*msegfun*/

Local double segfast(LINK)
struct LOC_seg *LINK;
{
  int g1, g2, g3, g4, i, j, k, l, f1, f2, s1, s2, slength;
  double val, temp1;
  double temp[maxchild][maxseg];
  double temp2[maxseg];
  int FORLIM, FORLIM1;
  thetavalues *WITH1, *WITH2;
  thisarray *WITH3;

  LINK->firstseg = LINK->fseg;
  slength = LINK->send - LINK->sstart + 1;
  for (k = 0; k < nchild; k++) {
    for (l = 0; l < slength; l++)
      temp[k][l] = 0.0;
  }
  WITH1 = LINK->firstsex;
  FORLIM = LINK->fend;
  for (i = LINK->fstart - 1; i < FORLIM; i++) {
    if (WITH1->segprob[i] == 0.0)
      LINK->firstseg++;
    else {
      temp1 = WITH1->segprob[i];
      f1 = invgenenum1[LINK->firstseg - 1];
      f2 = invgenenum2[LINK->firstseg - 1];
      LINK->secondseg = LINK->sseg;
      WITH2 = LINK->secondsex;
      if (f1 != f2) {
	FORLIM1 = LINK->send;
	for (j = LINK->sstart - 1; j < FORLIM1; j++) {
	  for (l = 0; l < slength; l++)
	    temp2[l] = temp1 * WITH2->segprob[j + l * slength];
	  s1 = invgenenum1[LINK->secondseg - 1];
	  s2 = invgenenum2[LINK->secondseg - 1];
	  if (s1 != s2) {
	    g1 = genenumber[f1 - 1][s1 - 1];
	    g2 = genenumber[f1 - 1][s2 - 1];
	    g3 = genenumber[f2 - 1][s1 - 1];
	    g4 = genenumber[f2 - 1][s2 - 1];
	    for (k = 0; k < nchild; k++) {
	      WITH3 = thischild[k];
	      val = WITH3->genarray[g1 - 1] + WITH3->genarray[g2 - 1] +
		    WITH3->genarray[g3 - 1] + WITH3->genarray[g4 - 1];
	      if (val != 0.0) {
		for (l = 0; l < slength; l++)
		  temp[k][l] += temp2[l] * val;
	      }
	    }
	  } else {
	    g1 = genenumber[f1 - 1][s1 - 1];
	    g3 = genenumber[f2 - 1][s1 - 1];
	    for (k = 0; k < nchild; k++) {
	      WITH3 = thischild[k];
	      val = 2.0 * (WITH3->genarray[g1 - 1] + WITH3->genarray[g3 - 1]);
	      if (val != 0.0) {
		for (l = 0; l < slength; l++)
		  temp[k][l] += temp2[l] * val;
	      }
	    }
	  }
	  LINK->secondseg++;
	}
      } else {
	FORLIM1 = LINK->send;
	for (j = LINK->sstart - 1; j < FORLIM1; j++) {
	  for (l = 0; l < slength; l++)
	    temp2[l] = temp1 * WITH2->segprob[j + l * slength];
	  s1 = invgenenum1[LINK->secondseg - 1];
	  s2 = invgenenum2[LINK->secondseg - 1];
	  if (s1 != s2) {
	    g1 = genenumber[f1 - 1][s1 - 1];
	    g2 = genenumber[f1 - 1][s2 - 1];
	    for (k = 0; k < nchild; k++) {
	      WITH3 = thischild[k];
	      val = 2.0 * (WITH3->genarray[g1 - 1] + WITH3->genarray[g2 - 1]);
	      if (val != 0.0) {
		for (l = 0; l < slength; l++)
		  temp[k][l] += temp2[l] * val;
	      }
	    }
	  } else {
	    g1 = genenumber[f1 - 1][s1 - 1];
	    for (k = 0; k < nchild; k++) {
	      WITH3 = thischild[k];
	      val = 4.0 * WITH3->genarray[g1 - 1];
	      if (val != 0.0) {
		for (l = 0; l < slength; l++)
		  temp[k][l] += temp2[l] * val;
	      }
	    }
	  }
	  LINK->secondseg++;
	}
      }
      LINK->firstseg++;
    }
  }
  temp1 = 0.0;
  for (l = 0; l < slength; l++) {
    val = 1.0;
    for (k = 0; k < nchild; k++)
      val *= temp[k][l];
    temp1 += val;
  }
  return temp1;
}  /*segfast*/

Void segsexctop(LINK)
struct LOC_seg *LINK;
{
  double segval;
  int first, second, sstop;
  boolean thiscensor, thisrare;
  censorrec *WITH;
  thisarray *WITH1;
  int FORLIM;
  thisarray *WITH2;

  initseg(LINK);
  WITH = censorstruct;
  if ((*LINK->p)->male) {
    WITH1 = (*LINK->p)->gen;
    FORLIM = LINK->nfirst;
    for (first = 0; first < FORLIM; first++) {
      if (WITH1->genarray[first] != 0.0) {
	segval = 0.0;
	LINK->fseg = first + 1;
	thisrare = rare[first];
	second = 1;
	if (thisrare) {
	  WITH2 = (*LINK->q)->gen;
	  do {
	    LINK->sstart = probstart[second - 1];
	    LINK->send = probend[second - 1];
	    LINK->sseg = segstart[second - 1];
	    sstop = second + LINK->send - LINK->sstart + 1;
	    if (thisc < maxcensor) {
	      thisc++;
	      thiscensor = WITH->censor[thisc - minint];
	    } else
	      thiscensor = (WITH2->genarray[second - 1] == 0.0) ||
			    (rare[second - 1] == true);
	    if (!thiscensor) {
	      if (mutsys != 0)
		segval += WITH2->genarray[second - 1] * msegsexf(LINK);
	      else
		segval += WITH2->genarray[second - 1] * segsexf(LINK);
	    }
	    second = sstop;
	  } while (second <= LINK->nsecond);
	} else {
	  WITH2 = (*LINK->q)->gen;
	  do {
	    LINK->sstart = probstart[second - 1];
	    LINK->send = probend[second - 1];
	    LINK->sseg = segstart[second - 1];
	    sstop = second + LINK->send - LINK->sstart + 1;
	    if (thisc < maxcensor) {
	      thisc++;
	      thiscensor = WITH->censor[thisc - minint];
	    } else
	      thiscensor = (WITH2->genarray[second - 1] == 0.0);
	    if (!thiscensor) {
	      if (mutsys != 0)
		segval += WITH2->genarray[second - 1] * msegsexf(LINK);
	      else
		segval += WITH2->genarray[second - 1] * segsexf(LINK);
	    }
	    second = sstop;
	  } while (second <= LINK->nsecond);
	}
	WITH1->genarray[first] *= segval * segscale;
      }
    }
  } else {
    WITH1 = (*LINK->p)->gen;
    FORLIM = LINK->nfirst;
    for (first = 0; first < FORLIM; first++) {
      if (WITH1->genarray[first] != 0.0) {
	segval = 0.0;
	LINK->fstart = probstart[first];
	LINK->fend = probend[first];
	LINK->fseg = segstart[first];
	thisrare = rare[first];
	second = 1;
	if (thisrare) {
	  WITH2 = (*LINK->q)->gen;
	  do {
	    LINK->sseg = second;
	    if (thisc < maxcensor) {
	      thisc++;
	      thiscensor = WITH->censor[thisc - minint];
	    } else
	      thiscensor = (WITH2->genarray[second - 1] == 0.0) ||
			   (rare[second - 1] == true);
	    if (!thiscensor) {
	      if (mutsys != 0)
		segval += WITH2->genarray[second - 1] * msegsex(LINK);
	      else
		segval += WITH2->genarray[second - 1] * segsex(LINK);
	    }
	    second++;
	  } while (second <= LINK->nsecond);
	} else {
	  WITH2 = (*LINK->q)->gen;
	  do {
	    LINK->sseg = second;
	    if (thisc < maxcensor) {
	      thisc++;
	      thiscensor = WITH->censor[thisc - minint];
	    } else
	      thiscensor = (WITH2->genarray[second - 1] == 0.0);
	    if (!thiscensor) {
	      if (mutsys != 0)
		segval += WITH2->genarray[second - 1] * msegsex(LINK);
	      else
		segval += WITH2->genarray[second - 1] * segsex(LINK);
	    }
	    second++;
	  } while (second <= LINK->nsecond);
	}
	WITH1->genarray[first] *= segval * segscale;
      }
    }
  }
  cleanup(LINK->q, LINK->LINK);
  exitseg(LINK);
}  /*segsexctop*/


Void segsextop(LINK)
struct LOC_seg *LINK;
{
  double segval, val;
  int first, second, sstop;
  censorrec *WITH;
  thisarray *WITH1;
  int FORLIM;
  thisarray *WITH2;

  initseg(LINK);
  WITH = censorstruct;
  if ((*LINK->p)->male) {
    WITH1 = (*LINK->p)->gen;
    FORLIM = LINK->nfirst;
    for (first = 1; first <= FORLIM; first++) {
      if (WITH1->genarray[first - 1] != 0.0) {
	segval = 0.0;
	LINK->fseg = first;
	second = 1;
	WITH2 = (*LINK->q)->gen;
	do {
	  LINK->sstart = probstart[second - 1];
	  LINK->send = probend[second - 1];
	  LINK->sseg = segstart[second - 1];
	  sstop = second + LINK->send - LINK->sstart + 1;
	  if (thisc <= maxcensor)
	    thisc++;
	  if (WITH2->genarray[second - 1] == 0.0) {
	    second = sstop;
	    if (thisc <= maxcensor)
	      WITH->censor[thisc - minint] = true;
	  } else {
	    if (mutsys != 0)
	      val = msegsexf(LINK);
	    else
	      val = segsexf(LINK);
	    if (thisc <= maxcensor) {
	      WITH->censor[thisc - minint] = (val==0.0);
	    }
	    segval += WITH2->genarray[second - 1] * val;
	    second = sstop;
	  }
	} while (second <= LINK->nsecond);
	WITH1->genarray[first - 1] *= segval * segscale;
      }
    }
  } else {
    WITH1 = (*LINK->p)->gen;
    FORLIM = LINK->nfirst;
    for (first = 0; first < FORLIM; first++) {
      if (WITH1->genarray[first] != 0.0) {
	segval = 0.0;
	LINK->fstart = probstart[first];
	LINK->fend = probend[first];
	LINK->fseg = segstart[first];
	second = 1;
	WITH2 = (*LINK->q)->gen;
	do {
	  LINK->sseg = second;
	  if (thisc <= maxcensor)
	    thisc++;
	  if (WITH2->genarray[second - 1] == 0.0) {
	    second++;
	    if (thisc <= maxcensor)
	      WITH->censor[thisc - minint] = true;
	  } else {
	    if (mutsys != 0)
	      val = msegsex(LINK);
	    else
	      val = segsex(LINK);
	    if (thisc <= maxcensor) {
	      WITH->censor[thisc - minint] = (val == 0.0);
	    }
	    segval += WITH2->genarray[second - 1] * val;
	    second++;
	  }
	} while (second <= LINK->nsecond);
	WITH1->genarray[first] *= segval * segscale;
      }
    }
  }
  cleanup(LINK->q, LINK->LINK);
  exitseg(LINK);
}  /*segsextop*/

Void segctop(LINK)
struct LOC_seg *LINK;
{
  double segval;
  int first, second, sstop;
  boolean thiscensor, thisrare;
  censorrec *WITH;
  thisarray *WITH1;
  thisarray *WITH2;

  initseg(LINK);
  WITH = censorstruct;
  WITH1 = (*LINK->p)->gen;
  for (first = 0; first < fgeno; first++) {
    if (WITH1->genarray[first] != 0.0) {
      segval = 0.0;
      LINK->fstart = probstart[first];
      LINK->fend = probend[first];
      LINK->fseg = segstart[first];
      thisrare = rare[first];
      second = 1;
      if (thisrare) {
	WITH2 = (*LINK->q)->gen;
	do {
	  LINK->sstart = probstart[second - 1];
	  LINK->send = probend[second - 1];
	  LINK->sseg = segstart[second - 1];
	  sstop = second + LINK->send - LINK->sstart + 1;
	  if (thisc < maxcensor) {
	    thisc++;
	    thiscensor = WITH->censor[thisc - minint];
	  } else
	    thiscensor = (WITH2->genarray[second - 1] == 0.0) ||
			 rare[second - 1];
	  if (thiscensor)
	    second = sstop;
	  else {
	    if (mutsys != 0)
	      segval += WITH2->genarray[second - 1] * msegfast(LINK);
	    else
	      segval += WITH2->genarray[second - 1] * segfast(LINK);
	    second = sstop;
	  }
	} while (second <= fgeno);
      } else {
	WITH2 = (*LINK->q)->gen;
	do {
	  LINK->sstart = probstart[second - 1];
	  LINK->send = probend[second - 1];
	  LINK->sseg = segstart[second - 1];
	  sstop = second + LINK->send - LINK->sstart + 1;
	  if (thisc < maxcensor) {
	    thisc++;
	    thiscensor = WITH->censor[thisc - minint];
	  } else
	    thiscensor = (WITH2->genarray[second - 1] == 0.0);
	  if (thiscensor)
	    second = sstop;
	  else {
	    if (mutsys != 0)
	      segval += WITH2->genarray[second - 1] * msegfast(LINK);
	    else
	      segval += WITH2->genarray[second - 1] * segfast(LINK);
	    second = sstop;
	  }
	} while (second <= fgeno);
      }
      WITH1->genarray[first] *= segval * segscale;
    }
  }
  cleanup(LINK->q, LINK->LINK);
  exitseg(LINK);
  if (approximate && firstapprox && ((*LINK->p)->pa==NULL))
    getapprox(LINK);
}  /*segctop*/

Void segtop(LINK)
struct LOC_seg *LINK;
{
  double segval, val;
  int first, second, sstop;
  boolean thisrare;
  censorrec *WITH;
  thisarray *WITH1;
  thisarray *WITH2;

  initseg(LINK);
  WITH = censorstruct;
  WITH1 = (*LINK->p)->gen;
  for (first = 0; first < fgeno; first++) {
    if (WITH1->genarray[first] != 0.0) {
      thisrare = rare[first];
      segval = 0.0;
      LINK->fstart = probstart[first];
      LINK->fend = probend[first];
      LINK->fseg = segstart[first];
      second = 1;
      if (thisrare) {
	WITH2 = (*LINK->q)->gen;
	do {
	  LINK->sstart = probstart[second - 1];
	  LINK->send = probend[second - 1];
	  LINK->sseg = segstart[second - 1];
	  sstop = second + LINK->send - LINK->sstart + 1;
	  if (thisc <= maxcensor)
	    thisc++;
	  if ((WITH2->genarray[second - 1] == 0.0) || (rare[second - 1])) {
	    second = sstop;
	    if (thisc <= maxcensor)
	      WITH->censor[thisc - minint] = true;
	  } else {
	    if (mutsys != 0)
	      val = msegfast(LINK);
	    else
	      val = segfast(LINK);
	    if (thisc <= maxcensor) {
	      WITH->censor[thisc - minint] = (val == 0.0);
	    }
	    segval += WITH2->genarray[second - 1] * val;
	    second = sstop;
	  }
	} while (second <= fgeno);
      } else {
	WITH2 = (*LINK->q)->gen;
	do {
	  LINK->sstart = probstart[second - 1];
	  LINK->send = probend[second - 1];
	  LINK->sseg = segstart[second - 1];
	  sstop = second + LINK->send - LINK->sstart + 1;
	  if (thisc <= maxcensor)
	    thisc++;
	  if (WITH2->genarray[second - 1] == 0.0) {
	    second = sstop;
	    if (thisc <= maxcensor)
	      WITH->censor[thisc - minint] = true;
	  } else {
	    if (mutsys != 0)
	      val = msegfast(LINK);
	    else
	      val = segfast(LINK);
	    if (thisc <= maxcensor) {
	      WITH->censor[thisc - minint] = (val == 0.0);
	    }
	    segval += WITH2->genarray[second - 1] * val;
	    second = sstop;
	  }
	} while (second <= fgeno);
      }
      WITH1->genarray[first] *= segval * segscale;
    }
  }
  cleanup(LINK->q, LINK->LINK);
  exitseg(LINK);
  if (approximate && firstapprox && ((*LINK->p)->pa==NULL))
    getapprox(LINK);
}  /*segtop*/


Void segcapprox(LINK)
struct LOC_seg *LINK;
{
  double segval;
  int first, second, sstop;
  boolean thiscensor;
  censorrec *WITH;
  thisarray *WITH2;
  thisarray *WITH3;

  initseg(LINK);
  WITH = censorstruct;
  WITH2 = (*LINK->p)->gen;
  for (first = 0; first < fgeno; first++) {
    if (WITH2->genarray[first] != 0.0) {
      segval = 0.0;
      LINK->fstart = probstart[first];
      LINK->fend = probend[first];
      LINK->fseg = segstart[first];
      second = 1;
      WITH3 = (*LINK->q)->gen;
      do {
	LINK->sstart = probstart[second - 1];
	LINK->send = probend[second - 1];
	LINK->sseg = segstart[second - 1];
	sstop = second + LINK->send - LINK->sstart + 1;
	if (thisc < maxcensor) {
	  thisc++;
	  thiscensor = WITH->censor[thisc - minint];
	} else
	  thiscensor = (WITH3->genarray[second - 1] == 0.0);
	if (thiscensor || 
	    (!(approxarray[LINK->LINK->thisped - 1][first])))
	  second = sstop;
	else {
	  if (mutsys != 0)
	    segval += WITH3->genarray[second - 1] * msegfast(LINK);
	  else
	    segval += WITH3->genarray[second - 1] * segfast(LINK);
	  second = sstop;
	}
      } while (second <= fgeno);
      WITH2->genarray[first] *= segval * segscale;
    }
  }
  cleanup(LINK->q, LINK->LINK);
  exitseg(LINK);
}  /*segcapprox*/

#include "oldsegup.c"

Void msegsexdown(LINK)
struct LOC_seg *LINK;
{
  int here;
  double val, temp2;
  int ms2, ms1, mf, j, first, second;
  censorrec *WITH1;
  thisarray *WITH2;
  int FORLIM1;
  thisarray *WITH3;
  int FORLIM2;
  thetavalues *WITH4;
  int FORLIM3;

  initseg(LINK);
  for (here = 0; here < fgeno; here++)
    gene[here] = 0.0;
  WITH1 = censorstruct;
  WITH2 = (*LINK->p)->gen;
  FORLIM1 = LINK->nfirst;
  for (first = 0; first < FORLIM1; first++) {
    if (WITH2->genarray[first] != 0.0) {
      LINK->fseg = first + 1;
      second = 1;
      WITH3 = (*LINK->q)->gen;
      FORLIM2 = LINK->nsecond;
      for (second = 0; second < FORLIM2; second++) {
	if (WITH3->genarray[second] != 0.0) {
	  val = WITH3->genarray[second] * (*LINK->p)->gen->genarray[first] *
		segscale;
	  LINK->sstart = probstart[second];
	  LINK->send = probend[second];
	  LINK->sseg = segstart[second];
	  if (nchild != 0)
	    val *= msegsex(LINK);
	  if (val != 0.0) {
	    mf = muthap[invgenenum1[LINK->fseg - 1] - 1];
	    LINK->secondseg = LINK->sseg;
	    WITH4 = femaletheta;
	    FORLIM3 = LINK->send;
	    for (j = LINK->sstart - 1; j < FORLIM3; j++) {
	      ms1 = muthap[invgenenum1[LINK->secondseg - 1] - 1];
	      ms2 = muthap[invgenenum2[LINK->secondseg - 1] - 1];
	      temp2 = WITH4->segprob[j];
	      if (temp2 != 0.0) {
		if ((*LINK->r)->male) {
		  here = invgenenum1[LINK->secondseg - 1];
		  gene[here - 1] += (1 - LINK->ps) * temp2 * val;
		  here = invgenenum2[LINK->secondseg - 1];
		  gene[here - 1] += (1 - LINK->ps) * temp2 * val;
		  here = ms1;
		  gene[here - 1] += LINK->ps * temp2 * val;
		  here = ms2;
		  gene[here - 1] += LINK->ps * temp2 * val;
		} else {
		  here = genenumber[invgenenum1[LINK->secondseg - 1] - 1]
		    [first];
		  gene[here - 1] += (1 - LINK->pf) * (1 - LINK->ps) * temp2 * val;
		  here = genenumber[invgenenum2[LINK->secondseg - 1] - 1]
		    [first];
		  gene[here - 1] += (1 - LINK->pf) * (1 - LINK->ps) * temp2 * val;
		  here = genenumber[invgenenum1[LINK->secondseg - 1] - 1]
		    [mf - 1];
		  gene[here - 1] += LINK->pf * (1 - LINK->ps) * temp2 * val;
		  here = genenumber[invgenenum2[LINK->secondseg - 1] - 1]
		    [mf - 1];
		  gene[here - 1] += LINK->pf * (1 - LINK->ps) * temp2 * val;
		  here = genenumber[ms1 - 1][first];
		  gene[here - 1] += (1 - LINK->pf) * LINK->ps * temp2 * val;
		  here = genenumber[ms2 - 1][first];
		  gene[here - 1] += (1 - LINK->pf) * LINK->ps * temp2 * val;
		  here = genenumber[ms1 - 1][mf - 1];
		  gene[here - 1] += LINK->pf * LINK->ps * temp2 * val;
		  here = genenumber[ms2 - 1][mf - 1];
		  gene[here - 1] += LINK->pf * LINK->ps * temp2 * val;
		}
	      }
	      LINK->secondseg++;
	    }
	  }
	}  /*second*/
      }
    }  /*first*/
  }
  memcpy((*LINK->p)->gen->genarray, (*LINK->r)->gen->genarray,
	 maxfemgen * sizeof(double));
  memcpy((*LINK->r)->gen->genarray, gene, maxfemgen * sizeof(double));
  WITH2 = (*LINK->r)->gen;
  for (first = 0; first < fgeno; first++)
    WITH2->genarray[first] *= (*LINK->p)->gen->genarray[first];
  cleanup(LINK->p, LINK->LINK);
  cleanup(LINK->q, LINK->LINK);
  exitseg(LINK);
}  /*msegsexdown*/


Void msegdown(LINK)
struct LOC_seg *LINK;
{
  int here;
  double val, temp, temp1, temp2;
  int i, j, first, second, f1, f2, s1, s2, mf1, mf2, ms1, ms2;
  censorrec *WITH1;
  thisarray *WITH2;
  thisarray *WITH3;
  thetavalues *WITH4;
  int FORLIM3;
  thetavalues *WITH5;
  int FORLIM4;

  initseg(LINK);
  for (here = 0; here < fgeno; here++)
    gene[here] = 0.0;
  WITH1 = censorstruct;
  WITH2 = (*LINK->p)->gen;
  for (first = 0; first < fgeno; first++) {
    if (WITH2->genarray[first] != 0.0) {
      LINK->fstart = probstart[first];
      LINK->fend = probend[first];
      LINK->fseg = segstart[first];
      WITH3 = (*LINK->q)->gen;
      for (second = 0; second < fgeno; second++) {
	if (WITH3->genarray[second] != 0.0) {
	  LINK->sstart = probstart[second];
	  LINK->send = probend[second];
	  LINK->sseg = segstart[second];
	  val = WITH3->genarray[second] * (*LINK->p)->gen->genarray[first] *
		segscale;
	  if (nchild != 0)
	    val = msegfun(LINK) * val;
	  if (val != 0.0) {
	    LINK->firstseg = LINK->fseg;
	    WITH4 = maletheta;
	    FORLIM3 = LINK->fend;
	    for (i = LINK->fstart - 1; i < FORLIM3; i++) {
	      temp1 = WITH4->segprob[i];
	      LINK->secondseg = LINK->sseg;
	      if (temp1 != 0.0) {
		WITH5 = femaletheta;
		FORLIM4 = LINK->send;
		for (j = LINK->sstart - 1; j < FORLIM4; j++) {
		  if (WITH5->segprob[j] == 0.0)
		    LINK->secondseg++;
		  else {
		    temp2 = WITH5->segprob[j];
		    f1 = invgenenum1[LINK->firstseg - 1];
		    f2 = invgenenum2[LINK->firstseg - 1];
		    s1 = invgenenum1[LINK->secondseg - 1];
		    s2 = invgenenum2[LINK->secondseg - 1];
		    temp = (1 - LINK->pf) * (1 - LINK->ps) * temp1 * temp2 * val;
		    here = genenumber[s1 - 1][f1 - 1];
		    gene[here - 1] += temp;
		    here = genenumber[s1 - 1][f2 - 1];
		    gene[here - 1] += temp;
		    here = genenumber[s2 - 1][f1 - 1];
		    gene[here - 1] += temp;
		    here = genenumber[s2 - 1][f2 - 1];
		    gene[here - 1] += temp;
		    ms1 = muthap[s1 - 1];
		    ms2 = muthap[s2 - 1];
		    temp = (1 - LINK->pf) * LINK->ps * temp1 * temp2 * val;
		    here = genenumber[ms1 - 1][f1 - 1];
		    gene[here - 1] += temp;
		    here = genenumber[ms1 - 1][f2 - 1];
		    gene[here - 1] += temp;
		    here = genenumber[ms2 - 1][f1 - 1];
		    gene[here - 1] += temp;
		    here = genenumber[ms2 - 1][f2 - 1];
		    gene[here - 1] += temp;
		    mf1 = muthap[f1 - 1];
		    mf2 = muthap[f2 - 1];
		    temp = LINK->pf * (1 - LINK->ps) * temp1 * temp2 * val;
		    here = genenumber[mf1 - 1][s1 - 1];
		    gene[here - 1] += temp;
		    here = genenumber[mf1 - 1][s2 - 1];
		    gene[here - 1] += temp;
		    here = genenumber[mf2 - 1][s1 - 1];
		    gene[here - 1] += temp;
		    here = genenumber[mf2 - 1][s2 - 1];
		    gene[here - 1] += temp;
		    temp = LINK->pf * LINK->ps * temp1 * temp2 * val;
		    here = genenumber[mf1 - 1][ms1 - 1];
		    gene[here - 1] += temp;
		    here = genenumber[mf1 - 1][ms2 - 1];
		    gene[here - 1] += temp;
		    here = genenumber[mf2 - 1][ms1 - 1];
		    gene[here - 1] += temp;
		    here = genenumber[mf2 - 1][ms2 - 1];
		    gene[here - 1] += temp;
		    LINK->secondseg++;
		  }
		}
	      }
	      LINK->firstseg++;
	    }
	  }
	}  /*second*/
      }
    }  /*first*/
  }
  memcpy((*LINK->p)->gen->genarray, (*LINK->r)->gen->genarray,
	 maxfemgen * sizeof(double));
  memcpy((*LINK->r)->gen->genarray, gene, maxfemgen * sizeof(double));
  WITH2 = (*LINK->r)->gen;
  for (first = 0; first < fgeno; first++)
    WITH2->genarray[first] *= (*LINK->p)->gen->genarray[first];
  cleanup(LINK->p, LINK->LINK);
  cleanup(LINK->q, LINK->LINK);
  exitseg(LINK);


}  /*msegdown*/




