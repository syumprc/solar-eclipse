/*************************************************************************/
/*  (C) 04.08.1997 Henrik Seidel (HS)                                    */
/*  <henrik@itb.biologie.hu-berlin.de>                                   */
/*************************************************************************/

/* $Id: acegr_np.h,v 1.1 1997/08/11 16:12:30 henrik Exp henrik $ */

#ifndef ACEGR_NPIPE_H_
#define ACEGR_NPIPE_H_

#ifndef EXIT_SUCCESS
#  define EXIT_SUCCESS 0
#endif

#ifndef EXIT_FAILURE
#  define EXIT_FAILURE -1
#endif

int ACEgrOpen (const int, char *xmgr_name);
int ACEgrClose (int terminated_by_user);
int ACEgrFlush (void);
int ACEgrPrintf (const char*, ...);
int ACEgrCommand (const char*);

#endif /* ACEGR_NPIPE_H */

