/*
 * tablefilecmd.cc provides a tcl interface to the TableFile class
 * Written by Charles Peterson May 26, 1998
 * Copyright (c) 1998 Southwest Foundation for Biomedical Research
 */

#include "solar.h"
#include "tablefile.h"

static TableFile *tablefiles[MAX_TABLEFILES];
static int tablefilecount = 0;

extern "C" int TableFileCmd (ClientData clientData, Tcl_Interp *interp,
		  int argc, char *argv[])
{
    int i;

    if (argc == 1)
    {
	for (i = 0; i < tablefilecount; i++)
	{
	    if (tablefiles[i])
	    {
		printf ("[%d]  %s\n", i, tablefiles[i]->filename() );
	    }
	}
	return TCL_OK;
    }

    if (argc == 2 && !Strcmp (argv[1], "fisher2"))
    {
	extern bool Fisher_2_Format_Allowed;
	Fisher_2_Format_Allowed = true;
	return TCL_OK;
    }
    
    if (argc < 3)
    {
	RESULT_LIT ("Invalid tablefile command");
	return TCL_ERROR;
    }
    
    const char *errmsg = 0;

    if (!StringCmp (argv[1], "open", case_ins))
    {
	TableFile *tf = TableFile::open (argv[2], &errmsg);
	if (errmsg)
	{
	    RESULT_BUF (errmsg);
	    return TCL_ERROR;
	}
	for (i = 0; i < tablefilecount; i++)
	{
	    if (!tablefiles[i]) break;
	}
	if (i >= MAX_TABLEFILES)
	{
	    RESULT_LIT ("Too many open tablefiles");
	    return TCL_ERROR;
	}
	tablefiles[i] = tf;
	if (tablefilecount <= i) tablefilecount = i+1;
	char buf[128];
	sprintf (buf, "%d", i);
	RESULT_BUF (buf);
	return TCL_OK;
    }

    int tablenum;
    char junk[1024];
    int count = sscanf (argv[1], "%d %s", &tablenum, junk);
    if (count != 1 || tablenum >= tablefilecount || !tablefiles[tablenum])
    {
	RESULT_LIT ("Invalid tablefile index");
	return TCL_ERROR;
    }
    TableFile *tf = tablefiles[tablenum];


    while (1)  // break on error
    {
	if (!StringCmp (argv[2], "names", case_ins))
	{
	    int count = 0;
	    const char **names = tf->names (&count, &errmsg);
	    if (errmsg) break;
	    for (i = 0; i < count; i++)
	    {
		Solar_AppendElement (interp, (char*) names[i]);
	    }
	    return TCL_OK;
	}
	if (!StringCmp (argv[2], "short_names", case_ins))
	{
	    const char **names = tf->short_names (&count, &errmsg);
	    if (errmsg) break;
	    for (i = 0; i < count; i++)
	    {
		Solar_AppendElement (interp, (char*) names[i]);
	    }
	    return TCL_OK;
	}
	if (!StringCmp (argv[2], "widths", case_ins))
	{
	    int *widths = tf->widths (&count, &errmsg);
	    if (errmsg) break;
	    for (i = 0; i < count; i++)
	    {
		char buf[64];
		sprintf (buf, "%d", widths[i]);
		Solar_AppendElement (interp, buf);
	    }
	    return TCL_OK;
	}
	if (!StringCmp (argv[2], "start_setup", case_ins))
	{
	    tf->start_setup (&errmsg);
	    if (errmsg) break;
	    return TCL_OK;
	}
	if (!StringCmp (argv[2], "setup", case_ins))
	{
	    if (argc < 4)
	    {
		errmsg = "Missing field name in tablefile setup command";
		break;
	    }
	    int ipos = tf->setup (argv[3], &errmsg);
	    if (errmsg) break;
	    char buf[128];
	    sprintf (buf, "%d", ipos);
	    RESULT_BUF (buf);
	    return TCL_OK;
	}
	if (!StringCmp (argv[2], "get", case_ins))
	{
	    char **record = tf->get (&errmsg);
	    if (errmsg && !strcmp (errmsg, "EOF")) return TCL_OK;
	    if (errmsg) break;
	    for (i = 0; record[i]; i++)
	    {
		Solar_AppendElement (interp, record[i]);
	    }
	    return TCL_OK;
	}
	if (!StringCmp (argv[2], "rewind", case_ins))
	{
	    tf->rewind (&errmsg);
	    if (errmsg) break;
	    return TCL_OK;
	}
	if (!StringCmp (argv[2], "close", case_ins))
	{
	    delete tf;
	    tablefiles[tablenum] = 0;
	    return TCL_OK;
	}
	if (!Strcmp (argv[2], "test_name"))
	{
	    int existed = tf->test_name (argv[3], &errmsg);
	    if (errmsg) break;
	    char buf[128];
	    sprintf (buf, "%d", existed);
	    RESULT_BUF (buf);
	    return TCL_OK;
	}
	if (!Strcmp (argv[2], "get_position"))
	{
	    int position = tf->get_position ();
	    char buf[128];
	    sprintf (buf, "%d", position);
	    RESULT_BUF (buf);
	    return TCL_OK;
	}
	if (!Strcmp (argv[2], "set_position"))
	{
 	    tf->set_position (atoi(argv[3]), &errmsg);
	    if (errmsg) break;
	    return TCL_OK;
	}
	RESULT_LIT ("Invalid tablefile command");
	return TCL_ERROR;
    }
    RESULT_BUF (errmsg);
    return TCL_ERROR;
}

		
	    

	
	
