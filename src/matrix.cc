/*
 * matrix.cc implements the matrix command and class
 * Written by Charles Peterson beginning on October 28, 1997
 * Copyright (c) 1997 Southwest Foundation for Biomedical Research
 */

// #define TR1 1

#ifdef TR1
#define STDPRE std::tr1
#else
#define STDPRE std
#endif

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <errno.h>
#include <string>

#ifdef TR1
#include <tr1/unordered_map>
#else
#include <unordered_map>
#endif

#include "solar.h"
#include "tablefile.h"
#include "pipeback.h"

int Matrix::count = 0;
Matrix* Matrix::Matrices[] = {0};
FD_Array<int> Matrix::Pedno(1024,1024);
FD_Array<int> Matrix::Pedsize(1024,1024);
int Matrix::last_pedno;
bool Matrix::Pedindex_Current = false;
int Matrix::Pedindex_Highest_Ibdid = 0;
bool Matrix::Famid_Present;
bool Matrix::Famid_Needed = false;
int Matrix::Missing_Ibdid = 0;

// Hash tables for ID->IBDID and ID,FAMID->IBDID
// Plan is to move this and the pedigree code here to pedigree.cc eventually

STDPRE::unordered_map<std::string, int>ID_ibdid;
STDPRE::unordered_map<std::string, int>IDFAM_ibdid;



//
// The most fun part is...getting and setting values.
//
// Since IBDID's start from 1, 1-based id's are always assumed.
//
// Matrix data is either stored in little matrices (one per pedigree)
// or one big matrix (for the entire set of pedigrees).  This is indicated
// by the boolean ids_within_peds
//
// Within the matrices, they matrices are compressed into the upper diagonal
// (?) half...the get_index() function takes care of that.  This is because
// all matrices are symmetric (?).
//
// This is all done because otherwise matrices can get so big they can
// consume all available memory.
//
// To avoid the boolean and conditional test, we could have a subclass
//    but then we would have method lookup, which would likely take more time
//
float Matrix::get (int id1, int id2)
{
    if (ids_within_peds)
    {
	int ped1 = Pedno[id1];
	int ped2 = Pedno[id2];
	if (ped1 != ped2) return 0.0;
	return pedmat[ped1].values[get_index(id1,id2,pedmat[ped1].start)];
    }
    return pedmat[0].values[get_index(id1,id2,1)];
}


int Matrix::set (int id1, int id2, float value)
{
    if (ids_within_peds)
    {
	int ped1 = Pedno[id1];
	int ped2 = Pedno[id2];
	if (ped1 == ped2) {
	    int index = get_index (id1,id2,pedmat[ped1].start);
	    pedmat[ped1].values[index] = value;
	}
	else
	{
//	    printf ("IDS NOT IN PEDS AS ASSUMED\n");
//  This forces break from loop, close, and new read
	    return -1;  // id's not in peds as assumed
	}
    }
    else
    {
	pedmat[0].values[get_index (id1,id2,1)] = value;
    }
    return 0;
}


Matrix::Matrix (const char *name)
{
    _name = Strdup (name);
    filename = Strdup ("");
    min = 1.0;
    max = 0.0;
    defaultable = false;
    second_matrix = 0;
    first_matrix = 0;
    _ibd = false;
    _d7 = false;
    pedmat = 0;
    pedmat_count = 0;
    highest_id = 0;
    sum = 0.0;
    ibdid_found = 0;
    load_option = 0;
}

Matrix::~Matrix ()
{
    remove ();
    free (_name);
    free (filename);
    delete [] pedmat;  // Initialized to 0 by constructor, allocated by new
    if (ibdid_found) free (ibdid_found);

    if (second_matrix)
    {
	delete second_matrix;
    }
}

// find this matrix and remove from Matrices array
// if not in Matrices array, nothing is done
void Matrix::remove ()
{
    int i;
    for (i = 0; i < count; i++)
    {
	if (Matrices[i] == this) break;
    }
    if (i >= count)
    {
	return;  // Wasn't added to array
    }
    for (i++; i < count; i++)
    {
	Matrices[i-1] = Matrices[i];
    }
    Matrices[i-1] = 0;
    count--;
}

	    
void Matrix::reset ()
{
    int i;
    for (i = count-1; i >= 0; i--)
    {
	delete Matrices[i];
    }
}

void Matrix::add ()
{
    Matrices[count++] = this;
}

Matrix *Matrix::find (const char *nam)
{
    Matrix *m;
    int i = count-1;
    for (; i >= 0; i--)
    {
	m = index (i);
	if (!m) continue;
	if (!StringCmp (nam, m->name(), case_ins))
	{
	    return m;
	}
	if (m->second_matrix)
	{
	    if (!StringCmp (nam, m->second_matrix->name(), case_ins))
	    {
		return m->second_matrix;
	    }
	}
    }
    return 0;
}

char *Matrix::command (char *buf)
{
    char option[64];
    option[0] = 0;
    if (load_option==1)
    {
	sprintf (option,"-allow ");
    }
    else if (load_option==-1)
    {
	sprintf (option,"-sample ");
    }
    if (second_matrix)
    {
	sprintf (buf, "matrix load %s%s %s %s", option, filename, name(), 
		 second_matrix->name()); // handle 2 matrices
    }
    else
    {
	sprintf (buf, "matrix load %s%s %s ", option, filename, name() );
    }
    return buf;
}

class MatrixPoint {
public:
    MatrixPoint () {ibdid1=-1;ibdid2=-1;value=0;}
    int ibdid1;
    int ibdid2;
    float value;
    bool valid() {return (ibdid1!=-1&&ibdid2!=-1);}
};

char *Matrix::describe (char *buf)
{
    MatrixPoint dmin;
    MatrixPoint dmax;
    MatrixPoint ndmin;
    MatrixPoint ndmax;
    MatrixPoint dmin2;
    MatrixPoint dmax2;
    MatrixPoint ndmin2;
    MatrixPoint ndmax2;

    Matrix* m2 = second_matrix;

    for (int i=1; i <= highest_id; i++)
    {
	for (int j=1; j <= highest_id; j++)
	{
	    float val = get(i,j);
	    float val2;
	    if (m2) val2 = m2->get(i,j);

	    if (i == j)
	    {
		if (!dmin.valid() || val < dmin.value)
		{
		    dmin.value = val;
		    dmin.ibdid1 = i;
		    dmin.ibdid2 = j;
		}
		if (!dmax.valid() || val > dmax.value)
		{
		    dmax.value = val;
		    dmax.ibdid1 = i;
		    dmax.ibdid2 = j;
		}
		if (m2)
		{
		    if (!dmin2.valid() || val2 < dmin2.value)
		    {
			dmin2.value = val2;
			dmin2.ibdid1 = i;
			dmin2.ibdid2 = j;
		    }
		    if (!dmax2.valid() || val2 > dmax2.value)
		    {
			dmax2.value = val2;
			dmax2.ibdid1 = i;
			dmax2.ibdid2 = j;
		    }
		}
	    }
	    else // i != j
	    {
		if (!ndmin.valid() || val < ndmin.value)
		{
		    ndmin.value = val;
		    ndmin.ibdid1 = i;
		    ndmin.ibdid2 = j;
		}
		if (!ndmax.valid() || val > ndmax.value)
		{
		    ndmax.value = val;
		    ndmax.ibdid1 = i;
		    ndmax.ibdid2 = j;
		}
		if (m2)
		{
		    if (!ndmin2.valid() || val2 < ndmin2.value)
		    {
			ndmin2.value = val2;
			ndmin2.ibdid1 = i;
			ndmin2.ibdid2 = j;
		    }
		    if (!ndmax2.valid() || val2 > ndmax2.value)
		    {
			ndmax2.value = val2;
			ndmax2.ibdid1 = i;
			ndmax2.ibdid2 = j;
		    }
		}
	    }
	}
    }

// create report

    if (!m2)
    {
	sprintf (buf, 
"matrix file=%s size=%d\n\
    name=%s sum=%-18.14g\n\
        ndmin(%d,%d)=%10.8f  ndmax(%d,%d)=%10.8f\n\
        dmin(%d,%d)=%10.8f   dmax(%d,%d)=%10.8f",
	    filename, highest_id, name(), sum,
	    ndmin.ibdid1, ndmin.ibdid2, ndmin.value,
	    ndmax.ibdid1, ndmax.ibdid2, ndmax.value,
	    dmin.ibdid1,dmin.ibdid2,dmin.value,
	    dmax.ibdid1,dmax.ibdid2,dmax.value);
    }
    else // report both matrixes
    {
	sprintf (buf, 
"matrix file=%s size=%d\n\
    name=%s sum=%-18.14g\n\
        ndmin(%d,%d)=%10.8f  ndmax(%d,%d)=%10.8f\n\
        dmin(%d,%d)=%10.8f   dmax(%d,%d)=%10.8f\n\
    name=%s sum=%-18.14g\n\
        ndmin(%d,%d)=%10.8f  ndmax(%d,%d)=%10.8f\n\
        dmin(%d,%d)=%10.8f   dmax(%d,%d)=%10.8f",
	    filename, highest_id, name(), sum,
	    ndmin.ibdid1, ndmin.ibdid2, ndmin.value,
	    ndmax.ibdid1, ndmax.ibdid2, ndmax.value,
	    dmin.ibdid1,dmin.ibdid2,dmin.value,
	    dmax.ibdid1,dmax.ibdid2,dmax.value,
		 
	    m2->name(), m2->sum,
	    ndmin2.ibdid1,
	    ndmin2.ibdid2,
	    ndmin2.value,
	    ndmax2.ibdid1,
	    ndmax2.ibdid2,
	    ndmax2.value,

	    dmin2.ibdid1,
	    dmin2.ibdid2,
	    dmin2.value,
	    dmax2.ibdid1,
	    dmax2.ibdid2,
	    dmax2.value);
    }
    return buf;
}

char *Matrix::describe_all (char *buf)
{
    int index = 0;
    Matrix *m;
    for (int i = 0; m = Matrix::index(i); i++)
    {
	m->describe (&buf[index]);
	index = strlen (buf);
	buf[index++] = '\n';
    }
    buf[index] = '\0';
    return buf;
}


// MATCHECK is fortran interface to Matrix::CheckIbdid
//    if load_option == 0 (the default)
//      errors reported, then flag is set to take error exit later
//    if load_option == 1 (-allow)
//      diagonals default to 1, no error reported
//    if load_option == -1 (-sample)
//      errors reported, notfound returned as 1 to exclude from sample

extern "C" void matcheck_ (int* ibdid, int* notfound)
{
//    printf ("CHECKING %d\n",*ibdid);
    *notfound = Matrix::CheckIbdid (*ibdid);
}

//CheckIbdid is public static function which checks each matrix
int Matrix::CheckIbdid (int ibdid)
{
    Matrix *m;
    int notfound = 0;
    for (int i = 0; m = Matrix::index(i); i++)
    {
	if (m->check_ibdid (ibdid))
	{
	    notfound = 1;
	}
    }
    return notfound;
}


// check_ibdid checks first matrix and second matrix if present
int Matrix::check_ibdid (int ibdid)
{
    int found = 0;
    int found2 = 1;

    if (ibdid <= Pedindex_Highest_Ibdid)
    {
	found = ibdid_found[ibdid];
    }
    if (!found)
    {
	if (load_option == 0)
	{
	    printf ("Matrix %s is missing person IBDID=%d\n",name(),ibdid);
	    if (load_option == 0)
	    {
		Missing_Ibdid = ibdid;
	    }
	}
    }

    if (second_matrix)
    {
	found2 = 0;
	if (ibdid <= Pedindex_Highest_Ibdid)
	{
	    found2 = second_matrix->ibdid_found[ibdid];
	}
	if (!found2)
	{
	    if (load_option == 0)
	    {
		printf ("Matrix %s is missing person IBDID=%d\n",
			second_matrix->name(),ibdid);
		if (load_option == 0)
		{
		    Missing_Ibdid = ibdid;
		}
	    }
	}
    }
    if (!found || !found2)
    {
	if (load_option == -1)
	{
	    return 1;
	}
    }
    return 0;
}

int Matrix::Check_Matrices ()
{
    return Missing_Ibdid;
}

char *Matrix::commands (char *buf)
{
    int index = 0;
    Matrix *m;
    for (int i = 0; m = Matrix::index(i); i++)
    {
	m->command (&buf[index]);
	index = strlen (buf);
	buf[index++] = '\n';
    }
    if (index > 0) index--;
    buf[index] = '\0';
    return buf;
}

int Matrix::return_all (Tcl_Interp* interp)
{
    char* buf = Strdup ("");
    char tbuf[1024];
    Matrix* m;

    for (int i = 0; m = Matrix::index(i); i++)
    {
	string__append (&buf,"{");
	m->command(tbuf);
	string__append (&buf,tbuf);
	string__append (&buf,"}\n");
    }
    RESULT_BUF(buf);
    free (buf);
    return TCL_OK;
}



void Matrix::write_commands (FILE *file)
{
    int index = 0;
    Matrix *m;
    for (int i = 0; m = Matrix::index(i); i++)
    {
	char buf[256];
	fprintf (file, "%s\n", m->command (&buf[index]));
    }
}

void Matrix::Changing_Pedigree ()
{
    if (Pedindex_Current)
    {
	Pedno.renew();
	Pedsize.renew();
    }
    Pedindex_Current = false;
}

const char* Matrix::load_pedigree ()
{
    if (!Pedindex_Current)
    {
//	printf ("Clearing existing map\n");
	ID_ibdid.clear();
	IDFAM_ibdid.clear();
	Famid_Needed = false;

//	printf ("Loading pedigree info to matrix...\n");
	const char *errmsg;
	Famid_Present = false;
	TableFile *pedindex = TableFile::open ("pedindex.out", &errmsg);
	if (errmsg)
	{
	    return "Pedigree must be loaded before matrices";
	}
	pedindex->start_setup (&errmsg);
	pedindex->setup ("IBDID", &errmsg);
	pedindex->setup ("PEDNO", &errmsg);
	pedindex->setup ("ID", &errmsg);
	if (pedindex->test_name ("FAMID", &errmsg))
	{
	    pedindex->setup ("FAMID", &errmsg);
	    Famid_Present = true;
//	    printf ("FAMID is present in pedindex\n");
	} else {
//	    printf ("FAMID is NOT present in pedindex\n");
	}

	if (errmsg)
	{
	    delete pedindex;
	    return "Something wrong with pedindex file";
	}
	Pedno.renew();
	Pedsize.renew();
	char** data;
	last_pedno = -1;
	int pedsize = 0;
	int ibdid;
	int pedno;
	int famid;

	while (0 != (data = pedindex->get (&errmsg)))
	{
	    ibdid = atoi (data[0]);
	    pedno = atoi (data[1]);
	    if (last_pedno != -1 && last_pedno != pedno)
	    {
		Pedsize.set (last_pedno, pedsize);
		pedsize = 0;
	    }
	    Pedno.set (ibdid, pedno);
	    last_pedno = pedno;
	    pedsize++;

// Now construct hash table with ID->IBDID mapping
// If FAMID present, also construct FAMID.ID->IBDID table
//   if ID->IBDID table fails, then we know FAMID is needed, set flag

	    char* id = data[2];

// Each individual should only occur once.
// Test for previous insertion

	    if (!Famid_Needed)
	    {
		STDPRE::unordered_map<std::string,int>::const_iterator found =
		    ID_ibdid.find(id);
		if (found == ID_ibdid.end())
		{
	          ID_ibdid.insert (std::make_pair<std::string,int>(id,ibdid));
		}
		else
		{
		    if (Famid_Present)
		    {
//			printf ("FAMID is going to be used\n");
			Famid_Needed = true;
		    }
		    else
		    {
			delete pedindex;
			return "Same ID found twice in pedindex.out...does matrix need FAMID?";
		    }
		}
	    }
	    if (Famid_Present)
	    {
		std::string idfamid = id;
		idfamid.append(".famid.");
		idfamid.append(data[3]);
		STDPRE::unordered_map<std::string,int>::const_iterator found =
		    IDFAM_ibdid.find(idfamid);
		if (found == IDFAM_ibdid.end())
		{
		    IDFAM_ibdid.insert (std::make_pair<std::string,int>(idfamid,ibdid));
		}
		else
		{
		    delete pedindex;
		    return 
		     "Same ID,FAMID found twice  in pedindex.out";
		}
	    }
	}
	Pedindex_Highest_Ibdid = ibdid;
	if (last_pedno != -1)
	{
	    Pedsize.set (last_pedno, pedsize);
	}
	if (errmsg && Strcmp (errmsg, "EOF")) {
	    Pedno.renew();
	    Pedsize.renew();
	    fprintf (stderr, "Errmsg: %s\n", errmsg);
	    delete pedindex;
	    return "Error reading pedigree for matrix";
 	}
	delete pedindex;
	Pedindex_Current = true;
    }
//  printf ("done loading pedigree info to matrix\n");
    return 0;
}

// May be new "load," or a "re-load" of same filename

const char* Matrix::load (const char *specified_filename)
{
// Remove this matrix from matrix array until done
    remove ();

// working variables and names
    int scount;
    Matrix *m1 = this;
    Matrix *m2 = m1->second_matrix;

// Update Pedno and Pedsize tables if necessary
    const char* errmsg;
    if (   (errmsg = load_pedigree()) )
    {
	if (specified_filename[0] == '-')
	{
	    return "Invalid matrix option or pedigree not loaded";
	}
	return errmsg;
    }


// Open matrix file to be sure it exists and is not empty
    char* loading_filename;
    if (specified_filename)
    {
	loading_filename = append_extension (specified_filename, ".gz");
    }
    else
    {
	loading_filename = Strdup (filename);
    }
    FILE *mfile = fopen (loading_filename, "r");
    if (!mfile)
    {
	return "Unable to open matrix file";
    }
    if (EOF == fgetc (mfile))
    {
	return "Matrix file is empty";
    }
    Fclose (mfile);
//
// Clear error file if it exists
//
    struct stat statbuf;
    if (!stat("matrix.load.err",&statbuf))
    {
	unlink ("matrix.load.err");
    }
    errno = 0;
    int errors_logged = 0;

//  printf ("Ready to load matrix file\n");

// NOW, READ THE MATRIX FILE !!!

    char buf[1024];

// set up arguments to read cksum
// which must be read outside of read loop

    int linenumber = 0;
    int first_ibdid1 = -1;
    int first_ibdid2 = -1;
    char savebuf[256];
    char errorsbuf[256];
    unsigned int matrix_cksum;
    bool got_matrix_cksum = false;
    bool got_possible_cksum = false;


// Set up arguments to unzip

    const char *pbarg[4];
    pbarg[0] = "gunzip";
    pbarg[1] = "-c";
    pbarg[2] = loading_filename;
    pbarg[3] = 0;

// 2013
// We start by assuming all id's occur within pedigrees
//   if not, we break out of outer read loop and start over
//
    ids_within_peds = true;
    if (m2) m2->ids_within_peds = true;

    for (;;)  // one time through, unless id's not within peds
    {
	bool must_retry = false;

// Clear out old matrix storage
//   From way before, or last time through for loop above


//    printf ("clearing out old matrix storage\n");
      if (pedmat_count) {
	delete [] pedmat;
	pedmat = 0;
	pedmat_count = 0;
	if (m2) {
	    delete [] m2->pedmat;
	    m2->pedmat = 0;
	    m2->pedmat_count = 0;
	}
      }
      if (ibdid_found) free (ibdid_found);
      ibdid_found = 0;
      if (m2)
      {
	  if (m2->ibdid_found) 
	  {
	      free (m2->ibdid_found);
	      m2->ibdid_found = 0;
	  }
      }
//
// Allocate ibdid_found vector used to ensure matrix completeness
//
      ibdid_found = (int*) Calloc (Pedindex_Highest_Ibdid+1, sizeof(int));
      if (m2)
      {
	  m2->ibdid_found = (int*) Calloc (Pedindex_Highest_Ibdid+1, sizeof(int));
      }

// Initialize diagonal elements to -1
//    OR, if option #1, default, set to 1.0
//    Maybe should make this the new standard in all cases...
      int initial = -1;
      if (load_option==1)
      {
	  initial = 1;
      }
//
// Allocate a single large matrix if id's are not all in peds
//   We allocate a "half matrix" since it's symmetric
//  
      if (!ids_within_peds)
      {
//	printf ("Allocating large matrix\n");
	int half_size = 1 + get_index (Pedindex_Highest_Ibdid, Pedindex_Highest_Ibdid, 1);
	pedmat_count = 1;
	pedmat = new PedMatrix[pedmat_count];
	pedmat->values = new float [half_size];
	memset ((void*) pedmat->values, 0, sizeof(float)*half_size);
	if (m2) {
	    m2->pedmat_count = pedmat_count;
	    m2->pedmat = new PedMatrix[pedmat_count];
	    m2->pedmat->values = new float[half_size];
	    memset ((void*) m2->pedmat->values, 0, sizeof(float)*half_size);
	}

	for (int i=1; i <= Pedindex_Highest_Ibdid; i++)
	{
	    m1->set (i,i,initial);
	    if (m2) m2->set (i,i,initial);
	}
      }
      else
      {

// Otherwise, 
//   Allocate small "half-matrices" for each pedigree

//	printf ("Allocating small matrix for each pedigree\n");
	pedmat_count = last_pedno;
	pedmat = new PedMatrix [ pedmat_count + 1 ];  // use 1-based index
	if (m2) {
	      m2->pedmat_count = pedmat_count;
	      m2->pedmat = new PedMatrix[ pedmat_count + 1 ];
	}
	int first_id = 1;
	for (int i=1; i <= last_pedno; i++)
	  {
	    int size = Pedsize[i];
	    int half_size = 1 + get_index (size, size, 1);
	    pedmat[i].values = new float [half_size];
	    memset ((void*) pedmat[i].values, 0, sizeof(float)*half_size);
	    pedmat[i].start = first_id;
	    for (int j = first_id; j < first_id+size; j++)
	    {
		set (j, j, initial);
	    }
	    if (m2) {
		m2->pedmat[i].values = new float[half_size];
		memset ((void*) m2->pedmat[i].values, 0, 
			sizeof(float)*half_size);
		m2->pedmat[i].start = first_id;
		for (int j = first_id; j < first_id+size; j++)
		{
		    m2->set (j, j, initial);
		}
	    }
	    first_id += size;
	}
      }

// declare tab format variables

      int first_len;
      int fileformat = 0;
      
// Initialize csv record indexes

      int id1pos = -1;
      int id2pos = -1;
      int matrix1pos = -1;
      int matrix2pos = -1;
      int famid1pos = -1;
      int famid2pos = -1;
      int maxposneeded = 0;

// Unzip matrix file and pipe back here


//    printf ("Opening matrix file through gunzip\n");
      mfile = pipeback_shell_open ("gunzip", pbarg);
      if (!mfile)
      {
	return "Unable to uncompress matrix file";
      }

      linenumber = 0;

// Read in file from pipe, parse, then store matrix values

      while (fgets (buf, 256, mfile))
      {
	int ibdid1;
	int ibdid2;
	float matrix1;
	float matrix2;

	linenumber++;

// fgets always includes newline in buffer, if there is one

	int lastchar = strlen(buf) - 1;
	if (buf[lastchar] == '\n')
	{
	    buf[lastchar] = '\0';
	    lastchar--;
	}
	if (buf[lastchar] == '\r')
	{
	    buf[lastchar] = '\0';
	    lastchar--;
	}

//	printf ("line in buffer is %s\n",buf);

	strncpy (errorsbuf, buf, 256);  // save for error message report

// If this is first line, check if original or csv format

	if (linenumber == 1)
	{
	  if (strchr (buf, ','))
	  {
	    fileformat = 1;

// Read CSV header now and set indexing positions

	    int headoff = 0;
	    int headind = 0;
	    while (strlen(buf+headoff) > 0)
	    {
		char* commaptr = strchr (&buf[headoff], ',');
		if (commaptr)
		{
		    *commaptr = '\0';
		}
//		printf ("identified field %s\n",&buf[headoff]);
		maxposneeded++;
		if (!Strcmp (&buf[headoff],"id1"))
		{
		    id1pos = headind;
		}
		else if (!Strcmp (&buf[headoff],"id2"))
		{
		    id2pos = headind;
		}
		else if (!Strcmp (&buf[headoff],"matrix1"))
		{
		    matrix1pos = headind;
		}
		else if (!Strcmp (&buf[headoff],"matrix2"))
		{
		    matrix2pos = headind;
		}
		else if (!Strcmp (&buf[headoff],"famid1"))
		{
		    famid1pos = headind;
		}
		else if (!Strcmp (&buf[headoff],"famid2"))
		{
		    famid2pos = headind;
		}
		else
		{
//		    printf ("unrecognized matrix field: %s\n", &buf[headoff]);
		}
		headind++;
		if (commaptr)
		{
		    headoff = commaptr+1 - buf;
		}
		else
		{
		    break;
		}
	    }


// Check for missing required fields

	    if (id1pos==-1 || id2pos ==-1 || matrix1pos==-1 || 
		(m2 && matrix2pos==-1) ||
		(Famid_Needed && ((famid1pos==-1) || (famid2pos==-1))))
	    {
		pipeback_shell_close (mfile);
		if (id1pos==-1)
		{
		    printf ("matrix line: %s\n",errorsbuf);
		    return "matrix file missing id1 field";
		}
		else if (id2pos==-1)
		{
		    printf ("matrix line: %s\n",errorsbuf);
		    return "matrix file missing id2 field";
		}
		else if (matrix1pos==-1)
		{
		    printf ("matrix line: %s\n",errorsbuf);
		    return "matrix file missing matrix1 field";
		}
		else if (m2 && matrix2pos==-1)
		{
		    printf ("matrix line: %s\n",errorsbuf);
		    return "matrix file missing requested matrix2 field";
		}
		else if (Famid_Needed && famid1pos==-1)
		{
		    printf ("matrix line: %s\n",errorsbuf);
		    return "csv matrix file lacking needed famid1";
		}
		else if (Famid_Needed && famid2pos==-1)
		{
		    printf ("matrix line: %s\n",errorsbuf);
		    return "csv matrix file lacking needed famid2";
		}
		else
		{
		    printf ("matrix line: %s\n",errorsbuf);
		    return "csv matrix file missing required field";
		}
	    }
	    continue;  // csv header processed now so read next line
	  }
	  else // no comma in first record
	  {

// Original tab format.  Find decimal position following original rules:
// Divide line(s) into two parts based on decimal point
//   First part is ID1 ID2 [Space] (both integers)
//   Second part is VAL1 [VAL2] (both floats)
//   Second part begins 1 or 2 characters to the left of decimal
//   If 2nd character left is digit, second part begins 1 character left
//   otherwise 2 characters left
//
// Note: The "load matrix" code in matrix.cc permits any fixed position
// for the beginning of data values.  However, because of the way the checksum
// is currently written by matcrc, starting in position 14, that requires
// the rest of the matrix to follow that precedent, with possible deviation
// of one character position (data values could start in column 13, though
// that is not recommended).  In future it may be required to make matcrc
// actually look at the rest of the file to allow for IBDID's higher than
// 99999, which would require making the starting data position higher.

	    char *decimal_ptr = strchr (buf, '.');
	    int dpos =  decimal_ptr - buf;  // compare pointers to get dpos
	    if (!decimal_ptr || dpos < 4)
	    {
		pipeback_shell_close (mfile);
		printf ("matrix line: %s\n",errorsbuf);
		return "Invalid tab matrix file format";
	    }

	    first_len = dpos - 2;
	    if (isdigit (buf[first_len])) first_len++;
	  }
	}

// Above we determined filetype and parsed header if csv
// Now parse actual records depending on format
// If csv, this section will start with the second line in the file

	if (fileformat==1)  // CSV
	{

// parse CSV format data records

	    char* rptr = buf;
	    char* eptr = buf;
	    char* eeptr = &buf[strlen(buf)];
	    char* id1 = 0;
	    char* id2 = 0;
	    char* famid1 = 0;
	    char* famid2 = 0;

	    for (int ifield=0; ifield < maxposneeded; ifield++)
	    {
		eptr = strchr (rptr, ',');
		if (!eptr)
		{
		    if (!strlen(rptr))
		    {
			pipeback_shell_close (mfile);
			printf ("matrix line: %s\n",errorsbuf);
			return "Matrix record has required last field blank\n";
		    }
		    eptr = eeptr;
		}
		else
		{
		    *eptr = '\0';
		}
		if (ifield == id1pos)
		{
		    if (!strlen(rptr))
		    {
			pipeback_shell_close (mfile);
			printf ("matrix line: %s\n",errorsbuf);
			return "matrix id1 value missing";
		    }
		    id1 = rptr;  // Note: user id's are strings
		}
		else if (ifield == id2pos)
		{
		    if (!strlen(rptr))
		    {
			pipeback_shell_close (mfile);
			printf ("matrix line: %s\n",errorsbuf);
			return "matrix id2 value missing";
		    }
		    id2 = rptr;  // Note: user id's are strings
		}
		else if (ifield == matrix1pos)
		{
		    char *endptr;
		    errno = 0;
//
// store literal value in case it is checksum
//   note: format spec doesn't require id's to come first, so can't assume we
//     have them yet...later we determine if ID1 is checksum
//
		    if (linenumber==2)
		    {
			char* cksum_pointer;
			if ((cksum_pointer = strchr (rptr, '.')))
			{
//		            printf ("scanning matrix value at %s for cksum\n",
//			    cksum_pointer);
			    cksum_pointer++;
		            if (sscanf (cksum_pointer, "%u", &matrix_cksum))
		            {
//			        printf ("got possible cksum\n");
			        got_possible_cksum = true;
		            }
			}
		    }		      
		    matrix1 = strtof (rptr, &endptr);
		    if (errno || endptr == rptr)
		    {
			pipeback_shell_close (mfile);
			if (!strlen(rptr))
			{
			    printf ("matrix line: %s\n",errorsbuf);
			    return "matrix record has blank matrix1 value";
			}
			printf ("matrix line: %s\n",errorsbuf);
			return "Matrix1 value is invalid";
		    }
		    while (*endptr != '\0')
		    {
			if (!isspace (*endptr++))
			{
			    pipeback_shell_close (mfile);
			    printf ("matrix line: %s\n",errorsbuf);
			    return "Matrix1 value has invalid suffix";
			}
		    }
		}
		else if (ifield == matrix2pos)
		{
		    char *endptr;
		    errno = 0;
		    matrix2 = strtof (rptr, &endptr);
		    if (errno || endptr == rptr)
		    {
			pipeback_shell_close (mfile);
			if (!strlen(rptr))
			{
			    printf ("matrix line: %s\n",errorsbuf);
			    return "matrix record has blank matrix2 value";
			}
			printf ("matrix line: %s\n",errorsbuf);
			return "Matrix2 value is invalid";
		    }
		    while (*endptr != '\0')
		    {
			if (!isspace (*endptr++))
			{
			    pipeback_shell_close (mfile);
			    printf ("matrix line: %s\n",errorsbuf);
			    return "Matrix2 value has invalid suffix";
			}
		    }
		}
		else if (ifield == famid1pos)
		{
		    if (!strlen(rptr))
		    {
			pipeback_shell_close (mfile);
			printf ("matrix line: %s\n",errorsbuf);
			return "matrix famid1 value missing";
		    }
		    famid1 = rptr;
		}
		else if (ifield == famid2pos)
		{
		    if (!strlen(rptr))
		    {
			pipeback_shell_close (mfile);
			printf ("matrix line: %s\n",errorsbuf);
			return "matrix famid2 value missing";
		    }
		    famid2 = rptr;
		}

// Move rptr to next field

		if (eptr+1 < eeptr)
		{
//		    printf ("moving to next field\n");
		    rptr = eptr+1;
		}
		else
		{
		    if (ifield+1 < maxposneeded)
		    {
			pipeback_shell_close (mfile);
			printf ("matrix line: %s\n",errorsbuf);
			return "Matrix record has last field blank";
		    }
		}
	    }
//
// Check for checksum (in which case can't validate ID's)
//
	    if (linenumber==2 && !Strcmp ("checksum",id1))
	    {
		if (!got_possible_cksum)
		{
		    pipeback_shell_close (mfile);
		    printf ("matrix line: %s\n",errorsbuf);
		    return "Matrix has incorrectly formatted checksum";
		}
		got_matrix_cksum = true;
		continue;  // advance to next line
	    }
//
//  Now determine IBDID's
//
	    if (!Famid_Needed)
	    {
//		printf ("id1 is %s and id2 is %s\n",id1,id2);
		STDPRE::unordered_map<std::string,int>::const_iterator got =
		    ID_ibdid.find(id1);
		if (got == ID_ibdid.end())
		{
		    printf ("Ignoring matrix ID %s not in pedigree\n", id1);
		    printf ("  matrix line: %s\n",errorsbuf);
		    FILE* errfile = fopen ("matrix.load.err","a");
		    if (errfile)
		    {
			fprintf (errfile,
			       "Ignoring matrix ID %s not in pedigree\n", id1);
			fprintf (errfile,
				 "  matrix line: %s\n",errorsbuf);
			fclose (errfile);
		    }
		    errno = 0;
		    errors_logged++;
		    continue;
		}
		ibdid1 =  got->second;  // second is the mapped value

		got = ID_ibdid.find(id2);
		if (got == ID_ibdid.end())
		{
		    printf ("Ignoring matrix ID %s not in pedigree\n", id2);
		    printf ("  matrix line: %s\n",errorsbuf);
		    FILE* errfile = fopen ("matrix.load.err","a");
		    if (errfile)
		    {
			fprintf (errfile,
				"Ignoring matrix ID %s not in pedigree\n", id2);
			fprintf (errfile,"  matrix line: %s\n",errorsbuf);
			fclose (errfile);
		    }
		    errno = 0;
		    errors_logged++;
		    continue;
		}
		ibdid2 = got->second;
//		printf ("ibdid1 is %d and ibdid2 is %d\n",ibdid1, ibdid2);
	    }
	    else // if famid is needed
	    {
		char idfamid[1024];
		sprintf (idfamid,"%s.famid.%s",id1,famid1);
		STDPRE::unordered_map<std::string,int>::const_iterator got =
		    IDFAM_ibdid.find(idfamid);
		if (got == IDFAM_ibdid.end())
		{
		    printf ("Ignoring matrix ID %s FAMID %s not in pedigree\n",
			    id1,famid1);
		    printf ("  matrix line: %s\n",errorsbuf);
		    FILE* errfile = fopen ("matrix.load.err","a");
		    if (errfile)
		    {
			fprintf (errfile,
			     "Ignoring matrix ID %s FAMID %s not in pedigree\n",
			    id1,famid1);
			fprintf (errfile,"  matrix line: %s\n",errorsbuf);
			fclose (errfile);
		    }
		    errno = 0;
		    errors_logged++;
		    continue;
		}
		ibdid1 = got->second;

		sprintf (idfamid,"%s.famid.%s",id2,famid2);
		got = IDFAM_ibdid.find(idfamid);
		if (got == IDFAM_ibdid.end())
		{
		    printf ("Ignoring matrix ID %s FAMID %s not in pedigree\n",
			    id2,famid2);
		    printf ("  matrix line: %s\n",errorsbuf);
		    FILE* errfile = fopen ("matrix.load.err","a");
		    if (errfile)
		    {
			fprintf (errfile,
			    "Ignoring matrix ID %s FAMID %s not in pedigree\n",
			    id2,famid2);
			fprintf (errfile,"  matrix line: %s\n",errorsbuf);
			fclose (errfile);
		    }
		    errno = 0;
		    errors_logged++;
		    continue;
		}
		ibdid2 = got->second;
	    }
	}
	else // parse record for tab format
	{
//	    printf ("first_len is %d\n",first_len);
	    char savech = buf[first_len];
	    char dummy[1024];
	    buf[first_len] = '\0';
	    scount = sscanf (buf, "%d %d %s", &ibdid1, &ibdid2, dummy);
	    if (scount != 2)
	    {
		pipeback_shell_close (mfile);
		printf ("matrix line: %s\n",errorsbuf);
		return "Error reading matrix file record";
	    }
	    if (ibdid1 > Pedindex_Highest_Ibdid ||
		ibdid2 > Pedindex_Highest_Ibdid)
	    {
		pipeback_shell_close (mfile);
		printf ("matrix line: %s\n",errorsbuf);
		return "Invalid ID found in matrix file";
	    }
	    buf[first_len] = savech;
	    if (!m2)
	    {
//		printf ("parsing first matrix1 value\n");
		scount = sscanf (&buf[first_len], "%f", &matrix1);
		if (scount != 1)
		{
		    pipeback_shell_close (mfile);
		    printf ("matrix line: %s\n",errorsbuf);
		    return "Error reading matrix1 value";
		}

	    } else { // include matrix2
		scount = sscanf (&buf[first_len], "%f %f", &matrix1,&matrix2);
		if (scount != 2)
		{
		    pipeback_shell_close (mfile);
		    printf ("matrix line: %s\n",errorsbuf);
		    return "Error reading matrix2 value";
		}
	    }

// Get cksum from file, if first two lines have same id pairs

	  if (linenumber == 1)
	  {
//	      printf ("checking for checksum\n");
	      first_ibdid1 = ibdid1;
	      first_ibdid2 = ibdid2;
	      strncpy (savebuf, buf, 256);
	  }
	  else if (linenumber == 2)
	  {
//	      printf ("checking for checksum 2\n");
	      if (first_ibdid1 == ibdid1 && first_ibdid2 == ibdid2)
	      {
//		  printf ("checking for decimal point\n");
		  char* cksum_pointer;
		  if ((cksum_pointer = strchr (savebuf, '.')))
		  {
//		      printf ("scanning matrix value at %s\n",cksum_pointer);
		      cksum_pointer++;
		      if (sscanf (cksum_pointer, "%u", &matrix_cksum))

		      {
//			  printf ("got matrix cksum\n");
			  got_matrix_cksum = true;
// Reset sum to zero since first line wasn't actual data
			  sum = 0;
			  if (m2) m2->sum = 0;
		      }
		  }
	      }
	  }
	}

// Now store matrix values requested
//   At this point, ibdid1 and ibdid2 should have been checked for range
//   However, m1->set will further validate the within-pedigree status

	if (-1 == m1->set (ibdid1, ibdid2, matrix1))
	{
//	    printf ("it appears ids are not within peds!\n");
	    ids_within_peds = false;
	    if (m2) m2->ids_within_peds = false;
	    must_retry = true;
	    break;
	}
	if (m1->max < matrix1) m1->max = matrix1;
	if (m1->min > matrix1) m1->min = matrix1;
	if (m1->_ibd && ibdid1 == ibdid2 && matrix1 == -1.0)
	{
	    for (int rc = 1; rc <= m1->Pedindex_Highest_Ibdid; rc++)
	    {
		m1->set (ibdid1, rc, -1.0);
	    }
	}
	if (ibdid1 > highest_id) highest_id = ibdid1;
	if (ibdid2 > highest_id) highest_id = ibdid2;
	sum = sum + matrix1;
	if (ibdid1==ibdid2)
	{
	    ibdid_found[ibdid1] = 1;
	}

	if (m2)
	{
	    m2->set (ibdid1, ibdid2, matrix2);
	    if (m2->max < matrix2) m2->max = matrix2;
	    if (m2->min > matrix2) m2->min = matrix2;
	    if (m2->_d7 && ibdid1 == ibdid2 && matrix2 == -1.0)
	    {
		for (int rc = 1; rc <= m2->Pedindex_Highest_Ibdid; rc++)
		{
		    m2->set (ibdid1, rc, -1.0);
		}
	    }
	    m2->highest_id = highest_id;
	    m2->sum = m2->sum + matrix2;
	    if (ibdid1==ibdid2)
	    {
		m2->ibdid_found[ibdid1] = 1;
	    }
	}

      } // end of file reading loop
      pipeback_shell_close (mfile);
//      printf ("Closed mfile\n");
      if (must_retry) continue;
//      printf ("breaking from read loop\n");
      break;

    } // end of try ids_within_peds true, then false, loop
    

    if (linenumber < 1)
    {
	return 	"Matrix empty or load failed for lack of memory";
    }

    if (got_matrix_cksum)
    {
//
// Get checksum from pedindex.out
//   Comment: This should be part of the pedigree reading procedure and
//   not have to be done all the time...
//   But right now, I don't trust "changing_pedigree" to actually
//   keep up with pedigree changes because people can simply change
//   directory, or change the file with a concurrent process.
//
//   So this remain here, for now, until I can fix the
//   changing directory issue and changing file issues.

	  FILE* pfile = fopen ("pedindex.out", "r");
	  if (!pfile)
	  {
	      return "Can't find pedindex.out";
	  }
	  fclose (pfile);

	  unsigned pedindex_cksum;
	  const char* carg[3];
	  carg[0] = "cksum";
	  carg[1] = "pedindex.out";
	  carg[2] = 0;
	  FILE* cfile = pipeback_shell_open ("cksum", carg);
	  if (!cfile)
	  {
	      return "Unable to run cksum on pedindex.out";
	  }

	  if (!fgets (buf, 256, cfile))
	  {
	      pipeback_shell_close (cfile);
	      return "Error reading checksum of pedindex.out";
	  }
	  if (1 != (sscanf (buf, "%u", &pedindex_cksum)))
	  {
	      pipeback_shell_close (cfile);
	      return "Error scanning checksum of pedindex.out";
	  }
	  pipeback_shell_close (cfile);
	  if (matrix_cksum != pedindex_cksum)
	  {
	      printf ("matrix cksum: %u\n",matrix_cksum);
	      printf ("pedindex cksum: %u\n",pedindex_cksum);
	      return "Checksum (cksum) in matrix doesn't match pedindex.out";
	  }
	  else
	  {
//	      printf ("Checksum verified\n");
	  }
      }

      free (filename);
      filename = loading_filename;
      add ();
      if (errors_logged)
      {
	  printf ("matrix has %d warnings written to matrix.load.err\n",
	      errors_logged);
	  printf ("matrix loaded incompletely\n");
      }
      return 0;
}


// If named matrix already exists, setup reloads it
// otherwise, it creates new matrix

const char* Matrix::setup (int option, const char *filename, const char *name1, 
		     const char *name2)
{
    Matrix* oldm = Matrix::find (name1);
    Matrix* m1;
    if (oldm)
    {
// Matrix with same name already exists.  Set up for re-load.
	m1 = oldm;
	if (name2) {
// Setup second matrix
	    if (m1->second_matrix) {
		free (m1->second_matrix->_name);
		m1->second_matrix->_name = Strdup (name2);
	    } else {
		m1->second_matrix = new Matrix (name2);
		if (m1->_ibd) m1->second_matrix->_d7 = true;
	    }
	} else {
// Delete previous second matrix, if any
	    if (m1->second_matrix) {
		delete m1->second_matrix;
		m1->second_matrix = 0;
	    }
	}
    }
    else
    {
// New primary matrix required
	m1 = new Matrix (name1);

// See if this is a twopoint (ibd,d7) or multipoint (mibd,d7) matrix
// Such matrices have follow the -1 convention in which a -1 on the diagonal
// Causes a (sub-)matrix to default to phi2/delta7

        char name3[128];
	strncpy (name3, name1, 3);
	name3[3] = '\0';

	char name4[128];
	strncpy (name4, name1, 4);
	name4[4] = '\0';
	
//	fprintf (stderr, "name3 is >%s< and name4 is >%s<\n", name3, name4);

	if (!Strcmp(name3, "ibd") || !Strcmp(name4, "mibd"))
	{
	    m1->_ibd = true;
	}
// Setup second matrix if required
	if (name2)
	{
	    Matrix *oldm2 = Matrix::find (name2);
	    if (oldm2)
	    {

// Previous matrix can be deleted if it is the second matrix of another
//   first matrix, or a first matrix having no second.  If it is the first
//   matrix of a matrix pair, it gets renamed.  If it is the second
//   matrix of another first matrix, be sure to zero that pointer.

		if (oldm2->first_matrix)
		{
//		    fprintf (stderr, "Setting backpointer to 0\n");
		    oldm2->first_matrix->second_matrix = 0;
		    delete oldm2;
		}
		else if (!oldm2->second_matrix)
		{
//		    fprintf (stderr, "Simply deleting it\n");
		    delete oldm2;
		}
		else
		{
//		    fprintf (stderr, "Renaming it\n");
		    char buf[1024];
		    strcpy (buf, "old_");
		    strncat (buf, oldm2->_name, 1024);
		    free (oldm2->_name);
		    oldm2->_name = Strdup (buf);
		}
	    }
	    m1->second_matrix = new Matrix (name2);
	    if (m1->_ibd) m1->second_matrix->_d7 = true;
	}
    }
    if (name2)
    {
	m1->second_matrix->first_matrix = m1;
    }
    m1->load_option = option;
    const char *message = m1->load (filename);
    if (message)
    {
	delete m1;
    }
    return message;
}

extern "C" int MathMatrixCmd (ClientData clientData, Tcl_Interp *interp,
			      int argc, char* argv[]);

extern "C" int MatrixCmd (ClientData clientData, Tcl_Interp *interp,
		  int argc, char *argv[])
{
    if (argc == 3 && !Strcmp ("new", argv[1]))
    {
	return MathMatrixCmd (clientData, interp, argc, argv);
    }
    
    if (argc == 3 && !Strcmp ("rows", argv[1]))
    {
	return MathMatrixCmd (clientData, interp, argc, argv);
    }
    
    if (argc == 3 && !Strcmp ("cols", argv[1]))
    {
	return MathMatrixCmd (clientData, interp, argc, argv);
    }
    
    if (argc == 3 && !Strcmp ("print", argv[1]))
    {
	return MathMatrixCmd (clientData, interp, argc, argv);
    }
    
    if (argc>1 && !Strcmp ("load", argv[1]))
    {

// This is a mathmatrix command if there is more than one non-hyphenated
// argument after all the hyphenated arguments.

	if (argc==3 || (argc>3 && (!Strcmp ("-noheader", argv[2]) ||
				   !Strcmp ("-cols", argv[2]) ||
				   !Strcmp ("-rows", argv[2]))))
	{
	    return MathMatrixCmd (clientData, interp, argc, argv);
	}
    }	

    if (argc>1 && argv[1][0] == '.')  // MathMatrixID as first argument
    {
	return MathMatrixCmd (clientData, interp, argc, argv);
    }

    if (argc==2 && !Strcmp ("reset", argv[1]))
    {
	return MathMatrixCmd (clientData, interp, argc, argv);
    }

    if (argc == 2 && !StringCmp ("help", argv[1], case_ins))
    {
	return Solar_Eval (interp, "help matrix");
    }

    if (argc == 1)
    {
	char buf[10000];
	sprintf (buf, "%s", Matrix::commands (buf));
	RESULT_BUF (buf);
	return TCL_OK;
    }

    if (argc == 2 && !StringCmp ("-return", argv[1], case_ins))
    {
	return  Matrix::return_all (interp);
    }

    if (argc == 2 && !StringCmp ("debug", argv[1], case_ins))
    {
	char buf[10000];
	printf ("%s", Matrix::describe_all (buf));
	return TCL_OK;
    }

    if (argc == 3 && !StringCmp (argv[1], "delete", case_ins))
    {
	if (0==strncmp( argv[2], ".mm.", 4))
	{
	    return MathMatrixCmd (clientData, interp, argc, argv);
	}

	Matrix *m = Matrix::find (argv[2]);
	if (!m)
	{
	    RESULT_LIT ("No such matrix");
	    return TCL_ERROR;
	}
	delete m;
	return TCL_OK;
    }

    if (argc == 2 && !StringCmp (argv[1], "delete_all", case_ins))
    {
	Matrix::reset();
	return TCL_OK;
    }

    if ((argc >= 4 && argc <= 6) && !StringCmp (argv[1], "load", case_ins))
    {
    // Setup new Matrices
	int option = 0;
	int offset = 0;
	if (!Strcmp (argv[2], "-allow"))
	{
	    if (argc==4)
	    {
		RESULT_LIT ("Incomplete matrix command");
		return TCL_ERROR;
	    }
	    offset = 1;
	    option = 1;
	}
	else if (!Strcmp (argv[2], "-sample"))
	{
	    if (argc==4)
	    {
		RESULT_LIT ("Incomplete matrix command");
		return TCL_ERROR;
	    }
	    offset = 1;
	    option = -1;
	}

	const char *message = 0;
	if (argc - offset == 4)
	{
	    message = Matrix::setup (option, argv[2+offset], argv[3+offset]);
	}
	else
	{
	    message = Matrix::setup (option, argv[2+offset], argv[3+offset], 
				     argv[4+offset]);
	}
	if (message)
	{
	    char buf[1024];
	    sprintf (buf, "%s:  %s", message, argv[2]);
	    RESULT_BUF (buf);
	    return TCL_ERROR;
	}
	return TCL_OK;
    }	
    RESULT_LIT ("Invalid matrix command");
    return TCL_ERROR;
}

int Matrix::bind (Tcl_Interp *interp)
{
// If pedigree changed, must reload all matrices.
// This is bad, but should be avoided by not re-loading same pedigree

    if (!Pedindex_Current && count>0)
    {
	int i;
	for (i=0; i < count; i++)
	{
	    {
		fprintf (stderr, "Pedigree changed; reloading matrix %d\n", i);
	    }
	    Matrix *m = Matrices[i];
	    m->load();
	}
    }

    Missing_Ibdid = 0;  // Determined during pinput

    return TCL_OK;
}

