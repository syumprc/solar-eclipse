/*
 * Copyright (c) 1998 Southwest Foundation for Biomedical Research
 *
 * fisherpedigree.cc implements the FisherPedigree class, which provides
 *   the pedigree information needed by the fisher component of SOLAR.
 * The pedigree information is read from the pedindex.out file created by
 *   the "load pedigree" command, but doled out to "preped" and "pinput"
 *   routines (originally from fisher) as if it were being read from a
 *   fisher-formatted file.
 *
 * Written by Charles Peterson beginning on June 4, 1998
 */

#include "solar.h"

class IntArray
{
    int *_array;
    int _count;
    int _size;
    int _position;
public:
    IntArray() {_position=0; _count=0; 
                _array=(int*) Calloc (sizeof(int), (_size=1));}
    ~IntArray() {free (_array);}
    void add (int number) {if (++_count>_size)
	_array=(int*)Realloc(_array,sizeof(int)*(++_size));
        _array[_count-1]=number;}
    void zero() {_count=0;_position=0;}
    int count () {return _count;}
    int index(int i) {return (i<_count) ? _array[i] : INT_MIN;}
    int change(int i, int newvalue) {return (i<_count) ? _array[i] = newvalue :
	INT_MIN;}
    void rewind() {_position=0;}
    int next() {return ((_position)<_count) ? _array[_position++] : INT_MIN;}
    IntArray *copy_to (IntArray *dest) {dest->zero (); int i=0; int value; 
        while (INT_MIN != (value = this->index (i++)))
	    dest->add(value); return dest;}
};

TableFile*   FisherPedigree::Ibd_Index = 0;
bool         FisherPedigree::Ibd_Famid = false;
int          FisherPedigree::_Highest_Pedcnt = 0;
int          FisherPedigree::last_grouping = 0;
int          FisherPedigree::last_merge_groups = 0;
IntArray*    FisherPedigree::_Pedcnts = new IntArray;
char*        FisherPedigree::_FirstIDs = 0;

/*
 * When household effects are being analyzed, pedigrees are merged into
 * "pedigree-household groups" which include all the members in all the
 * pedigrees that share households.
 *
 * The HouseholdGroup array maps pedigree numbers (index) to household group
 * number (value).  Multiple pedigrees may have the same household group
 * number, which is the same as the first pedigree number in that group
 * (When two pedigrees share a household, the pedigree-household group is
 * numbered by the lower pedigree number, i.e. pedigrees are "collapsed"
 * downward.)
 *
 * The HouseholdGroupCnts array indicates the number of individuals in each
 * household-pedigree group.  It is indexed by the household group number
 * (see preceeding paragraph) which is the number of the first pedigree in
 * the group. There may be gaps in this array corresponding to pedigrees
 * whose members were accumulated by a lower-numbered group.
 */
bool         FisherPedigree::_IfHouseholdGroups = false;
IntArray*    FisherPedigree::_HouseholdGroup = new IntArray;
IntArray*    FisherPedigree::_HouseholdGroupCnts = new IntArray;
int          FisherPedigree::_CurrentHouseholdGroup = 0;
int          FisherPedigree::_CurrentPedigree = 0;
int          FisherPedigree::_RemainingInCurrentPedigree = 0;
int          FisherPedigree::_TotalIndividuals = 0;
char         FisherPedigree::HouseholdGroupMessage[256];

int FisherPedigree::Number_of_Pedigrees ()
{
    return _Pedcnts->count();
}

char* FisherPedigree::Allocate_FirstIDs ()
{
    if (_FirstIDs) free (_FirstIDs);
    _FirstIDs = (char*) Calloc (sizeof(char[PERMID_LEN]), 
				     Number_of_Pedigrees());
    return _FirstIDs;
}

char* FisherPedigree::FirstID (int pedindex)
{
    if (_FirstIDs==0 || pedindex >= Number_of_Pedigrees())
    {
	throw Safe_Error_Return ("Invalid FirstID request\n");
    }
    return &_FirstIDs[pedindex * PERMID_LEN];
}
    
// FORTRAN interface

extern "C" void firstid_ (int *fortindex, char* firstid, int firstidlen)
{
    int cindex = *fortindex - 1;
    char *cptr = FisherPedigree::FirstID (cindex);
    int stlen = PERMID_LEN;
    if (stlen > firstidlen) stlen = firstidlen;
    int i = 0;
    for (; i < stlen; i++)
    {
	firstid[i] = cptr[i];
    }
    while (i < firstidlen)
    {
	firstid[i++] = ' ';
    }
}
    

void FisherPedigree::Changing_Pedigree ()
{
    if (Ibd_Index)
    {
	delete Ibd_Index;
	Ibd_Index = 0;
    }
    Ibd_Famid = false;
    _Pedcnts->zero ();
    _Highest_Pedcnt = 0;
    _HouseholdGroup->zero ();
    _HouseholdGroupCnts->zero ();
}

/*
 * Indexes into Ibd_Index record (not exactly same as in file)
 */
#define IBDID  0
#define FIBDID 1
#define MIBDID 2
#define SEX    3
#define MZTWIN 4
#define PEDNO  5
#define ID     6
#define FAMID  7

char FisherPedigree::_pedfilename[1024] = "";

int FisherPedigree::bind (Tcl_Interp *interp)
{
// Clear out last HouseholdGroup message

    HouseholdGroupMessage[0] = '\0';
/*
 * Get pedigree filename from pedigree.info
 */
    FILE *pinfo = fopen ("pedigree.info", "r");
    if (pinfo)
    {
	if (fgets (_pedfilename, 1024, pinfo))
	{
	    if (_pedfilename[strlen(_pedfilename)-1]=='\n') 
		_pedfilename[strlen(_pedfilename)-1]='\0';
	}
	else
	{
	    strcpy (_pedfilename, "");
	}
	fclose (pinfo);
    }
    Matrix *mat = Matrix::find ("house");
    if (mat && Parameter::find ("C2") && 
	Option::get_int ("MergeHousePeds"))
    {
	_IfHouseholdGroups = true;
    }
    else
    {
	_IfHouseholdGroups = false;
    }
/*
 * load reference info from pedindex.out
 * if it hasn't already been loaded (or if it has been reset)
 */
    const char *errmsg = 0;
    if (Ibd_Index && last_grouping == Option::get_int ("MergeAllPeds") &&
	!_IfHouseholdGroups && last_merge_groups==_IfHouseholdGroups)
    {
	return rewind (interp);
    }
    Changing_Pedigree ();
/*
 * Make arrays of pedigree sizes and family ID's
 * Also find biggest pedigree
 */
    last_merge_groups = _IfHouseholdGroups;
    last_grouping = Option::get_int ("MergeAllPeds");
    if (last_grouping)
    {
	_IfHouseholdGroups = false;  // Don't bother with hg's if merging all
    }


    Ibd_Index = TableFile::open ("pedindex.out", &errmsg);
    if (errmsg)
    {
	if (!strcmp (errmsg, "File not found"))
	{
	    RESULT_LIT ("Must load pedigree first (pedindex.out not found)");
	}
	else
	{
	    RESULT_LIT ("Invalid pedindex.out");
	}
	return TCL_ERROR;
    }
    int ncount;
    const char *fname = 0;
    switch (1)
    {
    default:
	Ibd_Index->short_names (&ncount, &errmsg);
	Ibd_Index->start_setup (&errmsg);
	if (errmsg) break;

	Ibd_Index->setup ("IBDID", &errmsg);
	fname = "IBDID";
	if (errmsg) break;

	Ibd_Index->setup ("FIBDID", &errmsg);
	fname = "FIBDID";
	if (errmsg) break;

	Ibd_Index->setup ("MIBDID", &errmsg);
	fname = "MIBDID";
	if (errmsg) break;

	Ibd_Index->setup ("SEX", &errmsg);
	fname = "SEX";
	if (errmsg) break;

	Ibd_Index->setup ("MZTWIN", &errmsg);
	fname = "MZTWIN";
	if (errmsg) break;

	Ibd_Index->setup ("PEDNO", &errmsg);
	fname = "PEDNO";
	if (errmsg) break;

	Ibd_Index->setup ("ID", &errmsg);
	fname = "ID";
	if (errmsg) break;

	fname = 0;
	if (Ibd_Index->test_name ("FAMID", &errmsg))
	{
	    if (errmsg) break;
	    Ibd_Famid = true;
	    Ibd_Index->setup ("FAMID", &errmsg);
	    fname = "FAMID";
	    if (errmsg) break;
	}
    }
    if (errmsg)
    {
	Changing_Pedigree ();
	if (fname)
	{
	    char buf[256];
	    sprintf (buf, "Missing field %s in pedindex.out", fname);
	    RESULT_BUF (buf);
	}
	else
	{
	    RESULT_LIT ("Error reading pedindex.out");
	}
	return TCL_ERROR;
    }

    int last_pedno = -1;
    int pedcnt = 0;
    _Highest_Pedcnt = 0;
    char *last_famid = 0;
    char **data;
    _TotalIndividuals = 0;
    IntArray *ipedno = new IntArray;
    while (0 != (data = Ibd_Index->get(&errmsg)))
    {
	_TotalIndividuals++;
	int pedno = atoi (data[PEDNO]);
	ipedno->add (pedno - 1);
	if (pedno != last_pedno && last_pedno != -1 && last_grouping == 0)
	{
	    _Pedcnts->add (pedcnt);
	    if (pedcnt > _Highest_Pedcnt)
	    {
		_Highest_Pedcnt = pedcnt;
	    }
	    last_pedno = pedno;
	    pedcnt = 1;
	}
	else
	{
	    last_pedno = pedno;
	    pedcnt++;
	}
    }
    if (last_pedno != -1)
    {
	_Pedcnts->add (pedcnt);
	if (pedcnt > _Highest_Pedcnt)
	{
	    _Highest_Pedcnt = pedcnt;
	}
    }
/*
 * Collapse pedigrees into household groups if household analysis
 *   (HouseholdGroups include one or more pedigrees tied by shared households)
 */
    if (_IfHouseholdGroups)
    {
	if (verbose ("ITERATE"))
	{
	    cfprintf (stderr,
		  "Merging pedigrees sharing household groups...\n","");
	}
	_CurrentHouseholdGroup = -1;
	_Pedcnts->copy_to (_HouseholdGroupCnts);

// pedindex pedno was 1 based, but it was converted to 0 based ipedno
// So every index here is 0 based

	int i,j,k;
	int number_of_household_groups = Number_of_Pedigrees();
	for (i = 0; i < number_of_household_groups; i++)
	{
	    _HouseholdGroup->add (i);
	}
	for (i = 0; i < _TotalIndividuals; i++)
	{
	    for (j = i; j < _TotalIndividuals; j++)
	    {
		if (mat->get (i+1, j+1) != 0.0)
		{
		    int igroup = _HouseholdGroup->index (ipedno->index (i));
		    int jgroup = _HouseholdGroup->index (ipedno->index (j));
		    if (igroup != jgroup)
		    {
			int newtotal;
			number_of_household_groups--;
			if (igroup < jgroup)
			{
			    for (k = jgroup; k < Number_of_Pedigrees(); k++)
			    {
				if (_HouseholdGroup->index(k) == jgroup)
				{
				    _HouseholdGroup->change (k, igroup);
				}
			    }
			    newtotal = _HouseholdGroupCnts->index (igroup) 
				+ _HouseholdGroupCnts->index (jgroup);
			    _HouseholdGroupCnts->change (igroup, newtotal);
			    _HouseholdGroupCnts->change (jgroup, 0);
			} else {
			    for (k = igroup; k < Number_of_Pedigrees(); k++)
			    {
				if (_HouseholdGroup->index(k) == igroup)
				{
				    _HouseholdGroup->change (k, jgroup);
				}
			    }
			    newtotal = _HouseholdGroupCnts->index (igroup) 
				+ _HouseholdGroupCnts->index (jgroup);
			    _HouseholdGroupCnts->change (jgroup, newtotal);
			    _HouseholdGroupCnts->change (igroup, 0);
			}
			if (newtotal > _Highest_Pedcnt)
			{
			    _Highest_Pedcnt = newtotal;
			}
		    }
		}
	    }
	}
	if (number_of_household_groups >= Number_of_Pedigrees())
	{
	    sprintf (HouseholdGroupMessage,
		      "No pedigrees shared households with other pedigrees\n",
		      "");
	    if (verbose ("ITERATE"))
	    {
		cfprintf (stderr, "%s\n", HouseholdGroupMessage);
	    }
	}
	else
	{
	    sprintf (HouseholdGroupMessage, 
		     "%d pedigrees merged into %d pedigree-household groups\n",
		     Number_of_Pedigrees(), number_of_household_groups);
	    if (verbose ("ITERATE"))
	    {
		cfprintf (stderr, "%s\n", HouseholdGroupMessage);
	    }

// Special testing and debugging section for household group merging

	    if (Option::get_int ("HouseGroupShow"))
	    {

// This part duplicates the logic in finding the size of each ped-house group
	    int current_household = -1;
	    int current_count;
	    _HouseholdGroupCnts->rewind ();
	    for (int pindex = 0; pindex < number_of_household_groups; pindex++)
	    {
		current_household++;
		while (0 == (current_count = _HouseholdGroupCnts->next()))
		    current_household++;
		if (current_count < 0)
		{
		    fprintf (stderr, 
          "Please report internal error gather merged household groups\n");
		    exit (-1);
		}
		fprintf (stderr, "Group %d [%d] includes pedigrees:\n",
			 pindex, current_household);
		
// This part finds all the pedigrees and invididuals in them within one
// pedigree household group.  The counts must match.

		int total_count = 0;
		for (int test_ped = 0; test_ped < Number_of_Pedigrees();
		     test_ped++)
		{
		    if (current_household == _HouseholdGroup->index(test_ped))
		    {
			fprintf (stderr, "%d ", test_ped);
			total_count += _Pedcnts->index (test_ped);
		    }
		}
		fprintf (stderr, "\nFound %d,  Expected %d\n", total_count,
			 current_count);
		if (total_count != current_count) exit (-1);

	    }
// End of special debugging session
	    }
	}
    }
    delete ipedno;
    return rewind (interp);
}

int FisherPedigree::rewind (Tcl_Interp *interp)
{
    _Pedcnts->rewind ();
    _HouseholdGroup->rewind ();
    _HouseholdGroupCnts->rewind ();

    const char *errmsg = 0;
// 7.5.0 Now checks if Ibd_Index is nonzero (trivial update)
    if (Ibd_Index)
    {
	Ibd_Index->rewind (&errmsg);
	if (errmsg)
	{
	    char buf[1024];
	    sprintf (buf, "Error reading ibdindex.dat: %s", errmsg);
	    RESULT_BUF (buf);
	    Changing_Pedigree ();
	    return TCL_ERROR;
	}
    }
    return TCL_OK;
}

/*
 * NEXTPEDINFO is called by preped to get the number of individuals in the
 * next pedigree so it can set up arrays prior to calling pinput to
 * actually get the data for the individuals.
 */
extern "C" int nextpedinfo_ (int *nptot)
{
    return FisherPedigree::next_pedinfo (nptot);
}

int FisherPedigree::next_pedinfo (int *nptot)
{
// Returns size of next pedigree, or
// If houshold analysis is in effect, returns size of next household group

    int next_nptot;

// Get size of current pedigree (or household group)
    if (_IfHouseholdGroups)
    {
// Skip over 0 length groups...members moved to other group
	_CurrentHouseholdGroup++;
	while (0 == (next_nptot = _HouseholdGroupCnts->next()))
	    _CurrentHouseholdGroup++;
// If next_nptot is INT_MIN, this is end of pedigrees (see below)
                     
// Set up to "select" members of current group from pedigree file
	_RemainingInCurrentPedigree = 0;  // Force check of first pedigree
	_CurrentPedigree = -1;
	_Pedcnts->rewind ();
	const char *errmsg = 0;
	Ibd_Index->rewind (&errmsg);
	if (errmsg) 
	{
	    fprintf (stderr, 
     "Please report internal error selecting household groups:\n%s\n", errmsg);
	    exit (-1);
	}
	    
// Return size of next pedigree
    } else {
	next_nptot = _Pedcnts->next ();
    }
    if (next_nptot == INT_MIN)
    {
	return 0;  // This signifies end of pedigrees
    }
    *nptot = next_nptot;
    return 1;
}

/*
 * NEXTPERSON is called by pinput to get pedigree-related data
 * It also loads the phenotypic data for access through GETPHENO,
 * which is in Class Phenotype.  NEXTPEDINFO must have been called
 * previously (preped does this) to get the total count of individuals
 * in this pedigree (or household-pedigree group).
 */
extern "C" int nextperson_ (double *egoid, double *faid, double *maid,
			     int *sex, double *group, char *permid, 
			     char *famid, int permidlen, int famidlen)
{
    return FisherPedigree::next_person (egoid, faid, maid, sex, group, permid,
					famid);
}

int FisherPedigree::next_person (double *egoid, double *faid, double *maid,
				 int *sex, double *group, char *permid, 
				 char *famid)
{
    const char *errmsg;
    char **data;

//  fprintf (stderr, "Getting next person\n");

// Select next person from current household group
    if (_IfHouseholdGroups)
    {
	if (_RemainingInCurrentPedigree==0)
	{
	    _CurrentPedigree++;
// _CurrentHouseholdGroup was set by next_pedinfo
	    while (_CurrentHouseholdGroup != 
		   _HouseholdGroup->index (_CurrentPedigree))
	    {
		int skip_over = _Pedcnts->next ();
		int i;
		for (i = 0; i < skip_over; i++)
		{
		    Ibd_Index->get (&errmsg);
		}
		_CurrentPedigree++;
	    }
	    _RemainingInCurrentPedigree = _Pedcnts->next ();
	}
	data = Ibd_Index->get (&errmsg);
	_RemainingInCurrentPedigree--;
    }

// Just read person from current pedigree
    else
    {
	data = Ibd_Index->get (&errmsg);
    }

    if (errmsg)
    {
	throw Safe_Error_Return 
	    ("Error getting pedigree data for next_person");
    }
    *egoid = atof (data[IBDID]);
    *faid = atof (data[FIBDID]);
    *maid = atof (data[MIBDID]);
    *sex = atoi (data[SEX]);
    if (!*sex)
    {
	*sex = (data[SEX][0]=='M'||data[SEX][0]=='m') ? 1 : 2;
    }
    Phenotypes::Current_Person_Sex = *sex;
    *group = atof (data[MZTWIN]);
    strncpy (permid, data[ID], PERMID_LEN);
    fortspad (permid, PERMID_LEN);

// Find the corresponding phenotype record and load it
// Also set famid if required

    try
    {
	if (Ibd_Famid && Phenotypes::found_famid())
	{
	    strncpy (famid, data[FAMID], FAMID_LEN);
	    fortspad (famid, FAMID_LEN);
	    Phenotypes::seek (data[ID], data[FAMID]);
	}
	else
	{
	    *famid = '\0';
	    fortspad (famid, FAMID_LEN);
	    Phenotypes::seek (data[ID]);
	}
    } 
    catch (Safe_Error_Return& ser)
    {
	char errmsg[1024];
	sprintf (errmsg, "Error reading phenotype record for ID: %s\n%s",
		 data[ID], ser.message() );
	throw Safe_Error_Return (errmsg);
    }
    return 0;
}


	
