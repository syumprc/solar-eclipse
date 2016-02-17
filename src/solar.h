/*
 * solar.h global declarations for solar
 *
 * Written by Charles Peterson beginning on October 27, 1997
 * Copyright (c) 1997 Southwest Foundation for Biomedical Research
 */

#ifndef SOLAR_H
#define SOLAR_H

// #define MATRIX_DEBUG 1

#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include <stdio.h>
#include <limits.h>

// Apple used to mess this up
// #ifdef __APPLE__
// #define isnan(D) __isnand((D))
// #endif

#define MAX_TRAITS 400
#define MAX_H2QINDEX 100
#define MAX_TABLEFILES 999

// Make sure Tcl_CmdProc is declared with "C" not "C++" type
extern "C" {
#include "tcl.h"
	   };

#include "expression.h"
#include "safelib.h"
#include "tablefile.h"

const int MAX_PHENOTYPE_FILES = 1000;
#define BIGK 40000
const int MAX_PARAMETERS = BIGK;
const int MAX_PHENOTYPES = BIGK;
const int MAX_PHEN_FILE_CHARS = 20*BIGK;


/* VNAMELEN is the maximimum phenotype name length */
/*  Needs to be preset because it defines interfaces between fortran and C */
/*  Used in ccsearch and zscore now (2010) */

#define VNAMELEN 96


/*
 * RESULT_LIT and RESULT_BUF make it easier to set the "interp->result"
 * now that we can no longer access it directly.  Use RESULT_LIT for string
 * literals (constant) and RESULT_BUF for local volatile buffers.  Interpreter
 * named "interp" is assumed.  Many messages fit into one line with these
 * macros instead of sprawling over 3 lines.
 */
#define RESULT_LIT(lit) Tcl_SetResult (interp, (char*)(lit), TCL_STATIC)
#define RESULT_BUF(buf) Tcl_SetResult (interp, (char*)(buf), TCL_VOLATILE)

/*
 * Create Solar_ equivalents to the Tcl functions which perform the required
 * casting to satisfy Workshop 5 compiler (better may be required later)
 */
#define Solar_CreateCommand(i,n,p,d,dp) \
     Tcl_CreateCommand ((i), (char*) (n), (p), (d), (dp))

#define Solar_EvalFile(i,fn) Tcl_EvalFile ((i), (char*) (fn))

#define Solar_AppendElement(i,b) Tcl_AppendElement ((i), (char*) (b))

#define Solar_AppendResult1(i,s1,s0) Tcl_AppendResult ((i), (char*) (s1), (s0))

#define Solar_AppendResult2(i,s1,s2,s0) \
    Tcl_AppendResult ((i), (char*)(s1), (char*)(s2), (s0))

#define Solar_AppendResult3(i,s1,s2,s3,s0) \
    Tcl_AppendResult ((i), (char*)(s1), (char*)(s2), (char*)(s3), (s0))

#define Solar_Eval(i,s) Tcl_Eval ((i), (char*) (s))

/*
 * Other misc stuff
 */

// FORTRAN callback to write on a FORTRAN unit
extern "C" void writestr_ (int *unitno, const char *string, int len);

#define MISSING_PHENOTYPE_STRING "-1.0E20"
const double MISSING_PHENOTYPE = atof (MISSING_PHENOTYPE_STRING);

const double MIN_VISIBLE = 1.0E-12;

const char *SolarVersion ();

void fortspad (char *string, int length);  // Replacement for buggy Spadf

extern "C" void error (const char *message);

// Exception classes

class Safe_Error_Return
{
    char *_message;
public:
    char *message () {return _message;}
    Safe_Error_Return (const char *message) {_message = Strdup (message);}
    Safe_Error_Return (const Safe_Error_Return& e) {
	_message = strdup (e._message);}
    ~Safe_Error_Return () {free(_message);}
};

class Duplicate_Parameter_Name {};

class Missing_Data {};

class Covariate_Parameter {};

// SolarFile adds field mapping to TableFile
//   standard encapsulation of major data formats supported by SOLAR:
//     PEDSYS, Comma Delimited, and Fisher 2 format
//   data is transparent to SOLAR regardless of file type

class SolarFile
{
    const char *file_description;
    char error_message[1024];
    const char *_errmsg;

    SolarFile (const char *desc) {file_description = desc; tf=0; _errmsg=0;}
public:
    TableFile *tf;
    ~SolarFile () {if (tf) delete tf;}
    static SolarFile* open (const char *desc, const char *filename, 
			    const char **errmsg);
    virtual int setup (const char *generic_name, const char **errmsg);
    bool test_name (const char *generic_name, const char **errmsg);
    const char* establish_name (const char *name, const char **errmsg);

// Most methods simply query tf object
//   also must check error status
#define checkerr if (0 != (*errmsg = _errmsg)) return
#define checkerr0 if (0 != (*errmsg = _errmsg)) return 0

    const char **names (int *count, const char **errmsg)
	{checkerr0; return tf->names (count, errmsg);}
    const char **short_names (int *count, const char **errmsg)
	{checkerr0; return tf->short_names (count, errmsg);}
    virtual void start_setup (const char **errmsg)
	{checkerr; tf->start_setup (errmsg);}
    virtual char **get (const char **errmsg)
	{checkerr0; return tf->get (errmsg);}
    virtual void rewind(const char **errmsg)
	{checkerr; tf->rewind (errmsg);}
    const char *filename () {return tf->filename();}
    long get_position () {return tf->get_position();}
    virtual void set_position (int pos, const char **errmsg)
	{checkerr; tf->set_position(pos, errmsg);}
    virtual int *widths (int *count, const char **errmsg)
	{checkerr0; return tf->widths (count, errmsg);}
};

// ScratchFile replaces the old scratch file system with dynamic memory

#ifndef MAX_SCRATCH_UNIT
#define MAX_SCRATCH_UNIT 64
#endif

enum ScratchType {Int, Real};

class ScratchFile
{
 private:
    static ScratchFile* ScratchFiles[MAX_SCRATCH_UNIT];  // Also shows in use
    static ScratchFile* Ptrs[MAX_SCRATCH_UNIT];
    static void CheckUnit (int unit) {if (unit < 0 && unit > MAX_SCRATCH_UNIT)
	{fprintf (stderr, "Invalid scratch unit: %d\n", unit); exit (10);}}
    static void CheckActive (int unit) {CheckUnit(unit);
        if (ScratchFiles[unit]==0) {fprintf (stderr, 
	    "Inactive scratch unit %d\n", unit); exit (10);}} // exit

    ScratchFile* next;
    int size;
    int index;
    double* rarray;
    int* iarray;
    
 public:
    ScratchFile () {next=0; rarray=0; iarray=0; index=0;}
    ~ScratchFile () {if (next) delete next; 
                         if (rarray) free (rarray); if (iarray) free (iarray);}
// Note: open is implied by the first write
    static void Write (int unit, int len, ScratchType stype,
		       int* ints, double* reals);
    static void Read (int unit, int len, ScratchType stype,
		      int* ints, double* reals);
    static void Rewind (int unit);
    static void Reset ();
};



// Main classes of solar

enum mandate {mandatory, optional};

class Field
{
    const char *_internal_name;
    const char *_user_name;
    const char *_description;
    mandate _mandate;
    Field *_next;
    static Field *Fields;
    Field (const char *name, const char *description, mandate man) {
	_internal_name=Strdup(name); _user_name=_internal_name;
	_description=description; _mandate=man; _next=0; _changed=false;}
    static Field* Find (const char *int_name);
    bool _changed;
public:
    Field *add (const char *name, const char *description, mandate man);
    static Field *Fields_Setup ();
    static const char *Change (const char *int_name, const char *user_name);
    static const char *Map (const char *int_name);
    static Field *All () {return Fields;}
    static int Info (Field **f, const char **name, const char **user_name, 
		 const char **description) 
	{if (!*f) return 0; *name=(*f)->_internal_name; 
	*user_name=(*f)->_user_name; *description=(*f)->_description; 
	*f=(*f)->_next; return 1;}
    static bool Changed (const char *int_name);
    static void Start ();
};

class Scale
{
 private:
    char* _name;
    double _adj;
    static Scale* Scales[MAX_PARAMETERS];
    static int Total_Count;
    static bool _all;
 public:
    Scale (const char* name, double adj) {_name=Strdup(name); _adj=adj;}
    ~Scale () {free (_name);}
    static void Write_Commands (FILE *fptr);
    static Scale* Find (const char *testname);
    static void Remove (const char *testname);
    static void Reset ();
    static void All();
    Scale* add ();
    double adj () {return _adj;}
    double adj (double newadj) {return _adj=newadj;}
};

class Parameter;

class Constraint;

class Definition
{
    char* _name;
    char* _definition;

    static Definition** Definitions;
    static int End;
    static int find_index (const char* name);
public:
    ~Definition ();
    const char* definition () {return _definition;}

    static int Delete (const char* name, Tcl_Interp* interp);
    static void Delete_All();
    static int Show_Names (Tcl_Interp* interp);
    static void Reset () {Delete_All();}
    static void Write_Commands (FILE *fptr);
    static Definition* Find (const char* name);
    static int Show (const char* name, Tcl_Interp* interp);
    static int Show_All (Tcl_Interp* interp);
    static int Rename (const char* oldname, const char* newname, Tcl_Interp* interp);
    static int Assign (const char* name, const char* expression, 
			   Tcl_Interp* interp);
};

/*
 * CovariateTerm represents one term in a covariate
 *   It is the building block of Covariate and CovariateList
 */
class CovariateTerm
{
private:
    CovariateTerm** end() {if (!next) return &next; return next->end();}
    const char* fullname_imp (char* buf, char* traitname);
public:
    Expression* expression;
    char* name;
    char* trait;
    int exponent;
    const char* fullname (char* buf) {return fullname_imp (buf, 0);}
;
    bool sex;
    int packed_number;
    double adj_value;
    bool adj_value_set;

    CovariateTerm *next;
    
    CovariateTerm (const char* vname, int exp, const char* qualifier);

    ~CovariateTerm ();

    void add (const char* newname, int exp, const char* qualifier);

    static bool compare (CovariateTerm *term1, CovariateTerm *term2) {
	if (Strcmp (term1->name, term2->name) || 
	    term1->exponent != term2->exponent ||
	    term1->trait && !term2->trait ||
	    Strcmp (term1->trait, term2->trait) ||
	    term1->next==0 && term2->next!= 0 ||
	    term1->next!=0 && term2->next==0) return false;
	if (term1->next) return compare(term1->next, term2->next);
	return true;
    }
};
	
class CovariateList;

class Covariate
{
    static Context* _Context;

    CovariateTerm* _terms;    // Each term has variable and exponent
    friend class Constraint;
    bool _active;             // active or inactive
    int _number;
    char* _mu_string;         // string with terms for mu equation
    Parameter* _beta[MAX_TRAITS];				     
    double _save_beta[MAX_TRAITS];       // Save beta value during suspension

    static int total_count;
    static Covariate* Covariates[];
    Covariate (CovariateTerm* terms);
    ~Covariate();
    static double _minrange; // minimum variable range
    static Covariate *find (CovariateTerm *terms);
    static int last_number_of_traits;
public:
    static void write_commands (FILE *file);
    static void add (CovariateList *cl);
    static void suspend (CovariateList *cl);
    static void remove (CovariateList *cl);
    static void restore (CovariateList *cl);
    static void zero_all ();
    static void new_traits ();

    int active () {return _active;}
    Parameter *beta (int t) {return _beta[t];}
    const char *fullname (char *buf) {return _terms->fullname (buf);}
    static void show_betanames (Tcl_Interp *interp);
    static void show_applicable (Tcl_Interp *interp, bool only_active);

    static Covariate *index (int i) 
	{return (i<total_count) ? Covariates[i] : 0;}
    static void show_all (Tcl_Interp *interp);
    static int bind (Tcl_Interp *interp);
    static int savemeans (int *nvar, double *vmean, double *vmin, 
			   double *vmax, int *output_unit);
    static double eval (int *malep, double *var, double *par, int trait);
    static int set_boundaries (int *PNUM, double *PAR, double *parmin,
			       double *parmax,
			       double *vvar, double *vmin, double *vmax);
    const char *mu_string (int trait);
    static void reset ();
    CovariateTerm *terms () {return _terms;}
    static int boundary_check (double factor, char *pnames, int plen);
};

class StringArray
{
    char **_array;
    int _count;
    int _size;
    int _position;
public:
    StringArray() {_position=0; _count=0; 
                _array=(char **) Calloc (sizeof(const char*), (_size=1));}
    ~StringArray() {for (int i=0; i<_count; i++) if (_array[i]) 
	free (_array[i]); free (_array);}
    void zero () {for (int i=0; i<_count; i++) if (_array[i]) 
	free (_array[i]); _count=0;}
    void add (const char* string);
// {if (++_count>_size)
//	_array=(char**)Realloc(_array,sizeof(char*)*(++_size));
//        _array[_count-1]=Strdup(string);}
    char* index(int i) {return (i<_count) ? _array[i] : 0;}
    int find (const char *target) {for (int i=0; i<_count; i++)
	if (!Strcmp (target, _array[i])) return i; return -1;}
    void rewind() {_position=0;}
    const char *next() 
	{return ((_position)<_count) ? _array[_position++] : "";}
};

class FileIndex;
class IntArray;

#define PERMID_LEN 18  // Must match definition in PINPUT.F
#define FAMID_LEN 18

class FisherPedigree
{
    static TableFile *Ibd_Index;
    static bool Ibd_Famid;
    static int _Highest_Pedcnt;
    static IntArray *_Pedcnts;
    static char _pedfilename[];
    static char* _FirstIDs;
    static int last_grouping;
    static int last_merge_groups;
    static int _TotalIndividuals;

// Household Groups are combinations of pedigree(s) including all members
// of common households
    static bool _IfHouseholdGroups;         // If household group analysis
    static IntArray* _HouseholdGroup;       // Maps ped to h. group
    static IntArray* _HouseholdGroupCnts;   // Size of h. group
    static int _CurrentHouseholdGroup;
    static int _RemainingInCurrentPedigree;
    static int _CurrentPedigree;
public:
    static int rewind (Tcl_Interp *interp);
    static int trait_index;
    static int bind (Tcl_Interp *interp);
    static int Number_of_Pedigrees ();
    static int close ();
    static int Highest_Pedcnt () {return _Highest_Pedcnt;}
    static int next_pedinfo (int *pedcnt);
    static int next_person (double *egoid, double *faid, double *maid,
		           int *sex, double *group, char *permid, char *famid);
    static void Changing_Pedigree ();
    static const char *pedfilename() {return _pedfilename;}
    static char* Allocate_FirstIDs ();
    static char* FirstID (int pedindex);
    static char HouseholdGroupMessage[256];
};

enum ResetMode {soft, hard};

class FilePosition
{
    int file_position[MAX_PHENOTYPE_FILES];
};


enum famid_status {not_known, not_found, found};

class Phenotypes
{
    static int Filecount;
    static const char **PhenoNames[MAX_PHENOTYPE_FILES];
    static int PhenoCount[MAX_PHENOTYPE_FILES];
    static int Phenotypes_In_File[MAX_PHENOTYPE_FILES];
    static time_t Last_Modified[MAX_PHENOTYPE_FILES];
    static Expression* Expressions[MAX_PHENOTYPES];
    
// These variables are used to dereference data and copy to Fortran
// For example, see setup and get_indexed_phenotypes

    static SolarFile *Sfile[MAX_PHENOTYPE_FILES];
    static char **Current_Person_Phenotypes[MAX_PHENOTYPE_FILES];
    static int File_Useage[MAX_PHENOTYPE_FILES];
    static int File_Index[MAX_PHENOTYPES];
    static int File_Pos[MAX_PHENOTYPES];
    static int Setup_Count;
    static StringArray *Setup_Names;
    static int _proband_pindex;
    static void position (int fpos);
    static int build_index (Tcl_Interp *interp);
public:
    static const char* get_var_name (int index);
    static char* maxphen_index (int index);
    static famid_status _found_famid;
    static void reset (ResetMode rm = hard);
    static void load_files (char** filenames, int nfiles, 
			   const char **errmsg);
    static void reload_files (char** filenames, int nfiles, 
			     const char **errmsg);
    static int describe (Tcl_Interp* interp, bool showall);
    static void write_commands (FILE *file);
    static bool available (const char *name);
    static int bind (Tcl_Interp *interp);
    static const char* get_indexed_phenotype (int pindex);
    static void seek (const char *id, const char *famid="");
    static void start_setup ();
    static void setup (const char *name);
    static void setup_proband (const char *name);
//  setup_virtual_position reserves VAR slot for "virtual" phenotypes
//  not in phenotype file such as expression trait and IBDID
    static void setup_virtual_position (const char* name);
    static int find_virtual_position (const char *name);  
    static bool found_famid () {return (_found_famid==found);}
    static const char *test_proband ();
    static const char *filenames ();
    static void start ();
    static int showfiles (Tcl_Interp *interp);
    static Tcl_Interp* Interpp;
// vectors for inverse normal, also used by Trait
    static double* INData;
    static char**  INNames;
    static int     INCount;
    static bool*   INIfclass;
    static int*    INClass;

    static void Initialize_Expressions();
    static void Setup_Expression (Expression *exp, int nvar);
    static void Covar_Of_Trait (int covindex, int traitindex);
    static int Current_Person_Sex;
};

double eval_phenotype (int eqindex);
double eval_sex (int eqindex);


class Zscore
{
public:
    char* pname;
    double mean;
    double sd;
    bool meanvalid;
    int vindex;  //index from EqVar to underlying phenotype

    static Zscore* Zscores[MAX_PARAMETERS];
    static int ZscoreCnt;
    static int Setup (char* pname, int vindex);
    static int SetupMean (char* pname, double mean, double sd, int vindex);
    static void Reset ();
    static Zscore* Get(int index);
    static bool Zneeded (char* phenname);
    static bool If_Ready ();
    int enable (bool validity);
};
double get_zscore (int vindex);

const char loglike_format[] = "%-9.6f";

class Parameter
{
    friend int Covariate::savemeans (int *nvar, double *vmean, double *vmin, 
			              double *vmax, int *output_unit);
    Parameter (const char *name);
    static Parameter *Parameters[];
    static int _count;
    int _number;
    ~Parameter ();
    char *_name;

public:
    static int return_all (Tcl_Interp* interp);
    const char *name () {return (const char*) _name;}
    static void return_names (Tcl_Interp *interp);
    static Parameter *add (const char *name);  // public "new" interface
    static Parameter *insert (const char *name);  
    void remove (); // public delete interface throws exception on error
    static Parameter *find (const char *name);
    static int Find0 (const char *name);
    static void display_all ();
    static void return_commands (Tcl_Interp *interp);
    char *to_string (char *buf);
    static void write_commands (FILE *file);

    void rename (const char *newname);

    double last;
    double start;
    double upper;
    double lower;
    double se;
    double score;
    Covariate *covariate;
    int trait;
    bool mod_start;
    bool mod_lower;
    bool mod_upper;
    bool fix_lower;
    bool fix_upper;

    void zero () {last=0.0;start=0.0;upper=0.0;lower=0.0;se=0.0;score=0.0;}
    int number () {return _number;}
    static int longest_name ();
    static Parameter *index (int i) 
	{return (i< _count) ? Parameters[i] : NULL;}
    static int count () {return _count;}
    static void set (int i, double start, double lower, double upper, 
		     double se, double score)
	{
	    if (0>i || i>=_count) return;
	    Parameters[i]->start = start;
	    Parameters[i]->lower = lower;
	    Parameters[i]->upper = upper;
	    Parameters[i]->se = se;
	    Parameters[i]->score = score;
	}
    static int bind (Tcl_Interp *interp);
    static void reset();
};

enum trait_type {unknown, quantitative, discrete};

class Mu;
class Omega;



class Trait
{
    Expression* _expression;
    char *_name;
    char *_stdev; /* Trait standard deviation (optional) */
    trait_type _type;
    int _index;

    static Context* _Context;
    static Trait* Traits[];
    static int Number_Of_Traits;
    static int _Number_Of_Expressions;
    static int RParameter (int index, const char* name, bool must_have);
public:
    Trait (char* name, char* sdev = 0);
    Trait (Expression* expression);
    ~Trait ();
    static int Number_Of () {return Number_Of_Traits;}
    static const char* Name (int i);
    static const char *Name_Stdev (int i);

// Return zero-based index to trait-specific parameter
//   return -1 if parameter not found

    static int Find_Parameter (int index, const char* name);

    static int One (char* trait_name, char* stdev_name);
    static int Add (char* trait_name, char* stdev_name);
    static int One (Expression* expression);
    static int Add (Expression* expression);
    static void Update_Parameters (Tcl_Interp* interp, char* oldtrait1, 
				   char* oldtrait2);
    static void Set_Discrete (int i);
    static void Set_Quantitative (int i);
    static trait_type Type (int i);
    static trait_type Maximization_Type ();
    static char *Describe_All ();  // Must free result
    static void Write_Commands (FILE *file);
    static int Bind (Tcl_Interp *t);
    static void Reset ();
    static bool Valid;
    static Expression* expression (int i);
    static int Number_Of_Expressions () {return _Number_Of_Expressions;}
};

class Omega
{
    static Expression *_expression;
    static Expression *first_expression;
    static Context *context;
public:
    static Expression *expression () {return _expression;}
    static void new_expression (const char *expr_string);
    static void write_commands (FILE *file);
    static char* free_string() 
	{return Strdup (Omega::expression()->string ());}
    static double eval () {return expression()->eval ();}
    static int bind (Tcl_Interp *interp);
    static void reset ();
    static int check_if_defined (Tcl_Interp* interp);
};

enum model_status {unoptimized, valid, failed};

class Loglike
{
    static model_status _status;
    static double _loglike;
    static double _quadratic;
    static int max_over_this_phen (Tcl_Interp *interp, const char *outfilename);
    static int max_inside_outdir (Tcl_Interp *interp, const char *outfilename);
public:
    static int status () {return (_status==valid)? 1 : 0;}
    static void check_status ();
    static char *get_string (char *buf);
    static char *get_quadratic_string (char *buf);
    static double loglike (double nl) {_status=valid; return _loglike=nl;}
    static double loglike () {check_status(); return _loglike;}
    static double quadratic (double q) {return _quadratic=q;}
    static double quadratic () {return _quadratic;}
    static int maximize (Tcl_Interp *interp, const char *outfilename);
    static void reset () {_status=unoptimized; _loglike=0.0; _quadratic=0.0;}
    friend int ccfisher (Tcl_Interp *interp);
};


// Needed for callbacks upstream of maximize
Tcl_Interp* get_maximize_interp ();



class Model
{
    static bool _loading;
public:
    static int write (FILE *file);
    static int load (const char *filename, Tcl_Interp *interp);
    static int renew (Tcl_Interp *interp);
    static void reset ();
    static bool loading () {return _loading;}
};

/*
 *Generic fast dynamic array with automatic expansion
 *   but no automatic element creation or deletion (takes time)
 *   (Not for classes requiring constructors or destructors)
 *   (Good for built-in types: int, float, ...)
 *
 * 0-based indexing, however, you can ignore that if you use only set and get
 */

class Out_of_Bounds {};
template<class A> class FD_Array
{
    A* _array;
    int _count;
    int _size;
    int _incr_size;
    int _starting_size;
    int _position;
    void expand(int newsize) {if (newsize>_size){_size=newsize+_incr_size;
    _array=(A*) Realloc (_array, sizeof(A)*_size);}}
    void allocate () {_size=_starting_size;
        _array=(A*) Calloc (sizeof(A), _size);}
public:
    FD_Array<A> () {FD_Array<A>(1,1);}
    FD_Array<A> (int starting_size, int incr_size) { _starting_size=starting_size;
        _incr_size=incr_size; zero(); allocate();}
    void renew() {zero(); free(_array); allocate();}
    void zero() {_count=0; _position=0;}
    ~FD_Array() {free (_array);}
    void add (A ele) {expand(_count);_array[_count++]=ele;}
    void set (int newindex, A ele) {expand(newindex+1);_array[newindex]=ele;
        if (_count<=newindex) _count=newindex+1;}
    int count () {return _count;}
    bool last () {return (_position<_count) ? 0 : 1;}
    A operator[] (int i) {if (i>=_count) throw Out_of_Bounds();
        return _array[i];}
    void rewind() {_position=0;}
    int next() {if (_position>=_count) throw Out_of_Bounds(); 
        return _array[_position++];}
};

class PedMatrix
{
    friend class Matrix;
    int start;
    float *values;
    PedMatrix () {start=1;values=0;}
    ~PedMatrix () {delete [] values;}
};


const int MAX_MATRICES = 100;

class Matrix
{
    friend int Omega::bind (Tcl_Interp *interp);
    friend double user_matrix (int key);
    friend double user_second_matrix (int key);
    char *_name;
    char *filename;
    bool _ibd;
    bool _d7;
    Matrix *second_matrix;
    Matrix *first_matrix;  // backpointer used only by second_matrix
    float min;  // min and max are for descriptive purposes
    float max;
    double sum;
    bool defaultable;
    Matrix *default_matrix;
    double *default_scalar;
    bool ids_within_peds;
    static bool Famid_Present;
    static bool Famid_Needed;
    PedMatrix *pedmat;
    int pedmat_count;
    int highest_id;
    void add ();
    void remove ();
    const char *load (const char *specified_filename=0);
    int set (int id1, int id2, float value);
    int* ibdid_found;
    int load_option;

    Matrix (const char *name);
    static Matrix *Matrices[MAX_MATRICES];
    static int count;
    static FD_Array<int> Pedno;   // ID -> Pedno
    static FD_Array<int> Pedsize; // Pedn -> size
    static int last_pedno;
    static const char* load_pedigree ();
    static bool Pedindex_Current;
    static int Pedindex_Highest_Ibdid;
    static int Missing_Ibdid;
    int get_index (int id1, int id2, int base) {
	return (id1 >= id2) ? ((id1+1-base)*(id1-base)/2)+(id2-base):
	                      ((id2+1-base)*(id2-base)/2)+(id1-base);}
    int check_ibdid (int ibdid);

public:    
    float get (int id1, int id2);
    static const char* setup (int option, const char *filename,
			      const char *name1, const char *name2=0);
    static void write_commands (FILE *file);
    static Matrix *index (int i) {return (i<count) ? Matrices[i] : 0;}
    const char *name () {return _name;}

    ~Matrix ();

    static Matrix *find (const char *varname);
    char *command (char *buf);
    char *describe (char *buf);
    static char *commands (char *buf);
    static char *describe_all (char *buf);
    static int return_all (Tcl_Interp* interp);
    static void reset();
    static int bind (Tcl_Interp *interp); 
    bool ibd () {return _ibd;}
    bool d7 () {return _d7;}
    static void Changing_Pedigree();
    static int CheckIbdid (int ibdid);
    static int Check_Matrices();
};

class Term
{
    double _factor;
    char *_name;
    Parameter *_parameter;
public:
    Term (double f, char *n) {_factor=f; _name=n?Strdup (n):0; _parameter=0;}
    ~Term () {if (_name) free (_name);}
    bool cmp (Term *t)
	{return (t->_factor==_factor && 
		 !Strcmp(t->_name,_name)) ? true : false;}
    const char *name () {return _name;}
    double factor() {return _factor;}
    Parameter *parameter() {return _parameter;}
    void add_to (double addend) {_factor += addend;}
    int bind (Tcl_Interp *interp);
    void negate () {_factor = - _factor;}
    static Term* scan (const char **string, const char **error_message);
};

const int MAX_TERMS = 999;
const int MAX_CONSTRAINTS = 999;

class Constraint
{
    static Constraint *Constraints[MAX_CONSTRAINTS];
    static int _count;
    int Nlterms;
    Term *lterms[MAX_TERMS];
    Term *rterm;
    char *string;
    Constraint (const char *str) {Nlterms = 0; rterm=0; string=Strdup(str);
                                  transposed=false;}
    int transpose ();
    bool transposed;
    const char* add ();
    int find_matching_left_side ();
    static Constraint* parse (const char *string, const char** message);
    char* principal();
public:
    ~Constraint ();
    static const char* parse_and_add 
	(const char* string);
    static void write_commands (FILE *file);
    static void display_all ();
    static int return_commands (Tcl_Interp *interp);
    static const char* remove (int number);
    static const char* remove_spec (const char *cstr);
    static int count ();
    static void apply (double *cnstr, double *cvalue);
    static int bind (Tcl_Interp *interp);
    static void reset();
    static int check_parameter (int pnum, double *START, double *LOWER,
				double *UPPER);
    static bool Check_by_name (const char *name);
    static void remove_if_needs_par (const char* pname);
    int gnu_noop () {return 0;}  // obviate stupid warning from gcc
    static void Not_Satisfied (int i);
};

const int MAX_OPTIONS = 256;
const int MAX_OPTION_NAME_LEN = 82;
class Option
{
    char *name;
    char *value;
    char *initial_value;
    int visible;
    int writeable;
    bool list_type;
    char** list;
    int list_len;
    bool session;

    Option (const char *name, const char *value);
    static void add (const char *name, const char *value);
    static void addi (const char *name, const char *value);  // invisible
    static void addnw (const char *name, const char *value);  // not writeable
    static void addl (const char *name, const char *value);  // list type
    static void addses (const char *name, const char *value);  // session
    static int count;
    static Option* Options[MAX_OPTIONS];
    static Option* index (int i) {return (i < count) ? Options[i] : 0;}
    void append (const char* value);
    void free_list ();
public:
    ~Option ();
    static void write_commands (FILE *file);
    static int return_commands (Tcl_Interp* interp);
    static void setup ();
    static int query (const char *name, Tcl_Interp *interp);
    static int get_int (const char *name);
    static double get_double (const char *name);
    static const char* get_string (const char *name);
    static const char* get_element (const char *name, int index);
    static int set (const char *name, const char *value, bool append,
		    Tcl_Interp *interp);
    static void show_all ();
    static void reset ();
    int gnu_noop () {return 0;}  // obviate stupid warning from gcc
    static bool check_list (const char* name, int ivalue);
};

// EqVar is used to ensure variables used in Mu and Omega equations
// and trait expression (if used)
// are included in analysis, and to get their indexes during evaluation

class EqVar
{
    static char *name[2*MAX_PHENOTYPES];
    static int packed_index[2*MAX_PHENOTYPES];
    static int count;
    static int _For_Trait[2*MAX_PHENOTYPES];
 public:
    static int bind (Tcl_Interp *interp); // Bind before Trait, Mu and Omega!
    static int add (const char *vname);  // Returns unpacked index
    static int add_for_trait (const char *vname, int trait);
    static int for_trait (int i);
    static const char* index (int i);
    static void set_packed_index (int i, int packed_index);
    static int get_packed_index (int i)
	{if (i>=count); return packed_index[i];}
};

class Mu
{
    static Expression* expression;  // other than default part
    static bool _default_part;
    static char *before_string;
    static char *after_string;
    static char *combined_string;
    static Context *context;
    static char *comment_string();
public:
    static void write_commands (FILE *file);
    static double *vmeans;  // Set by ccsearch, read by get_variable_means
    static double *vmins;
    static double *vmaxs;				   

    static bool default_part_used () {return _default_part;}
    static char* free_this_string ();
    static void new_expression (const char* string);
    static int bind (Tcl_Interp *interp);
    static double eval () {return (expression) ? expression->eval() : 0.;}
    static void add_term_before (const char *term_plus); // must include trailing '+'
    static void add_term_after (const char *plus_term);  // must include leading +/-
    static void reset ();
};

char *append_extension (const char *given_filename, const char *extension);
// extension must include "." or other delimiter

class EVD
{
public:
    static void Reset();
    static void Flush();
    static void Make_Matrix (int zindex, int* maxpeo, int* n, 
			     double* mu, double* vvar, int* nvar,
			     int* ierr);
};

int verbose (const char *keyword);

enum verbosity_level {min_verbosity=0, default_verbosity=255, 
		      plus_verbosity=0x3ff, max_verbosity=0xffff};

class Verbosity
{
    static unsigned int _verbosity;
public:
    static bool def () {return (_verbosity == default_verbosity) 
				? true : false;}
    static bool plus () {return (_verbosity == plus_verbosity) 
			    ? true : false;}
    static bool min () {return (_verbosity == min_verbosity) 
			    ? true : false;}
    static bool max () {return (_verbosity == max_verbosity) 
			    ? true : false;}
    static void set (int new_verbosity) {_verbosity = new_verbosity;}
    static int get () {return _verbosity;}
};

#define MAXLOC 10000
#define MMRKNM 20

struct Mrk {
    char name[MMRKNM+1];
    int ntyped;
    int nfoutyp;
};

class Marker
{
    char _filename[1024];
    char error_message[1024];
    SolarFile *Tfile;
    const char **_names;
    int *_widths;
    int _count;
    int _id_len;
    int _famid_len;
    int _gtype_len;
    int _nloci;
    struct Mrk *_mrk[MAXLOC];
public:
    Marker (const char*);
    ~Marker ();
    int load (const char, Tcl_Interp*);
    bool get_stats (const char**);
    int run_makeped (Tcl_Interp*);
    bool unload (bool, Tcl_Interp*);
    int discrep (const char*, Tcl_Interp*);
    void show_names (Tcl_Interp*);
    void show (char*, Tcl_Interp*);
    const char *filename () {return _filename;}
    int id_len () {return _id_len;}
    int famid_len () {return _famid_len;}
    int gtype_len () {return _gtype_len;}
    int nloci () {return _nloci;}
    const char *mrkname (int i) {return _mrk[i]->name;}
    int ntyped (int i) {return _mrk[i]->ntyped;}
    int nfoutyp (int i) {return _mrk[i]->nfoutyp;}
    int get_marker (const char*);
    void delete_Tfile ();
};

#define MAXALL 500
#define MALLNM 20

struct Allele {
    char name[MALLNM+1];
    double freq;
    double mle_se;
};

struct LocusFreq {
    char name[MMRKNM+1];
    char xlinked;
    char mle_status;
    char whence;
    double chi2_hwe;
    int nall;
    struct Allele *all[MAXALL];
};

class Freq
{
    char _filename[1024];
    char error_message[1024];
    int _nloci;
    struct LocusFreq *_locus[MAXLOC];
public:
    Freq (const char *fname) {strcpy(_filename,fname); _nloci=0;}
    ~Freq ();
    int load (char, Tcl_Interp*);
    bool read_info (const char**);
    int write_info (Tcl_Interp*);
    int get_freqs (const char*, char, Tcl_Interp*);
    void write_locfiles (struct LocusFreq*, bool);
    bool unload (Tcl_Interp*);
    int mlefreq (const char*, bool, bool, Tcl_Interp*);
    int save (const char*, Tcl_Interp*);
    const char *show (char*);
    const char *filename () {return _filename;}
    void change_filename (const char *fname) {strcpy(_filename,fname);}
    int nloci () {return _nloci;}
    const char *mrkname (int i) {return _locus[i]->name;}
    char xlinked (int i) {return _locus[i]->xlinked;}
    char mle_status (int i) {return _locus[i]->mle_status;}
    void set_mle_status (int i, char status) {_locus[i]->mle_status=status;}
    char whence (int i) {return _locus[i]->whence;}
    double chi2_hwe (int i) {return _locus[i]->chi2_hwe;}
    int nall (int i) {return _locus[i]->nall;}
    struct Allele *all (int i, int j) {return _locus[i]->all[j];}
    const char *name (int i, int j) {return _locus[i]->all[j]->name;}
    double freq (int i, int j) {return _locus[i]->all[j]->freq;}
    double mle_se (int i, int j) {return _locus[i]->all[j]->mle_se;}
    int get_marker (const char*);
    void remove_marker (int);
    void remove_marker (const char*);
};

struct LocusMap {
    char name[MMRKNM+1];
    double locn;
};

class Map
{
    char _filename[1024];
    char _chrnum[1024];
    char _mapfunc;
    int _nloci;
    struct LocusMap *_locus[MAXLOC];
public:
    Map (const char *fname) {strcpy(_filename,fname); _chrnum[0]=0; _nloci=0;}
    ~Map ();
    int load (Tcl_Interp*, char);
    void show_names (Tcl_Interp*);
    char *show (char*);
    char *filename () {return _filename;}
    char *chrnum () {return _chrnum;}
    char mapfunc () {return _mapfunc;}
    int nloci () {return _nloci;}
    char *mrkname (int i) {return _locus[i]->name;}
    double mrklocn (int i) {return _locus[i]->locn;}
    int get_marker (const char*);
};

struct Ped {
    int nfam;
    int nind;
    int nfou;
    int nlbrk;
    const char inbred;
};

class Pedigree
{
    char _filename[1024];
    char error_message[1024];
    SolarFile *Tfile;
    int *_widths;
    int _count;
    int _id_len;
    int _sex_len;
    int _mztwin_len;
    int _hhid_len;
    int _famid_len;
    int _nped;
    int _nfam;
    int _nind;
    int _nfou;
    struct Ped **_ped;
    Marker *_marker;
public:
    Pedigree (const char*);
    ~Pedigree ();
    int load (bool, Tcl_Interp*);
    char *show (char*);
    const char *filename () {return _filename;}
    int id_len () {return _id_len;}
    int sex_len () {return _sex_len;}
    int mztwin_len () {return _mztwin_len;}
    int hhid_len () {return _hhid_len;}
    int famid_len () {return _famid_len;}
    bool get_stats (const char**);
    int nped () {return _nped;}
    int nfam () {return _nfam;}
    int nfam (int i) {return _ped[i]->nfam;}
    int nind () {return _nind;}
    int nind (int i) {return _ped[i]->nind;}
    int nfou () {return _nfou;}
    int nfou (int i) {return _ped[i]->nfou;}
    int nlbrk (int i) {return _ped[i]->nlbrk;}
    char inbred (int i) {return _ped[i]->inbred;}
    Marker *marker ();
    void add_marker (Marker* m) {_marker=m;}
    void delete_marker ();
    void delete_Tfile ();

// Other classes can add hooks here for when pedigree is being changed
    static void Changing_Pedigree () {
	FisherPedigree::Changing_Pedigree();
	Matrix::Changing_Pedigree();
    }
};

extern Pedigree *currentPed;
extern Freq *currentFreq;
extern Map *currentMap;

extern bool loadedPed (void);
extern bool loadedFreq (void);

extern bool NoMLE;
extern bool XLinked;
extern bool MCarlo;
extern bool MMSibs;

class Voxel {
public:
    static bool _Valid;
    static int X;
    static int Y;
    static int Z;
    static int Show (Tcl_Interp *interp);
    static int Set (Tcl_Interp *interp, char *voxel_spec);
    static int Set (Tcl_Interp* interp, int x, int y, int z);
    static bool Valid () {return _Valid;}
    static void write_commands (FILE *file);
};



#endif
