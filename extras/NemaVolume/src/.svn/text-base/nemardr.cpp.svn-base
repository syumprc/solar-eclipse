/////////////////////////////////////////////////////////////////////
/// UTHSCSA - Research Imaging Center  ///
/////////////////////////////////////////////////////////////////////
//  Module Name: nemardr.C
//  Birth Date: Not sure, maybe 1992.
//  Programmer: Hunter Downs, Dan Nickerson
/////////////////////////////////////////////////////////////////////
/////////////////////////  Description  /////////////////////////////
//
/////////////////////////////////////////////////////////////////////
//////////////////////////  Revisions  //////////////////////////////
// Programmer   , Date     : Explanation  ///////////////////////////
/////////////////////////////////////////////////////////////////////
// Dan Nickerson, 3-25-97  : Added this comment block.
// Dan Nickerson, 3-25-97  : Fixed duplicate function - it didn't copy
//                           the image min and max.
// Bill Rogers 6-11-07 : removed all the minc stuff
// Bill Rogers 3-5-10 : Added wrapper to read and write multiple headers
/////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <ostream>
#include <ios>
#include <string.h>
#define NEMA_MOD
#include <ctype.h>
#include "dtypes.h"
#include "nemardr.h"

using namespace std;

// local function declarations
int xlate(char *, int, char **);
int nema_magic(char *);

int NEMA_get_totv(char **, int);
int NEMA_get_totv(char *bptr[], int n, int *keynum); // br
int NEMA_get_vol(char **, int);
int NEMA_get_vol(char *bptr[], int n, int *keynum); // br
int NEMA_get_tots(char **, int);
void NEMA_fill_slices(char **, int, char **, double *, long *);
int NEMA_get_type(char **, int);
int NEMA_get_bita(char **, int);
int NEMA_get_bits(char **, int);
int NEMA_get_bith(char **, int);
double NEMA_get_scale(char *bptr[], int n);
double NEMA_get_scale_offset(char *bptr[], int n);
int NEMA_get_rows(char **, int);
int NEMA_get_cols(char **, int);
int NEMA_check_slices(char **, int);
double NEMA_get_rowv(char **, int);
double NEMA_get_colv(char **, int);
double NEMA_get_slicev(char **, int);
double NEMA_get_value(char **, int, const char *which);
void NEMA_get_name(char **, int, char *);
void NEMA_get_number(char **, int, char *);
void NEMA_get_date(char **, int, char *);
void NEMA_get_unit(char **, int, char *);
void NEMA_get_xyorigin(char **, int, double *, double *);
void NEMA_get_orientation(char **bptr, int num_keys, char *orientation);
void NEMA_get_image_minmax(char **, int n, double *, double *, int*, int*);
void NEMA_get_xyzoffset(char **, int n, XYZPOINT *);
void NEMA_fill_slices2(char **, int, double *, double *, double *, int*, int*,
		int*);
void NEMA_get_version(char **bptr, int num_keys, char *sn_version,
		char *bp_version);

/*
 #define DEB	1
 */

/*!
 * Function to write the beginning of the header for a NEMA file that
 * may include multiple volumes. - br
 * @param fname - filename for NEMA .des file
 * @param nvol - number of volumes in file
 * @return - 1 on success, 0 on file error
 */
int NEMA_Start_Header(char *fname, int nvol)
{
	// Open the file
	fstream out;
	out.open(fname,ios::out);
	if (!out)
	{
		return 0;
	}

	out << NEMA_MAGIC << "\n";
	out << TOTAL_VOLUMES << "=" << nvol << "\n";
	out.close();
	return 1;
}

/*!
 * Function to read in entire header for a NEMA file that may include
 * multiple volumes. - br
 * @param NemaVol - pointer to array of NEMA headers
 * @param filename - name of NEMA file to open
 * @return - number of volumes in file
 */
int NEMA_READ(NEMAFILE **NemaVol, char *fname)
{
	ifstream infile;
	long int len; // TLA changed - previous nolonger worked
	char *buffer;
	char **bptr;
	int num_keys;
	// Verify that the filename is not too long
	int parsed;
	if ((strlen(fname)) > FILENAME_SIZE)
	{
		parsed = NEMA_FILENAME_TOO_LONG;
		return 0;
	}

	infile.open(fname, ios::in);
	if (!infile)
	{
		parsed = NEMA_FILE_DID_NOT_OPEN;
		return 0;
	}

	// Find out how long the file is
	infile.seekg(0, ios::end);
	len = infile.tellg();

	// Allocate a buffer large enough to hold the file
	buffer = new char[len];
	if (!buffer)
	{
		parsed = COULD_NOT_ALLOC_BUFFER;
		return 0;
	}

	// Read in the whole file
	infile.seekg(0, ios::beg);
	infile.read(buffer, (int) len);
	//  cout << buffer << '\n'; // TLA junk

	// Verify that this is a nema file
	if (!nema_magic(buffer))
	{
		parsed = NOT_NEMA_FILE;
		return 0;
	}
	else
	{
		// Allocate a buffer for parsed pointers
		bptr = new char *[(int) len];

		// Parse the buffer
		num_keys = xlate(buffer, (int) len, bptr);

#ifdef DEB
		std::cerr << num_keys << "\n";
		for (idx = 0; idx < num_keys;idx++)
		std::cerr << idx << ":  " << bptr[idx] << "\n";
#endif
		infile.close();
	}

	// find out how many volumes we have
	int totvolkey;	// key num for total volumes
	int nvol = NEMA_get_totv(bptr, num_keys,&totvolkey);

	// now that we know how many volumes there are, allocate memory
	NEMAFILE *NemaTemp;
	NemaTemp = new NEMAFILE[nvol];
	(*NemaVol) = NemaTemp;

	// read in each volume
	int volkey;	// key number to start reading in each volume
	volkey = totvolkey; 	// start reading just past total volumes
	num_keys -= totvolkey;		// now there are less keys to contend with
	for ( int n=0 ; n<nvol ; ++n)
	{
		// get the start of this volume in the list of key words
		int keyidx;
		NEMA_get_vol(&bptr[volkey], num_keys,&keyidx);
		volkey += (keyidx+1); 	// start reading just past volume keyword
		num_keys -= (keyidx+1);	// now there are less keys to contend with

		// to make our search easier, only read up to the next volume
		int last_key;
		if ( n < (nvol-1) ) // do all but last volume
			NEMA_get_vol(&bptr[volkey], num_keys,&last_key);
		else
			last_key = num_keys;

		// Get the total number of scans (a.k.a slices)
		int tot_scans = NEMA_get_tots(&bptr[volkey], last_key);
		if (!tot_scans)
			return 0;
		NemaTemp[n].init(tot_scans);

		// get the stuff for each individual volume

		// Copy the filename
		strcpy(NemaTemp[n].filename, fname);

		// Get the path from the file name, if no path, set to null.   --jls
		strcpy(NemaTemp[n].directory, fname);
		int i;
		for (i = strlen(NemaTemp[n].directory); ((NemaTemp[n].directory[i] != '/') && (i > 0)); i--)
		{
		}
		if (i > 0)
			NemaTemp[n].directory[i + 1] = '\0';
		else
			NemaTemp[n].directory[0] = '\0';

		// Get the parsed data
		NemaTemp[n].rows = NEMA_get_rows(&bptr[volkey], last_key);
		NemaTemp[n].cols = NEMA_get_cols(&bptr[volkey], last_key);
		NemaTemp[n].imagetype = NEMA_get_type(&bptr[volkey], last_key); // TLA
		NemaTemp[n].bitsall = NEMA_get_bita(&bptr[volkey], last_key);
		NemaTemp[n].bitstor = NEMA_get_bits(&bptr[volkey], last_key);
		NemaTemp[n].bithigh = NEMA_get_bith(&bptr[volkey], last_key);
		NemaTemp[n].image_scale = NEMA_get_scale(&bptr[volkey], last_key);
		NemaTemp[n].image_scale_offset = NEMA_get_scale_offset(&bptr[volkey], last_key);
		if (NemaTemp[n].imagetype == 4 && NemaTemp[n].bitsall == 64)
			NemaTemp[n].imagetype = 5;
		if (NemaTemp[n].imagetype == 5 && NemaTemp[n].bitsall == 32)
			NemaTemp[n].imagetype = 4;
		NemaTemp[n].pixel_bandwidth1 = NEMA_get_value(&bptr[volkey], last_key, PIXEL_BANDWIDTH1);
		NemaTemp[n].pixel_bandwidth2 = NEMA_get_value(&bptr[volkey], last_key, PIXEL_BANDWIDTH2);
		NemaTemp[n].repetition_time1 = (int) NEMA_get_value(&bptr[volkey], last_key,
				REPETITION_TIME1);
		NemaTemp[n].repetition_time2 = (int) NEMA_get_value(&bptr[volkey], last_key,
				REPETITION_TIME2);
		NemaTemp[n].echo1_time = (int) NEMA_get_value(&bptr[volkey], last_key, ECHO1_TIME);
		NemaTemp[n].echo2_time = (int) NEMA_get_value(&bptr[volkey], last_key, ECHO2_TIME);
		NEMA_get_name(&bptr[volkey], last_key, NemaTemp[n].name);
		NEMA_get_number(&bptr[volkey], last_key, NemaTemp[n].patient_number);
		NEMA_get_date(&bptr[volkey], last_key, NemaTemp[n].scan_date);
		NEMA_get_version(&bptr[volkey], last_key, NemaTemp[n].sn_version, NemaTemp[n].bp_version);
		NEMA_get_orientation(&bptr[volkey], last_key, NemaTemp[n].orientation);
		NEMA_get_unit(&bptr[volkey], last_key, NemaTemp[n].spatial_units);
		NemaTemp[n].rowvec = NEMA_get_rowv(&bptr[volkey], last_key);
		NemaTemp[n].colvec = NEMA_get_colv(&bptr[volkey], last_key);
		NemaTemp[n].slicevec = NEMA_get_slicev(&bptr[volkey], last_key);
		NEMA_get_xyorigin(&bptr[volkey], last_key, &NemaTemp[n].XYorigin[0],
				&NemaTemp[n].XYorigin[1]);
		NEMA_get_xyzoffset(&bptr[volkey], last_key, &NemaTemp[n].XYZoffset);
		NEMA_get_image_minmax(&bptr[volkey], last_key, &NemaTemp[n].image_min,
				&NemaTemp[n].image_max,&NemaTemp[n].imin_flag,
				&NemaTemp[n].imax_flag);
		NEMA_fill_slices(&bptr[volkey], last_key, NemaTemp[n].Slice_filename,
				NemaTemp[n].Position, NemaTemp[n].Offset);
		for (int idx = 0; idx < tot_scans; idx++)
			NemaTemp[n].data_scale[idx] = NemaTemp[n].image_scale;
		NEMA_fill_slices2(&bptr[volkey], last_key, NemaTemp[n].slice_min,
				NemaTemp[n].slice_max, NemaTemp[n].data_scale,
				&NemaTemp[n].smin_flag, &NemaTemp[n].smax_flag,
				&NemaTemp[n].dscale_flag);

		NemaTemp[n].load_position(&bptr[volkey], last_key); // TLA from diva
		parsed = PARSED_OK;
	}
	// Free allocated memory
	delete [] bptr;
	delete buffer;

	return nvol;
}

NEMAFILE::NEMAFILE()
{
	Position = NULL;
	Slice_filename = NULL;
	slice_min = NULL;
	slice_max = NULL;
	data_scale = NULL;
	Offset = NULL;
}
/*!
 * This initializes the class - added return value and default number
 * of scans - br
 * @param tot_scans
 * @return
 */
int NEMAFILE::init(int tot_scans)
{
	// Initialize everything
	total_scans = tot_scans;
	if ( total_scans == 0 )
	{
		Position = NULL;
		Slice_filename = NULL;
		slice_min = NULL;
		slice_max = NULL;
		data_scale = NULL;
		Offset = NULL;
	}
	else // allocate for these guys
	{
		// Allocate space for the filenames, offsets, and positions, slice_min, slice_max and data_scale
		Position = new double[tot_scans];
		Slice_filename = new char *[tot_scans];
		slice_min = new double[tot_scans];
		slice_max = new double[tot_scans];
		data_scale = new double[tot_scans];
		Offset = new long[tot_scans];

		if ((!Position) || (!Offset)
				|| (!Slice_filename) || !data_scale
				|| !slice_min || !slice_max)
		{
			parsed = BAD_ALLOC;
			return 0;
		}

		// allocate space for slice file names
		for (int idx = 0; idx < total_scans; idx++)
		{
			if ((Slice_filename[idx] = new char[FILENAME_SIZE]) == NULL)
			{
				parsed = BAD_ALLOC;
				return 0;
			}
		}
		// initialize min, max, and scale values
		for (int idx = 0; idx < total_scans; idx++)
		{
			slice_min[idx] = slice_max[idx] = 0.0;
			data_scale[idx] = 1.0;
		}
	}

	filename[0] = 0;
	parsed = NOT_PARSED;
	name[0] = 0;
	patient_number[0] = 0;
	scan_date[0] = 0;
	rows = 0;
	cols = 0;
	rowvec = 0.0;
	colvec = 0.0;
	slicevec = 0.0;
	spatial_units[0] = 0;
	bitsall = 0;
	bitstor = 0;
	bithigh = 0;
	image_scale = 1.0;
	imagetype = 0; // TLA changed type to imagetype
	image_min = 0;
	image_max = 0;
	XYorigin[0] = 0;
	XYorigin[1] = 0;
	orientation[0] = '\0';
	XYZoffset[0] = 0;
	XYZoffset[1] = 0;
	XYZoffset[2] = 0;
	smin_flag = 0;
	smax_flag = 0;
	imin_flag = 0;
	imax_flag = 0;
	dscale_flag = 0;
	pixel_bandwidth1 = pixel_bandwidth2 = 0.0L;
	echo1_time = echo2_time = 0L;
	repetition_time1 = repetition_time2 = 0L;

	return 1; // initialization ok
}

void NEMAFILE::operator()(char *file)
{
	parseit(file);
}

NEMAFILE::NEMAFILE(char *file)
{

	parseit(file);

	return;
}

NEMAFILE::~NEMAFILE()
{
	for (int idx = 0; idx < total_scans; idx++)
		delete Slice_filename[idx];
	delete[] Slice_filename;
	delete[] Position;
	delete[] Offset;
	delete[] data_scale;
	delete[] slice_min;
	delete[] slice_max;
}

void NEMAFILE::parseit(char *file)
{
	parsed = 0; // TLA added

	// cout << "in parseit - filename = " << file << " \n"; // TLA junke
	ifstream infile;
	// streamoff len ;
	long int len; // TLA changed - previous nolonger worked
	char *buffer;
	char **bptr;
	int num_keys;
	int idx, i;

	// Initialize everything
	this->init();

	//
	//
	// Verify that the filename is not too long
	if ((strlen(file)) > FILENAME_SIZE)
	{
		parsed = NEMA_FILENAME_TOO_LONG;
		return;
	}

	// Open the file
	//cout << "My file name is " << file << '\n'; // TLA test junk

	//	cout << "in parseit - after init before fail  = " << file << " \n"; // TLA junke
	infile.open(file, ios::in);
	if (!infile)
	{
		parsed = NEMA_FILE_DID_NOT_OPEN;
		return;
	}
	// Copy the filename
	strcpy(filename, file);

	// Get the path from the file name, if no path, set to null.   --jls
	strcpy(directory, file);
	for (i = strlen(directory); ((directory[i] != '/') && (i > 0)); i--)
	{
	}
	if (i > 0)
		directory[i + 1] = '\0';
	else
		directory[0] = '\0';

	// Find out how long the file is
	infile.seekg(0, ios::end);
	len = infile.tellg();

	// Allocate a buffer large enough to hold the file
	buffer = new char[len];
	if (!buffer)
	{
		parsed = COULD_NOT_ALLOC_BUFFER;
		return;
	}

	// Read in the whole file
	infile.seekg(0, ios::beg);
	infile.read(buffer, (int) len);
	//  cout << buffer << '\n'; // TLA junk

	// Verify that this is a nema file
	if (!nema_magic(buffer))
	{
		parsed = NOT_NEMA_FILE;
		return;
	}
	else
	{
		// Allocate a buffer for parsed pointers
		bptr = new char *[(int) len];

		// Parse the buffer
		num_keys = xlate(buffer, (int) len, bptr);
#ifdef DEB
		std::cerr << num_keys << "\n";
		for (idx = 0; idx < num_keys;idx++)
		std::cerr << idx << ":  " << bptr[idx] << "\n";
#endif
		// Get the total number of scans (a.k.a slices)
		total_scans = NEMA_get_tots(bptr, num_keys);
		if (!total_scans)
		{
			/*---- Verify that there are no images in this file...compensate for a bug in DIPS ----*/
			if ((total_scans = NEMA_check_slices(bptr, num_keys)) == 0)
			{
				parsed = BAD_TOTAL_SCANS;
				return;
			}
		}

		// Allocate space for the filenames, offsets, and positions, slice_min, slice_max and data_scale
		Slice_filename = new char *[total_scans];
		Position = new double[total_scans];
		Offset = new long[total_scans];
		data_scale = new double[total_scans];
		slice_min = new double[total_scans];
		slice_max = new double[total_scans];

		if ((!Position) || (!Offset) || (!Slice_filename) || !data_scale
				|| !slice_min || !slice_max)
		{
			parsed = BAD_ALLOC;
			return;
		}

		for (idx = 0; idx < total_scans; idx++)
			if ((Slice_filename[idx] = new char[FILENAME_SIZE]) == NULL)
			{
				parsed = BAD_ALLOC;
				return;
			}

		// Get the parsed data
		rows = NEMA_get_rows(bptr, num_keys);
		cols = NEMA_get_cols(bptr, num_keys);
		imagetype = NEMA_get_type(bptr, num_keys); // TLA
		bitsall = NEMA_get_bita(bptr, num_keys);
		bitstor = NEMA_get_bits(bptr, num_keys);
		bithigh = NEMA_get_bith(bptr, num_keys);
		image_scale = NEMA_get_scale(bptr, num_keys);
		if (imagetype == 4 && bitsall == 64)
			imagetype = 5;
		if (imagetype == 5 && bitsall == 32)
			imagetype = 4;
		pixel_bandwidth1 = NEMA_get_value(bptr, num_keys, PIXEL_BANDWIDTH1);
		pixel_bandwidth2 = NEMA_get_value(bptr, num_keys, PIXEL_BANDWIDTH2);
		repetition_time1 = (int) NEMA_get_value(bptr, num_keys,
				REPETITION_TIME1);
		repetition_time2 = (int) NEMA_get_value(bptr, num_keys,
				REPETITION_TIME2);
		echo1_time = (int) NEMA_get_value(bptr, num_keys, ECHO1_TIME);
		echo2_time = (int) NEMA_get_value(bptr, num_keys, ECHO2_TIME);
		NEMA_get_name(bptr, num_keys, name);
		NEMA_get_number(bptr, num_keys, patient_number);
		NEMA_get_date(bptr, num_keys, scan_date);
		NEMA_get_version(bptr, num_keys, sn_version, bp_version);
		NEMA_get_orientation(bptr, num_keys, orientation);
		NEMA_get_unit(bptr, num_keys, spatial_units);
		rowvec = NEMA_get_rowv(bptr, num_keys);
		colvec = NEMA_get_colv(bptr, num_keys);
		slicevec = NEMA_get_slicev(bptr, num_keys);
		NEMA_get_xyorigin(bptr, num_keys, &XYorigin[0], &XYorigin[1]);
		NEMA_get_xyzoffset(bptr, num_keys, &XYZoffset);
		NEMA_get_image_minmax(bptr, num_keys, &image_min, &image_max,
				&imin_flag, &imax_flag);
		NEMA_fill_slices(bptr, num_keys, Slice_filename, Position, Offset);
		for (idx = 0; idx < total_scans; idx++)
			data_scale[idx] = image_scale;
		NEMA_fill_slices2(bptr, num_keys, slice_min, slice_max, data_scale,
				&smin_flag, &smax_flag, &dscale_flag);

		load_position(bptr, num_keys); // TLA from diva
		parsed = PARSED_OK;
	}

	infile.close();
	// Free allocated memory
	delete bptr;
	delete buffer;

}

void NEMAFILE::set_orientation(char *orient)
{
	strcpy(orientation, orient);
}

void NEMAFILE::duplicate(NEMAFILE *orig)
{
	int idx;

	// Initialize everything
	strcpy(filename, orig->filename);
	parsed = orig->parsed;
	strcpy(name, orig->name);
	set_orientation(orig->get_orientation());
	rows = orig->rows;
	cols = orig->cols;
	rowvec = orig->rowvec;
	colvec = orig->colvec;
	slicevec = orig->slicevec;
	pixel_bandwidth1 = orig-> get_pixel_bandwidth2();
	pixel_bandwidth2 = orig-> get_pixel_bandwidth2();
	echo1_time = orig -> get_echo1_time();
	echo2_time = orig -> get_echo2_time();
	repetition_time1 = get_rep_time1();
	repetition_time2 = get_rep_time2();
	strcpy(spatial_units, orig->spatial_units);
	bitsall = orig->bitsall;
	bitstor = orig->bitstor;
	bithigh = orig->bithigh;
	imagetype = orig->imagetype;
	/*   imin_flag = orig->imin_flag;
	 imax_flag = orig->imax_flag;
	 smin_flag = orig->smin_flag;
	 smax_flag = orig->smax_flag;
	 dscale_flag = orig->dscale_flag;
	 */
	set_image_max(orig -> get_image_max());
	set_image_min(orig -> get_image_min());

	// Allocate space for the filenames, offsets, slice_min, slice_max,
	// data_scale, and positions. If we are already fully instantiated
	// then we need to delete the current space.
	for (idx = 0; idx < total_scans; idx++)
		delete Slice_filename[idx];
	//Now copy total_scans.
	total_scans = orig->total_scans;

	delete[] Slice_filename;
	Slice_filename = new char *[total_scans];
	delete[] Position;
	Position = new double[total_scans];
	delete[] Offset;
	Offset = new long[total_scans];
	delete[] data_scale;
	data_scale = new double[total_scans];
	delete[] slice_min;
	slice_min = new double[total_scans];
	delete[] slice_max;
	slice_max = new double[total_scans];

	if ((!Position) || (!Offset) || (!Slice_filename) || (!data_scale)
			|| (!slice_min) || (!slice_max))
	{
		parsed = BAD_ALLOC;
		return;
	}

	for (idx = 0; idx < total_scans; idx++)
	{
		if ((Slice_filename[idx] = new char[FILENAME_SIZE]) == NULL)
		{
			parsed = BAD_ALLOC;
			return;
		}
		strcpy(Slice_filename[idx], orig->Slice_filename[idx]);
		Offset[idx] = orig->Offset[idx];
		Position[idx] = orig->Position[idx];
		data_scale[idx] = orig->data_scale[idx];
		slice_min[idx] = orig->slice_min[idx];
		slice_max[idx] = orig->slice_max[idx];
	}

}

void NEMAFILE::write()
{
	ofstream out;
	int idx, j;
	const char* q = "\"";
	char temp[256];

	// Open the file
	out.open(filename);
	if (!out)
	{
		parsed = NEMA_OUTPUT_FILE_DID_NOT_OPEN;
		return;
	}

	out << NEMA_MAGIC << "\n";
	out << TOTAL_VOLUMES << "=1\n";
	out << VOLUME << "=1\n";
	out << TOTAL_SCANS << "=" << total_scans << "\n";
	out << PATIENT_NAME << "=" << '"' << name << '"' << "\n";
	out << PATIENT_NUMBER << "=" << '"' << patient_number << '"' << "\n";
	out << SCAN_DATE << "=" << '"' << scan_date << '"' << "\n";
	out << ROWS << "=" << rows << "\n";
	out << COLUMNS << "=" << cols << "\n";
	out << ROWVEC << "=" << rowvec << ",0.0,0.0\n";
	out << COLVEC << "=0.0," << colvec << ",0.0\n";
	out << SLICEVEC << "=0.0,0.0," << slicevec << "\n";
	if (strcmp(orientation, "NONE"))
		out << ORIENTATION << "=" << orientation << "\n";
	if (repetition_time1 != 0)
		out << REPETITION_TIME1 << "=" << repetition_time1 << "\n";
	if (echo1_time != 0)
		out << ECHO1_TIME << "=" << echo1_time << "\n";
	if (pixel_bandwidth1 != 0.0L)
		out << PIXEL_BANDWIDTH1 << "=" << pixel_bandwidth1 << "\n";
	if (repetition_time2 != 0)
		out << REPETITION_TIME2 << "=" << repetition_time2 << "\n";
	if (echo2_time != 0)
		out << ECHO2_TIME << "=" << echo2_time << "\n";
	if (pixel_bandwidth2 != 0)
		out << PIXEL_BANDWIDTH2 << "=" << pixel_bandwidth2 << "\n";

	out << XOFFSET << "=" << XYZoffset[0] << "\n";
	out << YOFFSET << "=" << XYZoffset[1] << "\n";
	out << ZOFFSET << "=" << XYZoffset[2] << "\n";
	if (imax_flag)
		out << IMAGE_MAX << "=" << image_max << "\n";
	if (imin_flag)
		out << IMAGE_MIN << "=" << image_min << "\n";
	if (spatial_units[0] != 0)
		out << SPATIAL_UNITS << "=" << spatial_units << "\n";
	out << BITS_ALLOCATED << "=" << bitsall << "\n";
	out << BITS_STORED << "=" << bitstor << "\n";
	out << HIGH_BIT << "=" << bithigh << "\n";
	out << PIXEL_REPRESENTATION << "=" << desouttypes[imagetype] << "\n";
	//  out << XYORIGIN << "=" << XYorigin[0] << "," << XYorigin[1] << "\n" ;
	out << SLICE << "=1\n";
	out << IMAGE_POSITION << "=0.0,0.0," << Position[0] << "\n";
	if (smin_flag)
		out << SLICE_MIN << "=" << slice_min[0] << "\n";
	if (smax_flag)
		out << SLICE_MAX << "=" << slice_max[0] << "\n";
	if (dscale_flag)
		out << DATA_SCALE << "=" << data_scale[0] << "\n";
	//  out << PADDING << "=0\n" ;

	if (Slice_filename[0][0] != '\"') // DIPS needs quotes around file name
	{
		for (j = strlen(Slice_filename[0]) - 1; j > 0; j--)
			if (Slice_filename[0][j] == '/')
			{
				j++;
				break;
			}
		sprintf(temp, "%s", (char *) &Slice_filename[0][j]);
		out << DATA << "=" << q << temp << q << "," << Offset[0] << "\n";
	}
	else
		out << DATA << "=" << temp << "," << Offset[0] << "\n";

	for (idx = 1; idx < total_scans; idx++)
	{
		out << SLICE << "=" << idx + 1 << "\n";
		out << IMAGE_POSITION << "=0.0,0.0," << Position[idx] << "\n";
		if (smin_flag)
			out << SLICE_MIN << "=" << slice_min[idx] << "\n";
		if (smax_flag)
			out << SLICE_MAX << "=" << slice_max[idx] << "\n";
		if (dscale_flag)
			out << DATA_SCALE << "=" << data_scale[idx] << "\n";
		if (Slice_filename[idx][0] != '\"') // DIPS needs quotes around file name
		{
			for (j = strlen(Slice_filename[idx]) - 1; j > 0; j--)
				if (Slice_filename[idx][j] == '/')
				{
					j++;
					break;
				}
			sprintf(temp, "%s", (char *) &Slice_filename[idx][j]);
			out << DATA << "=" << q << temp << q << "," << Offset[idx] << "\n";
		}
		else
			out << DATA << "=" << temp << "," << Offset[idx] << "\n";
	}

	out.close();
} // write ()

/*!
 * This routine is used to append volumes to a NEMA des header file - br
 * @param vnum - volume number
 */
void NEMAFILE::write_vol(int vnum)
{
	int idx, j;
	const char* q = "\"";
	char temp[256];

	// Open the file
	fstream out;
	out.open(filename,ios::out|ios::app);
	if (!out)
	{
		parsed = NEMA_OUTPUT_FILE_DID_NOT_OPEN;
		return;
	}

	//out << NEMA_MAGIC << "\n";
	//out << TOTAL_VOLUMES << "=1\n";
	out << VOLUME << "=" << vnum << "\n";
	out << TOTAL_SCANS << "=" << total_scans << "\n";
	out << PATIENT_NAME << "=" << '"' << name << '"' << "\n";
	out << PATIENT_NUMBER << "=" << '"' << patient_number << '"' << "\n";
	out << SCAN_DATE << "=" << '"' << scan_date << '"' << "\n";
	out << ROWS << "=" << rows << "\n";
	out << COLUMNS << "=" << cols << "\n";
	out << ROWVEC << "=" << rowvec << ",0.0,0.0\n";
	out << COLVEC << "=0.0," << colvec << ",0.0\n";
	out << SLICEVEC << "=0.0,0.0," << slicevec << "\n";
	if (strcmp(orientation, "NONE"))
		out << ORIENTATION << "=" << orientation << "\n";
	if (repetition_time1 != 0)
		out << REPETITION_TIME1 << "=" << repetition_time1 << "\n";
	if (echo1_time != 0)
		out << ECHO1_TIME << "=" << echo1_time << "\n";
	if (pixel_bandwidth1 != 0.0L)
		out << PIXEL_BANDWIDTH1 << "=" << pixel_bandwidth1 << "\n";
	if (repetition_time2 != 0)
		out << REPETITION_TIME2 << "=" << repetition_time2 << "\n";
	if (echo2_time != 0)
		out << ECHO2_TIME << "=" << echo2_time << "\n";
	if (pixel_bandwidth2 != 0)
		out << PIXEL_BANDWIDTH2 << "=" << pixel_bandwidth2 << "\n";

	out << XOFFSET << "=" << XYZoffset[0] << "\n";
	out << YOFFSET << "=" << XYZoffset[1] << "\n";
	out << ZOFFSET << "=" << XYZoffset[2] << "\n";
	if (imax_flag)
		out << IMAGE_MAX << "=" << image_max << "\n";
	if (imin_flag)
		out << IMAGE_MIN << "=" << image_min << "\n";
	if (spatial_units[0] != 0)
		out << SPATIAL_UNITS << "=" << spatial_units << "\n";
	out << BITS_ALLOCATED << "=" << bitsall << "\n";
	out << BITS_STORED << "=" << bitstor << "\n";
	out << HIGH_BIT << "=" << bithigh << "\n";
	out << PIXEL_REPRESENTATION << "=" << desouttypes[imagetype] << "\n";
	//  out << XYORIGIN << "=" << XYorigin[0] << "," << XYorigin[1] << "\n" ;
	out << SLICE << "=1\n";
	out << IMAGE_POSITION << "=0.0,0.0," << Position[0] << "\n";
	if (smin_flag)
		out << SLICE_MIN << "=" << slice_min[0] << "\n";
	if (smax_flag)
		out << SLICE_MAX << "=" << slice_max[0] << "\n";
	if (dscale_flag)
		out << DATA_SCALE << "=" << data_scale[0] << "\n";
	//  out << PADDING << "=0\n" ;

	if (Slice_filename[0][0] != '\"') // DIPS needs quotes around file name
	{
		for (j = strlen(Slice_filename[0]) - 1; j > 0; j--)
			if (Slice_filename[0][j] == '/')
			{
				j++;
				break;
			}
		sprintf(temp, "%s", (char *) &Slice_filename[0][j]);
		out << DATA << "=" << q << temp << q << "," << Offset[0] << "\n";
	}
	else
		out << DATA << "=" << temp << "," << Offset[0] << "\n";

	for (idx = 1; idx < total_scans; idx++)
	{
		out << SLICE << "=" << idx + 1 << "\n";
		out << IMAGE_POSITION << "=0.0,0.0," << Position[idx] << "\n";
		if (smin_flag)
			out << SLICE_MIN << "=" << slice_min[idx] << "\n";
		if (smax_flag)
			out << SLICE_MAX << "=" << slice_max[idx] << "\n";
		if (dscale_flag)
			out << DATA_SCALE << "=" << data_scale[idx] << "\n";
		if (Slice_filename[idx][0] != '\"') // DIPS needs quotes around file name
		{
			for (j = strlen(Slice_filename[idx]) - 1; j > 0; j--)
				if (Slice_filename[idx][j] == '/')
				{
					j++;
					break;
				}
			sprintf(temp, "%s", (char *) &Slice_filename[idx][j]);
			out << DATA << "=" << q << temp << q << "," << Offset[idx] << "\n";
		}
		else
			out << DATA << "=" << temp << "," << Offset[idx] << "\n";
	}

	out.close();
} // write ()


int NEMAFILE::num_rows() // Return the number of rows
{
	return (rows);
}

int NEMAFILE::num_cols() // Return the number of cols
{
	return (cols);
}

char *NEMAFILE::spatial_type() // Return the units of all sizes
{
	return (spatial_units);
}

double NEMAFILE::get_slicevec() // Return bits per pixel
{
	return (slicevec);
}

int NEMAFILE::bpp() // Return bits per pixel
{
	return (bitsall);
}

int NEMAFILE::bspp() // Return bits stored per pixel
{
	return (bitstor);
}

int NEMAFILE::hbp() // Return the position of the msbit
{
	return (bithigh);
}

int NEMAFILE::ptype() // Return the type of the pixel
{
	return (imagetype);
}

int NEMAFILE::tot_scans() // Return the number of scans (slices)
{
	return (total_scans);
}

double NEMAFILE::pos(int slicenum) // Return the position of slice slicenum
{
	return (Position[slicenum]);
}

char * NEMAFILE::fname(int slicenum, int rootflag) // Return the filename of slice slicenum
{
	if (!rootflag)
	{
		memset(fullname, '\0', 256);
		strcpy(fullname, directory);
		return (strcat(fullname, Slice_filename[slicenum]));
	}
	else
		return (Slice_filename[slicenum]);
}

long NEMAFILE::off(int slicenum) // Return the offset of slice slicenum
{
	if (slicenum < total_scans)
		return (Offset[slicenum]);
	return (0);
}

int NEMAFILE::status() // Return the parse status
{
	return (parsed);
}

char *NEMAFILE::desfname() // Return the filename of the descriptor file
{
	return (filename);
}

void NEMAFILE::xyorigin(double *x, double *y)
{
	*x = XYorigin[0];
	*y = XYorigin[1];
}

void NEMAFILE::set_desfname(char *fname) // Set the filename of the descriptor file
{
	int i;

	memset(filename, '\0', FILENAME_SIZE);
	strncpy(filename, fname, FILENAME_SIZE);
	// Get the path from the file name, if no path, set to null.   --jls
	memset(directory, '\0', 128);
	strcpy(directory, filename);
	for (i = strlen(directory); ((directory[i] != '/') && (i > 0)); i--)
	{
	}
	if (i > 0)
		directory[i + 1] = '\0';
	else
		directory[0] = '\0';
}

void NEMAFILE::set_name(char *patient_name)
{
	strcpy(name, patient_name);
}

void NEMAFILE::set_fname(char *fname, int slicenum) // Set the filename of the slice file
{
	if (slicenum < total_scans)
	{
		memset(Slice_filename[slicenum], '\0', FILENAME_SIZE);
		strcpy(Slice_filename[slicenum], fname);
	}
	else
		parsed = NEMA_SLICE_BEYOND_END;
}

void NEMAFILE::set_fnameallslices(char *fname) // Set same filename for all slices TLA added
{
	for (int idx = 0; idx < total_scans; idx++)
	{
		memset(Slice_filename[idx], '\0', FILENAME_SIZE); // do we need this???????
		strcpy(Slice_filename[idx], fname);
	}
}

void NEMAFILE::set_offset(long offset, int slicenum) // Set the offset of the slice within the file
{
	if (slicenum < total_scans)
		Offset[slicenum] = offset;
	else
		parsed = NEMA_SLICE_BEYOND_END;

}

void NEMAFILE::set_slice_min(double val, int slicenum)
{
	if (slicenum < total_scans)
	{
		smin_flag = 1;
		slice_min[slicenum] = val;
	}
	else
		parsed = NEMA_SLICE_BEYOND_END;

}

void NEMAFILE::set_image_min(double val)
{
	imin_flag = 1;
	image_min = val;
}

void NEMAFILE::set_slice_max(double val, int slicenum)
{
	if (slicenum < total_scans)
	{
		smax_flag = 1;
		slice_max[slicenum] = val;
	}
	else
		parsed = NEMA_SLICE_BEYOND_END;

}

void NEMAFILE::set_image_max(double val)
{
	imax_flag = 1;
	image_max = val;
}

void NEMAFILE::set_slice_data_scale(double val, int slicenum) // Set the data scale for a particular slice
{
	if (slicenum < total_scans)
	{
		data_scale[slicenum] = val;
		dscale_flag = 1;
	}
	else
		parsed = NEMA_SLICE_BEYOND_END;

}

void NEMAFILE::set_XYZoffset(XYZPOINT offset)
{
	XYZoffset[0] = offset[0];
	XYZoffset[1] = offset[1];
	XYZoffset[2] = offset[2];
	return;
}

double NEMAFILE::get_slice_min(int slicenum)
{
	if (slicenum < total_scans)
		return (slice_min[slicenum]);
	return (0);
}

double NEMAFILE::get_slice_max(int slicenum)
{
	if (slicenum < total_scans)
		return (slice_max[slicenum]);
	return (0);
}

double NEMAFILE::get_row_mm(void)
{
	return (rowvec);
}

double NEMAFILE::get_col_mm(void)
{
	return (colvec);
}

double NEMAFILE::get_image_min(void)
{
	return (image_min);
}

double NEMAFILE::get_image_max(void)
{
	return (image_max);
}

double NEMAFILE::get_slice_data_scale(int slicenum)
{
	if (slicenum < total_scans)
		return (data_scale[slicenum]);
	return (0);
}

void NEMAFILE::get_XYZoffset(XYZPOINT offset)
{
	offset[0] = XYZoffset[0];
	offset[1] = XYZoffset[1];
	offset[2] = XYZoffset[2];
	return;
}

void NEMAFILE::get_load_space(XYZPOINT load_space) // TLA added for diva
{
	load_space = Load_Space;
	return;
}

unsigned char NEMAFILE::get_load_position() // TLA added for diva
{
	return Load_Position;
}

void NEMAFILE::set_type(int nt, int nba, int nbs, int nbh)
{
	imagetype = nt;
	bitsall = nba;
	bitstor = nbs;
	bithigh = nbh;
}

void NEMAFILE::set_dims(int c, int r)
{
	cols = c;
	rows = r;
}

void NEMAFILE::set_sizes(double cs, double rs, double ss)
{
	int i;

	rowvec = rs; // Row element size
	colvec = cs; // Column element size
	slicevec = ss; // Slice element size
	for (i = 0; i < total_scans; i++)
		Position[i] = i * slicevec;

}

void NEMAFILE::set_spatial_type(char *ty)
{
	strncpy(spatial_units, ty, UNITS_SIZE);
}

void NEMAFILE::set_tot_scans(int scans)
{
	int idx;

	if (Position)
		delete Position;
	if (Slice_filename)
	{
		for (idx = 0; idx < total_scans; idx++)
			delete Slice_filename[idx];
		delete Slice_filename;
	}
	if (Offset)
		delete Offset;

	total_scans = scans;

	// Allocate space for the filenames, offsets, positions, slice_min/max, and data_scales
	Slice_filename = new char *[total_scans];
	Position = new double[total_scans];
	Offset = new long[total_scans];
	slice_min = new double[total_scans];
	slice_max = new double[total_scans];
	data_scale = new double[total_scans];

	if ((!Position) || (!Offset) || (!Slice_filename))
	{
		parsed = BAD_ALLOC;
		return;
	}

	for (idx = 0; idx < total_scans; idx++)
	{
		if ((Slice_filename[idx] = new char[FILENAME_SIZE]) == NULL)
		{
			parsed = BAD_ALLOC;
			return;
		}
		strcpy(Slice_filename[idx], "");
		Offset[idx] = 0;
		Position[idx] = 0.0;
		slice_min[idx] = 0.0;
		slice_max[idx] = 0.0;
		data_scale[idx] = 1.0;
	}

}

void NEMAFILE::set_xyorigin(double x, double y)
{
	XYorigin[0] = x;
	XYorigin[1] = y;
}

void NEMAFILE::set_pos(double pos, int slicenum) // Set the position of the slice
{
	if (slicenum < total_scans)
		Position[slicenum] = pos;
	else
		parsed = NEMA_SLICE_BEYOND_END;
}

int NEMAFILE::get_version()
{
	char buffer[10];
	int i, snv, dot_count = 0;
	float bpv;

	if (strlen(sn_version) > 0)
	{
		for (i = 0; i < (int) strlen(sn_version); i++)
			if (sn_version[i] == '.')
				break;
		strncpy(buffer, sn_version, i);
		snv = atoi(buffer);
		if (snv < 6)
			return 0;
		else
			return 1;
	}
	else
	{
		if (strlen(bp_version) > 0)
		{
			for (i = 0; i < (int) strlen(bp_version); i++)
			{
				if (bp_version[i] == '.')
				{
					dot_count++;
					if (dot_count > 1)
						break;
				}
			}
			strncpy(buffer, bp_version, i);
			bpv = atof(buffer);
			if (bpv < 3.599)
				return 0;
			else
				return 1;
		}
		else
			return -1;
	} //end else
}

/*---------------------------------------------------------------------------------*/
/*                                                                                 */
/*  MODULES CONTAINED IN THIS FILE:                                                */
/*                                                                                 */
/*      xlate  --  converts <CR> and other characters to null or space             */
/*      nema_magic  --  examines the file for the NEMA01 magic number line         */
/*      NEMA_get_cols --  extracts the number of columns from the description      */
/*      NEMA_get_rows  --  extracts the number of rows from the description        */
/*      NEMA_get_rowv  --  extracts the row vector value from the description      */
/*      NEMA_get_colv  --  extracts the column vector value from the description   */
/*      NEMA_get_slicev  --  extracts the slice vector value from the description  */
/*      NEMA_get_name  --  extracts the patient name from the description          */
/*      NEMA_get_unit  --  extracts the unit type from the description             */
/*      NEMA_get_type  --  extracts the pixel representation type from the data    */
/*      NEMA_get_tots --  extracts the number of slices from the description       */
/*      NEMA_get_bita  --  extracts the pixel depth from the data                  */
/*      NEMA_get_bits  --  extracts the number of significant bits                 */
/*      NEMA_get_bith  --  extracts the high bit number from the data              */
/*	NEMA_get_orientation -- extracts the orientation string			   */
/*      NEMA_fill_slices  -- extracts the filenames and offsets for each slice     */
/*	NEMA_fill_slices2 -- extracts the min, max, and data scales for each slice */
/*	NEMA_get_xyzoffset -- extracts xoffset, yoffset, and zoffset		   */
/*	NEMA_get_image_minmax -- extracts image_min and image_max		   */
/*	NEMA_get_version -- extracts the sn or mips version			   */
/*                                                                                 */
/***********************************************************************************/

/*---------------------------------------------------------------------------------*/
/*                                      xlate                                      */
/*---------------------------------------------------------------------------------*/
int xlate(char *buff, int n, char *bptr[])
{
	int i, z;
	char *tptr;

	/* initialize number of parsed items to zero */
	z = 0;

	/* set the temporary pointer at the start of the data buffer */
	tptr = buff;

	/* set the first buffer pointer at the "first" element */
	bptr[0] = tptr;

	/* check each character for delimiters */
	for (i = 0; i < n; i++)
		switch (buff[i])
		{
		case DELIMITER_1:
		case DELIMITER_2:
		case DELIMITER_3:
		case DELIMITER_4:
			buff[i] = 0x00;
			bptr[++z] = ++tptr;
			//    bptr[z] = ' ';  // JUNK TLA
			break;

		case SPEC_DELIMITER_1:
		case SPEC_DELIMITER_2:
			buff[i] = ' ';
			tptr++;
			break;

		default:
			tptr++;
		}

	/* return the number of items located and parsed in the data array */
	return (z);

} /* end xlate */

/*---------------------------------------------------------------------------------*/
/*                                    nema_magic                                   */
/*---------------------------------------------------------------------------------*/
int nema_magic(char *buff)
{
	const char *nemask = NEMA_MAGIC;
	int i;

	for (i = 0; i < 6; i++)
		if (buff[i] != nemask[i])
			return (0);
	return (1);

} /* end nema_magic */

/*---------------------------------------------------------------------------------*/
/*                                  NEMA_get_totv                                  */
/*---------------------------------------------------------------------------------*/
int NEMA_get_totv(char *bptr[], int n)
{
	int i;

	for (i = 0; i < n; i++)
		if (strcmp(bptr[i], TOTAL_VOLUMES) == 0)
			return (atoi(bptr[i + 1]));

	return (0);

} /* end NEMA_get_totv */

/*---------------------------------------------------------------------------------*/
/*                                  NEMA_get_totv                                  */
/*---------------------------------------------------------------------------------*/
int NEMA_get_totv(char *bptr[], int n, int *keynum)
{
	int i;

	for (i = 0; i < n; i++)
		if (strcmp(bptr[i], TOTAL_VOLUMES) == 0)
		{
			*keynum = i;
			return (atoi(bptr[i + 1]));
		}

	keynum = 0;
	return (0);

} /* end NEMA_get_totv */


/*---------------------------------------------------------------------------------*/
/*                                  NEMA_get_vol                                  */
/*---------------------------------------------------------------------------------*/
int NEMA_get_vol(char *bptr[], int n)
{
	int i;

	for (i = 0; i < n; i++)
		if (strcmp(bptr[i], VOLUME) == 0)
			return (atoi(bptr[i + 1]));

	return (0);

} /* end NEMA_get_vol */

/*---------------------------------------------------------------------------------*/
/*                                  NEMA_get_vol                                  */
/*---------------------------------------------------------------------------------*/
int NEMA_get_vol(char *bptr[], int n, int *keynum)
{
	int i;

	for (i = 0; i < n; i++)
		if (strcmp(bptr[i], VOLUME) == 0)
		{
			*keynum = i;
			return (atoi(bptr[i + 1]));
		}

	keynum = 0;
	return (0);

} /* end NEMA_get_vol */

/*---------------------------------------------------------------------------------*/
/*                                  NEMA_get_tots                                  */
/*---------------------------------------------------------------------------------*/
int NEMA_get_tots(char *bptr[], int n)
{
	int i;

	for (i = 0; i < n; i++)
		if (strcmp(bptr[i], TOTAL_SCANS) == 0)
			return (atoi(bptr[i + 1]));

	return (0);

} /* end NEMA_get_tots */

/*---------------------------------------------------------------------------------*/
/*                                 NEMA_fill_slices                                */
/*---------------------------------------------------------------------------------*/
void NEMA_fill_slices(char *bptr[], int n, char **filename, double *position,
		long *offset)
{
	int i, j; /* loop counters                               */
	int z; /* parsed slice number used for indexing       */
	double part1, part2, part3; /* 3-part vector value targets for scanf       */
	char str1[40], str2[40]; /* 2-part string for file name and offset data */

	/* examine each token in the parsed data array and look for slices */
	for (i = 0; i < n; i++)
	{
		if (strcmp(bptr[i], SLICE) == 0)
		{
			z = atoi(bptr[i + 1]) - 1;

			for (j = 1; i + j < n; j++)
			{ /* find the next image position in the list */
				if (strcmp(bptr[i + j], IMAGE_POSITION) == 0)
				{
					sscanf(bptr[i + j + 1], " %lf %lf %lf", &part1, &part2,
							&part3);
					position[z] = part3;
					j = n;
				}
			}

			for (j = 1; i + j < n; j++)
			{ /* find the next data line in the list */
				if (strcmp(bptr[i + j], DATA) == 0)
				{
					sscanf(bptr[i + j + 1], " %s %s", str1, str2);
					sprintf(filename[z], "%s", str1);
					offset[z] = atoi(str2);
					j = n;
				}
			}

		} /* end $SLICE target if block */

	} /* end for loop */

} /* end NEMA_fill_slices */

/*--------------------------------------------------------------------------------*/
/*				NEMA_fill_slices2				  */
/*--------------------------------------------------------------------------------*/
void NEMA_fill_slices2(char *bptr[], int n, double *slice_min,
		double *slice_max, double *data_scale, int *smin_flag, int *smax_flag,
		int *dscale_flag)
{
	int i, j, z;
	double val;

	for (i = 0; i < n; i++)
	{
		if (strcmp(bptr[i], SLICE) == 0) // SLICE is "$SLICE"
		{
			z = atoi(bptr[i + 1]) - 1;
			for (j = 1; j + i < n && strcmp(bptr[i + j], SLICE); j++)
			{
				if (strncmp(bptr[i + j], SLICE_MIN, 7) == 0)
				{
					sscanf(bptr[i + j + 1], "%lf", &val);
					if (strcmp(bptr[i + j], SLICE_MIN) == 0)
					{
						*smin_flag = 1;
						slice_min[z] = val;
					}
					else
					{
						*smax_flag = 1;
						slice_max[z] = val;
					}
				} // end if SLICE_M
				if (strcmp(bptr[i + j], DATA_SCALE) == 0)
				{
					sscanf(bptr[i + j + 1], "%lf", &data_scale[z]);
					*dscale_flag = 1;
				}
			} //end of j loop
		} //end if SLICE
	}//end of i loop
}

/*---------------------------------------------------------------------------------*/
/*                                 NEMA_check_slices                               */
/*---------------------------------------------------------------------------------*/
int NEMA_check_slices(char *bptr[], int n)
{
	int i; /* loop counters                               */
	int z; /* number of slices			           */

	/* examine each token in the parsed data array and look for slices */
	z = 0;
	for (i = 0; i < n; i++)
	{
		if (strcmp(bptr[i], SLICE) == 0)
			z++;
	} /* end for loop */

	return (z);

} /* end NEMA_check_slices */

/*---------------------------------------------------------------------------------*/
/*                                  NEMA_get_type                                  */
/*---------------------------------------------------------------------------------*/
int NEMA_get_type(char *bptr[], int n) // gets PIXEL_REPRESENTATION from des file
{
	int i, j;

	for (i = 0; i < n; i++)
		if (strcmp(bptr[i], PIXEL_REPRESENTATION) == 0)
			for (j = 1; j < 7; j++)
				if (strcmp(bptr[i + 1], destypes[j]) == 0)
					return (j);

	return (0);

} /* end NEMA_get_type */

/*---------------------------------------------------------------------------------*/
/*                                  NEMA_get_value                                 */
/*---------------------------------------------------------------------------------*/
double NEMA_get_value(char *bptr[], int n, const char *which)
{
	int i;

	for (i = 0; i < n; i++)
		if (strcmp(bptr[i], which) == 0)
			return (atof(bptr[i + 1]));

	return (0.0L);

} /* end NEMA_get_pbw */

/*---------------------------------------------------------------------------------*/
/*                                  NEMA_get_bita                                  */
/*---------------------------------------------------------------------------------*/
int NEMA_get_bita(char *bptr[], int n)
{
	int i;

	for (i = 0; i < n; i++)
		if (strcmp(bptr[i], BITS_ALLOCATED) == 0)
			return (atoi(bptr[i + 1]));

	return (0);

} /* end NEMA_get_bita */

/*---------------------------------------------------------------------------------*/
/*                                  NEMA_get_bits                                  */
/*---------------------------------------------------------------------------------*/
int NEMA_get_bits(char *bptr[], int n)
{
	int i;

	for (i = 0; i < n; i++)
		if (strcmp(bptr[i], BITS_STORED) == 0)
			return (atoi(bptr[i + 1]));

	return (0);

} /* end NEMA_get_bits */

/*---------------------------------------------------------------------------------*/
/*                                  NEMA_get_bith                                  */
/*---------------------------------------------------------------------------------*/
int NEMA_get_bith(char *bptr[], int n)
{
	int i;

	for (i = 0; i < n; i++)
		if (strcmp(bptr[i], HIGH_BIT) == 0)
			return (atoi(bptr[i + 1]));

	return (0);

} /* end NEMA_get_bith */

/*---------------------------------------------------------------------------------*/
/*                                  NEMA_get_scale                                 */
/*---------------------------------------------------------------------------------*/
double NEMA_get_scale(char *bptr[], int n)
{
	int i;

	for (i = 0; i < n; i++)
		if (strcmp(bptr[i], DATA_SCALE) == 0)
			return (atof(bptr[i + 1]));

	return (1); // default scale factor of 1 - br

} /* end NEMA_get_scale */

/*---------------------------------------------------------------------------------*/
/*                                NEMA_get_scale_offset - br                       */
/*---------------------------------------------------------------------------------*/
double NEMA_get_scale_offset(char *bptr[], int n)
{
	int i;

	for (i = 0; i < n; i++)
		if (strcmp(bptr[i], DATA_SCALE_OFFSET) == 0)
			return (atof(bptr[i + 1]));

	return (0); // default scale offset of 0 - br

} /* end NEMA_get_scale_offset */

/*---------------------------------------------------------------------------------*/
/*                                  NEMA_get_rows                                  */
/*---------------------------------------------------------------------------------*/
int NEMA_get_rows(char *bptr[], int n)
{
	int i;

	for (i = 0; i < n; i++)
		if (strcmp(bptr[i], ROWS) == 0)
			return (atoi(bptr[i + 1]));

	return (0);

} /* end NEMA_get_rows */

/*---------------------------------------------------------------------------------*/
/*                                  NEMA_get_cols                                  */
/*---------------------------------------------------------------------------------*/
int NEMA_get_cols(char *bptr[], int n)
{
	int i;

	for (i = 0; i < n; i++)
		if (strcmp(bptr[i], COLUMNS) == 0)
			return (atoi(bptr[i + 1]));

	return (0);

} /* end NEMA_get_cols */

/*---------------------------------------------------------------------------------*/
/*                                  NEMA_get_rowv                                  */
/*---------------------------------------------------------------------------------*/
double NEMA_get_rowv(char *bptr[], int n)
{
	int i;
	double part1, part2, part3;

	for (i = 0; i < n; i++)
	{
		if (strcmp(bptr[i], ROWVEC) == 0)
		{
			sscanf(bptr[i + 1], " %lf %lf %lf", &part1, &part2, &part3);
			return (part1);
		}
	}

	return (0.0);

} /* end NEMA_get_rowv */

/*---------------------------------------------------------------------------------*/
/*                                  NEMA_get_colv                                  */
/*---------------------------------------------------------------------------------*/
double NEMA_get_colv(char *bptr[], int n)
{
	int i;
	double part1, part2, part3;

	for (i = 0; i < n; i++)
	{
		if (strcmp(bptr[i], COLVEC) == 0)
		{
			sscanf(bptr[i + 1], " %lf %lf %lf", &part1, &part2, &part3);
			return (part2);
		}
	}

	return (0.0);

} /* end NEMA_get_colv */

/*---------------------------------------------------------------------------------*/
/*                                NEMA_get_slicev                                  */
/*---------------------------------------------------------------------------------*/
double NEMA_get_slicev(char *bptr[], int n)
{
	int i;
	double part1, part2, part3;

	for (i = 0; i < n; i++)
	{
		if (strcmp(bptr[i], SLICEVEC) == 0)
		{
			sscanf(bptr[i + 1], " %lf %lf %lf", &part1, &part2, &part3);
			return (part3);
		}
	}

	return (0.0);

} /* end NEMA_get_slicev */

/*---------------------------------------------------------------------------------*/
/*                                NEMA_get_xyorigin                                  */
/*---------------------------------------------------------------------------------*/
void NEMA_get_xyorigin(char *bptr[], int n, double *x, double *y)
{
	int i;

	for (i = 0; i < n; i++)
		if (strcmp(bptr[i], XYORIGIN) == 0)
			sscanf(bptr[i + 1], " %lf %lf", x, y);

}

/*---------------------------------------------------------------------------------*/
/*                                 NEMA_get_name                                   */
/*---------------------------------------------------------------------------------*/
void NEMA_get_name(char *bptr[], int n, char *name)
{
	int i;

	for (i = 0; i < n; i++)
	{
		if (strcmp(bptr[i], PATIENT_NAME) == 0)
		{
			sprintf(name, "%s", bptr[i + 1]);
			return;
		}
	}

	return;

} /* end NEMA_get_name */

/*---------------------------------------------------------------------------------*/
/*                                 NEMA_get_date                                   */
/*---------------------------------------------------------------------------------*/
void NEMA_get_date(char *bptr[], int n, char *date)
{
	int i;

	for (i = 0; i < n; i++)
	{
		if (strcmp(bptr[i], SCAN_DATE) == 0)
		{
			sprintf(date, "%s", bptr[i + 1]);
			return;
		}
	}

	return;

} /* end NEMA_get_date */

/*---------------------------------------------------------------------------------*/
/*                                 NEMA_get_number                                 */
/*---------------------------------------------------------------------------------*/
void NEMA_get_number(char *bptr[], int n, char *number)
{
	int i;

	for (i = 0; i < n; i++)
	{
		if (strcmp(bptr[i], PATIENT_NUMBER) == 0)
		{
			sprintf(number, "%s", bptr[i + 1]);
			return;
		}
	}

	return;

} /* end NEMA_get_number */

/*---------------------------------------------------------------------------------*/
/*                                 NEMA_get_orientation                            */
/*---------------------------------------------------------------------------------*/
void NEMA_get_orientation(char *bptr[], int n, char *orientation)
{
	int i, not_found = 1;

	for (i = 0; i < n; i++)
	{
		if (strcmp(bptr[i], ORIENTATION) == 0)
		{
			sprintf(orientation, "%s", bptr[i + 1]);
			for (i = 0; i < 3; i++)
				if (orientation[i] >= 'a')
					orientation[i] -= 'a' - 'A';
			not_found = 0;
			break;
		}
	}
	if (not_found)
		sprintf(orientation, "%s", "NONE");
	return;

} /* end NEMA_get_orientation */

/*---------------------------------------------------------------------------------*/
/*                                 NEMA_get_unit                                   */
/*---------------------------------------------------------------------------------*/
void NEMA_get_unit(char *bptr[], int n, char *name)
{
	int i;

	for (i = 0; i < n; i++)
	{
		if (strcmp(bptr[i], SPATIAL_UNITS) == 0)
		{
			sprintf(name, "%s", bptr[i + 1]);
			return;
		}
	}

	return;

} /* end NEMA_get_unit */

/*----------------------------------------------------------------------------------*/
/*				NEMA_get_xyzoffset				    */
/*----------------------------------------------------------------------------------*/

void NEMA_get_xyzoffset(char *bptr[], int n, XYZPOINT *XYZpoint) // XYZpoint is local
{
	int i, j = 0;

	for (i = 0; i < n && j < 3; i++)
	{
		if (strcmp(bptr[i], XOFFSET) == 0)
		{
			XYZpoint[0][0] = atof(bptr[i + 1]);
			j++;
		}
		if (strcmp(bptr[i], YOFFSET) == 0)
		{
			XYZpoint[0][1] = atof(bptr[i + 1]);
			j++;
		}
		if (strcmp(bptr[i], ZOFFSET) == 0)
		{
			XYZpoint[0][2] = atof(bptr[i + 1]);
			j++;
		}
	}
	return;
}

/*----------------------------------------------------------------------------------*/
/*				NEMA_get_image_minmax				    */
/*----------------------------------------------------------------------------------*/

void NEMA_get_image_minmax(char *bptr[], int n, double *image_min,
		double *image_max, int *imin_flag, int *imax_flag)
{
	int i, j = 0;

	for (i = 0; i < n && j < 2; i++)
	{
		if (strcmp(bptr[i], IMAGE_MIN) == 0)
		{
			sscanf(bptr[i + 1], "%lf", image_min);
			j++;
			*imin_flag = 1;
		}
		if (strcmp(bptr[i], IMAGE_MAX) == 0)
		{
			sscanf(bptr[i + 1], "%lf", image_max);
			j++;
			*imax_flag = 1;
		}
	}
	if (image_min == image_max)
	{
		imin_flag = imax_flag = 0;
	}
	return;
}

/*---------------------------------------------------------------------------------*/
/*                                 NEMA_get_version                                */
/*---------------------------------------------------------------------------------*/
void NEMA_get_version(char *bptr[], int n, char *sn_version, char *bp_version)
{
	int i;

	for (i = 0; i < n; i++)
	{
		if (strcmp(bptr[i], SN_VERSION) == 0)
		{
			sprintf(sn_version, "%s", bptr[i + 1]);
			sprintf(bp_version, "%s", "\0");
			return;
		}
		else
		{
			if (strcmp(bptr[i], BP_VERSION) == 0)
			{
				sprintf(bp_version, "%s", bptr[i + 1]);
				sprintf(sn_version, "%s", "\0");
				return;
			}
		}//end else
	}//end for

	sprintf(bp_version, "%s", "\0");
	sprintf(sn_version, "%s", "\0");
	return;

} /* end NEMA_get_version */

/*---------------------------------------------------------------------------------*/
/*                                 read_single_word                               */
/*---------------------------------------------------------------------------------*/

int NEMAFILE::read_single_keyword(char filename[], char key[], char result[])
{
	ifstream infile;
	// streamoff len ;
	long int len; // TLA changed - previous nolonger worked

	char *buffer;
	char **bptr;
	int i, num_keys;

	infile.open(filename);
	if (!infile)
	{
		result = NULL;
		return 0;
	}
	infile.seekg(0, ios::end);
	len = infile.tellg();
	// Allocate a buffer large enough to hold the file
	buffer = new char[len];
	if (!buffer)
	{
		result = NULL;
		return 0;
	}
	infile.seekg(0, ios::beg);
	infile.read(buffer, (int) len);
	bptr = new char *[(int) len];
	// Parse the buffer
	num_keys = xlate(buffer, (int) len, bptr);

	for (i = 0; i < num_keys; i++)
	{
		if (strcmp(bptr[i], key) == 0)
			strcpy(result, bptr[i + 1]);
	}
	delete[] buffer;
	delete[] bptr;
	return 1;

}

/* end read_single_word */

/*---------------------------------------------------------------------------------*/
/*                   LD_get_load_space no idea what this does  from diva           */
/*                   Also in ABM                                                   */
/*---------------------------------------------------------------------------------*/
void LD_get_load_space(char *bptr[], int n)
{
	int i;
	float part1, part2, part3;
	XYZPOINT Load_Space; // move to H file
	Load_Space[0] = 0.0;
	Load_Space[1] = 0.0;
	Load_Space[2] = 0.0;

	for (i = 0; i < n; i++)
		if (strcmp(bptr[i], "LOAD_SPACE") == 0)
			sscanf(bptr[i + 1], " %f %f %f", &part1, &part2, &part3);

	if (part1 > 0.0)
		Load_Space[0] = part1;
	if (part1 > 0.0)
		Load_Space[1] = part2;
	if (part1 > 0.0)
		Load_Space[2] = part3;

} /* end LD_get_load_space */

/*---------------------------------------------------------------------------------*/
/*            load_position no idea what this does  TLA copied from diva           */
/*            also used by ABM                                                     */
/*---------------------------------------------------------------------------------*/
void NEMAFILE::load_position(char *bptr[], int n)
{

#define LP_LR_LEFT       0x04
#define LP_LR_CENTER     0x08
#define LP_LR_RIGHT      0x0c
#define LP_AP_ANTERIOR   0x10
#define LP_AP_CENTER     0x20
#define LP_AP_POSTERIOR  0x30
#define LP_SI_SUPERIOR   0x40
#define LP_SI_CENTER     0x80
#define LP_SI_INFERIOR   0xc0

	int i;

	Load_Position = 0x00;

	for (i = 0; i < n; i++)
	{
		if (strcmp(bptr[i], "LOAD_POSITION") == 0)
		{
			if (strcmp(bptr[i + 1], "CENTER") == 0)
				Load_Position = 0x01;
			else if (strcmp(bptr[i + 1], "OFFSET") == 0)
				Load_Position = 0x02;
			else if (strcmp(bptr[i + 1], "ALIGNED") == 0)
			{
				Load_Position = 0x03;
				switch (bptr[i + 2][0])
				{
				case (unsigned char) 'L':
					Load_Position |= LP_LR_LEFT;
					break;
				case (unsigned char) 'C':
					Load_Position |= LP_LR_CENTER;
					break;
				case (unsigned char) 'R':
					Load_Position |= LP_LR_RIGHT;
					break;
				}
				switch (bptr[i + 2][1])
				{
				case (unsigned char) 'A':
					Load_Position |= LP_AP_ANTERIOR;
					break;
				case (unsigned char) 'C':
					Load_Position |= LP_AP_CENTER;
					break;
				case (unsigned char) 'P':
					Load_Position |= LP_AP_POSTERIOR;
					break;
				}
				switch (bptr[i + 2][2])
				{
				case (unsigned char) 'S':
					Load_Position |= LP_SI_SUPERIOR;
					break;
				case (unsigned char) 'C':
					Load_Position |= LP_SI_CENTER;
					break;
				case (unsigned char) 'I':
					Load_Position |= LP_SI_INFERIOR;
					break;
				}
			}
			else
				return;
		} /* end if */
	} /* end for */

} /* end load_position */

