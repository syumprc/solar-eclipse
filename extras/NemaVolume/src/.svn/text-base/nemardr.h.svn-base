/////////////////////////////////////////////////////////////////////
/// University of Texas at San Antonio - Research Imaging Center  ///
/////////////////////////////////////////////////////////////////////
//  Module Name: nemardr.H
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
// Dan Nickerson, 3-25-97  : Added inline functions check_slice_min,
//                           check_slice_max, check_image_min, and
//                           check_image_max.
// Dan Nickerson, 5-7-97   : Added function check_slice_data_scale.
// Bill Rogers 3-5-10 : Added wrapper to read and write multiple headers
/////////////////////////////////////////////////////////////////////
// TLA type changed to imagetype
// nickersd had many variables called type in different routines
//                    0        1         2        3       4         5                6
// char desout = {"UNKNOWN","SIGNED","UNSIGNED","ASCII","IEEE","IEEE_DOUBLE", "UNSIGNED_BYTE
//

#ifndef NEMARDR_H
#define NEMARDR_H
#include "dtypes.h"

// Definitions of ACR-NEMA Keywords
#define NEMA_MAGIC	"NEMA01"
#define TOTAL_VOLUMES	"TOTAL_VOLUMES"
#define VOLUME		"$VOLUME"
#define SLICE		"$SLICE"
#define IMAGE_POSITION	"IMAGE_POSITION"
#define DATA		"DATA"
#define PIXEL_REPRESENTATION	"PIXEL_REPRESENTATION"
#define TOTAL_SCANS	"TOTAL_SCANS"
#define BITS_ALLOCATED	"BITS_ALLOCATED"
#define BITS_STORED	"BITS_STORED"
#define HIGH_BIT	"HIGH_BIT"
#define ROWS		"ROWS"
#define COLUMNS		"COLUMNS"
#define ROWVEC		"ROWVEC"
#define COLVEC		"COLVEC"
#define SLICEVEC	"SLICEVEC"
#define PATIENT_NAME	"PATIENT_NAME"
#define PATIENT_NUMBER  "PATIENT_NUMBER"
#define SCAN_DATE       "SCANDATE"
#define SPATIAL_UNITS	"SPATIAL_UNITS"
#define PADDING		"PADDING"
#define XYORIGIN	"XYORIGIN"
#define ORIENTATION	"ORIENTATION"
#define XOFFSET		"XOFFSET"
#define YOFFSET		"YOFFSET"
#define ZOFFSET		"ZOFFSET"
#define IMAGE_MAX	"IMAGE_MAX"
#define IMAGE_MIN	"IMAGE_MIN"
#define SLICE_MIN	"SLICE_MIN"
#define SLICE_MAX	"SLICE_MAX"
#define DATA_SCALE	"DATA_SCALE"
#define DATA_SCALE_OFFSET "DATA_SCALE_OFFSET"
#define BP_VERSION	"BP_VERSION"
#define SN_VERSION	"SN_VERSION"
#define PIXEL_BANDWIDTH1 "PIXEL_BANDWIDTH1"
#define PIXEL_BANDWIDTH2 "PIXEL_BANDWIDTH2"
#define ECHO1_TIME      "ECHO1_TIME"
#define ECHO2_TIME      "ECHO2_TIME"
#define REPETITION_TIME1 "REPETITION_TIME1"
#define REPETITION_TIME2 "REPETITION_TIME2"

// NEMA Delimiters
#define DELIMITER_1	((char)0x0D)
#define DELIMITER_2	'='
#define DELIMITER_3	EOF
#define DELIMITER_4	((char)0x0A)
#define SPEC_DELIMITER_1	'"'
#define SPEC_DELIMITER_2	','

// NEMAFILE CLASS definitions
#define FILENAME_SIZE		1024
#define UNITS_SIZE		30
#define ORIENTATION_SIZE	30
#define VERSION_SIZE		30

// Parse state definitions
#define PARSED_OK			0
#define NOT_PARSED			-1
#define NEMA_FILENAME_TOO_LONG		-2
#define NEMA_FILE_DID_NOT_OPEN		-3
#define NOT_NEMA_FILE			-4
#define COULD_NOT_ALLOC_BUFFER		-5
#define BAD_ALLOC			-6
#define BAD_TOTAL_SCANS			-7
#define NEMA_OUTPUT_FILE_DID_NOT_OPEN	-8
#define NEMA_SLICE_BEYOND_END		-9

#define NEMA_UNKNOWN			0
#define NEMA_SIGNED			1
#define NEMA_UNSIGNED			2
#define NEMA_ASCII			3
#define NEMA_IEEE_FLOAT			4

#ifdef NEMA_MOD
const char *destypes[7] =
{	"UNKNOWN","SIGNED","UNSIGNED","ASCII","IEEE_FLOAT","IEEE", "UNSIGNED_BYTE"};
const char *desouttypes[7] =
{	"UNKNOWN","SIGNED","UNSIGNED","ASCII","IEEE","IEEE", "UNSIGNED"};
#endif


class NEMAFILE
{
public: // public just to get things rolling - br
	char filename[FILENAME_SIZE]; // The Description file name
	char directory[256];
	char fullname[256];
	int parsed; // Has the file been parsed
	char name[FILENAME_SIZE]; // Patient Name
	char patient_number[FILENAME_SIZE]; // Patient Number
	char scan_date[FILENAME_SIZE]; // Scan date
	int rows; // Number of rows in the image
	int cols; // Number of columns in the image
	int total_scans; // Total Scans
	double rowvec; // Row element size
	double colvec; // Column element size
	double slicevec; // Slice element size
	char spatial_units[UNITS_SIZE]; // Spatial units description
	int bitsall; // Number of bits per pixel
	int bitstor; // Number of bits stored
	int bithigh; // High Bit
	double image_scale; // Volume scale factor
	double image_scale_offset; // Volume scale factor offset
	int imagetype; // Image Data type (was type)
	double image_max; // maximum pixel value in the volume
	double image_min; // minimum pixel value in the volume
	double *Position; // Image Positions
	char **Slice_filename; // Array of Slice Filenames
	long *Offset; // Array of Slice offsets
	double *slice_min; // Array of slice minimum pixel values
	double *slice_max; // Array of slice maximum pixel values
	double *data_scale; // Array of slice data scales
	XYZPOINT XYZoffset; // The x, y, and z offsets of the volume
	double XYorigin[2]; // The world coordinates of 0,0 in the images
	char orientation[ORIENTATION_SIZE]; // The orientation of the image data
	char sn_version[VERSION_SIZE]; // The sn version
	char bp_version[VERSION_SIZE]; // The mips version
	int smin_flag; //These flags will be set if the nemafile read contained the
	int smax_flag; //corresponding data (data scale, image and slice min/max).
	int imin_flag;
	int minc_flag; // I am a minc file instead of a des TLA
	int imax_flag;
	int dscale_flag;
	double pixel_bandwidth1;
	double pixel_bandwidth2;
	int echo1_time;
	int echo2_time;
	int repetition_time1;
	int repetition_time2;
	void load_position(char *bptr[], int n); // TLA for diva
	unsigned char Load_Position; // used by routine added for diva
	XYZPOINT Load_Space; // used by routine added for diva

public:
	void write(); // Write the ACR-NEMA file
	void write_vol(int vnum); // write a volume header to stream - br
	NEMAFILE(); // Parameterless constructor
	NEMAFILE(char *); // Constructor with filename
	~NEMAFILE();
	int init(int tot_scans=0); // Initialize all variables
	int num_rows(); // Return the number of rows
	int num_cols(); // Return the number of cols
	// Return the pixel size along rows
	double row_sz(void) const
	{
		return (rowvec);
	}
	// Return the pixel size along cols
	double col_sz(void) const
	{
		return (colvec);
	}
	// Return the pixel size along slices
	double slice_sz(void) const
	{
		return (slicevec);
	}
	char *spatial_type(); // Return the units of all sizes
	int bpp(); // Return bits per pixel
	int bspp(); // Return bits stored per pixel
	int hbp(); // Return the position of the msbit
	int ptype(); // Return the type of the pixel
	int tot_scans(); // Return the number of scans (slices)
	int tot_volumes(); // Return the number of volumes
	double get_slice_min(int slicenum); // Get a slice's minimum value
	double get_slice_max(int slicenum); // Get a slice's maximum value
	double get_image_min();
	double get_image_max();
	double get_row_mm(); // added by TLA row/col millimeters
	double get_col_mm();
	double get_slicevec(); // TLA return slice vector
	double get_slice_data_scale(int slicenum); // Get a slice's data scale
	void get_XYZoffset(XYZPOINT offset); // Get the volume's offset (XOFFSET, YOFFSET, ZOFFSET)
	unsigned char get_load_position(); //  TLA for diva and ABM
	void get_load_space(XYZPOINT load_space);// TLA added for diva

	double pos(int slicenum); // Return the position of slice slicenum
	char * fname(int slicenum, int rootflag = 0); // Return the filename of slice slicenum w/ path if !rootflag
	long off(int slicenum); // Return the offset of slice slicenum
	void xyorigin(double *x, double *y); // Return the origin (world coordinates) of the image
	int status(); // Return the status of the parse
	char *desfname(); // Return the file name of the descriptor file
	// Get the orientation string
	char *get_orientation() const
	{
		return ((char*) orientation);
	}
	int get_version(); // returns -1, 0, or 1 for "don't know", "old version"(sn < 6.0
	// or mips < 3.6), and "new version" respectively.
	double get_pixel_bandwidth1(void)
	{
		return pixel_bandwidth1;
	}
	double get_pixel_bandwidth2(void)
	{
		return pixel_bandwidth2;
	}
	int get_echo1_time(void)
	{
		return echo1_time;
	}
	int get_echo2_time(void)
	{
		return echo2_time;
	}
	int get_rep_time1(void)
	{
		return repetition_time1;
	}
	int get_rep_time2(void)
	{
		return repetition_time2;
	}
	int check_slice_min(void)
	{
		return (smin_flag);
	}
	int check_slice_max(void)
	{
		return (smax_flag);
	}
	int check_image_min(void)
	{
		return (imin_flag);
	}
	int check_image_max(void)
	{
		return (imax_flag);
	}
	int is_minc(void)
	{
		return (minc_flag);
	} // TLA retuns true if a minc file
	int check_slice_data_scale(void)
	{
		return (dscale_flag);
	}

	void set_fname(char *fname, int slicenum); // Set the file of a particular slice
	void set_fnameallslices(char *fname);// Set same filename for all slices // TLA ADDED
	void set_name(char *patient_name); // Set the patient name
	void set_number(char *patient_number); // Set patient number
	void set_scan_date(char *scan_date); // Set scan date
	void set_offset(long off, int slicenum); // Set the offset of the slice within the file
	void set_slice_min(double val, int slicenum); // Set a slice's minimum value
	void set_slice_max(double val, int slicenum); // Set a slice's maximum value
	void set_image_min(double val);
	void set_image_max(double val);
	void set_slice_data_scale(double, int); // Set a slice's data scale
	void set_XYZoffset(XYZPOINT); // Set the volume's offset (XOFFSET, YOFFSET, ZOFFSET)
	void set_desfname(char *fname); // Set the filename of the descriptor file
	void set_type(int nt, int nba, int nbs, int bh);// Set the image type in the descriptor
	void set_xyorigin(double, double); // Set the origin (world coordinates) of the image
	void set_dims(int cols, int rows); // Set the dimensions of the images
	void set_sizes(double cs, double rs, double ss); // Set the size of the voxels
	void set_spatial_type(char *imagetype); // Set the units for the voxel // TLA changed type to imagetype
	void set_tot_scans(int scans); // Set the number of scans (slices)
	void set_pos(double pos, int slicenum); // Set the position of the slice
	void set_orientation(char *orient); // Set the image data orientation
	void set_pixel_bandwidth1(const double bw)
	{
		pixel_bandwidth1 = bw;
	}
	void set_pixel_bandwidth2(const double bw)
	{
		pixel_bandwidth2 = bw;
	}
	void set_echo1_time(const int et)
	{
		echo1_time = et;
	}
	void set_echo2_time(const int et)
	{
		echo2_time = et;
	}
	void set_rep_time1(const int rt)
	{
		repetition_time1 = rt;
	}
	void set_rep_time2(const int rt)
	{
		repetition_time2 = rt;
	}
	void parseit(char *); // Open and parse a new nemafile
	void operator()(char *); // Open and parse a new nemafile
	void duplicate(NEMAFILE *); // Duplicate an instance of a nema class

	//  The function to read a single keyword from the DES file.
	// returns 0 if not found and 1 if it is found;
	static int read_single_keyword(char* filename, char* key, char* result);
	// int NEMAFILE::read_single_keyword(char filename[], char key[], char result[])

};


/*!
 * Function to read in entire header for a NEMA file that may include
 * multiple volumes.
 * @param NemaVol - pointer to array of NEMA headers
 * @param filename - name of NEMA file to open
 * @return - number of volumes in file
 */
int NEMA_READ(NEMAFILE **NemaVol, char *filename);

/*!
 * Function to write the beginning of the header for multiple volume NEMA files
 * @param fname - file name to write to
 * @param nvol - number of volumes
 * @return
 */
int	NEMA_Start_Header(char *fname, int nvol);

#endif
