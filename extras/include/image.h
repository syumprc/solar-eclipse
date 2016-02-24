#ifndef IMAGE_H
#define IMAGE_H
/////////////////////////////////////////////////////////////////////
/// University of Texas at San Antonio - Research Imaging Center  ///
/////////////////////////////////////////////////////////////////////
//  Module Name: image.H
//  Birth Date: 
//  Programmer: Hunter Downs
/////////////////////////////////////////////////////////////////////
/////////////////////////  Description  /////////////////////////////
// Interface for image class.
/////////////////////////////////////////////////////////////////////
//////////////////////////  Revisions  //////////////////////////////
// Programmer   , Date     : Explanation  ///////////////////////////
/////////////////////////////////////////////////////////////////////
// Dan Nickerson, 4-10-97: Changed DATA operator [] to non-virtual for
//                         improve efficiency.  Also changed IMAGE getxy
//                         from return ((*this)[y * this -> dim_x + x])
//                         to return(operator [] (y*dim_x + x)).
// Dan Nickerson, 4-10-97: Changed getxy and setxy to inline functions.
// Dan Nickerson, 4-14-97: Added DATA member "float *fdata" for common
//                         internal representation that's directly indexable.
//                         Modified getxy and setxy to operate on fdata.
// Dan Nickerson, 4-29-97: Remove data_ptr function (bad design - defeats
//                         data hiding/abstraction principles).
/////////////////////////////////////////////////////////////////////
/* $Log: image.H,v $
 * Revision 1.11  2001/01/12 23:13:36  nickersd
 * Oops, for bad reasons, set_data_type needed to be public.
 *
 * Revision 1.10  2001/01/12 23:05:50  nickersd
 * Added new memeber function "set_data_type(UNS type)" to DATA class.
 *
 * Revision 1.9  1997/11/03 19:58:10  nickersd
 * Changed new_seedfile to use the STL stack class.
 *
 * Revision 1.8  1997/10/15 15:09:09  nickersd
 * Added functions new_seedfill, invert, and changed embed_line so it returns the
 * number of pixels set in the image.
 *
 * Revision 1.7  1997/06/20 18:29:43  nickersd
 * DSN 6-20-97: Fixed the overloaded operators.  They were returning and passing
 * IMAGE by value - very bad.
 *  RCS check in comments

*/


/*SCCS stuff */
//static char SCCSid_interface_image[]="$Id: image.H,v 1.11 2001/01/12 23:13:36 nickersd Exp tom $ date:$Date: 2001/01/12 23:13:36 $";
#ifndef DTYPES_H
#include <fstream>
#include "dtypes.h"		// Some convenient typedefs
#endif
// Pixel data types
#define UNDEF_PIXELS		0
#define	FLOAT_PIXELS		1	
#define INT_PIXELS		2
#define SHORT_PIXELS		3
#define CHAR_PIXELS		4
#define UN_INT_PIXELS		5
#define UN_SHORT_PIXELS		6
#define UN_CHAR_PIXELS		7
#define DOUBLE_PIXELS		8
#define BIT_PIXELS		9
#define BS_INT_PIXELS		10
#define BS_SHORT_PIXELS		11
#define BS_UN_INT_PIXELS	12
#define BS_UN_SHORT_PIXELS	13
#define BS_UN_12BIT_PIXELS      14
#define LAST_TYPE		15
#define BS_FLOAT_PIXELS         16      // TLA added
#define TYPE_MASK	0xF

// #define LINUXSWAP  // TLA if defined - reverse bytes compile in makefile with -DLINUXSWAP

// Image ERRORs
#define	IMG_OK				0
#define IMG_INVALID_DTYPE		1
#define IMG_ALLOC_FAILED		2
#define IMG_SECOND_IMAGE_INVALID	3
#define NUM_ERRORS			4

union DATA_PTR
{
  void *V ;
  float *F ;
  double *D ;
  int *I ;
  char *C ;
  short *S ;
  unsigned int *Ui ;
  unsigned short *Us ;
  unsigned char *Uc ;
} ;

class DATA
{
 public:
  // What is the data type (on disk) of this image
  UNS type(void) const;
  // Is this a valid type
  UNS valid_type(void) const;
  // How many bytes per datum
  UNS bytes_per_datum(void) const;
  // How many bits per datum
  UNS bpd(void) const;	
  // Set the type (we use float internally).
  void set_data_type(UNS type) 
    {data_type = type; bits_per_datum = bpd();}				
  int get_data_type(void)     // return the data type TLA added 
    {return (data_type);}				
  float *fdata;    // fdata seems to hold all the data for an image after conversion TLA moved to public at least for now
protected:
  // Only used in read function 
  float get(UNL idx) const;
  // Only used in write function
  void set(UNL idx, float value) ;
  // Common internal format 
  // Data type                                       
  int data_type ;
  // How many bits per datum
  UNS bits_per_datum ;
  union DATA_PTR Img ; // What does Img do? Unconverted data read in with void * pointer
  // copied to fdata and converted - then deleted
} ;

// Class Definitions for images
class IMAGE: public DATA
{
 // X dimension of the image
 UNS dim_x ;
 // Y dimension of the image
 UNS dim_y ;
 // Offset into data file
 UNL offset ;
 // 950202:  Added to help dealing with bit images
 UNL bytes_per_image ;
 UNS error ;
public:
 IMAGE();
 IMAGE(UNS x,UNS y,UNS type,UNL off ) ;
 ~IMAGE(void) ;
 // Has an error occurred
 UNS ok(void) ;
 // What is the error message
 char *error_msg(void) ;
 // What is the x dimension of this image
 UNS x_dim(void) ;
 // What is the y dimension of this image
 UNS y_dim(void) ;
 // 950202 How many bytes per image
 UNL bytes_per_img(void) ;
 // Function to get the pixel at location x,y
 double getxy(int x,int y) const;
 // Function to set the pixel at location x,y
 void setxy(int x,int y,double value) ;
 // Offset into the file (if any)
 UNL foffset(void) ;
 // Set the file offset
 void set_offset(UNL off) ;
 // How many pixels are in this image
 UNL pixels(void) ;
 // What is the maximum pixel in the image
 double max(void) ;
 // What is the minimum pixel in the image
 double min(void) ;
 // See if the given x,y are within the image
 int inrange(int x,int y) ;
 // Copy one image to another
 void copy(const IMAGE&) ;
 // Clear the image 
 void clear(void) ;
 // Initialize an empty image
 IMAGE& operator()(UNS,UNS,UNS,UNL) ;

    // Basic Arithmetic functions:
 // Add a floating point constant to the image
 IMAGE& operator+=(const double con) ;
 // Mult. a floating point constant to image
 IMAGE& operator*=(const double con) ;
 // Add an image to the image
 IMAGE& operator+=(const IMAGE &i2) ;
 // Multiply the image by another image
 IMAGE& operator*=(const IMAGE &i2) ;
 // Equate one image to another
 IMAGE& operator=(const IMAGE &i2) ;

		// Advanced functions:
 // This function is used to find the first point
 // along the specified line in the image that has
 // a value greater than or equal to the specified
 // surface value.
 int surf_lin_int(int ,int ,int ,int ,short *,short *,double) ;
 // Simple stack based seed filling function.
 void seedfill(int x, int y, double nv) ;
 void new_seedfill(int x, int y, double nv) ;
 // Embed a line in the image from x1,y1 to x2,y2 using pix_val. Returns
 // the number of pixels set pix_val while drawing the line.
 int  embed_line(int x1, int y1, int x2, int y2, double pix_val) ;
 // get the eight neighbors of a point
 void eightneigh(int x,int y,double *nbuff) ;
 // get the direct four neighbors of a point
 void dirfourneigh(int x,int y,double *nbuff) ;
 // get the diagonal four neighbors of a point
 void diagfourneigh(int x,int y,double *nbuff) ;
 // get the pixel values along x in the specified row
 void profile_x_x(int row,double *prof) ;
 // get the pixel values along y in the specified col
 void profile_x_y(int col,double *prof) ;
 // Make a mesh grid in an image
 void mesh(int spc_x,int spc_y,double value) ;
 // Return the value for coord x,y
 double interp2D(float x,float y,double reduceover,double reduceby,
		 int *wasreduced,int *ok) ;
 // Sets all values in the image that are equal to pix_val to zero and sets
 // all other values to pix_val - good for post-processing after an attempt
 // to seed_fill and ROI.  If the seed value was outside the ROI, then seed
 // fill will fill everywhere but the inside of the ROI.  Use this function
 // to "invert" the seed fill.
 void invert(const double pix_val);
 // I/O Functions
 std::ifstream &read(std::ifstream &istr);
 std::ofstream &write(std::ofstream &ostr);
} ;

// Generic Image file read 
inline std::ifstream &operator>>(std::ifstream &stream, IMAGE &i) 
{return(i.read(stream));}
// Generic Image file write
inline std::ofstream &operator<<(std::ofstream &stream, IMAGE &i) 
{return(i.write(stream));} 

inline double IMAGE :: getxy(int x, int y) const
{
 return((double) fdata[y * dim_x + x]); // fdata seems to hold all the data
}
inline void IMAGE :: setxy(int x, int y, double val)
{
 fdata[y * dim_x + x] = (float) val;
}

//}
#endif
