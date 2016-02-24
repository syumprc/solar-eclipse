/*
*
* Template Numerical Toolkit (TNT): Linear Algebra Module
*
*
*
*/


#ifndef TNT_H
#define TNT_H



//---------------------------------------------------------------------
// Define this macro if you want  TNT to track some of the out-of-bounds
// indexing. This can encur a small run-time overhead, but is recommended 
// while developing code.  It can be turned off for production runs.
// 
//       #define TNT_BOUNDS_CHECK
//---------------------------------------------------------------------
//

//#define TNT_BOUNDS_CHECK



#include "tnt_version.h"
#include "tnt_math_utils.h"
#include "tnt_array1d.h"
#include "tnt_array2d.h"
#include "tnt_array3d.h"
#include "tnt_array1d_utils.h"
#include "tnt_array2d_utils.h"
#include "tnt_array3d_utils.h"

#include "tnt_fortran_array1d.h"
#include "tnt_fortran_array2d.h"
#include "tnt_fortran_array3d.h"
#include "tnt_fortran_array1d_utils.h"
#include "tnt_fortran_array2d_utils.h"
#include "tnt_fortran_array3d_utils.h"

#include "tnt_sparse_matrix_csr.h"

#include "tnt_stopwatch.h"
#include "tnt_subscript.h"
#include "tnt_vec.h"
#include "tnt_cmat.h"


#endif
// TNT_H
