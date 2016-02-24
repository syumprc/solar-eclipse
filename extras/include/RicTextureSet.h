// -------------------------- RicTextureSet.h ---------------------------------
/*! \file 
Header file of the RicTextureSet class
Bill Rogers - November 2008
 */

#ifndef _RIC_TEXTURESET_H
#define _RIC_TEXTURESET_H

#include <string>
#include "RicTexture.h"

using namespace std;

/*!
This class deals with sets of textures. It can read and write Anatomist/AIMS
.tex texture files with with multiple time steps
*/
class RicTextureSet
{
	
public:
	
	int 		size;	///< number of elements in texture
	int			ntstep;	///< number of textures (time steps) in file
	RicTexture	*tex;	///< array of textures
	
	/// constructors
	RicTextureSet();
	RicTextureSet(int numnodes, int numtstep);
	RicTextureSet(string filename);
	~RicTextureSet();
	
	/// member functions
	void init(int nnodes, int nmode=1);
	int Read(string filename);
	int Write(string filename);
	int Write_Ascii(string filename);
};


#endif // _RIC_TEXTURESET_H
