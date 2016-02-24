// --------------------------- RicMeshSet.h ----------------------------
/*!
 * @file
 * This is the header file for the RicMeshSet class for sets of meshes.
 */

#ifndef RICMESHSET_H_
#define RICMESHSET_H_

#include <string>
#include "RicMesh.h"

using namespace std;

/*!
This class implements functions to read, write, and create Brain Visa -
Anatomist meshes. This version will handle polygons of dimension 2 and 3 - 
lines and triangles. Multiple meshes per file are allowed.
 */
class RicMeshSet
{
public:
	
	// variables in mesh file
	int  ntstep;	///< number of time steps in the mesh
	int  p_dim;		///< dimension of polygons (2 or 3)
	int v_size; 	///< number of vertices
	int n_size; 	///< number of normals
	int p_size; 	///< number of polygons
	RicMesh *mesh;	///< array of meshes
	
	// constructors
	RicMeshSet();
	RicMeshSet(int vs, int ns, int ps, int pd, int nt);
	RicMeshSet(string filename);
	~RicMeshSet();

	// methods
	int Init(int vs, int ns, int ps, int pd, int nt);
	int FindTimeStep(int ts);

	int Read(string file_in);
	int ReadText(string file_in);
	int ReadBinary(string file_in);
	
	int Write(string file_out);
	int WriteText(string file_out);
	int WriteBinary(string file_out);	
};

#endif /*RICMESHSET_H_*/
