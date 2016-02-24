// ------------------------ RicMesh.h --------------------------------
/*! \file
Header file for the RicMesh class. The class can open and write
binary and ascii BrainVisa/Anatomist file formats. This version will
handle polygons of dimension 2 and 3 - lines and triangles.
This is derived from the mesh class developed by Tom Arnow.

There are also class definitions here for vertex and triangle.

 */

#ifndef _RICMESH_H
#define _RICMESH_H

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include "RicUtil.h"

using namespace std;

/// macro returning square of x
#define sqr(x) (x)*(x)
#define FP_FUDGE 0.00001		///< used for non exact comparisons

class RicMeshSet; ///< just declaring the class here so we can reference it

/*!
 This class is a container to hold a vertex. A vertex is a 3D point along
 with a general purpose integer label. The class includes a host of operations
 that one can do with a vertex.
*/
class vertex
{
	public:

	Point pnt;	///< 3D point
	int  label;	///< label

	/// constructor - empty vertex
	vertex()
	{
		pnt.x=pnt.y=pnt.z=label=0;
	}

	/// constructor - new vertex from x,y,z values
	vertex(float x, float y, float z)
	{
		pnt.x=x;pnt.y=y;pnt.z=z;
		label = 0;
	}

	/// constructor - new vertex from point
	vertex( Point p)
	{
		pnt = p;
		label = 0;
	}

	/// assign x,y,z values to a vertex
	inline void assign(float x, float y, float z)
	{
		pnt.x=x;pnt.y=y;pnt.z=z;
	}

	/// assign point to a vertex
	inline void assign(Point p)
	{
		pnt = p;
	}

	/// set vertex point from passed vertex
	inline void set(vertex fu)
	{
		pnt = fu.pnt;
	}

	/// return the dot product of vertex with passed vertex
	inline float dot(vertex fu)
	{
		return  pnt.x*fu.pnt.x +pnt.y*fu.pnt.y+ pnt.z*fu.pnt.z;

	}

	/// add passed vertex to vertex
	inline void add(vertex v)
	{
		pnt.x+=v.pnt.x;
		pnt.y+=v.pnt.y;
		pnt.z+=v.pnt.z;
	}

	/// subtract passed vertex from vertex
	inline void minus(vertex v)
	{
		pnt.x-=v.pnt.x;
		pnt.y-=v.pnt.y;
		pnt.z-=v.pnt.z;
	}

	/// subtract passed vertex from vertex
	inline void subtract(vertex v)
	{
		pnt.x-=v.pnt.x;
		pnt.y-=v.pnt.y;
		pnt.z-=v.pnt.z;
	}

	/// divide vertex by value
	inline void divide(float val)
	{
		if ( val == 0 ) return;
		pnt.x/=val;
		pnt.y/=val;
		pnt.z/=val;
	}

	/// multiply vertex by value
	inline void multiply(float val)
	{
		pnt.x*=val;
		pnt.y*=val;
		pnt.z*=val;
	}

	/// average vertex with passed vertex
	inline void average(vertex fu)
	{
		pnt.x=(pnt.x+fu.pnt.x)*0.5;
		pnt.y=(pnt.y+fu.pnt.y)*0.5;
		pnt.z=(pnt.z+fu.pnt.z)*0.5;
	}

	/// average vertex with passed vertex
	inline vertex average2(vertex fu)
	{
		vertex v;
		v.pnt.x=(pnt.x+fu.pnt.x)*0.5;
		v.pnt.y=(pnt.y+fu.pnt.y)*0.5;
		v.pnt.z=(pnt.z+fu.pnt.z)*0.5;
		return v;
	}

	/// normalize vertex
	inline void normalize()
	{
		float norm=sqrt(pnt.x*pnt.x + pnt.y*pnt.y + pnt.z*pnt.z);
		pnt.x/=norm;
		pnt.y/=norm;
		pnt.z/=norm;
	}

	/// calculate eculidian distance to passed vertex
	inline float e_distance(vertex v)
	{return sqrt( sqr(pnt.x-v.pnt.x)+sqr(pnt.y-v.pnt.y)+sqr(pnt.z-v.pnt.z));}

	/// calculate the square of the distance to passed vertex
	inline float squ_distance(vertex v)
	{return sqr(pnt.x-v.pnt.x)+sqr(pnt.y-v.pnt.y)+sqr(pnt.z-v.pnt.z);}
};

// ---------------------------------------------------------------------------
/*!
The triangle class stores the indices for the triangle vertices along with
a couple of general purpose integer labels for the triangle. Note that only
indices for vertices, not the vertices themselves are stored by the class.
*/
class triangle
{
	public:

	int vidx[3];	///< vertex indices
	int labeled;	///< label
	int  id_labeled;///< keeps reference to the normal vectors

	/// constructor - empty triangle
	triangle()
	{
		vidx[0]=vidx[1]=vidx[2]=0;
		labeled=0;
		id_labeled=0;
	}

	/// constructor - create triangle from passed vertex indices
	triangle(int v1, int v2, int v3)
	{
		vidx[0]=v1; vidx[1]=v2; vidx[2]=v3;
		labeled=0;
		id_labeled=0;
	}

	/// assign vertex indices from passed indices
	inline void assign(int v1, int v2, int v3)
	{
		vidx[0]=v1; vidx[1]=v2; vidx[2]=v3;
		labeled=0;
		id_labeled=0;
	}

	/// checks to see if vertex id belongs to triangle and assigns side and norm_id if so
	inline int check_if_belong(int id, int side, int norm_id)
	{
		for (int i=0;i<3;i++)
		{
			if (vidx[i]==id)
			{
				labeled=side;
				id_labeled=norm_id;// keeps the id of the original normal array.
				return 1;
			}
		}
		return 0;
	}

};


// ---------------------------------------------------------------------------
/*!
This class implements functions to read, write, and manipulate Brain Visa -
Anatomist meshes. This version will handle polygons of dimension 2 and 3 -
lines and triangles.
 */
class RicMesh
{
	friend class RicMeshSet;

	public:

	// variables in mesh file
	int  p_dim;		///< dimension of polygons (2 or 3)
	int v_size; 	///< number of vertices
	int n_size; 	///< number of normals
	int p_size; 	///< number of polygons
	int	t_step;		///< time step
	vertex *vertices;	///< pointer to array of vertices
	triangle *polygons;	///< pointer to array of polygons
	vertex *normals;	///< pointer to array of normals

	// spatial limits of mesh
	float	xmin,ymin,zmin;	///< minimum limits of mesh
	float	xmax,ymax,zmax;	///< maximum limits of mesh

	// methods
	RicMesh ();
	RicMesh (int vsize, int nsize, int psize ,int dim);
	RicMesh (char *fname);
	~RicMesh();
	int Init(int vs, int ns, int ps, int pd);
	int Write(string outname);
	int Read(string inname);
	int init_vertices(int vertex_size);
	int init_normals(  int normals_size);
	int init_polygons(  int polygons_size);
	void assign_node(int index, float x,float y, float z);
	void assign_node(int index, vertex z);
	void assign_normal(int index, float x,float y, float z);
	void assign_normal(int index, vertex z);
	void assign_polygon(int index, int v1,int v2, int v3);
	void calc_limits(void);
	float volume(void);
	float area(void);
	int search_for_vertex(int start_index, int end_index, vertex tmp);
	RicMesh* super_sample(float min_area);
	RicMesh* super_sample2(float min_area);
	RicMesh* triangle_fix(float max_ratio);
	RicMesh* adjust_triangle_size(float min_area, float max_area);
	int	check_vertex(vertex v);
	int point_inside_mesh(Point p);
	RicMesh*	Merge(RicMesh* InMesh);

	// maintained for historical purposes
	int read_mesh_txt(char *file_in);
	int read_mesh_bin(char *file_in);
	void write_mesh_txt(char *file_out);
	void write_mesh_bin(char *file_out);

	/// return the normal at passed index as vertex
	inline vertex get_normal(int index){ return  normals[index];}

	/// return vertex at passed index
	inline vertex get_node(int index){return vertices[index];}

	/// return triangle at passed polygon index
	inline triangle get_triangle(int index){return polygons[index];}
};


#endif // _RICMESH_H
