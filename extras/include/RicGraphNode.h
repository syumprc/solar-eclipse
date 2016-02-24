// ---------------------------- RicGraphNode.h --------------------------------
/*!
@file
Header file for class to read and manipulate nodes from a BrainVisa Graph
file.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <RicUtil.h>

#ifndef RICGRAPHNODE_H_
#define RICGRAPHNODE_H_

using namespace std;
/*!
This class contains all the info from a node in a BrainVisa graph (.arg) file.
There is a member function to read a node from a file stream.
All of the member variables come from the names that are listed for a node in
a BrainVisa graph file. I don't have descriptions as to what they are.
*/
class RicGraphNode
{
public:
	string	type;
	int		number;
	int		bottom_point_number;
	string	name;
	string	Tmtktri_filename;
	int		point_number;
	float	size;
	int		skeleton_label;
	int		ss_point_number;
	Point	Tal_boundingbox_max;
	Point	Tal_boundingbox_min;
	int		Tmtktri_label;
	int		bottom_label;
	IPoint	boundingbox_max;
	IPoint	boundingbox_min;
	Point	gravity_center;
	int		index;
	string	label;
	float	maxdepth;
	float	mindepth;
	Point	normal;
	int		other_label;
	int		other_point_number;
	Point	refgravity_center;
	Point	refnormal;
	float	refsize;
	float	refsurface_area;
	int		rootsbassin;
	int		ss_label;
	float	surface_area;
	float	talcovar[9];

	// constructors
	RicGraphNode();
	RicGraphNode(ifstream *fp);
	~RicGraphNode();

	// member functions
	int Read(ifstream *fp);
};
#endif /*RICGRAPHNODE_H_*/
