// -------------------------- RicGraph.h ----------------------------------
/*!
@file
Header file for the RicGraph class.
*/

/*!
@mainpage
The RicGraph library contains the RicGraph and RicGraphNode classes. The
RicGraph class contains the header info from a BrainVisa graph (.arg) file
and has an array of nodes (RicGraphNode class) containing the nodes corresponding
to the sulci or gyri in the graph. The classes contain functions to read the
header and nodes from a file as well as functions to search for nodes by sulcal
or gyral label name.
*/


#ifndef RICGRAPH_H_
#define RICGRAPH_H_

#include <string>
#include <RicUtil.h>
#include "RicGraphNode.h"

#define MAXNODE 1000 ///< maximum number of graph nodes allowed

using namespace std;

// This routine breaks a string into tokens (strings) separated by delimiters.
int tokenize(const string& str, string *tokens,	const string& delimiters = " ");

/*!
Class for reading BrainVisa graph (.arg) files. This class parses the "GRAPH"
header and makes calls to the RicGraphNode class to read the nodes
corresponding to sulci or gyri in the file. A single sulcus or gyrus may be
spread across several nodes. This class has methods for searching for nodes by
sulcal label.

Most of the member variables come from the names that are listed for a node in
a BrainVisa graph file. Descriptions as to what they represent are not available.
*/
class RicGraph
{

public:
	// member variables from graph file - descriptions not available
	string	gname;
	IPoint 	boundingbox_max;
	IPoint	boundingbox_min;
	Point	voxel_size;
	string	CorticalFoldArg_VERSION;
	Point	Tal_boundingbox_max;
	Point	Tal_boundingbox_min;
	float	Talairach_rotation[8];
	Point	Talairach_scale;
	Point	Talairach_translation;
	IPoint	Tmtktri_label;
	IPoint	anterior_commissure;
	IPoint	bottom_label;
//	cortical.bck
//	cortical.global.bck
	string	datagraph_VERSION;
	string	datagraph_compatibility_model_VERSION;
	string	filename_base;
//	fold.bck
//	fold.global.bck
//	fold.global.tri
//	fold.tri
//	hull_junction.bck
//	hull_junction.global.bck
	IPoint	interhemi_point;
//	junction.bck
//	junction.global.bck
	IPoint	other_label;
//	plidepassage.bck
//	plidepassage.global.bck
	IPoint	posterior_commissure;
	IPoint	ss_label;
//	type.global.bck
//	type.global.tri

	//  variables for nodes
	int	nnodes;			///< number of nodes
	RicGraphNode *nodes;///< array of nodes from graph file
	int node_idx;		///< index to current node

	// constructors
	RicGraph();
	RicGraph(string fname);
	~RicGraph();

	// member functions
	int Read(string fname);
	int FindNodeLabel(string nodename);
	int FindNextNodeLabel(string nodename);
	int FindNodeLabelExact(string nodename);
	int FindNextNodeLabelExact(string nodename);
	int FindNodeName(string nodename);
	int FindNextNodeName(string nodename);
	int FindNodeNameExact(string nodename);
	int FindNextNodeNameExact(string nodename);

};
#endif /*RICGRAPH_H_*/
