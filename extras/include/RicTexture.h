// -------------------------- RicTexture.h ---------------------------------
/*! \file 
Header file of the RicTexture class
Bill Rogers - November 2008
 */

#ifndef _RIC_TEXTURE_H
#define _RIC_TEXTURE_H


class RicTextureSet; ///< just a declaration so we can reference it

/*! 
This was derived from the texture class developed by Tom Arnow but was 
renamed in order to avoid confusion with the rest of the world.
What a texture really is, is a per node value aligned with the vertices
of a mesh. The value can be any useful parameter such a curvature, 
thickness, etc. It reads and writes texture file formats from Anatomist/AIMS.
*/
class RicTexture 
{
	friend class RicTextureSet;
	
public:
	
	long 	size;	///< number of elements in texture
	float 	*nodes;	///< pointer to array of texture values
	double	min;	///< minimum texture value
	double	max;	///< maximum texture value
	double  avg;	///< average texture value
	double  med;	///< median texture value
	double	std_dev;///< standard deviation of texture values
	int		t_step;	///< time step for this mesh
	
	/// constructors
	RicTexture();
	RicTexture(int numnodes);
	~RicTexture();
	
	/// member functions
	void init(int nnodes);
	void set_node(int index,float value);
	float CalcMinMaxAvg(void);
	int write_texture(char *f_name);
	int read_texture (char *f_name);
};


#endif // _RIC_TEXTURE_H
