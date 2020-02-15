

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <string>
#include "defs.h"
#include "Geometry.h"
#include "BndGeom.h"

using std::string;

class CGeometry
{
public:
	static string fnameGrid,	// grids file
		          fbndGrid;		// boundary geometry definition file
	static int    NDummy;		// number of dummy layers
	int			  NBlade;		// number of blades, 动静子的叶片数不同
	REAL		  delta_th;
	static REAL   Len_ref;		// reference length
				  		
	int		IM, JM, KM;			// i(max index), j(max index), k(max index) of every block
	REAL	CellVol;			// total volume of every block

	TNode	***Node_coord;		// coordinate
	TFace	***FaceI,			// vector area
			***FaceJ,
			***FaceK;
	TCell	***Cell;


	// function
	CGeometry();
	~CGeometry();
	void  AllocateMemory();
	void  FaceVectorVolume();
	void  JacobiConvertCoeff();
private:
	CGeometry(const CGeometry &);					// copy constructor 
	CGeometry & operator = (const CGeometry &);		// assignment operator
};

#endif