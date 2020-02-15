

#include "Geometry.h"
#include <string>
//#include <iostream>
#include <fstream>

using namespace std;

string CGeometry::fnameGrid,
	   CGeometry::fbndGrid;	
int    CGeometry::NDummy   = 0;
REAL   CGeometry::Len_ref = 0.0;
	   

CGeometry::CGeometry()
{
	IM = 0;
	JM = 0;
	KM = 0;

	NBlade   = 0;
	delta_th = 0.0;
    CellVol  = 0.0;

	Node_coord = nullptr;
	FaceI	   = nullptr;
	FaceJ	   = nullptr;
	FaceK	   = nullptr;
	Cell	   = nullptr;
}

CGeometry::~CGeometry()
{
	int i, j;
	// release memory of node coord
	if(Node_coord != nullptr)
	{
		for (i=0;i<=IM;i++)
		{
			for (j=0;j<=JM;j++)
			{
				delete[] Node_coord[i][j];
			}
			delete[] Node_coord[i];
		}
		delete[] Node_coord;
		Node_coord = nullptr;
	}
	// release memory of FaceI geometry
	if(FaceI != nullptr)
	{
		for (i=0;i<=IM;i++)
		{
			for (j=0;j<JM;j++)
			{
				delete[] FaceI[i][j];
			}
			delete[] FaceI[i];
		}
		delete[] FaceI;
		FaceI = nullptr;
	}
	// release memory of FaceJ geometry
	if(FaceJ != nullptr)
	{
		for (i=0;i<IM;i++)
		{
			for (j=0;j<=JM;j++)
			{
				delete[] FaceJ[i][j];
			}
			delete[] FaceJ[i];
		}
		delete[] FaceJ;
		FaceJ = nullptr;
	}
	// release memory of FaceK geometry
	if(FaceK != nullptr)
	{
		for (i=0;i<IM;i++)
		{
			for (j=0;j<JM;j++)
			{
				delete[] FaceK[i][j];
			}
			delete[] FaceK[i];
		}
		delete[] FaceK;
		FaceK = nullptr;
	}
	// release memory of Cell geometry
	if(Cell != nullptr)
	{
		for (i=0;i<IM;i++)
		{
			for (j=0;j<JM;j++)
			{
				delete[] Cell[i][j];
			}
			delete[] Cell[i];
		}
		delete[] Cell;
		Cell = nullptr;
	}
}


void CGeometry::AllocateMemory()
{
	int i, j;
	Cell = new TCell**[IM];
	for(i=0;i<IM;i++)
		Cell[i] = new TCell*[JM];
	for(i=0;i<IM;i++)
	{
		for(j=0;j<JM;j++)
		{
			Cell[i][j] = new TCell[KM];
		}
	}

	FaceI = new TFace**[IM+1];
	for(i=0;i<=IM;i++)
		FaceI[i] = new TFace*[JM];
	for(i=0;i<=IM;i++)
	{
		for(j=0;j<JM;j++)
		{
			FaceI[i][j] = new TFace[KM];
		}
	}

	FaceJ = new TFace**[IM];
	for(i=0;i<IM;i++)
		FaceJ[i] = new TFace*[JM+1];
	for(i=0;i<IM;i++)
	{
		for(j=0;j<=JM;j++)
		{
			FaceJ[i][j] = new TFace[KM];
		}
	}

	FaceK = new TFace**[IM];
	for(i=0;i<IM;i++)
		FaceK[i] = new TFace*[JM];
	for(i=0;i<IM;i++)
	{
		for(j=0;j<JM;j++)
		{
			FaceK[i][j] = new TFace[KM+1];
		}
	}
}