

#ifndef CBLOCKS_H
#define CBLOCKS_H

#include "Geometry.h"
#include "BndConds.h"
#include "BndGeom.h"
#include "FluidProps.h"
#include "SpaceDiscr.h"
#include "TimeDiscr.h"


class CBlocks
{
public:
	int Block_ID;

	CGeometry	Geometry;
	CBndConds	BndCond;		//boundary condition
	CFluidProps FluidProp;
	CSpaceDiscr SpaceDiscr;
	CTimeDiscr  TimeDiscr;

	int  ***flag;				//use for BL model

	// function
	CBlocks();
	~CBlocks();

	void AllocateMemory();
	void Initial_Uniform();

	void  BL_Model(const CBndGeom &BndGeom);
	
private:
	CBlocks(const CBlocks &);					// copy constructor 
	CBlocks & operator = (const CBlocks &);		// assignment operator

	void CBlocks::BL_model_1d(int M, REAL *dens, REAL *Vvel, REAL *dis, REAL *Omega, REAL *mue_l, REAL *mue_t);
	//void CBlocks::BL_model_1d(int M, REAL *dens, REAL *Vvel, REAL *dis, REAL *Omega, REAL *mue_l, REAL *mue_t, REAL Vvel_wall);
};

#endif
