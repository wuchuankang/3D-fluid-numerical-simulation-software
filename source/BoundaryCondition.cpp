

#include "Mesh.h"
#include <iostream>
#include "Geometry.h"

using namespace std;

void CMesh::BoundaryCondition()
{
	int i, Block_ID, Block_ID_1, Block_ID_2;

	for(i=1;i<=BndGeom.nInlet;i++)
	{
		Block_ID = BndGeom.Inlet[i].Block_ID;

		Block[Block_ID].BndCond.BndCondInlet(BndGeom.Inlet[i], Block[Block_ID].Geometry, Block[Block_ID].FluidProp);
	}

	for(i=1;i<=BndGeom.nOutlet;i++)
	{
		Block_ID = BndGeom.Outlet[i].Block_ID;
		if(TurboType==TurboTypes::Axial)
			Block[Block_ID].BndCond.BndCondOutlet_Axial(BndGeom.Outlet[i], Block[Block_ID].Geometry, Block[Block_ID].FluidProp);
		else if(TurboType==TurboTypes::Centrifugal)
			Block[Block_ID].BndCond.BndCondOutlet_Centrifugal(BndGeom.Outlet[i], Block[Block_ID].Geometry, Block[Block_ID].FluidProp);
		else
			Block[Block_ID].BndCond.BndCondOutlet_General(BndGeom.Outlet[i], Block[Block_ID].Geometry, Block[Block_ID].FluidProp);
	}

	for(i=1;i<=BndGeom.nWall;i++)
	{
		Block_ID = BndGeom.Wall[i].Block_ID;
		if(BndGeom.Wall[i].rotate)
		{
			BndGeom.Wall[i].rot = Block[Block_ID].FluidProp.w_rot;
		}
		else
		{
			BndGeom.Wall[i].rot = 0.0;
		}
		Block[Block_ID].BndCond.BndCondWall_NoSlip(BndGeom.Wall[i], Block[Block_ID].Geometry, Block[Block_ID].FluidProp);
	}

	for(i=1;i<=BndGeom.nInterior;i++)
	{
		Block_ID_1 = BndGeom.Interior[i].Block_ID_1;
		Block_ID_2 = BndGeom.Interior[i].Block_ID_2;
		BndCondInterior(Block[Block_ID_1],Block[Block_ID_2], BndGeom.Interior[i]);							
	}

	for(i=1;i<=BndGeom.nPeriodic;i++)
	{
		Block_ID_1 = BndGeom.Periodic[i].Block_ID_1;
		Block_ID_2 = BndGeom.Periodic[i].Block_ID_2;
		BndCondPeriodic(Block[Block_ID_1],Block[Block_ID_2], BndGeom.Periodic[i]);							
	}

	for(i=1;i<=BndGeom.nMixingPlane;i++)
	{
		BndCondMixingPlane(BndGeom.MixingPlane[i]);
	}
}
