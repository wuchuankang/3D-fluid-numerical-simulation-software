


#include "BndGeom.h"


int CBndGeom::nInlet		= 0;
int CBndGeom::nOutlet		= 0;
int CBndGeom::nWall			= 0;
int CBndGeom::nPeriodic		= 0;
int CBndGeom::nInterior		= 0;
int CBndGeom::nMixingPlane  = 0;

CBndGeom::CBndGeom()
{
	Inlet		= nullptr;
	Outlet		= nullptr;
	Wall		= nullptr;
	Interior	= nullptr;
	Periodic	= nullptr;
	MixingPlane = nullptr;
}

CBndGeom::~CBndGeom()
{
	if(Inlet != nullptr)
	{
		delete[] Inlet;
		Inlet = nullptr;
	}
	if(Outlet!= nullptr)
	{
		delete[] Outlet;
		Outlet = nullptr;
	}
	if(Wall != nullptr)
	{
		delete[] Wall;
		Wall = nullptr;
	}
	if(Interior != nullptr)
	{
		delete[] Interior;
		Interior = nullptr;
	}
	if(Periodic != nullptr)
	{
		delete[] Periodic;
		Periodic = nullptr;
	}
	if(MixingPlane != nullptr)
	{
		delete[] MixingPlane;
		MixingPlane = nullptr;
	}
}
