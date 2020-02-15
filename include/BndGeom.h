

#ifndef BNDGEOM_H
#define BNDGEOM_H

#include "defs.h"

struct TInOutlet
{
public:
	int InOutlet_ID,
		Block_ID,
		Face_type,
		is, ie, js, je, ks, ke;
};

struct TWall
{
public:
	int Wall_ID,
		Block_ID,
		Face_type,
		is, ie, js, je, ks, ke;
	REAL rot;
	bool rotate;
};


struct TPeriodic
{
public:
	int Periodic_ID,
		Block_ID_1,
		Face_type_1,
		is_1, ie_1, js_1, je_1, ks_1, ke_1,
		orient,
		Block_ID_2,
		Face_type_2,
		is_2, ie_2, js_2, je_2, ks_2, ke_2;
	REAL delta_th;
};


struct TInterior
{
public:
	int Interior_ID,
		Block_ID_1,
		Face_type_1,
		is_1, ie_1, js_1, je_1, ks_1, ke_1,
		orient,
		Block_ID_2,
		Face_type_2,
		is_2, ie_2, js_2, je_2, ks_2, ke_2;
};

struct TMixingPlane
{
public:
	int MixingPlane_ID,
		Block_ID_1,
		Face_type_1,
		is_1, ie_1, js_1, je_1, ks_1, ke_1,
		orient,
		Block_ID_2,
		Face_type_2,
		is_2, ie_2, js_2, je_2, ks_2, ke_2;
};


class CBndGeom
{
public:
	static int  nInlet,
			    nOutlet,
			    nWall,
			    nPeriodic,
			    nInterior,
			    nMixingPlane;

	TInOutlet		*Inlet;
	TInOutlet		*Outlet;
	TWall			*Wall;
	TInterior		*Interior;
	TPeriodic		*Periodic;
	TMixingPlane	*MixingPlane;

	CBndGeom();
	~CBndGeom();
private:
	CBndGeom(const CBndGeom &BndGeom);					// override default copy constructor
	CBndGeom & operator = (const CBndGeom &BndGeom);	// and assignment operator
};

#endif