

#ifndef BNDCONDS_H
#define BNDCONDS_H

#include "defs.h"
#include "BndGeom.h"
#include "FluidProps.h"
#include "Geometry.h"

class CBndConds
{
public:
	static TNode InletAirAngle,
				 OutletAirAngle;
	static REAL Ptin,			//total pressure at inlet [Pa]
				Ttin,
				Pout,
				Pout_static,
				P21ratio,		/**< ratio of outlet to inlet static pressure (initial guess only) */
				ev_ratio;		/**< reference ev ratio, used to initialization  */


	// function
	CBndConds();
	~CBndConds();

	void BndCondInlet(const TInOutlet &Inlet, const CGeometry &Geometry, CFluidProps &FluidProp);
	void BndCondOutlet_Axial(const TInOutlet &Outlet, const CGeometry &Geometry, CFluidProps &FluidProp);
	void BndCondOutlet_Centrifugal(const TInOutlet &Outlet, const CGeometry &Geometry, CFluidProps &FluidProp);
	void BndCondOutlet_General(const TInOutlet &Outlet, const CGeometry &Geometry, CFluidProps &FluidProp);
	void BndCondWall_Slip(const TWall &Wall, const CGeometry &Geometry, CFluidProps &FluidProp);
	void BndCondWall_NoSlip(const TWall &Wall, const CGeometry &Geometry, CFluidProps &FluidProp);		// not include the interior, periodic, mixing plane

private:
	CBndConds(const CBndConds &);					// copy constructor 
	CBndConds & operator = (const CBndConds &);		// assignment operator
};


#endif