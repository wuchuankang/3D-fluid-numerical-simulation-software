


#include "SpaceDiscr.h"

using namespace std;

void CSpaceDiscr::ConvectiveFlux(const CGeometry &Geometry, const CFluidProps &FluidProp)
{
	if(spaceDiscrType==SpaceDiscrTypes::JST)
	{
		JST_Convective(Geometry, FluidProp);
	}
	else
	{
		Upwind_Convective(Geometry, FluidProp);
	}
}