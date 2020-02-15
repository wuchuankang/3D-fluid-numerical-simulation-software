


#include "BndConds.h"

using namespace std;

TNode  CBndConds::InletAirAngle;
TNode  CBndConds::OutletAirAngle;
REAL   CBndConds::Ptin = 0.0,		//total pressure at inlet [Pa]
	   CBndConds::Ttin = 0.0,
	   CBndConds::Pout = 0.0,
	   CBndConds::Pout_static = 0.0,
	   CBndConds::P21ratio = 0.0,		/**< ratio of inlet to outlet static pressure (initial guess only) */
	   CBndConds::ev_ratio = 0.0;		/**< reference ev ratio, used to initialization  */

CBndConds::CBndConds()
{
	
}


CBndConds::~CBndConds()
{
	
}

	
