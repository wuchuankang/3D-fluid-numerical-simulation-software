


#include "Mesh.h"

using namespace std;

void CMesh::Time_Marching()
{
	if(CTimeDiscr::timeDiscrType==TimeDiscrTypes::Runge_Kutta3)
		Runge_Kutta3_Marching();
	if(CTimeDiscr::timeDiscrType==TimeDiscrTypes::Runge_Kutta4)
		Runge_Kutta4_Marching();
	if(CTimeDiscr::timeDiscrType==TimeDiscrTypes::Runge_Kutta5)
		Runge_Kutta5_Marching();
	if(CTimeDiscr::timeDiscrType==TimeDiscrTypes::Runge_Kutta5_Hybrid)
		Runge_Kutta5_Hybrid_Marching();									//Hybrid used in central scheme£¬the other three used in upwind scheme
	if(CTimeDiscr::timeDiscrType==TimeDiscrTypes::LU_SGS)
		LU_SGS_Marching();
}