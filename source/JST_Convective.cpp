


#include "SpaceDiscr.h"

using namespace std;

void CSpaceDiscr::JST_Convective(const CGeometry &Geometry, const CFluidProps &FluidProp)
{
	int i, j, k, flag;
	TConsVars Fc;

	for(i=1;i<=IM;i++)
	{
		for(j=1;j<=JM-1;j++)
		{
			for(k=1;k<=KM-1;k++)
			{									
				flag = 100;
				JST_Scheme(i, j, k, flag, Geometry, FluidProp, Fc);
				FaceI_Fc[i][j][k] = Fc;
			}
		}
	}
		
	for(i=1;i<=IM-1;i++)
	{
		for(j=1;j<=JM;j++)		
		{
			for(k=1;k<=KM-1;k++)
			{
				flag = 200;
				JST_Scheme(i, j, k, flag, Geometry, FluidProp, Fc);
				FaceJ_Fc[i][j][k] = Fc;
			}
		}
	}
		
	for(i=1;i<=IM-1;i++)
	{
		for(j=1;j<=JM-1;j++)		
		{
			for(k=1;k<=KM;k++)
			{
				flag = 300;
				JST_Scheme(i, j, k, flag, Geometry, FluidProp, Fc);
				FaceK_Fc[i][j][k] = Fc;
			}
		}
	}
}

void CSpaceDiscr::JST_Scheme(int i, int j, int k, int flag, const CGeometry &Geometry, const CFluidProps &FluidProp, TConsVars &Fc)
{	
	REAL Gamma  = FluidProp.Gamma;
	REAL Vn_a, Vn_e, f_uvel_e, f_vvel_e;
	TNode S;
			
	TConsVars f_cv;
	TPrimVars f_pv;

	if(flag==100)
	{
		S = Geometry.FaceI[i][j][k].S;
		f_pv = 0.5f*(FluidProp.Cell_pv[i][j][k] + FluidProp.Cell_pv[i-1][j][k]);
		f_uvel_e = -1.0f*Geometry.FaceI[i][j][k].coord.y * FluidProp.w_rot;
		f_vvel_e = Geometry.FaceI[i][j][k].coord.x * FluidProp.w_rot;
	}
	else if(flag==200)
	{
		S = Geometry.FaceJ[i][j][k].S;
		f_pv = 0.5f*(FluidProp.Cell_pv[i][j][k] + FluidProp.Cell_pv[i][j-1][k]);
		f_uvel_e = -1.0f*Geometry.FaceJ[i][j][k].coord.y * FluidProp.w_rot;
		f_vvel_e = Geometry.FaceJ[i][j][k].coord.x * FluidProp.w_rot;
	}
	else
	{
		S = Geometry.FaceK[i][j][k].S;
		f_pv = 0.5f*(FluidProp.Cell_pv[i][j][k] + FluidProp.Cell_pv[i][j][k-1]);
		f_uvel_e = -1.0f*Geometry.FaceK[i][j][k].coord.y * FluidProp.w_rot;
		f_vvel_e = Geometry.FaceK[i][j][k].coord.x * FluidProp.w_rot;
	}

	f_cv.dens = f_pv.dens;
	f_cv.xmom = f_pv.dens * f_pv.uvel;
	f_cv.ymom = f_pv.dens * f_pv.vvel;
	f_cv.zmom = f_pv.dens * f_pv.wvel;
	Vn_a = f_pv.uvel*S.x + f_pv.vvel*S.y + f_pv.wvel*S.z;
	Vn_e = f_uvel_e*S.x + f_vvel_e*S.y + 0.0f*S.z;
	f_cv.ener = f_pv.press/(Gamma-1.0f) + f_pv.dens*(f_pv.uvel*f_pv.uvel+f_pv.vvel*f_pv.vvel+f_pv.wvel*f_pv.wvel)/2.0f;

	Fc.dens = f_cv.dens*(Vn_a-Vn_e);
	Fc.xmom = f_cv.xmom*(Vn_a-Vn_e) + f_pv.press*S.x;
	Fc.ymom = f_cv.ymom*(Vn_a-Vn_e) + f_pv.press*S.y;
	Fc.zmom = f_cv.zmom*(Vn_a-Vn_e) + f_pv.press*S.z;
	Fc.ener = (f_cv.ener + f_pv.press)*(Vn_a-Vn_e) + f_pv.press*Vn_e;	
}