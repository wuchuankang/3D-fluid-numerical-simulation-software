/*
一定要注意，在WallFlux中一定不可使用FliudProp.w_rot
因为壁面可能旋转，也可能不旋转，
直接用，造成所有的壁面都旋转，显然不对
*/

#include "SpaceDiscr.h"
#include <iostream>

void CSpaceDiscr::WallFlux(const TWall &Wall,const CGeometry & Geometry, const CFluidProps & FluidProp)
{
	int i, j, k;
	//REAL pressf;
	//TNode norm;
	//for(i=Wall.is;i<=Wall.ie;i++)
	//{
	//	for(j=Wall.js;j<=Wall.je;j++)
	//	{
	//		for(k=Wall.ks;k<=Wall.ke;k++)
	//		{
	//			if(Wall.Face_type==1 || Wall.Face_type==2)
	//			{
	//				pressf = 0.5f*(FluidProp.Cell_pv[i][j][k].press + FluidProp.Cell_pv[i-1][j][k].press);
	//				FaceI_Fc[i][j][k].dens = 0;
	//				FaceI_Fc[i][j][k].xmom = pressf*Geometry.FaceI[i][j][k].S.x;
	//				FaceI_Fc[i][j][k].ymom = pressf*Geometry.FaceI[i][j][k].S.y;
	//				FaceI_Fc[i][j][k].zmom = pressf*Geometry.FaceI[i][j][k].S.z;
	//				FaceI_Fc[i][j][k].ener = 0;
	//				if(CFluidProps::turbulenceType==TurbulenceTypes::SA)
	//					FaceI_Fc_tur[i][j][k].ev = 0;
	//				if(Wall.Face_type==1)
	//					FluidProp.Cell_tur[i-1][j][k].mue = -FluidProp.Cell_tur[i][j][k].mue;		// 壁面的湍流粘性系数==0
	//				else
	//					FluidProp.Cell_tur[i][j][k].mue = -FluidProp.Cell_tur[i-1][j][k].mue;
	//			}
	//			if(Wall.Face_type==3 || Wall.Face_type==4)
	//			{
	//				pressf = 0.5f*(FluidProp.Cell_pv[i][j][k].press + FluidProp.Cell_pv[i][j-1][k].press);
	//				FaceJ_Fc[i][j][k].dens = 0;
	//				FaceJ_Fc[i][j][k].xmom = pressf*Geometry.FaceJ[i][j][k].S.x;
	//				FaceJ_Fc[i][j][k].ymom = pressf*Geometry.FaceJ[i][j][k].S.y;
	//				FaceJ_Fc[i][j][k].zmom = pressf*Geometry.FaceJ[i][j][k].S.z;
	//				FaceJ_Fc[i][j][k].ener = 0;
	//				if(CFluidProps::turbulenceType==TurbulenceTypes::SA)
	//					FaceJ_Fc_tur[i][j][k].ev = 0;		
	//				if(Wall.Face_type==3)
	//					FluidProp.Cell_tur[i][j-1][k].mue = -FluidProp.Cell_tur[i][j][k].mue;		// 壁面的湍流粘性系数==0
	//				else
	//					FluidProp.Cell_tur[i][j][k].mue = -FluidProp.Cell_tur[i][j-1][k].mue;
	//			}
	//			if(Wall.Face_type==5 || Wall.Face_type==6)
	//			{
	//				pressf = 0.5f*(FluidProp.Cell_pv[i][j][k].press + FluidProp.Cell_pv[i][j][k-1].press);
	//				FaceK_Fc[i][j][k].dens = 0;
	//				FaceK_Fc[i][j][k].xmom = pressf*Geometry.FaceK[i][j][k].S.x;
	//				FaceK_Fc[i][j][k].ymom = pressf*Geometry.FaceK[i][j][k].S.y;
	//				FaceK_Fc[i][j][k].zmom = pressf*Geometry.FaceK[i][j][k].S.z;
	//				FaceK_Fc[i][j][k].ener = 0;
	//				if(CFluidProps::turbulenceType==TurbulenceTypes::SA)
	//					FaceK_Fc_tur[i][j][k].ev = 0;
	//				if(Wall.Face_type==5)
	//					FluidProp.Cell_tur[i][j][k-1].mue = -FluidProp.Cell_tur[i][j][k].mue;		// 壁面的湍流粘性系数==0
	//				else
	//					FluidProp.Cell_tur[i][j][k].mue = -FluidProp.Cell_tur[i][j][k-1].mue;
	//			}
	//		}
	//	}
	//}

	REAL Gamma = CFluidProps::Gamma;
	REAL Vn_a, Vn_e, f_uvel_e, f_vvel_e;
	TNode S;
	TConsVars f_cv;
	TPrimVars f_pv;
	for(i=Wall.is;i<=Wall.ie;i++)
	{
		for(j=Wall.js;j<=Wall.je;j++)
		{
			for(k=Wall.ks;k<=Wall.ke;k++)
			{
				if(Wall.Face_type==1 || Wall.Face_type==2)
				{
					S = Geometry.FaceI[i][j][k].S;

					f_pv = 0.5f*(FluidProp.Cell_pv[i][j][k] + FluidProp.Cell_pv[i-1][j][k]);
					f_uvel_e = -1.0f*Geometry.FaceI[i][j][k].coord.y * Wall.rot;
					f_vvel_e = Geometry.FaceI[i][j][k].coord.x * Wall.rot;
					//f_uvel_e = -1.0f*Geometry.FaceI[i][j][k].coord.y * FluidProp.w_rot;
					//f_vvel_e = Geometry.FaceI[i][j][k].coord.x * FluidProp.w_rot;

					f_cv.dens = f_pv.dens;
					f_cv.xmom = f_pv.dens * f_pv.uvel;
					f_cv.ymom = f_pv.dens * f_pv.vvel;
					f_cv.zmom = f_pv.dens * f_pv.wvel;
					Vn_a = f_pv.uvel*S.x + f_pv.vvel*S.y + f_pv.wvel*S.z;
					Vn_e = f_uvel_e*S.x + f_vvel_e*S.y + 0.0f*S.z;
					f_cv.ener = f_pv.press/(Gamma-1.0f) + f_pv.dens*(f_pv.uvel*f_pv.uvel+f_pv.vvel*f_pv.vvel+f_pv.wvel*f_pv.wvel)/2.0f;

					FaceI_Fc[i][j][k].dens = f_cv.dens*(Vn_a-Vn_e);
					FaceI_Fc[i][j][k].xmom = f_cv.xmom*(Vn_a-Vn_e) + f_pv.press*S.x;
					FaceI_Fc[i][j][k].ymom = f_cv.ymom*(Vn_a-Vn_e) + f_pv.press*S.y;
					FaceI_Fc[i][j][k].zmom = f_cv.zmom*(Vn_a-Vn_e) + f_pv.press*S.z;
					FaceI_Fc[i][j][k].ener = (f_cv.ener + f_pv.press)*(Vn_a-Vn_e) + f_pv.press*Vn_e;
					if(CFluidProps::turbulenceType==TurbulenceTypes::SA)
						FaceI_Fc_tur[i][j][k].ev = 0;
					if(Wall.Face_type==1)
						FluidProp.Cell_tur[i-1][j][k].mue = -FluidProp.Cell_tur[i][j][k].mue;		// 壁面的湍流粘性系数==0
					else
						FluidProp.Cell_tur[i][j][k].mue = -FluidProp.Cell_tur[i-1][j][k].mue;

				}
				if(Wall.Face_type==3 || Wall.Face_type==4)
				{
					S = Geometry.FaceJ[i][j][k].S;

					f_pv = 0.5f*(FluidProp.Cell_pv[i][j][k] + FluidProp.Cell_pv[i][j-1][k]);
					f_uvel_e = -1.0f*Geometry.FaceJ[i][j][k].coord.y * Wall.rot;
					f_vvel_e = Geometry.FaceJ[i][j][k].coord.x * Wall.rot;
										
					//f_uvel_e = -1.0f*Geometry.FaceJ[i][j][k].coord.y * FluidProp.w_rot;
					//f_vvel_e = Geometry.FaceJ[i][j][k].coord.x * FluidProp.w_rot;
					f_cv.dens = f_pv.dens;
					f_cv.xmom = f_pv.dens * f_pv.uvel;
					f_cv.ymom = f_pv.dens * f_pv.vvel;
					f_cv.zmom = f_pv.dens * f_pv.wvel;
					Vn_a = f_pv.uvel*S.x + f_pv.vvel*S.y + f_pv.wvel*S.z;
					Vn_e = f_uvel_e*S.x + f_vvel_e*S.y + 0.0f*S.z;
					f_cv.ener = f_pv.press/(Gamma-1.0f) + f_pv.dens*(f_pv.uvel*f_pv.uvel+f_pv.vvel*f_pv.vvel+f_pv.wvel*f_pv.wvel)/2.0f;

					FaceJ_Fc[i][j][k].dens = f_cv.dens*(Vn_a-Vn_e);
					FaceJ_Fc[i][j][k].xmom = f_cv.xmom*(Vn_a-Vn_e) + f_pv.press*S.x;
					FaceJ_Fc[i][j][k].ymom = f_cv.ymom*(Vn_a-Vn_e) + f_pv.press*S.y;
					FaceJ_Fc[i][j][k].zmom = f_cv.zmom*(Vn_a-Vn_e) + f_pv.press*S.z;
					FaceJ_Fc[i][j][k].ener = (f_cv.ener + f_pv.press)*(Vn_a-Vn_e) + f_pv.press*Vn_e;
					if(CFluidProps::turbulenceType==TurbulenceTypes::SA)
						FaceJ_Fc_tur[i][j][k].ev = 0;		
					if(Wall.Face_type==3)
						FluidProp.Cell_tur[i][j-1][k].mue = -FluidProp.Cell_tur[i][j][k].mue;		// 壁面的湍流粘性系数==0
					else
						FluidProp.Cell_tur[i][j][k].mue = -FluidProp.Cell_tur[i][j-1][k].mue;
				}
				if(Wall.Face_type==5 || Wall.Face_type==6)
				{

					S = Geometry.FaceK[i][j][k].S;

					f_pv = 0.5f*(FluidProp.Cell_pv[i][j][k] + FluidProp.Cell_pv[i][j][k-1]);
					f_uvel_e = -1.0f*Geometry.FaceK[i][j][k].coord.y * Wall.rot;
					f_vvel_e = Geometry.FaceK[i][j][k].coord.x * Wall.rot;
										
					//f_uvel_e = -1.0f*Geometry.FaceK[i][j][k].coord.y * FluidProp.w_rot;
					//f_vvel_e = Geometry.FaceK[i][j][k].coord.x * FluidProp.w_rot;
					f_cv.dens = f_pv.dens;
					f_cv.xmom = f_pv.dens * f_pv.uvel;
					f_cv.ymom = f_pv.dens * f_pv.vvel;
					f_cv.zmom = f_pv.dens * f_pv.wvel;
					Vn_a = f_pv.uvel*S.x + f_pv.vvel*S.y + f_pv.wvel*S.z;
					Vn_e = f_uvel_e*S.x + f_vvel_e*S.y + 0.0f*S.z;
					f_cv.ener = f_pv.press/(Gamma-1.0f) + f_pv.dens*(f_pv.uvel*f_pv.uvel+f_pv.vvel*f_pv.vvel+f_pv.wvel*f_pv.wvel)/2.0f;

					FaceK_Fc[i][j][k].dens = f_cv.dens*(Vn_a-Vn_e);
					FaceK_Fc[i][j][k].xmom = f_cv.xmom*(Vn_a-Vn_e) + f_pv.press*S.x;
					FaceK_Fc[i][j][k].ymom = f_cv.ymom*(Vn_a-Vn_e) + f_pv.press*S.y;
					FaceK_Fc[i][j][k].zmom = f_cv.zmom*(Vn_a-Vn_e) + f_pv.press*S.z;
					FaceK_Fc[i][j][k].ener = (f_cv.ener + f_pv.press)*(Vn_a-Vn_e) + f_pv.press*Vn_e;
					if(CFluidProps::turbulenceType==TurbulenceTypes::SA)
						FaceK_Fc_tur[i][j][k].ev = 0;
					if(Wall.Face_type==5)
						FluidProp.Cell_tur[i][j][k-1].mue = -FluidProp.Cell_tur[i][j][k].mue;		// 壁面的湍流粘性系数==0
					else
						FluidProp.Cell_tur[i][j][k].mue = -FluidProp.Cell_tur[i][j][k-1].mue;
				}
			}
		}
	}
}