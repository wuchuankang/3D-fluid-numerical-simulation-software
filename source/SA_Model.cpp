


#include "SpaceDiscr.h"
#include <fstream>
#include <iostream>
using namespace std;

void CSpaceDiscr::SA_Model(const CGeometry &Geometry, CFluidProps &FluidProp)
{
	int i, j, k;

	TCell ***Cell  = Geometry.Cell;
	TFace ***FaceI = Geometry.FaceI;
	TFace ***FaceJ = Geometry.FaceJ;
	TFace ***FaceK = Geometry.FaceK;
	
	TVisVars  ***Cell_tur  = FluidProp.Cell_tur;
	TPrimVars ***Cell_pv   = FluidProp.Cell_pv;
	TEddyVars ***Cell_ev   = FluidProp.Cell_ev;
	REAL ***gradx_I_ev	   = FluidProp.gradx_I_ev;
	REAL ***grady_I_ev	   = FluidProp.gradx_I_ev;
	REAL ***gradz_I_ev	   = FluidProp.gradx_I_ev;
						   
	REAL ***gradx_J_ev	   = FluidProp.gradx_J_ev;
	REAL ***gradz_J_ev	   = FluidProp.gradz_J_ev;
	REAL ***grady_J_ev	   = FluidProp.grady_J_ev;
						   
	REAL ***gradx_K_ev	   = FluidProp.gradx_K_ev;
	REAL ***grady_K_ev	   = FluidProp.gradz_K_ev;
	REAL ***gradz_K_ev	   = FluidProp.grady_K_ev;
						   
	REAL ***Omega		   = FluidProp.Omega;
	TNode norm;

	REAL densf, tempf, evf, ev_lam, mue_lam, tau_xx_t, tau_yy_t, tau_zz_t;
	REAL fv1, fv2, ft2, X_SA, fw, g_SA, r_SA, Sbar_SA;
	REAL Vn1, Vn2, Cell_Qc_tur;

	REAL Cb1, Cb2, Cv1, sigma, kavf;
	REAL Cw1, Cw2, Cw3, Ct1, Ct2, Ct3, Ct4;
	Cb1=0.1355f, Cb2=0.622f, Cv1=7.1f, sigma=2.0f/3.0f, kavf=0.41f;
	Cw1=Cb1/kavf/kavf+(1.0f+Cb2)/sigma;
	Cw2=0.3f, Cw3=2.0f, Ct1=1.0f, Ct2=2.0f, Ct3=1.2f, Ct4=0.5f;

	//湍流粘性
#pragma omp parallel for num_threads(NT) default(shared) private(i, j, k, ev_lam, X_SA, fv1) 
	for(i=1;i<=IM-1;i++)
	{
		for(j=1;j<=JM-1;j++)
		{
			for(k=1;k<=KM-1;k++)
			{
				ev_lam  = FluidProp.Re_ref * FluidProp.SutherLand(Cell_pv[i][j][k].temp) / Cell_pv[i][j][k].dens;		// Re!!!
				X_SA    = Cell_ev[i][j][k].ev / ev_lam;
				fv1     = POW(X_SA, 3.0f) / (POW(X_SA, 3.0f) + Cv1*Cv1*Cv1);

				Cell_tur[i][j][k].mue = fv1*Cell_pv[i][j][k].dens*Cell_ev[i][j][k].ev / FluidProp.Re_ref;		// Re!!!
			}			
		}
	}


	//对湍流粘性系数进行限制 (不能出现负值，不能超过层流粘性系数的MUT_MAX倍）
	REAL Max=200.0;
#pragma omp parallel for num_threads(NT) default(shared) private(i, j, k, mue_lam, ev_lam) 
	for(i=1;i<=IM-1;i++)
	{
		for(j=1;j<=JM-1;j++)
		{
			for(k=1;k<=KM-1;k++)
			{

				if(Cell_tur[i][j][k].mue<0.0)
					Cell_tur[i][j][k].mue = 0.0;

				if(Cell_ev[i][j][k].ev<0.0)
					Cell_ev[i][j][k].ev = 0.0;

				mue_lam = FluidProp.SutherLand(Cell_pv[i][j][k].temp);
				ev_lam  = FluidProp.Re_ref*mue_lam/Cell_pv[i][j][k].dens;
				if(mue_lam>0 && Cell_tur[i][j][k].mue>(Max*mue_lam))
					Cell_tur[i][j][k].mue = Max*mue_lam;

				if(ev_lam>0 && Cell_ev[i][j][k].ev>(Max*ev_lam))
					Cell_ev[i][j][k].ev = Max*ev_lam;
			}
		}
	}
	
#pragma omp parallel num_threads(NT) default(shared) private(i, j, k, densf, tempf, evf, ev_lam, tau_xx_t, tau_yy_t, tau_zz_t, norm, Vn1, Vn2,\
	X_SA, fv1, fv2, ft2, Sbar_SA, r_SA, g_SA, fw)
	{
		//对流项通量和粘性通量
		//I direction 
#pragma omp for nowait
		for(i=1;i<=IM;i++)
		{
			for(j=1;j<=JM-1;j++)
			{
				for(k=1;k<=KM-1;k++)
				{
					densf    = 0.5f*(Cell_pv[i][j][k].dens + Cell_pv[i-1][j][k].dens);
					tempf    = 0.5f*(Cell_pv[i][j][k].temp + Cell_pv[i-1][j][k].temp);
					evf      = 0.5f*(Cell_ev[i][j][k].ev   + Cell_ev[i-1][j][k].ev);
					ev_lam   = FluidProp.Re_ref*FluidProp.SutherLand(tempf)/densf;			// Re_ref 量纲转换
					tau_xx_t = (1.0f+Cb2)/sigma * (ev_lam+evf)*gradx_I_ev[i][j][k];
					tau_yy_t = (1.0f+Cb2)/sigma * (ev_lam+evf)*grady_I_ev[i][j][k];
					tau_zz_t = (1.0f+Cb2)/sigma * (ev_lam+evf)*gradz_I_ev[i][j][k];					

					norm.x = FaceI[i][j][k].S.x / FaceI[i][j][k].s;
					norm.y = FaceI[i][j][k].S.y / FaceI[i][j][k].s;
					norm.z = FaceI[i][j][k].S.z / FaceI[i][j][k].s;

					Vn1 = Cell_pv[i-1][j][k].uvel*norm.x + Cell_pv[i-1][j][k].vvel*norm.y + Cell_pv[i-1][j][k].wvel*norm.z;
					Vn2 = Cell_pv[i][j][k].uvel*norm.x   + Cell_pv[i][j][k].vvel*norm.y   + Cell_pv[i][j][k].wvel*norm.z;
					FaceI_Fc_tur[i][j][k].ev  = (0.5f*(Vn1+ABS(Vn1))*Cell_ev[i-1][j][k].ev + 0.5f*(Vn2-ABS(Vn2))*Cell_ev[i][j][k].ev)*FaceI[i][j][k].s;
					FaceI_Fv_tur[i][j][k].ev  = (tau_xx_t*FaceI[i][j][k].S.x + tau_yy_t*FaceI[i][j][k].S.y + tau_zz_t*FaceI[i][j][k].S.z) / FluidProp.Re_ref;		// Re_ref!!!
					FaceI_Fv1_tur[i][j][k].ev = gradx_I_ev[i][j][k]*FaceI[i][j][k].S.x + grady_I_ev[i][j][k]*FaceI[i][j][k].S.y + gradz_I_ev[i][j][k]*FaceI[i][j][k].S.z;
				}
			}
		}
		//J direction
#pragma omp for nowait
		for(i=1;i<=IM-1;i++)
		{
			for(j=1;j<=JM;j++)
			{
				for(k=1;k<=KM-1;k++)
				{
					densf    = 0.5f*(Cell_pv[i][j][k].dens + Cell_pv[i][j-1][k].dens);
					tempf    = 0.5f*(Cell_pv[i][j][k].temp + Cell_pv[i][j-1][k].temp);
					evf      = 0.5f*(Cell_ev[i][j][k].ev   + Cell_ev[i][j-1][k].ev);
					ev_lam   = FluidProp.Re_ref * FluidProp.SutherLand(tempf)/densf;			// Re_ref!!!
					tau_xx_t = (1.0f+Cb2)/sigma * (ev_lam+evf)*gradx_J_ev[i][j][k];
					tau_yy_t = (1.0f+Cb2)/sigma * (ev_lam+evf)*grady_J_ev[i][j][k];
					tau_zz_t = (1.0f+Cb2)/sigma * (ev_lam+evf)*gradz_J_ev[i][j][k];					

					norm.x = FaceJ[i][j][k].S.x / FaceJ[i][j][k].s;
					norm.y = FaceJ[i][j][k].S.y / FaceJ[i][j][k].s;
					norm.z = FaceJ[i][j][k].S.z / FaceJ[i][j][k].s;

					Vn1 = Cell_pv[i][j-1][k].uvel*norm.x + Cell_pv[i][j-1][k].vvel*norm.y + Cell_pv[i][j-1][k].wvel*norm.z;
					Vn2 = Cell_pv[i][j][k].uvel*norm.x   + Cell_pv[i][j][k].vvel*norm.y   + Cell_pv[i][j][k].wvel*norm.z;

					FaceJ_Fc_tur[i][j][k].ev  = (0.5f*(Vn1+ABS(Vn1))*Cell_ev[i][j-1][k].ev + 0.5f*(Vn2-ABS(Vn2))*Cell_ev[i][j][k].ev)*FaceJ[i][j][k].s;
					FaceJ_Fv_tur[i][j][k].ev  = (tau_xx_t*FaceJ[i][j][k].S.x + tau_yy_t*FaceJ[i][j][k].S.y + tau_zz_t*FaceJ[i][j][k].S.z) / FluidProp.Re_ref;		// Re_ref!!!
					FaceJ_Fv1_tur[i][j][k].ev = gradx_J_ev[i][j][k]*FaceJ[i][j][k].S.x + grady_J_ev[i][j][k]*FaceJ[i][j][k].S.y + gradz_J_ev[i][j][k]*FaceJ[i][j][k].S.z;
				}
			}
		}
		//K direction 
#pragma omp for
		for(i=1;i<=IM-1;i++)
		{
			for(j=1;j<=JM-1;j++)
			{
				for(k=1;k<=KM;k++)
				{
					densf    = 0.5f*(Cell_pv[i][j][k].dens + Cell_pv[i][j][k-1].dens);
					tempf    = 0.5f*(Cell_pv[i][j][k].temp + Cell_pv[i][j][k-1].temp);
					evf      = 0.5f*(Cell_ev[i][j][k].ev   + Cell_ev[i][j][k-1].ev);
					ev_lam   = FluidProp.Re_ref * FluidProp.SutherLand(tempf)/densf;				// Re_ref!!!
					tau_xx_t = (1.0f+Cb2)/sigma * (ev_lam+evf)*gradx_K_ev[i][j][k];
					tau_yy_t = (1.0f+Cb2)/sigma * (ev_lam+evf)*grady_K_ev[i][j][k];
					tau_zz_t = (1.0f+Cb2)/sigma * (ev_lam+evf)*gradz_K_ev[i][j][k];					

					norm.x = FaceK[i][j][k].S.x / FaceK[i][j][k].s;
					norm.y = FaceK[i][j][k].S.y / FaceK[i][j][k].s;
					norm.z = FaceK[i][j][k].S.z / FaceK[i][j][k].s;

					Vn1 = Cell_pv[i][j][k-1].uvel*norm.x + Cell_pv[i][j][k-1].vvel*norm.y + Cell_pv[i][j][k-1].wvel*norm.z;
					Vn2 = Cell_pv[i][j][k].uvel*norm.x   + Cell_pv[i][j][k].vvel*norm.y   + Cell_pv[i][j][k].wvel*norm.z;

					FaceK_Fc_tur[i][j][k].ev  = (0.5f*(Vn1+ABS(Vn1))*Cell_ev[i][j][k-1].ev + 0.5f*(Vn2-ABS(Vn2))*Cell_ev[i][j][k].ev)*FaceK[i][j][k].s;
					FaceK_Fv_tur[i][j][k].ev  = (tau_xx_t*FaceK[i][j][k].S.x + tau_yy_t*FaceK[i][j][k].S.y + tau_zz_t*FaceK[i][j][k].S.z) / FluidProp.Re_ref;		// Re_ref!!!
					FaceK_Fv1_tur[i][j][k].ev = gradx_K_ev[i][j][k]*FaceK[i][j][k].S.x + grady_K_ev[i][j][k]*FaceK[i][j][k].S.y + gradz_K_ev[i][j][k]*FaceK[i][j][k].S.z;
				}
			}
		}

#pragma omp for 
		//源项, 通量
		for(i=1;i<=IM-1;i++)
		{
			for(j=1;j<=JM-1;j++)
			{
				for(k=1;k<=KM-1;k++)
				{
					ev_lam  = FluidProp.Re_ref * FluidProp.SutherLand(Cell_pv[i][j][k].temp) / Cell_pv[i][j][k].dens;		// Re!!!
					X_SA    = Cell_ev[i][j][k].ev / ev_lam;
					fv1     = POW(X_SA,3.0f) / (POW(X_SA,3.0f)+Cv1*Cv1*Cv1);

					fv2     = 1.0f - X_SA/(1.0f+X_SA*fv1);
					ft2     = Ct3*EXP(0.0f-Ct4*X_SA*X_SA);																			// Re!!!
					Sbar_SA = MAX(Omega[i][j][k] + Cell_ev[i][j][k].ev*fv2/(POW(kavf*Cell[i][j][k].dis_Wall, 2.0f)*FluidProp.Re_ref), 0.3f*Omega[i][j][k]);
					r_SA    = MIN(Cell_ev[i][j][k].ev/(Sbar_SA*kavf*kavf*Cell[i][j][k].dis_Wall*Cell[i][j][k].dis_Wall*FluidProp.Re_ref), 10.0f);		// Re !!!
					g_SA    = r_SA + Cw2*(POW(r_SA,6.0f)-r_SA);
					fw      = g_SA*POW((1+POW(Cw3,6.0f))/(POW(g_SA,6.0f)+POW(Cw3,6.0f)), 1.0f/6.0f);

					//	Cell_Qc_tur = (Cb1*(1.0f-ft2)*Sbar_SA*Cell_ev[i][j][k].ev - (Cw1*fw-Cb1/(kavf*kavf)*ft2)*POW(Cell_ev[i][j][k].ev/Cell[i][j][k].dis_Wall,2.0)/FluidProp.Re_ref)*Cell[i][j][k].Vol;	// Re!!!
					Cell_Qc_tur = (Cb1*(1.0f-ft2)*Sbar_SA*Cell_ev[i][j][k].ev - (Cw1*fw-Cb1*((1-ft2)*fv2+ft2)/(kavf*kavf))*POW(Cell_ev[i][j][k].ev/Cell[i][j][k].dis_Wall,2.0)/FluidProp.Re_ref)*Cell[i][j][k].Vol;	// Re!!!

					Cell_Fv_tur[i][j][k].ev = Cell_Qc_tur + FaceI_Fv_tur[i+1][j][k].ev-FaceI_Fv_tur[i][j][k].ev + FaceJ_Fv_tur[i][j+1][k].ev-FaceJ_Fv_tur[i][j][k].ev + FaceK_Fv_tur[i][j][k+1].ev-FaceK_Fv_tur[i][j][k].ev - 
						Cb2*(ev_lam+Cell_ev[i][j][k].ev)/sigma/FluidProp.Re_ref*(FaceI_Fv1_tur[i+1][j][k].ev-FaceI_Fv1_tur[i][j][k].ev + FaceJ_Fv1_tur[i][j+1][k].ev-FaceJ_Fv1_tur[i][j][k].ev + FaceK_Fv1_tur[i][j][k+1].ev-FaceK_Fv1_tur[i][j][k].ev);
				}			// Re_ref !!!
			}
		}
	}
	
	for(j=1;j<=JM-1;j++)
	{
		for(k=1;k<=KM-1;k++)
		{
			ev_lam = FluidProp.Re_ref*FluidProp.SutherLand(Cell_pv[0][j][k].temp)/Cell_pv[0][j][k].dens;		// Re !!!
			X_SA   = Cell_ev[0][j][k].ev/ev_lam;
			fv1    = POW(X_SA,3.0f)/(POW(X_SA,3.0f)+Cv1*Cv1*Cv1);
			Cell_tur[0][j][k].mue = fv1*Cell_pv[0][j][k].dens*Cell_ev[0][j][k].ev/FluidProp.Re_ref;

			ev_lam = FluidProp.Re_ref*FluidProp.SutherLand(Cell_pv[IM][j][k].temp)/Cell_pv[IM][j][k].dens;
			X_SA   = Cell_ev[IM][j][k].ev/ev_lam;
			fv1    = POW(X_SA,3.0f)/(POW(X_SA,3.0f)+Cv1*Cv1*Cv1);
			Cell_tur[IM][j][k].mue = fv1*Cell_pv[IM][j][k].dens*Cell_ev[IM][j][k].ev/FluidProp.Re_ref;
		}
	}

	for(i=1;i<=IM-1;i++)
	{
		for(k=1;k<=KM-1;k++)
		{
			ev_lam = FluidProp.Re_ref*FluidProp.SutherLand(Cell_pv[i][0][k].temp)/Cell_pv[i][0][k].dens;
			X_SA   = Cell_ev[i][0][k].ev/ev_lam;
			fv1    = POW(X_SA,3.0f) / (POW(X_SA,3.0f)+Cv1*Cv1*Cv1);
			Cell_tur[i][0][k].mue = fv1*Cell_pv[i][0][k].dens*Cell_ev[i][0][k].ev/FluidProp.Re_ref;

			ev_lam = FluidProp.Re_ref*FluidProp.SutherLand(Cell_pv[i][JM][k].temp)/Cell_pv[i][JM][k].dens;
			X_SA   = Cell_ev[i][JM][k].ev/ev_lam;
			fv1    = POW(X_SA,3.0f)/(POW(X_SA,3.0f)+Cv1*Cv1*Cv1);
			Cell_tur[i][JM][k].mue = fv1*Cell_pv[i][JM][k].dens*Cell_ev[i][JM][k].ev/FluidProp.Re_ref;		
		}
	}

	for(i=1;i<=IM-1;i++)
	{
		for(j=1;j<=JM-1;j++)
		{
			ev_lam = FluidProp.Re_ref*FluidProp.SutherLand(Cell_pv[i][j][0].temp)/Cell_pv[i][j][0].dens;
			X_SA   = Cell_ev[i][j][0].ev/ev_lam;
			fv1    = POW(X_SA,3.0f)/(POW(X_SA,3.0f)+Cv1*Cv1*Cv1);
			Cell_tur[i][j][0].mue = fv1*Cell_pv[i][j][0].dens*Cell_ev[i][j][0].ev/FluidProp.Re_ref;

			ev_lam = FluidProp.Re_ref*FluidProp.SutherLand(Cell_pv[i][j][KM].temp)/Cell_pv[i][j][KM].dens;
			X_SA   = Cell_ev[i][j][KM].ev/ev_lam;
			fv1    = POW(X_SA,3.0f)/(POW(X_SA,3.0f)+Cv1*Cv1*Cv1);
			Cell_tur[i][j][KM].mue = fv1*Cell_pv[i][j][KM].dens*Cell_ev[i][j][KM].ev/FluidProp.Re_ref;
		}
	}
}