

#include "SpaceDiscr.h"
#include <fstream>
#include <iostream>
using namespace std;


void CSpaceDiscr::ViscousFlux(const CGeometry &Geometry, CFluidProps &FluidProp)
{
	int i, j, k;

	REAL Gamma  = FluidProp.Gamma;
	REAL Pr_lam = FluidProp.Pr_lam;
	REAL Pr_tur = FluidProp.Pr_tur;
	REAL Rgas   = FluidProp.Rgas;
	TPrimVars  ***gradx_I = FluidProp.gradx_I;
	TPrimVars  ***grady_I = FluidProp.grady_I;
	TPrimVars  ***gradz_I = FluidProp.gradz_I;
	TPrimVars  ***gradx_J = FluidProp.gradx_J;
	TPrimVars  ***grady_J = FluidProp.grady_J;
	TPrimVars  ***gradz_J = FluidProp.gradz_J;
	TPrimVars  ***gradx_K = FluidProp.gradx_K;
	TPrimVars  ***grady_K = FluidProp.grady_K;
	TPrimVars  ***gradz_K = FluidProp.gradz_K;

	REAL mue, mue_lam, mue_tur, cond, cond_lam, cond_tur;
	TNode S, F1, F2, F3, F4, F5;
	REAL tau_xx, tau_xy, tau_xz, tau_yy, tau_yz, tau_zz;
	REAL uvelf, vvelf ,wvelf, tempf;
	
#pragma omp parallel num_threads(NT) default(shared) firstprivate(i,j,k, mue, mue_lam, mue_tur, cond, cond_lam, cond_tur, S, F1, F2, F3, F4, F5, tau_xx, tau_xy, tau_xz, tau_yy, tau_yz, tau_zz, uvelf, vvelf ,wvelf, tempf)
	{
#pragma omp for nowait
		for(i=1;i<=IM;i++)
		{
			for(j=1;j<=JM-1;j++)
			{
				for(k=1;k<=KM-1;k++)
				{
					uvelf   = 0.5f*(FluidProp.Cell_pv[i][j][k].uvel + FluidProp.Cell_pv[i-1][j][k].uvel);
					vvelf   = 0.5f*(FluidProp.Cell_pv[i][j][k].vvel + FluidProp.Cell_pv[i-1][j][k].vvel);
					wvelf   = 0.5f*(FluidProp.Cell_pv[i][j][k].wvel + FluidProp.Cell_pv[i-1][j][k].wvel);
					tempf   = 0.5f*(FluidProp.Cell_pv[i][j][k].temp + FluidProp.Cell_pv[i-1][j][k].temp);

					if(tempf>=10000 || tempf<=0)
					{
						//tempf = 2E2;
						cout<<"ViscousFlux I"<<endl;
					}
					mue_lam  = FluidProp.SutherLand(tempf);
					cond_lam = mue_lam*Gamma*Rgas/(Gamma-1.0f)/Pr_lam;

					if (FluidProp.turbulenceType==TurbulenceTypes::SA)
					{
						mue_tur  = (FluidProp.Cell_tur[i-1][j][k].mue + FluidProp.Cell_tur[i][j][k].mue)/2.0f;
						cond_tur = mue_tur*Gamma*Rgas / (Gamma-1.0f) / Pr_tur;					
					}
					else if(FluidProp.turbulenceType==TurbulenceTypes::BL)
					{
						mue_tur  = (FluidProp.Cell_tur[i-1][j][k].mue + FluidProp.Cell_tur[i][j][k].mue)/2.0f;
						cond_tur = mue_tur*Gamma*Rgas / (Gamma-1.0f) / Pr_tur;
					}
					else if(FluidProp.turbulenceType==TurbulenceTypes::SST)
					{

					}
					else
					{
						mue_tur  = 0.0;
						cond_tur = 0.0;
					}
					mue  = mue_lam  + mue_tur;			// if we need plot the mue, we should make Cell_vs[][][].mue
					cond = cond_lam + cond_tur;

					tau_xx = mue * (4.0f/3.0f*gradx_I[i][j][k].uvel - 2.0f/3.0f*(grady_I[i][j][k].vvel+gradz_I[i][j][k].wvel));
					tau_xy = mue * (grady_I[i][j][k].uvel + gradx_I[i][j][k].vvel);
					tau_xz = mue * (gradz_I[i][j][k].uvel + gradx_I[i][j][k].wvel);
					tau_yy = mue * (4.0f/3.0f*grady_I[i][j][k].vvel - 2.0f/3.0f*(gradx_I[i][j][k].uvel+gradz_I[i][j][k].wvel));
					tau_yz = mue * (gradz_I[i][j][k].vvel + grady_I[i][j][k].wvel);
					tau_zz = mue * (4.0f/3.0f*gradz_I[i][j][k].wvel - 2.0f/3.0f*(gradx_I[i][j][k].uvel+grady_I[i][j][k].vvel));

					//				out<<tau_xx<<"    "<<tau_xy<<"   "<<tau_xz<<"   "<<tau_yy<<"   "<<tau_yz<<"   "<<tau_yy<<"   "<<tau_zz<<endl;

					F1.x = 0.0;
					F1.y = 0.0;
					F1.z = 0.0;
					F2.x = tau_xx;
					F2.y = tau_xy;
					F2.z = tau_xz;
					F3.x = tau_xy;
					F3.y = tau_yy;
					F3.z = tau_yz;
					F4.x = tau_xz;
					F4.y = tau_yz;
					F4.z = tau_zz;
					F5.x = uvelf*tau_xx + vvelf*tau_xy + wvelf*tau_xz + cond*gradx_I[i][j][k].temp;
					F5.y = uvelf*tau_xy + vvelf*tau_yy + wvelf*tau_yz + cond*grady_I[i][j][k].temp;
					F5.z = uvelf*tau_xz + vvelf*tau_yz + wvelf*tau_zz + cond*gradz_I[i][j][k].temp;

					S.x = Geometry.FaceI[i][j][k].S.x;
					S.y = Geometry.FaceI[i][j][k].S.y;
					S.z = Geometry.FaceI[i][j][k].S.z;

					FaceI_Fv[i][j][k].dens = F1*S;
					FaceI_Fv[i][j][k].xmom = F2*S;
					FaceI_Fv[i][j][k].ymom = F3*S;
					FaceI_Fv[i][j][k].zmom = F4*S;
					FaceI_Fv[i][j][k].ener = F5*S;
				}
			}
		}

#pragma omp for nowait
		for(i=1;i<=IM-1;i++)
		{
			for(j=1;j<=JM;j++)
			{
				for(k=1;k<=KM-1;k++)
				{
					uvelf = 0.5f*(FluidProp.Cell_pv[i][j][k].uvel + FluidProp.Cell_pv[i][j-1][k].uvel);
					vvelf = 0.5f*(FluidProp.Cell_pv[i][j][k].vvel + FluidProp.Cell_pv[i][j-1][k].vvel);
					wvelf = 0.5f*(FluidProp.Cell_pv[i][j][k].wvel + FluidProp.Cell_pv[i][j-1][k].wvel);
					tempf = 0.5f*(FluidProp.Cell_pv[i][j][k].temp + FluidProp.Cell_pv[i][j-1][k].temp);

					if(tempf>=10000 || tempf<=0)
					{
						//tempf = 2E2;
						cout<<"ViscousFlux J"<<endl;
					}
					mue_lam  = FluidProp.SutherLand(tempf);
					cond_lam = mue_lam*Gamma*Rgas/(Gamma-1.0f)/Pr_lam;

					if (FluidProp.turbulenceType==TurbulenceTypes::SA)
					{
						mue_tur  = (FluidProp.Cell_tur[i][j-1][k].mue + FluidProp.Cell_tur[i][j][k].mue)/2.0f;
						cond_tur = mue_tur*Gamma*Rgas / (Gamma-1.0f) / Pr_tur;					
					}
					else if(FluidProp.turbulenceType==TurbulenceTypes::BL)
					{
						mue_tur  = (FluidProp.Cell_tur[i][j-1][k].mue + FluidProp.Cell_tur[i][j][k].mue)/2.0f;
						cond_tur = mue_tur*Gamma*Rgas / (Gamma-1.0f) / Pr_tur;
					}
					else if(FluidProp.turbulenceType==TurbulenceTypes::SST)
					{

					}
					else
					{
						mue_tur  = 0.0;
						cond_tur = 0.0;
					}
					mue  = mue_lam  + mue_tur;			// if we need plot the mue, we should make Cell_vs[][][].mue
					cond = cond_lam + cond_tur;

					tau_xx = mue * (4.0f/3.0f*gradx_J[i][j][k].uvel - 2.0f/3.0f*(grady_J[i][j][k].vvel+gradz_J[i][j][k].wvel));
					tau_xy = mue * (grady_J[i][j][k].uvel + gradx_J[i][j][k].vvel);
					tau_xz = mue * (gradz_J[i][j][k].uvel + gradx_J[i][j][k].wvel);
					tau_yy = mue * (4.0f/3.0f*grady_J[i][j][k].vvel - 2.0f/3.0f*(gradx_J[i][j][k].uvel+gradz_J[i][j][k].wvel));
					tau_yz = mue * (gradz_J[i][j][k].vvel + grady_J[i][j][k].wvel);
					tau_zz = mue * (4.0f/3.0f*gradz_J[i][j][k].wvel - 2.0f/3.0f*(gradx_J[i][j][k].uvel+grady_J[i][j][k].vvel));

					F1.x = 0.0;
					F1.y = 0.0;
					F1.z = 0.0;
					F2.x = tau_xx;
					F2.y = tau_xy;
					F2.z = tau_xz;
					F3.x = tau_xy;
					F3.y = tau_yy;
					F3.z = tau_yz;
					F4.x = tau_xz;
					F4.y = tau_yz;
					F4.z = tau_zz;
					F5.x = uvelf*tau_xx + vvelf*tau_xy + wvelf*tau_xz + cond*gradx_J[i][j][k].temp;;
					F5.y = uvelf*tau_xy + vvelf*tau_yy + wvelf*tau_yz + cond*grady_J[i][j][k].temp;;
					F5.z = uvelf*tau_xz + vvelf*tau_yz + wvelf*tau_zz + cond*gradz_J[i][j][k].temp;;

					S.x = Geometry.FaceJ[i][j][k].S.x;
					S.y = Geometry.FaceJ[i][j][k].S.y;
					S.z = Geometry.FaceJ[i][j][k].S.z;

					FaceJ_Fv[i][j][k].dens = F1*S;
					FaceJ_Fv[i][j][k].xmom = F2*S;
					FaceJ_Fv[i][j][k].ymom = F3*S;
					FaceJ_Fv[i][j][k].zmom = F4*S;
					FaceJ_Fv[i][j][k].ener = F5*S;

				}
			}
		}

#pragma omp for
		for(i=1;i<=IM-1;i++)
		{
			for(j=1;j<=JM-1;j++)
			{
				for(k=1;k<=KM;k++)
				{
					uvelf = 0.5f*(FluidProp.Cell_pv[i][j][k].uvel + FluidProp.Cell_pv[i][j][k-1].uvel);
					vvelf = 0.5f*(FluidProp.Cell_pv[i][j][k].vvel + FluidProp.Cell_pv[i][j][k-1].vvel);
					wvelf = 0.5f*(FluidProp.Cell_pv[i][j][k].wvel + FluidProp.Cell_pv[i][j][k-1].wvel);
					tempf = 0.5f*(FluidProp.Cell_pv[i][j][k].temp + FluidProp.Cell_pv[i][j][k-1].temp);

					if(tempf>=10000 || tempf<=0)
					{
						//tempf = 2E2;
						cout<<"ViscousFlux K"<<endl;
					}
					mue_lam  = FluidProp.SutherLand(tempf);
					cond_lam = mue_lam*Gamma*Rgas/(Gamma-1.0f)/Pr_lam;

					if (FluidProp.turbulenceType==TurbulenceTypes::SA)
					{
						mue_tur  = (FluidProp.Cell_tur[i][j][k-1].mue + FluidProp.Cell_tur[i][j][k].mue)/2.0f;
						cond_tur = mue_tur*Gamma*Rgas / (Gamma-1.0f) / Pr_tur;					
					}
					else if(FluidProp.turbulenceType==TurbulenceTypes::BL)
					{
						mue_tur  = (FluidProp.Cell_tur[i][j][k-1].mue + FluidProp.Cell_tur[i][j][k].mue)/2.0f;
						cond_tur = mue_tur*Gamma*Rgas / (Gamma-1.0f) / Pr_tur;
					}
					else if(FluidProp.turbulenceType==TurbulenceTypes::SST)
					{

					}
					else
					{
						mue_tur  = 0.0;
						cond_tur = 0.0;
					}
					mue  = mue_lam  + mue_tur;			// if we need plot the mue, we should make Cell_vs[][][].mue
					cond = cond_lam + cond_tur;

					tau_xx = mue * (4.0f/3.0f*gradx_K[i][j][k].uvel - 2.0f/3.0f*(grady_K[i][j][k].vvel+gradz_K[i][j][k].wvel));
					tau_xy = mue * (grady_K[i][j][k].uvel + gradx_K[i][j][k].vvel);
					tau_xz = mue * (gradz_K[i][j][k].uvel + gradx_K[i][j][k].wvel);
					tau_yy = mue * (4.0f/3.0f*grady_K[i][j][k].vvel - 2.0f/3.0f*(gradx_K[i][j][k].uvel+gradz_K[i][j][k].wvel));
					tau_yz = mue * (gradz_K[i][j][k].vvel + grady_K[i][j][k].wvel);
					tau_zz = mue * (4.0f/3.0f*gradz_K[i][j][k].wvel - 2.0f/3.0f*(gradx_K[i][j][k].uvel+grady_K[i][j][k].vvel));

					F1.x = 0.0;
					F1.y = 0.0;
					F1.z = 0.0;
					F2.x = tau_xx;
					F2.y = tau_xy;
					F2.z = tau_xz;
					F3.x = tau_xy;
					F3.y = tau_yy;
					F3.z = tau_yz;
					F4.x = tau_xz;
					F4.y = tau_yz;
					F4.z = tau_zz;
					F5.x = uvelf*tau_xx + vvelf*tau_xy + wvelf*tau_xz + cond*gradx_K[i][j][k].temp;
					F5.y = uvelf*tau_xy + vvelf*tau_yy + wvelf*tau_yz + cond*grady_K[i][j][k].temp;
					F5.z = uvelf*tau_xz + vvelf*tau_yz + wvelf*tau_zz + cond*gradz_K[i][j][k].temp;

					S.x = Geometry.FaceK[i][j][k].S.x;
					S.y = Geometry.FaceK[i][j][k].S.y;
					S.z = Geometry.FaceK[i][j][k].S.z;

					FaceK_Fv[i][j][k].dens = F1*S;
					FaceK_Fv[i][j][k].xmom = F2*S;
					FaceK_Fv[i][j][k].ymom = F3*S;
					FaceK_Fv[i][j][k].zmom = F4*S;
					FaceK_Fv[i][j][k].ener = F5*S;
				}
			}
		}

#pragma omp for
		for(i=1;i<=IM-1;i++)
		{
			for(j=1;j<=JM-1;j++)		
			{
				for(k=1;k<=KM-1;k++)
				{
					Cell_Fv[i][j][k] = FaceI_Fv[i+1][j][k]-FaceI_Fv[i][j][k] + FaceJ_Fv[i][j+1][k]-FaceJ_Fv[i][j][k] + FaceK_Fv[i][j][k+1]-FaceK_Fv[i][j][k];
				}
			}
		}
	}
}
