


#include "Mesh.h"
#include <iostream>
#include <fstream>

using namespace std;

void CMesh::LU_SGS_Marching()
{
	int i, j, k, flag, nb, plane;

	REAL alpha[3];
	REAL sI, sJ, sK;
	REAL  relax_factor = CTimeDiscr::relax_factor;

	TConsVars delt_Fn_cv;
	TEddyVars delt_Fn_ev;
	TConsVars AWSI_cv, AWSJ_cv, AWSK_cv, VWSI_cv, VWSJ_cv, VWSK_cv;
	REAL	  AWSI_ev, AWSJ_ev, AWSK_ev, VWSI_ev, VWSJ_ev, VWSK_ev;
	TConsVars zero;
		
	zero.dens = 0.0;
	zero.xmom = 0.0;
	zero.ymom = 0.0;
	zero.zmom = 0.0;
	zero.ener = 0.0;

	for(nb=1;nb<=NB;nb++)
	{
#pragma omp parallel for num_threads(NT) private(i, j, k)
		for (i=1;i<=Block[nb].Geometry.IM-1;i++)
		{
			for (j=1;j<=Block[nb].Geometry.JM-1;j++)
			{
				for(k=1;k<=Block[nb].Geometry.KM-1;k++)
				{
					Block[nb].TimeDiscr.Cell_cv_old[i][j][k] = Block[nb].FluidProp.Cell_cv[i][j][k];	
					if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
					{
						Block[nb].TimeDiscr.Cell_ev_old[i][j][k] = Block[nb].FluidProp.Cell_ev[i][j][k];
					}
				}
			}
		}
	}	
		
	if(UseForceFunction == false)		//fine grid's UseForceFunction always is false, when coarse grid's is true, we needn't compute its' ResFlux 
	{
		ResFlux();
	}
	else					//UseForceFunction == true
	{
		ComputeForceFunction();
		AddForceFunction();
	}

	for(nb=1;nb<=NB;nb++)
	{
		//--------------------------------up scan
		for(plane=3;plane<=(Block[nb].Geometry.IM+Block[nb].Geometry.JM+Block[nb].Geometry.KM-3);plane++)
		{
#pragma omp parallel for num_threads(NT) default(shared) private(i,j,k,flag, alpha, AWSI_cv, AWSJ_cv, AWSK_cv, VWSI_cv,\
VWSJ_cv, VWSK_cv, sI, sJ, sK, AWSI_ev, AWSJ_ev, AWSK_ev, VWSI_ev, VWSJ_ev, VWSK_ev, delt_Fn_cv, delt_Fn_ev)
			for(i=1;i<=Block[nb].Geometry.IM-1;i++)
			{
				for(j=1;j<=Block[nb].Geometry.JM-1;j++)
				{
					k = plane-i-j;
					if(k<1 || k>Block[nb].Geometry.KM-1)
					{
						continue;
					}

					alpha[1] = Block[nb].Geometry.Cell[i][j][k].Vol/Block[nb].TimeDiscr.dt[i][j][k] + relax_factor*(Block[nb].FluidProp.lambda_cI[i][j][k] + Block[nb].FluidProp.lambda_cJ[i][j][k] + Block[nb].FluidProp.lambda_cK[i][j][k])+ 
						2.0f*(Block[nb].FluidProp.lambda_vI[i][j][k] + Block[nb].FluidProp.lambda_vJ[i][j][k] + Block[nb].FluidProp.lambda_vK[i][j][k]);

					if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
					{
						alpha[2] = Block[nb].Geometry.Cell[i][j][k].Vol/Block[nb].TimeDiscr.dt[i][j][k] + relax_factor*(Block[nb].FluidProp.lambda_cI[i][j][k] + Block[nb].FluidProp.lambda_cJ[i][j][k] + Block[nb].FluidProp.lambda_cK[i][j][k])+ 
							3.0f*(Block[nb].FluidProp.lambda_vI[i][j][k] + Block[nb].FluidProp.lambda_vJ[i][j][k] + Block[nb].FluidProp.lambda_vK[i][j][k]);
					}
					if(i != 1)
					{
						flag = 100;
						Block[nb].TimeDiscr.delt_Fn_up(i, j, k, flag, delt_Fn_cv, delt_Fn_ev, Block[nb].Geometry, Block[nb].FluidProp);
						sI = 0.5f*(Block[nb].Geometry.FaceI[i-1][j][k].s + Block[nb].Geometry.FaceI[i][j][k].s);
						AWSI_cv = 0.5f*(delt_Fn_cv*sI + relax_factor*Block[nb].FluidProp.lambda_cI[i-1][j][k]*Block[nb].TimeDiscr.delt_cv[i-1][j][k]);
						VWSI_cv = Block[nb].FluidProp.lambda_vI[i-1][j][k] * Block[nb].TimeDiscr.delt_cv[i-1][j][k];
						if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
						{
							AWSI_ev = 0.5f*(delt_Fn_ev.ev*sI + relax_factor* Block[nb].FluidProp.lambda_cI[i-1][j][k]*Block[nb].TimeDiscr.delt_ev[i-1][j][k].ev);
							VWSI_ev = Block[nb].FluidProp.lambda_vI[i-1][j][k] * Block[nb].TimeDiscr.delt_ev[i-1][j][k].ev;
						} 
					}
					else
					{
						AWSI_cv = zero;
						VWSI_cv = zero;
						if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
						{
							AWSI_ev = 0.0;
							VWSI_ev = 0.0;					
						}
					}

					if(j != 1)
					{
						flag = 200;
						Block[nb].TimeDiscr.delt_Fn_up(i, j, k, flag, delt_Fn_cv, delt_Fn_ev, Block[nb].Geometry, Block[nb].FluidProp);
						sJ = 0.5f*(Block[nb].Geometry.FaceJ[i][j-1][k].s + Block[nb].Geometry.FaceJ[i][j][k].s);
						AWSJ_cv = 0.5f*(delt_Fn_cv*sJ + relax_factor*Block[nb].FluidProp.lambda_cJ[i][j-1][k]*Block[nb].TimeDiscr.delt_cv[i][j-1][k]);					
						VWSJ_cv = Block[nb].FluidProp.lambda_vJ[i][j-1][k] * Block[nb].TimeDiscr.delt_cv[i][j-1][k];
						if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
						{
							AWSJ_ev = 0.5f*(delt_Fn_ev.ev*sJ + relax_factor*Block[nb].FluidProp.lambda_cJ[i][j-1][k]*Block[nb].TimeDiscr.delt_ev[i][j-1][k].ev);					
							VWSJ_ev = Block[nb].FluidProp.lambda_vJ[i][j-1][k] * Block[nb].TimeDiscr.delt_ev[i][j-1][k].ev;
						}
					}
					else
					{
						AWSJ_cv = zero;
						VWSJ_cv = zero;
						if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
						{
							AWSJ_ev = 0.0;
							VWSJ_ev = 0.0;						
						}
					}

					if(k != 1)
					{
						flag = 300;
						Block[nb].TimeDiscr.delt_Fn_up(i, j, k, flag, delt_Fn_cv, delt_Fn_ev, Block[nb].Geometry, Block[nb].FluidProp);
						sK = 0.5f*(Block[nb].Geometry.FaceK[i][j][k-1].s + Block[nb].Geometry.FaceK[i][j][k].s);
						AWSK_cv = 0.5f*(delt_Fn_cv*sK + relax_factor*Block[nb].FluidProp.lambda_cK[i][j][k-1]*Block[nb].TimeDiscr.delt_cv[i][j][k-1]);					
						VWSK_cv = Block[nb].FluidProp.lambda_vK[i][j][k-1] * Block[nb].TimeDiscr.delt_cv[i][j][k-1];
						if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
						{
							AWSK_ev = 0.5f*(delt_Fn_ev.ev*sK + relax_factor* Block[nb].FluidProp.lambda_cK[i][j][k-1]*Block[nb].TimeDiscr.delt_ev[i][j][k-1].ev);
							VWSK_ev = Block[nb].FluidProp.lambda_vK[i][j][k-1] * Block[nb].TimeDiscr.delt_ev[i][j][k-1].ev;
						}
					}
					else
					{
						AWSK_cv = zero;
						VWSK_cv = zero;
						if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
						{
							AWSK_ev = 0.0;
							VWSK_ev = 0.0;						
						}
					}
					Block[nb].TimeDiscr.delt_cv[i][j][k] = (-1.0f*Block[nb].SpaceDiscr.Res_cv[i][j][k] + AWSI_cv + AWSJ_cv + AWSK_cv + VWSI_cv + VWSJ_cv + VWSK_cv)/alpha[1];
					if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
						Block[nb].TimeDiscr.delt_ev[i][j][k].ev = (-1.0f*Block[nb].SpaceDiscr.Res_ev[i][j][k].ev + AWSI_ev + AWSJ_ev + AWSK_ev + VWSI_ev + VWSJ_ev + VWSK_ev)/alpha[2];
				}
			}
		}
		//------------------------------------------down scan
		for(plane=Block[nb].Geometry.IM+Block[nb].Geometry.JM+Block[nb].Geometry.KM-3;plane>=3;plane--)
		{
#pragma omp parallel for num_threads(NT) default(shared) private(i,j,k,flag, alpha, AWSI_cv, AWSJ_cv, AWSK_cv, VWSI_cv,\
	VWSJ_cv, VWSK_cv, sI, sJ, sK, AWSI_ev, AWSJ_ev, AWSK_ev, VWSI_ev, VWSJ_ev, VWSK_ev, delt_Fn_cv, delt_Fn_ev)
			for(i=Block[nb].Geometry.IM-1;i>=1;i--)
			{
				for(j=Block[nb].Geometry.JM-1;j>=1;j--)
				{
					k = plane-i-j;
					if(k<1 || k>Block[nb].Geometry.KM-1)
					{
						continue;
					}
					alpha[1] = Block[nb].Geometry.Cell[i][j][k].Vol/Block[nb].TimeDiscr.dt[i][j][k] + relax_factor*(Block[nb].FluidProp.lambda_cI[i][j][k] + Block[nb].FluidProp.lambda_cJ[i][j][k] + Block[nb].FluidProp.lambda_cK[i][j][k])+ 
						2.0f*(Block[nb].FluidProp.lambda_vI[i][j][k] + Block[nb].FluidProp.lambda_vJ[i][j][k] + Block[nb].FluidProp.lambda_vK[i][j][k]);

					if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
					{
						alpha[2] = Block[nb].Geometry.Cell[i][j][k].Vol/Block[nb].TimeDiscr.dt[i][j][k] + relax_factor*(Block[nb].FluidProp.lambda_cI[i][j][k] + Block[nb].FluidProp.lambda_cJ[i][j][k] + Block[nb].FluidProp.lambda_cK[i][j][k])+ 
							3.0f*(Block[nb].FluidProp.lambda_vI[i][j][k] + Block[nb].FluidProp.lambda_vJ[i][j][k] + Block[nb].FluidProp.lambda_vK[i][j][k]);
					}
					if(i != Block[nb].Geometry.IM-1)
					{
						flag = 100;
						Block[nb].TimeDiscr.delt_Fn_down(i, j, k, flag, delt_Fn_cv, delt_Fn_ev, Block[nb].Geometry, Block[nb].FluidProp);
						sI = 0.5f*(Block[nb].Geometry.FaceI[i+1][j][k].s + Block[nb].Geometry.FaceI[i+2][j][k].s);
						AWSI_cv = 0.5f*(delt_Fn_cv*sI - relax_factor*Block[nb].FluidProp.lambda_cI[i+1][j][k]*Block[nb].TimeDiscr.delt_cv[i+1][j][k]);					
						VWSI_cv = Block[nb].FluidProp.lambda_vI[i+1][j][k] * Block[nb].TimeDiscr.delt_cv[i+1][j][k];
						if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
						{
							AWSI_ev = 0.5f*(delt_Fn_ev.ev*sI - relax_factor* Block[nb].FluidProp.lambda_cI[i+1][j][k]*Block[nb].TimeDiscr.delt_ev[i+1][j][k].ev);
							VWSI_ev = Block[nb].FluidProp.lambda_vI[i+1][j][k] * Block[nb].TimeDiscr.delt_ev[i+1][j][k].ev;
						}
					}
					else
					{
						AWSI_cv = zero;
						VWSI_cv = zero;
						if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
						{
							AWSI_ev = 0.0;
							VWSI_ev = 0.0;					
						}
					}

					if(j != Block[nb].Geometry.JM-1)
					{
						flag = 200;
						Block[nb].TimeDiscr.delt_Fn_down(i, j, k, flag, delt_Fn_cv, delt_Fn_ev, Block[nb].Geometry, Block[nb].FluidProp);
						sJ = 0.5f*(Block[nb].Geometry.FaceJ[i][j+1][k].s + Block[nb].Geometry.FaceJ[i][j+2][k].s);
						AWSJ_cv = 0.5f*(delt_Fn_cv*sJ - relax_factor*Block[nb].FluidProp.lambda_cJ[i][j+1][k]*Block[nb].TimeDiscr.delt_cv[i][j+1][k]);					
						VWSJ_cv = Block[nb].FluidProp.lambda_vJ[i][j+1][k] * Block[nb].TimeDiscr.delt_cv[i][j+1][k];
						if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
						{
							AWSJ_ev = 0.5f*(delt_Fn_ev.ev*sJ - relax_factor* Block[nb].FluidProp.lambda_cJ[i][j+1][k]*Block[nb].TimeDiscr.delt_ev[i][j+1][k].ev);
							VWSJ_ev = Block[nb].FluidProp.lambda_vJ[i][j+1][k] * Block[nb].TimeDiscr.delt_ev[i][j+1][k].ev;
						}
					}
					else
					{
						AWSJ_cv = zero;
						VWSJ_cv = zero;
						if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
						{
							AWSJ_ev = 0.0;
							VWSJ_ev = 0.0;					
						}
					}

					if(k != Block[nb].Geometry.KM-1)
					{
						flag = 300;
						Block[nb].TimeDiscr.delt_Fn_down(i, j, k, flag, delt_Fn_cv, delt_Fn_ev, Block[nb].Geometry, Block[nb].FluidProp);
						sK = 0.5f*(Block[nb].Geometry.FaceK[i][j][k+1].s + Block[nb].Geometry.FaceK[i][j][k+2].s);
						AWSK_cv = 0.5f*(delt_Fn_cv*sK - relax_factor*Block[nb].FluidProp.lambda_cK[i][j][k+1]*Block[nb].TimeDiscr.delt_cv[i][j][k+1]);					
						VWSK_cv = Block[nb].FluidProp.lambda_vK[i][j][k+1] * Block[nb].TimeDiscr.delt_cv[i][j][k+1];
						if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
						{
							AWSK_ev = 0.5f*(delt_Fn_ev.ev*sK - relax_factor* Block[nb].FluidProp.lambda_cK[i][j][k+1]*Block[nb].TimeDiscr.delt_ev[i][j][k+1].ev);
							VWSK_ev = Block[nb].FluidProp.lambda_vK[i][j][k+1] * Block[nb].TimeDiscr.delt_ev[i][j][k+1].ev;
						}
					}
					else
					{
						AWSK_cv = zero;
						VWSK_cv = zero;
						if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
						{
							AWSK_ev = 0.0;
							VWSK_ev = 0.0;					
						}
					}

					Block[nb].TimeDiscr.delt_cv[i][j][k] = (alpha[1]*Block[nb].TimeDiscr.delt_cv[i][j][k] - AWSI_cv - AWSJ_cv - AWSK_cv + VWSI_cv + VWSJ_cv + VWSK_cv)/alpha[1];
					if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
						Block[nb].TimeDiscr.delt_ev[i][j][k].ev = (alpha[2]*Block[nb].TimeDiscr.delt_ev[i][j][k].ev - AWSI_ev - AWSJ_ev - AWSK_ev + VWSI_ev + VWSJ_ev + VWSK_ev)/alpha[2];
				}
			}
		}

		//----------------------------------------time stepping
#pragma omp parallel for num_threads(NT) private(i, j, k)
		for(i=1;i<=Block[nb].Geometry.IM-1;i++)
		{
			for(j=1;j<=Block[nb].Geometry.JM-1;j++)
			{
				for(k=1;k<=Block[nb].Geometry.KM-1;k++)
				{
					Block[nb].FluidProp.Cell_cv[i][j][k] = Block[nb].FluidProp.Cell_cv[i][j][k] + Block[nb].TimeDiscr.delt_cv[i][j][k];
					if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
						Block[nb].FluidProp.Cell_ev[i][j][k].ev = Block[nb].FluidProp.Cell_ev[i][j][k].ev + Block[nb].TimeDiscr.delt_ev[i][j][k].ev;

					//if(_finite(Block[nb].FluidProp.Cell_cv[i][j][k].dens)==0 || _finite(Block[nb].FluidProp.Cell_cv[i][j][k].xmom)==0 || _finite(Block[nb].FluidProp.Cell_cv[i][j][k].ymom)==0
					//	|| _finite(Block[nb].FluidProp.Cell_cv[i][j][k].zmom)==0 || _finite(Block[nb].FluidProp.Cell_cv[i][j][k].ener)==0)
					//{
					//	//Block[nb].FluidProp.Cell_cv[i][j][k] = Block[nb].TimeDiscr.Cell_cv_old[i][j][k];
					//	cout<<"LUSGS cv"<<endl;
					//}
					//if(_finite(Block[nb].FluidProp.Cell_ev[i][j][k].ev)==0)
					//	//Block[nb].FluidProp.Cell_ev[i][j][k].ev = Block[nb].TimeDiscr.Cell_ev_old[i][j][k].ev;
					//		cout<<"LUSGS ev"<<endl;

				}
			}
		}
	}
	PrimitiveVars();
	if(CFluidProps::limiter_vars==1)
		Limit_Vars();
}

void CTimeDiscr::delt_Fn_up(int i, int j, int k, int flag, TConsVars &delt_Fn_cv, TEddyVars &delt_Fn_ev, const CGeometry &Geometry,const CFluidProps &FluidProp)
{
	REAL dens, uvel, vvel, wvel, uvel_e, vvel_e, E;
	REAL Vn_old, Vn_new, Vn_e, press_new, press_old, ev_old, ev_new;
	TNode norm;
	TConsVars cv_old, cv_new;
	REAL Gamma = CFluidProps::Gamma;

	if(flag==100)
	{
		norm.x = (Geometry.FaceI[i][j][k].S.x/Geometry.FaceI[i][j][k].s + Geometry.FaceI[i-1][j][k].S.x/Geometry.FaceI[i-1][j][k].s)/2.0f;
		norm.y = (Geometry.FaceI[i][j][k].S.y/Geometry.FaceI[i][j][k].s + Geometry.FaceI[i-1][j][k].S.y/Geometry.FaceI[i-1][j][k].s)/2.0f;
		norm.z = (Geometry.FaceI[i][j][k].S.z/Geometry.FaceI[i][j][k].s + Geometry.FaceI[i-1][j][k].S.z/Geometry.FaceI[i-1][j][k].s)/2.0f;
		uvel_e = -1.0f*Geometry.Cell[i-1][j][k].coord.y * FluidProp.w_rot;
		vvel_e = Geometry.Cell[i-1][j][k].coord.x * FluidProp.w_rot;
		i = i-1;
	}
	if(flag==200)
	{
		norm.x = (Geometry.FaceJ[i][j][k].S.x/Geometry.FaceJ[i][j][k].s + Geometry.FaceJ[i][j-1][k].S.x/Geometry.FaceJ[i][j-1][k].s)/2.0f;
		norm.y = (Geometry.FaceJ[i][j][k].S.y/Geometry.FaceJ[i][j][k].s + Geometry.FaceJ[i][j-1][k].S.y/Geometry.FaceJ[i][j-1][k].s)/2.0f;
		norm.z = (Geometry.FaceJ[i][j][k].S.z/Geometry.FaceJ[i][j][k].s + Geometry.FaceJ[i][j-1][k].S.z/Geometry.FaceJ[i][j-1][k].s)/2.0f;
		uvel_e = -1.0f*Geometry.Cell[i][j-1][k].coord.y * FluidProp.w_rot;
		vvel_e = Geometry.Cell[i][j-1][k].coord.x * FluidProp.w_rot;
		j = j-1;
	}
	if(flag==300)
	{
		norm.x = (Geometry.FaceK[i][j][k].S.x/Geometry.FaceK[i][j][k].s + Geometry.FaceK[i][j][k-1].S.x/Geometry.FaceK[i][j][k-1].s)/2.0f;
		norm.y = (Geometry.FaceK[i][j][k].S.y/Geometry.FaceK[i][j][k].s + Geometry.FaceK[i][j][k-1].S.y/Geometry.FaceK[i][j][k-1].s)/2.0f;
		norm.z = (Geometry.FaceK[i][j][k].S.z/Geometry.FaceK[i][j][k].s + Geometry.FaceK[i][j][k-1].S.z/Geometry.FaceK[i][j][k-1].s)/2.0f;
		uvel_e = -1.0f*Geometry.Cell[i][j][k-1].coord.y * FluidProp.w_rot;
		vvel_e = Geometry.Cell[i][j][k-1].coord.x * FluidProp.w_rot;
		k = k-1;
	}

	cv_old = FluidProp.Cell_cv[i][j][k];
	cv_new = cv_old + delt_cv[i][j][k];

	dens = cv_new.dens;
	uvel = cv_new.xmom / cv_new.dens;
	vvel = cv_new.ymom / cv_new.dens;
	wvel = cv_new.zmom / cv_new.dens;
	E    = cv_new.ener / cv_new.dens;
	press_new  = (Gamma-1.0f)*dens*(E-0.5f*(uvel*uvel+vvel*vvel+wvel*wvel));
	press_old  =  FluidProp.Cell_pv[i][j][k].press;
	Vn_old = cv_old.xmom/cv_old.dens*norm.x + cv_old.ymom/cv_old.dens*norm.y + cv_old.zmom/cv_old.dens*norm.z;
	Vn_new = uvel*norm.x + vvel*norm.y + wvel*norm.z;
	Vn_e = uvel_e*norm.x + vvel_e*norm.y + 0.0f*norm.z;

	delt_Fn_cv.dens =  cv_new.dens*(Vn_new-Vn_e) - cv_old.dens*(Vn_old-Vn_e);
	delt_Fn_cv.xmom = (cv_new.xmom*(Vn_new-Vn_e) + press_new*norm.x) - (cv_old.xmom*(Vn_old-Vn_e) + press_old*norm.x);
	delt_Fn_cv.ymom = (cv_new.ymom*(Vn_new-Vn_e) + press_new*norm.y) - (cv_old.ymom*(Vn_old-Vn_e) + press_old*norm.y);
	delt_Fn_cv.zmom = (cv_new.zmom*(Vn_new-Vn_e) + press_new*norm.z) - (cv_old.zmom*(Vn_old-Vn_e) + press_old*norm.z);
	delt_Fn_cv.ener = ((cv_new.ener+press_new)*(Vn_new-Vn_e) + press_new*Vn_e) - ((cv_old.ener+press_old)*(Vn_old-Vn_e) + press_old*Vn_e);

	if(FluidProp.turbulenceType==TurbulenceTypes::SA)
	{
		ev_old = FluidProp.Cell_ev[i][j][k].ev;
		ev_new = ev_old + delt_ev[i][j][k].ev;
		delt_Fn_ev.ev = ev_new*Vn_new - ev_old*Vn_old;
	}
}

void CTimeDiscr::delt_Fn_down(int i, int j, int k, int flag, TConsVars &delt_Fn_cv, TEddyVars &delt_Fn_ev,const CGeometry &Geometry,const CFluidProps &FluidProp)
{
	REAL dens, uvel, vvel, wvel, uvel_e, vvel_e, E;
	REAL Vn_old, Vn_new, Vn_e, press_new, press_old, ev_old, ev_new;
	TNode norm;
	TConsVars cv_old, cv_new;
	REAL Gamma = CFluidProps::Gamma;

	if(flag==100)
	{
		norm.x = (Geometry.FaceI[i+2][j][k].S.x/Geometry.FaceI[i+2][j][k].s + Geometry.FaceI[i+1][j][k].S.x/Geometry.FaceI[i+1][j][k].s)/2.0f;
		norm.y = (Geometry.FaceI[i+2][j][k].S.y/Geometry.FaceI[i+2][j][k].s + Geometry.FaceI[i+1][j][k].S.y/Geometry.FaceI[i+1][j][k].s)/2.0f;
		norm.z = (Geometry.FaceI[i+2][j][k].S.z/Geometry.FaceI[i+2][j][k].s + Geometry.FaceI[i+1][j][k].S.z/Geometry.FaceI[i+1][j][k].s)/2.0f;
		uvel_e = -1.0f*Geometry.Cell[i+1][j][k].coord.y * FluidProp.w_rot;
		vvel_e = Geometry.Cell[i+1][j][k].coord.x * FluidProp.w_rot;
		i = i+1;
	}
	if(flag==200)
	{
		norm.x = (Geometry.FaceJ[i][j+2][k].S.x/Geometry.FaceJ[i][j+2][k].s + Geometry.FaceJ[i][j+1][k].S.x/Geometry.FaceJ[i][j+1][k].s)/2.0f;
		norm.y = (Geometry.FaceJ[i][j+2][k].S.y/Geometry.FaceJ[i][j+2][k].s + Geometry.FaceJ[i][j+1][k].S.y/Geometry.FaceJ[i][j+1][k].s)/2.0f;
		norm.z = (Geometry.FaceJ[i][j+2][k].S.z/Geometry.FaceJ[i][j+2][k].s + Geometry.FaceJ[i][j+1][k].S.z/Geometry.FaceJ[i][j+1][k].s)/2.0f;
		uvel_e = -1.0f*Geometry.Cell[i][j+1][k].coord.y * FluidProp.w_rot;
		vvel_e = Geometry.Cell[i][j+1][k].coord.x * FluidProp.w_rot;
		j = j+1;
	}
	if(flag==300)	// flag==300
	{
		norm.x = (Geometry.FaceK[i][j][k+2].S.x/Geometry.FaceK[i][j][k+2].s + Geometry.FaceK[i][j][k+1].S.x/Geometry.FaceK[i][j][k+1].s)/2.0f;
		norm.y = (Geometry.FaceK[i][j][k+2].S.y/Geometry.FaceK[i][j][k+2].s + Geometry.FaceK[i][j][k+1].S.y/Geometry.FaceK[i][j][k+1].s)/2.0f;
		norm.z = (Geometry.FaceK[i][j][k+2].S.z/Geometry.FaceK[i][j][k+2].s + Geometry.FaceK[i][j][k+1].S.z/Geometry.FaceK[i][j][k+1].s)/2.0f;
		uvel_e = -1.0f*Geometry.Cell[i][j][k+1].coord.y * FluidProp.w_rot;
		vvel_e = Geometry.Cell[i][j][k+1].coord.x * FluidProp.w_rot;
		k = k+1;
	}

	cv_old = FluidProp.Cell_cv[i][j][k];
	cv_new = cv_old + delt_cv[i][j][k];

	dens = cv_new.dens;
	uvel = cv_new.xmom / cv_new.dens;
	vvel = cv_new.ymom / cv_new.dens;
	wvel = cv_new.zmom / cv_new.dens;
	E    = cv_new.ener / cv_new.dens;
	press_new  = (Gamma-1.0f)*dens*(E-0.5f*(uvel*uvel+vvel*vvel+wvel*wvel));
	press_old  = FluidProp.Cell_pv[i][j][k].press;
	Vn_old = cv_old.xmom/cv_old.dens*norm.x + cv_old.ymom/cv_old.dens*norm.y + cv_old.zmom/cv_old.dens*norm.z;
	Vn_new = uvel*norm.x + vvel*norm.y + wvel*norm.z;
	Vn_e = uvel_e*norm.x + vvel_e*norm.y + 0.0f*norm.z;

	delt_Fn_cv.dens =  cv_new.dens*(Vn_new-Vn_e) - cv_old.dens*(Vn_old-Vn_e);
	delt_Fn_cv.xmom = (cv_new.xmom*(Vn_new-Vn_e) + press_new*norm.x) - (cv_old.xmom*(Vn_old-Vn_e) + press_old*norm.x);
	delt_Fn_cv.ymom = (cv_new.ymom*(Vn_new-Vn_e) + press_new*norm.y) - (cv_old.ymom*(Vn_old-Vn_e) + press_old*norm.y);
	delt_Fn_cv.zmom = (cv_new.zmom*(Vn_new-Vn_e) + press_new*norm.z) - (cv_old.zmom*(Vn_old-Vn_e) + press_old*norm.z);
	delt_Fn_cv.ener = ((cv_new.ener+press_new)*(Vn_new-Vn_e) + press_new*Vn_e) - ((cv_old.ener+press_old)*(Vn_old-Vn_e) + press_old*Vn_e);

	if(FluidProp.turbulenceType==TurbulenceTypes::SA)
	{
		ev_old = FluidProp.Cell_ev[i][j][k].ev;
		ev_new = ev_old + delt_ev[i][j][k].ev;
		delt_Fn_ev.ev = ev_new*Vn_new - ev_old*Vn_old;
	}
}
