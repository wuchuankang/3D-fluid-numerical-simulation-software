

#include "FluidProps.h"
#include <fstream>
#include "SpaceDiscr.h"

using namespace std;

void CFluidProps::SpectralRadii(const CGeometry &Geometry)
{
	int i, j, k;
	TNode S;
	REAL CS, csoun;
	REAL uvelr, vvelr, wvelr;
#pragma omp parallel for num_threads(NT) default(shared) private(i, j, k, uvelr, vvelr, wvelr, S, CS, csoun)
	for(i=1;i<=IM-1;i++)
	{
		for(j=1;j<=JM-1;j++)
		{
			for(k=1;k<=KM-1;k++)
			{
				uvelr = Cell_pv[i][j][k].uvel + w_rot*Geometry.Cell[i][j][k].coord.y;
				vvelr = Cell_pv[i][j][k].vvel - w_rot*Geometry.Cell[i][j][k].coord.x;
				wvelr = Cell_pv[i][j][k].wvel;

				csoun = SQRT(Gamma*Rgas*Cell_pv[i][j][k].temp);
				
				S  = 0.5f*(Geometry.FaceI[i][j][k].S + Geometry.FaceI[i+1][j][k].S);
				CS = csoun * SQRT(S.x*S.x+S.y*S.y+S.z*S.z);
				lambda_cI[i][j][k] = ABS(uvelr*S.x+vvelr*S.y+wvelr*S.z) + CS;

				S  = 0.5f*(Geometry.FaceJ[i][j][k].S + Geometry.FaceJ[i][j+1][k].S);
				CS = csoun * SQRT(S.x*S.x+S.y*S.y+S.z*S.z);
				lambda_cJ[i][j][k] = ABS(uvelr*S.x+vvelr*S.y+wvelr*S.z) + CS;

				S  = 0.5f*(Geometry.FaceK[i][j][k].S + Geometry.FaceK[i][j][k+1].S);
				CS = csoun * SQRT(S.x*S.x+S.y*S.y+S.z*S.z);
				lambda_cK[i][j][k] = ABS(uvelr*S.x+vvelr*S.y+wvelr*S.z) + CS;
			}
		}
	}
}

void CFluidProps::Curl(const CGeometry &Geometry)
{
	int i, j, k;

	REAL u_y, u_z, v_x, v_z, w_x, w_y;
	REAL u_kesi, u_yita, u_zita, v_kesi, v_yita, v_zita, w_kesi, w_yita, w_zita;

#pragma omp parallel for num_threads(NT) default(shared) private(i, j, k, u_y, u_z, v_x, v_z, w_x, w_y, u_kesi, u_yita, u_zita, v_kesi, v_yita, v_zita, w_kesi, w_yita, w_zita)
	for(k=1;k<=KM-1;k++)
	{
		for(j=1;j<=JM-1;j++)
		{
			for(i=1;i<=IM-1;i++)
			{
				u_kesi = 0.5f*(Cell_pv[i+1][j][k].uvel - Cell_pv[i-1][j][k].uvel);
				v_kesi = 0.5f*(Cell_pv[i+1][j][k].vvel - Cell_pv[i-1][j][k].vvel);
				w_kesi = 0.5f*(Cell_pv[i+1][j][k].wvel - Cell_pv[i-1][j][k].wvel);
					   	 
				u_yita = 0.5f*(Cell_pv[i][j+1][k].uvel - Cell_pv[i][j-1][k].uvel);
				v_yita = 0.5f*(Cell_pv[i][j+1][k].vvel - Cell_pv[i][j-1][k].vvel);
				w_yita = 0.5f*(Cell_pv[i][j+1][k].wvel - Cell_pv[i][j-1][k].wvel);
					   	 
				u_zita = 0.5f*(Cell_pv[i][j][k+1].uvel - Cell_pv[i][j][k-1].uvel);
				v_zita = 0.5f*(Cell_pv[i][j][k+1].vvel - Cell_pv[i][j][k-1].vvel);
				w_zita = 0.5f*(Cell_pv[i][j][k+1].wvel - Cell_pv[i][j][k-1].wvel);

				v_x = v_kesi*Geometry.Cell[i][j][k].gradx_kesi + v_yita*Geometry.Cell[i][j][k].gradx_yita + v_zita*Geometry.Cell[i][j][k].gradx_zita;
				w_x = w_kesi*Geometry.Cell[i][j][k].gradx_kesi + w_yita*Geometry.Cell[i][j][k].gradx_yita + w_zita*Geometry.Cell[i][j][k].gradx_zita;
				      										   	 										    
				u_y = u_kesi*Geometry.Cell[i][j][k].grady_kesi + u_yita*Geometry.Cell[i][j][k].grady_yita + u_zita*Geometry.Cell[i][j][k].grady_zita;
				w_y = w_kesi*Geometry.Cell[i][j][k].grady_kesi + w_yita*Geometry.Cell[i][j][k].grady_yita + w_zita*Geometry.Cell[i][j][k].grady_zita;
				      										   	 										    
				u_z = u_kesi*Geometry.Cell[i][j][k].gradz_kesi + u_yita*Geometry.Cell[i][j][k].gradz_yita + u_zita*Geometry.Cell[i][j][k].gradz_zita;
				v_z = v_kesi*Geometry.Cell[i][j][k].gradz_kesi + v_yita*Geometry.Cell[i][j][k].gradz_yita + v_zita*Geometry.Cell[i][j][k].gradz_zita;

				Omega[i][j][k] = SQRT(0.5f*((w_y-v_z)*(w_y-v_z) + (u_z-w_x)*(u_z-w_x) + (v_x-u_y)*(v_x-u_y)));
			}
		}
	}
}

void CFluidProps::Gradient(const CGeometry &Geometry)
{
	int i, j, k;
	int const IM = Geometry.IM;
	int const JM = Geometry.JM;
	int const KM = Geometry.KM;

	REAL u_kesi, v_kesi, w_kesi, T_kesi, ev_kesi,
		u_yita, v_yita, w_yita, T_yita, ev_yita, 
		u_zita, v_zita, w_zita, T_zita, ev_zita;

#pragma omp parallel num_threads(NT) default(shared) private(i, j, k, u_kesi, v_kesi, w_kesi, T_kesi, ev_kesi,\
	u_yita, v_yita, w_yita, T_yita, ev_yita, u_zita, v_zita, w_zita, T_zita, ev_zita)
	{
#pragma omp for nowait
		//I direction Viscous Flux
		for(i=1;i<=IM;i++)
		{
			for(j=1;j<=JM-1;j++)
			{
				for(k=1;k<=KM-1;k++)
				{	
					//kesi
					u_kesi = Cell_pv[i][j][k].uvel - Cell_pv[i-1][j][k].uvel;
					v_kesi = Cell_pv[i][j][k].vvel - Cell_pv[i-1][j][k].vvel;
					w_kesi = Cell_pv[i][j][k].wvel - Cell_pv[i-1][j][k].wvel;
					T_kesi = Cell_pv[i][j][k].temp - Cell_pv[i-1][j][k].temp;
					if (turbulenceType==TurbulenceTypes::SA)
						ev_kesi = Cell_ev[i][j][k].ev - Cell_ev[i-1][j][k].ev;

					if(i==1 && j==1)
					{
						u_yita = (Cell_pv[i][j][k].uvel+Cell_pv[i][j+1][k].uvel+Cell_pv[i-1][j+1][k].uvel+Cell_pv[i-1][j][k].uvel)/4.0f - (Cell_pv[i][j][k].uvel+Cell_pv[i-1][j][k].uvel+Cell_pv[i][j-1][k].uvel)/3.0f;
						v_yita = (Cell_pv[i][j][k].vvel+Cell_pv[i][j+1][k].vvel+Cell_pv[i-1][j+1][k].vvel+Cell_pv[i-1][j][k].vvel)/4.0f - (Cell_pv[i][j][k].vvel+Cell_pv[i-1][j][k].vvel+Cell_pv[i][j-1][k].vvel)/3.0f;
						w_yita = (Cell_pv[i][j][k].wvel+Cell_pv[i][j+1][k].wvel+Cell_pv[i-1][j+1][k].wvel+Cell_pv[i-1][j][k].wvel)/4.0f - (Cell_pv[i][j][k].wvel+Cell_pv[i-1][j][k].wvel+Cell_pv[i][j-1][k].wvel)/3.0f;
						T_yita = (Cell_pv[i][j][k].temp+Cell_pv[i][j+1][k].temp+Cell_pv[i-1][j+1][k].temp+Cell_pv[i-1][j][k].temp)/4.0f - (Cell_pv[i][j][k].temp+Cell_pv[i-1][j][k].temp+Cell_pv[i][j-1][k].temp)/3.0f;
						if (turbulenceType==TurbulenceTypes::SA)
							ev_yita = (Cell_ev[i][j][k].ev+Cell_ev[i][j+1][k].ev+Cell_ev[i-1][j+1][k].ev+Cell_ev[i-1][j][k].ev)/4.0f - (Cell_ev[i][j][k].ev+Cell_ev[i-1][j][k].ev+Cell_ev[i][j-1][k].ev)/3.0f;
					}
					else if(i==1 && j==(JM-1))
					{
						u_yita = (Cell_pv[i][j][k].uvel+Cell_pv[i][j+1][k].uvel+Cell_pv[i-1][j][k].uvel)/3.0f - (Cell_pv[i][j][k].uvel+Cell_pv[i-1][j][k].uvel+Cell_pv[i-1][j-1][k].uvel+Cell_pv[i][j-1][k].uvel)/4.0f;
						v_yita = (Cell_pv[i][j][k].vvel+Cell_pv[i][j+1][k].vvel+Cell_pv[i-1][j][k].vvel)/3.0f - (Cell_pv[i][j][k].vvel+Cell_pv[i-1][j][k].vvel+Cell_pv[i-1][j-1][k].vvel+Cell_pv[i][j-1][k].vvel)/4.0f;
						w_yita = (Cell_pv[i][j][k].wvel+Cell_pv[i][j+1][k].wvel+Cell_pv[i-1][j][k].wvel)/3.0f - (Cell_pv[i][j][k].wvel+Cell_pv[i-1][j][k].wvel+Cell_pv[i-1][j-1][k].wvel+Cell_pv[i][j-1][k].wvel)/4.0f;
						T_yita = (Cell_pv[i][j][k].temp+Cell_pv[i][j+1][k].temp+Cell_pv[i-1][j][k].temp)/3.0f - (Cell_pv[i][j][k].temp+Cell_pv[i-1][j][k].temp+Cell_pv[i-1][j-1][k].temp+Cell_pv[i][j-1][k].temp)/4.0f;
						if (turbulenceType==TurbulenceTypes::SA)
							ev_yita=(Cell_ev[i][j][k].ev+Cell_ev[i][j+1][k].ev+Cell_ev[i-1][j][k].ev)/3.0f-0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i-1][j][k].ev+Cell_ev[i-1][j-1][k].ev+Cell_ev[i][j-1][k].ev);
					}
					else if(i==IM && j==1)
					{
						u_yita = (Cell_pv[i][j][k].uvel+Cell_pv[i][j+1][k].uvel+Cell_pv[i-1][j+1][k].uvel+Cell_pv[i-1][j][k].uvel)/4.0f - (Cell_pv[i][j][k].uvel+Cell_pv[i-1][j][k].uvel+Cell_pv[i-1][j-1][k].uvel)/3.0f;
						v_yita = (Cell_pv[i][j][k].vvel+Cell_pv[i][j+1][k].vvel+Cell_pv[i-1][j+1][k].vvel+Cell_pv[i-1][j][k].vvel)/4.0f - (Cell_pv[i][j][k].vvel+Cell_pv[i-1][j][k].vvel+Cell_pv[i-1][j-1][k].vvel)/3.0f;
						w_yita = (Cell_pv[i][j][k].wvel+Cell_pv[i][j+1][k].wvel+Cell_pv[i-1][j+1][k].wvel+Cell_pv[i-1][j][k].wvel)/4.0f - (Cell_pv[i][j][k].wvel+Cell_pv[i-1][j][k].wvel+Cell_pv[i-1][j-1][k].wvel)/3.0f;
						T_yita = (Cell_pv[i][j][k].temp+Cell_pv[i][j+1][k].temp+Cell_pv[i-1][j+1][k].temp+Cell_pv[i-1][j][k].temp)/4.0f - (Cell_pv[i][j][k].temp+Cell_pv[i-1][j][k].temp+Cell_pv[i-1][j-1][k].temp)/3.0f;	
						if (turbulenceType==TurbulenceTypes::SA)
							ev_yita=0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i][j+1][k].ev+Cell_ev[i-1][j+1][k].ev+Cell_ev[i-1][j][k].ev)-(Cell_ev[i][j][k].ev+Cell_ev[i-1][j][k].ev+Cell_ev[i-1][j-1][k].ev)/3.0f;
					}
					else if(i==IM && j==(JM-1))
					{
						u_yita = (Cell_pv[i][j][k].uvel+Cell_pv[i-1][j+1][k].uvel+Cell_pv[i-1][j][k].uvel)/3.0f - (Cell_pv[i][j][k].uvel+Cell_pv[i-1][j][k].uvel+Cell_pv[i-1][j-1][k].uvel+Cell_pv[i][j-1][k].uvel)/4.0f;
						v_yita = (Cell_pv[i][j][k].vvel+Cell_pv[i-1][j+1][k].vvel+Cell_pv[i-1][j][k].vvel)/3.0f - (Cell_pv[i][j][k].vvel+Cell_pv[i-1][j][k].vvel+Cell_pv[i-1][j-1][k].vvel+Cell_pv[i][j-1][k].vvel)/4.0f;
						w_yita = (Cell_pv[i][j][k].wvel+Cell_pv[i-1][j+1][k].wvel+Cell_pv[i-1][j][k].wvel)/3.0f - (Cell_pv[i][j][k].wvel+Cell_pv[i-1][j][k].wvel+Cell_pv[i-1][j-1][k].wvel+Cell_pv[i][j-1][k].wvel)/4.0f;
						T_yita = (Cell_pv[i][j][k].temp+Cell_pv[i-1][j+1][k].temp+Cell_pv[i-1][j][k].temp)/3.0f - (Cell_pv[i][j][k].temp+Cell_pv[i-1][j][k].temp+Cell_pv[i-1][j-1][k].temp+Cell_pv[i][j-1][k].temp)/4.0f;		
						if (turbulenceType==TurbulenceTypes::SA)
							ev_yita=(Cell_ev[i][j][k].ev+Cell_ev[i-1][j+1][k].ev+Cell_ev[i-1][j][k].ev)/3.0f - (Cell_ev[i][j][k].ev+Cell_ev[i-1][j][k].ev+Cell_ev[i-1][j-1][k].ev+Cell_ev[i][j-1][k].ev)/4.0f;	
					}
					else
					{
						u_yita = (Cell_pv[i][j][k].uvel+Cell_pv[i][j+1][k].uvel+Cell_pv[i-1][j+1][k].uvel+Cell_pv[i-1][j][k].uvel)/4.0f - (Cell_pv[i][j][k].uvel+Cell_pv[i-1][j][k].uvel+Cell_pv[i-1][j-1][k].uvel+Cell_pv[i][j-1][k].uvel)/4.0f;
						v_yita = (Cell_pv[i][j][k].vvel+Cell_pv[i][j+1][k].vvel+Cell_pv[i-1][j+1][k].vvel+Cell_pv[i-1][j][k].vvel)/4.0f - (Cell_pv[i][j][k].vvel+Cell_pv[i-1][j][k].vvel+Cell_pv[i-1][j-1][k].vvel+Cell_pv[i][j-1][k].vvel)/4.0f;
						w_yita = (Cell_pv[i][j][k].wvel+Cell_pv[i][j+1][k].wvel+Cell_pv[i-1][j+1][k].wvel+Cell_pv[i-1][j][k].wvel)/4.0f - (Cell_pv[i][j][k].wvel+Cell_pv[i-1][j][k].wvel+Cell_pv[i-1][j-1][k].wvel+Cell_pv[i][j-1][k].wvel)/4.0f;
						T_yita = (Cell_pv[i][j][k].temp+Cell_pv[i][j+1][k].temp+Cell_pv[i-1][j+1][k].temp+Cell_pv[i-1][j][k].temp)/4.0f - (Cell_pv[i][j][k].temp+Cell_pv[i-1][j][k].temp+Cell_pv[i-1][j-1][k].temp+Cell_pv[i][j-1][k].temp)/4.0f;	
						if (turbulenceType==TurbulenceTypes::SA)
							ev_yita = (Cell_ev[i][j][k].ev+Cell_ev[i][j+1][k].ev+Cell_ev[i-1][j+1][k].ev+Cell_ev[i-1][j][k].ev)/4.0f - (Cell_ev[i][j][k].ev+Cell_ev[i-1][j][k].ev+Cell_ev[i-1][j-1][k].ev+Cell_ev[i][j-1][k].ev)/4.0f;
					}
					//zita
					if(i==1 && k==1)
					{
						u_zita = (Cell_pv[i][j][k].uvel+Cell_pv[i][j][k+1].uvel+Cell_pv[i-1][j][k+1].uvel+Cell_pv[i-1][j][k].uvel)/4.0f - (Cell_pv[i][j][k].uvel+Cell_pv[i-1][j][k].uvel+Cell_pv[i][j][k-1].uvel)/3.0f;
						v_zita = (Cell_pv[i][j][k].vvel+Cell_pv[i][j][k+1].vvel+Cell_pv[i-1][j][k+1].vvel+Cell_pv[i-1][j][k].vvel)/4.0f - (Cell_pv[i][j][k].vvel+Cell_pv[i-1][j][k].vvel+Cell_pv[i][j][k-1].vvel)/3.0f;
						w_zita = (Cell_pv[i][j][k].wvel+Cell_pv[i][j][k+1].wvel+Cell_pv[i-1][j][k+1].wvel+Cell_pv[i-1][j][k].wvel)/4.0f - (Cell_pv[i][j][k].wvel+Cell_pv[i-1][j][k].wvel+Cell_pv[i][j][k-1].wvel)/3.0f;
						T_zita = (Cell_pv[i][j][k].temp+Cell_pv[i][j][k+1].temp+Cell_pv[i-1][j][k+1].temp+Cell_pv[i-1][j][k].temp)/4.0f - (Cell_pv[i][j][k].temp+Cell_pv[i-1][j][k].temp+Cell_pv[i][j][k-1].temp)/3.0f;
						if (turbulenceType==TurbulenceTypes::SA)
							ev_zita=(Cell_ev[i][j][k].ev+Cell_ev[i][j][k+1].ev+Cell_ev[i-1][j][k+1].ev+Cell_ev[i-1][j][k].ev)/4.0f - (Cell_ev[i][j][k].ev+Cell_ev[i-1][j][k].ev+Cell_ev[i][j][k-1].ev)/3.0f;
					}
					else if(i==1 && k==(KM-1))
					{
						u_zita = (Cell_pv[i][j][k].uvel+Cell_pv[i][j][k+1].uvel+Cell_pv[i-1][j][k].uvel)/3.0f - (Cell_pv[i][j][k].uvel+Cell_pv[i-1][j][k].uvel+Cell_pv[i-1][j][k-1].uvel+Cell_pv[i][j][k-1].uvel)/4.0f;
						v_zita = (Cell_pv[i][j][k].vvel+Cell_pv[i][j][k+1].vvel+Cell_pv[i-1][j][k].vvel)/3.0f - (Cell_pv[i][j][k].vvel+Cell_pv[i-1][j][k].vvel+Cell_pv[i-1][j][k-1].vvel+Cell_pv[i][j][k-1].vvel)/4.0f;
						w_zita = (Cell_pv[i][j][k].wvel+Cell_pv[i][j][k+1].wvel+Cell_pv[i-1][j][k].wvel)/3.0f - (Cell_pv[i][j][k].wvel+Cell_pv[i-1][j][k].wvel+Cell_pv[i-1][j][k-1].wvel+Cell_pv[i][j][k-1].wvel)/4.0f;
						T_zita = (Cell_pv[i][j][k].temp+Cell_pv[i][j][k+1].temp+Cell_pv[i-1][j][k].temp)/3.0f - (Cell_pv[i][j][k].temp+Cell_pv[i-1][j][k].temp+Cell_pv[i-1][j][k-1].temp+Cell_pv[i][j][k-1].temp)/4.0f;
						if (turbulenceType==TurbulenceTypes::SA)
							ev_zita=(Cell_ev[i][j][k].ev+Cell_ev[i][j][k+1].ev+Cell_ev[i-1][j][k].ev)/3.0f-0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i-1][j][k].ev+Cell_ev[i-1][j][k-1].ev+Cell_ev[i][j][k-1].ev);
					}
					else if(i==IM && k==1)
					{
						u_zita = (Cell_pv[i][j][k].uvel+Cell_pv[i][j][k+1].uvel+Cell_pv[i-1][j][k+1].uvel+Cell_pv[i-1][j][k].uvel)/4.0f - (Cell_pv[i][j][k].uvel+Cell_pv[i-1][j][k].uvel+Cell_pv[i-1][j][k-1].uvel)/3.0f;
						v_zita = (Cell_pv[i][j][k].vvel+Cell_pv[i][j][k+1].vvel+Cell_pv[i-1][j][k+1].vvel+Cell_pv[i-1][j][k].vvel)/4.0f - (Cell_pv[i][j][k].vvel+Cell_pv[i-1][j][k].vvel+Cell_pv[i-1][j][k-1].vvel)/3.0f;
						w_zita = (Cell_pv[i][j][k].wvel+Cell_pv[i][j][k+1].wvel+Cell_pv[i-1][j][k+1].wvel+Cell_pv[i-1][j][k].wvel)/4.0f - (Cell_pv[i][j][k].wvel+Cell_pv[i-1][j][k].wvel+Cell_pv[i-1][j][k-1].wvel)/3.0f;
						T_zita = (Cell_pv[i][j][k].temp+Cell_pv[i][j][k+1].temp+Cell_pv[i-1][j][k+1].temp+Cell_pv[i-1][j][k].temp)/4.0f - (Cell_pv[i][j][k].temp+Cell_pv[i-1][j][k].temp+Cell_pv[i-1][j][k-1].temp)/3.0f;
						if (turbulenceType==TurbulenceTypes::SA)
							ev_zita=0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i][j][k+1].ev+Cell_ev[i-1][j][k+1].ev+Cell_ev[i-1][j][k].ev)-(Cell_ev[i][j][k].ev+Cell_ev[i-1][j][k].ev+Cell_ev[i-1][j][k-1].ev)/3.0f;
					}
					else if(i==IM && k==(KM-1))
					{
						u_zita = (Cell_pv[i][j][k].uvel+Cell_pv[i-1][j][k+1].uvel+Cell_pv[i-1][j][k].uvel)/3.0f - (Cell_pv[i][j][k].uvel+Cell_pv[i-1][j][k].uvel+Cell_pv[i-1][j][k-1].uvel+Cell_pv[i][j][k-1].uvel)/4.0f;
						v_zita = (Cell_pv[i][j][k].vvel+Cell_pv[i-1][j][k+1].vvel+Cell_pv[i-1][j][k].vvel)/3.0f - (Cell_pv[i][j][k].vvel+Cell_pv[i-1][j][k].vvel+Cell_pv[i-1][j][k-1].vvel+Cell_pv[i][j][k-1].vvel)/4.0f;
						w_zita = (Cell_pv[i][j][k].wvel+Cell_pv[i-1][j][k+1].wvel+Cell_pv[i-1][j][k].wvel)/3.0f - (Cell_pv[i][j][k].wvel+Cell_pv[i-1][j][k].wvel+Cell_pv[i-1][j][k-1].wvel+Cell_pv[i][j][k-1].wvel)/4.0f;
						T_zita = (Cell_pv[i][j][k].temp+Cell_pv[i-1][j][k+1].temp+Cell_pv[i-1][j][k].temp)/3.0f - (Cell_pv[i][j][k].temp+Cell_pv[i-1][j][k].temp+Cell_pv[i-1][j][k-1].temp+Cell_pv[i][j][k-1].temp)/4.0f;
						if (turbulenceType==TurbulenceTypes::SA)
							ev_zita=(Cell_ev[i][j][k].ev+Cell_ev[i-1][j][k+1].ev+Cell_ev[i-1][j][k].ev)/3.0f-0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i-1][j][k].ev+Cell_ev[i-1][j][k-1].ev+Cell_ev[i][j][k-1].ev);
					}
					else
					{
						u_zita = (Cell_pv[i][j][k].uvel+Cell_pv[i][j][k+1].uvel+Cell_pv[i-1][j][k+1].uvel+Cell_pv[i-1][j][k].uvel)/4.0f - (Cell_pv[i][j][k].uvel+Cell_pv[i-1][j][k].uvel+Cell_pv[i-1][j][k-1].uvel+Cell_pv[i][j][k-1].uvel)/4.0f;
						v_zita = (Cell_pv[i][j][k].vvel+Cell_pv[i][j][k+1].vvel+Cell_pv[i-1][j][k+1].vvel+Cell_pv[i-1][j][k].vvel)/4.0f - (Cell_pv[i][j][k].vvel+Cell_pv[i-1][j][k].vvel+Cell_pv[i-1][j][k-1].vvel+Cell_pv[i][j][k-1].vvel)/4.0f;
						w_zita = (Cell_pv[i][j][k].wvel+Cell_pv[i][j][k+1].wvel+Cell_pv[i-1][j][k+1].wvel+Cell_pv[i-1][j][k].wvel)/4.0f - (Cell_pv[i][j][k].wvel+Cell_pv[i-1][j][k].wvel+Cell_pv[i-1][j][k-1].wvel+Cell_pv[i][j][k-1].wvel)/4.0f;
						T_zita = (Cell_pv[i][j][k].temp+Cell_pv[i][j][k+1].temp+Cell_pv[i-1][j][k+1].temp+Cell_pv[i-1][j][k].temp)/4.0f - (Cell_pv[i][j][k].temp+Cell_pv[i-1][j][k].temp+Cell_pv[i-1][j][k-1].temp+Cell_pv[i][j][k-1].temp)/4.0f;
						if (turbulenceType==TurbulenceTypes::SA)
							ev_zita=0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i][j][k+1].ev+Cell_ev[i-1][j][k+1].ev+Cell_ev[i-1][j][k].ev)-0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i-1][j][k].ev+Cell_ev[i-1][j][k-1].ev+Cell_ev[i][j][k-1].ev);
					}

					gradx_I[i][j][k].uvel = u_kesi*Geometry.FaceI[i][j][k].gradx_kesi + u_yita*Geometry.FaceI[i][j][k].gradx_yita + u_zita*Geometry.FaceI[i][j][k].gradx_zita;
					gradx_I[i][j][k].vvel = v_kesi*Geometry.FaceI[i][j][k].gradx_kesi + v_yita*Geometry.FaceI[i][j][k].gradx_yita + v_zita*Geometry.FaceI[i][j][k].gradx_zita;
					gradx_I[i][j][k].wvel = w_kesi*Geometry.FaceI[i][j][k].gradx_kesi + w_yita*Geometry.FaceI[i][j][k].gradx_yita + w_zita*Geometry.FaceI[i][j][k].gradx_zita;
					gradx_I[i][j][k].temp = T_kesi*Geometry.FaceI[i][j][k].gradx_kesi + T_yita*Geometry.FaceI[i][j][k].gradx_yita + T_zita*Geometry.FaceI[i][j][k].gradx_zita;
					grady_I[i][j][k].uvel = u_kesi*Geometry.FaceI[i][j][k].grady_kesi + u_yita*Geometry.FaceI[i][j][k].grady_yita + u_zita*Geometry.FaceI[i][j][k].grady_zita;
					grady_I[i][j][k].vvel = v_kesi*Geometry.FaceI[i][j][k].grady_kesi + v_yita*Geometry.FaceI[i][j][k].grady_yita + v_zita*Geometry.FaceI[i][j][k].grady_zita;
					grady_I[i][j][k].wvel = w_kesi*Geometry.FaceI[i][j][k].grady_kesi + w_yita*Geometry.FaceI[i][j][k].grady_yita + w_zita*Geometry.FaceI[i][j][k].grady_zita;
					grady_I[i][j][k].temp = T_kesi*Geometry.FaceI[i][j][k].grady_kesi + T_yita*Geometry.FaceI[i][j][k].grady_yita + T_zita*Geometry.FaceI[i][j][k].grady_zita;
					gradz_I[i][j][k].uvel = u_kesi*Geometry.FaceI[i][j][k].gradz_kesi + u_yita*Geometry.FaceI[i][j][k].gradz_yita + u_zita*Geometry.FaceI[i][j][k].gradz_zita;
					gradz_I[i][j][k].vvel = v_kesi*Geometry.FaceI[i][j][k].gradz_kesi + v_yita*Geometry.FaceI[i][j][k].gradz_yita + v_zita*Geometry.FaceI[i][j][k].gradz_zita;
					gradz_I[i][j][k].wvel = w_kesi*Geometry.FaceI[i][j][k].gradz_kesi + w_yita*Geometry.FaceI[i][j][k].gradz_yita + w_zita*Geometry.FaceI[i][j][k].gradz_zita;
					gradz_I[i][j][k].temp = T_kesi*Geometry.FaceI[i][j][k].gradz_kesi + T_yita*Geometry.FaceI[i][j][k].gradz_yita + T_zita*Geometry.FaceI[i][j][k].gradz_zita;
					if (turbulenceType==TurbulenceTypes::SA)
					{
						gradx_I_ev[i][j][k] = ev_kesi*Geometry.FaceI[i][j][k].gradx_kesi + ev_yita*Geometry.FaceI[i][j][k].gradx_yita + ev_zita*Geometry.FaceI[i][j][k].gradx_zita;
						grady_I_ev[i][j][k] = ev_kesi*Geometry.FaceI[i][j][k].grady_kesi + ev_yita*Geometry.FaceI[i][j][k].grady_yita + ev_zita*Geometry.FaceI[i][j][k].grady_zita;
						gradz_I_ev[i][j][k] = ev_kesi*Geometry.FaceI[i][j][k].gradz_kesi + ev_yita*Geometry.FaceI[i][j][k].gradz_yita + ev_zita*Geometry.FaceI[i][j][k].gradz_zita;
					}
				}
			}
		}

#pragma omp for nowait
		// J direction
		for(i=1;i<=IM-1;i++)
		{
			for(j=1;j<=JM;j++)
			{
				for(k=1;k<=KM-1;k++)
				{
					//kesai
					if(j==1 && i==1)
					{
						u_kesi=0.25f*(Cell_pv[i][j][k].uvel+Cell_pv[i+1][j][k].uvel+Cell_pv[i+1][j-1][k].uvel+Cell_pv[i][j-1][k].uvel)-(Cell_pv[i][j][k].uvel+Cell_pv[i-1][j][k].uvel+Cell_pv[i][j-1][k].uvel)/3.0f;
						v_kesi=0.25f*(Cell_pv[i][j][k].vvel+Cell_pv[i+1][j][k].vvel+Cell_pv[i+1][j-1][k].vvel+Cell_pv[i][j-1][k].vvel)-(Cell_pv[i][j][k].vvel+Cell_pv[i-1][j][k].vvel+Cell_pv[i][j-1][k].vvel)/3.0f;
						w_kesi=0.25f*(Cell_pv[i][j][k].wvel+Cell_pv[i+1][j][k].wvel+Cell_pv[i+1][j-1][k].wvel+Cell_pv[i][j-1][k].wvel)-(Cell_pv[i][j][k].wvel+Cell_pv[i-1][j][k].wvel+Cell_pv[i][j-1][k].wvel)/3.0f;
						T_kesi=0.25f*(Cell_pv[i][j][k].temp+Cell_pv[i+1][j][k].temp+Cell_pv[i+1][j-1][k].temp+Cell_pv[i][j-1][k].temp)-(Cell_pv[i][j][k].temp+Cell_pv[i-1][j][k].temp+Cell_pv[i][j-1][k].temp)/3.0f;
						if (turbulenceType==TurbulenceTypes::SA)
							ev_kesi=0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i+1][j][k].ev+Cell_ev[i+1][j-1][k].ev+Cell_ev[i][j-1][k].ev)-(Cell_ev[i][j][k].ev+Cell_ev[i-1][j][k].ev+Cell_ev[i][j-1][k].ev)/3.0f;
					}
					else if(j==1 && i==(IM-1))
					{
						u_kesi=(Cell_pv[i][j][k].uvel+Cell_pv[i+1][j][k].uvel+Cell_pv[i][j-1][k].uvel)/3.0f-0.25f*(Cell_pv[i][j][k].uvel+Cell_pv[i-1][j][k].uvel+Cell_pv[i-1][j-1][k].uvel+Cell_pv[i][j-1][k].uvel);
						v_kesi=(Cell_pv[i][j][k].vvel+Cell_pv[i+1][j][k].vvel+Cell_pv[i][j-1][k].vvel)/3.0f-0.25f*(Cell_pv[i][j][k].vvel+Cell_pv[i-1][j][k].vvel+Cell_pv[i-1][j-1][k].vvel+Cell_pv[i][j-1][k].vvel);
						w_kesi=(Cell_pv[i][j][k].wvel+Cell_pv[i+1][j][k].wvel+Cell_pv[i][j-1][k].wvel)/3.0f-0.25f*(Cell_pv[i][j][k].wvel+Cell_pv[i-1][j][k].wvel+Cell_pv[i-1][j-1][k].wvel+Cell_pv[i][j-1][k].wvel);
						T_kesi=(Cell_pv[i][j][k].temp+Cell_pv[i+1][j][k].temp+Cell_pv[i][j-1][k].temp)/3.0f-0.25f*(Cell_pv[i][j][k].temp+Cell_pv[i-1][j][k].temp+Cell_pv[i-1][j-1][k].temp+Cell_pv[i][j-1][k].temp);
						if (turbulenceType==TurbulenceTypes::SA)
							ev_kesi=(Cell_ev[i][j][k].ev+Cell_ev[i+1][j][k].ev+Cell_ev[i][j-1][k].ev)/3.0f-0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i-1][j][k].ev+Cell_ev[i-1][j-1][k].ev+Cell_ev[i][j-1][k].ev);
					}
					else if(j==JM && i==1)
					{
						u_kesi=0.25f*(Cell_pv[i][j][k].uvel+Cell_pv[i+1][j][k].uvel+Cell_pv[i+1][j-1][k].uvel+Cell_pv[i][j-1][k].uvel)-(Cell_pv[i][j][k].uvel+Cell_pv[i-1][j-1][k].uvel+Cell_pv[i][j-1][k].uvel)/3.0f;
						v_kesi=0.25f*(Cell_pv[i][j][k].vvel+Cell_pv[i+1][j][k].vvel+Cell_pv[i+1][j-1][k].vvel+Cell_pv[i][j-1][k].vvel)-(Cell_pv[i][j][k].vvel+Cell_pv[i-1][j-1][k].vvel+Cell_pv[i][j-1][k].vvel)/3.0f;
						w_kesi=0.25f*(Cell_pv[i][j][k].wvel+Cell_pv[i+1][j][k].wvel+Cell_pv[i+1][j-1][k].wvel+Cell_pv[i][j-1][k].wvel)-(Cell_pv[i][j][k].wvel+Cell_pv[i-1][j-1][k].wvel+Cell_pv[i][j-1][k].wvel)/3.0f;
						T_kesi=0.25f*(Cell_pv[i][j][k].temp+Cell_pv[i+1][j][k].temp+Cell_pv[i+1][j-1][k].temp+Cell_pv[i][j-1][k].temp)-(Cell_pv[i][j][k].temp+Cell_pv[i-1][j-1][k].temp+Cell_pv[i][j-1][k].temp)/3.0f;
						if (turbulenceType==TurbulenceTypes::SA)
							ev_kesi=0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i+1][j][k].ev+Cell_ev[i+1][j-1][k].ev+Cell_ev[i][j-1][k].ev)-(Cell_ev[i][j][k].ev+Cell_ev[i-1][j-1][k].ev+Cell_ev[i][j-1][k].ev)/3.0f;
					}
					else if(j==JM && i==(IM-1))
					{
						u_kesi=(Cell_pv[i][j][k].uvel+Cell_pv[i+1][j-1][k].uvel+Cell_pv[i][j-1][k].uvel)/3.0f-0.25f*(Cell_pv[i][j][k].uvel+Cell_pv[i-1][j][k].uvel+Cell_pv[i-1][j-1][k].uvel+Cell_pv[i][j-1][k].uvel);
						v_kesi=(Cell_pv[i][j][k].vvel+Cell_pv[i+1][j-1][k].vvel+Cell_pv[i][j-1][k].vvel)/3.0f-0.25f*(Cell_pv[i][j][k].vvel+Cell_pv[i-1][j][k].vvel+Cell_pv[i-1][j-1][k].vvel+Cell_pv[i][j-1][k].vvel);
						w_kesi=(Cell_pv[i][j][k].wvel+Cell_pv[i+1][j-1][k].wvel+Cell_pv[i][j-1][k].wvel)/3.0f-0.25f*(Cell_pv[i][j][k].wvel+Cell_pv[i-1][j][k].wvel+Cell_pv[i-1][j-1][k].wvel+Cell_pv[i][j-1][k].wvel);
						T_kesi=(Cell_pv[i][j][k].temp+Cell_pv[i+1][j-1][k].temp+Cell_pv[i][j-1][k].temp)/3.0f-0.25f*(Cell_pv[i][j][k].temp+Cell_pv[i-1][j][k].temp+Cell_pv[i-1][j-1][k].temp+Cell_pv[i][j-1][k].temp);
						if (turbulenceType==TurbulenceTypes::SA)
							ev_kesi=(Cell_ev[i][j][k].ev+Cell_ev[i+1][j-1][k].ev+Cell_ev[i][j-1][k].ev)/3.0f-0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i-1][j][k].ev+Cell_ev[i-1][j-1][k].ev+Cell_ev[i][j-1][k].ev);
					}
					else
					{
						u_kesi=0.25f*(Cell_pv[i][j][k].uvel+Cell_pv[i+1][j][k].uvel+Cell_pv[i+1][j-1][k].uvel+Cell_pv[i][j-1][k].uvel)-0.25f*(Cell_pv[i][j][k].uvel+Cell_pv[i-1][j][k].uvel+Cell_pv[i-1][j-1][k].uvel+Cell_pv[i][j-1][k].uvel);
						v_kesi=0.25f*(Cell_pv[i][j][k].vvel+Cell_pv[i+1][j][k].vvel+Cell_pv[i+1][j-1][k].vvel+Cell_pv[i][j-1][k].vvel)-0.25f*(Cell_pv[i][j][k].vvel+Cell_pv[i-1][j][k].vvel+Cell_pv[i-1][j-1][k].vvel+Cell_pv[i][j-1][k].vvel);
						w_kesi=0.25f*(Cell_pv[i][j][k].wvel+Cell_pv[i+1][j][k].wvel+Cell_pv[i+1][j-1][k].wvel+Cell_pv[i][j-1][k].wvel)-0.25f*(Cell_pv[i][j][k].wvel+Cell_pv[i-1][j][k].wvel+Cell_pv[i-1][j-1][k].wvel+Cell_pv[i][j-1][k].wvel);
						T_kesi=0.25f*(Cell_pv[i][j][k].temp+Cell_pv[i+1][j][k].temp+Cell_pv[i+1][j-1][k].temp+Cell_pv[i][j-1][k].temp)-0.25f*(Cell_pv[i][j][k].temp+Cell_pv[i-1][j][k].temp+Cell_pv[i-1][j-1][k].temp+Cell_pv[i][j-1][k].temp);
						if (turbulenceType==TurbulenceTypes::SA)
							ev_kesi=0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i+1][j][k].ev+Cell_ev[i+1][j-1][k].ev+Cell_ev[i][j-1][k].ev)-0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i-1][j][k].ev+Cell_ev[i-1][j-1][k].ev+Cell_ev[i][j-1][k].ev);
					}
					//yita
					u_yita = Cell_pv[i][j][k].uvel - Cell_pv[i][j-1][k].uvel;
					v_yita = Cell_pv[i][j][k].vvel - Cell_pv[i][j-1][k].vvel;
					w_yita = Cell_pv[i][j][k].wvel - Cell_pv[i][j-1][k].wvel;
					T_yita = Cell_pv[i][j][k].temp - Cell_pv[i][j-1][k].temp;
					if (turbulenceType==TurbulenceTypes::SA)
						ev_yita = Cell_ev[i][j][k].ev - Cell_ev[i][j-1][k].ev;
					//zita
					if(j==1 && k==1)
					{
						u_zita=0.25f*(Cell_pv[i][j][k].uvel+Cell_pv[i][j][k+1].uvel+Cell_pv[i][j-1][k+1].uvel+Cell_pv[i][j-1][k].uvel)-(Cell_pv[i][j][k].uvel+Cell_pv[i][j-1][k].uvel+Cell_pv[i][j][k-1].uvel)/3.0f;
						v_zita=0.25f*(Cell_pv[i][j][k].vvel+Cell_pv[i][j][k+1].vvel+Cell_pv[i][j-1][k+1].vvel+Cell_pv[i][j-1][k].vvel)-(Cell_pv[i][j][k].vvel+Cell_pv[i][j-1][k].vvel+Cell_pv[i][j][k-1].vvel)/3.0f;
						w_zita=0.25f*(Cell_pv[i][j][k].wvel+Cell_pv[i][j][k+1].wvel+Cell_pv[i][j-1][k+1].wvel+Cell_pv[i][j-1][k].wvel)-(Cell_pv[i][j][k].wvel+Cell_pv[i][j-1][k].wvel+Cell_pv[i][j][k-1].wvel)/3.0f;
						T_zita=0.25f*(Cell_pv[i][j][k].temp+Cell_pv[i][j][k+1].temp+Cell_pv[i][j-1][k+1].temp+Cell_pv[i][j-1][k].temp)-(Cell_pv[i][j][k].temp+Cell_pv[i][j-1][k].temp+Cell_pv[i][j][k-1].temp)/3.0f;
						if (turbulenceType==TurbulenceTypes::SA)
							ev_zita=0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i][j][k+1].ev+Cell_ev[i][j-1][k+1].ev+Cell_ev[i][j-1][k].ev)-(Cell_ev[i][j][k].ev+Cell_ev[i][j-1][k].ev+Cell_ev[i][j][k-1].ev)/3.0f;
					}
					else if(j==1 && k==(KM-1))
					{
						u_zita=(Cell_pv[i][j][k].uvel+Cell_pv[i][j][k+1].uvel+Cell_pv[i][j-1][k].uvel)/3.0f-0.25f*(Cell_pv[i][j][k].uvel+Cell_pv[i][j-1][k].uvel+Cell_pv[i][j-1][k-1].uvel+Cell_pv[i][j][k-1].uvel);
						v_zita=(Cell_pv[i][j][k].vvel+Cell_pv[i][j][k+1].vvel+Cell_pv[i][j-1][k].vvel)/3.0f-0.25f*(Cell_pv[i][j][k].vvel+Cell_pv[i][j-1][k].vvel+Cell_pv[i][j-1][k-1].vvel+Cell_pv[i][j][k-1].vvel);
						w_zita=(Cell_pv[i][j][k].wvel+Cell_pv[i][j][k+1].wvel+Cell_pv[i][j-1][k].wvel)/3.0f-0.25f*(Cell_pv[i][j][k].wvel+Cell_pv[i][j-1][k].wvel+Cell_pv[i][j-1][k-1].wvel+Cell_pv[i][j][k-1].wvel);
						T_zita=(Cell_pv[i][j][k].temp+Cell_pv[i][j][k+1].temp+Cell_pv[i][j-1][k].temp)/3.0f-0.25f*(Cell_pv[i][j][k].temp+Cell_pv[i][j-1][k].temp+Cell_pv[i][j-1][k-1].temp+Cell_pv[i][j][k-1].temp);
						if (turbulenceType==TurbulenceTypes::SA)
							ev_zita=(Cell_ev[i][j][k].ev+Cell_ev[i][j][k+1].ev+Cell_ev[i][j-1][k].ev)/3.0f-0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i][j-1][k].ev+Cell_ev[i][j-1][k-1].ev+Cell_ev[i][j][k-1].ev);
					}
					else if(j==JM && k==1)
					{
						u_zita=0.25f*(Cell_pv[i][j][k].uvel+Cell_pv[i][j][k+1].uvel+Cell_pv[i][j-1][k+1].uvel+Cell_pv[i][j-1][k].uvel)-0.25f*(Cell_pv[i][j][k].uvel+Cell_pv[i][j-1][k].uvel+Cell_pv[i][j-1][k-1].uvel)/3.0f;
						v_zita=0.25f*(Cell_pv[i][j][k].vvel+Cell_pv[i][j][k+1].vvel+Cell_pv[i][j-1][k+1].vvel+Cell_pv[i][j-1][k].vvel)-0.25f*(Cell_pv[i][j][k].vvel+Cell_pv[i][j-1][k].vvel+Cell_pv[i][j-1][k-1].vvel)/3.0f;
						w_zita=0.25f*(Cell_pv[i][j][k].wvel+Cell_pv[i][j][k+1].wvel+Cell_pv[i][j-1][k+1].wvel+Cell_pv[i][j-1][k].wvel)-0.25f*(Cell_pv[i][j][k].wvel+Cell_pv[i][j-1][k].wvel+Cell_pv[i][j-1][k-1].wvel)/3.0f;
						T_zita=0.25f*(Cell_pv[i][j][k].temp+Cell_pv[i][j][k+1].temp+Cell_pv[i][j-1][k+1].temp+Cell_pv[i][j-1][k].temp)-0.25f*(Cell_pv[i][j][k].temp+Cell_pv[i][j-1][k].temp+Cell_pv[i][j-1][k-1].temp)/3.0f;
						if (turbulenceType==TurbulenceTypes::SA)
							ev_zita=0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i][j][k+1].ev+Cell_ev[i][j-1][k+1].ev+Cell_ev[i][j-1][k].ev)-0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i][j-1][k].ev+Cell_ev[i][j-1][k-1].ev)/3.0f;
					}
					else if(j==JM && k==(KM-1))
					{
						u_zita=(Cell_pv[i][j][k].uvel+Cell_pv[i][j-1][k+1].uvel+Cell_pv[i][j-1][k].uvel)/3.0f-0.25f*(Cell_pv[i][j][k].uvel+Cell_pv[i][j-1][k].uvel+Cell_pv[i][j-1][k-1].uvel+Cell_pv[i][j][k-1].uvel);
						v_zita=(Cell_pv[i][j][k].vvel+Cell_pv[i][j-1][k+1].vvel+Cell_pv[i][j-1][k].vvel)/3.0f-0.25f*(Cell_pv[i][j][k].vvel+Cell_pv[i][j-1][k].vvel+Cell_pv[i][j-1][k-1].vvel+Cell_pv[i][j][k-1].vvel);
						w_zita=(Cell_pv[i][j][k].wvel+Cell_pv[i][j-1][k+1].wvel+Cell_pv[i][j-1][k].wvel)/3.0f-0.25f*(Cell_pv[i][j][k].wvel+Cell_pv[i][j-1][k].wvel+Cell_pv[i][j-1][k-1].wvel+Cell_pv[i][j][k-1].wvel);
						T_zita=(Cell_pv[i][j][k].temp+Cell_pv[i][j-1][k+1].temp+Cell_pv[i][j-1][k].temp)/3.0f-0.25f*(Cell_pv[i][j][k].temp+Cell_pv[i][j-1][k].temp+Cell_pv[i][j-1][k-1].temp+Cell_pv[i][j][k-1].temp);
						if (turbulenceType==TurbulenceTypes::SA)
							ev_zita=(Cell_ev[i][j][k].ev+Cell_ev[i][j-1][k+1].ev+Cell_ev[i][j-1][k].ev)/3.0f-0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i][j-1][k].ev+Cell_ev[i][j-1][k-1].ev+Cell_ev[i][j][k-1].ev);
					}
					else
					{				
						u_zita=0.25f*(Cell_pv[i][j][k].uvel+Cell_pv[i][j][k+1].uvel+Cell_pv[i][j-1][k+1].uvel+Cell_pv[i][j-1][k].uvel)-0.25f*(Cell_pv[i][j][k].uvel+Cell_pv[i][j-1][k].uvel+Cell_pv[i][j-1][k-1].uvel+Cell_pv[i][j][k-1].uvel);
						v_zita=0.25f*(Cell_pv[i][j][k].vvel+Cell_pv[i][j][k+1].vvel+Cell_pv[i][j-1][k+1].vvel+Cell_pv[i][j-1][k].vvel)-0.25f*(Cell_pv[i][j][k].vvel+Cell_pv[i][j-1][k].vvel+Cell_pv[i][j-1][k-1].vvel+Cell_pv[i][j][k-1].vvel);
						w_zita=0.25f*(Cell_pv[i][j][k].wvel+Cell_pv[i][j][k+1].wvel+Cell_pv[i][j-1][k+1].wvel+Cell_pv[i][j-1][k].wvel)-0.25f*(Cell_pv[i][j][k].wvel+Cell_pv[i][j-1][k].wvel+Cell_pv[i][j-1][k-1].wvel+Cell_pv[i][j][k-1].wvel);
						T_zita=0.25f*(Cell_pv[i][j][k].temp+Cell_pv[i][j][k+1].temp+Cell_pv[i][j-1][k+1].temp+Cell_pv[i][j-1][k].temp)-0.25f*(Cell_pv[i][j][k].temp+Cell_pv[i][j-1][k].temp+Cell_pv[i][j-1][k-1].temp+Cell_pv[i][j][k-1].temp);
						if (turbulenceType==TurbulenceTypes::SA)
							ev_zita=0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i][j][k+1].ev+Cell_ev[i][j-1][k+1].ev+Cell_ev[i][j-1][k].ev)-0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i][j-1][k].ev+Cell_ev[i][j-1][k-1].ev+Cell_ev[i][j][k-1].ev);
					}

					gradx_J[i][j][k].uvel = u_kesi*Geometry.FaceJ[i][j][k].gradx_kesi + u_yita*Geometry.FaceJ[i][j][k].gradx_yita + u_zita*Geometry.FaceJ[i][j][k].gradx_zita;
					gradx_J[i][j][k].vvel = v_kesi*Geometry.FaceJ[i][j][k].gradx_kesi + v_yita*Geometry.FaceJ[i][j][k].gradx_yita + v_zita*Geometry.FaceJ[i][j][k].gradx_zita;
					gradx_J[i][j][k].wvel = w_kesi*Geometry.FaceJ[i][j][k].gradx_kesi + w_yita*Geometry.FaceJ[i][j][k].gradx_yita + w_zita*Geometry.FaceJ[i][j][k].gradx_zita;
					gradx_J[i][j][k].temp = T_kesi*Geometry.FaceJ[i][j][k].gradx_kesi + T_yita*Geometry.FaceJ[i][j][k].gradx_yita + T_zita*Geometry.FaceJ[i][j][k].gradx_zita;
					grady_J[i][j][k].uvel = u_kesi*Geometry.FaceJ[i][j][k].grady_kesi + u_yita*Geometry.FaceJ[i][j][k].grady_yita + u_zita*Geometry.FaceJ[i][j][k].grady_zita;
					grady_J[i][j][k].vvel = v_kesi*Geometry.FaceJ[i][j][k].grady_kesi + v_yita*Geometry.FaceJ[i][j][k].grady_yita + v_zita*Geometry.FaceJ[i][j][k].grady_zita;
					grady_J[i][j][k].wvel = w_kesi*Geometry.FaceJ[i][j][k].grady_kesi + w_yita*Geometry.FaceJ[i][j][k].grady_yita + w_zita*Geometry.FaceJ[i][j][k].grady_zita;
					grady_J[i][j][k].temp = T_kesi*Geometry.FaceJ[i][j][k].grady_kesi + T_yita*Geometry.FaceJ[i][j][k].grady_yita + T_zita*Geometry.FaceJ[i][j][k].grady_zita;
					gradz_J[i][j][k].uvel = u_kesi*Geometry.FaceJ[i][j][k].gradz_kesi + u_yita*Geometry.FaceJ[i][j][k].gradz_yita + u_zita*Geometry.FaceJ[i][j][k].gradz_zita;
					gradz_J[i][j][k].vvel = v_kesi*Geometry.FaceJ[i][j][k].gradz_kesi + v_yita*Geometry.FaceJ[i][j][k].gradz_yita + v_zita*Geometry.FaceJ[i][j][k].gradz_zita;
					gradz_J[i][j][k].wvel = w_kesi*Geometry.FaceJ[i][j][k].gradz_kesi + w_yita*Geometry.FaceJ[i][j][k].gradz_yita + w_zita*Geometry.FaceJ[i][j][k].gradz_zita;
					gradz_J[i][j][k].temp = T_kesi*Geometry.FaceJ[i][j][k].gradz_kesi + T_yita*Geometry.FaceJ[i][j][k].gradz_yita + T_zita*Geometry.FaceJ[i][j][k].gradz_zita;

					if (turbulenceType==TurbulenceTypes::SA)
					{
						gradx_J_ev[i][j][k] = ev_kesi*Geometry.FaceJ[i][j][k].gradx_kesi+ev_yita*Geometry.FaceJ[i][j][k].gradx_yita+ev_zita*Geometry.FaceJ[i][j][k].gradx_zita;
						grady_J_ev[i][j][k] = ev_kesi*Geometry.FaceJ[i][j][k].grady_kesi+ev_yita*Geometry.FaceJ[i][j][k].grady_yita+ev_zita*Geometry.FaceJ[i][j][k].grady_zita;
						gradz_J_ev[i][j][k] = ev_kesi*Geometry.FaceJ[i][j][k].gradz_kesi+ev_yita*Geometry.FaceJ[i][j][k].gradz_yita+ev_zita*Geometry.FaceJ[i][j][k].gradz_zita;
					}

				}
			}
		}

#pragma omp for nowait
		// K direction
		for(i=1;i<=IM-1;i++)
		{
			for(j=1;j<=JM-1;j++)
			{
				for(k=1;k<=KM;k++)
				{
					//kesi
					if(k==1 && i==1)
					{
						u_kesi=0.25f*(Cell_pv[i][j][k].uvel+Cell_pv[i+1][j][k].uvel+Cell_pv[i+1][j][k-1].uvel+Cell_pv[i][j][k-1].uvel)-(Cell_pv[i][j][k].uvel+Cell_pv[i-1][j][k].uvel+Cell_pv[i][j][k-1].uvel)/3.0f;
						v_kesi=0.25f*(Cell_pv[i][j][k].vvel+Cell_pv[i+1][j][k].vvel+Cell_pv[i+1][j][k-1].vvel+Cell_pv[i][j][k-1].vvel)-(Cell_pv[i][j][k].vvel+Cell_pv[i-1][j][k].vvel+Cell_pv[i][j][k-1].vvel)/3.0f;
						w_kesi=0.25f*(Cell_pv[i][j][k].wvel+Cell_pv[i+1][j][k].wvel+Cell_pv[i+1][j][k-1].wvel+Cell_pv[i][j][k-1].wvel)-(Cell_pv[i][j][k].wvel+Cell_pv[i-1][j][k].wvel+Cell_pv[i][j][k-1].wvel)/3.0f;
						T_kesi=0.25f*(Cell_pv[i][j][k].temp+Cell_pv[i+1][j][k].temp+Cell_pv[i+1][j][k-1].temp+Cell_pv[i][j][k-1].temp)-(Cell_pv[i][j][k].temp+Cell_pv[i-1][j][k].temp+Cell_pv[i][j][k-1].temp)/3.0f;
						if (turbulenceType==TurbulenceTypes::SA)
							ev_kesi=0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i+1][j][k].ev+Cell_ev[i+1][j][k-1].ev+Cell_ev[i][j][k-1].ev)-(Cell_ev[i][j][k].ev+Cell_ev[i-1][j][k].ev+Cell_ev[i][j][k-1].ev)/3.0f;
					}
					else if(k==1 && i==(IM-1))
					{
						u_kesi=(Cell_pv[i][j][k].uvel+Cell_pv[i+1][j][k].uvel+Cell_pv[i][j][k-1].uvel)/3.0f-0.25f*(Cell_pv[i][j][k].uvel+Cell_pv[i-1][j][k].uvel+Cell_pv[i-1][j][k-1].uvel+Cell_pv[i][j][k-1].uvel);
						v_kesi=(Cell_pv[i][j][k].vvel+Cell_pv[i+1][j][k].vvel+Cell_pv[i][j][k-1].vvel)/3.0f-0.25f*(Cell_pv[i][j][k].vvel+Cell_pv[i-1][j][k].vvel+Cell_pv[i-1][j][k-1].vvel+Cell_pv[i][j][k-1].vvel);
						w_kesi=(Cell_pv[i][j][k].wvel+Cell_pv[i+1][j][k].wvel+Cell_pv[i][j][k-1].wvel)/3.0f-0.25f*(Cell_pv[i][j][k].wvel+Cell_pv[i-1][j][k].wvel+Cell_pv[i-1][j][k-1].wvel+Cell_pv[i][j][k-1].wvel);
						T_kesi=(Cell_pv[i][j][k].temp+Cell_pv[i+1][j][k].temp+Cell_pv[i][j][k-1].temp)/3.0f-0.25f*(Cell_pv[i][j][k].temp+Cell_pv[i-1][j][k].temp+Cell_pv[i-1][j][k-1].temp+Cell_pv[i][j][k-1].temp);
						if (turbulenceType==TurbulenceTypes::SA)
							ev_kesi=(Cell_ev[i][j][k].ev+Cell_ev[i+1][j][k].ev+Cell_ev[i][j][k-1].ev)/3.0f-0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i-1][j][k].ev+Cell_ev[i-1][j][k-1].ev+Cell_ev[i][j][k-1].ev);
					}
					else if(k==KM && i==1)
					{
						u_kesi=0.25f*(Cell_pv[i][j][k].uvel+Cell_pv[i+1][j][k].uvel+Cell_pv[i+1][j][k-1].uvel+Cell_pv[i][j][k-1].uvel)-(Cell_pv[i][j][k].uvel+Cell_pv[i-1][j][k-1].uvel+Cell_pv[i][j][k-1].uvel)/3.0f;
						v_kesi=0.25f*(Cell_pv[i][j][k].vvel+Cell_pv[i+1][j][k].vvel+Cell_pv[i+1][j][k-1].vvel+Cell_pv[i][j][k-1].vvel)-(Cell_pv[i][j][k].vvel+Cell_pv[i-1][j][k-1].vvel+Cell_pv[i][j][k-1].vvel)/3.0f;
						w_kesi=0.25f*(Cell_pv[i][j][k].wvel+Cell_pv[i+1][j][k].wvel+Cell_pv[i+1][j][k-1].wvel+Cell_pv[i][j][k-1].wvel)-(Cell_pv[i][j][k].wvel+Cell_pv[i-1][j][k-1].wvel+Cell_pv[i][j][k-1].wvel)/3.0f;
						T_kesi=0.25f*(Cell_pv[i][j][k].temp+Cell_pv[i+1][j][k].temp+Cell_pv[i+1][j][k-1].temp+Cell_pv[i][j][k-1].temp)-(Cell_pv[i][j][k].temp+Cell_pv[i-1][j][k-1].temp+Cell_pv[i][j][k-1].temp)/3.0f;
						if (turbulenceType==TurbulenceTypes::SA)
							ev_kesi=0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i+1][j][k].ev+Cell_ev[i+1][j][k-1].ev+Cell_ev[i][j][k-1].ev)-(Cell_ev[i][j][k].ev+Cell_ev[i-1][j][k-1].ev+Cell_ev[i][j][k-1].ev)/3.0f;
					}
					else if(k==KM && i==(IM-1))
					{
						u_kesi=(Cell_pv[i][j][k].uvel+Cell_pv[i+1][j][k-1].uvel+Cell_pv[i][j][k-1].uvel)/3.0f-0.25f*(Cell_pv[i][j][k].uvel+Cell_pv[i-1][j][k].uvel+Cell_pv[i-1][j][k-1].uvel+Cell_pv[i][j][k-1].uvel);
						v_kesi=(Cell_pv[i][j][k].vvel+Cell_pv[i+1][j][k-1].vvel+Cell_pv[i][j][k-1].vvel)/3.0f-0.25f*(Cell_pv[i][j][k].vvel+Cell_pv[i-1][j][k].vvel+Cell_pv[i-1][j][k-1].vvel+Cell_pv[i][j][k-1].vvel);
						w_kesi=(Cell_pv[i][j][k].wvel+Cell_pv[i+1][j][k-1].wvel+Cell_pv[i][j][k-1].wvel)/3.0f-0.25f*(Cell_pv[i][j][k].wvel+Cell_pv[i-1][j][k].wvel+Cell_pv[i-1][j][k-1].wvel+Cell_pv[i][j][k-1].wvel);
						T_kesi=(Cell_pv[i][j][k].temp+Cell_pv[i+1][j][k-1].temp+Cell_pv[i][j][k-1].temp)/3.0f-0.25f*(Cell_pv[i][j][k].temp+Cell_pv[i-1][j][k].temp+Cell_pv[i-1][j][k-1].temp+Cell_pv[i][j][k-1].temp);
						if (turbulenceType==TurbulenceTypes::SA)
							ev_kesi=(Cell_ev[i][j][k].ev+Cell_ev[i+1][j][k-1].ev+Cell_ev[i][j][k-1].ev)/3.0f-0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i-1][j][k].ev+Cell_ev[i-1][j][k-1].ev+Cell_ev[i][j][k-1].ev);
					}
					else
					{
						u_kesi=0.25f*(Cell_pv[i][j][k].uvel+Cell_pv[i+1][j][k].uvel+Cell_pv[i+1][j][k-1].uvel+Cell_pv[i][j][k-1].uvel)-0.25f*(Cell_pv[i][j][k].uvel+Cell_pv[i-1][j][k].uvel+Cell_pv[i-1][j][k-1].uvel+Cell_pv[i][j][k-1].uvel);
						v_kesi=0.25f*(Cell_pv[i][j][k].vvel+Cell_pv[i+1][j][k].vvel+Cell_pv[i+1][j][k-1].vvel+Cell_pv[i][j][k-1].vvel)-0.25f*(Cell_pv[i][j][k].vvel+Cell_pv[i-1][j][k].vvel+Cell_pv[i-1][j][k-1].vvel+Cell_pv[i][j][k-1].vvel);
						w_kesi=0.25f*(Cell_pv[i][j][k].wvel+Cell_pv[i+1][j][k].wvel+Cell_pv[i+1][j][k-1].wvel+Cell_pv[i][j][k-1].wvel)-0.25f*(Cell_pv[i][j][k].wvel+Cell_pv[i-1][j][k].wvel+Cell_pv[i-1][j][k-1].wvel+Cell_pv[i][j][k-1].wvel);
						T_kesi=0.25f*(Cell_pv[i][j][k].temp+Cell_pv[i+1][j][k].temp+Cell_pv[i+1][j][k-1].temp+Cell_pv[i][j][k-1].temp)-0.25f*(Cell_pv[i][j][k].temp+Cell_pv[i-1][j][k].temp+Cell_pv[i-1][j][k-1].temp+Cell_pv[i][j][k-1].temp);
						if (turbulenceType==TurbulenceTypes::SA)
							ev_kesi=0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i+1][j][k].ev+Cell_ev[i+1][j][k-1].ev+Cell_ev[i][j][k-1].ev)-0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i-1][j][k].ev+Cell_ev[i-1][j][k-1].ev+Cell_ev[i][j][k-1].ev);
					}
					//yita
					if(k==1 && j==1)
					{
						u_yita=0.25f*(Cell_pv[i][j][k].uvel+Cell_pv[i][j+1][k].uvel+Cell_pv[i][j+1][k-1].uvel+Cell_pv[i][j][k-1].uvel)-(Cell_pv[i][j][k].uvel+Cell_pv[i][j-1][k].uvel+Cell_pv[i][j][k-1].uvel)/3.0f;
						v_yita=0.25f*(Cell_pv[i][j][k].vvel+Cell_pv[i][j+1][k].vvel+Cell_pv[i][j+1][k-1].vvel+Cell_pv[i][j][k-1].vvel)-(Cell_pv[i][j][k].vvel+Cell_pv[i][j-1][k].vvel+Cell_pv[i][j][k-1].vvel)/3.0f;
						w_yita=0.25f*(Cell_pv[i][j][k].wvel+Cell_pv[i][j+1][k].wvel+Cell_pv[i][j+1][k-1].wvel+Cell_pv[i][j][k-1].wvel)-(Cell_pv[i][j][k].wvel+Cell_pv[i][j-1][k].wvel+Cell_pv[i][j][k-1].wvel)/3.0f;
						T_yita=0.25f*(Cell_pv[i][j][k].temp+Cell_pv[i][j+1][k].temp+Cell_pv[i][j+1][k-1].temp+Cell_pv[i][j][k-1].temp)-(Cell_pv[i][j][k].temp+Cell_pv[i][j-1][k].temp+Cell_pv[i][j][k-1].temp)/3.0f;
						if (turbulenceType==TurbulenceTypes::SA)
							ev_yita=0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i][j+1][k].ev+Cell_ev[i][j+1][k-1].ev+Cell_ev[i][j][k-1].ev)-(Cell_ev[i][j][k].ev+Cell_ev[i][j-1][k].ev+Cell_ev[i][j][k-1].ev)/3.0f;
					}
					else if(k==1 && j==(JM-1))
					{
						u_yita=(Cell_pv[i][j][k].uvel+Cell_pv[i][j+1][k].uvel+Cell_pv[i][j][k-1].uvel)/3.0f-0.25f*(Cell_pv[i][j][k].uvel+Cell_pv[i][j-1][k].uvel+Cell_pv[i][j-1][k-1].uvel+Cell_pv[i][j][k-1].uvel);
						v_yita=(Cell_pv[i][j][k].vvel+Cell_pv[i][j+1][k].vvel+Cell_pv[i][j][k-1].vvel)/3.0f-0.25f*(Cell_pv[i][j][k].vvel+Cell_pv[i][j-1][k].vvel+Cell_pv[i][j-1][k-1].vvel+Cell_pv[i][j][k-1].vvel);
						w_yita=(Cell_pv[i][j][k].wvel+Cell_pv[i][j+1][k].wvel+Cell_pv[i][j][k-1].wvel)/3.0f-0.25f*(Cell_pv[i][j][k].wvel+Cell_pv[i][j-1][k].wvel+Cell_pv[i][j-1][k-1].wvel+Cell_pv[i][j][k-1].wvel);
						T_yita=(Cell_pv[i][j][k].temp+Cell_pv[i][j+1][k].temp+Cell_pv[i][j][k-1].temp)/3.0f-0.25f*(Cell_pv[i][j][k].temp+Cell_pv[i][j-1][k].temp+Cell_pv[i][j-1][k-1].temp+Cell_pv[i][j][k-1].temp);
						if (turbulenceType==TurbulenceTypes::SA)
							ev_yita=(Cell_ev[i][j][k].ev+Cell_ev[i][j+1][k].ev+Cell_ev[i][j][k-1].ev)/3.0f-0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i][j-1][k].ev+Cell_ev[i][j-1][k-1].ev+Cell_ev[i][j][k-1].ev);
					}
					else if(k==KM && j==1)
					{
						u_yita=0.25f*(Cell_pv[i][j][k].uvel+Cell_pv[i][j+1][k].uvel+Cell_pv[i][j+1][k-1].uvel+Cell_pv[i][j][k-1].uvel)-(Cell_pv[i][j][k].uvel+Cell_pv[i][j-1][k-1].uvel+Cell_pv[i][j][k-1].uvel)/3.0f;
						v_yita=0.25f*(Cell_pv[i][j][k].vvel+Cell_pv[i][j+1][k].vvel+Cell_pv[i][j+1][k-1].vvel+Cell_pv[i][j][k-1].vvel)-(Cell_pv[i][j][k].vvel+Cell_pv[i][j-1][k-1].vvel+Cell_pv[i][j][k-1].vvel)/3.0f;
						w_yita=0.25f*(Cell_pv[i][j][k].wvel+Cell_pv[i][j+1][k].wvel+Cell_pv[i][j+1][k-1].wvel+Cell_pv[i][j][k-1].wvel)-(Cell_pv[i][j][k].wvel+Cell_pv[i][j-1][k-1].wvel+Cell_pv[i][j][k-1].wvel)/3.0f;
						T_yita=0.25f*(Cell_pv[i][j][k].temp+Cell_pv[i][j+1][k].temp+Cell_pv[i][j+1][k-1].temp+Cell_pv[i][j][k-1].temp)-(Cell_pv[i][j][k].temp+Cell_pv[i][j-1][k-1].temp+Cell_pv[i][j][k-1].temp)/3.0f;
						if (turbulenceType==TurbulenceTypes::SA)
							ev_yita=0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i][j+1][k].ev+Cell_ev[i][j+1][k-1].ev+Cell_ev[i][j][k-1].ev)-(Cell_ev[i][j][k].ev+Cell_ev[i][j-1][k-1].ev+Cell_ev[i][j][k-1].ev)/3.0f;
					}
					else if(k==KM && j==(JM-1))
					{
						u_yita=(Cell_pv[i][j][k].uvel+Cell_pv[i][j+1][k-1].uvel+Cell_pv[i][j][k-1].uvel)/3.0f-0.25f*(Cell_pv[i][j][k].uvel+Cell_pv[i][j-1][k].uvel+Cell_pv[i][j-1][k-1].uvel+Cell_pv[i][j][k-1].uvel);
						v_yita=(Cell_pv[i][j][k].vvel+Cell_pv[i][j+1][k-1].vvel+Cell_pv[i][j][k-1].vvel)/3.0f-0.25f*(Cell_pv[i][j][k].vvel+Cell_pv[i][j-1][k].vvel+Cell_pv[i][j-1][k-1].vvel+Cell_pv[i][j][k-1].vvel);
						w_yita=(Cell_pv[i][j][k].wvel+Cell_pv[i][j+1][k-1].wvel+Cell_pv[i][j][k-1].wvel)/3.0f-0.25f*(Cell_pv[i][j][k].wvel+Cell_pv[i][j-1][k].wvel+Cell_pv[i][j-1][k-1].wvel+Cell_pv[i][j][k-1].wvel);
						T_yita=(Cell_pv[i][j][k].temp+Cell_pv[i][j+1][k-1].temp+Cell_pv[i][j][k-1].temp)/3.0f-0.25f*(Cell_pv[i][j][k].temp+Cell_pv[i][j-1][k].temp+Cell_pv[i][j-1][k-1].temp+Cell_pv[i][j][k-1].temp);
						if (turbulenceType==TurbulenceTypes::SA)
							ev_yita=(Cell_ev[i][j][k].ev+Cell_ev[i][j+1][k-1].ev+Cell_ev[i][j][k-1].ev)/3.0f-0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i][j-1][k].ev+Cell_ev[i][j-1][k-1].ev+Cell_ev[i][j][k-1].ev);
					}
					else
					{
						u_yita=0.25f*(Cell_pv[i][j][k].uvel+Cell_pv[i][j+1][k].uvel+Cell_pv[i][j+1][k-1].uvel+Cell_pv[i][j][k-1].uvel)-0.25f*(Cell_pv[i][j][k].uvel+Cell_pv[i][j-1][k].uvel+Cell_pv[i][j-1][k-1].uvel+Cell_pv[i][j][k-1].uvel);
						v_yita=0.25f*(Cell_pv[i][j][k].vvel+Cell_pv[i][j+1][k].vvel+Cell_pv[i][j+1][k-1].vvel+Cell_pv[i][j][k-1].vvel)-0.25f*(Cell_pv[i][j][k].vvel+Cell_pv[i][j-1][k].vvel+Cell_pv[i][j-1][k-1].vvel+Cell_pv[i][j][k-1].vvel);
						w_yita=0.25f*(Cell_pv[i][j][k].wvel+Cell_pv[i][j+1][k].wvel+Cell_pv[i][j+1][k-1].wvel+Cell_pv[i][j][k-1].wvel)-0.25f*(Cell_pv[i][j][k].wvel+Cell_pv[i][j-1][k].wvel+Cell_pv[i][j-1][k-1].wvel+Cell_pv[i][j][k-1].wvel);
						T_yita=0.25f*(Cell_pv[i][j][k].temp+Cell_pv[i][j+1][k].temp+Cell_pv[i][j+1][k-1].temp+Cell_pv[i][j][k-1].temp)-0.25f*(Cell_pv[i][j][k].temp+Cell_pv[i][j-1][k].temp+Cell_pv[i][j-1][k-1].temp+Cell_pv[i][j][k-1].temp);
						if (turbulenceType==TurbulenceTypes::SA)
							ev_yita=0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i][j+1][k].ev+Cell_ev[i][j+1][k-1].ev+Cell_ev[i][j][k-1].ev)-0.25f*(Cell_ev[i][j][k].ev+Cell_ev[i][j-1][k].ev+Cell_ev[i][j-1][k-1].ev+Cell_ev[i][j][k-1].ev);
					}
					u_zita=Cell_pv[i][j][k].uvel-Cell_pv[i][j][k-1].uvel;
					v_zita=Cell_pv[i][j][k].vvel-Cell_pv[i][j][k-1].vvel;
					w_zita=Cell_pv[i][j][k].wvel-Cell_pv[i][j][k-1].wvel;
					T_zita=Cell_pv[i][j][k].temp-Cell_pv[i][j][k-1].temp;
					if (turbulenceType==TurbulenceTypes::SA)
						ev_zita=Cell_ev[i][j][k].ev-Cell_ev[i][j][k-1].ev;

					gradx_K[i][j][k].uvel = u_kesi*Geometry.FaceK[i][j][k].gradx_kesi + u_yita*Geometry.FaceK[i][j][k].gradx_yita + u_zita*Geometry.FaceK[i][j][k].gradx_zita;
					gradx_K[i][j][k].vvel = v_kesi*Geometry.FaceK[i][j][k].gradx_kesi + v_yita*Geometry.FaceK[i][j][k].gradx_yita + v_zita*Geometry.FaceK[i][j][k].gradx_zita;
					gradx_K[i][j][k].wvel = w_kesi*Geometry.FaceK[i][j][k].gradx_kesi + w_yita*Geometry.FaceK[i][j][k].gradx_yita + w_zita*Geometry.FaceK[i][j][k].gradx_zita;
					gradx_K[i][j][k].temp = T_kesi*Geometry.FaceK[i][j][k].gradx_kesi + T_yita*Geometry.FaceK[i][j][k].gradx_yita + T_zita*Geometry.FaceK[i][j][k].gradx_zita;
					grady_K[i][j][k].uvel = u_kesi*Geometry.FaceK[i][j][k].grady_kesi + u_yita*Geometry.FaceK[i][j][k].grady_yita + u_zita*Geometry.FaceK[i][j][k].grady_zita;
					grady_K[i][j][k].vvel = v_kesi*Geometry.FaceK[i][j][k].grady_kesi + v_yita*Geometry.FaceK[i][j][k].grady_yita + v_zita*Geometry.FaceK[i][j][k].grady_zita;
					grady_K[i][j][k].wvel = w_kesi*Geometry.FaceK[i][j][k].grady_kesi + w_yita*Geometry.FaceK[i][j][k].grady_yita + w_zita*Geometry.FaceK[i][j][k].grady_zita;
					grady_K[i][j][k].temp = T_kesi*Geometry.FaceK[i][j][k].grady_kesi + T_yita*Geometry.FaceK[i][j][k].grady_yita + T_zita*Geometry.FaceK[i][j][k].grady_zita;
					gradz_K[i][j][k].uvel = u_kesi*Geometry.FaceK[i][j][k].gradz_kesi + u_yita*Geometry.FaceK[i][j][k].gradz_yita + u_zita*Geometry.FaceK[i][j][k].gradz_zita;
					gradz_K[i][j][k].vvel = v_kesi*Geometry.FaceK[i][j][k].gradz_kesi + v_yita*Geometry.FaceK[i][j][k].gradz_yita + v_zita*Geometry.FaceK[i][j][k].gradz_zita;
					gradz_K[i][j][k].wvel = w_kesi*Geometry.FaceK[i][j][k].gradz_kesi + w_yita*Geometry.FaceK[i][j][k].gradz_yita + w_zita*Geometry.FaceK[i][j][k].gradz_zita;
					gradz_K[i][j][k].temp = T_kesi*Geometry.FaceK[i][j][k].gradz_kesi + T_yita*Geometry.FaceK[i][j][k].gradz_yita + T_zita*Geometry.FaceK[i][j][k].gradz_zita;

					if (turbulenceType==TurbulenceTypes::SA)
					{
						gradx_K_ev[i][j][k] = ev_kesi*Geometry.FaceK[i][j][k].gradx_kesi + ev_yita*Geometry.FaceK[i][j][k].gradx_yita + ev_zita*Geometry.FaceK[i][j][k].gradx_zita;
						grady_K_ev[i][j][k] = ev_kesi*Geometry.FaceK[i][j][k].grady_kesi + ev_yita*Geometry.FaceK[i][j][k].grady_yita + ev_zita*Geometry.FaceK[i][j][k].grady_zita;
						gradz_K_ev[i][j][k] = ev_kesi*Geometry.FaceK[i][j][k].gradz_kesi + ev_yita*Geometry.FaceK[i][j][k].gradz_yita + ev_zita*Geometry.FaceK[i][j][k].gradz_zita;
					}

				}
			}
		}
	}
}

