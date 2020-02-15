
#include <iostream>
#include <fstream>
#include "FluidProps.h"
#include "Output.h"

using namespace std;

REAL COutput::MassFlowCaculation(int K, CBlocks &Block)
{
	int i, j, k;
	REAL MassFlow, uvel, vvel, wvel, dens;
	k = K;
	MassFlow = 0.0;
	for(i=1;i<=Block.Geometry.IM-1;i++)
	{
		for(j=1;j<=Block.Geometry.JM-1;j++)
		{
			dens = 0.5f*(Block.FluidProp.Cell_pv[i][j][k].dens+Block.FluidProp.Cell_pv[i][j][k-1].dens);
			uvel = 0.5f*(Block.FluidProp.Cell_pv[i][j][k].uvel+Block.FluidProp.Cell_pv[i][j][k-1].uvel);
			vvel = 0.5f*(Block.FluidProp.Cell_pv[i][j][k].vvel+Block.FluidProp.Cell_pv[i][j][k-1].vvel);
			wvel = 0.5f*(Block.FluidProp.Cell_pv[i][j][k].wvel+Block.FluidProp.Cell_pv[i][j][k-1].wvel);
			MassFlow = MassFlow+dens*(uvel*Block.Geometry.FaceK[i][j][k].S.x+vvel*Block.Geometry.FaceK[i][j][k].S.y+wvel*Block.Geometry.FaceK[i][j][k].S.z);
		}
	}
	MassFlow = MassFlow * Block.FluidProp.dens_ref * Block.FluidProp.Vel_ref;			//有量纲化
	return MassFlow;
}

REAL COutput::TotalPressureSecCal(int K, CBlocks &Block)
{
	int i, j, k=K;
	REAL MassFlow, MassFlow_cell;
	REAL Pt_cell, Pt_cellG, Pt_med, Pt_section;
	REAL uvel, vvel, wvel, dens, temp, Ma;
	REAL Gamma = Block.FluidProp.Gamma;
	REAL Rgas  = Block.FluidProp.Rgas;
	MassFlow=0.0;
	Pt_med=0.0;
	for(i=1;i<=Block.Geometry.IM-1;i++)
	{
		for(j=1;j<=Block.Geometry.JM-1;j++)
		{
			dens = 0.5f*(Block.FluidProp.Cell_pv[i][j][k].dens + Block.FluidProp.Cell_pv[i][j][k-1].dens);
			uvel = 0.5f*(Block.FluidProp.Cell_pv[i][j][k].uvel + Block.FluidProp.Cell_pv[i][j][k-1].uvel);
			vvel = 0.5f*(Block.FluidProp.Cell_pv[i][j][k].vvel + Block.FluidProp.Cell_pv[i][j][k-1].vvel);
			wvel = 0.5f*(Block.FluidProp.Cell_pv[i][j][k].wvel + Block.FluidProp.Cell_pv[i][j][k-1].wvel);
			MassFlow_cell = dens*(uvel*Block.Geometry.FaceK[i][j][k].S.x + vvel*Block.Geometry.FaceK[i][j][k].S.y + wvel*Block.Geometry.FaceK[i][j][k].S.z);
			temp = 0.5f*(Block.FluidProp.Cell_pv[i][j][k].temp + Block.FluidProp.Cell_pv[i][j][k-1].temp);
			Ma   = SQRT(uvel*uvel+vvel*vvel+wvel*wvel) / SQRT(Gamma*Rgas*temp);
			Pt_cell  = Block.FluidProp.Cell_pv[i][j][k].press * POW(1.0f+(Gamma-1.0f)/2.0f*Ma*Ma, Gamma/(Gamma-1.0f));
			Pt_cellG = POW(Pt_cell, (Gamma-1.0f)/Gamma) * MassFlow_cell;
			Pt_med   = Pt_med + Pt_cellG;
			MassFlow = MassFlow + MassFlow_cell;
		}
	}
	//MassFlow = MassFlow * Block.FluidProp.dens_ref * Block.FluidProp.Vel_ref;									//有量纲化
	//Pt_med   = Pt_med * Block.FluidProp.dens_ref * Block.FluidProp.Vel_ref * Block.FluidProp.Vel_ref;			//有量纲化
	Pt_section = POW(Pt_med/MassFlow, Gamma/(Gamma-1.0f));
	Pt_section = Pt_section * Block.FluidProp.dens_ref * Block.FluidProp.Vel_ref * Block.FluidProp.Vel_ref;		//有量纲化
	return Pt_section;
}

REAL COutput::TotalTemperatureSecCal(int K, CBlocks &Block)
{
	int i, j, k=K;
	REAL MassFlow, MassFlow_cell;
	REAL Tt_cell, Tt_med, Tt_section;
	REAL uvel, vvel, wvel, dens, temp, Ma;
	REAL Gamma = Block.FluidProp.Gamma;
	REAL Rgas  = Block.FluidProp.Rgas;
	MassFlow=0.0;
	Tt_med=0.0;
	for(i=1;i<=Block.Geometry.IM-1;i++)
	{
		for(j=1;j<=Block.Geometry.JM-1;j++)
		{
			dens = 0.5f*(Block.FluidProp.Cell_pv[i][j][k].dens + Block.FluidProp.Cell_pv[i][j][k-1].dens);
			uvel = 0.5f*(Block.FluidProp.Cell_pv[i][j][k].uvel + Block.FluidProp.Cell_pv[i][j][k-1].uvel);
			vvel = 0.5f*(Block.FluidProp.Cell_pv[i][j][k].vvel + Block.FluidProp.Cell_pv[i][j][k-1].vvel);
			wvel = 0.5f*(Block.FluidProp.Cell_pv[i][j][k].wvel + Block.FluidProp.Cell_pv[i][j][k-1].wvel);
			MassFlow_cell = dens*(uvel*Block.Geometry.FaceK[i][j][k].S.x+vvel*Block.Geometry.FaceK[i][j][k].S.y+wvel*Block.Geometry.FaceK[i][j][k].S.z);
			temp = 0.5f*(Block.FluidProp.Cell_pv[i][j][k].temp+Block.FluidProp.Cell_pv[i][j][k-1].temp);
			Ma   = SQRT(uvel*uvel+vvel*vvel+wvel*wvel)/SQRT(Gamma*Rgas*temp);
			Tt_cell  = Block.FluidProp.Cell_pv[i][j][k].temp*POW(1+(Gamma-1.0f)/2.0f*Ma*Ma,1);
			Tt_med   = Tt_med+Tt_cell*MassFlow_cell;
			MassFlow = MassFlow+MassFlow_cell;
		}
	}
	//MassFlow = MassFlow * Block.FluidProp.dens_ref * Block.FluidProp.Vel_ref;	//有量纲化
	//Tt_med   = Tt_med * Block.FluidProp.Temp_ref;								//有量纲化
	Tt_section = Tt_med/MassFlow;
	Tt_section = Tt_section * Block.FluidProp.Temp_ref;							//有量纲化
	return Tt_section;
}

REAL COutput::EfficiencyCal(REAL Pt_in, REAL Tt_in, REAL Pt_out, REAL Tt_out)
{
	REAL Gamma = CFluidProps::Gamma;
	REAL Rgas  = CFluidProps::Rgas;
	REAL PressureRatio, L_adk, L_k, Efficiency;
	PressureRatio = Pt_out/Pt_in;
	L_adk         = Gamma/(Gamma-1.0f)*Rgas*Tt_in * (POW(PressureRatio, (Gamma-1.0f)/Gamma)-1.0f);
	L_k           = Gamma/(Gamma-1.0f)*Rgas*Tt_in * (Tt_out/Tt_in-1.0f);
	Efficiency    = L_adk/L_k;
	return Efficiency;
}

void COutput::OutputSolution(CMesh &Mesh)
{
	int i, j, k, nb;

	TPrimVars pv;
	TNode coord;
	REAL  mue, Ma_re, Ma_abs, uvel_re, vvel_re, Pt;
	CBlocks *Block = Mesh.Block;
	REAL Gamma = CFluidProps::Gamma;
	REAL Rgas  = CFluidProps::Rgas;
	REAL x, y, z, uvel, vvel, w_rot;

	cout<<"Output ..."<<endl;
	ofstream out2("RelativeSolution.plt");
	out2<<"VARIABLES = x y z dens uvel vvel wvel temp press Ma Pt  miu_t "<<endl;
	for(nb=1;nb<=Mesh.NB;nb++)
	{
		w_rot = Block[nb].FluidProp.w_rot;
		out2<<"ZONE i="<<Block[nb].Geometry.IM<<" j="<<Block[nb].Geometry.JM<<" k="<<Block[nb].Geometry.KM<<" f=point"<<endl;
		for(k=1;k<=Block[nb].Geometry.KM;k++)
		{
			for(j=1;j<=Block[nb].Geometry.JM;j++)
			{
				for(i=1;i<=Block[nb].Geometry.IM;i++)
				{
					pv=0.125f*(Block[nb].FluidProp.Cell_pv[i][j][k] + Block[nb].FluidProp.Cell_pv[i-1][j][k] + Block[nb].FluidProp.Cell_pv[i-1][j][k-1] + Block[nb].FluidProp.Cell_pv[i][j][k-1] +
							 Block[nb].FluidProp.Cell_pv[i][j-1][k] + Block[nb].FluidProp.Cell_pv[i-1][j-1][k] + Block[nb].FluidProp.Cell_pv[i-1][j-1][k-1] + Block[nb].FluidProp.Cell_pv[i][j-1][k-1]);				
					mue=0.125f*(Block[nb].FluidProp.Cell_tur[i][j][k].mue+Block[nb].FluidProp.Cell_tur[i-1][j][k].mue+Block[nb].FluidProp.Cell_tur[i-1][j][k-1].mue+Block[nb].FluidProp.Cell_tur[i][j][k-1].mue+
						Block[nb].FluidProp.Cell_tur[i][j-1][k].mue+Block[nb].FluidProp.Cell_tur[i-1][j-1][k].mue+Block[nb].FluidProp.Cell_tur[i-1][j-1][k-1].mue+Block[nb].FluidProp.Cell_tur[i][j-1][k-1].mue);
					
					coord  = Block[nb].Geometry.Node_coord[i][j][k];										
					uvel_re = pv.uvel + w_rot*coord.y;
					vvel_re = pv.vvel - w_rot*coord.x;
					Ma_re   = SQRT(uvel_re*uvel_re+vvel_re*vvel_re+pv.wvel*pv.wvel) / SQRT(Gamma*Rgas*pv.temp);
					Ma_abs = SQRT(pv.uvel*pv.uvel+pv.vvel*pv.vvel+pv.wvel*pv.wvel) / SQRT(Gamma*Rgas*pv.temp);
					Pt		 = pv.press*POW(1+(Gamma-1)/2.0f*Ma_abs*Ma_abs, Gamma/(Gamma-1));

					// 有量纲化
					x = coord.x * CGeometry::Len_ref;
					y = coord.y * CGeometry::Len_ref;
					z = coord.z * CGeometry::Len_ref;
					pv.dens = pv.dens * CFluidProps::dens_ref;
					pv.uvel = pv.uvel * CFluidProps::Vel_ref;
					pv.vvel = pv.vvel * CFluidProps::Vel_ref;
					pv.wvel = pv.wvel * CFluidProps::Vel_ref;
					pv.temp = pv.temp * CFluidProps::Temp_ref;
					pv.press= pv.press * CFluidProps::dens_ref * CFluidProps::Vel_ref * CFluidProps::Vel_ref;										
					Pt  = Pt * CFluidProps::dens_ref * CFluidProps::Vel_ref * CFluidProps::Vel_ref;
					mue = mue * CFluidProps::mue_ref * CFluidProps::Re_ref;				// Re!!!
					out2<<x<<" "<<y<<" "<<z<<" "<<pv.dens<<" "<<uvel_re<<" "<<vvel_re<<" "<<pv.wvel<<" "<<pv.temp<<" "<<pv.press<<" "<<Ma_re<<" "<<Pt<<" "<<mue<<endl;
/*					if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
					{
						ev=0.125f*(Block[nb].FluidProp.Cell_ev[i][j][k].ev + Block[nb].FluidProp.Cell_ev[i-1][j][k].ev + Block[nb].FluidProp.Cell_ev[i-1][j][k-1].ev + Block[nb].FluidProp.Cell_ev[i][j][k-1].ev + 
							Block[nb].FluidProp.Cell_ev[i][j-1][k].ev + Block[nb].FluidProp.Cell_ev[i-1][j-1][k].ev + Block[nb].FluidProp.Cell_ev[i-1][j-1][k-1].ev + Block[nb].FluidProp.Cell_ev[i][j-1][k-1].ev);
						out2<<coord.x<<" "<<coord.y<<" "<<coord.z<<" "<<pv.dens<<" "<<pv.uvel<<" "<<pv.vvel<<" "<<pv.wvel<<" "<<pv.temp<<" "<<pv.press<<" "<<Ma_rel<<" "<<Pt<<" "<<ev<<" "<<mue<<endl;
					}
					else
					{
						out2<<coord.x<<" "<<coord.y<<" "<<coord.z<<" "<<pv.dens<<" "<<pv.uvel<<" "<<pv.vvel<<" "<<pv.wvel<<" "<<pv.temp<<" "<<pv.press<<" "<<Ma_rel<<" "<<Pt<<" "<<mue<<endl;
					}
*/
				}
			}
		}	
		out2<<"ZONE i="<<Block[nb].Geometry.IM<<" j="<<Block[nb].Geometry.JM<<" k="<<Block[nb].Geometry.KM<<" f=point"<<endl;
		for(k=1;k<=Block[nb].Geometry.KM;k++)
		{
			for(j=1;j<=Block[nb].Geometry.JM;j++)
			{
				for(i=1;i<=Block[nb].Geometry.IM;i++)
				{
					pv=0.125f*(Block[nb].FluidProp.Cell_pv[i][j][k]+Block[nb].FluidProp.Cell_pv[i-1][j][k]+Block[nb].FluidProp.Cell_pv[i-1][j][k-1]+Block[nb].FluidProp.Cell_pv[i][j][k-1] +
							Block[nb].FluidProp.Cell_pv[i][j-1][k]+Block[nb].FluidProp.Cell_pv[i-1][j-1][k]+Block[nb].FluidProp.Cell_pv[i-1][j-1][k-1]+Block[nb].FluidProp.Cell_pv[i][j-1][k-1]);				
					mue=0.125f*(Block[nb].FluidProp.Cell_tur[i][j][k].mue+Block[nb].FluidProp.Cell_tur[i-1][j][k].mue+Block[nb].FluidProp.Cell_tur[i-1][j][k-1].mue+Block[nb].FluidProp.Cell_tur[i][j][k-1].mue+
						Block[nb].FluidProp.Cell_tur[i][j-1][k].mue+Block[nb].FluidProp.Cell_tur[i-1][j-1][k].mue+Block[nb].FluidProp.Cell_tur[i-1][j-1][k-1].mue+Block[nb].FluidProp.Cell_tur[i][j-1][k-1].mue);

					coord  = Block[nb].Geometry.Node_coord[i][j][k];										
					uvel_re = pv.uvel + w_rot*coord.y;
					vvel_re = pv.vvel - w_rot*coord.x;
					Ma_re   = SQRT(uvel_re*uvel_re+vvel_re*vvel_re+pv.wvel*pv.wvel) / SQRT(Gamma*Rgas*pv.temp);
					Ma_abs = SQRT(pv.uvel*pv.uvel+pv.vvel*pv.vvel+pv.wvel*pv.wvel) / SQRT(Gamma*Rgas*pv.temp);
					Pt		 = pv.press*POW(1+(Gamma-1)/2.0f*Ma_abs*Ma_abs, Gamma/(Gamma-1));
/*					if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
					{
						ev=0.125f*(Block[nb].FluidProp.Cell_ev[i][j][k].ev + Block[nb].FluidProp.Cell_ev[i-1][j][k].ev + Block[nb].FluidProp.Cell_ev[i-1][j][k-1].ev + Block[nb].FluidProp.Cell_ev[i][j][k-1].ev + 
							Block[nb].FluidProp.Cell_ev[i][j-1][k].ev + Block[nb].FluidProp.Cell_ev[i-1][j-1][k].ev + Block[nb].FluidProp.Cell_ev[i-1][j-1][k-1].ev + Block[nb].FluidProp.Cell_ev[i][j-1][k-1].ev);
						out2<<coord.x*cos(delta_th)-coord.y*sin(delta_th)<<" "<<coord.y*cos(delta_th)+coord.x*sin(delta_th)<<" "<<coord.z<<" "<<pv.dens<<" "<<pv.uvel*cos(delta_th)-pv.vvel*sin(delta_th)
							<<" "<<pv.vvel*cos(delta_th)+pv.uvel*sin(delta_th)<<" "<<pv.wvel<<" "<<pv.temp<<" "<<pv.press<<" "<<Ma_rel<<" "<<Pt<<" "<<ev<<" "<<mue<<endl;
					}
					else
					{
						out2<<coord.x*cos(delta_th)-coord.y*sin(delta_th)<<" "<<coord.y*cos(delta_th)+coord.x*sin(delta_th)<<" "<<coord.z<<" "<<pv.dens<<" "<<pv.uvel*cos(delta_th)-pv.vvel*sin(delta_th)
							<<" "<<pv.vvel*cos(delta_th)+pv.uvel*sin(delta_th)<<" "<<pv.wvel<<" "<<pv.temp<<" "<<pv.press<<" "<<Ma_rel<<" "<<Pt<<" "<<mue<<endl;
					}
*/
					// 有量纲化
					x = (coord.x*COS(Block[nb].Geometry.delta_th)-coord.y*SIN(Block[nb].Geometry.delta_th)) * CGeometry::Len_ref;
					y = (coord.y*COS(Block[nb].Geometry.delta_th)+coord.x*SIN(Block[nb].Geometry.delta_th)) * CGeometry::Len_ref;
					z = coord.z * CGeometry::Len_ref;
					pv.dens = pv.dens * CFluidProps::dens_ref;
					uvel = (uvel_re*COS(Block[nb].Geometry.delta_th)-vvel_re*SIN(Block[nb].Geometry.delta_th)) * CFluidProps::Vel_ref;
					vvel = (vvel_re*COS(Block[nb].Geometry.delta_th)+uvel_re*SIN(Block[nb].Geometry.delta_th)) * CFluidProps::Vel_ref;
					pv.wvel = pv.wvel * CFluidProps::Vel_ref;
					pv.temp = pv.temp * CFluidProps::Temp_ref;
					pv.press = pv.press * CFluidProps::dens_ref * CFluidProps::Vel_ref * CFluidProps::Vel_ref;										
					Pt  = Pt * CFluidProps::dens_ref * CFluidProps::Vel_ref * CFluidProps::Vel_ref;
					mue = mue * CFluidProps::mue_ref * CFluidProps::Re_ref;				// Re!!!

					out2<<x<<" "<<y<<" "<<z<<" "<<pv.dens<<" "<<uvel<<" "<<vvel<<" "<<pv.wvel<<" "<<pv.temp<<" "<<pv.press<<" "<<Ma_re<<" "<<Pt<<" "<<mue<<endl;
				}
			}
		}	
	}
	out2.close();

	cout<<"Output Completed!"<<endl;
}

void COutput::OutputRestartSolution(CMesh &Mesh)
{
	int i, j, k, nb;
	TPrimVars pv;
	REAL ev, ener;
	CBlocks *Block = Mesh.Block;
	ofstream out;
	out.open("RestartSolution.dat");
	for(nb=1;nb<=Mesh.NB;nb++)
	{
		for(k=1;k<=Block[nb].Geometry.KM-1;k++)
		{
			for(j=1;j<=Block[nb].Geometry.JM-1;j++)
			{
				for(i=1;i<=Block[nb].Geometry.IM-1;i++)
				{
					pv = Block[nb].FluidProp.Cell_pv[i][j][k];
					ener = Block[nb].FluidProp.Cell_cv[i][j][k].ener;
					if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
					{
						ev = Block[nb].FluidProp.Cell_ev[i][j][k].ev;
						out<<nb<<" "<<i<<" "<<j<<" "<<k<<" "<<pv.press<<" "<<pv.uvel<<" "<<pv.vvel<<" "<<pv.wvel<<" "<<pv.temp<<" "<<pv.dens<<" "<<ener<<" "<<ev<<endl;								
					}
					else
					{
						out<<nb<<" "<<i<<" "<<j<<" "<<k<<" "<<pv.press<<" "<<pv.uvel<<" "<<pv.vvel<<" "<<pv.wvel<<" "<<pv.temp<<" "<<pv.dens<<" "<<ener<<endl;								
					}
				}
			}
		}
	}
	out.close();
}

void COutput::OutputGrid(CMesh &Mesh)			//只是细网格输出
{
	int i, j, k, nb;
	ofstream out1("grid_out.plt");
	for(nb=1;nb<=Mesh.NB;nb++)
	{
		out1<<"zone i="<<Mesh.Block[nb].Geometry.IM<<" j="<<Mesh.Block[nb].Geometry.JM<<" k="<<Mesh.Block[nb].Geometry.KM<<endl;
		for(k=1;k<=Mesh.Block[nb].Geometry.KM;k++)
		{
			for(j=1;j<=Mesh.Block[nb].Geometry.JM;j++)
			{
				for(i=1;i<=Mesh.Block[nb].Geometry.IM;i++)
				{
					out1<<Mesh.Block[nb].Geometry.Node_coord[i][j][k].x<<" "<<Mesh.Block[nb].Geometry.Node_coord[i][j][k].y<<" "<<Mesh.Block[nb].Geometry.Node_coord[i][j][k].z<<endl;
				}
			}
		}
	}
	out1.close();
}

void COutput::OutputGeneralParameters(REAL &MassFlow_in, REAL &MassFlow_out, REAL &Efficiency, REAL &MassFlow_ratio, REAL &Press_ratio, CMesh &Mesh)
{
	int Block_ID, ks;	
	REAL Pt_in, Tt_in, Pt_out, Tt_out;
	Block_ID	= Mesh.BndGeom.Inlet[1].Block_ID;
	ks			= Mesh.BndGeom.Inlet[1].ks;
	MassFlow_in = MassFlowCaculation(ks, Mesh.Block[Block_ID])*Mesh.Block[Block_ID].Geometry.NBlade;	
	Pt_in		= TotalPressureSecCal(ks, Mesh.Block[Block_ID]);
	Tt_in		= TotalTemperatureSecCal(ks, Mesh.Block[Block_ID]);

	Block_ID	= Mesh.BndGeom.Outlet[1].Block_ID;
	ks			= Mesh.BndGeom.Outlet[1].ks;
	MassFlow_out= MassFlowCaculation(ks, Mesh.Block[Block_ID])*Mesh.Block[Block_ID].Geometry.NBlade;		
	Pt_out		= TotalPressureSecCal(ks, Mesh.Block[Block_ID]);
	Tt_out		= TotalTemperatureSecCal(ks, Mesh.Block[Block_ID]);
	Efficiency	= EfficiencyCal(Pt_in,Tt_in,Pt_out,Tt_out);
	MassFlow_ratio = ABS((MassFlow_in-MassFlow_out)/MassFlow_in);
	Press_ratio = Pt_out/Pt_in;
}