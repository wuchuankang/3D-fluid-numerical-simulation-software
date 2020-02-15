

#include "Mesh.h"
#include <iostream>
#include <fstream>
#include "UserInput.h"

using namespace std;

void CMesh::AllocateMemory()
{
	int nb;
	for(nb=1;nb<=NB;nb++)
	{
		Block[nb].Geometry.AllocateMemory();
		Block[nb].FluidProp.AllocateMemory(Block[nb].Geometry);		
		Block[nb].SpaceDiscr.AllocateMemory(Block[nb].Geometry);
		Block[nb].TimeDiscr.AllocateMemory(Block[nb].Geometry);
		Block[nb].AllocateMemory();

		Block[nb].SpaceDiscr.Mesh_ID = this->Mesh_ID;
	}
}

void CMesh::Geometry()
{
	int nb;
	TotalVol = 0.0;
	for(nb=1;nb<=NB;nb++)
	{
		Block[nb].Geometry.FaceVectorVolume();
		TotalVol = TotalVol + Block[nb].Geometry.CellVol;
	}
	for(nb=1;nb<=NB;nb++)
		Block[nb].Geometry.JacobiConvertCoeff();
}

void CMesh::Initial()
{
	int nb;
	if(InitialType == InitialTypes::Uniform)
	{
		for(nb=1;nb<=NB;nb++)
			Block[nb].Initial_Uniform();
	}
	else
		Initial_Interpolate(Block);
}

void CMesh::Spectral_curl_grad()
{
	int nb;
	for(nb=1;nb<=NB;nb++)
		Block[nb].FluidProp.SpectralRadii(Block[nb].Geometry);	

	if(CFluidProps::equationType==EquationTypes::NavierStokes)
	{
		for(nb=1;nb<=NB;nb++)
		{
			Block[nb].FluidProp.Curl(Block[nb].Geometry);
			Block[nb].FluidProp.Gradient(Block[nb].Geometry);
		}
	}
}

void CMesh::TimeStep()
{
	int nb;
	for(nb=1;nb<=NB;nb++)
		Block[nb].TimeDiscr.TimeStep(Block[nb].Geometry, Block[nb].FluidProp);
}

void CMesh::ConvectiveFlux()
{
	int i, j, k, nb;
		
	for(nb=1;nb<=NB;nb++)
	{
		Block[nb].SpaceDiscr.ConvectiveFlux(Block[nb].Geometry, Block[nb].FluidProp);
	}

	WallFlux();

	//convective = convective + source 
	REAL uvel, vvel, coord_x, coord_y, w_rot;
	TConsVars  Cell_Qc;
	for(nb=1;nb<=NB;nb++)
	{
		w_rot = Block[nb].FluidProp.w_rot;
#pragma omp parallel for num_threads(NT) default(shared) private(i,j,k, uvel, vvel, coord_x, coord_y, Cell_Qc)
		for(i=1;i<=Block[nb].Geometry.IM-1;i++)
		{
			for(j=1;j<=Block[nb].Geometry.JM-1;j++)		
			{
				for(k=1;k<=Block[nb].Geometry.KM-1;k++)
				{
					uvel    = Block[nb].FluidProp.Cell_pv[i][j][k].uvel;
					vvel    = Block[nb].FluidProp.Cell_pv[i][j][k].vvel;
					coord_x = Block[nb].Geometry.Cell[i][j][k].coord.x;
					coord_y = Block[nb].Geometry.Cell[i][j][k].coord.y;
					Cell_Qc.dens = 0.0;
					Cell_Qc.xmom = Block[nb].FluidProp.Cell_pv[i][j][k].dens * (w_rot*vvel) * Block[nb].Geometry.Cell[i][j][k].Vol;
					Cell_Qc.ymom = -1.0f*Block[nb].FluidProp.Cell_pv[i][j][k].dens * (w_rot*uvel) * Block[nb].Geometry.Cell[i][j][k].Vol;
					Cell_Qc.zmom = 0.0;
					Cell_Qc.ener = 0.0;

					Block[nb].SpaceDiscr.Cell_Fc[i][j][k] = Block[nb].SpaceDiscr.FaceI_Fc[i+1][j][k] - Block[nb].SpaceDiscr.FaceI_Fc[i][j][k] + Block[nb].SpaceDiscr.FaceJ_Fc[i][j+1][k] - 
														    Block[nb].SpaceDiscr.FaceJ_Fc[i][j][k] + Block[nb].SpaceDiscr.FaceK_Fc[i][j][k+1] - Block[nb].SpaceDiscr.FaceK_Fc[i][j][k] - Cell_Qc;																							
					if(CFluidProps::turbulenceType==TurbulenceTypes::SA)
					{
						Block[nb].SpaceDiscr.Cell_Fc_tur[i][j][k].ev = Block[nb].SpaceDiscr.FaceI_Fc_tur[i+1][j][k].ev - Block[nb].SpaceDiscr.FaceI_Fc_tur[i][j][k].ev + 
																	   Block[nb].SpaceDiscr.FaceJ_Fc_tur[i][j+1][k].ev - Block[nb].SpaceDiscr.FaceJ_Fc_tur[i][j][k].ev + 
																	   Block[nb].SpaceDiscr.FaceK_Fc_tur[i][j][k+1].ev - Block[nb].SpaceDiscr.FaceK_Fc_tur[i][j][k].ev;
					}
				}
			}
		}
	}
}

void CMesh::WallFlux()
{
	int i, Block_ID;
	for(i=1;i<=BndGeom.nWall;i++)
	{
		Block_ID = BndGeom.Wall[i].Block_ID;
		Block[Block_ID].SpaceDiscr.WallFlux(BndGeom.Wall[i], Block[Block_ID].Geometry, Block[Block_ID].FluidProp);
	}
}

void CMesh::ViscousFlux()
{
	int i, j, k, nb;
	if(CFluidProps::equationType==EquationTypes::NavierStokes)
	{
		for(nb=1;nb<=NB;nb++)
			Block[nb].SpaceDiscr.ViscousFlux(Block[nb].Geometry, Block[nb].FluidProp);
	}
	else
	{
		for(nb=1;nb<=NB;nb++)
		{
#pragma omp parallel for num_threads(NT) private(i,j,k)
			for(i=1;i<=Block[nb].Geometry.IM-1;i++)
			{
				for(j=1;j<=Block[nb].Geometry.JM-1;j++)
				{
					for(k=1;k<=Block[nb].Geometry.KM-1;k++)
					{
						Block[nb].SpaceDiscr.Cell_Fv[i][j][k].dens = 0.0;
						Block[nb].SpaceDiscr.Cell_Fv[i][j][k].xmom = 0.0;
						Block[nb].SpaceDiscr.Cell_Fv[i][j][k].ymom = 0.0;
						Block[nb].SpaceDiscr.Cell_Fv[i][j][k].zmom = 0.0;
						Block[nb].SpaceDiscr.Cell_Fv[i][j][k].ener = 0.0;
					}
				}
			}
		}
	}
}

void CMesh::Dissipation()
{
	int i, j, k, nb;
	if(CSpaceDiscr::spaceDiscrType==SpaceDiscrTypes::JST)
	{
		for(nb=1;nb<=NB;nb++)
			Block[nb].SpaceDiscr.JST_Dissip(Block[nb].FluidProp);
	}
	else
	{
		for(nb=1;nb<=NB;nb++)
		{
#pragma omp parallel for num_threads(NT) private(i,j,k)
			for(i=1;i<=Block[nb].Geometry.IM-1;i++)
			{
				for(j=1;j<=Block[nb].Geometry.JM-1;j++)
				{
					for(k=1;k<=Block[nb].Geometry.KM-1;k++)
					{
						Block[nb].SpaceDiscr.Cell_Fd[i][j][k].dens = 0.0;
						Block[nb].SpaceDiscr.Cell_Fd[i][j][k].xmom = 0.0;
						Block[nb].SpaceDiscr.Cell_Fd[i][j][k].ymom = 0.0;
						Block[nb].SpaceDiscr.Cell_Fd[i][j][k].zmom = 0.0;
						Block[nb].SpaceDiscr.Cell_Fd[i][j][k].ener = 0.0;
					}
				}
			}
		}
	}
}

void CMesh::PrimitiveVars()
{
	int nb;
	for(nb=1;nb<=NB;nb++)
		Block[nb].FluidProp.PrimitiveVars();
}

void CMesh::Turbulence_Model()
{
	int nb;
	if(CFluidProps::turbulenceType==TurbulenceTypes::None)
	{

	}
	if(CFluidProps::turbulenceType==TurbulenceTypes::BL)
	{
		for(nb=1;nb<=NB;nb++)
			Block[nb].BL_Model(BndGeom);
	}
	if(CFluidProps::turbulenceType==TurbulenceTypes::SA)
	{
		for(nb=1;nb<=NB;nb++)
			Block[nb].SpaceDiscr.SA_Model(Block[nb].Geometry, Block[nb].FluidProp);
	}
	if(CFluidProps::turbulenceType==TurbulenceTypes::SST)
	{

	}
}

void CMesh::ResFlux()
{
	int i, j, k, nb;

	BoundaryCondition();
	Spectral_curl_grad();
	if(Mesh_ID == 1)						// turbulence model is only used at fine mesh 
		Turbulence_Model();					// turbulence model also need treat the wallflux when it is SA model

	TimeStep();								// because timestep need the turbulence viscous, so it is reasonable after the Turbulence_Model()
	ConvectiveFlux();
	ViscousFlux();
	Dissipation();

	for(nb=1;nb<=NB;nb++)
	{
#pragma omp parallel for num_threads(NT) private(i,j,k)
		for (i=1;i<=Block[nb].Geometry.IM-1;i++)
		{
			for (j=1;j<=Block[nb].Geometry.JM-1;j++)
			{
				for(k=1;k<=Block[nb].Geometry.KM-1;k++)
				{
					Block[nb].SpaceDiscr.Res_cv[i][j][k] = Block[nb].SpaceDiscr.Cell_Fc[i][j][k] - Block[nb].SpaceDiscr.Cell_Fv[i][j][k] - Block[nb].SpaceDiscr.Cell_Fd[i][j][k];	
					if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
					{
						Block[nb].SpaceDiscr.Res_ev[i][j][k].ev = Block[nb].SpaceDiscr.Cell_Fc_tur[i][j][k].ev - Block[nb].SpaceDiscr.Cell_Fv_tur[i][j][k].ev;
					}
				}
			}
		}
	}
}

void CMesh::Limit_Vars()
{
	int nb;
#pragma omp parallel for num_threads(NT)
	for(nb=1;nb<=NB;nb++)
	{
		Block[nb].FluidProp.Limit_Vars();
		if(CFluidProps::turbulenceType==TurbulenceTypes::SA)
			Block[nb].FluidProp.Limit_SA();
	}
}

void CMesh::ComputeForceFunction()
{
	int i, j, k, nb;
	for(nb=1;nb<=NB;nb++)
	{
#pragma omp parallel for num_threads(NT) private(i, j, k)
		for(i=1;i<=Block[nb].Geometry.IM-1;i++)
		{
			for(j=1;j<=Block[nb].Geometry.JM-1;j++)
			{
				for(k=1;k<=Block[nb].Geometry.KM-1;k++)
				{
					Qf_2h[nb][i][j][k] = Res_2h[nb][i][j][k] -  Block[nb].SpaceDiscr.Res_cv[i][j][k];
				}
			}
		}
	}
}

void CMesh::AddForceFunction()
{
	int i, j, k, nb;
	// 细网格上不用强迫函数
	if(Mesh_ID == 2)
	{
		for(nb=1;nb<=NB;nb++)
		{
#pragma omp parallel for num_threads(NT) private(i,j,k)
			for (i=1;i<=Block[nb].Geometry.IM-1;i++)
			{
				for (j=1;j<=Block[nb].Geometry.JM-1;j++)
				{
					for(k=1;k<=Block[nb].Geometry.KM-1;k++)
					{
						Block[nb].SpaceDiscr.Res_cv[i][j][k] = Block[nb].SpaceDiscr.Res_cv[i][j][k] + Qf_2h[nb][i][j][k];	
					}
				}
			}
		}		
	}
}

void CMesh::IRS()
{
	int nb;
	for(nb=1;nb<=NB;nb++)
		Block[nb].TimeDiscr.UpwindIRS(Block[nb].Geometry, Block[nb].FluidProp, Block[nb].SpaceDiscr);
}

void CMesh::Res_Calculate(TConsVars &Res, REAL &Res_ev)
{
	int i, j, k, nb;		
	TConsVars r_med;
	REAL	  r_med_ev, res_dens, res_xmom, res_ymom, res_zmom, res_ener, res_ev;

	res_dens = 0.0;
	res_xmom = 0.0;
	res_ymom = 0.0;
	res_zmom = 0.0;
	res_ener = 0.0;
	res_ev	 = 0.0;

	for(nb=1;nb<=NB;nb++)	//residual
	{
#pragma omp parallel for num_threads(NT) private(i,j,k, r_med, r_med_ev) reduction(+:res_dens, res_xmom, res_ymom, res_zmom, res_ener, res_ev)
		for(i=1;i<=Block[nb].Geometry.IM-1;i++)	
		{
			for(j=1;j<=Block[nb].Geometry.JM-1;j++)
			{
				for(k=1;k<=Block[nb].Geometry.KM-1;k++)
				{	
					r_med.dens = ABS((Block[nb].TimeDiscr.Cell_cv_old[i][j][k].dens-Block[nb].FluidProp.Cell_cv[i][j][k].dens)/Block[nb].TimeDiscr.Cell_cv_old[i][j][k].dens);
					r_med.xmom = ABS((Block[nb].TimeDiscr.Cell_cv_old[i][j][k].xmom-Block[nb].FluidProp.Cell_cv[i][j][k].xmom)/Block[nb].TimeDiscr.Cell_cv_old[i][j][k].xmom);
					r_med.ymom = ABS((Block[nb].TimeDiscr.Cell_cv_old[i][j][k].ymom-Block[nb].FluidProp.Cell_cv[i][j][k].ymom)/Block[nb].TimeDiscr.Cell_cv_old[i][j][k].ymom);
					r_med.zmom = ABS((Block[nb].TimeDiscr.Cell_cv_old[i][j][k].zmom-Block[nb].FluidProp.Cell_cv[i][j][k].zmom)/Block[nb].TimeDiscr.Cell_cv_old[i][j][k].zmom);
					r_med.ener = ABS((Block[nb].TimeDiscr.Cell_cv_old[i][j][k].ener-Block[nb].FluidProp.Cell_cv[i][j][k].ener)/Block[nb].TimeDiscr.Cell_cv_old[i][j][k].ener);

					r_med = r_med*Block[nb].Geometry.Cell[i][j][k].Vol;
					res_dens = res_dens + r_med.dens;
					res_xmom = res_xmom + r_med.xmom;
					res_ymom = res_ymom + r_med.ymom;
					res_zmom = res_zmom + r_med.zmom;
					res_ener = res_ener + r_med.ener;
					if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
					{
						r_med_ev = ABS((Block[nb].TimeDiscr.Cell_ev_old[i][j][k].ev-Block[nb].FluidProp.Cell_ev[i][j][k].ev)/Block[nb].TimeDiscr.Cell_ev_old[i][j][k].ev);	
						r_med_ev = r_med_ev*Block[nb].Geometry.Cell[i][j][k].Vol;
						res_ev    = res_ev + r_med_ev;
					}
				}
			}
		}
	}

	Res.dens = res_dens/TotalVol;
	Res.xmom = res_xmom/TotalVol;
	Res.ymom = res_ymom/TotalVol;
	Res.zmom = res_zmom/TotalVol;
	Res.ener = res_ener/TotalVol;
	if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
	{
		Res_ev = res_ev/TotalVol;
	}
}

void CMesh::CornerTreatment()
{
	int i, j, k, nb;
	int IM, JM, KM;
	//12条边
	for(nb=1;nb<=NB;nb++)
	{
		IM = Block[nb].Geometry.IM;
		JM = Block[nb].Geometry.JM;
		KM = Block[nb].Geometry.KM;
		for(i=1;i<=IM-1;i++)
		{
			Block[nb].FluidProp.Cell_cv[i][0][0]   = 0.5f*(Block[nb].FluidProp.Cell_cv[i][1][0]		+ Block[nb].FluidProp.Cell_cv[i][0][1]);
			Block[nb].FluidProp.Cell_cv[i][JM][0]  = 0.5f*(Block[nb].FluidProp.Cell_cv[i][JM-1][0]  + Block[nb].FluidProp.Cell_cv[i][JM][1]);
			Block[nb].FluidProp.Cell_cv[i][0][KM]  = 0.5f*(Block[nb].FluidProp.Cell_cv[i][0][KM-1]  + Block[nb].FluidProp.Cell_cv[i][1][KM]);
			Block[nb].FluidProp.Cell_cv[i][JM][KM] = 0.5f*(Block[nb].FluidProp.Cell_cv[i][JM-1][KM] + Block[nb].FluidProp.Cell_cv[i][JM][KM-1]);
					
			delt_cv[nb][i][0][0]   = 0.5f*(delt_cv[nb][i][1][0]		+ delt_cv[nb][i][0][1]);
			delt_cv[nb][i][JM][0]  = 0.5f*(delt_cv[nb][i][JM-1][0]  + delt_cv[nb][i][JM][1]);
			delt_cv[nb][i][0][KM]  = 0.5f*(delt_cv[nb][i][0][KM-1]  + delt_cv[nb][i][1][KM]);
			delt_cv[nb][i][JM][KM] = 0.5f*(delt_cv[nb][i][JM-1][KM] + delt_cv[nb][i][JM][KM-1]);
		}
		for(j=1;j<=JM-1;j++)
		{
			Block[nb].FluidProp.Cell_cv[0][j][0]   = 0.5f*(Block[nb].FluidProp.Cell_cv[1][j][0]     + Block[nb].FluidProp.Cell_cv[0][j][1]);
			Block[nb].FluidProp.Cell_cv[0][j][KM]  = 0.5f*(Block[nb].FluidProp.Cell_cv[1][j][KM]    + Block[nb].FluidProp.Cell_cv[0][j][KM-1]);
			Block[nb].FluidProp.Cell_cv[IM][j][KM] = 0.5f*(Block[nb].FluidProp.Cell_cv[IM-1][j][KM] + Block[nb].FluidProp.Cell_cv[IM][j][KM-1]);
			Block[nb].FluidProp.Cell_cv[IM][j][0]  = 0.5f*(Block[nb].FluidProp.Cell_cv[IM-1][j][0]  + Block[nb].FluidProp.Cell_cv[IM][j][1]);
					
			delt_cv[nb][0][j][0]   = 0.5f*(delt_cv[nb][1][j][0]     + delt_cv[nb][0][j][1]);
			delt_cv[nb][0][j][KM]  = 0.5f*(delt_cv[nb][1][j][KM]    + delt_cv[nb][0][j][KM-1]);
			delt_cv[nb][IM][j][KM] = 0.5f*(delt_cv[nb][IM-1][j][KM] + delt_cv[nb][IM][j][KM-1]);
			delt_cv[nb][IM][j][0]  = 0.5f*(delt_cv[nb][IM-1][j][0]  + delt_cv[nb][IM][j][1]);
		}
		for(k=1;k<=KM-1;k++)
		{
			Block[nb].FluidProp.Cell_cv[0][0][k]   = 0.5f*(Block[nb].FluidProp.Cell_cv[1][0][k]     + Block[nb].FluidProp.Cell_cv[0][1][k]);
			Block[nb].FluidProp.Cell_cv[IM][0][k]  = 0.5f*(Block[nb].FluidProp.Cell_cv[IM][1][k]    + Block[nb].FluidProp.Cell_cv[IM-1][0][k]);
			Block[nb].FluidProp.Cell_cv[0][JM][k]  = 0.5f*(Block[nb].FluidProp.Cell_cv[1][JM][k]    + Block[nb].FluidProp.Cell_cv[0][JM-1][k]);
			Block[nb].FluidProp.Cell_cv[IM][JM][k] = 0.5f*(Block[nb].FluidProp.Cell_cv[IM-1][JM][k] + Block[nb].FluidProp.Cell_cv[IM][JM-1][k]);
				
			delt_cv[nb][0][0][k]   = 0.5f*(delt_cv[nb][1][0][k]     + delt_cv[nb][0][1][k]);
			delt_cv[nb][IM][0][k]  = 0.5f*(delt_cv[nb][IM][1][k]    + delt_cv[nb][IM-1][0][k]);
			delt_cv[nb][0][JM][k]  = 0.5f*(delt_cv[nb][1][JM][k]    + delt_cv[nb][0][JM-1][k]);
			delt_cv[nb][IM][JM][k] = 0.5f*(delt_cv[nb][IM-1][JM][k] + delt_cv[nb][IM][JM-1][k]);
		}

		Block[nb].FluidProp.Cell_cv[0][0][0]    = (2.0f*(Block[nb].FluidProp.Cell_cv[1][0][0]      + Block[nb].FluidProp.Cell_cv[0][1][0]      + Block[nb].FluidProp.Cell_cv[0][0][1])      - (Block[nb].FluidProp.Cell_cv[1][1][0]        + Block[nb].FluidProp.Cell_cv[1][0][1]        + Block[nb].FluidProp.Cell_cv[0][1][1]))/3.0f;
		Block[nb].FluidProp.Cell_cv[IM][0][0]   = (2.0f*(Block[nb].FluidProp.Cell_cv[IM-1][0][0]   + Block[nb].FluidProp.Cell_cv[IM][1][0]     + Block[nb].FluidProp.Cell_cv[IM][0][1])     - (Block[nb].FluidProp.Cell_cv[IM-1][1][0]     + Block[nb].FluidProp.Cell_cv[IM-1][0][1]     + Block[nb].FluidProp.Cell_cv[IM][1][1]))/3.0f;
		Block[nb].FluidProp.Cell_cv[0][JM][0]   = (2.0f*(Block[nb].FluidProp.Cell_cv[1][JM][0]     + Block[nb].FluidProp.Cell_cv[0][JM-1][0]   + Block[nb].FluidProp.Cell_cv[0][JM][1])     - (Block[nb].FluidProp.Cell_cv[1][JM-1][0]     + Block[nb].FluidProp.Cell_cv[1][JM][1]       + Block[nb].FluidProp.Cell_cv[0][JM-1][1]))/3.0f;
		Block[nb].FluidProp.Cell_cv[IM][JM][0]  = (2.0f*(Block[nb].FluidProp.Cell_cv[IM-1][JM][0]  + Block[nb].FluidProp.Cell_cv[IM][JM-1][0]  + Block[nb].FluidProp.Cell_cv[IM][JM][1])    - (Block[nb].FluidProp.Cell_cv[IM-1][JM-1][0]  + Block[nb].FluidProp.Cell_cv[IM-1][JM][1]    + Block[nb].FluidProp.Cell_cv[IM][JM-1][1]))/3.0f;
		Block[nb].FluidProp.Cell_cv[IM][0][KM]  = (2.0f*(Block[nb].FluidProp.Cell_cv[IM-1][0][KM]  + Block[nb].FluidProp.Cell_cv[IM][1][KM]	   + Block[nb].FluidProp.Cell_cv[IM][0][KM-1])  - (Block[nb].FluidProp.Cell_cv[IM-1][1][KM]    + Block[nb].FluidProp.Cell_cv[IM-1][0][KM-1]  + Block[nb].FluidProp.Cell_cv[IM][1][KM-1]))/3.0f;
		Block[nb].FluidProp.Cell_cv[0][0][KM]   = (2.0f*(Block[nb].FluidProp.Cell_cv[1][0][KM]     + Block[nb].FluidProp.Cell_cv[0][1][KM]	   + Block[nb].FluidProp.Cell_cv[0][0][KM-1])   - (Block[nb].FluidProp.Cell_cv[1][1][KM]       + Block[nb].FluidProp.Cell_cv[1][0][KM-1]     + Block[nb].FluidProp.Cell_cv[0][1][KM-1]))/3.0f;
		Block[nb].FluidProp.Cell_cv[0][JM][KM]  = (2.0f*(Block[nb].FluidProp.Cell_cv[1][JM][KM]    + Block[nb].FluidProp.Cell_cv[0][JM-1][KM]  + Block[nb].FluidProp.Cell_cv[0][JM][KM-1])  - (Block[nb].FluidProp.Cell_cv[1][JM-1][KM]    + Block[nb].FluidProp.Cell_cv[1][JM][KM-1]    + Block[nb].FluidProp.Cell_cv[0][JM-1][KM-1]))/3.0f;
		Block[nb].FluidProp.Cell_cv[IM][JM][KM] = (2.0f*(Block[nb].FluidProp.Cell_cv[IM-1][JM][KM] + Block[nb].FluidProp.Cell_cv[IM][JM-1][KM] + Block[nb].FluidProp.Cell_cv[IM][JM][KM-1]) - (Block[nb].FluidProp.Cell_cv[IM-1][JM-1][KM] + Block[nb].FluidProp.Cell_cv[IM-1][JM][KM-1] + Block[nb].FluidProp.Cell_cv[IM][JM-1][KM-1]))/3.0f;
			
		delt_cv[nb][0][0][0]    = (2.0f*(delt_cv[nb][1][0][0]      + delt_cv[nb][0][1][0]      + delt_cv[nb][0][0][1])      - (delt_cv[nb][1][1][0]        + delt_cv[nb][1][0][1]        + delt_cv[nb][0][1][1]))/3.0f;
		delt_cv[nb][IM][0][0]   = (2.0f*(delt_cv[nb][IM-1][0][0]   + delt_cv[nb][IM][1][0]     + delt_cv[nb][IM][0][1])     - (delt_cv[nb][IM-1][1][0]     + delt_cv[nb][IM-1][0][1]     + delt_cv[nb][IM][1][1]))/3.0f;
		delt_cv[nb][0][JM][0]   = (2.0f*(delt_cv[nb][1][JM][0]     + delt_cv[nb][0][JM-1][0]   + delt_cv[nb][0][JM][1])     - (delt_cv[nb][1][JM-1][0]     + delt_cv[nb][1][JM][1]       + delt_cv[nb][0][JM-1][1]))/3.0f;
		delt_cv[nb][IM][JM][0]  = (2.0f*(delt_cv[nb][IM-1][JM][0]  + delt_cv[nb][IM][JM-1][0]  + delt_cv[nb][IM][JM][1])    - (delt_cv[nb][IM-1][JM-1][0]  + delt_cv[nb][IM-1][JM][1]    + delt_cv[nb][IM][JM-1][1]))/3.0f;
		delt_cv[nb][IM][0][KM]  = (2.0f*(delt_cv[nb][IM-1][0][KM]  + delt_cv[nb][IM][1][KM]	   + delt_cv[nb][IM][0][KM-1])  - (delt_cv[nb][IM-1][1][KM]    + delt_cv[nb][IM-1][0][KM-1]  + delt_cv[nb][IM][1][KM-1]))/3.0f;
		delt_cv[nb][0][0][KM]   = (2.0f*(delt_cv[nb][1][0][KM]     + delt_cv[nb][0][1][KM]	   + delt_cv[nb][0][0][KM-1])   - (delt_cv[nb][1][1][KM]       + delt_cv[nb][1][0][KM-1]     + delt_cv[nb][0][1][KM-1]))/3.0f;
		delt_cv[nb][0][JM][KM]  = (2.0f*(delt_cv[nb][1][JM][KM]    + delt_cv[nb][0][JM-1][KM]  + delt_cv[nb][0][JM][KM-1])  - (delt_cv[nb][1][JM-1][KM]    + delt_cv[nb][1][JM][KM-1]    + delt_cv[nb][0][JM-1][KM-1]))/3.0f;
		delt_cv[nb][IM][JM][KM] = (2.0f*(delt_cv[nb][IM-1][JM][KM] + delt_cv[nb][IM][JM-1][KM] + delt_cv[nb][IM][JM][KM-1]) - (delt_cv[nb][IM-1][JM-1][KM] + delt_cv[nb][IM-1][JM][KM-1] + delt_cv[nb][IM][JM-1][KM-1]))/3.0f;
	}
}



