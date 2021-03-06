

#include "TimeDiscr.h"

using namespace std;

void CTimeDiscr::TimeStep(const CGeometry &Geometry, CFluidProps &FluidProp)
{
	int i, j, k;

	REAL Pr_lam = FluidProp.Pr_lam;
	REAL Pr_tur = FluidProp.Pr_tur;
	REAL Gamma  = FluidProp.Gamma;
	int  Flag_low_dt = FluidProp.Flag_low_dt;
	
	TNode SI, SJ, SK;

	REAL eq_SI, eq_SJ, eq_SK, mue_lam, mue_tur, lambda_c_all, lambda_v_all, dens;

	REAL dt_factor;
	// use the limit flow variables
	if(FluidProp.Flag_low_dt == 0)
		dt_factor = 1.0f;
	else
		dt_factor = 1.0f/(1.0f+FluidProp.Flag_low_dt);

#pragma omp parallel for num_threads(NT) default(shared) private(i, j, k, SI, SJ, SK, eq_SI, eq_SJ, eq_SK, mue_lam, mue_tur, lambda_c_all, lambda_v_all, dens)
	for(i=1;i<=IM-1;i++)
	{
		for(j=1;j<=JM-1;j++)
		{
			for(k=1;k<=KM-1;k++)
			{	
				SI  = 0.5*(Geometry.FaceI[i][j][k].S + Geometry.FaceI[i+1][j][k].S);
				SJ  = 0.5*(Geometry.FaceJ[i][j][k].S + Geometry.FaceJ[i][j+1][k].S);
				SK  = 0.5*(Geometry.FaceK[i][j][k].S + Geometry.FaceK[i][j][k+1].S);
	
				eq_SI = SI.x*SI.x + SI.y*SI.y + SI.z*SI.z;
				eq_SJ = SJ.x*SJ.x + SJ.y*SJ.y + SJ.z*SJ.z;
				eq_SK = SK.x*SK.x + SK.y*SK.y + SK.z*SK.z;

				if(FluidProp.equationType==EquationTypes::NavierStokes)
				{
					mue_lam = FluidProp.SutherLand(FluidProp.Cell_pv[i][j][k].temp);
					mue_tur = FluidProp.Cell_tur[i][j][k].mue;

					dens   = FluidProp.Cell_pv[i][j][k].dens;
					FluidProp.lambda_vI[i][j][k] = MAX(4.0f/3.0f/dens, Gamma/dens) * (mue_lam/Pr_lam+mue_tur/Pr_tur)*eq_SI/Geometry.Cell[i][j][k].Vol;
					FluidProp.lambda_vJ[i][j][k] = MAX(4.0f/3.0f/dens, Gamma/dens) * (mue_lam/Pr_lam+mue_tur/Pr_tur)*eq_SJ/Geometry.Cell[i][j][k].Vol;
					FluidProp.lambda_vK[i][j][k] = MAX(4.0f/3.0f/dens, Gamma/dens) * (mue_lam/Pr_lam+mue_tur/Pr_tur)*eq_SK/Geometry.Cell[i][j][k].Vol;
					lambda_c_all = FluidProp.lambda_cI[i][j][k] + FluidProp.lambda_cJ[i][j][k] + FluidProp.lambda_cK[i][j][k];
					lambda_v_all = FluidProp.lambda_vI[i][j][k] + FluidProp.lambda_vJ[i][j][k] + FluidProp.lambda_vK[i][j][k];
				}
				else  // non viscous flow
				{
					lambda_c_all = FluidProp.lambda_cI[i][j][k] + FluidProp.lambda_cJ[i][j][k] + FluidProp.lambda_cK[i][j][k];
					lambda_v_all = 0.0;
				}
				dt[i][j][k]  = dt_factor*CFL*Geometry.Cell[i][j][k].Vol/(lambda_c_all + C_timestep*lambda_v_all);
			}
		}
	}

	if(timestepType==TimestepTypes::Global)
	{
		REAL dt_min = 1E10;
		for(i=1;i<=IM-1;i++)
		{
			for(j=1;j<=JM-1;j++)
			{
				for(k=1;k<=KM-1;k++)
				{
					dt_min = MIN(dt_min, dt[i][j][k]);
				}
			}
		}
		for(i=1;i<=IM-1;i++)
		{
			for(j=1;j<=JM-1;j++)
			{
				for(k=1;k<=KM-1;k++)
				{
					dt[i][j][k] = dt_min;
				}
			}
		}
	}
	
	if(FluidProp.Flag_low_dt > 0) 
	{
		FluidProp.Flag_low_dt = FluidProp.Flag_low_dt - 1;		//降低时间步长的标志，今后Flag_low_dt步内，还将降低时间步长
	}
}