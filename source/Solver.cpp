


#include "Solver.h"

using namespace std;

bool   CSolver::UseFlowMassImposed = false;

CSolver::CSolver()
{
	reUseSolution     = false;
	reUseWallDistance = false;
	UseMultigrid = false;
	MaxIter		= 0;
	OutStep		= 0;
	iter		= 0;
	Coarse_iter = 0;
	ConvTol     = 0.0;
	NG			= 4;
	CFL_res		= 1.0E-3f;

	MassFlow_in		= 0.0;
	MassFlow_out	= 0.0;
	Efficiency		= 0.0;
	MassFlow_ratio	= 0.0;
	Press_ratio		= 0.0;

	MassFlow_Imposed = 0.0;
	relax_press = 0.0;
	resc_mass_imp = 0.0;
	iter_mass_imp = 0; 

	Res_ev = 0.0;
}

CSolver::~CSolver()
{
	if(Mesh != nullptr)
	{
		delete[] Mesh;
		Mesh = nullptr;
	}
}


void CSolver::Solver()
{
	// CFL amplify
	REAL CFL_amplify=1.2f, CFL_MAX = 10.0f;
	if(iter>1)
	{
		if(Res.dens<=CFL_res)
		{
			CTimeDiscr::CFL = CTimeDiscr::CFL * CFL_amplify;
			CFL_res = CFL_res / 10.0f;
		}
		if(CTimeDiscr::CFL > CFL_MAX)
			CTimeDiscr::CFL = CFL_MAX; 
	}

	if(NG==1)
	{
		Mesh[1].Time_Marching();
		Mesh[1].Res_Calculate(Res, Res_ev);
	}
	else
	{
		if(reUseSolution==false)
		{
			if(iter < Coarse_iter)
			{
				Mesh[2].Time_Marching();
				Mesh[2].Res_Calculate(Res, Res_ev);
			}
			else if(iter == Coarse_iter)
			{
				Mesh[2].Time_Marching();
				Mesh[2].Res_Calculate(Res, Res_ev);
				Mesh[2].CornerTreatment();
				Interpolate(Mesh[1], Mesh[2]);
				Mesh[1].PrimitiveVars();
				if(CFluidProps::limiter_vars==1)
				{				
					Mesh[1].Limit_Vars();
				}
			}
			else
			{
				Mesh[2].UseForceFunction  = true;

				Mesh[1].Time_Marching();
				Mesh[1].Res_Calculate(Res, Res_ev);

				Mesh[1].ResFlux();
				Prolong(Mesh[1], Mesh[2]);
				Mesh[2].PrimitiveVars();							
				if(CFluidProps::limiter_vars==1)
				{
					Mesh[2].Limit_Vars();
				}

				Mesh[2].Time_Marching();
				ComputeDeltConsVars(Mesh[2]);
				Mesh[2].CornerTreatment();
				Interpolate(Mesh[1], Mesh[2]);
				Mesh[1].PrimitiveVars();				
				if(CFluidProps::limiter_vars==1)
				{					
					Mesh[1].Limit_Vars();
				}
			}
		}
		else
		{
			Mesh[2].UseForceFunction  = true;

			Mesh[1].Time_Marching();
			Mesh[1].Res_Calculate(Res, Res_ev);

			Mesh[1].ResFlux();
			Prolong(Mesh[1], Mesh[2]);
			Mesh[2].PrimitiveVars();				
			if(CFluidProps::limiter_vars==1)
			{				
				Mesh[2].Limit_Vars();
			}

			Mesh[2].Time_Marching();
			ComputeDeltConsVars(Mesh[2]);
			Mesh[2].CornerTreatment();
			Interpolate(Mesh[1], Mesh[2]);
			Mesh[1].PrimitiveVars();
			if(CFluidProps::limiter_vars==1)
			{				
				Mesh[1].Limit_Vars();
			}
		}
	}
}

void CSolver::MassFlowImposed()
{
	CBndConds::Pout = CBndConds::Pout + relax_press*CBndConds::Pout*(MassFlow_out-MassFlow_Imposed)/MassFlow_Imposed;
}


