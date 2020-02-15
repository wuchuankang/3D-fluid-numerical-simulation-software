


#include "FluidProps.h"
#include "SpaceDiscr.h"

bool			CFluidProps::limiter_vars;
EquationTypes	CFluidProps::equationType;
TurbulenceTypes CFluidProps::turbulenceType;
string			CFluidProps::equationName;
string			CFluidProps::turbulenceName;

int		CFluidProps::turbulenceScheme = 0;				
REAL	CFluidProps::Gamma    = 0.0,
    	CFluidProps::Pr_lam   = 0.0,
    	CFluidProps::Pr_tur   = 0.0,
    	CFluidProps::Rgas	  = 0.0,
		CFluidProps::Ma0	  = 0.0,
    	CFluidProps::Vel_ref  = 0.0,
		CFluidProps::Temp_ref = 0.0,
		CFluidProps::Press_ref= 0.0,
    	CFluidProps::dens_ref = 0.0,
		CFluidProps::mue_ref  = 0.0,
    	CFluidProps::Ma_ref   = 0.0,
		CFluidProps::Re_ref   = 0.0;

CFluidProps::CFluidProps()
{
	IM = 0;
	JM = 0;
	KM = 0;
	NDummy = 0;
	
	Flag_low_dt = 0;
	LdensMin  = 0.0;
	LdensMax  = 0.0;
	LpressMax = 0.0;
	LpressMin = 0.0;
	LVvelMax  = 0.0;

	rpm	  = 0.0;
	w_rot = 0.0;

	Cell_cv		= nullptr;
	Cell_pv		= nullptr;
	Cell_ev	    = nullptr;
	lambda_cI   = nullptr;
	lambda_cJ   = nullptr;
	lambda_cK   = nullptr;	
	lambda_vI   = nullptr;
	lambda_vJ   = nullptr;
	lambda_vK   = nullptr;

	Omega		= nullptr;

	DummyIS_cv	= nullptr;
	DummyIP_cv	= nullptr;
	DummyJS_cv	= nullptr;
	DummyJP_cv	= nullptr;
	DummyKS_cv	= nullptr;
	DummyKP_cv	= nullptr;
	DummyIS_pv	= nullptr;
	DummyIP_pv	= nullptr;
	DummyJS_pv	= nullptr;
	DummyJP_pv	= nullptr;
	DummyKS_pv	= nullptr;
	DummyKP_pv	= nullptr;

	Cell_tur    = nullptr;	// turbulence heat conduction ,dynamic viscosity 

	gradx_I		= nullptr;
	grady_I		= nullptr;
	gradz_I		= nullptr;
	gradx_J		= nullptr;
	grady_J		= nullptr;
	gradz_J		= nullptr;
	gradx_K		= nullptr;
	grady_K		= nullptr;
	gradz_K		= nullptr;

	gradx_I_ev  = nullptr;
	grady_I_ev  = nullptr;
	gradz_I_ev  = nullptr;
	gradx_J_ev  = nullptr;
	grady_J_ev  = nullptr;
	gradz_J_ev  = nullptr;
	gradx_K_ev  = nullptr;
	grady_K_ev  = nullptr;
	gradz_K_ev  = nullptr;
}

CFluidProps::~CFluidProps()
{
	int i, j;

	if(Cell_cv != nullptr)
	{
		for (i=0;i<=IM;i++)	
		{
			for (j=0;j<=JM;j++)	
			{
				delete[] Cell_cv[i][j];
			}
			delete[] Cell_cv[i];
		}
		delete[] Cell_cv;
		Cell_cv = nullptr;
	}
	if(Cell_pv != nullptr)
	{
		for (i=0;i<=IM;i++)	
		{
			for (j=0;j<=JM;j++)	
			{
				delete[] Cell_pv[i][j];
			}
			delete[] Cell_pv[i];
		}
		delete[] Cell_pv;
		Cell_pv = nullptr;
	}
	if(Cell_tur != nullptr)
	{
		for (i=0;i<=IM;i++)	
		{
			for (j=0;j<=JM;j++)	
			{
				delete[] Cell_tur[i][j];
			}
			delete[] Cell_tur[i];
		}
		delete[] Cell_tur;
		Cell_tur = nullptr;
	}
	//------------------------------------
	if(DummyIS_cv != nullptr)
	{
		for (i=0;i<=NDummy;i++)	
		{
			for (j=0;j<JM;j++)	
			{
				delete[] DummyIS_cv[i][j];		
			}
			delete[] DummyIS_cv[i];
		}
		delete[] DummyIS_cv;
		DummyIS_cv = nullptr;
	}
	if(DummyIP_cv != nullptr)
	{
		for (i=0;i<=NDummy;i++)	
		{
			for (j=0;j<JM;j++)	
			{
				delete[] DummyIP_cv[i][j];		
			}
			delete[] DummyIP_cv[i];
		}
		delete[] DummyIP_cv;
		DummyIP_cv = nullptr;
	}
	if(DummyIS_pv != nullptr)
	{
		for (i=0;i<=NDummy;i++)	
		{
			for (j=0;j<JM;j++)	
			{
				delete[] DummyIS_pv[i][j];			
			}
			delete[] DummyIS_pv[i];
		}
		delete[] DummyIS_pv;
		DummyIS_pv = nullptr;
	}
	if(DummyIP_pv != nullptr)
	{
		for (i=0;i<=NDummy;i++)	
		{
			for (j=0;j<JM;j++)	
			{
				delete[] DummyIP_pv[i][j];			
			}
			delete[] DummyIP_pv[i];
		}
		delete[] DummyIP_pv;
		DummyIP_pv = nullptr;
	}

	if(DummyJS_cv != nullptr)
	{
		for (i=0;i<IM;i++)	
		{
			for (j=0;j<=NDummy;j++)	
			{
				delete[] DummyJS_cv[i][j];		
			}
			delete[] DummyJS_cv[i];
		}
		delete[] DummyJS_cv;
		DummyJS_cv = nullptr;
	}
	if(DummyJP_cv != nullptr)
	{
		for (i=0;i<IM;i++)	
		{
			for (j=0;j<=NDummy;j++)	
			{				
				delete[] DummyJP_cv[i][j];		
			}
			delete[] DummyJP_cv[i];
		}
		delete[] DummyJP_cv;
		DummyJP_cv = nullptr;
	}
	if(DummyJS_pv != nullptr)
	{
		for (i=0;i<IM;i++)	
		{
			for (j=0;j<=NDummy;j++)	
			{
				delete[] DummyJS_pv[i][j];			
			}
			delete[] DummyJS_pv[i];
		}
		delete[] DummyJS_pv;
		DummyJS_pv = nullptr;
	}
	if(DummyJP_pv != nullptr)
	{
		for (i=0;i<IM;i++)	
		{
			for (j=0;j<=NDummy;j++)	
			{
				delete[] DummyJP_pv[i][j];			
			}
			delete[] DummyJP_pv[i];
		}
		delete[] DummyJP_pv;
		DummyJP_pv = nullptr;
	}

	if(DummyKS_cv != nullptr)
	{
		for (i=0;i<IM;i++)	
		{
			for (j=0;j<JM;j++)	
			{
				delete[] DummyKS_cv[i][j];		
			}
			delete[] DummyKS_cv[i];
		}
		delete[] DummyKS_cv;
		DummyKS_cv = nullptr;
	}
	if(DummyKP_cv != nullptr)
	{
		for (i=0;i<IM;i++)	
		{
			for (j=0;j<JM;j++)	
			{
				delete[] DummyKP_cv[i][j];		
			}
			delete[] DummyKP_cv[i];
		}
		delete[] DummyKP_cv;
		DummyKP_cv = nullptr;
	}
	if(DummyKS_pv != nullptr)
	{
		for (i=0;i<IM;i++)	
		{
			for (j=0;j<JM;j++)	
			{
				delete[] DummyKS_pv[i][j];			
			}
			delete[] DummyKS_pv[i];
		}
		delete[] DummyKS_pv;
		DummyKS_pv = nullptr;
	}
	if(DummyKP_pv != nullptr)
	{
		for (i=0;i<IM;i++)	
		{
			for (j=0;j<JM;j++)	
			{	
				delete[] DummyKP_pv[i][j];			
			}
			delete[] DummyKP_pv[i];
		}
		delete[] DummyKP_pv;
		DummyKP_pv = nullptr;
	}
	//-----------------------------------------------
	if(lambda_cI != nullptr)
	{
		for (i=0;i<IM;i++)	
		{
			for (j=0;j<JM;j++)	
			{
				delete[] lambda_cI[i][j];
			}
			delete[] lambda_cI[i];
		}
		delete[] lambda_cI;
		lambda_cI = nullptr;
	}
	if(lambda_cJ != nullptr)
	{
		for (i=0;i<IM;i++)	
		{
			for (j=0;j<JM;j++)	
			{
				delete[] lambda_cJ[i][j];
			}
			delete[] lambda_cJ[i];
		}
		delete[] lambda_cJ;
		lambda_cJ = nullptr;
	}
	if(lambda_cK != nullptr)
	{
		for (i=0;i<IM;i++)	
		{
			for (j=0;j<JM;j++)	
			{
				delete[] lambda_cK[i][j];
			}
			delete[] lambda_cK[i];
		}
		delete[] lambda_cK;
		lambda_cK = nullptr;
	}

	//---------------------------------
	if (equationType==EquationTypes::NavierStokes)
	{
		if(Omega != nullptr)
		{
			for (i=0;i<=IM;i++)	
			{
				for (j=0;j<=JM;j++)	
				{
					delete[] Omega[i][j];
				}
				delete[] Omega[i];
			}
			delete[] Omega;
			Omega = nullptr;
		}

		if(lambda_vI != nullptr)
		{
			for (i=0;i<IM;i++)	
			{
				for (j=0;j<JM;j++)	
				{
					delete[] lambda_vI[i][j];
				}
				delete[] lambda_vI[i];
			}
			delete[] lambda_vI;
			lambda_vI = nullptr;
		}
		if(lambda_vJ != nullptr)
		{
			for (i=0;i<IM;i++)	
			{
				for (j=0;j<JM;j++)	
				{
					delete[] lambda_vJ[i][j];
				}
				delete[] lambda_vJ[i];
			}
			delete[] lambda_vJ;
			lambda_vJ = nullptr;
		}
		if(lambda_vK != nullptr)
		{
			for (i=0;i<IM;i++)	
			{
				for (j=0;j<JM;j++)	
				{
					delete[] lambda_vK[i][j];
				}
				delete[] lambda_vK[i];
			}
			delete[] lambda_vK;
			lambda_vK = nullptr;
		}

		if(gradx_I != nullptr)
		{
			for (i=0;i<=IM;i++)	
			{
				for (j=0;j<JM;j++)	
				{
					delete[] gradx_I[i][j];
				}
				delete[] gradx_I[i];
			}
			delete[] gradx_I;
			gradx_I = nullptr;
		}
		if(grady_I != nullptr)
		{
			for (i=0;i<=IM;i++)	
			{
				for (j=0;j<JM;j++)	
				{
					delete[] grady_I[i][j];
				}
				delete[] grady_I[i];
			}
			delete[] grady_I;
			grady_I = nullptr;
		}
		if(gradz_I != nullptr)
		{
			for (i=0;i<=IM;i++)	
			{
				for (j=0;j<JM;j++)	
				{
					delete[] gradz_I[i][j];
				}
				delete[] gradz_I[i];
			}
			delete[] gradz_I;
			gradz_I = nullptr;
		}

		if(gradx_J != nullptr)
		{
			for (i=0;i<IM;i++)	
			{
				for (j=0;j<=JM;j++)	
				{
					delete[] gradx_J[i][j];
				}
				delete[] gradx_J[i];
			}
			delete[] gradx_J;
			gradx_J = nullptr;
		}
		if(grady_J != nullptr)
		{
			for (i=0;i<IM;i++)	
			{
				for (j=0;j<=JM;j++)	
				{
					delete[] grady_J[i][j];
				}
				delete[] grady_J[i];
			}
			delete[] grady_J;
			grady_J = nullptr;
		}
		if(gradz_J != nullptr)
		{
			for (i=0;i<IM;i++)	
			{
				for (j=0;j<=JM;j++)	
				{
					delete[] gradz_J[i][j];
				}
				delete[] gradz_J[i];
			}
			delete[] gradz_J;
			gradz_J = nullptr;
		}

		if(gradx_K != nullptr)
		{
			for (i=0;i<IM;i++)	
			{
				for (j=0;j<JM;j++)	
				{
					delete[] gradx_K[i][j];
				}
				delete[] gradx_K[i];
			}
			delete[] gradx_K;
			gradx_K = nullptr;
		}

		if(grady_K != nullptr)
		{
			for (i=0;i<IM;i++)	
			{
				for (j=0;j<JM;j++)	
				{
					delete[] grady_K[i][j];
				}
				delete[] grady_K[i];
			}
			delete[] grady_K;
			grady_K = nullptr;
		}
		if(gradz_K != nullptr)
		{
			for (i=0;i<IM;i++)	
			{
				for (j=0;j<JM;j++)	
				{
					delete[] gradz_K[i][j];
				}
				delete[] gradz_K[i];
			}
			delete[] gradz_K;
			gradz_K = nullptr;
		}

		if (turbulenceType==TurbulenceTypes::SA)
		{
			if(Cell_ev != nullptr)
			{
				for (i=0;i<=IM;i++)	
				{
					for (j=0;j<=JM;j++)	
					{
						delete[] Cell_ev[i][j];
					}
					delete[] Cell_ev[i];
				}
				delete[] Cell_ev;
				Cell_ev = nullptr;
			}

			if(gradx_I_ev != nullptr)
			{
				for (i=0;i<=IM;i++)	
				{
					for (j=0;j<JM;j++)	
					{
						delete[] gradx_I_ev[i][j];
					}
					delete[] gradx_I_ev[i];
				}
				delete[] gradx_I_ev;
				gradx_I_ev = nullptr;
			}
			if(grady_I_ev != nullptr)
			{
				for (i=0;i<=IM;i++)	
				{
					for (j=0;j<JM;j++)	
					{
						delete[] grady_I_ev[i][j];
					}
					delete[] grady_I_ev[i];
				}
				delete[] grady_I_ev;
				grady_I_ev = nullptr;
			}
			if(gradz_I_ev != nullptr)
			{
				for (i=0;i<=IM;i++)	
				{
					for (j=0;j<JM;j++)	
					{
						delete[] gradz_I_ev[i][j];
					}
					delete[] gradz_I_ev[i];
				}
				delete[] gradz_I_ev;
				gradz_I_ev = nullptr;
			}
			if(gradx_J_ev != nullptr)
			{
				for (i=0;i<IM;i++)	
				{
					for (j=0;j<=JM;j++)	
					{
						delete[] gradx_J_ev[i][j];
					}
					delete[] gradx_J_ev[i];
				}
				delete[] gradx_J_ev;
				gradx_J_ev = nullptr;
			}
			if(grady_J_ev != nullptr)
			{
				for (i=0;i<IM;i++)	
				{
					for (j=0;j<=JM;j++)	
					{
						delete[] grady_J_ev[i][j];
					}
					delete[] grady_J_ev[i];
				}
				delete[] grady_J_ev;
				grady_J_ev = nullptr;
			}
			if(gradz_J_ev != nullptr)
			{
				for (i=0;i<IM;i++)	
				{
					for (j=0;j<=JM;j++)	
					{
						delete[] gradz_J_ev[i][j];
					}
					delete[] gradz_J_ev[i];
				}
				delete[] gradz_J_ev;
				gradz_J_ev = nullptr;
			}
			if(gradx_K_ev != nullptr)
			{
				for (i=0;i<IM;i++)	
				{
					for (j=0;j<JM;j++)	
					{
						delete[] gradx_K_ev[i][j];
					}
					delete[] gradx_K_ev[i];
				}
				delete[] gradx_K_ev;
				gradx_K_ev = nullptr;
			}
			if(grady_K_ev != nullptr)
			{
				for (i=0;i<IM;i++)	
				{
					for (j=0;j<JM;j++)	
					{
						delete[] grady_K_ev[i][j];
					}
					delete[] grady_K_ev[i];
				}
				delete[] grady_K_ev;
				grady_K_ev = nullptr;
			}
			if(gradz_K_ev != nullptr)
			{
				for (i=0;i<IM;i++)	
				{
					for (j=0;j<JM;j++)	
					{
						delete[] gradz_K_ev[i][j];
					}
					delete[] gradz_K_ev[i];
				}
				delete[] gradz_K_ev;
				gradz_K_ev = nullptr;
			}
		}
	}
}

void CFluidProps::AllocateMemory(const CGeometry &Geometry)
{
	int i, j;
	IM     = Geometry.IM;
	JM     = Geometry.JM;
	KM     = Geometry.KM;
	NDummy = Geometry.NDummy;
	// cell memory
	Cell_cv  = new TConsVars **[IM+1];
	Cell_pv  = new TPrimVars **[IM+1];
	Cell_tur = new TVisVars **[IM+1];
	for (i=0;i<=IM;i++)
	{
		Cell_cv[i]  = new TConsVars *[JM+1];
		Cell_pv[i]  = new TPrimVars *[JM+1];
		Cell_tur[i] = new TVisVars *[JM+1];
	}
	for (i=0;i<=IM;i++)	
	{
		for (j=0;j<=JM;j++)	
		{
			Cell_cv[i][j]  = new TConsVars[KM+1];
			Cell_pv[i][j]  = new TPrimVars[KM+1];
			Cell_tur[i][j] = new TVisVars[KM+1];
		}
	}

	// dummy memory
	// I direction
	DummyIS_cv = new TConsVars **[NDummy+1];
	DummyIP_cv = new TConsVars **[NDummy+1];
	DummyIS_pv = new TPrimVars  **[NDummy+1];
	DummyIP_pv = new TPrimVars  **[NDummy+1];
	for (i=0;i<=NDummy;i++)
	{
		DummyIS_cv[i] = new TConsVars *[JM];
		DummyIP_cv[i] = new TConsVars *[JM];
		DummyIS_pv[i] = new TPrimVars  *[JM];
		DummyIP_pv[i] = new TPrimVars  *[JM];
	}
	for (i=0;i<=NDummy;i++)	
	{
		for (j=0;j<JM;j++)	
		{
			DummyIS_cv[i][j] = new TConsVars[KM];
			DummyIP_cv[i][j] = new TConsVars[KM];
			DummyIS_pv[i][j] = new TPrimVars[KM];
			DummyIP_pv[i][j] = new TPrimVars[KM];			
		}
	}

	// J direction
	DummyJS_cv = new TConsVars **[IM];
	DummyJP_cv = new TConsVars **[IM];
	DummyJS_pv = new TPrimVars  **[IM];
	DummyJP_pv = new TPrimVars  **[IM];
	for (i=0;i<IM;i++)
	{
		DummyJS_cv[i] = new TConsVars *[NDummy+1];
		DummyJP_cv[i] = new TConsVars *[NDummy+1];
		DummyJS_pv[i] = new TPrimVars *[NDummy+1];
		DummyJP_pv[i] = new TPrimVars *[NDummy+1];
	}
	for (i=0;i<IM;i++)	
	{
		for (j=0;j<=NDummy;j++)	
		{
			DummyJS_cv[i][j] = new TConsVars[KM];
			DummyJP_cv[i][j] = new TConsVars[KM];
			DummyJS_pv[i][j] = new TPrimVars[KM];
			DummyJP_pv[i][j] = new TPrimVars[KM];			
		}
	} 

	// K direction
	DummyKS_cv = new TConsVars **[IM];
	DummyKP_cv = new TConsVars **[IM];
	DummyKS_pv = new TPrimVars  **[IM];
	DummyKP_pv = new TPrimVars  **[IM];
	for (i=0;i<IM;i++)
	{
		DummyKS_cv[i] = new TConsVars *[JM];
		DummyKP_cv[i] = new TConsVars *[JM];
		DummyKS_pv[i] = new TPrimVars  *[JM];
		DummyKP_pv[i] = new TPrimVars  *[JM];
	}
	for (i=0;i<IM;i++)	
	{
		for (j=0;j<JM;j++)	
		{
			DummyKS_cv[i][j] = new TConsVars[NDummy+1];
			DummyKP_cv[i][j] = new TConsVars[NDummy+1];
			DummyKS_pv[i][j] = new TPrimVars[NDummy+1];	
			DummyKP_pv[i][j] = new TPrimVars[NDummy+1];			
		}
	} 

	lambda_cI = new REAL **[IM];
	lambda_cJ = new REAL **[IM];
	lambda_cK = new REAL **[IM];
	for (i=0;i<IM;i++)
	{
		lambda_cI[i] = new REAL *[JM];
		lambda_cJ[i] = new REAL *[JM];
		lambda_cK[i] = new REAL *[JM];
	}
	for (i=0;i<IM;i++)	
	{
		for (j=0;j<JM;j++)	
		{
			lambda_cI[i][j] = new REAL[KM];
			lambda_cJ[i][j] = new REAL[KM];
			lambda_cK[i][j] = new REAL[KM];
		}
	}

	// viscous dependent and viscous flux memory
	if (equationType==EquationTypes::NavierStokes)
	{
		lambda_vI = new REAL **[IM];
		lambda_vJ = new REAL **[IM];
		lambda_vK = new REAL **[IM];
		for (i=0;i<IM;i++)
		{
			lambda_vI[i] = new REAL *[JM];
			lambda_vJ[i] = new REAL *[JM];
			lambda_vK[i] = new REAL *[JM];
		}
		for (i=0;i<IM;i++)	
		{
			for (j=0;j<JM;j++)	
			{
				lambda_vI[i][j] = new REAL[KM];
				lambda_vJ[i][j] = new REAL[KM];
				lambda_vK[i][j] = new REAL[KM];
			}
		}
		// gradient
		gradx_I = new TPrimVars **[IM+1];
		grady_I = new TPrimVars **[IM+1];
		gradz_I = new TPrimVars **[IM+1];
		for (i=0;i<=IM;i++)
		{
			gradx_I[i] = new TPrimVars *[JM];
			grady_I[i] = new TPrimVars *[JM];
			gradz_I[i] = new TPrimVars *[JM];
		}
		for (i=0;i<=IM;i++)	
		{
			for (j=0;j<JM;j++)	
			{
				gradx_I[i][j] = new TPrimVars[KM];
				grady_I[i][j] = new TPrimVars[KM];
				gradz_I[i][j] = new TPrimVars[KM];
			}
		}

		gradx_J = new TPrimVars **[IM];
		grady_J = new TPrimVars **[IM];
		gradz_J = new TPrimVars **[IM];
		for (i=0;i<IM;i++)
		{
			gradx_J[i] = new TPrimVars *[JM+1];
			grady_J[i] = new TPrimVars *[JM+1];
			gradz_J[i] = new TPrimVars *[JM+1];
		}
		for (i=0;i<IM;i++)	
		{
			for (j=0;j<=JM;j++)	
			{
				gradx_J[i][j] = new TPrimVars[KM];
				grady_J[i][j] = new TPrimVars[KM];
				gradz_J[i][j] = new TPrimVars[KM];
			}
		}

		gradx_K = new TPrimVars **[IM];
		grady_K = new TPrimVars **[IM];
		gradz_K = new TPrimVars **[IM];
		for (i=0;i<IM;i++)
		{
			gradx_K[i] = new TPrimVars *[JM];
			grady_K[i] = new TPrimVars *[JM];
			gradz_K[i] = new TPrimVars *[JM];
		}
		for (i=0;i<IM;i++)	
		{
			for (j=0;j<JM;j++)	
			{
				gradx_K[i][j] = new TPrimVars[KM+1];
				grady_K[i][j] = new TPrimVars[KM+1];
				gradz_K[i][j] = new TPrimVars[KM+1];
			}
		}

		// turbulence viscous and curl at cell
		Omega	 = new REAL  **[IM+1];
		for (i=0;i<=IM;i++)
			Omega[i]    = new REAL	*[JM+1];
		for (i=0;i<=IM;i++)	
		{
			for (j=0;j<=JM;j++)	
			{
				Omega[i][j]	   = new REAL[KM+1];
			}
		}

		if (turbulenceType==TurbulenceTypes::SA)
		{
			// cell viscous memory
			Cell_ev = new TEddyVars **[IM+1];
			for (i=0;i<=IM;i++)
				Cell_ev[i] = new TEddyVars *[JM+1];
			for (i=0;i<=IM;i++)	
			{
				for (j=0;j<=JM;j++)	
				{
					Cell_ev[i][j] = new TEddyVars[KM+1];
				}
			}

			gradx_I_ev = new REAL **[IM+1];
			grady_I_ev = new REAL **[IM+1];
			gradz_I_ev = new REAL **[IM+1];
			for (i=0;i<=IM;i++)
			{
				gradx_I_ev[i] = new REAL *[JM];
				grady_I_ev[i] = new REAL *[JM];
				gradz_I_ev[i] = new REAL *[JM];
			}
			for (i=0;i<=IM;i++)	
			{
				for (j=0;j<JM;j++)	
				{
					gradx_I_ev[i][j] = new REAL[KM];
					grady_I_ev[i][j] = new REAL[KM];
					gradz_I_ev[i][j] = new REAL[KM];
				}
			}

			gradx_J_ev = new REAL **[IM];
			grady_J_ev = new REAL **[IM];
			gradz_J_ev = new REAL **[IM];
			for (i=0;i<IM;i++)
			{
				gradx_J_ev[i] = new REAL *[JM+1];
				grady_J_ev[i] = new REAL *[JM+1];
				gradz_J_ev[i] = new REAL *[JM+1];
			}
			for (i=0;i<IM;i++)	
			{
				for (j=0;j<=JM;j++)	
				{
					gradx_J_ev[i][j] = new REAL[KM];
					grady_J_ev[i][j] = new REAL[KM];
					gradz_J_ev[i][j] = new REAL[KM];
				}
			}

			gradx_K_ev = new REAL **[IM];
			grady_K_ev = new REAL **[IM];
			gradz_K_ev = new REAL **[IM];
			for (i=0;i<IM;i++)
			{
				gradx_K_ev[i] = new REAL *[JM];
				grady_K_ev[i] = new REAL *[JM];
				gradz_K_ev[i] = new REAL *[JM];
			}
			for (i=0;i<IM;i++)	
			{
				for (j=0;j<JM;j++)	
				{
					gradx_K_ev[i][j] = new REAL[KM+1];
					grady_K_ev[i][j] = new REAL[KM+1];
					gradz_K_ev[i][j] = new REAL[KM+1];
				}
			}
		}
	}	
}

