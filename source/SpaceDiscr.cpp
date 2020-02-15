


#include <iostream>
#include "SpaceDiscr.h"


using namespace std;

SpaceDiscrTypes  CSpaceDiscr::spaceDiscrType;
ReconstructTypes CSpaceDiscr::reconstructType;
ReconstructVars  CSpaceDiscr::reconstructVar;
string			 CSpaceDiscr::spaceDiscrName;
string			 CSpaceDiscr::reconstructTypeName;
string			 CSpaceDiscr::reconstructVarName;

REAL CSpaceDiscr::epsEntr; 
REAL CSpaceDiscr::k2;
REAL CSpaceDiscr::k4;
int  CSpaceDiscr::DiscrScheme;
int  CSpaceDiscr::ReconstructScheme;
int  CSpaceDiscr::ReconstructVar;


CSpaceDiscr::CSpaceDiscr()
{
	IM = 0;
	JM = 0;
	KM = 0;
	Mesh_ID = 0;

	Res_cv   = nullptr;
	Res_ev   = nullptr;

	Cell_Fc	 = nullptr;
	Cell_Fv	 = nullptr;
	Cell_Fd	 = nullptr;
			 
	FaceI_Fc = nullptr;
	FaceJ_Fc = nullptr;
	FaceK_Fc = nullptr;
			 
	FaceI_Fv = nullptr;
	FaceJ_Fv = nullptr;
	FaceK_Fv = nullptr;
			 
	FaceI_Fd = nullptr;
	FaceJ_Fd = nullptr;
	FaceK_Fd = nullptr;

	Cell_Fc_tur  = nullptr;
	Cell_Fv_tur  = nullptr;
	FaceI_Fc_tur = nullptr;
	FaceJ_Fc_tur = nullptr;
	FaceK_Fc_tur = nullptr;
	FaceI_Fv_tur = nullptr;
	FaceJ_Fv_tur = nullptr;
	FaceK_Fv_tur = nullptr;
	FaceI_Fv1_tur = nullptr;
	FaceJ_Fv1_tur = nullptr;
	FaceK_Fv1_tur = nullptr;
}

CSpaceDiscr::~CSpaceDiscr()
{
	int i, j;

	if(Res_cv != nullptr)
	{
		for (i=0;i<IM;i++)
		{
			for (j=0;j<JM;j++)
			{
				delete[] Res_cv[i][j];
			}
			delete[] Res_cv[i];
		}
		delete[] Res_cv;
		Res_cv = nullptr;
	}

	if(Cell_Fc != nullptr)
	{
		for (i=0;i<IM;i++)	
		{
			for (j=0;j<JM;j++)	
			{
				delete[] Cell_Fc[i][j];	
			}
			delete[] Cell_Fc[i];
		}
		delete[] Cell_Fc;
		Cell_Fc = nullptr;
	}

	if(FaceI_Fc != nullptr)
	{
		for (i=0;i<=IM;i++)	
		{
			for (j=0;j<JM;j++)	
			{
				delete[] FaceI_Fc[i][j];
			}
			delete[] FaceI_Fc[i];
		}
		delete[] FaceI_Fc;
		FaceI_Fc = nullptr;
	}
	if(FaceJ_Fc != nullptr)
	{
		for (i=0;i<IM;i++)	
		{
			for (j=0;j<=JM;j++)	
			{
				delete[] FaceJ_Fc[i][j];
			}
			delete[] FaceJ_Fc[i];
		}
		delete[] FaceJ_Fc;
		FaceJ_Fc = nullptr;
	}
	if(FaceK_Fc != nullptr)
	{
		for (i=0;i<IM;i++)	
		{
			for (j=0;j<JM;j++)	
			{
				delete[] FaceK_Fc[i][j];	
			}
			delete[] FaceK_Fc[i];	
		}
		delete[] FaceK_Fc;	
		FaceK_Fc = nullptr;
	}

	if(Cell_Fv != nullptr)
	{
		for (i=0;i<IM;i++)	
		{
			for (j=0;j<JM;j++)	
			{
				delete[] Cell_Fv[i][j];	
			}
			delete[] Cell_Fv[i];	
		}
		delete[] Cell_Fv;	
		Cell_Fv = nullptr;
	}

	if (CFluidProps::equationType==EquationTypes::NavierStokes)
	{
		if(FaceI_Fv != nullptr)
		{
			for (i=0;i<=IM;i++)	
			{
				for (j=0;j<JM;j++)	
				{
					delete[] FaceI_Fv[i][j];
				}
				delete[] FaceI_Fv[i];
			}
			delete[] FaceI_Fv;
			FaceI_Fv = nullptr;
		}
		if(FaceJ_Fv != nullptr)
		{
			for (i=0;i<IM;i++)	
			{
				for (j=0;j<=JM;j++)	
				{
					delete[] FaceJ_Fv[i][j];
				}
				delete[] FaceJ_Fv[i];
			}
			delete[] FaceJ_Fv;
			FaceJ_Fv = nullptr;
		}
		if(FaceK_Fv != nullptr)
		{
			for (i=0;i<IM;i++)	
			{
				for (j=0;j<JM;j++)	
				{
					delete[] FaceK_Fv[i][j];
				}
				delete[] FaceK_Fv[i];
			}
			delete[] FaceK_Fv;
			FaceK_Fv = nullptr;
		}

		if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
		{
			if(Res_ev != nullptr)
			{
				for (i=0;i<IM;i++)
				{
					for (j=0;j<JM;j++)
					{
						delete[] Res_ev[i][j];
					}
					delete[] Res_ev[i];
				}
				delete[] Res_ev;
				Res_ev = nullptr;
			}

			if(Cell_Fc_tur != nullptr)
			{
				for (i=0;i<IM;i++)	
				{
					for (j=0;j<JM;j++)	
					{
						delete[] Cell_Fc_tur[i][j];
					}
					delete[] Cell_Fc_tur[i];
				}
				delete[] Cell_Fc_tur;
				Cell_Fc_tur = nullptr;
			}
			if(Cell_Fv_tur != nullptr)
			{
				for (i=0;i<IM;i++)	
				{
					for (j=0;j<JM;j++)	
					{
						delete[] Cell_Fv_tur[i][j];
					}
					delete[] Cell_Fv_tur[i];
				}
				delete[] Cell_Fv_tur;
				Cell_Fv_tur = nullptr;
			}

			if( FaceI_Fc_tur != nullptr)
			{
				for (i=0;i<=IM;i++)	
				{
					for (j=0;j<JM;j++)	
					{
						delete[] FaceI_Fc_tur[i][j];
					}
					delete[] FaceI_Fc_tur[i];
				}
				delete[] FaceI_Fc_tur;
				FaceI_Fc_tur = nullptr;
			}
			if(FaceI_Fv_tur != nullptr)
			{
				for (i=0;i<=IM;i++)	
				{
					for (j=0;j<JM;j++)	
					{
						delete[] FaceI_Fv_tur[i][j];
					}
					delete[] FaceI_Fv_tur[i];
				}
				delete[] FaceI_Fv_tur;
				FaceI_Fv_tur = nullptr;
			}
			if(FaceI_Fv1_tur != nullptr)
			{
				for (i=0;i<=IM;i++)	
				{
					for (j=0;j<JM;j++)	
					{
						delete[] FaceI_Fv1_tur[i][j];
					}
					delete[] FaceI_Fv1_tur[i];
				}
				delete[] FaceI_Fv1_tur;
				FaceI_Fv1_tur = nullptr;
			}
			if(FaceJ_Fc_tur != nullptr)
			{
				for (i=0;i<IM;i++)	
				{
					for (j=0;j<=JM;j++)	
					{
						delete[] FaceJ_Fc_tur[i][j];;
					}
					delete[] FaceJ_Fc_tur[i];
				}
				delete[] FaceJ_Fc_tur;
				FaceJ_Fc_tur = nullptr;
			}
			if(FaceJ_Fv_tur != nullptr)
			{
				for (i=0;i<IM;i++)	
				{
					for (j=0;j<=JM;j++)	
					{
						delete[] FaceJ_Fv_tur[i][j];
					}
					delete[] FaceJ_Fv_tur[i];
				}
				delete[] FaceJ_Fv_tur;
				FaceJ_Fv_tur = nullptr;
			}
			if(FaceJ_Fv1_tur != nullptr)
			{
				for (i=0;i<IM;i++)	
				{
					for (j=0;j<=JM;j++)	
					{
						delete[] FaceJ_Fv1_tur[i][j];
					}
					delete[] FaceJ_Fv1_tur[i];
				}
				delete[] FaceJ_Fv1_tur;
				FaceJ_Fv1_tur = nullptr;
			}
			if(FaceK_Fc_tur != nullptr)
			{
				for (i=0;i<IM;i++)	
				{
					for (j=0;j<JM;j++)	
					{
						delete[] FaceK_Fc_tur[i][j];
					}
					delete[] FaceK_Fc_tur[i];
				}
				delete[] FaceK_Fc_tur;
				FaceK_Fc_tur = nullptr;
			}
			if(FaceK_Fv_tur != nullptr)
			{
				for (i=0;i<IM;i++)	
				{
					for (j=0;j<JM;j++)	
					{
						delete[] FaceK_Fv_tur[i][j];
					}
					delete[] FaceK_Fv_tur[i];
				}
				delete[] FaceK_Fv_tur;
				FaceK_Fv_tur = nullptr;
			}
			if(FaceK_Fv1_tur != nullptr)
			{
				for (i=0;i<IM;i++)	
				{
					for (j=0;j<JM;j++)	
					{
						delete[] FaceK_Fv1_tur[i][j];
					}
					delete[] FaceK_Fv1_tur[i];
				}
				delete[] FaceK_Fv1_tur;
				FaceK_Fv1_tur = nullptr;
			}
		}
	}

	if(Cell_Fd != nullptr)
	{
		for (i=0;i<IM;i++)	
		{
			for (j=0;j<JM;j++)	
			{
				delete[] Cell_Fd[i][j];	
			}
			delete[] Cell_Fd[i];	
		}
		delete[] Cell_Fd;	
	}

	if(spaceDiscrType==SpaceDiscrTypes::JST)
	{
		if(FaceI_Fd != nullptr)
		{
			for (i=0;i<=IM;i++)	
			{
				for (j=0;j<JM;j++)	
				{
					delete[] FaceI_Fd[i][j];
				}
				delete[] FaceI_Fd[i];
			}
			delete[] FaceI_Fd;
		}
		if(FaceJ_Fd != nullptr)
		{
			for (i=0;i<IM;i++)	
			{
				for (j=0;j<=JM;j++)	
				{
					delete[] FaceJ_Fd[i][j];
				}
				delete[] FaceJ_Fd[i];
			}
			delete[] FaceJ_Fd;
		}
		if(FaceK_Fd != nullptr)
		{
			for (i=0;i<IM;i++)	
			{
				for (j=0;j<JM;j++)	
				{
					delete[] FaceK_Fd[i][j];
				}
				delete[] FaceK_Fd[i];
			}
			delete[] FaceK_Fd;
		}
	}
}

void CSpaceDiscr::AllocateMemory(const CGeometry &Geometry)
{
	int i, j;
	IM = Geometry.IM;
	JM = Geometry.JM;
	KM = Geometry.KM;

	// cell convective flux memory
	Res_cv  = new TConsVars **[IM];
	Cell_Fc = new TConsVars **[IM];
	for (i=0;i<IM;i++)
	{
		Res_cv[i]  = new TConsVars *[JM];
		Cell_Fc[i] = new TConsVars *[JM];
	}
	for (i=0;i<IM;i++)	
	{
		for (j=0;j<JM;j++)	
		{
			Res_cv[i][j]  = new TConsVars[KM];
			Cell_Fc[i][j] = new TConsVars[KM];	
		}
	}

	// Face convective flux memory
	FaceI_Fc = new TConsVars **[IM+1];
	for (i=0;i<=IM;i++)
	{
		FaceI_Fc[i] = new TConsVars *[JM];
	}
	for (i=0;i<=IM;i++)	
	{
		for (j=0;j<JM;j++)	
		{
			FaceI_Fc[i][j] = new TConsVars[KM];
		}
	}

	FaceJ_Fc = new TConsVars **[IM];
	for (i=0;i<IM;i++)
	{
		FaceJ_Fc[i] = new TConsVars *[JM+1];
	}
	for (i=0;i<IM;i++)	
	{
		for (j=0;j<=JM;j++)	
		{
			FaceJ_Fc[i][j] = new TConsVars[KM];
		}
	}

	FaceK_Fc = new TConsVars **[IM];
	for (i=0;i<IM;i++)
	{
		FaceK_Fc[i] = new TConsVars *[JM];
	}
	for (i=0;i<IM;i++)	
	{
		for (j=0;j<JM;j++)	
		{
			FaceK_Fc[i][j] = new TConsVars[KM+1];
		}
	}

	// viscous dependent and viscous flux memory
	Cell_Fv = new TConsVars **[IM];
	for (i=0;i<IM;i++)
		Cell_Fv[i] = new TConsVars *[JM];
	for (i=0;i<IM;i++)	
	{
		for (j=0;j<JM;j++)	
		{
			Cell_Fv[i][j] = new TConsVars[KM];	
		}
	}

	if (CFluidProps::equationType==EquationTypes::NavierStokes)
	{
		FaceI_Fv = new TConsVars **[IM+1];
		for (i=0;i<=IM;i++)
			FaceI_Fv[i] = new TConsVars *[JM];
		for (i=0;i<=IM;i++)	
		{
			for (j=0;j<JM;j++)	
			{
				FaceI_Fv[i][j] = new TConsVars[KM];
			}
		}

		FaceJ_Fv = new TConsVars **[IM];
		for (i=0;i<IM;i++)
			FaceJ_Fv[i] = new TConsVars *[JM+1];
		for (i=0;i<IM;i++)	
		{
			for (j=0;j<=JM;j++)	
			{
				FaceJ_Fv[i][j] = new TConsVars[KM];
			}
		}

		FaceK_Fv = new TConsVars **[IM];
		for (i=0;i<IM;i++)
			FaceK_Fv[i] = new TConsVars *[JM];
		for (i=0;i<IM;i++)	
		{
			for (j=0;j<JM;j++)	
			{
				FaceK_Fv[i][j] = new TConsVars[KM+1];
			}
		}

		if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
		{
			Res_ev      = new  TEddyVars **[IM];
			Cell_Fc_tur = new  TEddyVars **[IM];
			Cell_Fv_tur = new  TEddyVars **[IM];
			for (i=0;i<IM;i++)
			{
				Res_ev[i]	   = new  TEddyVars *[JM];
				Cell_Fc_tur[i] = new  TEddyVars *[JM];
				Cell_Fv_tur[i] = new  TEddyVars *[JM];
			}
			for (i=0;i<IM;i++)	
			{
				for (j=0;j<JM;j++)	
				{
					Res_ev[i][j]	  = new  TEddyVars[KM];
					Cell_Fc_tur[i][j] = new  TEddyVars[KM];
					Cell_Fv_tur[i][j] = new  TEddyVars[KM];
				}
			}

			FaceI_Fc_tur = new TEddyVars **[IM+1];
			FaceI_Fv_tur = new TEddyVars **[IM+1];
			FaceI_Fv1_tur = new TEddyVars **[IM+1];
			for (i=0;i<=IM;i++)
			{
				FaceI_Fc_tur[i] = new TEddyVars *[JM];
				FaceI_Fv_tur[i] = new TEddyVars *[JM];
				FaceI_Fv1_tur[i] = new TEddyVars *[JM];
			}
			for (i=0;i<=IM;i++)	
			{
				for (j=0;j<JM;j++)	
				{
					FaceI_Fc_tur[i][j] = new TEddyVars[KM];
					FaceI_Fv_tur[i][j] = new TEddyVars[KM];
					FaceI_Fv1_tur[i][j] = new TEddyVars[KM];
				}
			}

			FaceJ_Fc_tur = new TEddyVars **[IM];
			FaceJ_Fv_tur = new TEddyVars **[IM];
			FaceJ_Fv1_tur = new TEddyVars **[IM];
			for (i=0;i<IM;i++)
			{
				FaceJ_Fc_tur[i] = new TEddyVars *[JM+1];
				FaceJ_Fv_tur[i] = new TEddyVars *[JM+1];
				FaceJ_Fv1_tur[i] = new TEddyVars *[JM+1];
			}
			for (i=0;i<IM;i++)	
			{
				for (j=0;j<=JM;j++)	
				{
					FaceJ_Fc_tur[i][j] = new TEddyVars[KM];
					FaceJ_Fv_tur[i][j] = new TEddyVars[KM];
					FaceJ_Fv1_tur[i][j] = new TEddyVars[KM];
				}
			}

			FaceK_Fc_tur = new TEddyVars **[IM];
			FaceK_Fv_tur = new TEddyVars **[IM];
			FaceK_Fv1_tur = new TEddyVars **[IM];
			for (i=0;i<IM;i++)
			{
				FaceK_Fc_tur[i] = new TEddyVars *[JM];
				FaceK_Fv_tur[i] = new TEddyVars *[JM];
				FaceK_Fv1_tur[i] = new TEddyVars *[JM];
			}
			for (i=0;i<IM;i++)	
			{
				for (j=0;j<JM;j++)	
				{
					FaceK_Fc_tur[i][j] = new TEddyVars[KM+1];
					FaceK_Fv_tur[i][j] = new TEddyVars[KM+1];
					FaceK_Fv1_tur[i][j] = new TEddyVars[KM+1];
				}
			}
		}
	}

	// used in JST scheme's dissipation	
	Cell_Fd = new TConsVars **[IM];
	for (i=0;i<IM;i++)
		Cell_Fd[i] = new TConsVars *[JM];
	for (i=0;i<IM;i++)	
	{
		for (j=0;j<JM;j++)	
		{
			Cell_Fd[i][j] = new TConsVars[KM];	
		}
	}
	if(spaceDiscrType==SpaceDiscrTypes::JST)
	{
		FaceI_Fd = new TConsVars **[IM+1];
		for (i=0;i<=IM;i++)
			FaceI_Fd[i] = new TConsVars *[JM];
		for (i=0;i<=IM;i++)	
		{
			for (j=0;j<JM;j++)	
			{
				FaceI_Fd[i][j] = new TConsVars[KM];
			}
		}

		FaceJ_Fd = new TConsVars **[IM];
		for (i=0;i<IM;i++)
			FaceJ_Fd[i] = new TConsVars *[JM+1];
		for (i=0;i<IM;i++)	
		{
			for (j=0;j<=JM;j++)	
			{
				FaceJ_Fd[i][j] = new TConsVars[KM];
			}
		}

		FaceK_Fd = new TConsVars **[IM];
		for (i=0;i<IM;i++)
			FaceK_Fd[i] = new TConsVars *[JM];
		for (i=0;i<IM;i++)	
		{
			for (j=0;j<JM;j++)	
			{
				FaceK_Fd[i][j] = new TConsVars[KM+1];
			}
		}
	}
}