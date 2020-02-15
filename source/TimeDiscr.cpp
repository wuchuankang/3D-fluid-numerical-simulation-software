


#include "TimeDiscr.h"

 TimestepTypes  CTimeDiscr::timestepType;
 TimeDiscrTypes CTimeDiscr::timeDiscrType;
 string			CTimeDiscr::timestepName;
 string			CTimeDiscr::timeDiscrName;
 string			CTimeDiscr::timeDiscrGeneralName;

 int			CTimeDiscr::TimeDiscrScheme;
 REAL			CTimeDiscr::CFL;
 REAL			CTimeDiscr::relax_factor;	
 REAL			CTimeDiscr::C_timestep;	
 bool			CTimeDiscr::IRS = false;

CTimeDiscr::CTimeDiscr()
{
	IM = 0;
	JM = 0;
	KM = 0;

	Cell_cv_old = nullptr;
	Cell_Fv_old = nullptr;
	Cell_Fd_old = nullptr;
	Cell_ev_old = nullptr;
	dt			= nullptr;
	delt_cv		= nullptr;
	delt_ev		= nullptr;

	epsI   = nullptr;
	epsJ   = nullptr;
	epsK   = nullptr;
	a_IRS = nullptr;
	b_IRS = nullptr;
	c_IRS = nullptr;
	R_IRS = nullptr;
	Rs_IRS = nullptr;
}

CTimeDiscr::~CTimeDiscr()
{
	int i, j;

	if(dt != nullptr)
	{
		for (i=0;i<IM;i++)
		{
			for (j=0;j<JM;j++)
			{
				delete[] dt[i][j];
			}
			delete[] dt[i];
		}
		delete[] dt;
		dt = nullptr;
	}

	if(Cell_cv_old != nullptr)
	{
		for (i=0;i<IM;i++)
		{
			for (j=0;j<JM;j++)
			{
				delete[] Cell_cv_old[i][j];
			}
			delete[] Cell_cv_old[i];
		}
		delete[] Cell_cv_old;
		Cell_cv_old = nullptr;
	}

	if(string(timeDiscrGeneralName)==string("Runge_Kutta"))
	{
		if(Cell_Fv_old != nullptr)
		{
			for (i=0;i<IM;i++)
			{
				for (j=0;j<JM;j++)
				{
					delete[] Cell_Fv_old[i][j];
				}
				delete[] Cell_Fv_old[i];
			}
			delete[] Cell_Fv_old;
			Cell_Fv_old = nullptr;
		}
		if(Cell_Fd_old != nullptr)
		{
			for (i=0;i<IM;i++)
			{
				for (j=0;j<JM;j++)
				{
					delete[] Cell_Fd_old[i][j];
				}
				delete[] Cell_Fd_old[i];
			}
			delete[] Cell_Fd_old;
			Cell_Fd_old = nullptr;
		}
	}
	if(IRS == 1)
	{		
		if(epsI != nullptr)
		{
			delete[] epsI;
			epsI = nullptr;
		}
		if(epsJ != nullptr)
		{
			delete[] epsJ;
			epsJ = nullptr;
		}
		if(epsK != nullptr)
		{
			delete[] epsK;
			epsK = nullptr;
		}
		if(a_IRS != nullptr)
		{
			delete[] a_IRS;
			a_IRS = nullptr;
		}

		if(b_IRS != nullptr)
		{
			delete[] b_IRS;
			b_IRS = nullptr;
		};
		if(c_IRS != nullptr)
		{
			delete[] c_IRS;
			c_IRS = nullptr;
		}
		if(R_IRS != nullptr)
		{
			delete[] R_IRS;
			R_IRS = nullptr;
		}
		if(Rs_IRS != nullptr)
		{
			delete[] Rs_IRS;
			Rs_IRS = nullptr;
		}
	}
	if(timeDiscrType==TimeDiscrTypes::LU_SGS)
	{
		if(delt_cv != nullptr)
		{
			for (i=0;i<IM;i++)
			{
				for (j=0;j<JM;j++)
				{
					delete[] delt_cv[i][j];
				}
				delete[] delt_cv[i];
			}
			delete[] delt_cv;
			delt_cv = nullptr;
		}
		if(delt_ev != nullptr)
		{
			for (i=0;i<IM;i++)
			{
				for (j=0;j<JM;j++)
				{
					delete[] delt_ev[i][j];
				}
				delete[] delt_ev[i];
			}
			delete[] delt_ev;
			delt_ev = nullptr;
		}
	}

	if(CFluidProps::turbulenceType==TurbulenceTypes::SA)
	{
		if(Cell_ev_old != nullptr)
		{
			for (i=0;i<IM;i++)
			{
				for (j=0;j<JM;j++)
				{
					delete[] Cell_ev_old[i][j];
				}
				delete[] Cell_ev_old[i];
			}
			delete[] Cell_ev_old;
			Cell_ev_old = nullptr;
		}
	}
}

void CTimeDiscr::AllocateMemory(const CGeometry &Geometry)
{
	int i, j;
	IM = Geometry.IM;
	JM = Geometry.JM;
	KM = Geometry.KM;

	Cell_cv_old = new TConsVars **[IM];
	dt			= new REAL  **[IM];
	for (i=0;i<IM;i++)
	{
		Cell_cv_old[i] = new TConsVars *[JM];
		dt[i]		   = new REAL  *[JM];
	}
	for (i=0;i<IM;i++)
	{
		for (j=0;j<JM;j++)
		{
			Cell_cv_old[i][j] = new TConsVars [KM];
			dt[i][j]          = new REAL [KM];
		}
	}
	if(string(timeDiscrGeneralName)==string("Runge_Kutta"))
	{	
		Cell_Fv_old = new TConsVars **[IM];
		Cell_Fd_old = new TConsVars **[IM];
		for (i=0;i<IM;i++)
		{		
			Cell_Fv_old[i] = new TConsVars *[JM];
			Cell_Fd_old[i] = new TConsVars *[JM];
		}
		for (i=0;i<IM;i++)
		{
			for (j=0;j<JM;j++)
			{			
				Cell_Fv_old[i][j] = new TConsVars [KM];
				Cell_Fd_old[i][j] = new TConsVars [KM];
			}
		}			
	}
	if(IRS ==1)
	{
		int CMax;
		epsI = new REAL [IM+1];
		epsJ = new REAL [JM+1];
		epsK = new REAL [KM+1];
		CMax = int(MAX3(REAL(IM),REAL(JM),REAL(KM)));
		a_IRS = new REAL [CMax+1];
		b_IRS = new REAL [CMax+1];
		c_IRS = new REAL [CMax+1];
		R_IRS = new TConsVars [CMax+1];
		Rs_IRS = new TConsVars [CMax+1];
	}
	if(timeDiscrType==TimeDiscrTypes::LU_SGS)
	{
		delt_cv	= new TConsVars **[IM];
		delt_ev	= new TEddyVars **[IM];
		for (i=0;i<IM;i++)
		{
			delt_cv[i]  = new TConsVars *[JM];
			delt_ev[i]  = new TEddyVars *[JM];
		}
		for (i=0;i<IM;i++)
		{
			for (j=0;j<JM;j++)
			{
				delt_cv[i][j]  = new TConsVars [KM];
				delt_ev[i][j]  = new TEddyVars [KM];
			}
		}
	}
	if(CFluidProps::turbulenceType==TurbulenceTypes::SA)
	{
		Cell_ev_old = new TEddyVars **[IM];
		for (i=0;i<IM;i++)
		{
			Cell_ev_old[i] = new TEddyVars *[JM];
		}
		for (i=0;i<IM;i++)
		{
			for (j=0;j<JM;j++)
			{
				Cell_ev_old[i][j] = new TEddyVars[KM];
			}
		}
	}
}