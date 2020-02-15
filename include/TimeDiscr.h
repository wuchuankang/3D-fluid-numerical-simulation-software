

#ifndef TIMEDISCR_H
#define TIMEDISCR_H

#include "defs.h"
#include <string>
#include "Geometry.h"
#include "FluidProps.h"
#include "SpaceDiscr.h"


using namespace std;

class CTimeDiscr
{
public:
	static TimestepTypes  timestepType;
	static TimeDiscrTypes timeDiscrType;
	static int TimeDiscrScheme;
	static string  timestepName, timeDiscrName, timeDiscrGeneralName;
	static REAL CFL, C_timestep;		// time step coefficient	
	static REAL relax_factor;
	static bool IRS;

	int IM, JM, KM;

	TConsVars	***Cell_cv_old;
	TConsVars	***delt_cv;
	TConsVars	***Cell_Fv_old;
	TConsVars	***Cell_Fd_old;
	TEddyVars	***Cell_ev_old;
	TEddyVars	***delt_ev;
	REAL		***dt;

	REAL		*epsI, *epsJ, *epsK, *a_IRS, // 用于求隐式残差光顺
				*b_IRS, *c_IRS;
	TConsVars	*R_IRS, *Rs_IRS;	// 用于求隐式残差光顺

	// function
	CTimeDiscr();
	~CTimeDiscr();
	void AllocateMemory(const CGeometry &Geometry);
	void TimeStep(const CGeometry &Geometry, CFluidProps &FluidProp);
	void CentralIRS(const CFluidProps &FluidProp, CSpaceDiscr &SpaceDiscr);
	void UpwindIRS(const CGeometry &Geometry, const CFluidProps &FluidProp, CSpaceDiscr &SpaceDiscr);

	void  delt_Fn_up(int i, int j, int k, int flag, TConsVars &delt_Fn_cv, TEddyVars &delt_Fn_ev,const CGeometry &Geometry, const CFluidProps &FluidProp);		// used for LU_SGS 
	void  delt_Fn_down(int i, int j, int k, int flag, TConsVars &delt_Fn_cv, TEddyVars &delt_Fn_ev,const CGeometry &Geometry, const CFluidProps &FluidProp);

private:
	CTimeDiscr(const CTimeDiscr &);						// copy constructor 
	CTimeDiscr & operator = (const CTimeDiscr &);		// assignment operator

	void Thomas(int TMAX, REAL *a, REAL *b, REAL *c, TConsVars *d, TConsVars *x); // 用于求隐式残差光顺
};


#endif