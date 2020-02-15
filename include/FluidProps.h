

#ifndef FLUIDPROPS_H
#define FLUIDPROPS_H

#include "defs.h"
#include <string>
#include "Geometry.h"


class CFluidProps
{
public:
	static bool			   limiter_vars;
	static int             turbulenceScheme;
	static EquationTypes   equationType;
	static TurbulenceTypes turbulenceType;
	static string		   equationName,  turbulenceName;

	static REAL  Gamma,			/**< ratio of specific heat coefficients */
				 Pr_lam,		/**< laminar Prandtl number */
				 Pr_tur,		/**< turbulence Prandtl number */
				 Rgas,			/**< Gas_Constant */
				 Ma0;			/**< use for initial*/

	static REAL  Vel_ref,		/**  reference V */
				 dens_ref,
				 Temp_ref,
				 Press_ref,
				 mue_ref,
				 Ma_ref, 		/**< reference Ma, used to initialization  */	 
				 Re_ref;			/**< reylod number  */
				 
	 REAL  rpm,			/**< rpm---revolutions per minute */
		   w_rot;

	int IM, JM, KM, NDummy;

	int Flag_low_dt;												// use to limit flow variables
	REAL LdensMin, LdensMax, LpressMin, LpressMax, LVvelMax;		// use to limit flow variables
				 
	TConsVars	 ***Cell_cv,		// conservation variables
				 ***DummyIS_cv,
				 ***DummyIP_cv,
				 ***DummyJS_cv,
				 ***DummyJP_cv,
				 ***DummyKS_cv,
				 ***DummyKP_cv;

	TPrimVars    ***Cell_pv,		// primitive variables
				 ***DummyIS_pv,
				 ***DummyIP_pv,
				 ***DummyJS_pv,
				 ***DummyJP_pv,
				 ***DummyKS_pv,
				 ***DummyKP_pv;


	TPrimVars    ***gradx_I,			// gradient
				 ***grady_I,
				 ***gradz_I,
				 ***gradx_J,
				 ***grady_J,
				 ***gradz_J,
			     ***gradx_K,			     
			     ***grady_K,		    
			     ***gradz_K;
	REAL		 ***gradx_I_ev,
				 ***grady_I_ev,
				 ***gradz_I_ev,
				 ***gradx_J_ev,
				 ***grady_J_ev,
				 ***gradz_J_ev,
				 ***gradx_K_ev,
				 ***grady_K_ev,
				 ***gradz_K_ev;

	TVisVars	 ***Cell_tur;

	TEddyVars	 ***Cell_ev;	 // turbulence eddy viscostity

	REAL		 ***lambda_cI,	// convective spectral radius
				 ***lambda_cJ,
				 ***lambda_cK,
				 ***lambda_vI,	// viscous spectral radius
				 ***lambda_vJ,
				 ***lambda_vK;
			     
	REAL		 ***Omega;		// curl
	//REAL		 ***uvel_e,
	//			 ***vvel_e;		// embroil velocity
	// function
	CFluidProps();
	~CFluidProps();
	// void AllocateMemory(int IM, int JM, int KM);
	void AllocateMemory(const CGeometry &Geometry);
	void PrimitiveVars();
	void SpectralRadii(const CGeometry &Geometry);
	void Gradient(const CGeometry &Geometry);
	void Curl(const CGeometry &Geometry);
	//void Bulk_Velocity(const CGeometry &Geometry);		

	void Limit_Vars();
	void Limit_SA();


	REAL SutherLand(REAL temp)
	  {
		REAL mue;
		mue = 1.0f/Re_ref*POW(temp, 3.0f/2.0f) * (1.0f + 110.4f/288.15f)/(temp + 110.4f/288.15f);
		return mue;
	  }

private:
	CFluidProps(const CFluidProps &);					// copy constructor 
	CFluidProps & operator = (const CFluidProps &);		// assignment operator

	REAL sign(REAL a, REAL b)
	{
		REAL c;
		if(b<0)
			return c = -1.0f*ABS(a);
		else
			return c = ABS(a);
	}
};


#endif

 