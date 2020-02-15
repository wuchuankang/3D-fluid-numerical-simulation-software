
#ifndef MESH_H
#define MESH_H

#include "Blocks.h"
#include "BndGeom.h"
#include "defs.h"

class CMesh
{
public:
	static int NB;					  /**< number of blocks **/
	int Mesh_ID;
	static TurboTypes TurboType;
	static InitialTypes InitialType;
	bool UseForceFunction;

	REAL TotalVol;
						      /**< which the iteration process is stopped) */
	CBlocks		*Block;
	CBndGeom	BndGeom;

	TConsVars ****cv_2h,		//use for multigrids
			  ****delt_cv,
			  ****Res_2h,
			  ****Qf_2h;

		// function
	CMesh();
	~CMesh();
	void AllocateMemory();
	void AllocateMemory_Pro_Inter_Vars();
	void Geometry();
	void Initial();
	void PrimitiveVars();
	void Time_Marching();
	void ResFlux();
	void ComputeForceFunction();
	void AddForceFunction();
	void Res_Calculate(TConsVars &res, REAL &res_ev);
	void CornerTreatment();
	void WallDistanceSA();		/**< wall distance for SA model  */

	void Limit_Vars();

	void Initial_Interpolate(CBlocks *);

private:


	void Spectral_curl_grad();

	void TimeStep();

	void ConvectiveFlux();

	void WallFlux();

	void Dissipation();

	void ViscousFlux();

	void IRS();

	void Turbulence_Model();

	void Runge_Kutta3_Marching();
	void Runge_Kutta4_Marching();
	void Runge_Kutta5_Marching();
	void Runge_Kutta5_Hybrid_Marching();

	void LU_SGS_Marching();
	// BoundaryCondition() and WallDistanceSA() have been computed already, so we need not compute this in the SolveAllBlockVars.cpp
	void BoundaryCondition();

	void BndCondInterior(CBlocks &A, CBlocks &B, const TInterior &Interior);
	void BndCondPeriodic(CBlocks &A, CBlocks &B, const TPeriodic &Periodic);
								  
	void BndCondMixingPlane(const TMixingPlane &MixingPlane);
								  
	void ShadowOut(CBlocks &A, int fc, int nd, int is, int ie, int js, int je, int ks, int ke);
	void Trans(int ori, REAL dr);
	void ShadowIn(CBlocks &B, int fc, int nd, int is, int ie, int js, int je, int ks, int ke);

	void LinearFit(int, int, REAL*, REAL*, REAL*, REAL*);

		
	CMesh(const CMesh &);					// copy constructor 
	CMesh& operator = (const CMesh &);		// assignment operator
};


#endif

