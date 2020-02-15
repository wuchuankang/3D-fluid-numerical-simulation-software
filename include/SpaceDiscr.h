

#ifndef SPACEDISCR_H
#define SPACEDISCR_H

#include "defs.h"
#include <string>
#include "Geometry.h"
#include "BndConds.h"
#include "FluidProps.h"

class CSpaceDiscr
{
public:
	static SpaceDiscrTypes spaceDiscrType;
	static ReconstructTypes reconstructType;
	static ReconstructVars  reconstructVar;
	static string  spaceDiscrName,  reconstructTypeName, reconstructVarName;
	static REAL epsEntr, k2, k4;	/**< entropy correction coefficient */ 
									// k2, k4 used for artificial Dissipation
	static int  DiscrScheme,		/**< JST----------------100
										 Steger-warming-----200
										 Roe----------------300  */
				ReconstructScheme,  /**< None---------------0 do not use the reconstruction
										 NND----------------100
										 MUSCL--------------200  */
				ReconstructVar;	/**<conservation--------100
										primitive-----------200
										charactristic-------300	 */			
	int IM, JM, KM;
	int Mesh_ID;
   
	TConsVars  ***Res_cv,
			   ***Cell_Fc,	// convective flux at cell
			   ***Cell_Fv,	// viscous flux at cell
			   ***Cell_Fd,	// artificial dissipation flux at cell
			   ***FaceI_Fc,	// convective flux at face
			   ***FaceJ_Fc,
			   ***FaceK_Fc,
			   ***FaceI_Fv,
			   ***FaceJ_Fv,
			   ***FaceK_Fv,
			   ***FaceI_Fd,
			   ***FaceJ_Fd,
			   ***FaceK_Fd;
			   
  TEddyVars	   ***Res_ev,
			   ***Cell_Fc_tur,
			   ***Cell_Fv_tur,
			   ***FaceI_Fc_tur,	// convective flux at face
			   ***FaceJ_Fc_tur,
			   ***FaceK_Fc_tur,
			   ***FaceI_Fv_tur,
			   ***FaceJ_Fv_tur,
			   ***FaceK_Fv_tur,
			   ***FaceI_Fv1_tur,
			   ***FaceJ_Fv1_tur,
			   ***FaceK_Fv1_tur;
	// function
	CSpaceDiscr();
	~CSpaceDiscr();
	void AllocateMemory(const CGeometry &Geometry);
	void ConvectiveFlux(const CGeometry &Geometry, const CFluidProps &FluidProp);
	void ViscousFlux(const CGeometry &Geometry, CFluidProps &FluidProp);
	void JST_Dissip(const CFluidProps &FluidProp);

	void WallFlux(const TWall &Wall, const CGeometry & Geometry, const CFluidProps & FluidProp);

	void SA_Model(const CGeometry &Geometry, CFluidProps &FluidProp);

private:
	void JST_Convective(const CGeometry &Geometry, const CFluidProps &FluidProp);
	void Upwind_Convective(const CGeometry &Geometry, const CFluidProps &FluidProp);
	 	
	REAL sensorI(int i, int j, int k, const CFluidProps &FluidProp);	// used for JST_dissip()
	REAL sensorJ(int i, int j, int k, const CFluidProps &FluidProp);
	REAL sensorK(int i, int j, int k, const CFluidProps &FluidProp);

	// entropy revised function used in Roe
	REAL EntropyCorr( REAL z, REAL d )
	{
		if (z > d)
			return z;
		else
			return 0.5f*(z*z+d*d)/d;
	}

	// JST scheme
	void JST_Scheme(int, int, int, int, const CGeometry &, const CFluidProps &, TConsVars &);
	//void RoeScheme(int, int, int, int, const CGeometry &,TConsVars &, TConsVars &, TConsVars &, TConsVars &, TConsVars &, TConsVars &, TConsVars &);
	void RoeScheme(int , int , int , int , const CGeometry &, const CFluidProps &, TNode &, TNode &, TNode &,TConsVars &, TConsVars &, TConsVars &);
	// StegerWarimg
	void StegerWarimgScheme(int, int, int, int, const CGeometry &, const CFluidProps &, TNode &, TNode &, TNode &, TConsVars &, TConsVars &, TConsVars &);
	// L_F_scheme
	void L_F_Scheme(int, int, int, int, const CGeometry &, const CFluidProps &, TNode &, TNode &, TNode &, TConsVars &, TConsVars &, TConsVars &);
	// Van leer
	void  Van_Leer_Scheme(int, int, int, int, const CGeometry &, const CFluidProps &, TNode &, TNode &, TNode &, TConsVars &, TConsVars &, TConsVars &);
	// AUSM_plus
	void AUSM_Plus_Scheme(int, int, int, int, const CGeometry &, const CFluidProps &, TNode &, TNode &, TNode &, TConsVars &, TConsVars &, TConsVars &);
	// AUSM_PW
	void AUSM_PW_Scheme(int, int, int, int, const CGeometry &, const CFluidProps &, TNode &, TNode &, TNode &, TConsVars &, TConsVars &, TConsVars &);
	// HLLC
	void HLLC_Scheme(int, int, int, int, const CGeometry &, const CFluidProps &, TNode &, TNode &, TNode &, TConsVars &, TConsVars &, TConsVars &);
	// conservation reconstruction
	void ConservationReconstruct(TConsVars &, TConsVars &, TConsVars &, TConsVars &, TConsVars &, TConsVars &);
	// primitive reconstruction
	void PrimitiveReconstruct(TPrimVars &, TPrimVars &, TPrimVars &, TPrimVars &, TConsVars &, TConsVars &);
	// characteristic reconstruction
	//void CharacterReconstruct(TConsVars &, TConsVars &, TConsVars &, TConsVars &, TConsVars &, TConsVars &);
	void CharacterReconstruct(TPrimVars &, TPrimVars &, TConsVars &, TConsVars &, TConsVars &, TConsVars &, TConsVars &, TConsVars &);

	// coordinate rotation
	void Coord_Rotate(int, int, int, int, const CGeometry &, const CFluidProps &, TNode &, TNode &, TNode &, TConsVars &, TConsVars &, TConsVars &, TConsVars &, REAL&);
	// minmod() function, used for NND scheme
	TConsVars minmod(const TConsVars &A, const TConsVars &B);
	// override minmod()
	TPrimVars minmod(const TPrimVars &A, const TPrimVars &B);

	CSpaceDiscr(const CSpaceDiscr &);					// copy constructor 
	CSpaceDiscr & operator = (const CSpaceDiscr &);		// assignment operator
};

#endif