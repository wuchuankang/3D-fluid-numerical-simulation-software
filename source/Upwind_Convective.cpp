
#include "defs.h"
#include "SpaceDiscr.h"

void CSpaceDiscr::Upwind_Convective(const CGeometry &Geometry, const CFluidProps &FluidProp)
{
	int i, j, k, flag;	  

	TConsVars ***Cell_cv	= FluidProp.Cell_cv;
	TConsVars ***DummyIS_cv = FluidProp.DummyIS_cv;
	TConsVars ***DummyIP_cv = FluidProp.DummyIP_cv;
	TConsVars ***DummyJS_cv = FluidProp.DummyJS_cv;
	TConsVars ***DummyJP_cv = FluidProp.DummyJP_cv;
	TConsVars ***DummyKS_cv = FluidProp.DummyKS_cv;
	TConsVars ***DummyKP_cv = FluidProp.DummyKP_cv;

	TPrimVars ***Cell_pv	= FluidProp.Cell_pv;
	TPrimVars ***DummyIS_pv = FluidProp.DummyIS_pv;
	TPrimVars ***DummyIP_pv = FluidProp.DummyIP_pv;
	TPrimVars ***DummyJS_pv = FluidProp.DummyJS_pv;
	TPrimVars ***DummyJP_pv = FluidProp.DummyJP_pv;
	TPrimVars ***DummyKS_pv = FluidProp.DummyKS_pv;
	TPrimVars ***DummyKP_pv = FluidProp.DummyKP_pv;

	TConsVars CellF_cv, CellB_cv, CellFF_cv, CellBB_cv, cv_L, cv_R, Fc;
	TPrimVars CellF_pv, CellB_pv, CellFF_pv, CellBB_pv;
	TNode norm, tau1, tau2;
	//----------------------------------------------------------------------------------
#pragma omp parallel num_threads(NT) default(shared) firstprivate(i,j,k, CellF_cv, CellB_cv, CellFF_cv, CellBB_cv, cv_L, cv_R, Fc, CellF_pv, CellB_pv, CellFF_pv, CellBB_pv, norm, tau1, tau2)
	{
#pragma omp for nowait
		for(i=1;i<=IM;i++)
		{
			for(j=1;j<=JM-1;j++)
			{
				for(k=1;k<=KM-1;k++)
				{
					if(i==1)
					{
						CellFF_cv = DummyIS_cv[2][j][k];
						CellF_cv  = DummyIS_cv[1][j][k];
						CellB_cv  = Cell_cv[i][j][k];
						CellBB_cv = Cell_cv[i+1][j][k];

						CellFF_pv = DummyIS_pv[2][j][k];
						CellF_pv  = DummyIS_pv[1][j][k];
						CellB_pv  = Cell_pv[i][j][k];
						CellBB_pv = Cell_pv[i+1][j][k];
					}
					else if(i==2)
					{
						CellFF_cv = DummyIS_cv[1][j][k];
						CellF_cv  = Cell_cv[i-1][j][k];
						CellB_cv  = Cell_cv[i][j][k];
						CellBB_cv = Cell_cv[i+1][j][k];

						CellFF_pv = DummyIS_pv[1][j][k];
						CellF_pv  = Cell_pv[i-1][j][k];
						CellB_pv  = Cell_pv[i][j][k];
						CellBB_pv = Cell_pv[i+1][j][k];
					}
					else if(i==IM-1)
					{
						CellFF_cv = Cell_cv[i-2][j][k];
						CellF_cv  = Cell_cv[i-1][j][k];
						CellB_cv  = Cell_cv[i][j][k];
						CellBB_cv = DummyIP_cv[1][j][k];

						CellFF_pv = Cell_pv[i-2][j][k];
						CellF_pv  = Cell_pv[i-1][j][k];
						CellB_pv  = Cell_pv[i][j][k];
						CellBB_pv = DummyIP_pv[1][j][k];
					}
					else if(i==IM)
					{
						CellFF_cv = Cell_cv[i-2][j][k];
						CellF_cv  = Cell_cv[i-1][j][k];
						CellB_cv  = DummyIP_cv[1][j][k];
						CellBB_cv = DummyIP_cv[2][j][k];

						CellFF_pv = Cell_pv[i-2][j][k];
						CellF_pv  = Cell_pv[i-1][j][k];
						CellB_pv  = DummyIP_pv[1][j][k];
						CellBB_pv = DummyIP_pv[2][j][k];
					}				
					else
					{
						CellFF_cv = Cell_cv[i-2][j][k];
						CellF_cv  = Cell_cv[i-1][j][k];
						CellB_cv  = Cell_cv[i][j][k];
						CellBB_cv = Cell_cv[i+1][j][k];

						CellFF_pv = Cell_pv[i-2][j][k];
						CellF_pv  = Cell_pv[i-1][j][k];
						CellB_pv  = Cell_pv[i][j][k];
						CellBB_pv = Cell_pv[i+1][j][k];
					}
					//NND scheme
					if(reconstructVar == ReconstructVars::ConsVars)
						ConservationReconstruct(CellFF_cv, CellF_cv, CellB_cv, CellBB_cv, cv_L, cv_R);

					if(reconstructVar == ReconstructVars::PrimVars)
						PrimitiveReconstruct(CellFF_pv, CellF_pv, CellB_pv, CellBB_pv, cv_L, cv_R);

					if(reconstructVar == ReconstructVars::CharactVars)
						CharacterReconstruct(CellF_pv, CellB_pv, CellFF_cv, CellF_cv, CellB_cv, CellBB_cv, cv_L, cv_R);

					//-----------------------------------------------------
					flag = 100;

					if(spaceDiscrType==SpaceDiscrTypes::StegerWarming)
						StegerWarimgScheme(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Fc);

					if(spaceDiscrType==SpaceDiscrTypes::Roe)
						RoeScheme(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Fc);

					if(spaceDiscrType==SpaceDiscrTypes::Van_Leer)
						Van_Leer_Scheme(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Fc);

					if(spaceDiscrType==SpaceDiscrTypes::AUSM_Plus)
						AUSM_Plus_Scheme(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Fc);

					if(spaceDiscrType==SpaceDiscrTypes::AUSM_PW)
						AUSM_PW_Scheme(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Fc);

					if(spaceDiscrType==SpaceDiscrTypes::HLLC)
						HLLC_Scheme(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Fc);

					if(spaceDiscrType==SpaceDiscrTypes::L_F)
						L_F_Scheme(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Fc);

					FaceI_Fc[i][j][k].dens = Fc.dens * Geometry.FaceI[i][j][k].s;
					FaceI_Fc[i][j][k].xmom = (Fc.xmom*norm.x + Fc.ymom*tau1.x +Fc.zmom*tau2.x) * Geometry.FaceI[i][j][k].s;
					FaceI_Fc[i][j][k].ymom = (Fc.xmom*norm.y + Fc.ymom*tau1.y +Fc.zmom*tau2.y) * Geometry.FaceI[i][j][k].s;
					FaceI_Fc[i][j][k].zmom = (Fc.xmom*norm.z + Fc.ymom*tau1.z +Fc.zmom*tau2.z) * Geometry.FaceI[i][j][k].s;
					FaceI_Fc[i][j][k].ener = Fc.ener * Geometry.FaceI[i][j][k].s;
				}
			}
		}

#pragma omp for nowait
		for(i=1;i<=IM-1;i++)
		{
			for(j=1;j<=JM;j++)		
			{
				for(k=1;k<=KM-1;k++)
				{
					if(j==1)
					{
						CellFF_cv = DummyJS_cv[i][2][k];
						CellF_cv  = DummyJS_cv[i][1][k];
						CellB_cv  = Cell_cv[i][j][k];
						CellBB_cv = Cell_cv[i][j+1][k];

						CellFF_pv = DummyJS_pv[i][2][k];
						CellF_pv  = DummyJS_pv[i][1][k];
						CellB_pv  = Cell_pv[i][j][k];
						CellBB_pv = Cell_pv[i][j+1][k];
					}
					else if(j==2)
					{
						CellFF_cv = DummyJS_cv[i][1][k];
						CellF_cv  = Cell_cv[i][j-1][k];
						CellB_cv  = Cell_cv[i][j][k];
						CellBB_cv = Cell_cv[i][j+1][k];

						CellFF_pv = DummyJS_pv[i][1][k];
						CellF_pv  = Cell_pv[i][j-1][k];
						CellB_pv  = Cell_pv[i][j][k];
						CellBB_pv = Cell_pv[i][j+1][k];
					}
					else if(j==JM-1)
					{
						CellFF_cv = Cell_cv[i][j-2][k];
						CellF_cv  = Cell_cv[i][j-1][k];
						CellB_cv  = Cell_cv[i][j][k];
						CellBB_cv = DummyJP_cv[i][1][k];

						CellFF_pv = Cell_pv[i][j-2][k];
						CellF_pv  = Cell_pv[i][j-1][k];
						CellB_pv  = Cell_pv[i][j][k];
						CellBB_pv = DummyJP_pv[i][1][k];
					}
					else if(j==JM)
					{
						CellFF_cv = Cell_cv[i][j-2][k];
						CellF_cv  = Cell_cv[i][j-1][k];
						CellB_cv  = DummyJP_cv[i][1][k];
						CellBB_cv = DummyJP_cv[i][2][k];

						CellFF_pv = Cell_pv[i][j-2][k];
						CellF_pv  = Cell_pv[i][j-1][k];
						CellB_pv  = DummyJP_pv[i][1][k];
						CellBB_pv = DummyJP_pv[i][2][k];
					}				
					else
					{
						CellFF_cv = Cell_cv[i][j-2][k];
						CellF_cv  = Cell_cv[i][j-1][k];
						CellB_cv  = Cell_cv[i][j][k];
						CellBB_cv = Cell_cv[i][j+1][k];

						CellFF_pv = Cell_pv[i][j-2][k];
						CellF_pv  = Cell_pv[i][j-1][k];
						CellB_pv  = Cell_pv[i][j][k];
						CellBB_pv = Cell_pv[i][j+1][k];
					}

					if(reconstructVar == ReconstructVars::ConsVars)
						ConservationReconstruct(CellFF_cv, CellF_cv, CellB_cv, CellBB_cv, cv_L, cv_R);

					if(reconstructVar == ReconstructVars::PrimVars)
						PrimitiveReconstruct(CellFF_pv, CellF_pv, CellB_pv, CellBB_pv, cv_L, cv_R);

					if(reconstructVar == ReconstructVars::CharactVars)
						CharacterReconstruct(CellF_pv, CellB_pv, CellFF_cv, CellF_cv, CellB_cv, CellBB_cv, cv_L, cv_R);

					//-----------------------------------------------------
					flag = 200;							
					if(spaceDiscrType==SpaceDiscrTypes::StegerWarming)
						StegerWarimgScheme(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Fc);

					if(spaceDiscrType==SpaceDiscrTypes::Roe)
						RoeScheme(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Fc);

					if(spaceDiscrType==SpaceDiscrTypes::Van_Leer)
						Van_Leer_Scheme(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Fc);

					if(spaceDiscrType==SpaceDiscrTypes::AUSM_Plus)
						AUSM_Plus_Scheme(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Fc);

					if(spaceDiscrType==SpaceDiscrTypes::AUSM_PW)
						AUSM_PW_Scheme(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Fc);

					if(spaceDiscrType==SpaceDiscrTypes::HLLC)
						HLLC_Scheme(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Fc);

					if(spaceDiscrType==SpaceDiscrTypes::L_F)
						L_F_Scheme(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Fc);

					FaceJ_Fc[i][j][k].dens = Fc.dens * Geometry.FaceJ[i][j][k].s;
					FaceJ_Fc[i][j][k].xmom = (Fc.xmom*norm.x + Fc.ymom*tau1.x +Fc.zmom*tau2.x) * Geometry.FaceJ[i][j][k].s;
					FaceJ_Fc[i][j][k].ymom = (Fc.xmom*norm.y + Fc.ymom*tau1.y +Fc.zmom*tau2.y) * Geometry.FaceJ[i][j][k].s;
					FaceJ_Fc[i][j][k].zmom = (Fc.xmom*norm.z + Fc.ymom*tau1.z +Fc.zmom*tau2.z) * Geometry.FaceJ[i][j][k].s;
					FaceJ_Fc[i][j][k].ener = Fc.ener * Geometry.FaceJ[i][j][k].s;
				}
			}
		}

#pragma omp for nowait
		for(i=1;i<=IM-1;i++)
		{
			for(j=1;j<=JM-1;j++)		
			{
				for(k=1;k<=KM;k++)
				{
					if(k==1)
					{
						CellFF_cv = DummyKS_cv[i][j][2];
						CellF_cv  = DummyKS_cv[i][j][1];
						CellB_cv  = Cell_cv[i][j][k];
						CellBB_cv = Cell_cv[i][j][k+1];

						CellFF_pv = DummyKS_pv[i][j][2];
						CellF_pv  = DummyKS_pv[i][j][1];
						CellB_pv  = Cell_pv[i][j][k];
						CellBB_pv = Cell_pv[i][j][k+1];
					}
					else if(k==2)
					{
						CellFF_cv = DummyKS_cv[i][j][1];
						CellF_cv  = Cell_cv[i][j][k-1];
						CellB_cv  = Cell_cv[i][j][k];
						CellBB_cv = Cell_cv[i][j][k+1];

						CellFF_pv = DummyKS_pv[i][j][1];
						CellF_pv  = Cell_pv[i][j][k-1];
						CellB_pv  = Cell_pv[i][j][k];
						CellBB_pv = Cell_pv[i][j][k+1];
					}
					else if(k==KM-1)
					{
						CellFF_cv = Cell_cv[i][j][k-2];
						CellF_cv  = Cell_cv[i][j][k-1];
						CellB_cv  = Cell_cv[i][j][k];
						CellBB_cv = DummyKP_cv[i][j][1];

						CellFF_pv = Cell_pv[i][j][k-2];
						CellF_pv  = Cell_pv[i][j][k-1];
						CellB_pv  = Cell_pv[i][j][k];
						CellBB_pv = DummyKP_pv[i][j][1];
					}
					else if(k==KM)
					{
						CellFF_cv = Cell_cv[i][j][k-2];
						CellF_cv  = Cell_cv[i][j][k-1];
						CellB_cv  = DummyKP_cv[i][j][1];
						CellBB_cv = DummyKP_cv[i][j][2];

						CellFF_pv = Cell_pv[i][j][k-2];
						CellF_pv  = Cell_pv[i][j][k-1];
						CellB_pv  = DummyKP_pv[i][j][1];
						CellBB_pv = DummyKP_pv[i][j][2];
					}				
					else
					{
						CellFF_cv = Cell_cv[i][j][k-2];
						CellF_cv  = Cell_cv[i][j][k-1];
						CellB_cv  = Cell_cv[i][j][k];
						CellBB_cv = Cell_cv[i][j][k+1];

						CellFF_pv = Cell_pv[i][j][k-2];
						CellF_pv  = Cell_pv[i][j][k-1];
						CellB_pv  = Cell_pv[i][j][k];
						CellBB_pv = Cell_pv[i][j][k+1];
					}

					if(reconstructVar == ReconstructVars::ConsVars)
						ConservationReconstruct(CellFF_cv, CellF_cv, CellB_cv, CellBB_cv, cv_L, cv_R);

					if(reconstructVar == ReconstructVars::PrimVars)
						PrimitiveReconstruct(CellFF_pv, CellF_pv, CellB_pv, CellBB_pv, cv_L, cv_R);

					if(reconstructVar == ReconstructVars::CharactVars)
						CharacterReconstruct(CellF_pv, CellB_pv, CellFF_cv, CellF_cv, CellB_cv, CellBB_cv, cv_L, cv_R);

					//-----------------------------------------------------
					flag = 300;						
					if(spaceDiscrType==SpaceDiscrTypes::StegerWarming)
						StegerWarimgScheme(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Fc);

					if(spaceDiscrType==SpaceDiscrTypes::Roe)
						RoeScheme(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Fc);

					if(spaceDiscrType==SpaceDiscrTypes::Van_Leer)
						Van_Leer_Scheme(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Fc);

					if(spaceDiscrType==SpaceDiscrTypes::AUSM_Plus)
						AUSM_Plus_Scheme(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Fc);

					if(spaceDiscrType==SpaceDiscrTypes::AUSM_PW)
						AUSM_PW_Scheme(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Fc);

					if(spaceDiscrType==SpaceDiscrTypes::HLLC)
						HLLC_Scheme(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Fc);

					if(spaceDiscrType==SpaceDiscrTypes::L_F)
						L_F_Scheme(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Fc);


					FaceK_Fc[i][j][k].dens = Fc.dens * Geometry.FaceK[i][j][k].s;
					FaceK_Fc[i][j][k].xmom = (Fc.xmom*norm.x + Fc.ymom*tau1.x +Fc.zmom*tau2.x) * Geometry.FaceK[i][j][k].s;
					FaceK_Fc[i][j][k].ymom = (Fc.xmom*norm.y + Fc.ymom*tau1.y +Fc.zmom*tau2.y) * Geometry.FaceK[i][j][k].s;
					FaceK_Fc[i][j][k].zmom = (Fc.xmom*norm.z + Fc.ymom*tau1.z +Fc.zmom*tau2.z) * Geometry.FaceK[i][j][k].s;
					FaceK_Fc[i][j][k].ener = Fc.ener * Geometry.FaceK[i][j][k].s;
				}
			}
		}
	}
}


