
#include "SpaceDiscr.h"

using namespace std;


void CSpaceDiscr::ConservationReconstruct(TConsVars &CellFF_cv, TConsVars &CellF_cv, TConsVars &CellB_cv, TConsVars &CellBB_cv, TConsVars &cv_L, TConsVars &cv_R)
{
	REAL epsilon=1E-6f, eta=0.05f;			//eta是MUSCL2o中的耗散系数 >0, 值越大，耗散越大
	REAL k2c=1.0f, k2u=-1.0f;
	TConsVars r1, r2, f1, f;
	TConsVars SL_limiter, SR_limiter;
	REAL Gamma = CFluidProps::Gamma;
	REAL lim_zero=1E-20f;
	REAL pressL, pressR;

	if(reconstructType == ReconstructTypes::NND)
	{
		cv_L = CellF_cv + 0.5f*minmod(CellB_cv-CellF_cv, CellF_cv-CellFF_cv);
		cv_R = CellB_cv - 0.5f*minmod(CellB_cv-CellF_cv, CellBB_cv-CellB_cv);
	}
	if(reconstructType == ReconstructTypes::MUSCL2u)
	{
		cv_L = CellF_cv + 0.25*((1.0f-k2u)*minmod(2.0f*(CellB_cv-CellF_cv), CellF_cv-CellFF_cv) + (1.0f+k2u)*minmod(2.0f*(CellF_cv-CellFF_cv), CellB_cv-CellF_cv));
		cv_R = CellB_cv - 0.25*((1.0f-k2u)*minmod(2.0f*(CellB_cv-CellF_cv), CellBB_cv-CellB_cv) + (1.0f+k2u)*minmod(2.0f*(CellBB_cv-CellB_cv), CellB_cv-CellF_cv));
	}
	if(reconstructType == ReconstructTypes::MUSCL2c)
	{							   
		cv_L = CellF_cv + 0.25*((1.0f-k2c)*minmod(2.0f*(CellB_cv-CellF_cv), CellF_cv-CellFF_cv) + (1.0f+k2c)*minmod(CellB_cv-CellF_cv, 2.0f*(CellF_cv-CellFF_cv)));
		cv_R = CellB_cv - 0.25*((1.0f-k2c)*minmod(2.0f*(CellB_cv-CellF_cv), CellBB_cv-CellB_cv) + (1.0f+k2c)*minmod(CellB_cv-CellF_cv, 2.0f*(CellBB_cv-CellB_cv)));
	}
	if(reconstructType == ReconstructTypes::MUSCL2o)
	{
		r1 = (CellF_cv-CellFF_cv+epsilon) / (CellB_cv-CellF_cv+epsilon);
		r2 = (CellB_cv-CellF_cv+epsilon) / (CellBB_cv-CellB_cv+epsilon);
		f1 = 1.0f-eta + (eta-0.55f)/2.0f/r2 + (eta+0.55f)/2.0f*r1;
		f.dens = MAX(0.0f, MIN3(2.0f, f1.dens, 2.0f*r1.dens));
		f.xmom = MAX(0.0f, MIN3(2.0f, f1.xmom, 2.0f*r1.xmom));
		f.ymom = MAX(0.0f, MIN3(2.0f, f1.ymom, 2.0f*r1.ymom));
		f.zmom = MAX(0.0f, MIN3(2.0f, f1.zmom, 2.0f*r1.zmom));
		f.ener = MAX(0.0f, MIN3(2.0f, f1.ener, 2.0f*r1.ener));
		cv_L = CellF_cv + 0.5f*f*(CellB_cv-CellF_cv);

		r1 = (CellBB_cv-CellB_cv+epsilon) / (CellB_cv-CellF_cv+epsilon);
		r2 = (CellB_cv-CellF_cv+epsilon) / (CellF_cv-CellFF_cv+epsilon);
		f1 = 1.0f-eta + (eta-0.55f)/2.0f/r2 + (eta+0.55f)/2.0f*r1;
		f.dens = MAX(0.0f, MIN3(2.0f, f1.dens, 2.0f*r1.dens));
		f.xmom = MAX(0.0f, MIN3(2.0f, f1.xmom, 2.0f*r1.xmom));
		f.ymom = MAX(0.0f, MIN3(2.0f, f1.ymom, 2.0f*r1.ymom));
		f.zmom = MAX(0.0f, MIN3(2.0f, f1.zmom, 2.0f*r1.zmom));
		f.ener = MAX(0.0f, MIN3(2.0f, f1.ener, 2.0f*r1.ener));
		cv_R = CellB_cv - 0.5f*f*(CellB_cv-CellF_cv);
	}
	if(reconstructType == ReconstructTypes::MUSCL3)
	{
		SL_limiter = (2.0f*(CellF_cv - CellFF_cv)*(CellB_cv - CellF_cv) + epsilon) / (POW(CellF_cv-CellFF_cv, 2.0f) + POW(CellB_cv-CellF_cv, 2.0f) + epsilon);
		SR_limiter = (2.0f*(CellBB_cv - CellB_cv)*(CellB_cv - CellF_cv) + epsilon) / (POW(CellBB_cv-CellB_cv, 2.0f) + POW(CellB_cv-CellF_cv, 2.0f) + epsilon);
		cv_L   = CellF_cv + SL_limiter/4.0f*((1.0f-SL_limiter/3.0f)*(CellF_cv-CellFF_cv) + (1.0f+SL_limiter/3.0f)*(CellB_cv-CellF_cv));	
		cv_R   = CellB_cv - SR_limiter/4.0f*((1.0f-SR_limiter/3.0f)*(CellBB_cv-CellB_cv) + (1.0f+SR_limiter/3.0f)*(CellB_cv-CellF_cv));
	}
	// 重构完成之后，要检查press 和 density 是否为负，如果是，则用一阶格式
	pressL = (Gamma-1.0f)*(cv_L.ener - 0.5f*(cv_L.xmom*cv_L.xmom+cv_L.ymom*cv_L.ymom+cv_L.zmom*cv_L.zmom)/cv_L.dens);
	pressR = (Gamma-1.0f)*(cv_R.ener - 0.5f*(cv_R.xmom*cv_R.xmom+cv_R.ymom*cv_R.ymom+cv_R.zmom*cv_R.zmom)/cv_R.dens);
	if(cv_L.dens<=lim_zero || pressL<=lim_zero)	
		cv_L = CellF_cv;
	if(cv_R.dens<=lim_zero || pressR<=lim_zero)	
		cv_R = CellB_cv;
}

void CSpaceDiscr::PrimitiveReconstruct(TPrimVars &CellFF_pv, TPrimVars &CellF_pv, TPrimVars &CellB_pv, TPrimVars &CellBB_pv, TConsVars &cv_L, TConsVars &cv_R)
{
	REAL epsilon=1E-6f, eta=0.2f;			//eta是MUSCL2o中的耗散系数 >0, 值越大，耗散越大
	REAL k2c=1.0, k2u=-1.0;
	TPrimVars pv_L, pv_R;
	TPrimVars r1, r2, f1, f;
	TPrimVars SL_limiter, SR_limiter;
	TPrimVars aL, bL, deltaL, aR, bR, deltaR;
	REAL Gamma = CFluidProps::Gamma;
	REAL lim_zero=1E-20f;

	if(reconstructType == ReconstructTypes::NND)
	{
		pv_L = CellF_pv + 0.5*minmod(CellB_pv-CellF_pv, CellF_pv-CellFF_pv);
		pv_R = CellB_pv - 0.5*minmod(CellB_pv-CellF_pv, CellBB_pv-CellB_pv);
	}
	if(reconstructType == ReconstructTypes::MUSCL2u)
	{
		pv_L = CellF_pv + 0.25f*((1.0f-k2u)*minmod(2.0f*(CellB_pv-CellF_pv), CellF_pv-CellFF_pv) + (1.0f+k2u)*minmod(2.0f*(CellF_pv-CellFF_pv), CellB_pv-CellF_pv));
		pv_R = CellB_pv - 0.25f*((1.0f-k2u)*minmod(2.0f*(CellB_pv-CellF_pv), CellBB_pv-CellB_pv) + (1.0f+k2u)*minmod(2.0f*(CellBB_pv-CellB_pv), CellB_pv-CellF_pv));
	}
	if(reconstructType == ReconstructTypes::MUSCL2c)
	{							   
		pv_L = CellF_pv + 0.25f*((1.0f-k2c)*minmod(2.0f*(CellB_pv-CellF_pv), CellF_pv-CellFF_pv) + (1.0f+k2c)*minmod(CellB_pv-CellF_pv, 2.0*(CellF_pv-CellFF_pv)));
		pv_R = CellB_pv - 0.25f*((1.0f-k2c)*minmod(2.0f*(CellB_pv-CellF_pv), CellBB_pv-CellB_pv) + (1.0f+k2c)*minmod(CellB_pv-CellF_pv, 2.0*(CellBB_pv-CellB_pv)));
	}						  
	if(reconstructType == ReconstructTypes::MUSCL2o)
	{
		r1 = (CellF_pv-CellFF_pv+epsilon) / (CellB_pv-CellF_pv+epsilon);
		r2 = (CellB_pv-CellF_pv+epsilon) / (CellBB_pv-CellB_pv+epsilon);
		f1 = 1.0f-eta + (eta-0.55f)/2.0f/r2 + (eta+0.55f)/2.0f*r1;
		f.dens = MAX(0.0f, MIN3(2.0f, f1.dens, 2.0f*r1.dens));
		f.uvel = MAX(0.0f, MIN3(2.0f, f1.uvel, 2.0f*r1.uvel));
		f.vvel = MAX(0.0f, MIN3(2.0f, f1.vvel, 2.0f*r1.vvel));
		f.wvel = MAX(0.0f, MIN3(2.0f, f1.wvel, 2.0f*r1.wvel));
		f.press = MAX(0.0f, MIN3(2.0f, f1.press, 2.0f*r1.press));
		pv_L = CellF_pv + 0.5f*f*(CellB_pv-CellF_pv);

		r1 = (CellBB_pv-CellB_pv+epsilon) / (CellB_pv-CellF_pv+epsilon);
		r2 = (CellB_pv-CellF_pv+epsilon) / (CellF_pv-CellFF_pv+epsilon);
		f1 = 1.0f-eta + (eta-0.55f)/2.0f/r2 + (eta+0.55f)/2.0f*r1;
		f.dens = MAX(0.0f, MIN3(2.0f, f1.dens, 2.0f*r1.dens));
		f.uvel = MAX(0.0f, MIN3(2.0f, f1.uvel, 2.0f*r1.uvel));
		f.vvel = MAX(0.0f, MIN3(2.0f, f1.vvel, 2.0f*r1.vvel));
		f.wvel = MAX(0.0f, MIN3(2.0f, f1.wvel, 2.0f*r1.wvel));
		f.press = MAX(0.0f, MIN3(2.0f, f1.press, 2.0f*r1.press));
		pv_R = CellB_pv - 0.5f*f*(CellB_pv-CellF_pv);
	}
	//if(reconstructType == ReconstructTypes::MUSCL2a)	// use Albada limiter
	//{
	//	eps = POW(Vol, 1.0/3.0);
	//	aL = CellB_cv-CellF_cv;		bL = CellF_cv-CellFF_cv;
	//	aR = CellBB_cv-CellB_cv;	bR = CellB_cv-CellF_cv;
	//	deltaL = (aL*(bL*bL+eps) + bL*(aL*aL+eps)) / (aL*aL + bL*bL + 2*eps);
	//	deltaR = (aR*(bR*bR+eps) + bR*(aR*aR+eps)) / (aR*aR + bR*bR + 2*eps);
	//	cv_L = CellF_cv + deltaL/2.0;
	//	cv_R = CellB_cv - deltaR/2.0;
	//}

	//if(reconstructType == ReconstructTypes::MUSCL2a2)	// use Albada limiter
	//{
	//	TConsVars rL, rR, psiL, psiR;

	//	rL = (CellB_cv-CellF_cv) / (CellF_cv-CellFF_cv);
	//	rR = (CellB_cv-CellF_cv) / (CellBB_cv-CellB_cv);
	//	if(rL.dens>0.0)
	//		psiL.dens = (rL.dens + rL.dens*rL.dens) / (1.0 + rL.dens*rL.dens);
	//	else
	//		psiL.dens = 0.0;
	//	if(rL.xmom>0.0)
	//		psiL.xmom = (rL.xmom + rL.xmom*rL.xmom) / (1.0 + rL.xmom*rL.xmom);
	//	else
	//		psiL.xmom = 0.0;
	//	if(rL.ymom>0.0)
	//		psiL.ymom = (rL.ymom + rL.ymom*rL.ymom) / (1.0 + rL.ymom*rL.ymom);
	//	else
	//		psiL.ymom = 0.0;
	//	if(rL.zmom>0.0)
	//		psiL.zmom = (rL.zmom + rL.zmom*rL.zmom) / (1.0 + rL.zmom*rL.zmom);
	//	else
	//		psiL.zmom = 0.0;
	//	if(rL.ener>0.0)
	//		psiL.ener = (rL.ener + rL.ener*rL.ener) / (1.0 + rL.ener*rL.ener);
	//	else
	//		psiL.ener = 0.0;

	//	if(rR.dens>0.0)
	//		psiR.dens = (rR.dens + rR.dens*rR.dens) / (1.0 + rR.dens*rR.dens);
	//	else
	//		psiR.dens = 0.0;
	//	if(rR.xmom>0.0)
	//		psiR.xmom = (rR.xmom + rR.xmom*rR.xmom) / (1.0 + rR.xmom*rR.xmom);
	//	else
	//		psiR.xmom = 0.0;
	//	if(rR.ymom>0.0)
	//		psiR.ymom = (rR.ymom + rR.ymom*rR.ymom) / (1.0 + rR.ymom*rR.ymom);
	//	else
	//		psiR.ymom = 0.0;
	//	if(rR.zmom>0.0)
	//		psiR.zmom = (rR.zmom + rR.zmom*rR.zmom) / (1.0 + rR.zmom*rR.zmom);
	//	else
	//		psiR.zmom = 0.0;
	//	if(rR.ener>0.0)
	//		psiR.ener = (rR.ener + rR.ener*rR.ener) / (1.0 + rR.ener*rR.ener);
	//	else
	//		psiR.ener = 0.0;

	//	cv_L = CellF_cv + psiL*(CellF_cv-CellFF_cv)/2.0;
	//	cv_R = CellB_cv - psiR*(CellBB_cv-CellB_cv)/2.0;
	//}

	//if(reconstructType == ReconstructTypes::MUSCL2mod)	
	//{
	//	TConsVars uint, rL, rR, psiL, psiR;
	//	
	//	rL = (CellB_cv-CellF_cv) / (CellF_cv-CellFF_cv);
	//	rR = (CellB_cv-CellF_cv) / (CellBB_cv-CellB_cv);

	//	uint.dens = 1.0;
	//	uint.xmom = 1.0;
	//	uint.ymom = 1.0;
	//	uint.zmom = 1.0;
	//	uint.ener = 1.0;
	//	psiL = minmod(uint, rL);
	//	psiR = minmod(uint, rR);

	//	cv_L = CellF_cv + psiL*(CellF_cv-CellFF_cv)/2.0;
	//	cv_R = CellB_cv - psiR*(CellBB_cv-CellB_cv)/2.0;
	//}

	//if(reconstructType == ReconstructTypes::MUSCL2r)	
	//{
	//	TConsVars rL, rR, psiL, psiR;

	//	rL = (CellB_cv-CellF_cv) / (CellF_cv-CellFF_cv);
	//	rR = (CellB_cv-CellF_cv) / (CellBB_cv-CellB_cv);
	//	psiL.dens = (rL.dens+ABS(rL.dens)) / (1.0+rL.dens);
	//	psiL.xmom = (rL.dens+ABS(rL.xmom)) / (1.0+rL.xmom);
	//	psiL.ymom = (rL.dens+ABS(rL.ymom)) / (1.0+rL.ymom);
	//	psiL.zmom = (rL.dens+ABS(rL.zmom)) / (1.0+rL.zmom);
	//	psiL.ener = (rL.dens+ABS(rL.ener)) / (1.0+rL.ener);
	//	psiR.dens = (rR.dens+ABS(rR.dens)) / (1.0+rR.dens);
	//	psiR.xmom = (rR.dens+ABS(rR.xmom)) / (1.0+rR.xmom);
	//	psiR.ymom = (rR.dens+ABS(rR.ymom)) / (1.0+rR.ymom);
	//	psiR.zmom = (rR.dens+ABS(rR.zmom)) / (1.0+rR.zmom);
	//	psiR.ener = (rR.dens+ABS(rR.ener)) / (1.0+rR.ener);
	//	cv_L = CellF_cv + psiL*(CellF_cv-CellFF_cv)/2.0;
	//	cv_R = CellB_cv - psiR*(CellBB_cv-CellB_cv)/2.0;
	//}
	//if(reconstructType == ReconstructTypes::MUSCL2max)	
	//{
	//	TConsVars rL, rR, psiL, psiR;

	//	rL = (CellB_cv-CellF_cv) / (CellF_cv-CellFF_cv);
	//	rR = (CellB_cv-CellF_cv) / (CellBB_cv-CellB_cv);
	//	psiL.dens = MAX3(0.0, MIN(1.0, 2.0*rL.dens), MIN(2.0, rL.dens));
	//	psiL.xmom = MAX3(0.0, MIN(1.0, 2.0*rL.xmom), MIN(2.0, rL.xmom));
	//	psiL.ymom = MAX3(0.0, MIN(1.0, 2.0*rL.ymom), MIN(2.0, rL.ymom));
	//	psiL.zmom = MAX3(0.0, MIN(1.0, 2.0*rL.zmom), MIN(2.0, rL.zmom));
	//	psiL.ener = MAX3(0.0, MIN(1.0, 2.0*rL.ener), MIN(2.0, rL.ener));
	//	psiR.dens = MAX3(0.0, MIN(1.0, 2.0*rR.dens), MIN(2.0, rR.dens));
	//	psiR.xmom = MAX3(0.0, MIN(1.0, 2.0*rR.xmom), MIN(2.0, rR.xmom));
	//	psiR.ymom = MAX3(0.0, MIN(1.0, 2.0*rR.ymom), MIN(2.0, rR.ymom));
	//	psiR.zmom = MAX3(0.0, MIN(1.0, 2.0*rR.zmom), MIN(2.0, rR.zmom));
	//	psiR.ener = MAX3(0.0, MIN(1.0, 2.0*rR.ener), MIN(2.0, rR.ener));
	//			
	//	cv_L = CellF_cv + psiL*(CellF_cv-CellFF_cv)/2.0;
	//	cv_R = CellB_cv - psiR*(CellBB_cv-CellB_cv)/2.0;
	//}
	if(reconstructType == ReconstructTypes::MUSCL3)
	{
		SL_limiter = (2.0*(CellF_pv - CellFF_pv)*(CellB_pv - CellF_pv) + epsilon) / (POW(CellF_pv-CellFF_pv, 2.0) + POW(CellB_pv-CellF_pv, 2.0) + epsilon);
		SR_limiter = (2.0*(CellBB_pv - CellB_pv)*(CellB_pv - CellF_pv) + epsilon) / (POW(CellBB_pv-CellB_pv, 2.0) + POW(CellB_pv-CellF_pv, 2.0) + epsilon);
		pv_L   = CellF_pv + SL_limiter/4.0*((1.0-SL_limiter/3.0)*(CellF_pv-CellFF_pv) + (1.0+SL_limiter/3.0)*(CellB_pv-CellF_pv));	
		pv_R   = CellB_pv - SR_limiter/4.0*((1.0-SR_limiter/3.0)*(CellBB_pv-CellB_pv) + (1.0+SR_limiter/3.0)*(CellB_pv-CellF_pv));
	}
	// 重构完成之后，要检查press 和 density 是否为负，如果是，则用一阶格式
	if(pv_L.dens<=lim_zero || pv_L.press<=lim_zero)	
		pv_L = CellF_pv;
	if(pv_R.dens<=lim_zero || pv_R.press<=lim_zero)	
		pv_R = CellB_pv;

	cv_L.dens = pv_L.dens;
	cv_L.xmom = pv_L.dens*pv_L.uvel;
	cv_L.ymom = pv_L.dens*pv_L.vvel;
	cv_L.zmom = pv_L.dens*pv_L.wvel;
	cv_L.ener = pv_L.press/(Gamma-1.0f) + pv_L.dens*(pv_L.uvel*pv_L.uvel+pv_L.vvel*pv_L.vvel+pv_L.wvel*pv_L.wvel)/2.0f;
		
	cv_R.dens = pv_R.dens;
	cv_R.xmom = pv_R.dens*pv_R.uvel;
	cv_R.ymom = pv_R.dens*pv_R.vvel;
	cv_R.zmom = pv_R.dens*pv_R.wvel;
	cv_R.ener = pv_R.press/(Gamma-1.0f) + pv_R.dens*(pv_R.uvel*pv_R.uvel+pv_R.vvel*pv_R.vvel+pv_R.wvel*pv_R.wvel)/2.0f;


}

//void CSpaceDiscr::CharacterReconstruct(TConsVars &CellFF_cv, TConsVars &CellF_cv, TConsVars &CellB_cv, TConsVars &CellBB_cv, TConsVars &cv_L, TConsVars &cv_R)
//{
//	TConsVars CellF_ch, CellB_ch, CellFF_ch, CellBB_ch;
//	TConsVars SL_limiter, SR_limiter;
//	TConsVars cv_F, ch_L, ch_R;
//	REAL S[6][6], S1[6][6]; 
//	REAL Gamma =  CFluidProps::Gamma;
//	REAL Rgas =   CFluidProps::Rgas;
//	REAL uvel, vvel, wvel, temp, press, V2, C2, V_bar, C;
//	REAL Med0, Med1, Med2, Med3, Med4;
//	REAL epsilon=1E-6f, eta=0.05f;		// Van Albada limiter's infinitesimal
//	REAL k2c=1.0f, k2u=-1.0f;
//	TConsVars r1, r2, f1, f;
//
//	cv_F = 0.5f*(CellF_cv + CellB_cv);
//	uvel = cv_F.xmom / cv_F.dens;
//	vvel = cv_F.ymom / cv_F.dens;
//	wvel = cv_F.zmom / cv_F.dens;
//	V2   = 0.5f*(uvel*uvel + vvel*vvel + wvel*wvel);
//	press= (Gamma-1.0f)*(CellF_cv.ener - cv_F.dens*V2);
//	temp = press/cv_F.dens/Rgas;
//	C = SQRT(Gamma*Rgas*temp);
//	C2 = C*C;
//	Med0 = C2/(Gamma-1.0f); 
//	Med1 = (Gamma-1.0f)/C*SQRT(3.0f);
//	Med2 = (Gamma-1.0f)/C2;
//	Med3 = 1.0f/Med0;
//	Med4 = 0.5f/C/SQRT(3.0f);
//	V_bar= (uvel+vvel+wvel);
//
//	S[1][1] = V2-Med0;			S[1][2] = -uvel;			   S[1][3] = -vvel;				 S[1][4] = -wvel;						S[1][5] = 1.0f;
//	S[2][1] = uvel-vvel;		S[2][2] = -1.0f;			   S[2][3] = 1.0f;				 S[2][4] = 0.0;							S[2][5] = 0.0;
//	S[3][1] = uvel-wvel;		S[3][2] = -1.0f;			   S[3][3] = 0.0;				 S[3][4] = 1.0f;						S[3][5] = 0.0;
//	S[4][1] = -V_bar-V2*Med1;	S[4][2] = 1.0f+Med1*uvel;	   S[4][3] = 1.0f+Med1*vvel;	 S[4][4] = 1.0f+Med1*wvel;				S[4][5] = -Med1;
//	S[5][1] = -V_bar+V2*Med1;	S[5][2] = 1.0f-Med1*uvel;	   S[5][3] = 1.0f-Med1*vvel;	 S[5][4] = 1.0f-Med1*wvel;				S[5][5] = Med1;
//
//	S1[1][1] = -Med3;			S1[1][2] = 0.0;				   S1[1][3] = 0.0;				 S1[1][4] = -Med4;						S1[1][5] = Med4;
//	S1[2][1] = -Med3*uvel;		S1[2][2] = -1.0f/3.0f;		   S1[2][3] = -1.0f/3.0f;		 S1[2][4] = 1.0f/6.0f-Med4*uvel;		S1[2][5] = 1.0f/6.0f+Med4*uvel;
//	S1[3][1] = -Med3*vvel;		S1[3][2] = 2.0f/3.0f;		   S1[3][3] = -1.0f/3.0f;		 S1[3][4] = 1.0f/6.0f-Med4*vvel;		S1[3][5] = 1.0f/6.0f+Med4*vvel;
//	S1[4][1] = -Med3*wvel;		S1[4][2] = -1.0f/3.0f;		   S1[4][3] = 2.0f/3.0f;		 S1[4][4] = 1.0f/6.0f-Med4*wvel;		S1[4][5] = 1.0f/6.0f+Med4*wvel;
//	S1[5][1] = -Med3*V2;;		S1[5][2] = vvel-V_bar/3.0f;	   S1[5][3] = wvel-V_bar/3.0f;	 S1[5][4] = V_bar/6.0f-Med4*(V2+Med0);	S1[5][5] = V_bar/6.0f+Med4*(V2+Med0);
//
//	CellFF_ch = S*CellFF_cv;
//	CellF_ch  = S*CellF_cv;
//	CellB_ch  = S*CellB_cv;
//	CellBB_ch = S*CellBB_cv;
//
//	if(reconstructType == ReconstructTypes::NND)
//	{
//		ch_L = CellF_ch + 0.5f*minmod(CellB_ch-CellF_ch, CellF_ch-CellFF_ch);
//		ch_R = CellB_ch - 0.5f*minmod(CellB_ch-CellF_ch, CellBB_ch-CellB_ch);
//	}
//	if(reconstructType == ReconstructTypes::MUSCL2u)
//	{
//		ch_L = CellF_ch + 0.25*((1.0f-k2u)*minmod(2.0f*(CellB_ch-CellF_ch), CellF_ch-CellFF_ch) + (1.0f+k2u)*minmod(2.0f*(CellF_ch-CellFF_ch), CellB_ch-CellF_ch));
//		ch_R = CellB_ch - 0.25*((1.0f-k2u)*minmod(2.0f*(CellB_ch-CellF_ch), CellBB_ch-CellB_ch) + (1.0f+k2u)*minmod(2.0f*(CellBB_ch-CellB_ch), CellB_ch-CellF_ch));
//	}
//	if(reconstructType == ReconstructTypes::MUSCL2c)
//	{							   
//		ch_L = CellF_ch + 0.25*((1.0f-k2c)*minmod(2.0f*(CellB_ch-CellF_ch), CellF_ch-CellFF_ch) + (1.0f+k2c)*minmod(CellB_ch-CellF_ch, 2.0f*(CellF_ch-CellFF_ch)));
//		ch_R = CellB_ch - 0.25*((1.0f-k2c)*minmod(2.0f*(CellB_ch-CellF_ch), CellBB_ch-CellB_ch) + (1.0f+k2c)*minmod(CellB_ch-CellF_ch, 2.0f*(CellBB_ch-CellB_ch)));
//	}
//	if(reconstructType == ReconstructTypes::MUSCL2o)
//	{
//		r1 = (CellF_ch-CellFF_ch+epsilon) / (CellB_ch-CellF_ch+epsilon);
//		r2 = (CellB_ch-CellF_ch+epsilon) / (CellBB_ch-CellB_ch+epsilon);
//		f1 = 1.0f-eta + (eta-0.55f)/2.0f/r2 + (eta+0.55f)/2.0f*r1;
//		f.dens = MAX(0.0f, MIN3(2.0f, f1.dens, 2.0f*r1.dens));
//		f.xmom = MAX(0.0f, MIN3(2.0f, f1.xmom, 2.0f*r1.xmom));
//		f.ymom = MAX(0.0f, MIN3(2.0f, f1.ymom, 2.0f*r1.ymom));
//		f.zmom = MAX(0.0f, MIN3(2.0f, f1.zmom, 2.0f*r1.zmom));
//		f.ener = MAX(0.0f, MIN3(2.0f, f1.ener, 2.0f*r1.ener));
//		ch_L = CellF_ch + 0.5f*f*(CellB_ch-CellF_ch);
//
//		r1 = (CellBB_ch-CellB_ch+epsilon) / (CellB_ch-CellF_ch+epsilon);
//		r2 = (CellB_ch-CellF_ch+epsilon) / (CellF_ch-CellFF_ch+epsilon);
//		f1 = 1.0f-eta + (eta-0.55f)/2.0f/r2 + (eta+0.55f)/2.0f*r1;
//		f.dens = MAX(0.0f, MIN3(2.0f, f1.dens, 2.0f*r1.dens));
//		f.xmom = MAX(0.0f, MIN3(2.0f, f1.xmom, 2.0f*r1.xmom));
//		f.ymom = MAX(0.0f, MIN3(2.0f, f1.ymom, 2.0f*r1.ymom));
//		f.zmom = MAX(0.0f, MIN3(2.0f, f1.zmom, 2.0f*r1.zmom));
//		f.ener = MAX(0.0f, MIN3(2.0f, f1.ener, 2.0f*r1.ener));
//		ch_R = CellB_ch - 0.5f*f*(CellB_ch-CellF_ch);
//	}
//	if(reconstructType == ReconstructTypes::MUSCL3)
//	{
//		SL_limiter = (2.0f*(CellF_ch - CellFF_ch)*(CellB_ch - CellF_ch) + epsilon) / (POW(CellF_ch-CellFF_ch, 2.0f) + POW(CellB_ch-CellF_ch, 2.0f) + epsilon);
//		SR_limiter = (2.0f*(CellBB_ch - CellB_ch)*(CellB_ch - CellF_ch) + epsilon) / (POW(CellBB_ch-CellB_ch, 2.0f) + POW(CellB_ch-CellF_ch, 2.0f) + epsilon);
//		ch_L   = CellF_ch + SL_limiter/4.0f*((1.0f-SL_limiter/3.0f)*(CellF_ch-CellFF_ch) + (1.0f+SL_limiter/3.0f)*(CellB_ch-CellF_ch));	
//		ch_R   = CellB_ch - SR_limiter/4.0f*((1.0f-SR_limiter/3.0f)*(CellBB_ch-CellB_ch) + (1.0f+SR_limiter/3.0f)*(CellB_ch-CellF_ch));
//	}
//
//	cv_L   = S1*ch_L; 
//	cv_R   = S1*ch_R; 
//}

void CSpaceDiscr::CharacterReconstruct(TPrimVars &CellF_pv, TPrimVars &CellB_pv, TConsVars &CellFF_cv, TConsVars &CellF_cv, TConsVars &CellB_cv, TConsVars &CellBB_cv, TConsVars &cv_L, TConsVars &cv_R)
{
	TConsVars CellF_ch, CellB_ch, CellFF_ch, CellBB_ch;
	TConsVars SL_limiter, SR_limiter;
	TConsVars ch_L, ch_R;
	TPrimVars pv_F;
	REAL S[6][6], S1[6][6]; 
	REAL Gamma =  CFluidProps::Gamma;
	REAL Rgas =   CFluidProps::Rgas;
	REAL uvel, vvel, wvel, temp, press, V2, C2, V_bar, C;
	REAL Med0, Med1, Med2, Med3, Med4;
	REAL epsilon=1E-6f, eta=0.05f;		// Van Albada limiter's infinitesimal
	REAL k2c=1.0f, k2u=-1.0f;
	TConsVars r1, r2, f1, f;
	REAL lim_zero=1E-20f;
	REAL pressL, pressR;

	pv_F = 0.5f*(CellF_pv + CellB_pv);
	uvel = pv_F.uvel;
	vvel = pv_F.vvel;
	wvel = pv_F.wvel;
	press= pv_F.press;
	V2   = 0.5f*(uvel*uvel + vvel*vvel + wvel*wvel);
	temp = press/pv_F.dens/Rgas;
	C = SQRT(Gamma*Rgas*temp);
	C2 = C*C;
	Med0 = C2/(Gamma-1.0f); 
	Med1 = (Gamma-1.0f)/C*SQRT(3.0f);
	Med2 = (Gamma-1.0f)/C2;
	Med3 = 1.0f/Med0;
	Med4 = 0.5f/C/SQRT(3.0f);
	V_bar= (uvel+vvel+wvel);

	S[1][1] = V2-Med0;			S[1][2] = -uvel;			   S[1][3] = -vvel;				 S[1][4] = -wvel;						S[1][5] = 1.0f;
	S[2][1] = uvel-vvel;		S[2][2] = -1.0f;			   S[2][3] = 1.0f;				 S[2][4] = 0.0;							S[2][5] = 0.0;
	S[3][1] = uvel-wvel;		S[3][2] = -1.0f;			   S[3][3] = 0.0;				 S[3][4] = 1.0f;						S[3][5] = 0.0;
	S[4][1] = -V_bar-V2*Med1;	S[4][2] = 1.0f+Med1*uvel;	   S[4][3] = 1.0f+Med1*vvel;	 S[4][4] = 1.0f+Med1*wvel;				S[4][5] = -Med1;
	S[5][1] = -V_bar+V2*Med1;	S[5][2] = 1.0f-Med1*uvel;	   S[5][3] = 1.0f-Med1*vvel;	 S[5][4] = 1.0f-Med1*wvel;				S[5][5] = Med1;

	S1[1][1] = -Med3;			S1[1][2] = 0.0;				   S1[1][3] = 0.0;				 S1[1][4] = -Med4;						S1[1][5] = Med4;
	S1[2][1] = -Med3*uvel;		S1[2][2] = -1.0f/3.0f;		   S1[2][3] = -1.0f/3.0f;		 S1[2][4] = 1.0f/6.0f-Med4*uvel;		S1[2][5] = 1.0f/6.0f+Med4*uvel;
	S1[3][1] = -Med3*vvel;		S1[3][2] = 2.0f/3.0f;		   S1[3][3] = -1.0f/3.0f;		 S1[3][4] = 1.0f/6.0f-Med4*vvel;		S1[3][5] = 1.0f/6.0f+Med4*vvel;
	S1[4][1] = -Med3*wvel;		S1[4][2] = -1.0f/3.0f;		   S1[4][3] = 2.0f/3.0f;		 S1[4][4] = 1.0f/6.0f-Med4*wvel;		S1[4][5] = 1.0f/6.0f+Med4*wvel;
	S1[5][1] = -Med3*V2;;		S1[5][2] = vvel-V_bar/3.0f;	   S1[5][3] = wvel-V_bar/3.0f;	 S1[5][4] = V_bar/6.0f-Med4*(V2+Med0);	S1[5][5] = V_bar/6.0f+Med4*(V2+Med0);

	CellFF_ch = S*CellFF_cv;
	CellF_ch  = S*CellF_cv;
	CellB_ch  = S*CellB_cv;
	CellBB_ch = S*CellBB_cv;

	if(reconstructType == ReconstructTypes::NND)
	{
		ch_L = CellF_ch + 0.5f*minmod(CellB_ch-CellF_ch, CellF_ch-CellFF_ch);
		ch_R = CellB_ch - 0.5f*minmod(CellB_ch-CellF_ch, CellBB_ch-CellB_ch);
	}
	if(reconstructType == ReconstructTypes::MUSCL2u)
	{
		ch_L = CellF_ch + 0.25*((1.0f-k2u)*minmod(2.0f*(CellB_ch-CellF_ch), CellF_ch-CellFF_ch) + (1.0f+k2u)*minmod(2.0f*(CellF_ch-CellFF_ch), CellB_ch-CellF_ch));
		ch_R = CellB_ch - 0.25*((1.0f-k2u)*minmod(2.0f*(CellB_ch-CellF_ch), CellBB_ch-CellB_ch) + (1.0f+k2u)*minmod(2.0f*(CellBB_ch-CellB_ch), CellB_ch-CellF_ch));
	}
	if(reconstructType == ReconstructTypes::MUSCL2c)
	{							   
		ch_L = CellF_ch + 0.25*((1.0f-k2c)*minmod(2.0f*(CellB_ch-CellF_ch), CellF_ch-CellFF_ch) + (1.0f+k2c)*minmod(CellB_ch-CellF_ch, 2.0f*(CellF_ch-CellFF_ch)));
		ch_R = CellB_ch - 0.25*((1.0f-k2c)*minmod(2.0f*(CellB_ch-CellF_ch), CellBB_ch-CellB_ch) + (1.0f+k2c)*minmod(CellB_ch-CellF_ch, 2.0f*(CellBB_ch-CellB_ch)));
	}
	if(reconstructType == ReconstructTypes::MUSCL2o)
	{
		r1 = (CellF_ch-CellFF_ch+epsilon) / (CellB_ch-CellF_ch+epsilon);
		r2 = (CellB_ch-CellF_ch+epsilon) / (CellBB_ch-CellB_ch+epsilon);
		f1 = 1.0f-eta + (eta-0.55f)/2.0f/r2 + (eta+0.55f)/2.0f*r1;
		f.dens = MAX(0.0f, MIN3(2.0f, f1.dens, 2.0f*r1.dens));
		f.xmom = MAX(0.0f, MIN3(2.0f, f1.xmom, 2.0f*r1.xmom));
		f.ymom = MAX(0.0f, MIN3(2.0f, f1.ymom, 2.0f*r1.ymom));
		f.zmom = MAX(0.0f, MIN3(2.0f, f1.zmom, 2.0f*r1.zmom));
		f.ener = MAX(0.0f, MIN3(2.0f, f1.ener, 2.0f*r1.ener));
		ch_L = CellF_ch + 0.5f*f*(CellB_ch-CellF_ch);

		r1 = (CellBB_ch-CellB_ch+epsilon) / (CellB_ch-CellF_ch+epsilon);
		r2 = (CellB_ch-CellF_ch+epsilon) / (CellF_ch-CellFF_ch+epsilon);
		f1 = 1.0f-eta + (eta-0.55f)/2.0f/r2 + (eta+0.55f)/2.0f*r1;
		f.dens = MAX(0.0f, MIN3(2.0f, f1.dens, 2.0f*r1.dens));
		f.xmom = MAX(0.0f, MIN3(2.0f, f1.xmom, 2.0f*r1.xmom));
		f.ymom = MAX(0.0f, MIN3(2.0f, f1.ymom, 2.0f*r1.ymom));
		f.zmom = MAX(0.0f, MIN3(2.0f, f1.zmom, 2.0f*r1.zmom));
		f.ener = MAX(0.0f, MIN3(2.0f, f1.ener, 2.0f*r1.ener));
		ch_R = CellB_ch - 0.5f*f*(CellB_ch-CellF_ch);
	}
	if(reconstructType == ReconstructTypes::MUSCL3)
	{
		SL_limiter = (2.0f*(CellF_ch - CellFF_ch)*(CellB_ch - CellF_ch) + epsilon) / (POW(CellF_ch-CellFF_ch, 2.0f) + POW(CellB_ch-CellF_ch, 2.0f) + epsilon);
		SR_limiter = (2.0f*(CellBB_ch - CellB_ch)*(CellB_ch - CellF_ch) + epsilon) / (POW(CellBB_ch-CellB_ch, 2.0f) + POW(CellB_ch-CellF_ch, 2.0f) + epsilon);
		ch_L   = CellF_ch + SL_limiter/4.0f*((1.0f-SL_limiter/3.0f)*(CellF_ch-CellFF_ch) + (1.0f+SL_limiter/3.0f)*(CellB_ch-CellF_ch));	
		ch_R   = CellB_ch - SR_limiter/4.0f*((1.0f-SR_limiter/3.0f)*(CellBB_ch-CellB_ch) + (1.0f+SR_limiter/3.0f)*(CellB_ch-CellF_ch));
	}

	cv_L   = S1*ch_L; 
	cv_R   = S1*ch_R; 
		
	// 重构完成之后，要检查press 和 density 是否为负，如果是，则用一阶格式
	pressL = (Gamma-1.0f)*(cv_L.ener - 0.5f*(cv_L.xmom*cv_L.xmom+cv_L.ymom*cv_L.ymom+cv_L.zmom*cv_L.zmom)/cv_L.dens);
	pressR = (Gamma-1.0f)*(cv_R.ener - 0.5f*(cv_R.xmom*cv_R.xmom+cv_R.ymom*cv_R.ymom+cv_R.zmom*cv_R.zmom)/cv_R.dens);
	if(cv_L.dens<=lim_zero || pressL<=lim_zero)	
		cv_L = CellF_cv;
	if(cv_R.dens<=lim_zero || pressR<=lim_zero)	
		cv_R = CellB_cv;
}


// minmod() function, used for NND scheme
TConsVars CSpaceDiscr::minmod(const TConsVars &A, const TConsVars &B)	//it is meanningless to set minmod to be inline, because there is "if" clause
{
	TConsVars C;
	if (A.dens*B.dens < 0.0)
		C.dens = 0.0;
	else if (ABS(A.dens) <= ABS(B.dens))
		C.dens = A.dens;
	else
		C.dens = B.dens;

	if (A.xmom*B.xmom < 0.0)
		C.xmom = 0.0;
	else if (ABS(A.xmom) <= ABS(B.xmom))
		C.xmom = A.xmom;
	else
		C.xmom = B.xmom;

	if (A.ymom*B.ymom < 0.0)
		C.ymom = 0.0;
	else if (ABS(A.ymom) <= ABS(B.ymom))
		C.ymom = A.ymom;
	else
		C.ymom = B.ymom;

	if (A.zmom*B.zmom < 0.0)
		C.zmom = 0.0;
	else if (ABS(A.zmom) <= ABS(B.zmom))
		C.zmom = A.zmom;
	else
		C.zmom = B.zmom;

	if (A.ener*B.ener < 0.0)
		C.ener = 0.0;
	else if (ABS(A.ener) <= ABS(B.ener))
		C.ener = A.ener;
	else
		C.ener = B.ener;

	return C;
}

// override minmod() function with TPrimVars, used for NND scheme
TPrimVars CSpaceDiscr::minmod(const TPrimVars &A, const TPrimVars &B)	//it is meanningless to set minmod to be inline, because there is "if" clause
{
	TPrimVars C;
	if (A.dens*B.dens < 0.0)
		C.dens = 0.0;
	else if (ABS(A.dens) <= ABS(B.dens))
		C.dens = A.dens;
	else
		C.dens = B.dens;

	if (A.uvel*B.uvel < 0.0)
		C.uvel = 0.0;
	else if (ABS(A.uvel) <= ABS(B.uvel))
		C.uvel = A.uvel;
	else
		C.uvel = B.uvel;

	if (A.vvel*B.vvel < 0.0)
		C.vvel = 0.0;
	else if (ABS(A.vvel) <= ABS(B.vvel))
		C.vvel = A.vvel;
	else
		C.vvel = B.vvel;

	if (A.wvel*B.wvel < 0.0)
		C.wvel = 0.0;
	else if (ABS(A.wvel) <= ABS(B.wvel))
		C.wvel = A.wvel;
	else
		C.wvel = B.wvel;

	if (A.temp*B.temp < 0.0)
		C.temp = 0.0;
	else if (ABS(A.temp) <= ABS(B.temp))
		C.temp = A.temp;
	else
		C.temp = B.temp;

	if (A.press*B.press < 0.0)
		C.press = 0.0;
	else if (ABS(A.press) <= ABS(B.press))
		C.press = A.press;
	else
		C.press = B.press;

	return C;
}