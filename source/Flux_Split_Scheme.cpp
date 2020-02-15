
#include "defs.h"
#include "SpaceDiscr.h"


void CSpaceDiscr::L_F_Scheme(int i, int j, int k, int flag, const CGeometry &Geometry, const CFluidProps &FluidProp, TNode &norm, TNode &tau1, TNode &tau2, TConsVars &cv_L, TConsVars &cv_R, TConsVars &Fc)
{
	TConsVars Q_L, Q_R;
	REAL densL, uvelL, vvelL, wvelL, enerL, pL, csounL, CL2;
	REAL densR, uvelR, vvelR, wvelR, enerR, pR, csounR, CR2;
	REAL lambda_L, lambda_R;
	REAL Gamma = CFluidProps::Gamma;
	REAL Q_uvel_e = 0.0f;

	TConsVars Diss_L, Diss_R, Fc_L, Fc_R, Fc_plus, Fc_minus; 

	Coord_Rotate(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Q_L, Q_R, Q_uvel_e);

	densL = Q_L.dens;
	uvelL = Q_L.xmom / Q_L.dens;
	vvelL = Q_L.ymom / Q_L.dens;
	wvelL = Q_L.zmom / Q_L.dens;
	enerL = Q_L.ener / Q_L.dens;
	pL	  = (Gamma-1.0f)*densL*(enerL - 0.5f*(uvelL*uvelL + vvelL*vvelL + wvelL*wvelL));
	CL2	  = Gamma*(Gamma-1.0f)*(enerL - 0.5f*(uvelL*uvelL + vvelL*vvelL + wvelL*wvelL));
	if(CL2<0)
	{
		CL2 = 1E-5f;
	}
	csounL = SQRT(CL2);
	lambda_L = ABS(uvelL-Q_uvel_e) + csounL;

	densR = Q_R.dens;
	uvelR = Q_R.xmom / Q_R.dens;
	vvelR = Q_R.ymom / Q_R.dens;
	wvelR = Q_R.zmom / Q_R.dens;
	enerR = Q_R.ener / Q_R.dens;
	pR	  = (Gamma-1.0f)*densR*(enerR - 0.5f*(uvelR*uvelR + vvelR*vvelR + wvelR*wvelR));
	CR2	  = Gamma*(Gamma-1.0f)*(enerR - 0.5f*(uvelR*uvelR + vvelR*vvelR + wvelR*wvelR));
	if(CR2<0)
	{
		CR2 = 1E-5f;
	}
	csounR = SQRT(CR2);
	lambda_R = ABS(uvelR-Q_uvel_e) + csounR;

	// dissipation
	Diss_L = lambda_L * Q_L;
	Diss_R = lambda_R * Q_R;

	// convection
	Fc_L.dens = Q_L.dens*uvelL;
	Fc_L.xmom = Q_L.xmom*uvelL + pL;
	Fc_L.ymom = Q_L.ymom*uvelL;
	Fc_L.zmom = Q_L.zmom*uvelL;
	Fc_L.ener = (Q_L.ener+pL) * uvelL;

	Fc_R.dens = Q_R.dens*uvelR;
	Fc_R.xmom = Q_R.xmom*uvelR + pR;
	Fc_R.ymom = Q_R.ymom*uvelR;
	Fc_R.zmom = Q_R.zmom*uvelR;
	Fc_R.ener = (Q_R.ener+pR) * uvelR;

	// plus and minus flux
	Fc_plus.dens = (Fc_L.dens + Diss_L.dens)/2.0f;
	Fc_plus.xmom = (Fc_L.xmom + Diss_L.xmom)/2.0f;
	Fc_plus.ymom = (Fc_L.ymom + Diss_L.ymom)/2.0f;
	Fc_plus.zmom = (Fc_L.zmom + Diss_L.zmom)/2.0f;
	Fc_plus.ener = (Fc_L.ener + Diss_L.ener)/2.0f;
							    
	Fc_minus.dens = (Fc_R.dens - Diss_R.dens)/2.0f;
	Fc_minus.xmom = (Fc_R.xmom - Diss_R.xmom)/2.0f;
	Fc_minus.ymom = (Fc_R.ymom - Diss_R.ymom)/2.0f;
	Fc_minus.zmom = (Fc_R.zmom - Diss_R.zmom)/2.0f;
	Fc_minus.ener = (Fc_R.ener - Diss_R.ener)/2.0f;

	// flux
	Fc = Fc_plus + Fc_minus;
}


void CSpaceDiscr::StegerWarimgScheme(int i, int j, int k, int flag, const CGeometry &Geometry, const CFluidProps &FluidProp, TNode &norm, TNode &tau1, TNode &tau2, TConsVars &cv_L, TConsVars &cv_R, TConsVars &Fc)
{
	TConsVars Q_L, Q_R;
	REAL fplus[6], fminus[6], lambda_plus[6], lambda_minus[6];
	REAL A, B, W2, epsilon=1E-6f;
	REAL Gamma = CFluidProps::Gamma;
	REAL densL, uvelL, vvelL, wvelL, enerL, csounL, tmppL;
	REAL densR, uvelR, vvelR, wvelR, enerR, csounR, tmppR;
	REAL Q_uvel_e = 0.0f;

	Coord_Rotate(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Q_L, Q_R, Q_uvel_e);

	densL = Q_L.dens;
	uvelL = Q_L.xmom / Q_L.dens;
	vvelL = Q_L.ymom / Q_L.dens;
	wvelL = Q_L.zmom / Q_L.dens;
	enerL = Q_L.ener / Q_L.dens;
	tmppL = Gamma*(Gamma-1.0f)*(enerL - 0.5f*(uvelL*uvelL + vvelL*vvelL + wvelL*wvelL));
	tmppL = SQRT(tmppL*tmppL);
	csounL = SQRT(tmppL);

	lambda_plus[1] = 0.5f*((uvelL-Q_uvel_e) + SQRT((uvelL-Q_uvel_e)*(uvelL-Q_uvel_e) + epsilon*epsilon));
	lambda_plus[2] = 0.5f*((uvelL-Q_uvel_e) + SQRT((uvelL-Q_uvel_e)*(uvelL-Q_uvel_e) + epsilon*epsilon));
	lambda_plus[3] = 0.5f*((uvelL-Q_uvel_e) + SQRT((uvelL-Q_uvel_e)*(uvelL-Q_uvel_e) + epsilon*epsilon));
	lambda_plus[4] = 0.5f*(((uvelL-Q_uvel_e)-csounL) + SQRT(((uvelL-Q_uvel_e)-csounL)*((uvelL-Q_uvel_e)-csounL)+epsilon*epsilon));
	lambda_plus[5] = 0.5f*(((uvelL-Q_uvel_e)+csounL) + SQRT(((uvelL-Q_uvel_e)+csounL)*((uvelL-Q_uvel_e)+csounL)+epsilon*epsilon));

	W2 = (3.0f-Gamma)*(lambda_plus[4]+lambda_plus[5])*csounL*csounL/2.0f/(Gamma-1.0f);
	A = (uvelL-csounL)*(uvelL-csounL) + vvelL*vvelL +wvelL*wvelL;
	B = (uvelL+csounL)*(uvelL+csounL) + vvelL*vvelL+ wvelL*wvelL;
	fplus[1] = densL/2.0f/Gamma * (2.0f*(Gamma-1.0f)*lambda_plus[1]    + lambda_plus[4] + lambda_plus[5]);
	fplus[2] = densL/2.0f/Gamma * (2.0f*(Gamma-1.0f)*lambda_plus[1]*uvelL + lambda_plus[4]*(uvelL-csounL) + lambda_plus[5]*(uvelL+csounL));
	fplus[3] = densL/2.0f/Gamma * (2.0f*(Gamma-1.0f)*lambda_plus[1]*vvelL + lambda_plus[4]*vvelL + lambda_plus[5]*vvelL);
	fplus[4] = densL/2.0f/Gamma * (2.0f*(Gamma-1.0f)*lambda_plus[1]*wvelL + lambda_plus[4]*wvelL + lambda_plus[5]*wvelL);
	fplus[5] = densL/2.0f/Gamma * ((Gamma-1.0f)*lambda_plus[1]*(uvelL*uvelL+vvelL*vvelL +wvelL*wvelL) + 0.5f*lambda_plus[4]*A + 0.5f*lambda_plus[5]*B + W2);
	//----------------------------------------------------------

	densR = Q_R.dens;
	uvelR = Q_R.xmom / Q_R.dens;
	vvelR = Q_R.ymom / Q_R.dens;
	wvelR = Q_R.zmom / Q_R.dens;
	enerR = Q_R.ener / Q_R.dens;
	tmppR = Gamma*(Gamma-1.0f)*(enerR - 0.5f*(uvelR*uvelR + vvelR*vvelR + wvelR*wvelR));
	tmppR = SQRT(tmppR*tmppR);
	csounR = SQRT(tmppR);


	lambda_minus[1] = 0.5f*((uvelR-Q_uvel_e) - SQRT((uvelR-Q_uvel_e)*(uvelR-Q_uvel_e) + epsilon*epsilon));
	lambda_minus[2] = 0.5f*((uvelR-Q_uvel_e) - SQRT((uvelR-Q_uvel_e)*(uvelR-Q_uvel_e) + epsilon*epsilon));
	lambda_minus[3] = 0.5f*((uvelR-Q_uvel_e) - SQRT((uvelR-Q_uvel_e)*(uvelR-Q_uvel_e) + epsilon*epsilon));
	lambda_minus[4] = 0.5f*(((uvelR-Q_uvel_e)-csounR) - SQRT(((uvelR-Q_uvel_e)-csounR)*((uvelR-Q_uvel_e)-csounR)+epsilon*epsilon));
	lambda_minus[5] = 0.5f*(((uvelR-Q_uvel_e)+csounR) - SQRT(((uvelR-Q_uvel_e)+csounR)*((uvelR-Q_uvel_e)+csounR)+epsilon*epsilon));

	W2 = (3.0f-Gamma)*(lambda_minus[4]+lambda_minus[5])*csounR*csounR/2.0f/(Gamma-1.0f);
	A = (uvelR-csounR)*(uvelR-csounR) + vvelR*vvelR +wvelR*wvelR;
	B = (uvelR+csounR)*(uvelR+csounR) + vvelR*vvelR+ wvelR*wvelR;
	fminus[1] = densR/2.0f/Gamma * (2.0f*(Gamma-1.0f)*lambda_minus[1]    + lambda_minus[4] + lambda_minus[5]);
	fminus[2] = densR/2.0f/Gamma * (2.0f*(Gamma-1.0f)*lambda_minus[1]*uvelR + lambda_minus[4]*(uvelR-csounR) + lambda_minus[5]*(uvelR+csounR));
	fminus[3] = densR/2.0f/Gamma * (2.0f*(Gamma-1.0f)*lambda_minus[1]*vvelR + lambda_minus[4]*vvelR + lambda_minus[5]*vvelR);
	fminus[4] = densR/2.0f/Gamma * (2.0f*(Gamma-1.0f)*lambda_minus[1]*wvelR + lambda_minus[4]*wvelR + lambda_minus[5]*wvelR);
	fminus[5] = densR/2.0f/Gamma * ((Gamma-1.0f)*lambda_minus[1]*(uvelR*uvelR+vvelR*vvelR +wvelR*wvelR) + 0.5f*lambda_minus[4]*A + 0.5f*lambda_minus[5]*B + W2);

	Fc.dens = fplus[1] + fminus[1];
	Fc.xmom = fplus[2] + fminus[2];
	Fc.ymom = fplus[3] + fminus[3];
	Fc.zmom = fplus[4] + fminus[4];
	Fc.ener = fplus[5] + fminus[5];
}


void CSpaceDiscr::Van_Leer_Scheme(int i, int j, int k, int flag, const CGeometry &Geometry, const CFluidProps &FluidProp, TNode &norm, TNode &tau1, TNode &tau2, TConsVars &cv_L, TConsVars &cv_R, TConsVars &Fc)
{
	TConsVars Q_L, Q_R;
	REAL Gamma = CFluidProps::Gamma;
	REAL densL, uvelL, vvelL, wvelL, enerL, csounL, pL, MaL;
	REAL densR, uvelR, vvelR, wvelR, enerR, csounR, pR, MaR;
	REAL tmp, Map, Mam;
	REAL fp[6], fm[6];
	REAL Q_uvel_e = 0.0f;

	Coord_Rotate(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Q_L, Q_R, Q_uvel_e);

	densL = Q_L.dens;
	uvelL = Q_L.xmom / Q_L.dens;
	vvelL = Q_L.ymom / Q_L.dens;
	wvelL = Q_L.zmom / Q_L.dens;
	enerL = Q_L.ener / Q_L.dens;
	pL	  = (Gamma-1.0f)*densL*(enerL - 0.5f*(uvelL*uvelL + vvelL*vvelL + wvelL*wvelL));

	densR = Q_R.dens;
	uvelR = Q_R.xmom / Q_R.dens;
	vvelR = Q_R.ymom / Q_R.dens;
	wvelR = Q_R.zmom / Q_R.dens;
	enerR = Q_R.ener / Q_R.dens;
	pR	  = (Gamma-1.0f)*densR*(enerR - 0.5f*(uvelR*uvelR + vvelR*vvelR + wvelR*wvelR));
	
	csounL = SQRT(Gamma*pL/densL);
	csounR = SQRT(Gamma*pR/densR);
	MaL = uvelL/csounL;
	MaR = uvelR/csounR;

	if(MaL>1.0f)
	{
		fp[1] = densL*uvelL;
		fp[2] = densL*uvelL*uvelL + pL;
		fp[3] = densL*uvelL*vvelL;
		fp[4] = densL*uvelL*wvelL;
		fp[5] = uvelL*(Gamma*pL/(Gamma-1.0f) + 0.5f*densL*(uvelL*uvelL+vvelL*vvelL+wvelL*wvelL));
	}
	else if(ABS(MaL)<1.0f)	
	{
		Map = 0.25f*(1.0f+MaL)*(1.0f+MaL);
		tmp = densL*csounL*Map;
		fp[1] = tmp;
		fp[2] = tmp*((Gamma-1.0f)*uvelL+2.0f*csounL)/Gamma;
		fp[3] = tmp*vvelL;
		fp[4] = tmp*wvelL;
		fp[5] = tmp*(POW((Gamma-1.0f)*uvelL+2.0f*csounL, 2.0f)*0.5f/(Gamma*Gamma-1.0f) + 0.5f*(uvelL*uvelL+vvelL*vvelL+wvelL*wvelL));
	}
	else
	{
		fp[1] = 0.0f;
		fp[2] = 0.0f;
		fp[3] = 0.0f;
		fp[4] = 0.0f;
		fp[5] = 0.0f;
	}

	if(MaR>1.0f)
	{
		fm[1] = 0.0f;
		fm[2] = 0.0f;
		fm[3] = 0.0f;
		fm[4] = 0.0f;
		fm[5] = 0.0f;
	}
	else if(ABS(MaR)<1.0f)	
	{
		Mam = -0.25f*(1.0f-MaR)*(1.0f-MaR);
		tmp = densR*csounR*Mam;
		fm[1] = tmp;
		fm[2] = tmp*((Gamma-1.0f)*uvelR-2.0f*csounR)/Gamma;
		fm[3] = tmp*vvelR;
		fm[4] = tmp*wvelR;
		fm[5] = tmp*(POW((Gamma-1.0f)*uvelR-2.0f*csounR, 2.0f)*0.5f/(Gamma*Gamma-1.0f) + 0.5f*(uvelR*uvelR+vvelR*vvelR+wvelR*wvelR));
	}
	else
	{
		fm[1] = densR*uvelR;
		fm[2] = densR*uvelR*uvelR + pR;
		fm[3] = densR*uvelR*vvelR;
		fm[4] = densR*uvelR*wvelR;
		fm[5] = uvelR*(Gamma*pR/(Gamma-1.0f) + 0.5f*densR*(uvelR*uvelR+vvelR*vvelR+wvelR*wvelR));
	}

	 Fc.dens = fp[1]+fm[1];
	 Fc.xmom = fp[2]+fm[2];
	 Fc.ymom = fp[3]+fm[3];
	 Fc.zmom = fp[4]+fm[4];
	 Fc.ener = fp[5]+fm[5];
}


void CSpaceDiscr::AUSM_Plus_Scheme(int i, int j, int k, int flag, const CGeometry &Geometry, const CFluidProps &FluidProp, TNode &norm, TNode &tau1, TNode &tau2, TConsVars &cv_L, TConsVars &cv_R, TConsVars &Fc)
{
	TConsVars Q_L, Q_R;
	REAL Gamma = CFluidProps::Gamma;
	REAL densL, uvelL, vvelL, wvelL, enerL, csounL, pL, HL, MaL, fL, ppL;
	REAL densR, uvelR, vvelR, wvelR, enerR, csounR, pR, HR, MaR, fR, ppR;
	REAL csoun, Ma, Ma_abs, Ma_plus, Ma_minus, p_plus, p_minus, p;
	REAL fp[6], fm[6];
	REAL Q_uvel_e = 0.0f, delta=0.2f;

	Coord_Rotate(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Q_L, Q_R, Q_uvel_e);

	densL = Q_L.dens;
	uvelL = Q_L.xmom / Q_L.dens;
	vvelL = Q_L.ymom / Q_L.dens;
	wvelL = Q_L.zmom / Q_L.dens;
	enerL = Q_L.ener / Q_L.dens;
	pL	  = (Gamma-1.0f)*densL*(enerL - 0.5f*(uvelL*uvelL + vvelL*vvelL + wvelL*wvelL));
	HL    = enerL + pL/densL;

	densR = Q_R.dens;
	uvelR = Q_R.xmom / Q_R.dens;
	vvelR = Q_R.ymom / Q_R.dens;
	wvelR = Q_R.zmom / Q_R.dens;
	enerR = Q_R.ener / Q_R.dens;
	pR	  = (Gamma-1.0f)*densR*(enerR - 0.5f*(uvelR*uvelR + vvelR*vvelR + wvelR*wvelR));
	HR    = enerR + pR/densR;

	csounL = SQRT(Gamma*pL/densL);
	csounR = SQRT(Gamma*pR/densR);
	csoun = 0.5f*(csounL+csounR);

	MaL = (uvelL-Q_uvel_e)/csoun;
	MaR = (uvelR-Q_uvel_e)/csoun;
		
	if((MIN(pL/pR, pR/pL)>=0.75) && (MIN(pL/pR, pR/pL)<1.0))
	{
		ppL = 4.0f*MIN(pL/pR, pR/pL) - 3.0f;
		ppR = 4.0f*MIN(pL/pR, pR/pL) - 3.0f;
	}
	else
	{
		ppL = 0.0f;
		ppR = 0.0f;
	}
		
	if(ABS(MaL)>=1.0f)
		p_plus = 0.5f*(ABS(MaL)+MaL)/MaL;
	else
		//p_plus = 0.25f*(MaL+1.0f)*(MaL+1.0f)*(2.0f-MaL) + MaL*(MaL*MaL-1.0f)*(MaL*MaL-1.0f)*3.0f/16.0f;
		p_plus = 0.5f*(MaL+1.0f);		//用此公式更好

	if(ABS(MaR)>=1.0f)
		p_minus = 0.5f*(MaR-ABS(MaR))/MaR;
	else
		//p_minus = 0.25f*(MaR-1.0f)*(MaR-1.0f)*(2.0f+MaR) - MaR*(MaR*MaR-1.0f)*(MaR*MaR-1.0f)*3.0f/16.0f;
		p_minus = 0.5f*(1.0f-MaR);
	
	p = pL*p_plus + pR*p_minus;

	if(ABS(MaL)>=1.0f)
	{
		Ma_plus = 0.5f*(MaL+ABS(MaL));
		fL = 0.0f;
	}
	else
	{
		Ma_plus = 0.25f*(MaL+1.0f)*(MaL+1.0f) + 0.125f*(MaL*MaL-1.0f)*(MaL*MaL-1.0f);
		fL = (pL/p-1.0f)*ppL*ABS(0.25f*(MaL+1.0f)*(MaL+1.0f))*MIN(1.0f, POW(SQRT(uvelL*uvelL+vvelL*vvelL+wvelL*wvelL)/csoun, 0.25f));
	}
	if(ABS(MaR)>=1.0f)
	{
		Ma_minus = 0.5f*(MaR-ABS(MaR));
		fR = 0.0f;
	}
	else
	{
		Ma_minus = -0.25f*(MaR-1.0f)*(MaR-1.0f) - 0.125f*(MaR*MaR-1.0f)*(MaR*MaR-1.0f);
		fR = (pR/p-1.0f)*ppR*ABS(-0.25f*(MaR-1.0f)*(MaR-1.0f))*MIN(1.0f, POW(SQRT(uvelR*uvelR+vvelR*vvelR+wvelR*wvelR)/csoun, 0.25f));
	}

	Ma = (1.0f+fL)*Ma_plus + (1.0f+fR)*Ma_minus;
	Ma_abs = ABS(Ma);
	if(Ma_abs<=delta)
		Ma_abs = 0.5f*(Ma*Ma + delta*delta)/delta;

	fp[1] = csoun*densL;
	fp[2] = csoun*densL*uvelL;
	fp[3] = csoun*densL*vvelL;
	fp[4] = csoun*densL*wvelL;
	fp[5] = csoun*densL*HL;

	fm[1] = csoun*densR;
	fm[2] = csoun*densR*uvelR;
	fm[3] = csoun*densR*vvelR;
	fm[4] = csoun*densR*wvelR;
	fm[5] = csoun*densR*HR;

	Fc.dens = 0.5f*Ma*(fp[1]+fm[1]) - 0.5f*Ma_abs*(fm[1]-fp[1]);
	Fc.xmom = 0.5f*Ma*(fp[2]+fm[2]) - 0.5f*Ma_abs*(fm[2]-fp[2]) + p;
	Fc.ymom = 0.5f*Ma*(fp[3]+fm[3]) - 0.5f*Ma_abs*(fm[3]-fp[3]); 
	Fc.zmom = 0.5f*Ma*(fp[4]+fm[4]) - 0.5f*Ma_abs*(fm[4]-fp[4]); 
	Fc.ener = 0.5f*Ma*(fp[5]+fm[5]) - 0.5f*Ma_abs*(fm[5]-fp[5]) + p*Q_uvel_e;

}


void CSpaceDiscr::AUSM_PW_Scheme(int i, int j, int k, int flag, const CGeometry &Geometry, const CFluidProps &FluidProp, TNode &norm, TNode &tau1, TNode &tau2, TConsVars &cv_L, TConsVars &cv_R, TConsVars &Fc)
{
	TConsVars Q_L, Q_R;
	TConsVars TL, TR, GL, GR;
	REAL Gamma = CFluidProps::Gamma;
	REAL densL, uvelL, vvelL, wvelL, enerL, csounL, pL, HL, MaL, fL, ppL;
	REAL densR, uvelR, vvelR, wvelR, enerR, csounR, pR, HR, MaR, fR, ppR;
	REAL csoun, Ma, Ma_plus, Ma_minus, p_plus, p_minus, p, wp, mass, V2;
	REAL Q_uvel_e = 0.0f;

	Coord_Rotate(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Q_L, Q_R, Q_uvel_e);

	densL = Q_L.dens;
	uvelL = Q_L.xmom / Q_L.dens;
	vvelL = Q_L.ymom / Q_L.dens;
	wvelL = Q_L.zmom / Q_L.dens;
	enerL = Q_L.ener / Q_L.dens;
	V2    = uvelL*uvelL + vvelL*vvelL + wvelL*wvelL;
	pL	  = (Gamma-1.0f)*densL*(enerL - 0.5f*V2);
	HL    = enerL + pL/densL;

	densR = Q_R.dens;
	uvelR = Q_R.xmom / Q_R.dens;
	vvelR = Q_R.ymom / Q_R.dens;
	wvelR = Q_R.zmom / Q_R.dens;
	enerR = Q_R.ener / Q_R.dens;
	V2    = uvelR*uvelR + vvelR*vvelR + wvelR*wvelR;
	pR	  = (Gamma-1.0f)*densR*(enerR - 0.5f*V2);
	HR    = enerR + pR/densR;

	//HtL	  = (Gamma*pL/densL)/(Gamma-1.0f) + 0.5f*V2;
	//csounCriticalL = SQRT(2.0f*(Gamma-1.0f)/(Gamma+1.0f)*HtL);
	//csounL = csounCriticalL*csounCriticalL / MAX(ABS(uvelL-Q_uvel_e), csounCriticalL);
	//HtR	  = (Gamma*pR/densR)/(Gamma-1.0f) + 0.5f*V2;
	//csounCriticalR = SQRT(2.0f*(Gamma-1.0f)/(Gamma+1.0f)*HtR);
	//csounR = csounCriticalR*csounCriticalR / MAX(ABS(uvelR-Q_uvel_e), csounCriticalR);

	//csoun  = MIN(csounL, csounR);
		
	csounL = SQRT(Gamma*pL/densL);
	csounR = SQRT(Gamma*pR/densR);
	csoun = 0.5f*(csounL+csounR);
	MaL = (uvelL-Q_uvel_e)/csoun;
	MaR = (uvelR-Q_uvel_e)/csoun;

	if((MIN(pL/pR, pR/pL)>=0.75) && (MIN(pL/pR, pR/pL)<1.0))
	{
		ppL = 4.0f*MIN(pL/pR, pR/pL) - 3.0f;
		ppR = 4.0f*MIN(pL/pR, pR/pL) - 3.0f;
	}
	else
	{
		ppL = 0.0f;
		ppR = 0.0f;
	}

	if(ABS(MaL)>1.0f)
		p_plus = 0.5f*(ABS(MaL)+MaL)/MaL;
	else
		p_plus = 0.25f*(MaL+1.0f)*(MaL+1.0f)*(2.0f-MaL) + 3.0f/16.0f*MaL*(MaL*MaL-1.0f)*(MaL*MaL-1.0f);
		//p_plus = 0.5f*(MaL+1.0f);


	if(ABS(MaR)>1.0f)
		p_minus = 0.5f*(MaR-ABS(MaR))/MaR;
	else
		p_minus = 0.25f*(MaR-1.0f)*(MaR-1.0f)*(2.0f+MaR) - 3.0f/16.0f*MaR*(MaR*MaR-1.0f)*(MaR*MaR-1.0f);	//用此公式更好
	   //p_minus = 0.5f*(1.0f-MaR);

	p = pL*p_plus + pR*p_minus;

	if(ABS(MaL)>=1.0f)
	{
		Ma_plus = 0.5f*(MaL+ABS(MaL));
		fL = 0.0f;
	}
	else
	{
		Ma_plus = 0.25f*(MaL+1.0f)*(MaL+1.0f) + 0.125f*(MaL*MaL-1.0f)*(MaL*MaL-1.0f);
		fL = (pL/p-1.0f)*ppL*ABS(0.25f*(MaL+1.0f)*(MaL+1.0f))*MIN(1.0f, POW(SQRT(uvelL*uvelL+vvelL*vvelL+wvelL*wvelL)/csoun, 0.25f));
	}
	if(ABS(MaR)>=1.0f)
	{
		Ma_minus = 0.5f*(MaR-ABS(MaR));
		fR = 0.0f;
	}
	else
	{
		Ma_minus = -0.25f*(MaR-1.0f)*(MaR-1.0f) - 0.125f*(MaR*MaR-1.0f)*(MaR*MaR-1.0f);
		fR = (pR/p-1.0f)*ppR*ABS(-0.25f*(MaR-1.0f)*(MaR-1.0f))*MIN(1.0f, POW(SQRT(uvelR*uvelR+vvelR*vvelR+wvelR*wvelR)/csoun, 0.25f));
	}

	Ma = Ma_plus + Ma_minus;
	wp = 1.0f - POW(MIN(pL/pR, pR/pL), 3.0f);
	//wp = 0.0f;

	GL.dens = (1.0f-wp)*Q_L.dens + wp*Q_R.dens;
	GL.xmom = (1.0f-wp)*Q_L.xmom + wp*Q_R.xmom;
	GL.ymom = (1.0f-wp)*Q_L.ymom + wp*Q_R.ymom;
	GL.zmom = (1.0f-wp)*Q_L.zmom + wp*Q_R.zmom;
	GL.ener = ((1.0f-wp)*densL + wp*densR)*HL;
		
	GR.dens = (1.0f-wp)*Q_R.dens + wp*Q_L.dens;
	GR.xmom = (1.0f-wp)*Q_R.xmom + wp*Q_L.xmom;
	GR.ymom = (1.0f-wp)*Q_R.ymom + wp*Q_L.ymom;
	GR.zmom = (1.0f-wp)*Q_R.zmom + wp*Q_L.zmom;
	GR.ener = ((1.0f-wp)*densR + wp*densL)*HR;

	TL.dens = Q_L.dens;
	TL.xmom = Q_L.xmom;
	TL.ymom = Q_L.ymom;
	TL.zmom = Q_L.zmom;
	TL.ener = densL*HL;
		
	TR.dens = Q_R.dens;
	TR.xmom = Q_R.xmom;
	TR.ymom = Q_R.ymom;
	TR.zmom = Q_R.zmom;
	TR.ener = densR*HR;

	mass = densL*Ma_plus + densR*Ma_minus;
	if(mass>=0.0)
	{
		Fc.dens = (1.0f+fL)*Ma_plus*csoun*TL.dens + (1.0f+fR)*Ma_minus*csoun*GL.dens;
		Fc.xmom = (1.0f+fL)*Ma_plus*csoun*TL.xmom + (1.0f+fR)*Ma_minus*csoun*GL.xmom + p;
		Fc.ymom = (1.0f+fL)*Ma_plus*csoun*TL.ymom + (1.0f+fR)*Ma_minus*csoun*GL.ymom;
		Fc.zmom = (1.0f+fL)*Ma_plus*csoun*TL.zmom + (1.0f+fR)*Ma_minus*csoun*GL.zmom;
		Fc.ener = (1.0f+fL)*Ma_plus*csoun*TL.ener + (1.0f+fR)*Ma_minus*csoun*GL.ener + p*Q_uvel_e;
	}
	else
	{
		Fc.dens = (1.0f+fL)*Ma_plus*csoun*GR.dens + (1.0f+fR)*Ma_minus*csoun*TR.dens;
		Fc.xmom = (1.0f+fL)*Ma_plus*csoun*GR.xmom + (1.0f+fR)*Ma_minus*csoun*TR.xmom + p;
		Fc.ymom = (1.0f+fL)*Ma_plus*csoun*GR.ymom + (1.0f+fR)*Ma_minus*csoun*TR.ymom;
		Fc.zmom = (1.0f+fL)*Ma_plus*csoun*GR.zmom + (1.0f+fR)*Ma_minus*csoun*TR.zmom;
		Fc.ener = (1.0f+fL)*Ma_plus*csoun*GR.ener + (1.0f+fR)*Ma_minus*csoun*TR.ener + p*Q_uvel_e;
	}
}


void CSpaceDiscr::HLLC_Scheme(int i, int j, int k, int flag, const CGeometry &Geometry, const CFluidProps &FluidProp, TNode &norm, TNode &tau1, TNode &tau2, TConsVars &cv_L, TConsVars &cv_R, TConsVars &Fc)
{
	TConsVars Q_L, Q_R;
	REAL Gamma = CFluidProps::Gamma;
	REAL densL, uvelL, vvelL, wvelL, enerL, csounL, pL;
	REAL densR, uvelR, vvelR, wvelR, enerR, csounR, pR;
	REAL  p_pvrs,  p_star, qqL, qqR, SL, SR, S_star, tmpL, tmpR, tmp;
	TConsVars FcL, FcR, Fc_starL, Fc_starR;
	REAL Q_uvel_e = 0.0f;

	Coord_Rotate(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Q_L, Q_R, Q_uvel_e);

	densL = Q_L.dens;
	uvelL = Q_L.xmom / Q_L.dens;
	vvelL = Q_L.ymom / Q_L.dens;
	wvelL = Q_L.zmom / Q_L.dens;
	enerL = Q_L.ener / Q_L.dens;
	pL	  = (Gamma-1.0f)*densL*(enerL - 0.5f*(uvelL*uvelL + vvelL*vvelL + wvelL*wvelL));
	csounL = SQRT(Gamma*pL/densL);

	densR = Q_R.dens;
	uvelR = Q_R.xmom / Q_R.dens;
	vvelR = Q_R.ymom / Q_R.dens;
	wvelR = Q_R.zmom / Q_R.dens;
	enerR = Q_R.ener / Q_R.dens;
	pR	  = (Gamma-1.0f)*densR*(enerR - 0.5f*(uvelR*uvelR + vvelR*vvelR + wvelR*wvelR));
	csounR = SQRT(Gamma*pR/densR);

	p_pvrs = 0.5f*(pL+pR) - (uvelR-uvelL)*(densL+densR)*(csounL+csounR)*0.125f;
    p_star = MAX(0.0f, p_pvrs);  

	if(p_star<=pL)
		qqL = 1.0f;
	else
		qqL = SQRT(1.0f + ((Gamma+1.0f)/(2.0f*Gamma)) * (p_star/pL-1.0f));
	if(p_star<=pR)
		qqR = 1.0f;
	else
		qqR = SQRT(1.0f + ((Gamma+1.0f)/(2.0f*Gamma)) * (p_star/pL-1.0f));

	SL = uvelL - csounL*qqL;			//speed of the left and the right shockwaves
	SR = uvelR + csounR*qqR;

	FcL.dens = Q_L.xmom;	FcL.xmom = Q_L.xmom*uvelL+pL;	FcL.ymom = Q_L.ymom*uvelL;	FcL.zmom = Q_L.zmom*uvelL;	FcL.ener = (Q_L.ener+pL)*uvelL;
	FcR.dens = Q_R.xmom;	FcR.xmom = Q_R.xmom*uvelR+pR;	FcR.ymom = Q_R.ymom*uvelR;	FcR.zmom = Q_R.zmom*uvelR;	FcR.ener = (Q_R.ener+pR)*uvelR;

    S_star = (pR-pL + densL*uvelL*(SL-uvelL) - densR*uvelR*(SR-uvelR)) / (densL*(SL-uvelL) - densR*(SR-uvelR));  
	tmpL = densL*(SL-uvelL)/(SL-S_star);
	tmpR = densR*(SR-uvelR)/(SR-S_star);

	Fc_starL.dens = FcL.dens + SL*(tmpL-Q_L.dens); 
	Fc_starL.xmom = FcL.xmom + SL*(tmpL*S_star-Q_L.xmom); 
	Fc_starL.ymom = FcL.ymom + SL*(tmpL*vvelL-Q_L.ymom); 
	Fc_starL.zmom = FcL.zmom + SL*(tmpL*wvelL-Q_L.zmom); 
	tmp = Q_L.ener/densL + (S_star-uvelL)*(S_star+pL/(densL*(SL-uvelL)));
	Fc_starL.ener = FcL.ener + SL*(tmpL*tmp-Q_L.ener); 
		
	Fc_starR.dens = FcR.dens + SR*(tmpR-Q_R.dens); 
	Fc_starR.xmom = FcR.xmom + SR*(tmpR*S_star-Q_R.xmom); 
	Fc_starR.ymom = FcR.ymom + SR*(tmpR*vvelR-Q_R.ymom); 
	Fc_starR.zmom = FcR.zmom + SR*(tmpR*wvelR-Q_R.zmom); 
	tmp = Q_R.ener/densR + (S_star-uvelR)*(S_star+pR/(densR*(SR-uvelR)));
	Fc_starR.ener = FcR.ener + SR*(tmpR*tmp-Q_R.ener); 

	if(SL>=0.0)
		Fc = FcL;
	else if(SR<=0.0)
		Fc = FcR;
	else if(S_star>=0)
		Fc = Fc_starL;
	else
		Fc = Fc_starR;	     
}


//坐标旋转，但是用的是李新亮老师的方法
void CSpaceDiscr::RoeScheme(int i, int j, int k, int flag, const CGeometry &Geometry, const CFluidProps &FluidProp, TNode &norm, TNode &tau1, TNode &tau2,TConsVars &cv_L, TConsVars &cv_R, TConsVars &Fc)
{

	REAL Gamma = CFluidProps::Gamma;
	REAL densL, uvelL, vvelL, wvelL, EL, pL, HL;
	REAL densR, uvelR, vvelR, wvelR, ER, pR, HR; 
	REAL dens_av, uvel_av, vvel_av, wvel_av, H_av, c_av;	//roe平均量
	REAL delta, lamda1, lamda2, lamda5, tmpp;
	REAL arr[6];
	TConsVars Q_L, Q_R, Fc_L, Fc_R;
	REAL Q_uvel_e = 0.0f;

	Coord_Rotate(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Q_L, Q_R, Q_uvel_e);

	densL = Q_L.dens;
	uvelL = Q_L.xmom / Q_L.dens;
	vvelL = Q_L.ymom / Q_L.dens;
	wvelL = Q_L.zmom / Q_L.dens;
	EL	  = Q_L.ener / Q_L.dens;
	pL	  = (Gamma-1.0f)*densL*(EL - 0.5f*(uvelL*uvelL + vvelL*vvelL + wvelL*wvelL));
	HL    = EL + pL/densL;


	densR = Q_R.dens;
	uvelR = Q_R.xmom / Q_R.dens;
	vvelR = Q_R.ymom / Q_R.dens;
	wvelR = Q_R.zmom / Q_R.dens;
	ER	  = Q_R.ener / Q_R.dens;
	pR	  = (Gamma-1.0f)*densR*(ER - 0.5f*(uvelR*uvelR + vvelR*vvelR + wvelR*wvelR));
	HR    = ER + pR/densR;

	Fc_L.dens = Q_L.dens*(uvelL-Q_uvel_e);		Fc_L.xmom = Q_L.xmom*(uvelL-Q_uvel_e) + pL;		Fc_L.ymom = Q_L.ymom*(uvelL-Q_uvel_e);		Fc_L.zmom = Q_L.zmom*(uvelL-Q_uvel_e);		Fc_L.ener = (Q_L.ener+pL)*(uvelL-Q_uvel_e) + Q_uvel_e*pL;
	Fc_R.dens = Q_R.dens*(uvelR-Q_uvel_e);		Fc_R.xmom = Q_R.xmom*(uvelR-Q_uvel_e) + pR;		Fc_R.ymom = Q_R.ymom*(uvelR-Q_uvel_e);		Fc_R.zmom = Q_R.zmom*(uvelR-Q_uvel_e);		Fc_R.ener = (Q_R.ener+pR)*(uvelR-Q_uvel_e) + Q_uvel_e*pR;

	//Roe平均
	tmpp = SQRT(densR/densL);
	dens_av = 0.25f*densL*(1.0f+tmpp)*(1.0f+tmpp);
	uvel_av = (uvelL+tmpp*uvelR) / (1.0f+tmpp);
	vvel_av = (vvelL+tmpp*vvelR) / (1.0f+tmpp);
	wvel_av = (wvelL+tmpp*wvelR) / (1.0f+tmpp);
	H_av = (HL+tmpp*HR) / (1.0f+tmpp);
	c_av = SQRT((Gamma-1.0f)*(H_av-0.5f*(uvel_av*uvel_av+vvel_av*vvel_av+wvel_av*wvel_av)));

	lamda1 = ABS((uvel_av-Q_uvel_e)-c_av);
	lamda2 = ABS((uvel_av-Q_uvel_e));
	lamda5 = ABS((uvel_av-Q_uvel_e)+c_av);
	//delta = epsEntr*(ABS((uvel_av-Q_uvel_e))+c_av);	//也许可行，在试验中
	delta = 0.5f;

	lamda1 = EntropyCorr(lamda1, delta);
	lamda2 = EntropyCorr(lamda2, delta);
	lamda5 = EntropyCorr(lamda5, delta);

	arr[2] = (Q_R.ymom-Q_L.ymom) / (vvel_av+1E-8f) - (densR-densL);
	//arr[3] = (Q_R.zmom-Q_L.zmom) / (wvel_av) - (Q_R.xmom-Q_L.xmom);		//李老师的程序中的bug
	arr[3] = (Q_R.zmom-Q_L.zmom) / (wvel_av+1E-8f) - (densR-densL);
	arr[4] = (Gamma-1.0f)*((densR-densL)*(H_av-uvel_av*uvel_av-vvel_av*vvel_av-wvel_av*wvel_av) + uvel_av*(Q_R.xmom-Q_L.xmom) + 
			 vvel_av*(Q_R.ymom-Q_L.ymom) + wvel_av*(Q_R.zmom-Q_L.zmom) -(Q_R.ener-Q_L.ener)) / (c_av*c_av);
	arr[1] = ((uvel_av+c_av)*(densR-densL) - (Q_R.xmom-Q_L.xmom) - c_av*arr[4]) / (2.0f*c_av);
	arr[5] = densR-densL - (arr[1]+arr[4]);

	Fc.dens = 0.5f*(Fc_L.dens+Fc_R.dens) - 0.5f*(lamda1*arr[1]+lamda2*arr[4]+lamda5*arr[5]);
	Fc.xmom = 0.5f*(Fc_L.xmom+Fc_R.xmom) - 0.5f*(lamda1*arr[1]*(uvel_av-c_av)+lamda2*arr[4]*uvel_av+lamda5*arr[5]*(uvel_av+c_av));
	Fc.ymom = 0.5f*(Fc_L.ymom+Fc_R.ymom) - 0.5f*(lamda1*arr[1]*vvel_av+lamda2*arr[2]*vvel_av+lamda2*arr[4]*vvel_av+lamda5*arr[5]*vvel_av);
	Fc.zmom = 0.5f*(Fc_L.zmom+Fc_R.zmom) - 0.5f*(lamda1*arr[1]*wvel_av+lamda2*arr[3]*wvel_av+lamda2*arr[4]*wvel_av+lamda5*arr[5]*wvel_av);
	Fc.ener = 0.5f*(Fc_L.ener+Fc_R.ener) - 0.5f*(lamda1*arr[1]*(H_av-uvel_av*c_av)+lamda2*arr[2]*vvel_av*vvel_av+lamda2*arr[3]*wvel_av*wvel_av + 
			lamda2*arr[4]*0.5f*(uvel_av*uvel_av+vvel_av*vvel_av+wvel_av*wvel_av)+lamda5*arr[5]*(H_av+uvel_av*c_av));	 
}

//坐标旋转，用的是自己的方法
//void CSpaceDiscr::RoeScheme(int i, int j, int k, int flag, const CGeometry &Geometry, const CFluidProps &FluidProp, TNode &norm, TNode &tau1, TNode &tau2,TConsVars &cv_L, TConsVars &cv_R, TConsVars &Fc)
//{
//
//	REAL Gamma = CFluidProps::Gamma;
//	REAL densL, uvelL, vvelL, wvelL, EL, pL, HL;
//	REAL densR, uvelR, vvelR, wvelR, ER, pR, HR; 
//	REAL dens_av, uvel_av, vvel_av, wvel_av, H_av, c_av, q_av;	//roe平均量
//	REAL ddens, duvel, dvvel, dwvel, dp;
//	REAL delta, lamda1, lamda2, lamda5, med;
//	REAL h1, h2, h5;
//	TConsVars Q_L, Q_R, Fc_L, Fc_R, F1, F2, F5;
//	REAL uvel_e = FluidProp.uvel_e[i][j][k];
//	REAL vvel_e = FluidProp.vvel_e[i][j][k];
//	REAL Q_uvel_e;
//
//	Coord_Rotate(i, j, k, flag, Geometry, FluidProp, norm, tau1, tau2, cv_L, cv_R, Q_L, Q_R, Q_uvel_e);
//
//	densL = Q_L.dens;
//	uvelL = Q_L.xmom / Q_L.dens;
//	vvelL = Q_L.ymom / Q_L.dens;
//	wvelL = Q_L.zmom / Q_L.dens;
//	EL	  = Q_L.ener / Q_L.dens;
//	pL	  = (Gamma-1.0f)*densL*(EL - 0.5f*(uvelL*uvelL + vvelL*vvelL + wvelL*wvelL));
//	HL    = EL + pL/densL;
//
//
//	densR = Q_R.dens;
//	uvelR = Q_R.xmom / Q_R.dens;
//	vvelR = Q_R.ymom / Q_R.dens;
//	wvelR = Q_R.zmom / Q_R.dens;
//	ER	  = Q_R.ener / Q_R.dens;
//	pR	  = (Gamma-1.0f)*densR*(ER - 0.5f*(uvelR*uvelR + vvelR*vvelR + wvelR*wvelR));
//	HR    = ER + pR/densR;
//
//	Fc_L.dens = Q_L.dens*(uvelL-Q_uvel_e);		Fc_L.xmom = Q_L.xmom*(uvelL-Q_uvel_e) + pL;		Fc_L.ymom = Q_L.ymom*(uvelL-Q_uvel_e);		Fc_L.zmom = Q_L.zmom*(uvelL-Q_uvel_e);		Fc_L.ener = (Q_L.ener+pL)*(uvelL-Q_uvel_e) + Q_uvel_e*pL;
//	Fc_R.dens = Q_R.dens*(uvelR-Q_uvel_e);		Fc_R.xmom = Q_R.xmom*(uvelR-Q_uvel_e) + pR;		Fc_R.ymom = Q_R.ymom*(uvelR-Q_uvel_e);		Fc_R.zmom = Q_R.zmom*(uvelR-Q_uvel_e);		Fc_R.ener = (Q_R.ener+pR)*(uvelR-Q_uvel_e) + Q_uvel_e*pR;
//
//	med = SQRT(densR/densL);
//	dens_av = densL*med;
//	uvel_av = (uvelL + uvelR*med) / (1.0f + med);
//	vvel_av = (vvelL + vvelR*med) / (1.0f + med);
//	wvel_av = (wvelL + wvelR*med) / (1.0f + med);
//	H_av    = (HL + HR*med) / (1.0f + med);
//
//	q_av    = uvel_av*uvel_av + vvel_av*vvel_av + wvel_av*wvel_av;
//	c_av    = SQRT((Gamma-1.0f) * (H_av-q_av/2.0f));
//
//	lamda1 = ABS((uvel_av-Q_uvel_e)-c_av);
//	lamda2 = ABS((uvel_av-Q_uvel_e));
//	lamda5 = ABS((uvel_av-Q_uvel_e)+c_av);
//	//delta = epsEntr*lamda5;	//也可以试验一下
//	delta = epsEntr*(ABS((uvel_av-Q_uvel_e))+c_av);
//
//	lamda1 = EntropyCorr(lamda1, delta);
//	lamda2 = EntropyCorr(lamda2, delta);
//	lamda5 = EntropyCorr(lamda5, delta);
//
//	ddens = densR - densL;
//	duvel = uvelR - uvelL;
//	dvvel = vvelR - vvelL;
//	dwvel = wvelR - wvelL;
//	dp	  = pR - pL;
//
//
//	h1 = (dp - dens_av*c_av*duvel) /2.0f/c_av/c_av;
//	h2 = ddens - dp/c_av/c_av;
//	h5 = (dp + dens_av*c_av*duvel) /2.0f/c_av/c_av;
//
//	F1.dens = lamda1 * h1;
//	F1.xmom = lamda1 * h1 * (uvel_av - c_av);
//	F1.ymom = lamda1 * h1 * (vvel_av);
//	F1.zmom = lamda1 * h1 * (wvel_av);
//	F1.ener = lamda1 * h1 * (H_av    - c_av*uvel_av);
//
//	F2.dens = lamda2 * h2;
//	F2.xmom = lamda2 * (h2*uvel_av);
//	F2.ymom = lamda2 * (h2*vvel_av + dens_av*dvvel);
//	F2.zmom = lamda2 * (h2*wvel_av + dens_av*dwvel);
//	F2.ener = lamda2 * (h2*q_av/2.0f + dens_av*(vvel_av*dvvel + wvel_av*dwvel));
//
//	F5.dens = lamda5 * h5;
//	F5.xmom = lamda5 * h5 * (uvel_av + c_av);
//	F5.ymom = lamda5 * h5 * (vvel_av);
//	F5.zmom = lamda5 * h5 * (wvel_av);
//	F5.ener = lamda5 * h5 * (H_av + c_av*uvel_av);
//
//	Fc=0.5f*(Fc_L+Fc_R - (F1+F2+F5)); 
//}


void CSpaceDiscr::Coord_Rotate(int i, int j, int k, int flag, const CGeometry &Geometry, const CFluidProps &FluidProp, TNode &norm, TNode &tau1, TNode &tau2, TConsVars &cv_L, TConsVars &cv_R, TConsVars &Q_L, TConsVars &Q_R, REAL &Q_uvel_e)
{
	TNode r;
	REAL f_uvel_e, f_vvel_e;

	if(flag==100)
	{
		norm = Geometry.FaceI[i][j][k].S / Geometry.FaceI[i][j][k].s;
		r    = Geometry.Node_coord[i][j+1][k+1] - Geometry.Node_coord[i][j][k];
		f_uvel_e = -1.0f*Geometry.FaceI[i][j][k].coord.y * FluidProp.w_rot;
		f_vvel_e = Geometry.FaceI[i][j][k].coord.x * FluidProp.w_rot;
	}
	else if(flag==200)
	{
		norm = Geometry.FaceJ[i][j][k].S / Geometry.FaceJ[i][j][k].s;
		r    = Geometry.Node_coord[i+1][j][k+1] - Geometry.Node_coord[i][j][k];
		f_uvel_e = -1.0f*Geometry.FaceJ[i][j][k].coord.y * FluidProp.w_rot;
		f_vvel_e = Geometry.FaceJ[i][j][k].coord.x * FluidProp.w_rot;
	}
	else	//flag==300
	{
		norm = Geometry.FaceK[i][j][k].S / Geometry.FaceK[i][j][k].s;
		r    = Geometry.Node_coord[i+1][j+1][k] - Geometry.Node_coord[i][j][k];
		f_uvel_e = -1.0f*Geometry.FaceK[i][j][k].coord.y * FluidProp.w_rot;
		f_vvel_e = Geometry.FaceK[i][j][k].coord.x * FluidProp.w_rot;
	}

	tau1 = r / SQRT(r.x*r.x + r.y*r.y + r.z*r.z);
	r.x = norm.y*tau1.z - norm.z*tau1.y;				//单位化！！！
	r.y = norm.z*tau1.x - norm.x*tau1.z;
	r.z = norm.x*tau1.y - norm.y*tau1.x;
	tau2 = r / SQRT(r.x*r.x + r.y*r.y + r.z*r.z);

	// left
	Q_L.dens = cv_L.dens;			
	Q_L.xmom = cv_L.xmom*norm.x + cv_L.ymom*norm.y + cv_L.zmom*norm.z;
	Q_L.ymom = cv_L.xmom*tau1.x + cv_L.ymom*tau1.y + cv_L.zmom*tau1.z;
	Q_L.zmom = cv_L.xmom*tau2.x + cv_L.ymom*tau2.y + cv_L.zmom*tau2.z;
	Q_L.ener = cv_L.ener;	

	// right
	Q_R.dens = cv_R.dens;			
	Q_R.xmom = cv_R.xmom*norm.x + cv_R.ymom*norm.y + cv_R.zmom*norm.z;
	Q_R.ymom = cv_R.xmom*tau1.x + cv_R.ymom*tau1.y + cv_R.zmom*tau1.z;
	Q_R.zmom = cv_R.xmom*tau2.x + cv_R.ymom*tau2.y + cv_R.zmom*tau2.z;
	Q_R.ener = cv_R.ener;

	Q_uvel_e = f_uvel_e*norm.x + f_vvel_e*norm.y + 0.0f*norm.z;
}


