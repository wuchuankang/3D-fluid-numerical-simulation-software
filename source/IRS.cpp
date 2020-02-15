

#include "TimeDiscr.h"

using namespace std;

//void CTimeDiscr::CentralIRS(const CFluidProps &FluidProp, CSpaceDiscr &SpaceDiscr)
//{
//	int i, j, k;
//	REAL phi_CIRS, sxtos, rJI, rKI, rKJ, rIJ, rIK, rJK; 
//	
//	phi_CIRS = 0.0625;
//	sxtos    = 2.0f;
//
//
//	//I direction IRS
//	for(k=1;k<=KM-1;k++)
//	{
//		for(j=1;j<=JM-1;j++)
//		{
//			for(i=1;i<=IM-1;i++)
//			{
//				rJI = FluidProp.lambda_cJ[i][j][k]/FluidProp.lambda_cI[i][j][k];
//				rKI = FluidProp.lambda_cK[i][j][k]/FluidProp.lambda_cI[i][j][k];
//				epsI[i]   = MAX(0.25f*(POW(sxtos/(1.0f+phi_CIRS*(rJI+rKI)), 2.0f)-1.0f),0.0f);
//				a_IRS[i] = 0.0f-epsI[i];
//				b_IRS[i] = (1.0f+2.0f*epsI[i]);
//				c_IRS[i] = 0.0f-epsI[i];
//			}
//
//			for(i=1;i<=IM-1;i++)
//				d_IRS[i] = SpaceDiscr.Res_cv[i][j][k];
//			Thomas(IM-1, a_IRS, b_IRS, c_IRS, d_IRS, x_IRS);
//			for(i=1;i<=IM-1;i++)
//				SpaceDiscr.Res_cv[i][j][k] = x_IRS[i];
//		}
//	}
//
//	//J direction IRS
//	for(k=1;k<=KM-1;k++)
//	{
//		for(i=1;i<=IM-1;i++)
//		{
//			for(j=1;j<=JM-1;j++)
//			{
//				rIJ = FluidProp.lambda_cI[i][j][k]/FluidProp.lambda_cJ[i][j][k];
//				rKJ = FluidProp.lambda_cK[i][j][k]/FluidProp.lambda_cJ[i][j][k];
//				epsJ[j]   = MAX(0.25f*(POW(sxtos/(1.0f+phi_CIRS*(rKJ+rIJ)),2.0f)-1.0f),0.0f);
//				a_IRS[j] = 0.0f-epsJ[j];
//				b_IRS[j] = (1.0f+2.0f*epsJ[j]);
//				c_IRS[j] = 0.0f-epsJ[j];
//			}
//
//			for(j=1;j<=JM-1;j++)
//				d_IRS[j] = SpaceDiscr.Res_cv[i][j][k];
//			Thomas(JM-1, a_IRS, b_IRS, c_IRS, d_IRS, x_IRS);
//			for(j=1;j<=JM-1;j++)
//				SpaceDiscr.Res_cv[i][j][k] = x_IRS[j];
//		}
//	}
//
//	//K direction IRS
//	for(i=1;i<=IM-1;i++)
//	{
//		for(j=1;j<=JM-1;j++)
//		{
//			for(k=1;k<=KM-1;k++)
//			{
//				rIK = FluidProp.lambda_cI[i][j][k]/FluidProp.lambda_cK[i][j][k];
//				rJK = FluidProp.lambda_cJ[i][j][k]/FluidProp.lambda_cK[i][j][k];
//				epsK[k]   =MAX(0.25f*(POW(sxtos/(1.0f+phi_CIRS*(rIK+rJK)),2.0f)-1.0f),0.0f);
//				a_IRS[k] = 0.0f-epsK[k];
//				b_IRS[k] = (1.0f+2.0f*epsK[k]);
//				c_IRS[k] = 0.0f-epsK[k];
//			}			  
//			for(k=1;k<=KM-1;k++)
//				d_IRS[k] = SpaceDiscr.Res_cv[i][j][k];
//			Thomas(KM-1, a_IRS, b_IRS, c_IRS, d_IRS, x_IRS);
//			for(k=1;k<=KM-1;k++)
//				SpaceDiscr.Res_cv[i][j][k] = x_IRS[k];
//		}
//	}
//
//}

//void CTimeDiscr::Thomas(int TMAX, REAL *a, REAL *b, REAL *c, TConsVars *d, TConsVars *x)
//{
//	int i;
//	TConsVars *y;
//	REAL *u;
//	REAL *l;
//	y = new TConsVars [TMAX+1];
//	u = new REAL [TMAX+1];
//	l = new REAL [TMAX+1];
//
//	for(i=1;i<=TMAX;i++)
//	{
//		if(i==1)
//		{
//			u[1]=b[1];
//			y[1]=d[1];
//		}
//		else
//		{
//			l[i]=a[i]/u[i-1];
//			u[i]=b[i]-l[i]*c[i-1];
//			y[i]=d[i]-l[i]*y[i-1];
//		}
//	}
//	for(i=TMAX;i>=1;i--)
//	{
//		if(i==TMAX)
//		{
//			x[i]=y[i]/u[i];
//		}
//		else
//		{
//			x[i]=1.0f/u[i]*(y[i]-c[i]*x[i+1]);
//		}
//	}
//
//	delete []y;
//	delete []u;
//	delete []l;
//}

void CTimeDiscr::UpwindIRS(const CGeometry &Geometry, const CFluidProps &FluidProp, CSpaceDiscr &SpaceDiscr)
{
	int i, j, k;
	REAL rJI, rKI, rKJ, rIJ, rIK, rJK, eI, eJ, eK;
	REAL Vn,  csound, Man, eps=1.0f; 
	TNode norm;
#pragma omp parallel num_threads(NT) default(shared) private(i, j, k, rJI, rKI, rKJ, rIJ, rIK, rJK, eI, eJ, eK, Vn, csound, Man, norm)
	{
#pragma omp for nowait
		//I direction IRS
		for(k=1;k<=KM-1;k++)
		{
			for(j=1;j<=JM-1;j++)
			{
				for(i=1;i<=IM-1;i++)
				{
					rIJ = FluidProp.lambda_cI[i][j][k]/FluidProp.lambda_cJ[i][j][k];
					rIK = FluidProp.lambda_cI[i][j][k]/FluidProp.lambda_cK[i][j][k];
					eI = eps*MIN3(1.0f, rIJ, rIK);
					norm.x = (Geometry.FaceI[i][j][k].S.x/Geometry.FaceI[i][j][k].s + Geometry.FaceI[i+1][j][k].S.x/Geometry.FaceI[i+1][j][k].s)/2.0f;
					norm.y = (Geometry.FaceI[i][j][k].S.y/Geometry.FaceI[i][j][k].s + Geometry.FaceI[i+1][j][k].S.y/Geometry.FaceI[i+1][j][k].s)/2.0f;
					norm.z = (Geometry.FaceI[i][j][k].S.z/Geometry.FaceI[i][j][k].s + Geometry.FaceI[i+1][j][k].S.z/Geometry.FaceI[i+1][j][k].s)/2.0f;
					Vn = FluidProp.Cell_pv[i][j][k].uvel*norm.x + FluidProp.Cell_pv[i][j][k].vvel*norm.y + FluidProp.Cell_pv[i][j][k].wvel*norm.z;
					csound = SQRT(CFluidProps::Gamma*CFluidProps::Rgas*FluidProp.Cell_pv[i][j][k].temp);
					Man = Vn/csound;
					if(Man > 1.0f)
					{
						a_IRS[i] = -eI;
						b_IRS[i] = 1.0f + eI;
						c_IRS[i] = 0.0f;
					}
					else if(Man < -1.0f)	
					{
						a_IRS[i] = 0.0f;
						b_IRS[i] = 1.0f + eI;
						c_IRS[i] = -eI;
					}
					else
					{
						a_IRS[i] = -eI;
						b_IRS[i] = 1.0f + 2.0f*eI;
						c_IRS[i] = -eI;
					}
				}
				for(i=1;i<=IM-1;i++)
					R_IRS[i] = SpaceDiscr.Res_cv[i][j][k];
				Thomas(IM-1, a_IRS, b_IRS, c_IRS, R_IRS, Rs_IRS);
				for(i=1;i<=IM-1;i++)
					SpaceDiscr.Res_cv[i][j][k] = Rs_IRS[i];
			}
		}

#pragma omp for nowait
		//J direction IRS
		for(k=1;k<=KM-1;k++)
		{
			for(i=1;i<=IM-1;i++)
			{
				for(j=1;j<=JM-1;j++)
				{
					rJI = FluidProp.lambda_cJ[i][j][k]/FluidProp.lambda_cI[i][j][k];
					rJK = FluidProp.lambda_cJ[i][j][k]/FluidProp.lambda_cK[i][j][k];
					eJ = eps*MIN3(1.0f, rJI, rJK);
					norm.x = (Geometry.FaceJ[i][j][k].S.x/Geometry.FaceJ[i][j][k].s + Geometry.FaceJ[i][j+1][k].S.x/Geometry.FaceJ[i][j+1][k].s)/2.0f;
					norm.y = (Geometry.FaceJ[i][j][k].S.y/Geometry.FaceJ[i][j][k].s + Geometry.FaceJ[i][j+1][k].S.y/Geometry.FaceJ[i][j+1][k].s)/2.0f;
					norm.z = (Geometry.FaceJ[i][j][k].S.z/Geometry.FaceJ[i][j][k].s + Geometry.FaceJ[i][j+1][k].S.z/Geometry.FaceJ[i][j+1][k].s)/2.0f;
					Vn = FluidProp.Cell_pv[i][j][k].uvel*norm.x + FluidProp.Cell_pv[i][j][k].vvel*norm.y + FluidProp.Cell_pv[i][j][k].wvel*norm.z;
					csound = SQRT(CFluidProps::Gamma*CFluidProps::Rgas*FluidProp.Cell_pv[i][j][k].temp);
					Man = Vn/csound;
					if(Man > 1.0f)
					{
						a_IRS[j] = -eJ;
						b_IRS[j] = 1.0f + eJ;
						c_IRS[j] = 0.0f;
					}
					else if(Man < -1.0f)	
					{
						a_IRS[j] = 0.0f;
						b_IRS[j] = 1.0f + eJ;
						c_IRS[j] = -eJ;
					}
					else
					{
						a_IRS[j] = -eJ;
						b_IRS[j] = 1.0f + 2.0f*eJ;
						c_IRS[j] = -eJ;
					}
				}
				for(j=1;j<=JM-1;j++)
					R_IRS[j] = SpaceDiscr.Res_cv[i][j][k];
				Thomas(JM-1, a_IRS, b_IRS, c_IRS, R_IRS, Rs_IRS);
				for(j=1;j<=JM-1;j++)
					SpaceDiscr.Res_cv[i][j][k] = Rs_IRS[j];
			}
		}

#pragma omp for nowait
		//k direction IRS
		for(j=1;j<=JM-1;j++)
		{
			for(i=1;i<=IM-1;i++)
			{
				for(k=1;k<=KM-1;k++)
				{
					rKI = FluidProp.lambda_cK[i][j][k]/FluidProp.lambda_cI[i][j][k];
					rKJ = FluidProp.lambda_cK[i][j][k]/FluidProp.lambda_cJ[i][j][k];
					eK = eps*MIN3(1.0f, rKI, rKJ);
					norm.x = (Geometry.FaceK[i][j][k].S.x/Geometry.FaceK[i][j][k].s + Geometry.FaceK[i][j][k+1].S.x/Geometry.FaceK[i][j][k+1].s)/2.0f;
					norm.y = (Geometry.FaceK[i][j][k].S.y/Geometry.FaceK[i][j][k].s + Geometry.FaceK[i][j][k+1].S.y/Geometry.FaceK[i][j][k+1].s)/2.0f;
					norm.z = (Geometry.FaceK[i][j][k].S.z/Geometry.FaceK[i][j][k].s + Geometry.FaceK[i][j][k+1].S.z/Geometry.FaceK[i][j][k+1].s)/2.0f;
					Vn = FluidProp.Cell_pv[i][j][k].uvel*norm.x + FluidProp.Cell_pv[i][j][k].vvel*norm.y + FluidProp.Cell_pv[i][j][k].wvel*norm.z;
					csound = SQRT(CFluidProps::Gamma*CFluidProps::Rgas*FluidProp.Cell_pv[i][j][k].temp);
					Man = Vn/csound;
					if(Man > 1.0f)
					{
						a_IRS[k] = -eK;
						b_IRS[k] = 1.0f + eK;
						c_IRS[k] = 0.0f;
					}
					else if(Man < -1.0f)	
					{
						a_IRS[k] = 0.0f;
						b_IRS[k] = 1.0f + eK;
						c_IRS[k] = -eK;
					}
					else
					{
						a_IRS[k] = -eK;
						b_IRS[k] = 1.0f + 2.0f*eK;
						c_IRS[k] = -eK;
					}
				}
				for(k=1;k<=KM-1;k++)
					R_IRS[k] = SpaceDiscr.Res_cv[i][j][k];
				Thomas(KM-1, a_IRS, b_IRS, c_IRS, R_IRS, Rs_IRS);
				for(k=1;k<=KM-1;k++)
					SpaceDiscr.Res_cv[i][j][k] = Rs_IRS[k];
			}
		}
	}
}


void CTimeDiscr::Thomas(int n, REAL *a, REAL *b, REAL *c, TConsVars *R, TConsVars *Rs)
{
	int i;
	TConsVars *Q;
	REAL *P, tmp;
	P = new REAL[n+2];
	Q = new TConsVars[n+2];
	P[1] = 0.0;
	P[n] = 0.0;
	Q[1] = R[1];
	Q[n] = R[n];
	for(i=2;i<=n;i++)
	{
		tmp = 1.0f/(a[i]*P[i-1] + b[i]);
		P[i] = -c[i]*tmp;
		Q[i] = (R[i] - a[i]*Q[i-1]) * tmp;
	}
	Rs[n] = R[n];
	Rs[1] = R[1];
	for(i=n-1;i>=2;i--)
	{
		Rs[i] = P[i]*Rs[i+1] + Q[i];
	}

	delete[] P;
	delete[] Q;
}
