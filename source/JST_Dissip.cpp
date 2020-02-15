


#include "SpaceDiscr.h"

using namespace std;

REAL CSpaceDiscr::sensorI(int i, int j, int k, const CFluidProps &FluidProp)
{
	REAL vv, PF, PB, PC;

	if(j<1 || j>JM || k<1 || k>KM)
	{
		vv = 0.0;
	}
	else
	{
		if(i==1)
		{
			PF = FluidProp.Cell_pv[2][j][k].press;
			PC = FluidProp.Cell_pv[1][j][k].press;
			PB = FluidProp.DummyIS_pv[1][j][k].press;
			vv = ABS((PF-2.0f*PC+PB)/(PF+2.0f*PC+PB));
		}
		else if(i==0)
		{
			PF = FluidProp.Cell_pv[1][j][k].press;
			PC = FluidProp.DummyIS_pv[1][j][k].press;
			PB = FluidProp.DummyIS_pv[2][j][k].press;
			vv = ABS((PF-2.0f*PC+PB)/(PF+2.0f*PC+PB));
		}
		else if(i==(-1))
		{
			PF = FluidProp.DummyIS_pv[1][j][k].press;
			PC = FluidProp.DummyIS_pv[2][j][k].press;
			PB = FluidProp.DummyIS_pv[3][j][k].press;
			vv = ABS((PF-2.0f*PC+PB)/(PF+2.0f*PC+PB));
		}
		else if(i==(IM-1))
		{
			PF = FluidProp.DummyIP_pv[1][j][k].press;
			PC = FluidProp.Cell_pv[IM-1][j][k].press;
			PB = FluidProp.Cell_pv[IM-2][j][k].press;
			vv = ABS((PF-2.0f*PC+PB)/(PF+2.0f*PC+PB));
		}
		else if(i==IM)
		{
			PF = FluidProp.DummyIP_pv[2][j][k].press;
			PC = FluidProp.DummyIP_pv[1][j][k].press;
			PB = FluidProp.Cell_pv[IM-1][j][k].press;
			vv = ABS((PF-2.0f*PC+PB)/(PF+2.0f*PC+PB));
		}
		else if(i==(IM+1))
		{
			PF = FluidProp.DummyIP_pv[3][j][k].press;
			PC = FluidProp.DummyIP_pv[2][j][k].press;
			PB = FluidProp.DummyIP_pv[1][j][k].press;
			vv = ABS((PF-2.0f*PC+PB)/(PF+2.0f*PC+PB));
		}
		else
		{
			PF = FluidProp.Cell_pv[i+1][j][k].press;
			PC = FluidProp.Cell_pv[i][j][k].press;
			PB = FluidProp.Cell_pv[i-1][j][k].press;
			vv = ABS((PF-2.0f*PC+PB)/(PF+2.0f*PC+PB));
		}		
	}
	return(vv);
}

REAL CSpaceDiscr::sensorJ(int i, int j, int k, const CFluidProps &FluidProp)
{
	REAL vv, PF, PB, PC;

	if(i<1 || i>IM || k<1 || k>KM)
	{
		vv = 0.0;
	}
	else
	{
		if(j==1)
		{
			PF = FluidProp.Cell_pv[i][2][k].press;
			PC = FluidProp.Cell_pv[i][1][k].press;
			PB = FluidProp.DummyJS_pv[i][1][k].press;
			vv = ABS((PF-2.0f*PC+PB)/(PF+2.0f*PC+PB));
		}
		else if(j==0)
		{
			PF = FluidProp.Cell_pv[i][1][k].press;
			PC = FluidProp.DummyJS_pv[i][1][k].press;
			PB = FluidProp.DummyJS_pv[i][2][k].press;
			vv = ABS((PF-2.0f*PC+PB)/(PF+2.0f*PC+PB));
		}
		else if(j==(-1))
		{
			PF = FluidProp.DummyJS_pv[i][1][k].press;
			PC = FluidProp.DummyJS_pv[i][2][k].press;
			PB = FluidProp.DummyJS_pv[i][3][k].press;
			vv = ABS((PF-2.0f*PC+PB)/(PF+2.0f*PC+PB));
		}
		else if(j==(JM-1))
		{
			PF = FluidProp.DummyJP_pv[i][1][k].press;
			PC = FluidProp.Cell_pv[i][JM-1][k].press;
			PB = FluidProp.Cell_pv[i][JM-2][k].press;
			vv = ABS((PF-2.0f*PC+PB)/(PF+2.0f*PC+PB));
		}
		else if(j==JM)
		{
			PF = FluidProp.DummyJP_pv[i][2][k].press;
			PC = FluidProp.DummyJP_pv[i][1][k].press;
			PB = FluidProp.Cell_pv[i][JM-1][k].press;
			vv = ABS((PF-2.0f*PC+PB)/(PF+2.0f*PC+PB));
		}
		else if(j==(JM+1))
		{
			PF = FluidProp.DummyJP_pv[i][3][k].press;
			PC = FluidProp.DummyJP_pv[i][2][k].press;
			PB = FluidProp.DummyJP_pv[i][1][k].press;
			vv = ABS((PF-2.0f*PC+PB)/(PF+2.0f*PC+PB));
		}
		else
		{
			PF = FluidProp.Cell_pv[i][j+1][k].press;
			PC = FluidProp.Cell_pv[i][j][k].press;
			PB = FluidProp.Cell_pv[i][j-1][k].press;
			vv = ABS((PF-2.0f*PC+PB)/(PF+2.0f*PC+PB));
		}		
	}
	return(vv);
}

REAL CSpaceDiscr::sensorK(int i, int j, int k, const CFluidProps &FluidProp)
{
	REAL vv, PF, PB, PC;

	if(i<1 || i>IM || j<1 || j>KM)
	{
		vv = 0.0;
	}
	else
	{
		if(k==1)
		{
			PF = FluidProp.Cell_pv[i][j][2].press;
			PC = FluidProp.Cell_pv[i][j][1].press;
			PB = FluidProp.DummyKS_pv[i][j][1].press;
			vv = ABS((PF-2.0f*PC+PB)/(PF+2.0f*PC+PB));
		}
		else if(k==0)
		{
			PF = FluidProp.Cell_pv[i][j][1].press;
			PC = FluidProp.DummyKS_pv[i][j][1].press;
			PB = FluidProp.DummyKS_pv[i][j][2].press;
			vv = ABS((PF-2.0f*PC+PB)/(PF+2.0f*PC+PB));
		}
		else if(k==(-1))
		{
			PF = FluidProp.DummyKS_pv[i][j][1].press;
			PC = FluidProp.DummyKS_pv[i][j][2].press;
			PB = FluidProp.DummyKS_pv[i][j][3].press;
			vv = ABS((PF-2.0f*PC+PB)/(PF+2.0f*PC+PB));
		}
		else if(k==(KM-1))
		{
			PF = FluidProp.DummyKP_pv[i][j][1].press;
			PC = FluidProp.Cell_pv[i][j][KM-1].press;
			PB = FluidProp.Cell_pv[i][j][KM-2].press;
			vv = ABS((PF-2.0f*PC+PB)/(PF+2.0f*PC+PB));
		}
		else if(k==KM)
		{
			PF = FluidProp.DummyKP_pv[i][j][2].press;
			PC = FluidProp.DummyKP_pv[i][j][1].press;
			PB = FluidProp.Cell_pv[i][j][KM-1].press;
			vv = ABS((PF-2.0f*PC+PB)/(PF+2.0f*PC+PB));
		}
		else if(k==(KM+1))
		{
			PF = FluidProp.DummyKP_pv[i][j][3].press;
			PC = FluidProp.DummyKP_pv[i][j][2].press;
			PB = FluidProp.DummyKP_pv[i][j][1].press;
			vv = ABS((PF-2.0f*PC+PB)/(PF+2.0f*PC+PB));
		}
		else
		{
			PF = FluidProp.Cell_pv[i][j][k+1].press;
			PC = FluidProp.Cell_pv[i][j][k].press;
			PB = FluidProp.Cell_pv[i][j][k-1].press;
			vv = ABS((PF-2.0f*PC+PB)/(PF+2.0f*PC+PB));
		}		
	}
	return(vv);
}

void CSpaceDiscr::JST_Dissip( const CFluidProps &FluidProp)
{
	int i, j, k;

	REAL ***lambda_cI = FluidProp.lambda_cI;
	REAL ***lambda_cJ = FluidProp.lambda_cJ;
	REAL ***lambda_cK = FluidProp.lambda_cK;

	TConsVars cv_1, cv_2, cv_3, cv_4;
	REAL alpha, phi_1, phi_2, VA_1, VA_2, epsilon2, epsilon4;
	
	 //i方向耗散项计算
	for(i=1;i<=IM;i++)
	{
		for(j=1;j<=JM-1;j++)
		{
			for(k=1;k<=KM-1;k++)
			{
				if(i==1)
				{
					phi_1 = 1.0f + POW(lambda_cJ[i][j][k]/lambda_cI[i][j][k],0.4f) + POW(lambda_cK[i][j][k]/lambda_cI[i][j][k],0.4f);
					VA_1  = phi_1*lambda_cI[i][j][k];
					alpha = VA_1;
				}
				else if(i==IM)
				{
					phi_1 = 1.0f + POW(lambda_cJ[i-1][j][k]/lambda_cI[i-1][j][k],0.4f) + POW(lambda_cK[i-1][j][k]/lambda_cI[i-1][j][k],0.4f);
					VA_1  = phi_1*lambda_cI[i-1][j][k];
					alpha = VA_1;
				}
				else
				{
					phi_1 = 1.0f + POW(lambda_cJ[i][j][k]/lambda_cI[i][j][k],0.4f) + POW(lambda_cK[i][j][k]/lambda_cI[i][j][k],0.4f);
					phi_2 = 1.0f + POW(lambda_cJ[i-1][j][k]/lambda_cI[i-1][j][k],0.4f) + POW(lambda_cK[i-1][j][k]/lambda_cI[i-1][j][k],0.4f);
					VA_1  = phi_1*lambda_cI[i][j][k];
					VA_2  = phi_2*lambda_cI[i-1][j][k];
					alpha = 0.5f*(VA_1+VA_2);
				}
				epsilon2 = k2*MAX4(sensorI(i+1,j,k,FluidProp),sensorI(i,j,k,FluidProp),sensorI(i-1,j,k,FluidProp),sensorI(i-2,j,k,FluidProp));
				epsilon4 = MAX(0.0f, k4-epsilon2);
				if(i==1)
				{
					cv_1 = FluidProp.Cell_cv[i+1][j][k];
					cv_2 = FluidProp.Cell_cv[i][j][k];
					cv_3 = FluidProp.DummyIS_cv[1][j][k];
					cv_4 = FluidProp.DummyIS_cv[2][j][k];
				}
				else if(i==2)
				{
					cv_1 = FluidProp.Cell_cv[i+1][j][k];
					cv_2 = FluidProp.Cell_cv[i][j][k];
					cv_3 = FluidProp.Cell_cv[i-1][j][k];
					cv_4 = FluidProp.DummyIS_cv[1][j][k];
				}
				else if(i==IM-1)
				{
					cv_1 = FluidProp.DummyIP_cv[1][j][k];
					cv_2 = FluidProp.Cell_cv[i][j][k];
					cv_3 = FluidProp.Cell_cv[i-1][j][k];
					cv_4 = FluidProp.Cell_cv[i-2][j][k];
				}
				else if(i==IM)
				{
					cv_1 = FluidProp.DummyIP_cv[2][j][k];
					cv_2 = FluidProp.DummyIP_cv[1][j][k];
					cv_3 = FluidProp.Cell_cv[i-1][j][k];
					cv_4 = FluidProp.Cell_cv[i-2][j][k];
				}
				else
				{
					cv_1 = FluidProp.Cell_cv[i+1][j][k];
					cv_2 = FluidProp.Cell_cv[i][j][k];
					cv_3 = FluidProp.Cell_cv[i-1][j][k];
					cv_4 = FluidProp.Cell_cv[i-2][j][k];
				}

				FaceI_Fd[i][j][k] = alpha*(epsilon2*(cv_2-cv_3) - epsilon4*(cv_1-3.0*cv_2 + 3.0*cv_3-cv_4));
			}
		}
	}

	//j方向耗散项计算
	for(i=1;i<=IM-1;i++)
	{
		for(j=1;j<=JM;j++)
		{
			for(k=1;k<=KM-1;k++)
			{
				if(j==1)
				{
					phi_1 = 1.0f + POW(lambda_cI[i][j][k]/lambda_cJ[i][j][k],0.4f) + POW(lambda_cK[i][j][k]/lambda_cJ[i][j][k],0.4f);
					VA_1 = phi_1*lambda_cJ[i][j][k];
					alpha = VA_1;
				}
				else if(j==JM)
				{
					phi_1 = 1.0f + POW(lambda_cI[i][j-1][k]/lambda_cJ[i][j-1][k],0.4f) + POW(lambda_cK[i][j-1][k]/lambda_cJ[i][j-1][k],0.4f);
					VA_1 = phi_1*lambda_cJ[i][j-1][k];
					alpha = VA_1;
				}
				else
				{
					phi_1 = 1.0f + POW(lambda_cI[i][j][k]/lambda_cJ[i][j][k],0.4f) + POW(lambda_cK[i][j][k]/lambda_cJ[i][j][k],0.4f);
					phi_2 = 1.0f + POW(lambda_cI[i][j-1][k]/lambda_cJ[i][j-1][k],0.4f) + POW(lambda_cK[i][j-1][k]/lambda_cJ[i][j-1][k],0.4f);
					VA_1 = phi_1*lambda_cJ[i][j][k];
					VA_2 = phi_2*lambda_cJ[i][j-1][k];
					alpha = 0.5f*(VA_1 + VA_2);
				}
				
				epsilon2 = k2*MAX4(sensorJ(i,j+1,k,FluidProp),sensorJ(i,j,k,FluidProp),sensorJ(i,j-1,k,FluidProp),sensorJ(i,j-2,k,FluidProp));
				epsilon4 = MAX(0.0f, k4-epsilon2);

				if(j==1)
				{
					cv_1 = FluidProp.Cell_cv[i][j+1][k];
					cv_2 = FluidProp.Cell_cv[i][j][k];
					cv_3 = FluidProp.DummyJS_cv[i][1][k];
					cv_4 = FluidProp.DummyJS_cv[i][2][k];
				}
				else if(j==2)
				{
					cv_1 = FluidProp.Cell_cv[i][j+1][k];
					cv_2 = FluidProp.Cell_cv[i][j][k];
					cv_3 = FluidProp.Cell_cv[i][j-1][k];
					cv_4 = FluidProp.DummyJS_cv[i][1][k];							  
				}
				else if(j==JM-1)
				{
					cv_1 = FluidProp.DummyJP_cv[i][1][k];
					cv_2 = FluidProp.Cell_cv[i][j][k];
					cv_3 = FluidProp.Cell_cv[i][j-1][k];
					cv_4 = FluidProp.Cell_cv[i][j-2][k];
				}
				else if(j==JM)
				{
					cv_1 = FluidProp.DummyJP_cv[i][2][k];
					cv_2 = FluidProp.DummyJP_cv[i][1][k];
					cv_3 = FluidProp.Cell_cv[i][j-1][k];
					cv_4 = FluidProp.Cell_cv[i][j-2][k];
				}
				else
				{
					cv_1 = FluidProp.Cell_cv[i][j+1][k];
					cv_2 = FluidProp.Cell_cv[i][j][k];
					cv_3 = FluidProp.Cell_cv[i][j-1][k];
					cv_4 = FluidProp.Cell_cv[i][j-2][k];
				}

				FaceJ_Fd[i][j][k] = alpha*(epsilon2*(cv_2-cv_3) - epsilon4*(cv_1-3.0*cv_2 + 3.0*cv_3-cv_4));
			}
		}
	}

	//k方向耗散项计算
	for(i=1;i<=IM-1;i++)
	{
		for(j=1;j<=JM-1;j++)
		{
			for(k=1;k<=KM;k++)
			{
				if(k==1)
				{
					phi_1 = 1.0f + POW(lambda_cJ[i][j][k]/lambda_cK[i][j][k],0.4f) + POW(lambda_cI[i][j][k]/lambda_cK[i][j][k],0.4f);
					VA_1  = phi_1*lambda_cK[i][j][k];
					alpha = VA_1;
				}
				else if(k==KM)
				{
					phi_1 = 1.0f + POW(lambda_cJ[i][j][k-1]/lambda_cK[i][j][k-1],0.4f) + POW(lambda_cI[i][j][k-1]/lambda_cK[i][j][k-1],0.4f);
					VA_1  = phi_1*lambda_cK[i][j][k-1];
					alpha = VA_1;
				}
				else
				{
					phi_1 = 1.0f + POW(lambda_cJ[i][j][k]/lambda_cK[i][j][k],0.4f) + POW(lambda_cI[i][j][k]/lambda_cK[i][j][k],0.4f);
					phi_2 = 1.0f + POW(lambda_cJ[i][j][k-1]/lambda_cK[i][j][k-1],0.4f) + POW(lambda_cI[i][j][k-1]/lambda_cK[i][j][k-1],0.4f);
					VA_1  = phi_1*lambda_cK[i][j][k];
					VA_2  = phi_2*lambda_cK[i][j][k-1];
					alpha = 0.5f*(VA_1+VA_2);
				}
				epsilon2 = k2*MAX4(sensorK(i,j,k+1,FluidProp),sensorK(i,j,k,FluidProp),sensorK(i,j,k-1,FluidProp),sensorK(i,j,k-2,FluidProp));
				epsilon4 = MAX(0.0f, k4-epsilon2);

				if(k==1)
				{
					cv_1 = FluidProp.Cell_cv[i][j][k+1];
					cv_2 = FluidProp.Cell_cv[i][j][k];
					cv_3 = FluidProp.DummyKS_cv[i][j][1];
					cv_4 = FluidProp.DummyKS_cv[i][j][2];		
				}
				else if(k==2)
				{
					cv_1 = FluidProp.Cell_cv[i][j][k+1];
					cv_2 = FluidProp.Cell_cv[i][j][k];
					cv_3 = FluidProp.Cell_cv[i][j][k-1];
					cv_4 = FluidProp.DummyKS_cv[i][j][1];		
				}
				else if(k==KM-1)
				{
					cv_1 = FluidProp.DummyKP_cv[i][j][1];
					cv_2 = FluidProp.Cell_cv[i][j][k];
					cv_3 = FluidProp.Cell_cv[i][j][k-1];
					cv_4 = FluidProp.Cell_cv[i][j][k-2];
				}
				else if(k==KM)
				{
					cv_1 = FluidProp.DummyKP_cv[i][j][2];
					cv_2 = FluidProp.DummyKP_cv[i][j][1];
					cv_3 = FluidProp.Cell_cv[i][j][k-1];
					cv_4 = FluidProp.Cell_cv[i][j][k-2];
				}
				else
				{
					cv_1 = FluidProp.Cell_cv[i][j][k+1];
					cv_2 = FluidProp.Cell_cv[i][j][k];
					cv_3 = FluidProp.Cell_cv[i][j][k-1];
					cv_4 = FluidProp.Cell_cv[i][j][k-2];
				}

				FaceK_Fd[i][j][k] = alpha*(epsilon2*(cv_2-cv_3) - epsilon4*(cv_1-3.0*cv_2 + 3.0*cv_3-cv_4));
			}
		}
	}

	//cell dissipation 
	for(i=1;i<=IM-1;i++)
	{
		for(j=1;j<=JM-1;j++)
		{
			for(k=1;k<=KM-1;k++)
			{
				Cell_Fd[i][j][k] = FaceI_Fd[i+1][j][k]-FaceI_Fd[i][j][k] + FaceJ_Fd[i][j+1][k]-FaceJ_Fd[i][j][k] + FaceK_Fd[i][j][k+1]-FaceK_Fd[i][j][k];
			}
		}
	}
}