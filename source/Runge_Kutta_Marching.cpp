

#include "Mesh.h"
#include <iostream>
#include <omp.h>
using namespace std;


void CMesh::Runge_Kutta3_Marching()		//用于迎风格式
{
	int i, j, k, l, nb;
	REAL alpha[4];

	alpha[1]=0.1918f;	
	alpha[2]=0.4929f;	
	alpha[3]=1.0f;	
	
	for(nb=1;nb<=NB;nb++)
	{
#pragma omp parallel for num_threads(NT) private(i, j, k)
		for(i=1;i<=Block[nb].Geometry.IM-1;i++)
		{
			for(j=1;j<=Block[nb].Geometry.JM-1;j++)
			{
				for(k=1;k<=Block[nb].Geometry.KM-1;k++)
				{
					Block[nb].TimeDiscr.Cell_cv_old[i][j][k] = Block[nb].FluidProp.Cell_cv[i][j][k];
					if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
					{
						Block[nb].TimeDiscr.Cell_ev_old[i][j][k] = Block[nb].FluidProp.Cell_ev[i][j][k];
					}
				}
			}
		}
	}

	for(l=1;l<=3;l++)
	{	
		ResFlux();
		if((UseForceFunction == true) && (l==1))
		{
			ComputeForceFunction();
		}
		AddForceFunction();

		for(nb=1;nb<=NB;nb++)
		{
#pragma omp parallel for num_threads(NT) default(shared) private(i, j, k)
			for(i=1;i<=Block[nb].Geometry.IM-1;i++)
			{
				for(j=1;j<=Block[nb].Geometry.JM-1;j++)
				{
					for(k=1;k<=Block[nb].Geometry.KM-1;k++)
					{	
						Block[nb].SpaceDiscr.Res_cv[i][j][k] = alpha[l]*Block[nb].TimeDiscr.dt[i][j][k]*Block[nb].SpaceDiscr.Res_cv[i][j][k]/Block[nb].Geometry.Cell[i][j][k].Vol;
						if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
						{
							Block[nb].SpaceDiscr.Res_ev[i][j][k].ev = alpha[l]*Block[nb].TimeDiscr.dt[i][j][k]*Block[nb].SpaceDiscr.Res_ev[i][j][k].ev/Block[nb].Geometry.Cell[i][j][k].Vol;
						}
					}
				}
			}
		}
		if(CTimeDiscr::IRS==1)
			IRS();

		for(nb=1;nb<=NB;nb++)
		{
#pragma omp parallel for num_threads(NT) default(shared) private(i, j, k)
			for(i=1;i<=Block[nb].Geometry.IM-1;i++)
			{
				for(j=1;j<=Block[nb].Geometry.JM-1;j++)
				{
					for(k=1;k<=Block[nb].Geometry.KM-1;k++)
					{				
						Block[nb].FluidProp.Cell_cv[i][j][k] = Block[nb].TimeDiscr.Cell_cv_old[i][j][k] - Block[nb].SpaceDiscr.Res_cv[i][j][k];
						if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
						{
							Block[nb].FluidProp.Cell_ev[i][j][k].ev = Block[nb].TimeDiscr.Cell_ev_old[i][j][k].ev - Block[nb].SpaceDiscr.Res_ev[i][j][k].ev;
						}
					}
				}
			}	
		}
		PrimitiveVars();
		if(CFluidProps::limiter_vars==1)
			Limit_Vars();
	}
}

void CMesh::Runge_Kutta4_Marching()		//用于迎风格式
{
	int i, j, k, l, nb;
	REAL alpha[5];

	alpha[1]=0.1084f;	
	alpha[2]=0.2602f;	
	alpha[3]=0.5052f;	
	alpha[4]=1.0f;	
	
	for(nb=1;nb<=NB;nb++)
	{
#pragma omp parallel for num_threads(NT) default(shared) private(i, j, k)
		for(i=1;i<=Block[nb].Geometry.IM-1;i++)
		{
			for(j=1;j<=Block[nb].Geometry.JM-1;j++)
			{
				for(k=1;k<=Block[nb].Geometry.KM-1;k++)
				{
					Block[nb].TimeDiscr.Cell_cv_old[i][j][k] = Block[nb].FluidProp.Cell_cv[i][j][k];
					if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
					{
						Block[nb].TimeDiscr.Cell_ev_old[i][j][k] = Block[nb].FluidProp.Cell_ev[i][j][k];
					}
				}
			}
		}
	}

	for(l=1;l<=4;l++)
	{	
		ResFlux();
		if((UseForceFunction == true) && (l==1))
		{
			ComputeForceFunction();
		}
		AddForceFunction();

		for(nb=1;nb<=NB;nb++)
		{
#pragma omp parallel for num_threads(NT) default(shared) private(i, j, k)
			for(i=1;i<=Block[nb].Geometry.IM-1;i++)
			{
				for(j=1;j<=Block[nb].Geometry.JM-1;j++)
				{
					for(k=1;k<=Block[nb].Geometry.KM-1;k++)
					{	
						Block[nb].SpaceDiscr.Res_cv[i][j][k] = alpha[l]*Block[nb].TimeDiscr.dt[i][j][k]*Block[nb].SpaceDiscr.Res_cv[i][j][k]/Block[nb].Geometry.Cell[i][j][k].Vol;
						if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
						{
							Block[nb].SpaceDiscr.Res_ev[i][j][k].ev = alpha[l]*Block[nb].TimeDiscr.dt[i][j][k]*Block[nb].SpaceDiscr.Res_ev[i][j][k].ev/Block[nb].Geometry.Cell[i][j][k].Vol;
						}
					}
				}
			}
		}
		if(CTimeDiscr::IRS==1)
			IRS();
		for(nb=1;nb<=NB;nb++)
		{
#pragma omp parallel for num_threads(NT) default(shared) private(i, j, k)
			for(i=1;i<=Block[nb].Geometry.IM-1;i++)
			{
				for(j=1;j<=Block[nb].Geometry.JM-1;j++)
				{
					for(k=1;k<=Block[nb].Geometry.KM-1;k++)
					{				
						Block[nb].FluidProp.Cell_cv[i][j][k] = Block[nb].TimeDiscr.Cell_cv_old[i][j][k] - Block[nb].SpaceDiscr.Res_cv[i][j][k];
						if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
						{
							Block[nb].FluidProp.Cell_ev[i][j][k].ev = Block[nb].TimeDiscr.Cell_ev_old[i][j][k].ev - Block[nb].SpaceDiscr.Res_ev[i][j][k].ev;
						}
					}
				}
			}	
		}
		PrimitiveVars();
		if(CFluidProps::limiter_vars==1)
			Limit_Vars();
	}
}

void CMesh::Runge_Kutta5_Marching()		//用于迎风格式
{
	int i, j, k, l, nb;
	REAL alpha[6];

	alpha[1]=0.06954f;	
	alpha[2]=0.1602f;	
	alpha[3]=0.2898f;	
	alpha[4]=0.5056f;	
	alpha[5]=1.0f;
	
	for(nb=1;nb<=NB;nb++)
	{
#pragma omp parallel for num_threads(NT) default(shared) private(i, j, k)
		for(i=1;i<=Block[nb].Geometry.IM-1;i++)
		{
			for(j=1;j<=Block[nb].Geometry.JM-1;j++)
			{
				for(k=1;k<=Block[nb].Geometry.KM-1;k++)
				{
					Block[nb].TimeDiscr.Cell_cv_old[i][j][k] = Block[nb].FluidProp.Cell_cv[i][j][k];
					if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
					{
						Block[nb].TimeDiscr.Cell_ev_old[i][j][k] = Block[nb].FluidProp.Cell_ev[i][j][k];
					}
				}
			}
		}
	}

	for(l=1;l<=5;l++)
	{	
		ResFlux();
		if((UseForceFunction == true) && (l==1))
		{
			ComputeForceFunction();
		}
		AddForceFunction();

		for(nb=1;nb<=NB;nb++)
		{
#pragma omp parallel for num_threads(NT) default(shared) private(i, j, k)
			for(i=1;i<=Block[nb].Geometry.IM-1;i++)
			{
				for(j=1;j<=Block[nb].Geometry.JM-1;j++)
				{
					for(k=1;k<=Block[nb].Geometry.KM-1;k++)
					{	
						Block[nb].SpaceDiscr.Res_cv[i][j][k] = alpha[l]*Block[nb].TimeDiscr.dt[i][j][k]*Block[nb].SpaceDiscr.Res_cv[i][j][k]/Block[nb].Geometry.Cell[i][j][k].Vol;
						if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
						{
							Block[nb].SpaceDiscr.Res_ev[i][j][k].ev = alpha[l]*Block[nb].TimeDiscr.dt[i][j][k]*Block[nb].SpaceDiscr.Res_ev[i][j][k].ev/Block[nb].Geometry.Cell[i][j][k].Vol;
						}
					}
				}
			}
		}
		if(CTimeDiscr::IRS==1)
			IRS();
		for(nb=1;nb<=NB;nb++)
		{
#pragma omp parallel for num_threads(NT) default(shared) private(i, j, k)
			for(i=1;i<=Block[nb].Geometry.IM-1;i++)
			{
				for(j=1;j<=Block[nb].Geometry.JM-1;j++)
				{
					for(k=1;k<=Block[nb].Geometry.KM-1;k++)
					{				
						Block[nb].FluidProp.Cell_cv[i][j][k] = Block[nb].TimeDiscr.Cell_cv_old[i][j][k] - Block[nb].SpaceDiscr.Res_cv[i][j][k];
						if(Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
						{
							Block[nb].FluidProp.Cell_ev[i][j][k].ev = Block[nb].TimeDiscr.Cell_ev_old[i][j][k].ev - Block[nb].SpaceDiscr.Res_ev[i][j][k].ev;
						}
					}
				}
			}	
		}
		PrimitiveVars();
		if(CFluidProps::limiter_vars==1)
			Limit_Vars();
	}
}

void CMesh::Runge_Kutta5_Hybrid_Marching()		//主要用于central，因为耗散项只要在奇数项求，也可用于迎风格式
{
	int i, j, k, l, nb;
	REAL alpha[6], beta[6];
	TConsVars zero;

	alpha[1]=1.0f/4.0f;			beta[1]=1.0f;
	alpha[2]=1.0f/6.0f;			beta[2]=0.0;
	alpha[3]=3.0f/8.0f;			beta[3]=0.56f;
	alpha[4]=1.0f/2.0f;			beta[4]=0.0;
	alpha[5]=1.0f;				beta[5]=0.44f;

	zero.dens = 0.0;
	zero.xmom = 0.0;
	zero.ymom = 0.0;
	zero.zmom = 0.0;
	zero.ener = 0.0;
	
	for(nb=1;nb<=NB;nb++)
	{
#pragma omp parallel for num_threads(NT) default(shared) private(i, j, k)
		for(i=1;i<=Block[nb].Geometry.IM-1;i++)
		{
			for(j=1;j<=Block[nb].Geometry.JM-1;j++)
			{
				for(k=1;k<=Block[nb].Geometry.KM-1;k++)
				{
					Block[nb].TimeDiscr.Cell_cv_old[i][j][k] = Block[nb].FluidProp.Cell_cv[i][j][k];
				}
			}
		}
	}

	for(l=1;l<=5;l++)
	{	
		BoundaryCondition();
		Spectral_curl_grad();
		TimeStep();
		ConvectiveFlux();
		if(l==1 || l==3 || l==5)		
		{
			Turbulence_Model();
			Dissipation();			
			ViscousFlux();
			//#pragma omp parallel for num_threads(number_of_thread)
			for(nb=1;nb<=NB;nb++)
			{
#pragma omp parallel for num_threads(NT) default(shared) private(i, j, k)
				for(i=1;i<=Block[nb].Geometry.IM-1;i++)
				{
					for(j=1;j<=Block[nb].Geometry.JM-1;j++)
					{
						for(k=1;k<=Block[nb].Geometry.KM-1;k++)
						{
							if(l==1)
							{
								Block[nb].TimeDiscr.Cell_Fd_old[i][j][k] = zero;	  				
								Block[nb].TimeDiscr.Cell_Fv_old[i][j][k] = zero;
							}
							
							Block[nb].SpaceDiscr.Cell_Fd[i][j][k] = beta[l]*Block[nb].SpaceDiscr.Cell_Fd[i][j][k] + (1.0f-beta[l])*Block[nb].TimeDiscr.Cell_Fd_old[i][j][k];

							Block[nb].TimeDiscr.Cell_Fd_old[i][j][k]  = Block[nb].SpaceDiscr.Cell_Fd[i][j][k];

							Block[nb].SpaceDiscr.Cell_Fv[i][j][k] = beta[l]*Block[nb].SpaceDiscr.Cell_Fv[i][j][k] + (1.0f-beta[l])*Block[nb].TimeDiscr.Cell_Fv_old[i][j][k];

							Block[nb].TimeDiscr.Cell_Fv_old[i][j][k]  = Block[nb].SpaceDiscr.Cell_Fv[i][j][k];
						}
					}
				}
			}
		}

		for(nb=1;nb<=NB;nb++)
		{
#pragma omp parallel for num_threads(NT) default(shared) private(i, j, k)
			for(i=1;i<=Block[nb].Geometry.IM-1;i++)
			{
				for(j=1;j<=Block[nb].Geometry.JM-1;j++)
				{
					for(k=1;k<=Block[nb].Geometry.KM-1;k++)
					{	
						Block[nb].SpaceDiscr.Res_cv[i][j][k] = alpha[l]*Block[nb].TimeDiscr.dt[i][j][k]*(Block[nb].SpaceDiscr.Cell_Fc[i][j][k] - Block[nb].SpaceDiscr.Cell_Fv[i][j][k] - Block[nb].SpaceDiscr.Cell_Fd[i][j][k])/Block[nb].Geometry.Cell[i][j][k].Vol;
					}
				}
			}
		}

		if(CTimeDiscr::IRS==1)
			IRS();

		for(nb=1;nb<=NB;nb++)
		{
#pragma omp parallel for num_threads(NT) default(shared) private(i, j, k)
			for(i=1;i<=Block[nb].Geometry.IM-1;i++)
			{
				for(j=1;j<=Block[nb].Geometry.JM-1;j++)
				{
					for(k=1;k<=Block[nb].Geometry.KM-1;k++)
					{				
						Block[nb].FluidProp.Cell_cv[i][j][k] = Block[nb].TimeDiscr.Cell_cv_old[i][j][k] - Block[nb].SpaceDiscr.Res_cv[i][j][k];				
						//if(_finite(Block[nb].FluidProp.Cell_cv[i][j][k].dens) == 0 || Block[nb].FluidProp.Cell_cv[i][j][k].dens>=1000 || Block[nb].FluidProp.Cell_cv[i][j][k].dens<=0.0)
						//{
						//	cout<<"Density="<<Block[nb].FluidProp.Cell_cv[i][j][k].dens<<" broken in: Block"<<nb<<" I="<<i<<" J="<<j<<" K="<<k<<endl;
						//	system("pause");
						//}
					}
				}
			}	
		}
		PrimitiveVars();
		if(CFluidProps::limiter_vars==1)
			Limit_Vars();
	}
}