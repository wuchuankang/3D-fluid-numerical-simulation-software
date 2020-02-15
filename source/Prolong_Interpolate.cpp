

#include "Solver.h"

using namespace std;

void CSolver::Prolong(CMesh &Mesh_1, CMesh &Mesh_2)			//fine grid prolong to coarse grid
{
	int i, j, k, nb;
	int i1, j1, k1, i2, j2, k2;
	TConsVars cv111, cv112, cv121, cv122, cv211, cv212, cv221, cv222;
	REAL Vol2;
		
	for(nb=1;nb<=CMesh::NB;nb++)
	{
#pragma omp parallel for num_threads(NT) default(shared) private(i, j, k, i1, j1, k1, i2, j2, k2, cv111, cv112, cv121, cv122, cv211, cv212, cv221, cv222, Vol2)
		for(i=1;i<=Mesh_2.Block[nb].Geometry.IM-1;i++)
		{
			for(j=1;j<=Mesh_2.Block[nb].Geometry.JM-1;j++)
			{
				for(k=1;k<=Mesh_2.Block[nb].Geometry.KM-1;k++)
				{
					i1=2*i-1; 		
					j1=2*j-1;
					k1=2*k-1;
					i2=2*i;
					j2=2*j;
					k2=2*k;

					cv111 = Mesh_1.Block[nb].FluidProp.Cell_cv[i1][j1][k1] * Mesh_1.Block[nb].Geometry.Cell[i1][j1][k1].Vol;
					cv112 = Mesh_1.Block[nb].FluidProp.Cell_cv[i1][j1][k2] * Mesh_1.Block[nb].Geometry.Cell[i1][j1][k2].Vol;
					cv121 = Mesh_1.Block[nb].FluidProp.Cell_cv[i1][j2][k1] * Mesh_1.Block[nb].Geometry.Cell[i1][j2][k1].Vol;
					cv122 = Mesh_1.Block[nb].FluidProp.Cell_cv[i1][j2][k2] * Mesh_1.Block[nb].Geometry.Cell[i1][j2][k2].Vol;
					cv211 = Mesh_1.Block[nb].FluidProp.Cell_cv[i2][j1][k1] * Mesh_1.Block[nb].Geometry.Cell[i2][j1][k1].Vol;
					cv212 = Mesh_1.Block[nb].FluidProp.Cell_cv[i2][j1][k2] * Mesh_1.Block[nb].Geometry.Cell[i2][j1][k2].Vol;
					cv221 = Mesh_1.Block[nb].FluidProp.Cell_cv[i2][j2][k1] * Mesh_1.Block[nb].Geometry.Cell[i2][j2][k1].Vol;
					cv222 = Mesh_1.Block[nb].FluidProp.Cell_cv[i2][j2][k2] * Mesh_1.Block[nb].Geometry.Cell[i2][j2][k2].Vol;

					Vol2  = Mesh_2.Block[nb].Geometry.Cell[i][j][k].Vol;
					Mesh_2.cv_2h[nb][i][j][k] = (cv111 + cv112 + cv121 + cv122 + cv211 + cv212 + cv221 + cv222) / Vol2;
										
					Mesh_2.Res_2h[nb][i][j][k] = Mesh_1.Block[nb].SpaceDiscr.Res_cv[i1][j1][k1] + Mesh_1.Block[nb].SpaceDiscr.Res_cv[i1][j1][k2] + Mesh_1.Block[nb].SpaceDiscr.Res_cv[i1][j2][k1] + Mesh_1.Block[nb].SpaceDiscr.Res_cv[i1][j2][k2] + 
												 Mesh_1.Block[nb].SpaceDiscr.Res_cv[i2][j1][k1] + Mesh_1.Block[nb].SpaceDiscr.Res_cv[i2][j1][k2] + Mesh_1.Block[nb].SpaceDiscr.Res_cv[i2][j2][k1] + Mesh_1.Block[nb].SpaceDiscr.Res_cv[i2][j2][k2];

					Mesh_2.Block[nb].FluidProp.Cell_cv[i][j][k] = Mesh_2.cv_2h[nb][i][j][k];
				}
			}
		}
	}
}

void CSolver::Interpolate(CMesh &Mesh_1, CMesh &Mesh_2)		//coarse grid prolong to fine grid
{
	int i, j, k, nb, IMmax, JMmax, KMmax;
	TConsVars Temp1, Temp2,  Temp3, Temp4, Delt_cv;
	REAL a1=27.0f/64.0f, a2=9.0f/64.0f, a3=3.0f/64.0f, a4=1.0f/64.0f;

	int **ia, **ja, **ka;
	ia = new int *[3];
	ja = new int *[3];
	ka = new int *[3];
	//ia(1,i) 是距离i点最近的粗网格点的下标；ia(2,i)是次近点的下标

	IMmax = Mesh_1.Block[1].Geometry.IM;
	JMmax = Mesh_1.Block[1].Geometry.JM;
	KMmax = Mesh_1.Block[1].Geometry.KM;
	for(nb=2;nb<=CMesh::NB;nb++)
	{
		IMmax = MAX(IMmax, Mesh_1.Block[nb].Geometry.IM);
		JMmax = MAX(JMmax, Mesh_1.Block[nb].Geometry.JM);
		KMmax = MAX(KMmax, Mesh_1.Block[nb].Geometry.KM);
	}

	for(i=0;i<=2;i++)
	{
		ia[i] = new int [IMmax];
		ja[i] = new int [JMmax];
		ka[i] = new int [KMmax];
	}

	for(nb=1;nb<=CMesh::NB;nb++)
	{
#pragma omp parallel for num_threads(NT) default(shared) private(i, j, k, Temp1, Temp2,  Temp3, Temp4, Delt_cv)
		for(i=1;i<=Mesh_1.Block[nb].Geometry.IM-1;i++)
		{
			for(j=1;j<=Mesh_1.Block[nb].Geometry.JM-1;j++)
			{
				for(k=1;k<=Mesh_1.Block[nb].Geometry.KM-1;k++)
				{
					if(i%2 ==0)
					{
						ia[1][i] = i/2;
						ia[2][i] = i/2 + 1;
					}
					else
					{
						ia[1][i] = i/2 + 1;
						ia[2][i] = i/2;
					}
					if(j%2 ==0)
					{
						ja[1][j] = j/2;
						ja[2][j] = j/2 + 1;
					}
					else
					{
						ja[1][j] = j/2 + 1;
						ja[2][j] = j/2;
					}
					if(k%2 ==0)
					{
						ka[1][k] = k/2;
						ka[2][k] = k/2 + 1;
					}
					else
					{
						ka[1][k] = k/2 + 1;
						ka[2][k] = k/2;
					}
					if(reUseSolution==false)
					{
						if(iter == Coarse_iter)
						{
							Temp1 = a1 * Mesh_2.Block[nb].FluidProp.Cell_cv[ia[1][i]][ja[1][j]][ka[1][k]];
							Temp2 = a2 * (Mesh_2.Block[nb].FluidProp.Cell_cv[ia[1][i]][ja[1][j]][ka[2][k]] + Mesh_2.Block[nb].FluidProp.Cell_cv[ia[1][i]][ja[2][j]][ka[1][k]] + Mesh_2.Block[nb].FluidProp.Cell_cv[ia[2][i]][ja[1][j]][ka[1][k]]);
							Temp3 = a3 * (Mesh_2.Block[nb].FluidProp.Cell_cv[ia[2][i]][ja[2][j]][ka[1][k]] + Mesh_2.Block[nb].FluidProp.Cell_cv[ia[1][i]][ja[2][j]][ka[2][k]] + Mesh_2.Block[nb].FluidProp.Cell_cv[ia[2][i]][ja[1][j]][ka[2][k]]);
							Temp4 = a4 * Mesh_2.Block[nb].FluidProp.Cell_cv[ia[2][i]][ja[2][j]][ka[2][k]];
							Mesh_1.Block[nb].FluidProp.Cell_cv[i][j][k] = Temp1 + Temp2 + Temp3 + Temp4;
						}
						else
						{
							Temp1 = a1 * Mesh_2.delt_cv[nb][ia[1][i]][ja[1][j]][ka[1][k]];
							Temp2 = a2 * (Mesh_2.delt_cv[nb][ia[1][i]][ja[1][j]][ka[2][k]] + Mesh_2.delt_cv[nb][ia[1][i]][ja[2][j]][ka[1][k]] + Mesh_2.delt_cv[nb][ia[2][i]][ja[1][j]][ka[1][k]]);
							Temp3 = a3 * (Mesh_2.delt_cv[nb][ia[2][i]][ja[2][j]][ka[1][k]] + Mesh_2.delt_cv[nb][ia[1][i]][ja[2][j]][ka[2][k]] + Mesh_2.delt_cv[nb][ia[2][i]][ja[1][j]][ka[2][k]]);
							Temp4 = a4 * Mesh_2.delt_cv[nb][ia[2][i]][ja[2][j]][ka[2][k]];
							Delt_cv = Temp1 + Temp2 + Temp3 + Temp4;
							Mesh_1.Block[nb].FluidProp.Cell_cv[i][j][k] = Mesh_1.Block[nb].FluidProp.Cell_cv[i][j][k] + Delt_cv;	//更新细网格上的守恒量
						}
					}
					else
					{
						Temp1 = a1 * Mesh_2.delt_cv[nb][ia[1][i]][ja[1][j]][ka[1][k]];
						Temp2 = a2 * (Mesh_2.delt_cv[nb][ia[1][i]][ja[1][j]][ka[2][k]] + Mesh_2.delt_cv[nb][ia[1][i]][ja[2][j]][ka[1][k]] + Mesh_2.delt_cv[nb][ia[2][i]][ja[1][j]][ka[1][k]]);
						Temp3 = a3 * (Mesh_2.delt_cv[nb][ia[2][i]][ja[2][j]][ka[1][k]] + Mesh_2.delt_cv[nb][ia[1][i]][ja[2][j]][ka[2][k]] + Mesh_2.delt_cv[nb][ia[2][i]][ja[1][j]][ka[2][k]]);
						Temp4 = a4 * Mesh_2.delt_cv[nb][ia[2][i]][ja[2][j]][ka[2][k]];
						Delt_cv = Temp1 + Temp2 + Temp3 + Temp4;
						Mesh_1.Block[nb].FluidProp.Cell_cv[i][j][k] = Mesh_1.Block[nb].FluidProp.Cell_cv[i][j][k] + Delt_cv;	//更新细网格上的守恒量
					}					
				}
			}
		}
	}	
	delete []ia;
	delete []ja;
	delete []ka;              
}


void CSolver::ComputeDeltConsVars(CMesh &Mesh_2)
{
	int i, j, k, nb;
	for(nb=1;nb<=CMesh::NB;nb++)
	{
#pragma omp parallel for num_threads(NT) private(i, j, k)
		for(i=1;i<=Mesh_2.Block[nb].Geometry.IM-1;i++)
		{
			for(j=1;j<=Mesh_2.Block[nb].Geometry.JM-1;j++)
			{
				for(k=1;k<=Mesh_2.Block[nb].Geometry.KM-1;k++)
				{
					Mesh_2.delt_cv[nb][i][j][k] = Mesh_2.Block[nb].FluidProp.Cell_cv[i][j][k] -  Mesh_2.cv_2h[nb][i][j][k];
				}
			}
		}
	}
}