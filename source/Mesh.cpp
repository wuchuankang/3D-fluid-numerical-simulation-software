

#include "Mesh.h"

using namespace std;

int  CMesh::NB  = 0;
TurboTypes CMesh::TurboType;
InitialTypes CMesh::InitialType;

CMesh::CMesh()
{
	Mesh_ID = 1;
	TotalVol = 0.0;
	UseForceFunction = false;

	Block    = nullptr;
	cv_2h    = nullptr;
	delt_cv  = nullptr;
	Res_2h	 = nullptr;
	Qf_2h	 = nullptr;
}

CMesh::~CMesh()
{
	int i, j, nb;
	if(Block != nullptr)
	{
		delete[] Block;
		Block = nullptr;
	}

	if(Mesh_ID == 2)
	{
		if(cv_2h != nullptr)
		{
			for(nb=0;nb<=CMesh::NB;nb++)
			{
				for(i=0;i<=Block[nb].Geometry.IM;i++)
				{
					for(j=0;j<=Block[nb].Geometry.JM;j++)
					{
						delete [] cv_2h[nb][i][j];
					}
					delete [] cv_2h[nb][i];
				}
				delete [] cv_2h[nb];
			}
			delete [] cv_2h;
			cv_2h = nullptr;
		}
				
		if(delt_cv != nullptr)
		{
			for(nb=0;nb<=CMesh::NB;nb++)
			{
				for(i=0;i<=Block[nb].Geometry.IM;i++)
				{
					for(j=0;j<=Block[nb].Geometry.JM;j++)
					{
						delete [] delt_cv[nb][i][j];
					}
					delete [] delt_cv[nb][i];
				}
				delete [] delt_cv[nb];
			}
			delete [] delt_cv;
			delt_cv = nullptr;
		}
				
		if(Res_2h != nullptr)
		{
			for(nb=0;nb<=CMesh::NB;nb++)
			{
				for(i=0;i<=Block[nb].Geometry.IM;i++)
				{
					for(j=0;j<=Block[nb].Geometry.JM;j++)
					{
						delete [] Res_2h[nb][i][j];
					}
					delete [] Res_2h[nb][i];
				}
				delete [] Res_2h[nb];
			}
			delete [] Res_2h;
			Res_2h = nullptr;
		}
				
		if(Qf_2h != nullptr)
		{
			for(nb=0;nb<=CMesh::NB;nb++)
			{
				for(i=0;i<=Block[nb].Geometry.IM;i++)
				{
					for(j=0;j<=Block[nb].Geometry.JM;j++)
					{
						delete [] Qf_2h[nb][i][j];
					}
					delete [] Qf_2h[nb][i];
				}
				delete [] Qf_2h[nb];
			}
			delete [] Qf_2h;
			Qf_2h = nullptr;
		}
	}
}

void CMesh::AllocateMemory_Pro_Inter_Vars()
{
	int i, j, nb;
	if(Mesh_ID == 2)
	{
		cv_2h    = new TConsVars ***[CMesh::NB+1];
		delt_cv  = new TConsVars ***[CMesh::NB+1];
		Res_2h	 = new TConsVars ***[CMesh::NB+1];
		Qf_2h	 = new TConsVars ***[CMesh::NB+1];
		for(nb=0;nb<=CMesh::NB;nb++)
		{
			cv_2h[nb]    = new TConsVars **[Block[nb].Geometry.IM+1];
			delt_cv[nb]  = new TConsVars **[Block[nb].Geometry.IM+1];
			Res_2h[nb]	 = new TConsVars **[Block[nb].Geometry.IM+1];
			Qf_2h[nb]	 = new TConsVars **[Block[nb].Geometry.IM+1];
		}
				
		for(nb=0;nb<=CMesh::NB;nb++)
		{
			for(i=0;i<=Block[nb].Geometry.IM;i++)
			{
				cv_2h[nb][i]    = new TConsVars *[Block[nb].Geometry.JM+1];
				delt_cv[nb][i]  = new TConsVars *[Block[nb].Geometry.JM+1];
				Res_2h[nb][i]	= new TConsVars *[Block[nb].Geometry.JM+1];
				Qf_2h[nb][i]	= new TConsVars *[Block[nb].Geometry.JM+1];
			}
		}
				
		for(nb=0;nb<=CMesh::NB;nb++)
		{
			for(i=0;i<=Block[nb].Geometry.IM;i++)
			{
				for(j=0;j<=Block[nb].Geometry.JM;j++)
				{
					cv_2h[nb][i][j]    = new TConsVars [Block[nb].Geometry.KM+1];
					delt_cv[nb][i][j]  = new TConsVars [Block[nb].Geometry.KM+1];
					Res_2h[nb][i][j]   = new TConsVars [Block[nb].Geometry.KM+1];
					Qf_2h[nb][i][j]    = new TConsVars [Block[nb].Geometry.KM+1];
				}
			}
		}
	}
}