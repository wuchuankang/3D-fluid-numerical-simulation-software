
/*********************************
说明：这里的壁面上的虚元值是和壁面上的值是相同的，
这是为了重构时，避免壁面上的物理量与周边的量差距大，初始迭代会造成发散
这只是对于迎风格式而言的，对于中心格式则不需要如此处理
虚元的原始变量也是按照壁面条件给定的；

上面方法证明不适用 2016/4/3
现在重新思考：也许还得那么用，再次试一试  2016/4/13

但是先用一阶Roe格式，算到残差降到一定水平后，再用二阶Roe格式算，此时，无需对壁面虚元采用以上处理方式
***********************************/


#include "BndConds.h"
#include "SpaceDiscr.h"
#include "defs.h"


void CBndConds::BndCondWall_NoSlip(const TWall &Wall, const CGeometry &Geometry, CFluidProps &FluidProp)
{
	int i, j, k, iv, jv, kv, nd;

	int is = Wall.is;
	int ie = Wall.ie;
	int js = Wall.js;
	int je = Wall.je;
	int ks = Wall.ks;
	int ke = Wall.ke;
	int IM = Geometry.IM;
	int JM = Geometry.JM;
	int KM = Geometry.KM;
	int NDummy = Geometry.NDummy;
	TFace	    ***FaceI  = Geometry.FaceI;
	TFace	    ***FaceJ  = Geometry.FaceJ;
	TFace	    ***FaceK  = Geometry.FaceK;
	TConsVars   ***Cell_cv	  = FluidProp.Cell_cv; 
	TPrimVars   ***Cell_pv	  = FluidProp.Cell_pv; 
	TEddyVars	***Cell_ev	  = FluidProp.Cell_ev; 
	TConsVars   ***DummyIS_cv = FluidProp.DummyIS_cv; 
	TConsVars   ***DummyIP_cv = FluidProp.DummyIP_cv; 
	TPrimVars   ***DummyIS_pv = FluidProp.DummyIS_pv; 
	TPrimVars   ***DummyIP_pv = FluidProp.DummyIP_pv; 
	TConsVars   ***DummyJS_cv = FluidProp.DummyJS_cv; 
	TConsVars   ***DummyJP_cv = FluidProp.DummyJP_cv; 
	TPrimVars   ***DummyJS_pv = FluidProp.DummyJS_pv; 
	TPrimVars   ***DummyJP_pv = FluidProp.DummyJP_pv; 
	TConsVars   ***DummyKS_cv = FluidProp.DummyKS_cv; 
	TConsVars   ***DummyKP_cv = FluidProp.DummyKP_cv; 
	TPrimVars   ***DummyKS_pv = FluidProp.DummyKS_pv; 
	TPrimVars   ***DummyKP_pv = FluidProp.DummyKP_pv; 

	REAL  Rgas   = FluidProp.Rgas;
	REAL  Gamma  = FluidProp.Gamma;

	TNode Vin, Vout;
	REAL  x, y, Vwall_x, Vwall_y;

	if(is==ie)
	{
		for(k=ks;k<=ke;k++)
		{		
			for(j=js;j<=je;j++)
			{
				if(is==1)
				{
					i =1;
					iv=0;
					x = FaceI[i][j][k].coord.x;
					y = FaceI[i][j][k].coord.y;
				}
				else if(is==IM)
				{
					i =IM-1;
					iv=IM;
					x = FaceI[iv][j][k].coord.x;
					y = FaceI[iv][j][k].coord.y;
				}			
				//iv Dummy Cell
				Cell_pv[iv][j][k].press  = Cell_pv[i][j][k].press;
				Cell_pv[iv][j][k].temp   = Cell_pv[i][j][k].temp;
				Cell_pv[iv][j][k].dens   = Cell_pv[i][j][k].dens;

				Vwall_x = -Wall.rot*y;		//轮毂、叶片以Wall.rot旋转，机匣不动
				Vwall_y = Wall.rot*x;

				Vin.x = Cell_pv[i][j][k].uvel - Vwall_x;
				Vin.y = Cell_pv[i][j][k].vvel - Vwall_y;
				Vin.z = Cell_pv[i][j][k].wvel - 0.0f;

				Vout.x = 0.0f - Vin.x;
				Vout.y = 0.0f - Vin.y;
				Vout.z = 0.0f - Vin.z;
				
				Cell_pv[iv][j][k].uvel = Vout.x + Vwall_x;
				Cell_pv[iv][j][k].vvel = Vout.y + Vwall_y;
				Cell_pv[iv][j][k].wvel = Vout.z + 0.0f;    

				Cell_cv[iv][j][k].dens = Cell_cv[i][j][k].dens;
				Cell_cv[iv][j][k].xmom = Cell_pv[iv][j][k].uvel * Cell_cv[iv][j][k].dens;
				Cell_cv[iv][j][k].ymom = Cell_pv[iv][j][k].vvel * Cell_cv[iv][j][k].dens;
				Cell_cv[iv][j][k].zmom = Cell_pv[iv][j][k].wvel * Cell_cv[iv][j][k].dens;
				Cell_cv[iv][j][k].ener = Cell_cv[i][j][k].ener;
				if (FluidProp.turbulenceType==TurbulenceTypes::SA)
				{
					//Cell_ev[iv][j][k].ev = 0.0f-Cell_ev[i][j][k].ev;
					Cell_ev[iv][j][k].ev = 0.0f;
				}

				if(is==1)
				{
					for(nd=1;nd<=NDummy;nd++)
					{
						DummyIS_pv[nd][j][k] = 0.5*(Cell_pv[iv][j][k]+Cell_pv[1][j][k]);
						DummyIS_cv[nd][j][k] = 0.5*(Cell_cv[iv][j][k]+Cell_cv[1][j][k]);
						//DummyIS_pv[nd][j][k] = Cell_pv[iv][j][k];
						//DummyIS_cv[nd][j][k] = Cell_cv[iv][j][k];
					}
				}
				else if(is==IM)
				{
					for(nd=1;nd<=NDummy;nd++)
					{
						DummyIP_pv[nd][j][k] = 0.5*(Cell_pv[iv][j][k]+Cell_pv[IM-1][j][k]);
						DummyIP_cv[nd][j][k] = 0.5*(Cell_cv[iv][j][k]+Cell_cv[IM-1][j][k]);
						//DummyIP_pv[nd][j][k] = Cell_pv[iv][j][k];
						//DummyIP_cv[nd][j][k] = Cell_cv[iv][j][k];
					}
				}	
			}
		}
	}
	else if(js==je)
	{	
		for(k=ks;k<=ke;k++)
		{		
			for(i=is;i<=ie;i++)		
			{
				if(js==1)
				{
					j=1;
					jv=0;
					x  = FaceJ[i][j][k].coord.x;
					y  = FaceJ[i][j][k].coord.y;
				}
				else if(js==JM)
				{
					j=JM-1;
					jv=JM;
					x  = FaceJ[i][jv][k].coord.x;
					y  = FaceJ[i][jv][k].coord.y;
				}			
				//jv Dummy Cell
				Cell_pv[i][jv][k].press  = Cell_pv[i][j][k].press;
				Cell_pv[i][jv][k].temp   = Cell_pv[i][j][k].temp;
				Cell_pv[i][jv][k].dens   = Cell_pv[i][j][k].dens;

				Vwall_x = -Wall.rot*y;
				Vwall_y = Wall.rot*x;

				Vin.x = Cell_pv[i][j][k].uvel - Vwall_x;
				Vin.y = Cell_pv[i][j][k].vvel - Vwall_y;
				Vin.z = Cell_pv[i][j][k].wvel - 0.0f;

				Vout.x = 0.0f - Vin.x;
				Vout.y = 0.0f - Vin.y;
				Vout.z = 0.0f - Vin.z;

				Cell_pv[i][jv][k].uvel = Vout.x + Vwall_x;
				Cell_pv[i][jv][k].vvel = Vout.y + Vwall_y;
				Cell_pv[i][jv][k].wvel = Vout.z + 0.0f;    

				Cell_cv[i][jv][k].dens = Cell_cv[i][j][k].dens;
				Cell_cv[i][jv][k].xmom = Cell_pv[i][jv][k].uvel * Cell_cv[i][jv][k].dens;
				Cell_cv[i][jv][k].ymom = Cell_pv[i][jv][k].vvel * Cell_cv[i][jv][k].dens;
				Cell_cv[i][jv][k].zmom = Cell_pv[i][jv][k].wvel * Cell_cv[i][jv][k].dens;
				Cell_cv[i][jv][k].ener = Cell_cv[i][j][k].ener;
				if (FluidProp.turbulenceType==TurbulenceTypes::SA)
				{
					//Cell_ev[i][jv][k].ev = 0.0f-Cell_ev[i][j][k].ev;
					Cell_ev[i][jv][k].ev = 0.0f;
				}


				if(js==1)
				{
					for(nd=1;nd<=NDummy;nd++)
					{
						DummyJS_pv[i][nd][k] = 0.5*(Cell_pv[i][jv][k]+Cell_pv[i][1][k]);
						DummyJS_cv[i][nd][k] = 0.5*(Cell_cv[i][jv][k]+Cell_cv[i][1][k]);	
						//DummyJS_pv[i][nd][k] = Cell_pv[i][jv][k];
						//DummyJS_cv[i][nd][k] = Cell_cv[i][jv][k];
					}
				}
				else if(js==JM)
				{
					for(nd=1;nd<=NDummy;nd++)
					{
						DummyJP_pv[i][nd][k] = 0.5*(Cell_pv[i][jv][k]+Cell_pv[i][JM-1][k]);
						DummyJP_cv[i][nd][k] = 0.5*(Cell_cv[i][jv][k]+Cell_cv[i][JM-1][k]);
						//DummyJP_pv[i][nd][k] = Cell_pv[i][jv][k];
						//DummyJP_cv[i][nd][k] = Cell_cv[i][jv][k];
					}
				}	
			}		
		}
	}
	else if(ks==ke)
	{
		for(j=js;j<=je;j++)
		{		
			for(i=is;i<=ie;i++)
			{
				if(ks==1)
				{
					k=1;
					kv=0;
					x  = FaceK[i][j][k].coord.x;
					y  = FaceK[i][j][k].coord.y;
				}
				else if(ks==KM)
				{
					k=KM-1;
					kv=KM;
					x  = FaceK[i][j][kv].coord.x;
					y  = FaceK[i][j][kv].coord.y;
				}			
				//kv Dummy Cell
				Cell_pv[i][j][kv].press  = Cell_pv[i][j][k].press;
				Cell_pv[i][j][kv].temp   = Cell_pv[i][j][k].temp;
				Cell_pv[i][j][kv].dens   = Cell_pv[i][j][k].dens;

				Vwall_x = -Wall.rot*y;
				Vwall_y = Wall.rot*x;

				Vin.x = Cell_pv[i][j][k].uvel - Vwall_x;
				Vin.y = Cell_pv[i][j][k].vvel - Vwall_y;
				Vin.z = Cell_pv[i][j][k].wvel - 0.0f;

				Vout.x = 0.0f - Vin.x;
				Vout.y = 0.0f - Vin.y;
				Vout.z = 0.0f - Vin.z;

				Cell_pv[i][j][kv].uvel = Vout.x + Vwall_x; 
				Cell_pv[i][j][kv].vvel = Vout.y + Vwall_y; 
				Cell_pv[i][j][kv].wvel = Vout.z + 0.0f; 

				Cell_cv[i][j][kv].dens = Cell_cv[i][j][k].dens;
				Cell_cv[i][j][kv].xmom = Cell_pv[i][j][kv].uvel * Cell_cv[i][j][kv].dens;
				Cell_cv[i][j][kv].ymom = Cell_pv[i][j][kv].vvel * Cell_cv[i][j][kv].dens;
				Cell_cv[i][j][kv].zmom = Cell_pv[i][j][kv].wvel * Cell_cv[i][j][kv].dens;
				Cell_cv[i][j][kv].ener = Cell_cv[i][j][k].ener;

				if (FluidProp.turbulenceType==TurbulenceTypes::SA)
				{
					//Cell_ev[i][j][kv].ev = 0.0f-Cell_ev[i][j][k].ev;
					Cell_ev[i][j][kv].ev = 0.0f;
				}
	

				if(ks==1)
				{
					for(nd=1;nd<=NDummy;nd++)
					{
						DummyKS_pv[i][j][nd] = 0.5*(Cell_pv[i][j][kv]+Cell_pv[i][j][1]);
						DummyKS_cv[i][j][nd] = 0.5*(Cell_cv[i][j][kv]+Cell_cv[i][j][1]);
						//DummyKS_pv[i][j][nd] = Cell_pv[i][j][kv];
						//DummyKS_cv[i][j][nd] = Cell_cv[i][j][kv];
					}
				}
				else if(ks==KM)
				{
					for(nd=1;nd<=NDummy;nd++)
					{
						DummyKP_pv[i][j][nd] = 0.5*(Cell_pv[i][j][kv]+Cell_pv[i][j][KM-1]);
						DummyKP_cv[i][j][nd] = 0.5*(Cell_cv[i][j][kv]+Cell_cv[i][j][KM-1]);
						//DummyKP_pv[i][j][nd] = Cell_pv[i][j][kv];
						//DummyKP_cv[i][j][nd] = Cell_cv[i][j][kv];
					}
				}
			}		
		}
	}	
}


void CBndConds::BndCondWall_Slip(const TWall &Wall, const CGeometry &Geometry, CFluidProps & FluidProp)
{
	
}