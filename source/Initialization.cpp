// initialize the flow variables of every Block

/*---------------------------------------------------
初始化方法：
方法1. 线性插值：压力p， 气流角，  需要人为给定：压比， 出口气流角
		此方法的函数是在整个区域上的，所以在solver.h中声明，定义

方法2. 均匀值

----------------------------------------------------*/

#include "Mesh.h"
#include "Blocks.h"
#include "defs.h"
#include <fstream>

void CBlocks::Initial_Uniform()
{
	int  i, j, k;

	TNode InletAirAngle = BndCond.InletAirAngle;
	REAL  Ma0	   = FluidProp.Ma0;
	REAL  Gamma	   = FluidProp.Gamma;	//γ specific heat ratio
	REAL  Rgas	   = FluidProp.Rgas;	
	REAL  Ptin     = BndCond.Ptin;
	REAL  Ttin     = BndCond.Ttin;
	REAL  Pout     = BndCond.Pout;
	REAL  ev_ratio = BndCond.ev_ratio;


	REAL  press, temp, dens, Vel, n_m, uvel, vvel, wvel, ener, V2, ev_initial;
		  

	press = Ptin * POW(1.0f + (Gamma-1.0f)/2.0f *Ma0*Ma0, 0.0f-Gamma/(Gamma-1.0f));
	temp  = Ttin * POW(1.0f + (Gamma-1.0f)/2.0f *Ma0*Ma0, 0.0f-1.0f);
	dens  = press/(Rgas*temp);
	Vel	  = POW(Gamma*Rgas*temp, 0.5f) * Ma0;
	n_m   = POW(InletAirAngle.x*InletAirAngle.x + InletAirAngle.y*InletAirAngle.y + InletAirAngle.z*InletAirAngle.z,0.5);	

	ev_initial = ev_ratio * 1.0f;		// see UserInput
	//ev_initial = ev_ratio * FluidProp.SutherLand(temp)/dens * FluidProp.Re_ref;		//Re!!!

#pragma omp parallel for num_threads(NT) default(shared) private(i, j, k, uvel, vvel, wvel, V2, ener)
	for(i=1;i<=Geometry.IM-1;i++)
	{
		for(j=1;j<=Geometry.JM-1;j++)
		{
			for(k=1;k<=Geometry.KM-1;k++)
			{
				// axial
				uvel  = Vel * InletAirAngle.x / n_m;
				vvel  = Vel * InletAirAngle.y / n_m;
				wvel  = Vel * InletAirAngle.z / n_m;
			
				FluidProp.Cell_pv[i][j][k].uvel  = uvel;
				FluidProp.Cell_pv[i][j][k].vvel  = vvel;
				FluidProp.Cell_pv[i][j][k].wvel  = wvel;
				FluidProp.Cell_pv[i][j][k].dens  = dens;
				FluidProp.Cell_pv[i][j][k].temp  = temp;
				FluidProp.Cell_pv[i][j][k].press = press;
				
				FluidProp.Cell_cv[i][j][k].dens = dens;
				FluidProp.Cell_cv[i][j][k].xmom = dens * uvel;
				FluidProp.Cell_cv[i][j][k].ymom = dens * vvel;
				FluidProp.Cell_cv[i][j][k].zmom = dens * wvel;
				V2	  = uvel*uvel + vvel*vvel + wvel*wvel;
				ener  = press/(dens * (Gamma-1.0f)) + V2/2.0f;
				FluidProp.Cell_cv[i][j][k].ener = dens * ener;
				
				if (FluidProp.turbulenceType==TurbulenceTypes::SA)
					FluidProp.Cell_ev[i][j][k].ev = ev_initial;		
			}			
		}	
	}	

}

void CMesh::Initial_Interpolate(CBlocks Block[])
{
	int  i, j, k, nb;
	int Block_ID_in  =  BndGeom.Inlet[1].Block_ID;		//只能有一个进口和出口,
	int Block_ID_out =  BndGeom.Outlet[1].Block_ID;

	TNode InletAirAngle  = CBndConds::InletAirAngle;
	TNode OutletAirAngle = CBndConds::OutletAirAngle;
	TNode dBeta, Beta;
	REAL  Ptin     = CBndConds::Ptin;
	REAL  Ttin     = CBndConds::Ttin;
	REAL  Pout     = CBndConds::Pout;
	REAL  P21ratio = CBndConds::P21ratio;
	REAL  ev_ratio = CBndConds::ev_ratio;
	REAL  Gamma	   = CFluidProps::Gamma;	
	REAL  Rgas	   = CFluidProps::Rgas;	

	REAL Pin, press, temp, dens, Vel, n_m, uvel, vvel, wvel, cs, ener, Ma, V2, ev_initial;	
	REAL zmin, zmax, dz, dp;


	zmin = 1E+32f;
	zmax = -1E+32f;			//假定进口坐标小于出口坐标
	for(i=1;i<=Block[Block_ID_in].Geometry.IM-1;i++)
	{
		for(j=1;j<=Block[Block_ID_in].Geometry.JM-1;j++)
		{
			zmin = MIN(Block[Block_ID_in].Geometry.FaceK[i][j][BndGeom.Inlet[1].ks].coord.z, zmin);		//限定了只能是K面
		}
	}
	for(i=1;i<=Block[Block_ID_out].Geometry.IM-1;i++)
	{
		for(j=1;j<=Block[Block_ID_out].Geometry.JM-1;j++)
		{
			zmax = MAX(Block[Block_ID_out].Geometry.FaceK[i][j][BndGeom.Outlet[1].ke].coord.z, zmax);
		}
	}
	dz = zmax - zmin;

	Pin = Pout / P21ratio;
	if(Pin>=Ptin)
		Pin = 0.99f*Ptin;
	dp = Pout - Pin;

    dBeta = OutletAirAngle - InletAirAngle;

    temp = Ttin * POW((Pin/Ptin), (Gamma-1.0f)/Gamma);
	dens = Pin / (Rgas*temp);
    cs   = SQRT(Gamma*Pin/dens);
    Ma   = SQRT(2.0f*((Ttin/temp)-1.0f)/(Gamma-1.0f));
    Vel  = Ma*cs;

	ev_initial = ev_ratio * 1.0f;
	//ev_initial = ev_ratio * SutherLand(temp)/dens;
		
	for(nb=1;nb<=NB;nb++)
	{
#pragma omp parallel for num_threads(NT) default(shared) private(i, j, k, Beta, press, dens, n_m, uvel, vvel, wvel, V2, ener)
		for(i=1;i<=Block[nb].Geometry.IM-1;i++)
		{
			for(j=1;j<=Block[nb].Geometry.JM-1;j++)
			{
				for(k=1;k<=Block[nb].Geometry.KM-1;k++)
				{
					Beta  = InletAirAngle + dBeta*(Block[nb].Geometry.Cell[i][j][k].coord.z-zmin)/dz;
					press = Pin + dp*(Block[nb].Geometry.Cell[i][j][k].coord.z-zmin)/dz;
					dens  = press/Rgas/temp;
					n_m = SQRT(Beta.x*Beta.x + Beta.y*Beta.y + Beta.z*Beta.z);

					uvel  = Vel * Beta.x / n_m;
					vvel  = Vel * Beta.y / n_m;
					wvel  = Vel * Beta.z / n_m;
					Block[nb].FluidProp.Cell_pv[i][j][k].uvel  = uvel;
					Block[nb].FluidProp.Cell_pv[i][j][k].vvel  = vvel;
					Block[nb].FluidProp.Cell_pv[i][j][k].wvel  = wvel;
					Block[nb].FluidProp.Cell_pv[i][j][k].dens  = dens;
					Block[nb].FluidProp.Cell_pv[i][j][k].temp  = temp;
					Block[nb].FluidProp.Cell_pv[i][j][k].press = press;

					Block[nb].FluidProp.Cell_cv[i][j][k].dens = dens;
					Block[nb].FluidProp.Cell_cv[i][j][k].xmom = dens * uvel;
					Block[nb].FluidProp.Cell_cv[i][j][k].ymom = dens * vvel;
					Block[nb].FluidProp.Cell_cv[i][j][k].zmom = dens * wvel;
					V2	  = uvel*uvel + vvel*vvel + wvel*wvel;
					ener  = press/(Gamma-1.0f) + dens*V2/2.0f;
					Block[nb].FluidProp.Cell_cv[i][j][k].ener = ener;

					if (Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
						Block[nb].FluidProp.Cell_ev[i][j][k].ev = ev_initial;		
				}			
			}	
		}	
	}
}