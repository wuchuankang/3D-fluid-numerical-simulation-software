


#include "BndConds.h"
#include "defs.h"
#include <iostream>
#include "SpaceDiscr.h"

using namespace std;



void CBndConds::BndCondOutlet_Axial(const TInOutlet &Outlet, const CGeometry &Geometry, CFluidProps & FluidProp)
{
	int i, j, k, kv, nd, iref;
	int IM = Geometry.IM;
	int JM = Geometry.JM;
	int KM = Geometry.KM;
	int NDummy = Geometry.NDummy;
	int ks = Outlet.ks;
	int ke = Outlet.ke;

	TConsVars    ***Cell_cv	   = FluidProp.Cell_cv; 
	TPrimVars    ***Cell_pv    = FluidProp.Cell_pv; 
	TFace	     ***FaceK      = Geometry.FaceK;
	REAL  Rgas   = FluidProp.Rgas;
	REAL  Gamma  = FluidProp.Gamma;
	REAL  w_rot  = FluidProp.w_rot;

	REAL  uvel, vvel, wvel, V2_med;
	REAL  th, cu_temp, S_temp, rn, rn_1, dp_dr;
	REAL  cu[200];

	if(ks==ke)
	{
		if(ke==KM)
		{
			kv = KM;
			k  = KM-1;
		}
		for(i=1;i<=IM-1;i++)
		{
			cu[i]  = 0.0;
			S_temp = 0.0;
			for(j=1;j<=JM-1;j++)
			{
				uvel = Cell_pv[i][j][k].uvel;		// absolute velocity
				vvel = Cell_pv[i][j][k].vvel;

				th   = FaceK[i][j][KM].th;
				cu_temp = -1.0f*uvel*SIN(th) + vvel*COS(th);
				cu[i]   = cu[i]  + cu_temp*FaceK[i][j][KM].s;
				S_temp  = S_temp + FaceK[i][j][KM].s;
			}
			cu[i] = cu[i] / S_temp;
		}
				
		iref = int(0.5*IM);
		for(j=1;j<=JM-1;j++)
		{
			//center
			Cell_pv[iref][j][KM].press = Pout;
			//down
			for(i=iref-1;i>=1;i--)
			{
				rn    = FaceK[i][j][KM].r;
				rn_1  = FaceK[i+1][j][KM].r;
				dp_dr = Cell_cv[i][j][KM].dens*cu[i]*cu[i]/rn;
				Cell_pv[i][j][KM].press = Cell_pv[i+1][j][KM].press - dp_dr*(rn_1-rn);		
			}
			//up
			for(i=iref+1;i<=IM-1;i++)
			{
				rn_1 = FaceK[i][j][k].r;
				rn   = FaceK[i-1][j][k].r;
				dp_dr = Cell_cv[i][j][KM].dens*cu[i]*cu[i]/rn;	
				Cell_pv[i][j][KM].press = Cell_pv[i-1][j][KM].press + dp_dr*(rn_1-rn);							
			}
		}
				
		for(i=1;i<=IM-1;i++)
		{
			for(j=1;j<=JM-1;j++)
			{
				Cell_cv[i][j][KM].dens = Cell_cv[i][j][KM-1].dens;		// interpolate
				Cell_cv[i][j][KM].xmom = Cell_cv[i][j][KM-1].xmom;
				Cell_cv[i][j][KM].ymom = Cell_cv[i][j][KM-1].ymom;
				Cell_cv[i][j][KM].zmom = Cell_cv[i][j][KM-1].zmom;
									
				Cell_pv[i][j][KM].dens = Cell_cv[i][j][KM-1].dens;
				Cell_pv[i][j][KM].uvel = Cell_pv[i][j][KM-1].uvel;
				Cell_pv[i][j][KM].vvel = Cell_pv[i][j][KM-1].vvel;
				Cell_pv[i][j][KM].wvel = Cell_pv[i][j][KM-1].wvel;
				Cell_pv[i][j][KM].temp = Cell_pv[i][j][KM].press/Rgas/Cell_cv[i][j][KM].dens;

				uvel = Cell_pv[i][j][KM].uvel;
				vvel = Cell_pv[i][j][KM].vvel;
				wvel = Cell_pv[i][j][KM].wvel;
				V2_med = uvel*uvel + vvel*vvel + wvel*wvel;
				Cell_cv[i][j][KM].ener = Cell_pv[i][j][KM].press/(Gamma-1.0f) + Cell_pv[i][j][KM].dens*V2_med/2.0f;
								
				if (FluidProp.turbulenceType==TurbulenceTypes::SA)
					FluidProp.Cell_ev[i][j][KM].ev = FluidProp.Cell_ev[i][j][KM-1].ev;

				for(nd=1;nd<=NDummy;nd++)
				{
					FluidProp.DummyKP_cv[i][j][nd] = Cell_cv[i][j][KM];
					FluidProp.DummyKP_pv[i][j][nd] = Cell_pv[i][j][KM];
				}
			}
		}
	}
}


void CBndConds::BndCondOutlet_Centrifugal(const TInOutlet &Outlet, const CGeometry &Geometry, CFluidProps & FluidProp)
{
	int i, j, k, kv, nd;

	int is = Outlet.is;
	int ie = Outlet.ie;
	int js = Outlet.js;
	int je = Outlet.je;
	int ks = Outlet.ks;
	int ke = Outlet.ke;
	int KM = Geometry.KM;
	int NDummy = Geometry.NDummy;

	TConsVars    ***Cell_cv	   = FluidProp.Cell_cv;  
	TPrimVars    ***Cell_pv    = FluidProp.Cell_pv; 
	TEddyVars	 ***Cell_ev	   = FluidProp.Cell_ev; 
	TConsVars    ***DummyKS_cv = FluidProp.DummyKS_cv; 
	TConsVars    ***DummyKP_cv = FluidProp.DummyKP_cv; 
	TPrimVars    ***DummyKS_pv = FluidProp.DummyKS_pv; 
	TPrimVars    ***DummyKP_pv = FluidProp.DummyKP_pv; 
	TFace	     ***FaceK      = Geometry.FaceK;
	TNode		 norm;
	REAL  Rgas   = FluidProp.Rgas;
	REAL  Gamma  = FluidProp.Gamma;
	REAL dens0, c0,	V2_d,  ener;
	REAL  w_rot  = FluidProp.w_rot;

	if(ks==ke)
	{
		for(i=is;i<=ie;i++)
		{
			for(j=js;j<=je;j++)
			{
				if(ks==1)
				{
					k  = 1;
					kv = 0;
					norm = 	FaceK[i][j][k].S / FaceK[i][j][k].s;
				}
				else if(ks==KM)
				{
					k  = KM-1;
					kv = KM;
					norm = 	FaceK[i][j][kv].S / FaceK[i][j][kv].s;
				}
				dens0  = Cell_pv[i][j][k].dens;
				c0     = SQRT(Gamma * Rgas * Cell_pv[i][j][k].temp);
				Cell_pv[i][j][kv].press= Pout;
				Cell_pv[i][j][kv].dens = Cell_pv[i][j][k].dens + (Pout-Cell_pv[i][j][k].press)/c0/c0;
				Cell_pv[i][j][kv].uvel = Cell_pv[i][j][k].uvel + norm.x*(Cell_pv[i][j][k].press-Pout)/dens0/c0;
				Cell_pv[i][j][kv].vvel = Cell_pv[i][j][k].vvel + norm.y*(Cell_pv[i][j][k].press-Pout)/dens0/c0;
				Cell_pv[i][j][kv].wvel = Cell_pv[i][j][k].wvel + norm.z*(Cell_pv[i][j][k].press-Pout)/dens0/c0;
				Cell_pv[i][j][kv].temp = Cell_pv[i][j][kv].press/Cell_pv[i][j][kv].dens/Rgas;


				V2_d = Cell_pv[i][j][kv].uvel*Cell_pv[i][j][kv].uvel + Cell_pv[i][j][kv].vvel*Cell_pv[i][j][kv].vvel + Cell_pv[i][j][kv].wvel*Cell_pv[i][j][kv].wvel;
				ener = Cell_pv[i][j][kv].press/Cell_pv[i][j][kv].dens/(Gamma-1.0f) + V2_d/2.0f;

				Cell_cv[i][j][kv].dens = Cell_pv[i][j][kv].dens;
				Cell_cv[i][j][kv].xmom = Cell_pv[i][j][kv].dens * Cell_pv[i][j][kv].uvel;
				Cell_cv[i][j][kv].ymom = Cell_pv[i][j][kv].dens * Cell_pv[i][j][kv].vvel;
				Cell_cv[i][j][kv].zmom = Cell_pv[i][j][kv].dens * Cell_pv[i][j][kv].wvel;
				Cell_cv[i][j][kv].ener = Cell_pv[i][j][kv].dens * ener;

				if (FluidProp.turbulenceType==TurbulenceTypes::SA)
					Cell_ev[i][j][kv].ev = Cell_ev[i][j][k].ev;
				

				if(ks==1)
				{
					for(nd=1;nd<=NDummy;nd++)
					{
						DummyKS_cv[i][j][nd] = Cell_cv[i][j][kv];
						DummyKS_pv[i][j][nd] = Cell_pv[i][j][kv];
					}				
				}
				else if(ks==KM)
				{
					for(nd=1;nd<=NDummy;nd++)
					{
						DummyKP_cv[i][j][nd] = Cell_cv[i][j][kv];
						DummyKP_pv[i][j][nd] = Cell_pv[i][j][kv];
					}
				}
			}
		}
	}
}

// 判断是否是亚声速或是超声速，不同对待
void CBndConds::BndCondOutlet_General(const TInOutlet &Outlet, const CGeometry &Geometry, CFluidProps & FluidProp)
{
	int i, j, k, kv, nd;

	int is = Outlet.is;
	int ie = Outlet.ie;
	int js = Outlet.js;
	int je = Outlet.je;
	int ks = Outlet.ks;
	int ke = Outlet.ke;
	int KM = Geometry.KM;
	int NDummy = Geometry.NDummy;

	TConsVars    ***Cell_cv	   = FluidProp.Cell_cv;  
	TPrimVars    ***Cell_pv    = FluidProp.Cell_pv; 
	TEddyVars	 ***Cell_ev	   = FluidProp.Cell_ev; 
	TConsVars    ***DummyKS_cv = FluidProp.DummyKS_cv; 
	TConsVars    ***DummyKP_cv = FluidProp.DummyKP_cv; 
	TPrimVars    ***DummyKS_pv = FluidProp.DummyKS_pv; 
	TPrimVars    ***DummyKP_pv = FluidProp.DummyKP_pv; 
	TFace	     ***FaceK      = Geometry.FaceK;
	TNode		 norm;
	REAL  Rgas   = FluidProp.Rgas;
	REAL  Gamma  = FluidProp.Gamma;
	REAL dens0, c0,	V2_d,  ener;
	REAL  Vn,  Ma_n;

	if(ks==ke)
	{
		for(i=is;i<=ie;i++)
		{
			for(j=js;j<=je;j++)
			{
				if(ks==1)
				{
					k  = 1;
					kv = 0;
					norm = 	FaceK[i][j][k].S / FaceK[i][j][k].s;
				}
				else if(ks==KM)
				{
					k  = KM-1;
					kv = KM;
					norm = 	FaceK[i][j][kv].S / FaceK[i][j][kv].s;
				}
				dens0  = Cell_pv[i][j][k].dens;
				c0     = SQRT(Gamma * Rgas * Cell_pv[i][j][k].temp);
				Vn = Cell_pv[i][j][k].uvel*norm.x +  Cell_pv[i][j][k].vvel*norm.y + Cell_pv[i][j][k].wvel*norm.z;
				Ma_n = Vn / c0;
				if(Ma_n<=1.0)
				{
					Cell_pv[i][j][kv].press= Pout;
					Cell_pv[i][j][kv].dens = Cell_pv[i][j][k].dens + (Pout-Cell_pv[i][j][k].press)/c0/c0;
					Cell_pv[i][j][kv].uvel = Cell_pv[i][j][k].uvel + norm.x*(Cell_pv[i][j][k].press-Pout)/dens0/c0;
					Cell_pv[i][j][kv].vvel = Cell_pv[i][j][k].vvel + norm.y*(Cell_pv[i][j][k].press-Pout)/dens0/c0;
					Cell_pv[i][j][kv].wvel = Cell_pv[i][j][k].wvel + norm.z*(Cell_pv[i][j][k].press-Pout)/dens0/c0;
				}
				else  //超声速
				{
					Cell_pv[i][j][kv].press= Cell_pv[i][j][k].press;
					Cell_pv[i][j][kv].dens = Cell_pv[i][j][k].dens;
					Cell_pv[i][j][kv].uvel = Cell_pv[i][j][k].uvel;
					Cell_pv[i][j][kv].vvel = Cell_pv[i][j][k].vvel;
					Cell_pv[i][j][kv].wvel = Cell_pv[i][j][k].wvel;
				}

				Cell_pv[i][j][kv].temp = Cell_pv[i][j][kv].press/Cell_pv[i][j][kv].dens/Rgas;
				V2_d = Cell_pv[i][j][kv].uvel*Cell_pv[i][j][kv].uvel + Cell_pv[i][j][kv].vvel*Cell_pv[i][j][kv].vvel + Cell_pv[i][j][kv].wvel*Cell_pv[i][j][kv].wvel;
				ener = Cell_pv[i][j][kv].press/Cell_pv[i][j][kv].dens/(Gamma-1.0f) + V2_d/2.0f;

				Cell_cv[i][j][kv].dens = Cell_pv[i][j][kv].dens;
				Cell_cv[i][j][kv].xmom = Cell_pv[i][j][kv].dens * Cell_pv[i][j][kv].uvel;
				Cell_cv[i][j][kv].ymom = Cell_pv[i][j][kv].dens * Cell_pv[i][j][kv].vvel;
				Cell_cv[i][j][kv].zmom = Cell_pv[i][j][kv].dens * Cell_pv[i][j][kv].wvel;
				Cell_cv[i][j][kv].ener = Cell_pv[i][j][kv].dens * ener;

				if (FluidProp.turbulenceType==TurbulenceTypes::SA)
					Cell_ev[i][j][kv].ev = Cell_ev[i][j][k].ev;
				

				if(ks==1)
				{
					for(nd=1;nd<=NDummy;nd++)
					{
						DummyKS_cv[i][j][nd] = Cell_cv[i][j][kv];
						DummyKS_pv[i][j][nd] = Cell_pv[i][j][kv];
					}				
				}
				else if(ks==KM)
				{
					for(nd=1;nd<=NDummy;nd++)
					{
						DummyKP_cv[i][j][nd] = Cell_cv[i][j][kv];
						DummyKP_pv[i][j][nd] = Cell_pv[i][j][kv];
					}
				}
			}
		}
	}
}

 // 压力非均匀出口
//void CBndConds::BndCondOutlet_General(const TInOutlet &Outlet, const CGeometry &Geometry, CFluidProps & FluidProp)
//{
//	int i, j, k, kv, nd;
//	int IM = Geometry.IM;
//	int JM = Geometry.JM;
//	int KM = Geometry.KM;
//	int NDummy = Geometry.NDummy;
//	int ks = Outlet.ks;
//	int ke = Outlet.ke;
//
//	TConsVars    ***Cell_cv	   = FluidProp.Cell_cv; 
//	TPrimVars    ***Cell_pv    = FluidProp.Cell_pv; 
//	TFace	     ***FaceK      = Geometry.FaceK;
//	REAL Rgas   = FluidProp.Rgas;
//	REAL Gamma  = FluidProp.Gamma;
//	REAL w_rot  = FluidProp.w_rot;
//	REAL gradP_x, gradP_y, gradP_z;
//	REAL P_x[200][200], P_y[200][200], P_z[200][200], P[200][200];
//	REAL gradPn, dx, dy, dz;
//	REAL Pa=0.0f;
//
//	REAL  uvel, vvel, wvel, V2_med;
//
//	if(ks==ke)
//	{
//		if(ke==KM)
//		{
//			kv = KM;
//			k  = KM-1;
//		}
//				
//		for(i=1;i<=IM-1;i++)
//		{
//			for(j=1;j<=JM-1;j++)
//			{
//				gradP_x = 2.0f*w_rot*Cell_pv[i][j][k].dens*Cell_pv[i][j][k].vvel + Cell_pv[i][j][k].dens*w_rot*w_rot*FaceK[i][j][kv].coord.x;
//				gradP_y = -2.0f*w_rot*Cell_pv[i][j][k].dens*Cell_pv[i][j][k].uvel + Cell_pv[i][j][k].dens*w_rot*w_rot*FaceK[i][j][kv].coord.y;
//				gradP_z = 0.0;
//				gradPn = (gradP_x*FaceK[i][j][kv].S.x + gradP_y*FaceK[i][j][kv].S.y + gradP_z*FaceK[i][j][kv].S.z) / FaceK[i][j][kv].s;
//				//P_x[i][j] = gradP_x - gradPn*FaceK[i][j][kv].S.x/FaceK[i][j][kv].s;
//				//P_y[i][j] = gradP_y - gradPn*FaceK[i][j][kv].S.y/FaceK[i][j][kv].s;
//				//P_z[i][j] = gradP_z - gradPn*FaceK[i][j][kv].S.z/FaceK[i][j][kv].s;li
//								
//				P_x[i][j] = gradPn*FaceK[i][j][kv].S.x/FaceK[i][j][kv].s;
//				P_y[i][j] = gradPn*FaceK[i][j][kv].S.y/FaceK[i][j][kv].s;
//				P_z[i][j] = gradPn*FaceK[i][j][kv].S.z/FaceK[i][j][kv].s;
//
//			}
//		}
//			
//		P[1][1] = Pa;
//		for(i=2;i<=IM-1;i++)
//		{
//			dx = (FaceK[i][1][kv].coord.x-FaceK[i-1][1][kv].coord.x);
//			dy = (FaceK[i][1][kv].coord.y-FaceK[i-1][1][kv].coord.y);
//			dz = (FaceK[i][1][kv].coord.z-FaceK[i-1][1][kv].coord.z);
//			P[i][1] = P[i-1][1] + dx*P_x[i-1][1] + dy*P_y[i-1][1] + dz*P_z[i-1][1];
//		}
//		for(j=2;j<=JM-1;j++)
//		{
//			dx = (FaceK[1][j][kv].coord.x-FaceK[1][j-1][kv].coord.x);
//			dy = (FaceK[1][j][kv].coord.y-FaceK[1][j-1][kv].coord.y);
//			dz = (FaceK[1][j][kv].coord.z-FaceK[1][j-1][kv].coord.z);
//			P[1][j] = P[1][j-1] + dx*P_x[1][j-1] + dy*P_y[1][j-1] + dz*P_z[1][j-1];
//		}
//		for(i=2;i<=IM-1;i++)
//		{
//			for(j=2;j<=JM-1;j++)
//			{
//				dx = (FaceK[i][j][kv].coord.x-FaceK[i-1][j][kv].coord.x)+(FaceK[i][j][kv].coord.x-FaceK[i][j-1][kv].coord.x) + (FaceK[i][j][kv].coord.x-FaceK[i-1][j-1][kv].coord.x);
//				dy = (FaceK[i][j][kv].coord.y-FaceK[i-1][j][kv].coord.y)+(FaceK[i][j][kv].coord.y-FaceK[i][j-1][kv].coord.y) + (FaceK[i][j][kv].coord.y-FaceK[i-1][j-1][kv].coord.y);
//				dz = (FaceK[i][j][kv].coord.z-FaceK[i-1][j][kv].coord.z)+(FaceK[i][j][kv].coord.z-FaceK[i][j-1][kv].coord.z) + (FaceK[i][j][kv].coord.z-FaceK[i-1][j-1][kv].coord.z);
//				P[i][j] = 1.0f/3.0f*(P[i-1][j]+P[i][j-1]+P[i-1][j-1]) - 1.0f/3.0f*P_x[i][j]*dx - 1.0f/3.0f*P_y[i][j]*dy - 1.0f/3.0f*P_z[i][j]*dz;
//			}
//		}
//		for(i=1;i<=IM-1;i++)
//		{
//			for(j=1;j<=JM-1;j++)
//			{
//				Cell_pv[i][j][kv].press = Pout +  P[i][j];
//							
//				Cell_cv[i][j][kv].dens = Cell_cv[i][j][k].dens;		// interpolate
//				Cell_cv[i][j][kv].xmom = Cell_cv[i][j][k].xmom;
//				Cell_cv[i][j][kv].ymom = Cell_cv[i][j][k].ymom;
//				Cell_cv[i][j][kv].zmom = Cell_cv[i][j][k].zmom;
//													  
//				Cell_pv[i][j][kv].dens = Cell_cv[i][j][k].dens;
//				Cell_pv[i][j][kv].uvel = Cell_pv[i][j][k].uvel;
//				Cell_pv[i][j][kv].vvel = Cell_pv[i][j][k].vvel;
//				Cell_pv[i][j][kv].wvel = Cell_pv[i][j][k].wvel;
//				Cell_pv[i][j][kv].temp = Cell_pv[i][j][kv].press/Rgas/Cell_cv[i][j][kv].dens;
//
//				uvel = Cell_pv[i][j][kv].uvel;
//				vvel = Cell_pv[i][j][kv].vvel;
//				wvel = Cell_pv[i][j][kv].wvel;
//				V2_med = uvel*uvel + vvel*vvel + wvel*wvel;
//				Cell_cv[i][j][kv].ener = Cell_pv[i][j][kv].press/(Gamma-1.0f) + Cell_pv[i][j][kv].dens*V2_med/2.0f;
//								
//				if (FluidProp.turbulenceType==TurbulenceTypes::SA)
//					FluidProp.Cell_ev[i][j][kv].ev = FluidProp.Cell_ev[i][j][k].ev;
//
//				for(nd=1;nd<=NDummy;nd++)
//				{
//					FluidProp.DummyKP_cv[i][j][nd] = Cell_cv[i][j][kv];
//					FluidProp.DummyKP_pv[i][j][nd] = Cell_pv[i][j][kv];
//				}
//			}
//		}
//	}
//}