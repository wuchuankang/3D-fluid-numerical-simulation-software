

#include "BndConds.h"
#include "defs.h"
#include "SpaceDiscr.h"

void CBndConds::BndCondInlet(const TInOutlet &Inlet, const CGeometry &Geometry, CFluidProps & FluidProp)
{
	int i, j, k, kv, nd;

	int is = Inlet.is;
	int ie = Inlet.ie;
	int js = Inlet.js;
	int je = Inlet.je;
	int ks = Inlet.ks;
	int ke = Inlet.ke;
	int KM = Geometry.KM;
	int NDummy = Geometry.NDummy;

	TConsVars    ***Cell_cv	   = FluidProp.Cell_cv;  
	TPrimVars    ***Cell_pv    = FluidProp.Cell_pv; 
	TNode		 S;
	REAL  Rgas   = FluidProp.Rgas;
	REAL  Gamma  = FluidProp.Gamma;
	REAL  w_rot  = FluidProp.w_rot;

	REAL  TempBC, PressBC, DensBC, VBC;
	REAL  uvel, vvel, wvel, ener;
	REAL  meduA, meduB, Vn_m, V2_d, R_minus, c_d, c_b, c_02, costh, n_m, vdn;



	Vn_m = SQRT(InletAirAngle.x*InletAirAngle.x + InletAirAngle.y*InletAirAngle.y + InletAirAngle.z*InletAirAngle.z);
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
					S	= Geometry.FaceK[i][j][k].S;
					n_m = SQRT(S.x*S.x + S.y*S.y + S.z*S.z);
	
					uvel = Cell_pv[i][j][k].uvel;	
					vvel = Cell_pv[i][j][k].vvel; 
					wvel = Cell_pv[i][j][k].wvel;					
					vdn	= -1.0f*(uvel*S.x/n_m + vvel*S.y/n_m + wvel*S.z/n_m);	// normal vector pointe to the out of inlet;
				}
				else if(ks==KM)
				{
					k  = KM-1;
					kv = KM;
					S	= Geometry.FaceK[i][j][kv].S;
					n_m = SQRT(S.x*S.x + S.y*S.y + S.z*S.z);

					uvel = Cell_pv[i][j][k].uvel;
					vvel = Cell_pv[i][j][k].vvel;
					wvel = Cell_pv[i][j][k].wvel;			
					vdn	= (uvel*S.x/n_m + vvel*S.y/n_m + wvel*S.z/n_m);	// normal vector pointe to the out of inlet;
				}

				c_d     = SQRT(Gamma * Rgas * Cell_pv[i][j][k].temp);
				V2_d	= uvel*uvel + vvel*vvel + wvel*wvel;
				costh	= -1.0f*vdn / SQRT(V2_d);
				c_02	= c_d*c_d + (Gamma-1.0f)*V2_d/2.0f;
				R_minus = vdn - 2.0f*c_d/(Gamma-1.0f);
				meduA   = (-R_minus*(Gamma-1.0f)) / ((Gamma-1.0f) *costh*costh + 2.0f);
				meduB   = SQRT(c_02 * ((Gamma-1.0f)*costh*costh+2.0f) / (Gamma-1.0f) / (R_minus*R_minus) - (Gamma-1.0f)/2.0f);
				c_b     = meduA * (1.0f + costh*meduB);

				TempBC  = Ttin * (c_b*c_b / c_02);
				if(TempBC>Ttin)
					TempBC = 0.99f*Ttin;
				PressBC = Ptin * POW(TempBC/Ttin, Gamma/(Gamma-1.0f));
				

				DensBC	= PressBC / Rgas / TempBC;
				VBC      = SQRT(2.0f*Gamma/(Gamma-1) * Rgas * (Ttin-TempBC));		

				Cell_pv[i][j][kv].uvel = VBC * InletAirAngle.x/Vn_m;
				Cell_pv[i][j][kv].vvel = VBC * InletAirAngle.y/Vn_m;
				Cell_pv[i][j][kv].wvel = VBC * InletAirAngle.z/Vn_m;

				uvel = Cell_pv[i][j][kv].uvel;
				vvel = Cell_pv[i][j][kv].vvel;
				wvel = Cell_pv[i][j][kv].wvel;
				V2_d = uvel*uvel + vvel*vvel + wvel*wvel;
				ener = PressBC/DensBC/(Gamma-1.0f) + V2_d/2.0f;

				Cell_pv[i][j][kv].press = PressBC;
				Cell_pv[i][j][kv].temp  = TempBC;
				Cell_pv[i][j][kv].dens  = DensBC;

				Cell_cv[i][j][kv].dens = DensBC;
				Cell_cv[i][j][kv].xmom = DensBC * uvel;
				Cell_cv[i][j][kv].ymom = DensBC * vvel;
				Cell_cv[i][j][kv].zmom = DensBC * wvel;
				Cell_cv[i][j][kv].ener = DensBC * ener;

			
				if (FluidProp.turbulenceType==TurbulenceTypes::SA)
					FluidProp.Cell_ev[i][j][kv].ev = ev_ratio * 1.0f;		// see UeserInput		
					//Cell_ev[i][j][kv].ev  = ev_ratio * FluidProp.SutherLand(TempBC)/DensBC * FluidProp.Re_ref;		// Re

				if(ks==1)
				{
					for(nd=1;nd<=NDummy;nd++)
					{
						FluidProp.DummyKS_cv[i][j][nd] = Cell_cv[i][j][kv];
						FluidProp.DummyKS_pv[i][j][nd] = Cell_pv[i][j][kv];
					}				
				}
				else if(ks==KM)
				{
					for(nd=1;nd<=NDummy;nd++)
					{
						FluidProp.DummyKP_cv[i][j][nd] = Cell_cv[i][j][kv];
						FluidProp.DummyKP_pv[i][j][nd] = Cell_pv[i][j][kv];
					}
				}
			}
		}
	}
}
