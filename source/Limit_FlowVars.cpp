

#include "FluidProps.h"
#include <iostream>

using namespace  std;

void CFluidProps::Limit_Vars()
{
	int i, j, k, Nneg = 0;
	int kn, ia, ja, ka, i1, j1, k1;
	int Kneg[4][11];
	REAL dens, uvel, vvel, wvel, press;
	REAL sn1;

	LdensMin  = 1.0E-6f;
	LdensMax  = 1000.0f;
	LpressMin = 1.0E-6f;
	LpressMax = 1000.0f;
	LVvelMax  = 1000.0f;

	for(i=1;i<=IM-1;i++)
	{
		for(j=1;j<=JM-1;j++)
		{
			for(k=1;k<=KM-1;k++)
			{
				if(Cell_pv[i][j][k].dens<LdensMin || Cell_pv[i][j][k].dens>LdensMax || Cell_pv[i][j][k].press<LpressMin || Cell_pv[i][j][k].press>LpressMax ||
				   Cell_pv[i][j][k].uvel>LVvelMax || Cell_pv[i][j][k].vvel>LVvelMax || Cell_pv[i][j][k].wvel>LVvelMax)
				{
					Nneg = Nneg + 1;
					if(Nneg <=10)			// The location of Limited point
					{
						Kneg[1][Nneg] = i;
						Kneg[2][Nneg] = j;
						Kneg[3][Nneg] = k;
					}
					kn = 0;
					dens = 0.0; uvel = 0.0; vvel = 0.0; wvel = 0.0; press = 0.0;

					for(ia=-1;ia<=1;ia++)
					{
						for(ja=-1;ja<=1;ja++)
						{
							for(ka=-1;ka<=1;ka++)
							{
								k1 = k + ka; 									; 
								j1 = j + ja; 
								i1 = i + ia;
								if(i1>0 && i1<IM && j1>0 && j1<JM && k1>0 && k1<KM)
								{
									if(false==(Cell_pv[i1][j1][k1].dens<LdensMin || Cell_pv[i1][j1][k1].dens>LdensMax || Cell_pv[i1][j1][k1].press<LpressMin || Cell_pv[i1][j1][k1].press>LpressMax ||
											   Cell_pv[i1][j1][k1].uvel>LVvelMax || Cell_pv[i1][j1][k1].vvel>LVvelMax || Cell_pv[i1][j1][k1].wvel>LVvelMax))
									{
										dens  = dens + Cell_pv[i1][j1][k1].dens;
										uvel  = uvel + Cell_pv[i1][j1][k1].uvel;
										vvel  = vvel + Cell_pv[i1][j1][k1].vvel;
										wvel  = wvel + Cell_pv[i1][j1][k1].wvel;
										press = press + Cell_pv[i1][j1][k1].press;
										kn = kn + 1;
									}
								}
							}
						}
					}
					if(kn ==0)		// 周围全部为“坏点”
					{
						dens  = Cell_pv[i1][j1][k1].dens;
						uvel  = Cell_pv[i1][j1][k1].uvel;
						vvel  = Cell_pv[i1][j1][k1].vvel;
						wvel  = Cell_pv[i1][j1][k1].wvel;
						press = Cell_pv[i1][j1][k1].press;

						if(dens < LdensMin)		 dens = LdensMin;
						if(dens > LdensMax)		 dens = LdensMax;
						if(press < LpressMin)	 press = LpressMin;
						if(press > LpressMax)	 press = LpressMax;
						if(ABS(uvel) > LVvelMax) uvel = sign(LVvelMax, uvel);
						if(ABS(vvel) > LVvelMax) vvel = sign(LVvelMax, vvel);
						if(ABS(wvel) > LVvelMax) wvel = sign(LVvelMax, wvel);
					}
					else
					{
						sn1 = 1.0f / kn;
						dens  = dens * sn1;  
						uvel  = uvel * sn1;  
						vvel  = vvel * sn1;  
						wvel  = wvel * sn1;  
						press = press * sn1;				
					}
					Cell_pv[i][j][k].dens  = dens;
					Cell_pv[i][j][k].uvel  = uvel ;
					Cell_pv[i][j][k].vvel  = vvel ;
					Cell_pv[i][j][k].wvel  = wvel ;
					Cell_pv[i][j][k].press = press;
					Cell_pv[i][j][k].temp = press/dens/Rgas;

					Cell_cv[i][j][k].dens = dens;
					Cell_cv[i][j][k].xmom = dens * uvel;
					Cell_cv[i][j][k].ymom = dens * vvel;
					Cell_cv[i][j][k].zmom = dens * wvel;
					Cell_cv[i][j][k].ener = press/(Gamma-1.0f) + dens*(uvel*uvel + vvel*vvel + wvel*wvel)/2.0f;
				}
			}
		}
	}

	if(Nneg > 0)
	{
		cout<<"Warning::limit flow variables..."<<endl;
		cout<<"LdensMin = "<<LdensMin<<"  "<<"LdensMax = "<<LdensMax<<"  "<<"LpressMin = "<<LpressMin<<"  "<<"LpressMax = "<<LpressMax<<"  "<<"LVvelMax = "<<LVvelMax<<endl;
		for(i=1;i<=MIN(Nneg, 10);i++)
		{
			i1=Kneg[1][i]; j1=Kneg[2][i]; k1=Kneg[3][i];
			cout<<"i="<<i1<<"  "<<"j="<<j1<<"  "<<"k="<<k1<<endl;
		}
		Flag_low_dt = 10;		// 今后10步内，本块降低时间步
	}
}


void CFluidProps::Limit_SA()
{
	int i, j, k, kn, ia, ja, ka, i1, j1, k1;
	REAL SAmin = 1.0E-8f, LSAmax=1000.0f;
	REAL s0;

	for(i=1;i<=IM-1;i++)
	{
		for(j=1;j<=JM-1;j++)
		{
			for(k=1;k<=KM-1;k++)
			{
				if(Cell_ev[i][j][k].ev < SAmin)
					Cell_ev[i][j][k].ev = SAmin;
				if(Cell_ev[i][j][k].ev > LSAmax)
				{
					s0 = 0.0;
					kn = 0;
					for(ka=-1;ka<=1;ka++)
					{
						for(ja=-1;ja<=1;ja++)
						{
							for(ia=-1;ia<=1;ia++)
							{
								k1 = k + ka; 									; 
								j1 = j + ja; 
								i1 = i + ia;
								if(i1>0 && i1<IM && j1>0 && j1<JM && k1>0 && k1<KM)
								{
									if((Cell_ev[i1][j1][k1].ev>=SAmin) && (Cell_ev[i1][j1][k1].ev<=LSAmax))
									{
										s0  = s0 + Cell_ev[i1][j1][k1].ev;
										kn = kn + 1;
									}
								}
							}
						}
					}
					if(kn ==0)								 //周围全部为“坏点”
						Cell_ev[i][j][k].ev = LSAmax;			 
					else
						Cell_ev[i][j][k].ev = s0/kn;
				}
			}
		}
	}
}
  