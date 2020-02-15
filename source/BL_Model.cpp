
#include "../include/Blocks.h"
#include <iostream>

using namespace std;


void CBlocks::BL_Model(const CBndGeom &BndGeom)		//李老师的
{
	int i, j, k, nb, i1, j1, k1, i0, j0, k0;
	int  IM = Geometry.IM;
	int  JM = Geometry.JM;
	int  KM = Geometry.KM;
	

	TVisVars  ***Cell_tur  = FluidProp.Cell_tur;
	REAL uvel, vvel, wvel;
	REAL *dens, *Vvel, *Omega, *dis, *mue_l, *mue_t;
	dens  = nullptr;
	Vvel  = nullptr;
	dis	  = nullptr;
	Omega = nullptr;
	mue_l = nullptr;
	mue_t = nullptr;

	int MM = MAX3(IM, JM, KM);
	dens   = new REAL [MM];
	Vvel   = new REAL [MM];
	dis    = new REAL [MM];
	Omega  = new REAL [MM];
	mue_l  = new REAL [MM];
	mue_t  = new REAL [MM];

	for(i=1;i<=IM-1;i++)
	{
		for(j=1;j<=JM-1;j++)
		{
			for(k=1;k<=KM-1;k++)
			{
				flag[i][j][k] = 0;
			}
		}
	}
		
	//I direction turbulence viscous			
	for(nb=1;nb<=BndGeom.nWall;nb++)
	{
		if(this->Block_ID==BndGeom.Wall[nb].Block_ID && (BndGeom.Wall[nb].Face_type==1 || BndGeom.Wall[nb].Face_type==2))
		{
#pragma omp parallel for  num_threads(NT) default(shared) firstprivate(i,j,k, i1, i0, uvel, vvel, wvel, dis, dens, Vvel, Omega, mue_l, mue_t)
			for(k=BndGeom.Wall[nb].ks;k<=BndGeom.Wall[nb].ke;k++)
			{
				for(j=BndGeom.Wall[nb].js;j<=BndGeom.Wall[nb].je;j++)
				{
					for(i=1;i<=IM-1;i++)
					{
						if(BndGeom.Wall[nb].Face_type==1)
						{
							i1 = 1;
							i0 = i;
						}
						else
						{
							i1 = IM;
							i0 = IM-i;
						}

						dis[i0] = SQRT(POW((Geometry.Cell[i][j][k].coord.x - Geometry.FaceI[i1][j][k].coord.x),2.0) + POW((Geometry.Cell[i][j][k].coord.y - Geometry.FaceI[i1][j][k].coord.y),2.0)+
								  POW((Geometry.Cell[i][j][k].coord.z - Geometry.FaceI[i1][j][k].coord.z),2.0));
						dens[i0] = FluidProp.Cell_pv[i][j][k].dens;
						uvel	 = FluidProp.Cell_pv[i][j][k].uvel;
						vvel	 = FluidProp.Cell_pv[i][j][k].vvel;
						wvel	 = FluidProp.Cell_pv[i][j][k].wvel;
						Vvel[i0] = SQRT(uvel*uvel + vvel*vvel + wvel*wvel);
						Omega[i0] = FluidProp.Omega[i][j][k];
						mue_l[i0] = FluidProp.SutherLand(FluidProp.Cell_pv[i][j][k].temp);	
					}

					BL_model_1d(IM, dens, Vvel, dis, Omega, mue_l, mue_t);

					//如果一个点处于多条线上，取粘性系数最小的值
					for(i=1;i<=IM-1;i++)
					{
						if(BndGeom.Wall[nb].Face_type==1)
							i0 = i;
						else
							i0 = IM-i;

						if(flag[i][j][k] ==0)
						{
							flag[i][j][k] = 1;
							Cell_tur[i][j][k].mue = mue_t[i0];
						}
						else
						{
							Cell_tur[i][j][k].mue = MIN(Cell_tur[i][j][k].mue, mue_t[i0]);
						}
					}
				}
			}
			//设置壁面第1层网格上的湍流粘性系数为0,以保证壁面上mut=0 
			for(j=BndGeom.Wall[nb].js;j<=BndGeom.Wall[nb].je;j++)
			{
				for(k=BndGeom.Wall[nb].ks;k<=BndGeom.Wall[nb].ke;k++)
				{
					if(BndGeom.Wall[nb].Face_type==1)
					{
						Cell_tur[1][j][k].mue = 0.0;
						Cell_tur[0][j][k].mue = 0.0;
						//Cell_tur[0][j][k].mue = -Cell_tur[1][j][k].mue;		//此种做法，在网格较为稀疏或是质量不好的情况下，会使得计算出现负密度
					}
					else
					{
						Cell_tur[IM-1][j][k].mue = 0.0;
						Cell_tur[IM][j][k].mue = 0.0;
						//Cell_tur[IM][j][k].mue = -Cell_tur[IM-1][j][k].mue;
					}
				}
			}
		}

		if(this->Block_ID==BndGeom.Wall[nb].Block_ID && (BndGeom.Wall[nb].Face_type==3 || BndGeom.Wall[nb].Face_type==4))
		{
#pragma omp parallel for  num_threads(NT) default(shared) firstprivate(i,j,k, j1, j0, uvel, vvel, wvel, dis, dens, Vvel, Omega, mue_l, mue_t)
			for(k=BndGeom.Wall[nb].ks;k<=BndGeom.Wall[nb].ke;k++)
			{
				for(i=BndGeom.Wall[nb].is;i<=BndGeom.Wall[nb].ie;i++)
				{
					for(j=1;j<=JM-1;j++)
					{
						if(BndGeom.Wall[nb].Face_type==3)
						{
							j1 = 1;
							j0 = j;
						}
						else
						{
							j1 = JM;
							j0 = JM-j;
						}

						dis[j0] = SQRT(POW((Geometry.Cell[i][j][k].coord.x - Geometry.FaceJ[i][j1][k].coord.x),2.0) + POW((Geometry.Cell[i][j][k].coord.y - Geometry.FaceJ[i][j1][k].coord.y),2.0)+
								  POW((Geometry.Cell[i][j][k].coord.z - Geometry.FaceJ[i][j1][k].coord.z),2.0));
						dens[j0] = FluidProp.Cell_pv[i][j][k].dens;
						uvel	 = FluidProp.Cell_pv[i][j][k].uvel;
						vvel	 = FluidProp.Cell_pv[i][j][k].vvel;
						wvel	 = FluidProp.Cell_pv[i][j][k].wvel;
						Vvel[j0] = SQRT(uvel*uvel + vvel*vvel + wvel*wvel);
						Omega[j0] = FluidProp.Omega[i][j][k];
						mue_l[j0] =  FluidProp.SutherLand(FluidProp.Cell_pv[i][j][k].temp);	
					}

					BL_model_1d(JM, dens, Vvel, dis, Omega, mue_l, mue_t);

					//如果一个点处于多条线上，取粘性系数最小的值
					for(j=1;j<=JM-1;j++)
					{
						if(BndGeom.Wall[nb].Face_type==3)
							j0 = j;
						else
							j0 = JM-j;

						if(flag[i][j][k] ==0)
						{
							flag[i][j][k] = 1;
							Cell_tur[i][j][k].mue = mue_t[j0];
						}
						else
						{
							Cell_tur[i][j][k].mue = MIN(Cell_tur[i][j][k].mue, mue_t[j0]);
						}
					}
				}
			}

			//设置壁面第1层网格上的湍流粘性系数为0,以保证壁面上mut=0 
			for(k=BndGeom.Wall[nb].ks;k<=BndGeom.Wall[nb].ke;k++)
			{
				for(i=BndGeom.Wall[nb].is;i<=BndGeom.Wall[nb].ie;i++)
				{
					if(BndGeom.Wall[nb].Face_type==3)
					{
						Cell_tur[i][1][k].mue = 0.0;
						Cell_tur[i][0][k].mue = 0.0;
						//Cell_tur[i][0][k].mue = -Cell_tur[i][1][k].mue;
					}
					else
					{
						Cell_tur[i][JM-1][k].mue = 0.0;
						Cell_tur[i][JM][k].mue = 0.0;
						//Cell_tur[i][JM][k].mue = -Cell_tur[i][JM-1][k].mue;
					}
				}
			}
		}

		if(this->Block_ID==BndGeom.Wall[nb].Block_ID && (BndGeom.Wall[nb].Face_type==5 || BndGeom.Wall[nb].Face_type==6))
		{
#pragma omp parallel for  num_threads(NT) default(shared) firstprivate(i,j,k, k1, k0, uvel, vvel, wvel, dis, dens, Vvel, Omega, mue_l, mue_t)

			for(j=BndGeom.Wall[nb].js;j<=BndGeom.Wall[nb].je;j++)
			{
				for(i=BndGeom.Wall[nb].is;i<=BndGeom.Wall[nb].ie;i++)
				{
					for(k=1;k<=KM-1;k++)
					{
						if(BndGeom.Wall[nb].Face_type==5)
						{
							k1 = 1;
							k0 = k;
						}
						else
						{
							k1 = KM;
							k0 = KM-k;
						}

						dis[k0] = SQRT(POW((Geometry.Cell[i][j][k].coord.x - Geometry.FaceK[i][j][k1].coord.x),2.0) + POW((Geometry.Cell[i][j][k].coord.y - Geometry.FaceK[i][j][k1].coord.y),2.0)+
								  POW((Geometry.Cell[i][j][k].coord.z - Geometry.FaceK[i][j][k1].coord.z),2.0));
						dens[k0] = FluidProp.Cell_pv[i][j][k].dens;
						uvel	 = FluidProp.Cell_pv[i][j][k].uvel;
						vvel	 = FluidProp.Cell_pv[i][j][k].vvel;
						wvel	 = FluidProp.Cell_pv[i][j][k].wvel;
						Vvel[k0] = SQRT(uvel*uvel + vvel*vvel + wvel*wvel);
						Omega[k0] = FluidProp.Omega[i][j][k];
						mue_l[k0] =  FluidProp.SutherLand(FluidProp.Cell_pv[i][j][k].temp);	
					}

					BL_model_1d(JM, dens, Vvel, dis, Omega, mue_l, mue_t);

					//如果一个点处于多条线上，取粘性系数最小的值
					for(k=1;k<=KM-1;k++)
					{
						if(BndGeom.Wall[nb].Face_type==5)
							k0 = k;
						else
							k0 = KM-k;

						if(flag[i][j][k] ==0)
						{
							flag[i][j][k] = 1;
							Cell_tur[i][j][k].mue = mue_t[k0];
						}
						else
						{
							Cell_tur[i][j][k].mue = MIN(Cell_tur[i][j][k].mue, mue_t[k0]);
						}
					}
				}
			}

			//设置壁面第1层网格上的湍流粘性系数为0,以保证壁面上mut=0 
			for(j=BndGeom.Wall[nb].js;j<=BndGeom.Wall[nb].je;j++)
			{
				for(i=BndGeom.Wall[nb].is;i<=BndGeom.Wall[nb].ie;i++)
				{
					if(BndGeom.Wall[nb].Face_type==5)
					{
						Cell_tur[i][j][1].mue = 0.0;
						Cell_tur[i][j][0].mue = 0.0;
						//Cell_tur[i][j][0].mue = -Cell_tur[i][j][1].mue;
					}
					else
					{
						Cell_tur[i][j][KM-1].mue = 0.0;
						Cell_tur[i][j][KM].mue = 0.0;
						//Cell_tur[i][j][KM].mue = -Cell_tur[i][j][KM-1].mue;
					}
				}
			}
		}
	}
	delete[] dens;		dens  = nullptr;
	delete[] Vvel;	    Vvel  = nullptr;
	delete[] dis;	    dis	  = nullptr;
	delete[] Omega;		Omega = nullptr;
	delete[] mue_l;		mue_l = nullptr;
	delete[] mue_t;		mue_t = nullptr;
}



void CBlocks::BL_model_1d(int M, REAL dens[], REAL Vvel[], REAL dis[], REAL Omega[], REAL mue_l[], REAL mue_t[])	//李老师的
{
	int i;

	REAL Tw, Ret, F, Fmax, y_max, Udif, y_plus, bL;
	REAL Fwake, Fkleb, mue_in, mue_out;
	REAL A_plus=26.0f, Cwk=1.0f, k=0.4f, Ckleb=0.3f, alpha=0.0168f, Ccp=1.6f;

	Tw = ABS(Omega[1]*mue_l[1]);
	for(i=1;i<M;i++)
	{
		if(ABS(Omega[i]*mue_l[i]) > Tw)
			Tw = ABS(Omega[i]*mue_l[i]);
	}
	Ret = SQRT(dens[1]*Tw)/mue_l[1];

	Fmax=0.0 ;   y_max=0.0 ;     Udif=0.0;	

	//find the Fmax, y_max
	for(i=1;i<M;i++)
	{
		if(Vvel[i] > Udif)
			Udif = Vvel[i];
		y_plus = dis[i]*Ret;
		F = dis[i]*ABS(Omega[i])*(1.0f-EXP(-y_plus/A_plus));
		if(F > Fmax)
		{
			Fmax = F;
			y_max = dis[i];
		}
	}

	Fwake = MIN(y_max*Fmax, Cwk*y_max*Udif*Udif/Fmax);
	for(i=1;i<M;i++)
	{
		y_plus = Ret*dis[i];
		bL = k*dis[i]*(1.0f-EXP(-y_plus/A_plus));
		mue_in = dens[i]*bL*bL*Omega[i];
		Fkleb = 1.0f/(1.0f + 5.5f*POW((Ckleb*dis[i]/y_max), 6.0));
		mue_out = alpha*Ccp*dens[i]*Fwake*Fkleb;

		if(ABS(mue_in) > ABS(mue_out))
			mue_t[i] = mue_out;
		else
			mue_t[i] = mue_in;
	}
}


//自己修正的并不好
//void CBlocks::BL_Model(const CBndGeom &BndGeom)
//{
//	int i, j, k, nb, i1, j1, k1, i0, j0, k0;
//	int  IM = Geometry.IM;
//	int  JM = Geometry.JM;
//	int  KM = Geometry.KM;
//	
//
//	TVisVars  ***Cell_tur  = FluidProp.Cell_tur;
//	REAL uvel, vvel, wvel, Vvel_wall;
//	REAL *dens, *Vvel, *Omega, *dis, *mue_l, *mue_t;
//	dens  = nullptr;
//	Vvel  = nullptr;
//	dis	  = nullptr;
//	Omega = nullptr;
//	mue_l = nullptr;
//	mue_t = nullptr;
//
//	
//	for(i=1;i<=IM-1;i++)
//	{
//		for(j=1;j<=JM-1;j++)
//		{
//			for(k=1;k<=KM-1;k++)
//			{
//				flag[i][j][k] = 0;
//			}
//		}
//	}
//
//	//I direction turbulence viscous			
//	for(nb=1;nb<=BndGeom.nWall;nb++)
//	{
//		if(this->Block_ID==BndGeom.Wall[nb].Block_ID && (BndGeom.Wall[nb].Face_type==1 || BndGeom.Wall[nb].Face_type==2))
//		{
//			dens   = new REAL [IM];
//			Vvel   = new REAL [IM];
//			dis    = new REAL [IM];
//			Omega  = new REAL [IM];
//			mue_l  = new REAL [IM];
//			mue_t  = new REAL [IM];
//			for(k=BndGeom.Wall[nb].ks;k<=BndGeom.Wall[nb].ke;k++)
//			{
//				for(j=BndGeom.Wall[nb].js;j<=BndGeom.Wall[nb].je;j++)
//				{
//					for(i=1;i<=IM-1;i++)
//					{
//						if(BndGeom.Wall[nb].Face_type==1)
//						{
//							i1 = 1;
//							i0 = i;
//						}
//						else
//						{
//							i1 = IM;
//							i0 = IM-i;
//						}
//
//						dis[i0] = SQRT(POW((Geometry.Cell[i][j][k].coord.x - Geometry.FaceI[i1][j][k].coord.x),2.0) + POW((Geometry.Cell[i][j][k].coord.y - Geometry.FaceI[i1][j][k].coord.y),2.0)+
//								  POW((Geometry.Cell[i][j][k].coord.z - Geometry.FaceI[i1][j][k].coord.z),2.0));
//						dens[i0] = FluidProp.Cell_pv[i][j][k].dens;
//						uvel	 = FluidProp.Cell_pv[i][j][k].uvel;
//						vvel	 = FluidProp.Cell_pv[i][j][k].vvel;
//						wvel	 = FluidProp.Cell_pv[i][j][k].wvel;
//						Vvel[i0] = SQRT(uvel*uvel + vvel*vvel + wvel*wvel);
//						Omega[i0] = FluidProp.Omega[i][j][k];
//						mue_l[i0] = SutherLand(FluidProp.Cell_pv[i][j][k].temp);	
//					}
//					Vvel_wall = SQRT(POW(FluidProp.Cell_pv[i1-1][j][k].uvel+FluidProp.Cell_pv[i1][j][k].uvel, 2.0f) + POW(FluidProp.Cell_pv[i1-1][j][k].vvel+FluidProp.Cell_pv[i1][j][k].vvel, 2.0f) + POW(FluidProp.Cell_pv[i1-1][j][k].wvel+FluidProp.Cell_pv[i1][j][k].wvel, 2.0f));
//					BL_model_1d(IM, dens, Vvel, dis, Omega, mue_l, mue_t, Vvel_wall);
//					
//					//如果一个点处于多条线上，取粘性系数最小的值
//					for(i=1;i<=IM-1;i++)
//					{
//						if(BndGeom.Wall[nb].Face_type==1)
//							i0 = i;
//						else
//							i0 = IM-i;
//
//						if(flag[i][j][k] ==0)
//						{
//							flag[i][j][k] = 1;
//							Cell_tur[i][j][k].mue = mue_t[i0];
//						}
//						else
//						{
//							Cell_tur[i][j][k].mue = MIN(Cell_tur[i][j][k].mue, mue_t[i0]);
//						}
//					}
//				}
//			}
//			//设置壁面第1层网格上的湍流粘性系数为0,以保证壁面上mut=0 
//			for(j=BndGeom.Wall[nb].js;j<=BndGeom.Wall[nb].je;j++)
//			{
//				for(k=BndGeom.Wall[nb].ks;k<=BndGeom.Wall[nb].ke;k++)
//				{
//					if(BndGeom.Wall[nb].Face_type==1)
//					{
//						Cell_tur[1][j][k].mue = 0.0;
//						Cell_tur[0][j][k].mue = 0.0;
//						//Cell_tur[0][j][k].mue = -Cell_tur[1][j][k].mue;		//此种做法，在网格较为稀疏或是质量不好的情况下，会使得计算出现负密度
//					}
//					else
//					{
//						Cell_tur[IM-1][j][k].mue = 0.0;
//						Cell_tur[IM][j][k].mue = 0.0;
//						//Cell_tur[IM][j][k].mue = -Cell_tur[IM-1][j][k].mue;
//					}
//				}
//			}
//			delete[] dens;		dens  = nullptr;
//			delete[] Vvel;	    Vvel  = nullptr;
//			delete[] dis;	    dis	  = nullptr;
//			delete[] Omega;		Omega = nullptr;
//			delete[] mue_l;		mue_l = nullptr;
//			delete[] mue_t;		mue_t = nullptr;
//
//		}
//
//		if(this->Block_ID==BndGeom.Wall[nb].Block_ID && (BndGeom.Wall[nb].Face_type==3 || BndGeom.Wall[nb].Face_type==4))
//		{
//			dens   = new REAL [JM];
//			Vvel   = new REAL [JM];
//			dis    = new REAL [JM];
//			Omega  = new REAL [JM];
//			mue_l  = new REAL [JM];
//			mue_t  = new REAL [JM];
//			for(k=BndGeom.Wall[nb].ks;k<=BndGeom.Wall[nb].ke;k++)
//			{
//				for(i=BndGeom.Wall[nb].is;i<=BndGeom.Wall[nb].ie;i++)
//				{
//					for(j=1;j<=JM-1;j++)
//					{
//						if(BndGeom.Wall[nb].Face_type==3)
//						{
//							j1 = 1;
//							j0 = j;
//						}
//						else
//						{
//							j1 = JM;
//							j0 = JM-j;
//						}
//
//						dis[j0] = SQRT(POW((Geometry.Cell[i][j][k].coord.x - Geometry.FaceJ[i][j1][k].coord.x),2.0) + POW((Geometry.Cell[i][j][k].coord.y - Geometry.FaceJ[i][j1][k].coord.y),2.0)+
//								  POW((Geometry.Cell[i][j][k].coord.z - Geometry.FaceJ[i][j1][k].coord.z),2.0));
//						dens[j0] = FluidProp.Cell_pv[i][j][k].dens;
//						uvel	 = FluidProp.Cell_pv[i][j][k].uvel;
//						vvel	 = FluidProp.Cell_pv[i][j][k].vvel;
//						wvel	 = FluidProp.Cell_pv[i][j][k].wvel;
//						Vvel[j0] = SQRT(uvel*uvel + vvel*vvel + wvel*wvel);
//						Omega[j0] = FluidProp.Omega[i][j][k];
//						mue_l[j0] = SutherLand(FluidProp.Cell_pv[i][j][k].temp);	
//					}
//					Vvel_wall = SQRT(POW(FluidProp.Cell_pv[i][j1-1][k].uvel+FluidProp.Cell_pv[i][j1][k].uvel, 2.0f) + POW(FluidProp.Cell_pv[i][j1-1][k].vvel+FluidProp.Cell_pv[i][j1][k].vvel, 2.0f) + POW(FluidProp.Cell_pv[i][j1-1][k].wvel+FluidProp.Cell_pv[i][j1][k].wvel, 2.0f));
//					BL_model_1d(JM, dens, Vvel, dis, Omega, mue_l, mue_t, Vvel_wall);
//
//					//如果一个点处于多条线上，取粘性系数最小的值
//					for(j=1;j<=JM-1;j++)
//					{
//						if(BndGeom.Wall[nb].Face_type==3)
//							j0 = j;
//						else
//							j0 = JM-j;
//
//						if(flag[i][j][k] ==0)
//						{
//							flag[i][j][k] = 1;
//							Cell_tur[i][j][k].mue = mue_t[j0];
//						}
//						else
//						{
//							Cell_tur[i][j][k].mue = MIN(Cell_tur[i][j][k].mue, mue_t[j0]);
//						}
//					}
//				}
//			}
//
//			//设置壁面第1层网格上的湍流粘性系数为0,以保证壁面上mut=0 
//			for(k=BndGeom.Wall[nb].ks;k<=BndGeom.Wall[nb].ke;k++)
//			{
//				for(i=BndGeom.Wall[nb].is;i<=BndGeom.Wall[nb].ie;i++)
//				{
//					if(BndGeom.Wall[nb].Face_type==3)
//					{
//						Cell_tur[i][1][k].mue = 0.0;
//						Cell_tur[i][0][k].mue = 0.0;
//						//Cell_tur[i][0][k].mue = -Cell_tur[i][1][k].mue;
//					}
//					else
//					{
//						Cell_tur[i][JM-1][k].mue = 0.0;
//						Cell_tur[i][JM][k].mue = 0.0;
//						//Cell_tur[i][JM][k].mue = -Cell_tur[i][JM-1][k].mue;
//					}
//				}
//			}
//			delete[] dens;		dens  = nullptr;
//			delete[] Vvel;	    Vvel  = nullptr;
//			delete[] dis;	    dis	  = nullptr;
//			delete[] Omega;		Omega = nullptr;
//			delete[] mue_l;		mue_l = nullptr;
//			delete[] mue_t;		mue_t = nullptr;
//		}
//
//
//
//		if(this->Block_ID==BndGeom.Wall[nb].Block_ID && (BndGeom.Wall[nb].Face_type==5 || BndGeom.Wall[nb].Face_type==6))
//		{
//			dens   = new REAL [KM];
//			Vvel   = new REAL [KM];
//			dis    = new REAL [KM];
//			Omega  = new REAL [KM];
//			mue_l  = new REAL [KM];
//			mue_t  = new REAL [KM];
//			for(j=BndGeom.Wall[nb].js;j<=BndGeom.Wall[nb].je;j++)
//			{
//				for(i=BndGeom.Wall[nb].is;i<=BndGeom.Wall[nb].ie;i++)
//				{
//					for(k=1;k<=KM-1;k++)
//					{
//						if(BndGeom.Wall[nb].Face_type==5)
//						{
//							k1 = 1;
//							k0 = k;
//						}
//						else
//						{
//							k1 = KM;
//							k0 = KM-k;
//						}
//
//						dis[k0] = SQRT(POW((Geometry.Cell[i][j][k].coord.x - Geometry.FaceK[i][j][k1].coord.x),2.0) + POW((Geometry.Cell[i][j][k].coord.y - Geometry.FaceK[i][j][k1].coord.y),2.0)+
//								  POW((Geometry.Cell[i][j][k].coord.z - Geometry.FaceK[i][j][k1].coord.z),2.0));
//						dens[k0] = FluidProp.Cell_pv[i][j][k].dens;
//						uvel	 = FluidProp.Cell_pv[i][j][k].uvel;
//						vvel	 = FluidProp.Cell_pv[i][j][k].vvel;
//						wvel	 = FluidProp.Cell_pv[i][j][k].wvel;
//						Vvel[k0] = SQRT(uvel*uvel + vvel*vvel + wvel*wvel);
//						Omega[k0] = FluidProp.Omega[i][j][k];
//						mue_l[k0] = SutherLand(FluidProp.Cell_pv[i][j][k].temp);	
//					}
//
//					Vvel_wall = SQRT(POW(FluidProp.Cell_pv[i][j][k1-1].uvel+FluidProp.Cell_pv[i][j][k1].uvel, 2.0f) + POW(FluidProp.Cell_pv[i][j][k1-1].vvel+FluidProp.Cell_pv[i][j][k1].vvel, 2.0f) + POW(FluidProp.Cell_pv[i][j][k1-1].wvel+FluidProp.Cell_pv[i][j][k1].wvel, 2.0f));
//					BL_model_1d(KM, dens, Vvel, dis, Omega, mue_l, mue_t, Vvel_wall);
//
//					//如果一个点处于多条线上，取粘性系数最小的值
//					for(k=1;k<=KM-1;k++)
//					{
//						if(BndGeom.Wall[nb].Face_type==5)
//							k0 = k;
//						else
//							k0 = KM-k;
//
//						if(flag[i][j][k] ==0)
//						{
//							flag[i][j][k] = 1;
//							Cell_tur[i][j][k].mue = mue_t[k0];
//						}
//						else
//						{
//							Cell_tur[i][j][k].mue = MIN(Cell_tur[i][j][k].mue, mue_t[k0]);
//						}
//					}
//				}
//			}
//
//			//设置壁面第1层网格上的湍流粘性系数为0,以保证壁面上mut=0 
//			for(j=BndGeom.Wall[nb].js;j<=BndGeom.Wall[nb].je;j++)
//			{
//				for(i=BndGeom.Wall[nb].is;i<=BndGeom.Wall[nb].ie;i++)
//				{
//					if(BndGeom.Wall[nb].Face_type==5)
//					{
//						Cell_tur[i][j][1].mue = 0.0;
//						Cell_tur[i][j][0].mue = 0.0;
//						//Cell_tur[i][j][0].mue = -Cell_tur[i][j][1].mue;
//					}
//					else
//					{
//						Cell_tur[i][j][KM-1].mue = 0.0;
//						Cell_tur[i][j][KM].mue = 0.0;
//						//Cell_tur[i][j][KM].mue = -Cell_tur[i][j][KM-1].mue;
//					}
//				}
//			}
//
//			delete[] dens;		dens  = nullptr;
//			delete[] Vvel;	    Vvel  = nullptr;
//			delete[] dis;	    dis	  = nullptr;
//			delete[] Omega;		Omega = nullptr;
//			delete[] mue_l;		mue_l = nullptr;
//			delete[] mue_t;		mue_t = nullptr;
//
//		}
//	}
//}
//
//void CBlocks::BL_model_1d(int M, REAL dens[], REAL Vvel[], REAL dis[], REAL Omega[], REAL mue_l[], REAL mue_t[], REAL Vvel_wall)	//自己修正
//{
//	int i;
//
//	REAL Tw, Ret, F, Fmax, y_max, U_max, U_y_max, Udif, y_plus, bL;
//	REAL Fwake, Fkleb, mue_in, mue_out;
//	REAL A_plus=26.0f, Cwk=1.0f, k=0.4f, Ckleb=0.3f, alpha=0.0168f, Ccp=1.6f;
//	REAL  Vvel_grad_wall;
//
//	//Tw = ABS(Omega[1]*mue_l[1]);
//	//for(i=1;i<M;i++)
//	//{
//	//	if(ABS(Omega[i]*mue_l[i]) > Tw)
//	//		Tw = ABS(Omega[i]*mue_l[i]);
//	//}
//	//Ret = SQRT(dens[1]*Tw)/mue_l[1];
//	
//	Vvel_grad_wall = (Vvel[1]-Vvel_wall)/dis[1];
//	Tw = mue_l[1]*Vvel_grad_wall;
//	Ret = SQRT(dens[1]*Tw)/mue_l[1];
//
//
//	Fmax=0.0 ;   y_max=0.0 ;     Udif=0.0;	U_max=0.0;
//
//	//find the Fmax, y_max
//	for(i=1;i<M;i++)
//	{
//		//if(Vvel[i] > Udif)
//		//	Udif = Vvel[i];
//		//y_plus = dis[i]*Ret;
//		//F = dis[i]*ABS(Omega[i])*(1.0f-EXP(-y_plus/A_plus));
//		//if(F > Fmax)
//		//{
//		//	Fmax = F;
//		//	y_max = dis[i];
//		//}
//				
//		if(Vvel[i] > U_max)
//			U_max = Vvel[i];
//		y_plus = dis[i]*Ret;
//		F = dis[i]*ABS(Omega[i])*(1.0f-EXP(-y_plus/A_plus));
//		if(F > Fmax)
//		{
//			Fmax = F;
//			y_max = dis[i];
//			U_y_max = Vvel[i];
//		}
//	}
//	Udif = U_max - U_y_max;
//
//	Fwake = MIN(y_max*Fmax, Cwk*y_max*Udif*Udif/Fmax);
//	for(i=1;i<M;i++)
//	{
//		y_plus = Ret*dis[i];
//		bL = k*dis[i]*(1.0f-EXP(-y_plus/A_plus));
//		mue_in = dens[i]*bL*bL*Omega[i];
//		Fkleb = 1.0f/(1.0f + 5.5f*POW((Ckleb*dis[i]/y_max), 6.0));
//		mue_out = alpha*Ccp*dens[i]*Fwake*Fkleb;
//
//		if(ABS(mue_in) > ABS(mue_out))
//			mue_t[i] = mue_out;
//		else
//			mue_t[i] = mue_in;
//	}
//}
