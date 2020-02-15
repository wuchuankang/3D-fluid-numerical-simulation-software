

#include "defs.h"
#include "Mesh.h"
#include "Geometry.h"
#include <iostream>

using namespace std;

//void CMesh::BndCondMixingPlane(const TMixingPlane &MixingPlane)							 
//{
//	int i, j, KM, nd;
//	int NDummy = CGeometry::NDummy;
//
//	int B_ID_1 = MixingPlane.Block_ID_1;
//	int Ft_1 = MixingPlane.Face_type_1;
//	int is_1 = MixingPlane.is_1; 
//	int ie_1 = MixingPlane.ie_1; 
//	int js_1 = MixingPlane.js_1; 
//	int je_1 = MixingPlane.je_1; 
//	int ks_1 = MixingPlane.ks_1; 
//	int ke_1 = MixingPlane.ke_1; 
//		
//	int B_ID_2 = MixingPlane.Block_ID_2;	
//	int Ft_2 = MixingPlane.Face_type_2;
//	int is_2 = MixingPlane.is_2; 
//	int ie_2 = MixingPlane.ie_2; 
//	int js_2 = MixingPlane.js_2; 
//	int je_2 = MixingPlane.je_2; 
//	int ks_2 = MixingPlane.ks_2; 
//	int ke_2 = MixingPlane.ke_2; 
//
//	REAL *Iu_1, *Iu_2, *Iu_3, *Iu_4, *Iu_5, *c_u, *Ma_u;
//	REAL *Id_1, *Id_2, *Id_3, *Id_4, *Id_5, *c_d, *Ma_d;
//	TNode *un_u, *un_d;
//	TPrimVars *pv_u, *pv_d;
//
//	Iu_1  = new REAL [ie_1-is_1+2];
//	Iu_2  = new REAL [ie_1-is_1+2];
//	Iu_3  = new REAL [ie_1-is_1+2];
//	Iu_4  = new REAL [ie_1-is_1+2];
//	Iu_5  = new REAL [ie_1-is_1+2];
//	Id_1  = new REAL [ie_2-is_2+2];
//	Id_2  = new REAL [ie_2-is_2+2];
//	Id_3  = new REAL [ie_2-is_2+2];
//	Id_4  = new REAL [ie_2-is_2+2];
//	Id_5  = new REAL [ie_2-is_2+2];
//
//	c_u	  = new REAL [ie_1-is_1+2];
//	c_d	  = new REAL [ie_2-is_2+2];
//	Ma_u  = new REAL [ie_1-is_1+2];
//	Ma_d  = new REAL [ie_2-is_2+2];
//
//	un_u  = new TNode [ie_1-is_1+2];
//	un_d  = new TNode [ie_2-is_2+2];
//	pv_u  = new TPrimVars [ie_1-is_1+2];
//	pv_d  = new TPrimVars [ie_1-is_1+2];
//	
//	TNode un, un_sum;
//	REAL temp_s, s_sum, Vn, Ma, Ma_switch, temp_A, temp_B, temp_C, V2;
//	REAL Gamma = CFluidProps::Gamma;
//	REAL Rgas = CFluidProps::Rgas;
//	REAL C1ex_u, C2ex_u, C3ex_u, C4ex_u, C5ex_u;
//	REAL C1ex_d, C2ex_d, C3ex_d, C4ex_d, C5ex_d;
//	//Average upstream
//	for(i=is_1;i<=ie_1;i++)
//	{
//		s_sum=0.0;
//		un_sum.x = 0.0;
//		un_sum.y = 0.0;
//		un_sum.z = 0.0;
//		Iu_1[i]=0.0;
//		Iu_2[i]=0.0;
//		Iu_3[i]=0.0;
//		Iu_4[i]=0.0;
//		Iu_5[i]=0.0;
//		Ma_switch=0.0;
//		KM = Block[B_ID_1].Geometry.KM;
//		for(j=js_1;j<=je_1;j++)
//		{
//			un.x = Block[B_ID_1].Geometry.FaceK[i][j][KM].S.x / Block[B_ID_1].Geometry.FaceK[i][j][KM].s;
//			un.y = Block[B_ID_1].Geometry.FaceK[i][j][KM].S.y / Block[B_ID_1].Geometry.FaceK[i][j][KM].s;
//			un.z = Block[B_ID_1].Geometry.FaceK[i][j][KM].S.z / Block[B_ID_1].Geometry.FaceK[i][j][KM].s;
//			temp_s = Block[B_ID_1].Geometry.FaceK[i][j][KM].s;				  
//			un_sum.x = un_sum.x + un.x*temp_s;
//			un_sum.y = un_sum.y + un.y*temp_s;
//			un_sum.z = un_sum.z + un.z*temp_s;
//			s_sum = s_sum + temp_s;
//			Vn = Block[B_ID_1].FluidProp.Cell_pv[i][j][KM-1].uvel*un.x + Block[B_ID_1].FluidProp.Cell_pv[i][j][KM-1].vvel*un.y + Block[B_ID_1].FluidProp.Cell_pv[i][j][KM-1].wvel*un.z;
//			Iu_1[i]=Iu_1[i] +  Block[B_ID_1].FluidProp.Cell_cv[i][j][KM-1].dens*Vn*temp_s;
//			Iu_2[i]=Iu_2[i] + (Block[B_ID_1].FluidProp.Cell_cv[i][j][KM-1].xmom*Vn + un.x*Block[B_ID_1].FluidProp.Cell_pv[i][j][KM-1].press)*temp_s;
//			Iu_3[i]=Iu_3[i] + (Block[B_ID_1].FluidProp.Cell_cv[i][j][KM-1].ymom*Vn + un.y*Block[B_ID_1].FluidProp.Cell_pv[i][j][KM-1].press)*temp_s;
//			Iu_4[i]=Iu_4[i] + (Block[B_ID_1].FluidProp.Cell_cv[i][j][KM-1].zmom*Vn + un.z*Block[B_ID_1].FluidProp.Cell_pv[i][j][KM-1].press)*temp_s;
//			Iu_5[i]=Iu_5[i] + (Block[B_ID_1].FluidProp.Cell_cv[i][j][KM-1].ener + Block[B_ID_1].FluidProp.Cell_pv[i][j][KM-1].press)*Vn*temp_s;
//			Ma = SQRT(POW(Block[B_ID_1].FluidProp.Cell_pv[i][j][KM-1].uvel, 2.0f) + POW(Block[B_ID_1].FluidProp.Cell_pv[i][j][KM-1].vvel, 2.0f) + POW(Block[B_ID_1].FluidProp.Cell_pv[i][j][KM-1].wvel, 2.0f)) / (Gamma*Rgas*Block[B_ID_1].FluidProp.Cell_pv[i][j][KM-1].temp);
//			Ma_switch = Ma_switch + Ma*temp_s;
//		}
//		Iu_1[i] = Iu_1[i] / s_sum;
//		Iu_2[i] = Iu_2[i] / s_sum;
//		Iu_3[i] = Iu_3[i] / s_sum;
//		Iu_4[i] = Iu_4[i] / s_sum;
//		Iu_5[i] = Iu_5[i] / s_sum;
//		un_u[i] = un_sum / s_sum;
//		Ma_switch = Ma_switch / s_sum;
//				  
//		temp_A = Iu_2[i]*un_u[i].x + Iu_3[i]*un_u[i].y + Iu_4[i]*un_u[i].z;
//		temp_B = 2.0f*Iu_1[i]*Iu_5[i] - (Iu_2[i]*Iu_2[i]+Iu_3[i]*Iu_3[i]+Iu_4[i]*Iu_4[i]);
//		temp_C = un_u[i].x*un_u[i].x + un_u[i].y*un_u[i].y + un_u[i].z*un_u[i].z;
//		pv_u[i].press = (temp_A/temp_C + SQRT(temp_A*temp_A/temp_C/temp_C - (Gamma*Gamma-1.0f)*temp_B/temp_C)) / (Gamma+1.0f);		
//		pv_u[i].uvel = (Iu_2[i] - pv_u[i].press*un_u[i].x) / Iu_1[i];
//		pv_u[i].vvel = (Iu_3[i] - pv_u[i].press*un_u[i].y) / Iu_1[i];
//		pv_u[i].wvel = (Iu_4[i] - pv_u[i].press*un_u[i].z) / Iu_1[i];
//		pv_u[i].dens = Iu_1[i] / (pv_u[i].uvel*un_u[i].x+pv_u[i].vvel*un_u[i].y+pv_u[i].wvel*un_u[i].z);
//		pv_u[i].temp = pv_u[i].press/Rgas/pv_u[i].dens;
//		c_u[i]  = SQRT(Gamma*Rgas*pv_u[i].temp);
//		Ma_u[i] = SQRT(pv_u[i].uvel*pv_u[i].uvel+pv_u[i].vvel*pv_u[i].vvel+pv_u[i].wvel*pv_u[i].wvel) / c_u[i];
//		if(_isnan(pv_u[i].dens))
//		{
//			cout<<"Mixing plane upstream average density broken in block:"<<B_ID_1<<" I="<<i<<" Iu_1="<<Iu_1[i]<<endl;
//			system("pause");
//		}
//		if(pv_u[i].dens<0.0)
//		{
//			cout<<"Mixing plane upstream average density negetive in block:"<<B_ID_1<<" I="<<i<<" Iu_1="<<Iu_1[i]<<endl;
//			system("pause");
//		}
//	}
//		
//	//Average downstream
//	for(i=is_2;i<=ie_2;i++)
//	{
//		s_sum=0.0;
//		un_sum.x = 0.0;
//		un_sum.y = 0.0;
//		un_sum.z = 0.0;
//		Id_1[i]=0.0;
//		Id_2[i]=0.0;
//		Id_3[i]=0.0;
//		Id_4[i]=0.0;
//		Id_5[i]=0.0;
//		Ma_switch=0.0;
//		for(j=js_2;j<=je_2;j++)
//		{
//			un.x = Block[B_ID_2].Geometry.FaceK[i][j][1].S.x / Block[B_ID_2].Geometry.FaceK[i][j][1].s;
//			un.y = Block[B_ID_2].Geometry.FaceK[i][j][1].S.y / Block[B_ID_2].Geometry.FaceK[i][j][1].s;
//			un.z = Block[B_ID_2].Geometry.FaceK[i][j][1].S.z / Block[B_ID_2].Geometry.FaceK[i][j][1].s;
//			temp_s = Block[B_ID_2].Geometry.FaceK[i][j][1].s;
//			un_sum.x = un_sum.x + un.x*temp_s;
//			un_sum.y = un_sum.y + un.y*temp_s;
//			un_sum.z = un_sum.z + un.z*temp_s;
//			s_sum = s_sum + temp_s;
//			Vn = Block[B_ID_2].FluidProp.Cell_pv[i][j][1].uvel*un.x + Block[B_ID_2].FluidProp.Cell_pv[i][j][1].vvel*un.y + Block[B_ID_2].FluidProp.Cell_pv[i][j][1].wvel*un.z;
//			Id_1[i]=Id_1[i] +  Block[B_ID_2].FluidProp.Cell_cv[i][j][1].dens*Vn*temp_s;
//			Id_2[i]=Id_2[i] + (Block[B_ID_2].FluidProp.Cell_cv[i][j][1].xmom*Vn + un.x*Block[B_ID_2].FluidProp.Cell_pv[i][j][1].press)*temp_s;
//			Id_3[i]=Id_3[i] + (Block[B_ID_2].FluidProp.Cell_cv[i][j][1].ymom*Vn + un.y*Block[B_ID_2].FluidProp.Cell_pv[i][j][1].press)*temp_s;
//			Id_4[i]=Id_4[i] + (Block[B_ID_2].FluidProp.Cell_cv[i][j][1].zmom*Vn + un.z*Block[B_ID_2].FluidProp.Cell_pv[i][j][1].press)*temp_s;
//			Id_5[i]=Id_5[i] + (Block[B_ID_2].FluidProp.Cell_cv[i][j][1].ener + Block[B_ID_2].FluidProp.Cell_pv[i][j][1].press)*Vn*temp_s;
//		}
//
//	    Id_1[i] = Id_1[i] / s_sum;
//		Id_2[i] = Id_2[i] / s_sum;
//		Id_3[i] = Id_3[i] / s_sum;
//		Id_4[i] = Id_4[i] / s_sum;
//		Id_5[i] = Id_5[i] / s_sum;
//		un_d[i] = un_sum / s_sum;
//		
//		temp_A = Id_2[i]*un_d[i].x + Id_3[i]*un_d[i].y + Id_4[i]*un_d[i].z;
//		temp_B = 2.0f*Id_1[i]*Id_5[i] - (Id_2[i]*Id_2[i]+Id_3[i]*Id_3[i]+Id_4[i]*Id_4[i]);
//		temp_C = un_d[i].x*un_d[i].x + un_d[i].y*un_d[i].y + un_d[i].z*un_d[i].z;
//		pv_d[i].press = (temp_A/temp_C + SQRT(temp_A*temp_A/temp_C/temp_C - (Gamma*Gamma-1.0f)*temp_B/temp_C)) / (Gamma+1.0f);		
//		pv_d[i].uvel  = (Id_2[i]-pv_d[i].press*un_d[i].x) / Id_1[i];
//		pv_d[i].vvel  = (Id_3[i]-pv_d[i].press*un_d[i].y) / Id_1[i];
//		pv_d[i].wvel  = (Id_4[i]-pv_d[i].press*un_d[i].z) / Id_1[i];
//		pv_d[i].dens  = Id_1[i] / (pv_d[i].uvel*un_d[i].x+pv_d[i].vvel*un_d[i].y+pv_d[i].wvel*un_d[i].z);
//		pv_d[i].temp  = pv_d[i].press/Rgas/pv_d[i].dens;
//		c_d[i]   = SQRT(Gamma*Rgas*pv_d[i].temp);
//		Ma_d[i]  = SQRT(pv_d[i].uvel*pv_d[i].uvel+pv_d[i].vvel*pv_d[i].vvel+pv_d[i].wvel*pv_d[i].wvel) / c_d[i];
//				
//		if(_isnan(pv_d[i].press))
//		{
//			cout<<"Mixing plane downstream pressure broken in block:"<<B_ID_2<<" I="<<i<<" K="<<1<<endl;
//			system("pause");
//		}
//		if(_isnan(pv_d[i].press))
//		{
//			cout<<"Mixing plane downstream average density broken in block:"<<B_ID_2<<" I="<<i<<" Id_1="<<Id_1[i]<<endl;
//			system("pause");
//		}
//		if(pv_d[i].dens<0.0)
//		{
//			cout<<"Mixing plane downstream average density negetive in block:"<<B_ID_2<<" I="<<i<<" Id_1="<<Id_1[i]<<endl;
//			system("pause");
//		}
//	}
//
//	//Mixing upstream
//	for(i=is_1;i<=ie_1;i++)
//	{
//		for(j=js_1;j<=je_1;j++)
//		{			
//			KM = Block[B_ID_1].Geometry.KM;
//			C1ex_u = -c_u[i]*c_u[i]*(Block[B_ID_1].FluidProp.Cell_pv[i][j][ke_1-1].dens-pv_u[i].dens) + (Block[B_ID_1].FluidProp.Cell_pv[i][j][ke_1-1].press-pv_u[i].press);
//			C2ex_u = pv_u[i].dens*c_u[i]*(Block[B_ID_1].FluidProp.Cell_pv[i][j][ke_1-1].wvel-pv_u[i].wvel)  + (Block[B_ID_1].FluidProp.Cell_pv[i][j][ke_1-1].press-pv_u[i].press);
//			C3ex_u = pv_u[i].dens*c_u[i]*(Block[B_ID_1].FluidProp.Cell_pv[i][j][ke_1-1].vvel-pv_u[i].vvel);
//			C4ex_u = pv_u[i].dens*c_u[i]*(Block[B_ID_1].FluidProp.Cell_pv[i][j][ke_1-1].uvel-pv_u[i].uvel);
//			C5ex_u = 0.0;
//						
//			Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].dens = pv_d[i].dens + (-C1ex_u + (C2ex_u+C5ex_u)/2.0f) / (c_d[i]*c_d[i]);
//			Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].wvel = pv_d[i].wvel + (C2ex_u-C5ex_u) / (2.0f*pv_d[i].dens*c_d[i]);
//			Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].vvel = Block[B_ID_1].FluidProp.Cell_pv[i][j][KM-1].vvel;
//			Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].uvel = Block[B_ID_1].FluidProp.Cell_pv[i][j][KM-1].uvel;
//			Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].press = pv_d[i].press + (C2ex_u+C5ex_u)/2.0f;
//
//			Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].temp  = Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].press/Rgas/Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].dens;
//			
//			V2 = POW(Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].uvel, 2.0f)+ POW(Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].vvel, 2.0f) + POW(Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].wvel, 2.0f);
//
//			Block[B_ID_1].FluidProp.Cell_cv[i][j][KM].dens = Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].dens;
//			Block[B_ID_1].FluidProp.Cell_cv[i][j][KM].xmom = Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].dens * Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].uvel;
//			Block[B_ID_1].FluidProp.Cell_cv[i][j][KM].ymom = Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].dens * Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].vvel;
//			Block[B_ID_1].FluidProp.Cell_cv[i][j][KM].zmom = Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].dens * Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].wvel;
//			Block[B_ID_1].FluidProp.Cell_cv[i][j][KM].ener = Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].press/(Gamma-1.0f) + Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].dens*V2/2.0f;							
//
//			for(nd=1;nd<=NDummy;nd++)
//			{
//				Block[B_ID_1].FluidProp.DummyKP_pv[i][j][nd] = Block[B_ID_1].FluidProp.Cell_pv[i][j][KM];
//				Block[B_ID_1].FluidProp.DummyKP_cv[i][j][nd] = Block[B_ID_1].FluidProp.Cell_cv[i][j][KM];
//				// evÔÝÊ±Ã»ÓÐ¿¼ÂÇ
//			}
//		}
//	}
//
//	//Mixing downstream
//	for(i=is_2;i<=ie_2;i++)
//	{	
//		for(j=js_2;j<=je_2;j++)
//		{		
//			C1ex_d=0.0;
//			C2ex_d=0.0;
//			C3ex_d=0.0;
//			C4ex_d=0.0;
//			C5ex_d=(0.0f-pv_d[i].dens*c_d[i]*(Block[B_ID_2].FluidProp.Cell_pv[i][j][1].wvel-pv_d[i].wvel)+(Block[B_ID_2].FluidProp.Cell_pv[i][j][1].press-pv_d[i].press));
//						
//			Block[B_ID_2].FluidProp.Cell_pv[i][j][0].dens  = pv_u[i].dens + (-C1ex_d + (C2ex_d+C5ex_d)/2.0f)/(c_u[i]*c_u[i]);
//			Block[B_ID_2].FluidProp.Cell_pv[i][j][0].wvel  = pv_u[i].wvel - C5ex_d/(2.0f*pv_u[i].dens*c_u[i]);
//			Block[B_ID_2].FluidProp.Cell_pv[i][j][0].vvel  = pv_u[i].vvel + C3ex_d/(pv_u[i].dens*c_u[i]);
//			Block[B_ID_2].FluidProp.Cell_pv[i][j][0].uvel  = pv_u[i].uvel + C4ex_d/(pv_u[i].dens*c_u[i]);
//			Block[B_ID_2].FluidProp.Cell_pv[i][j][0].press = pv_u[i].press + (C2ex_d+C5ex_d)/2.0f;
//
//			Block[B_ID_2].FluidProp.Cell_pv[i][j][0].temp  = Block[B_ID_2].FluidProp.Cell_pv[i][j][0].press/Rgas/Block[B_ID_2].FluidProp.Cell_pv[i][j][0].dens;
//			
//			V2 = POW(Block[B_ID_2].FluidProp.Cell_pv[i][j][0].uvel, 2.0f)+ POW(Block[B_ID_2].FluidProp.Cell_pv[i][j][0].vvel, 2.0f) + POW(Block[B_ID_2].FluidProp.Cell_pv[i][j][0].wvel, 2.0f);
//
//			Block[B_ID_2].FluidProp.Cell_cv[i][j][0].dens = Block[B_ID_2].FluidProp.Cell_pv[i][j][0].dens; 
//			Block[B_ID_2].FluidProp.Cell_cv[i][j][0].xmom = Block[B_ID_2].FluidProp.Cell_pv[i][j][0].dens * Block[B_ID_2].FluidProp.Cell_pv[i][j][0].uvel; 
//			Block[B_ID_2].FluidProp.Cell_cv[i][j][0].ymom = Block[B_ID_2].FluidProp.Cell_pv[i][j][0].dens * Block[B_ID_2].FluidProp.Cell_pv[i][j][0].vvel; 
//			Block[B_ID_2].FluidProp.Cell_cv[i][j][0].zmom = Block[B_ID_2].FluidProp.Cell_pv[i][j][0].dens * Block[B_ID_2].FluidProp.Cell_pv[i][j][0].wvel; 
//			Block[B_ID_2].FluidProp.Cell_cv[i][j][0].ener = Block[B_ID_2].FluidProp.Cell_pv[i][j][0].press/(Gamma-1.0f) + Block[B_ID_2].FluidProp.Cell_pv[i][j][0].dens*V2/2.0f;	
//			for(nd=1;nd<=NDummy;nd++)
//			{
//				Block[B_ID_2].FluidProp.DummyKS_pv[i][j][nd] = Block[B_ID_2].FluidProp.Cell_pv[i][j][0];
//				Block[B_ID_2].FluidProp.DummyKS_cv[i][j][nd] = Block[B_ID_2].FluidProp.Cell_cv[i][j][0];
//			}
//		}
//	}
//	delete [] Iu_1;
//	delete [] Iu_2;
//	delete [] Iu_3;
//	delete [] Iu_4;
//	delete [] Iu_5;
//	delete [] Id_1;
//	delete [] Id_2;
//	delete [] Id_3;
//	delete [] Id_4;
//	delete [] Id_5;
//
//	delete [] c_u;
//	delete [] c_d;
//	delete [] Ma_u;
//	delete [] Ma_d;
//	delete [] un_u;
//	delete [] un_d;
//	delete [] pv_u;
//	delete [] pv_d;
//}

void CMesh::BndCondMixingPlane(const TMixingPlane &MixingPlane)							 
{
	int i, j, KM, nd;
	int NDummy = CGeometry::NDummy;

	int B_ID_1 = MixingPlane.Block_ID_1;
	int Ft_1 = MixingPlane.Face_type_1;
	int is_1 = MixingPlane.is_1; 
	int ie_1 = MixingPlane.ie_1; 
	int js_1 = MixingPlane.js_1; 
	int je_1 = MixingPlane.je_1; 
	int ks_1 = MixingPlane.ks_1; 
	int ke_1 = MixingPlane.ke_1; 

	int B_ID_2 = MixingPlane.Block_ID_2;	
	int Ft_2 = MixingPlane.Face_type_2;
	int is_2 = MixingPlane.is_2; 
	int ie_2 = MixingPlane.ie_2; 
	int js_2 = MixingPlane.js_2; 
	int je_2 = MixingPlane.je_2; 
	int ks_2 = MixingPlane.ks_2; 
	int ke_2 = MixingPlane.ke_2; 

	REAL Gamma = CFluidProps::Gamma;
	REAL Rgas = CFluidProps::Rgas;
	REAL Iu[7], Id[7], c_u, c_d, Ma_u, Ma_d, ev_u, ev_d;
	TNode n_u, n_d, n, n_sum;
	TPrimVars pv_u, pv_d;
	REAL temp_s, temp_A, temp_B, temp_C, s_sum, Vn, Ma, Ma_switch, V2;
	REAL C1ex_u, C2ex_u, C3ex_u, C4ex_u, C5ex_u;
	REAL C1ex_d, C2ex_d, C3ex_d, C4ex_d, C5ex_d;

	for(i=is_1;i<=ie_1;i++)
	{	
		//Average upstream
		s_sum=0.0;
		n_sum.x = 0.0;
		n_sum.y = 0.0;
		n_sum.z = 0.0;
		Iu[1] = 0.0;
		Iu[2] = 0.0;
		Iu[3] = 0.0;
		Iu[4] = 0.0;
		Iu[5] = 0.0;
		Iu[6] = 0.0;
		Ma_switch = 0.0;
		KM = Block[B_ID_1].Geometry.KM;
		for(j=js_1;j<=je_1;j++)
		{
			n.x = Block[B_ID_1].Geometry.FaceK[i][j][KM].S.x / Block[B_ID_1].Geometry.FaceK[i][j][KM].s;
			n.y = Block[B_ID_1].Geometry.FaceK[i][j][KM].S.y / Block[B_ID_1].Geometry.FaceK[i][j][KM].s;
			n.z = Block[B_ID_1].Geometry.FaceK[i][j][KM].S.z / Block[B_ID_1].Geometry.FaceK[i][j][KM].s;
			temp_s = Block[B_ID_1].Geometry.FaceK[i][j][KM].s;				  
			n_sum.x = n_sum.x + n.x*temp_s;
			n_sum.y = n_sum.y + n.y*temp_s;
			n_sum.z = n_sum.z + n.z*temp_s;
			s_sum = s_sum + temp_s;
			Vn = Block[B_ID_1].FluidProp.Cell_pv[i][j][KM-1].uvel*n.x + Block[B_ID_1].FluidProp.Cell_pv[i][j][KM-1].vvel*n.y + Block[B_ID_1].FluidProp.Cell_pv[i][j][KM-1].wvel*n.z;
			Iu[1] = Iu[1] +  Block[B_ID_1].FluidProp.Cell_cv[i][j][KM-1].dens*Vn*temp_s;
			Iu[2] = Iu[2] + (Block[B_ID_1].FluidProp.Cell_cv[i][j][KM-1].xmom*Vn + n.x*Block[B_ID_1].FluidProp.Cell_pv[i][j][KM-1].press)*temp_s;
			Iu[3] = Iu[3] + (Block[B_ID_1].FluidProp.Cell_cv[i][j][KM-1].ymom*Vn + n.y*Block[B_ID_1].FluidProp.Cell_pv[i][j][KM-1].press)*temp_s;
			Iu[4] = Iu[4] + (Block[B_ID_1].FluidProp.Cell_cv[i][j][KM-1].zmom*Vn + n.z*Block[B_ID_1].FluidProp.Cell_pv[i][j][KM-1].press)*temp_s;
			Iu[5] = Iu[5] + (Block[B_ID_1].FluidProp.Cell_cv[i][j][KM-1].ener + Block[B_ID_1].FluidProp.Cell_pv[i][j][KM-1].press)*Vn*temp_s;
			if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
				Iu[6] = Iu[6] + Block[B_ID_1].FluidProp.Cell_cv[i][j][KM-1].dens*Block[B_ID_1].FluidProp.Cell_ev[i][j][KM-1].ev*Vn*temp_s;
			V2 = POW(Block[B_ID_1].FluidProp.Cell_pv[i][j][KM-1].uvel, 2.0f) + POW(Block[B_ID_1].FluidProp.Cell_pv[i][j][KM-1].vvel, 2.0f) + POW(Block[B_ID_1].FluidProp.Cell_pv[i][j][KM-1].wvel, 2.0f);
			Ma = SQRT(V2) / (Gamma*Rgas*Block[B_ID_1].FluidProp.Cell_pv[i][j][KM-1].temp);
			Ma_switch = Ma_switch + Ma*temp_s;
		}
		Iu[1] = Iu[1] / s_sum;
		Iu[2] = Iu[2] / s_sum;
		Iu[3] = Iu[3] / s_sum;
		Iu[4] = Iu[4] / s_sum;
		Iu[5] = Iu[5] / s_sum;
		if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
		{
			Iu[6] = Iu[6] / s_sum;
			ev_u = Iu[6] / Iu[1];
		}
		n_u = n_sum / s_sum;
		Ma_switch = Ma_switch / s_sum;

		temp_A = Iu[2]*n_u.x + Iu[3]*n_u.y + Iu[4]*n_u.z;
		temp_B = 2.0f*Iu[1]*Iu[5] - (Iu[2]*Iu[2]+Iu[3]*Iu[3]+Iu[4]*Iu[4]);
		temp_C = n_u.x*n_u.x + n_u.y*n_u.y + n_u.z*n_u.z;
		pv_u.press = (temp_A/temp_C + SQRT(temp_A*temp_A/temp_C/temp_C - (Gamma*Gamma-1.0f)*temp_B/temp_C)) / (Gamma+1.0f);		
		pv_u.uvel = (Iu[2] - pv_u.press*n_u.x) / Iu[1];
		pv_u.vvel = (Iu[3] - pv_u.press*n_u.y) / Iu[1];
		pv_u.wvel = (Iu[4] - pv_u.press*n_u.z) / Iu[1];
		pv_u.dens = Iu[1] / (pv_u.uvel*n_u.x + pv_u.vvel*n_u.y + pv_u.wvel*n_u.z);
		pv_u.temp = pv_u.press/Rgas/pv_u.dens;	
		c_u  = SQRT(Gamma*Rgas*pv_u.temp);
		Ma_u = SQRT(pv_u.uvel*pv_u.uvel+pv_u.vvel*pv_u.vvel+pv_u.wvel*pv_u.wvel) / c_u;
		if(_isnan(pv_u.dens))
		{
			cout<<"Mixing plane upstream average density broken in block:"<<B_ID_1<<" I="<<i<<" Iu_1="<<Iu[1]<<endl;
			system("pause");
		}
		if(pv_u.dens<0.0)
		{
			cout<<"Mixing plane upstream average density negetive in block:"<<B_ID_1<<" I="<<i<<" Iu_1="<<Iu[1]<<endl;
			system("pause");
		}


		//Average downstream
		s_sum=0.0;
		n_sum.x = 0.0;
		n_sum.y = 0.0;
		n_sum.z = 0.0;
		Id[1] = 0.0;
		Id[2] = 0.0;
		Id[3] = 0.0;
		Id[4] = 0.0;
		Id[5] = 0.0;
		Id[6] = 0.0;
		Ma_switch = 0.0;
		for(j=js_2;j<=je_2;j++)
		{
			n.x = Block[B_ID_2].Geometry.FaceK[i][j][1].S.x / Block[B_ID_2].Geometry.FaceK[i][j][1].s;
			n.y = Block[B_ID_2].Geometry.FaceK[i][j][1].S.y / Block[B_ID_2].Geometry.FaceK[i][j][1].s;
			n.z = Block[B_ID_2].Geometry.FaceK[i][j][1].S.z / Block[B_ID_2].Geometry.FaceK[i][j][1].s;
			temp_s = Block[B_ID_2].Geometry.FaceK[i][j][1].s;
			n_sum.x = n_sum.x + n.x*temp_s;
			n_sum.y = n_sum.y + n.y*temp_s;
			n_sum.z = n_sum.z + n.z*temp_s;
			s_sum = s_sum + temp_s;
			Vn = Block[B_ID_2].FluidProp.Cell_pv[i][j][1].uvel*n.x + Block[B_ID_2].FluidProp.Cell_pv[i][j][1].vvel*n.y + Block[B_ID_2].FluidProp.Cell_pv[i][j][1].wvel*n.z;
			Id[1] = Id[1] +  Block[B_ID_2].FluidProp.Cell_cv[i][j][1].dens*Vn*temp_s;
			Id[2] = Id[2] + (Block[B_ID_2].FluidProp.Cell_cv[i][j][1].xmom*Vn + n.x*Block[B_ID_2].FluidProp.Cell_pv[i][j][1].press)*temp_s;
			Id[3] = Id[3] + (Block[B_ID_2].FluidProp.Cell_cv[i][j][1].ymom*Vn + n.y*Block[B_ID_2].FluidProp.Cell_pv[i][j][1].press)*temp_s;
			Id[4] = Id[4] + (Block[B_ID_2].FluidProp.Cell_cv[i][j][1].zmom*Vn + n.z*Block[B_ID_2].FluidProp.Cell_pv[i][j][1].press)*temp_s;
			Id[5] = Id[5] + (Block[B_ID_2].FluidProp.Cell_cv[i][j][1].ener + Block[B_ID_2].FluidProp.Cell_pv[i][j][1].press)*Vn*temp_s;
			if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
				Id[6] = Id[6] + Block[B_ID_2].FluidProp.Cell_cv[i][j][1].dens*Block[B_ID_2].FluidProp.Cell_ev[i][j][1].ev*Vn*temp_s;
		}

		Id[1] = Id[1] / s_sum;
		Id[2] = Id[2] / s_sum;
		Id[3] = Id[3] / s_sum;
		Id[4] = Id[4] / s_sum;
		Id[5] = Id[5] / s_sum;
		if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
		{
			Id[6] = Id[6] / s_sum;
			ev_d = Id[6] / Id[1];
		}
		n_d = n_sum / s_sum;

		temp_A = Id[2]*n_d.x + Id[3]*n_d.y + Id[4]*n_d.z;
		temp_B = 2.0f*Id[1]*Id[5] - (Id[2]*Id[2]+Id[3]*Id[3]+Id[4]*Id[4]);
		temp_C = n_d.x*n_d.x + n_d.y*n_d.y + n_d.z*n_d.z;
		pv_d.press = (temp_A/temp_C + SQRT(temp_A*temp_A/temp_C/temp_C - (Gamma*Gamma-1.0f)*temp_B/temp_C)) / (Gamma+1.0f);		
		pv_d.uvel  = (Id[2]-pv_d.press*n_d.x) / Id[1];
		pv_d.vvel  = (Id[3]-pv_d.press*n_d.y) / Id[1];
		pv_d.wvel  = (Id[4]-pv_d.press*n_d.z) / Id[1];
		pv_d.dens  = Id[1] / (pv_d.uvel*n_d.x+pv_d.vvel*n_d.y+pv_d.wvel*n_d.z);
		pv_d.temp  = pv_d.press/Rgas/pv_d.dens;
		c_d   = SQRT(Gamma*Rgas*pv_d.temp);
		Ma_d  = SQRT(pv_d.uvel*pv_d.uvel+pv_d.vvel*pv_d.vvel+pv_d.wvel*pv_d.wvel) / c_d;

		if(_isnan(pv_d.press))
		{
			cout<<"Mixing plane downstream pressure broken in block:"<<B_ID_2<<" I="<<i<<" K="<<1<<endl;
			system("pause");
		}
		if(pv_d.dens<0.0)
		{
			cout<<"Mixing plane downstream average density negetive in block:"<<B_ID_2<<" I="<<i<<endl;
			system("pause");
		}

		//Mixing upstream
		for(j=js_1;j<=je_1;j++)
		{			
			KM = Block[B_ID_1].Geometry.KM;
			C1ex_u = -c_u*c_u*(Block[B_ID_1].FluidProp.Cell_pv[i][j][ke_1-1].dens-pv_u.dens) + (Block[B_ID_1].FluidProp.Cell_pv[i][j][ke_1-1].press-pv_u.press);
			C2ex_u = pv_u.dens*c_u*(Block[B_ID_1].FluidProp.Cell_pv[i][j][ke_1-1].wvel-pv_u.wvel)  + (Block[B_ID_1].FluidProp.Cell_pv[i][j][ke_1-1].press-pv_u.press);
			C3ex_u = pv_u.dens*c_u*(Block[B_ID_1].FluidProp.Cell_pv[i][j][ke_1-1].vvel-pv_u.vvel);
			C4ex_u = pv_u.dens*c_u*(Block[B_ID_1].FluidProp.Cell_pv[i][j][ke_1-1].uvel-pv_u.uvel);
			C5ex_u = 0.0;

			Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].dens = pv_d.dens + (-C1ex_u + (C2ex_u+C5ex_u)/2.0f) / (c_d*c_d);
			Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].wvel = pv_d.wvel + (C2ex_u-C5ex_u) / (2.0f*pv_d.dens*c_d);
			Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].vvel = Block[B_ID_1].FluidProp.Cell_pv[i][j][KM-1].vvel;
			Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].uvel = Block[B_ID_1].FluidProp.Cell_pv[i][j][KM-1].uvel;
			Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].press = pv_d.press + (C2ex_u+C5ex_u)/2.0f;
			Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].temp  = Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].press/Rgas/Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].dens;

			V2 = POW(Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].uvel, 2.0f)+ POW(Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].vvel, 2.0f) + POW(Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].wvel, 2.0f);
			Block[B_ID_1].FluidProp.Cell_cv[i][j][KM].dens = Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].dens;
			Block[B_ID_1].FluidProp.Cell_cv[i][j][KM].xmom = Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].dens * Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].uvel;
			Block[B_ID_1].FluidProp.Cell_cv[i][j][KM].ymom = Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].dens * Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].vvel;
			Block[B_ID_1].FluidProp.Cell_cv[i][j][KM].zmom = Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].dens * Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].wvel;
			Block[B_ID_1].FluidProp.Cell_cv[i][j][KM].ener = Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].press/(Gamma-1.0f) + Block[B_ID_1].FluidProp.Cell_pv[i][j][KM].dens*V2/2.0f;							
			if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
				Block[B_ID_1].FluidProp.Cell_ev[i][j][KM].ev = ev_u;
			for(nd=1;nd<=NDummy;nd++)
			{
				Block[B_ID_1].FluidProp.DummyKP_pv[i][j][nd] = Block[B_ID_1].FluidProp.Cell_pv[i][j][KM];
				Block[B_ID_1].FluidProp.DummyKP_cv[i][j][nd] = Block[B_ID_1].FluidProp.Cell_cv[i][j][KM];
			}
		}

		//Mixing downstream
		for(j=js_2;j<=je_2;j++)
		{		
			C1ex_d=0.0;
			C2ex_d=0.0;
			C3ex_d=0.0;
			C4ex_d=0.0;
			C5ex_d=(0.0f-pv_d.dens*c_d*(Block[B_ID_2].FluidProp.Cell_pv[i][j][1].wvel-pv_d.wvel)+(Block[B_ID_2].FluidProp.Cell_pv[i][j][1].press-pv_d.press));

			Block[B_ID_2].FluidProp.Cell_pv[i][j][0].dens  = pv_u.dens + (-C1ex_d + (C2ex_d+C5ex_d)/2.0f)/(c_u*c_u);
			Block[B_ID_2].FluidProp.Cell_pv[i][j][0].wvel  = pv_u.wvel - C5ex_d/(2.0f*pv_u.dens*c_u);
			Block[B_ID_2].FluidProp.Cell_pv[i][j][0].vvel  = pv_u.vvel + C3ex_d/(pv_u.dens*c_u);
			Block[B_ID_2].FluidProp.Cell_pv[i][j][0].uvel  = pv_u.uvel + C4ex_d/(pv_u.dens*c_u);
			Block[B_ID_2].FluidProp.Cell_pv[i][j][0].press = pv_u.press + (C2ex_d+C5ex_d)/2.0f;
			Block[B_ID_2].FluidProp.Cell_pv[i][j][0].temp  = Block[B_ID_2].FluidProp.Cell_pv[i][j][0].press/Rgas/Block[B_ID_2].FluidProp.Cell_pv[i][j][0].dens;

			V2 = POW(Block[B_ID_2].FluidProp.Cell_pv[i][j][0].uvel, 2.0f)+ POW(Block[B_ID_2].FluidProp.Cell_pv[i][j][0].vvel, 2.0f) + POW(Block[B_ID_2].FluidProp.Cell_pv[i][j][0].wvel, 2.0f);
			Block[B_ID_2].FluidProp.Cell_cv[i][j][0].dens = Block[B_ID_2].FluidProp.Cell_pv[i][j][0].dens; 
			Block[B_ID_2].FluidProp.Cell_cv[i][j][0].xmom = Block[B_ID_2].FluidProp.Cell_pv[i][j][0].dens * Block[B_ID_2].FluidProp.Cell_pv[i][j][0].uvel; 
			Block[B_ID_2].FluidProp.Cell_cv[i][j][0].ymom = Block[B_ID_2].FluidProp.Cell_pv[i][j][0].dens * Block[B_ID_2].FluidProp.Cell_pv[i][j][0].vvel; 
			Block[B_ID_2].FluidProp.Cell_cv[i][j][0].zmom = Block[B_ID_2].FluidProp.Cell_pv[i][j][0].dens * Block[B_ID_2].FluidProp.Cell_pv[i][j][0].wvel; 
			Block[B_ID_2].FluidProp.Cell_cv[i][j][0].ener = Block[B_ID_2].FluidProp.Cell_pv[i][j][0].press/(Gamma-1.0f) + Block[B_ID_2].FluidProp.Cell_pv[i][j][0].dens*V2/2.0f;
			if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
				Block[B_ID_2].FluidProp.Cell_ev[i][j][0].ev = ev_u;
			for(nd=1;nd<=NDummy;nd++)
			{
				Block[B_ID_2].FluidProp.DummyKS_pv[i][j][nd] = Block[B_ID_2].FluidProp.Cell_pv[i][j][0];
				Block[B_ID_2].FluidProp.DummyKS_cv[i][j][nd] = Block[B_ID_2].FluidProp.Cell_cv[i][j][0];
			}
		}
	}
}

void CMesh::LinearFit(int n_point_base, int n_point_out, REAL* x_ori, REAL* y_ori, REAL* x_res, REAL* y_res)
{
	const int MAX = 250;
	REAL x_base[MAX], y_base[MAX];
	REAL x_out[MAX], y_out[MAX];
	REAL a_lt[MAX], b_lt[MAX];

	int i, j;
	for(i = 1; i <= n_point_base; i++)
	{
		x_base[i]=x_ori[i];
		y_base[i]=y_ori[i];
	}
	
    for(i = 1; i <= n_point_base-1; i++)	
	{
		a_lt[i]=(y_base[i]-y_base[i+1])/(x_base[i]-x_base[i+1]);
		b_lt[i]=y_base[i]-a_lt[i]*x_base[i];
	}

	for(j=1;j<=n_point_out;j++)
	{
		x_out[j]=x_res[j];
	}

	for(j=1;j<=n_point_out;j++)
	{
		for(i = 1; i <= n_point_base-1; i++)
		{
			if(x_out[j]<=x_base[i+1] && x_out[j]>=x_base[i])
			{
				y_out[j]=a_lt[i]*x_out[j]+b_lt[i];
				break;
			}
		}
	}

	for(j=1;j<=n_point_out;j++)
	{
		y_res[j]=y_out[j];		
	}
}