


#include "Mesh.h"
#include "FluidProps.h"

using namespace std;

int l1_max, m1_max, l2_max, m2_max;
TConsVars  ShadowCell_A_cv[300][300], ShadowCell_B_cv[300][300];
TPrimVars  ShadowCell_A_pv[300][300], ShadowCell_B_pv[300][300];
TEddyVars  ShadowCell_A_ev[300][300], ShadowCell_B_ev[300][300];
		

void CMesh::ShadowOut(CBlocks &A, int fc, int nd, int is, int ie, int js, int je, int ks, int ke)
{
	int i, j, k, l, m;
	switch(fc)
	{
	case 1:	//i-
		l1_max = abs(ke-ks)+1;
		m1_max = abs(je-js)+1;
		 
		l = 1;
		k = ks;
		do
		{
			m = 1;
			j = js;
			do
			{
				ShadowCell_A_cv[l][m] = A.FluidProp.Cell_cv[nd][j][k];
				ShadowCell_A_pv[l][m] = A.FluidProp.Cell_pv[nd][j][k];
				if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
					ShadowCell_A_ev[l][m] = A.FluidProp.Cell_ev[nd][j][k];
				m = m+1;
				j = j+1;
			}while(j<=je);
			l = l+1;
			k = k+1;
		}while(k<=ke);
		break;

	case 2:	//i+
		l1_max = abs(je-js)+1;
		m1_max = abs(ke-ks)+1;

		l = 1;
		j = js;
		do
		{
			m = 1;
			k = ks;
			do
			{
				ShadowCell_A_cv[l][m] = A.FluidProp.Cell_cv[A.Geometry.IM-nd][j][k];
				ShadowCell_A_pv[l][m] = A.FluidProp.Cell_pv[A.Geometry.IM-nd][j][k];
				if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
					ShadowCell_A_ev[l][m] = A.FluidProp.Cell_ev[A.Geometry.IM-nd][j][k];
				m = m+1;
				k = k+1;			
			}while(k<=ke);
			l = l+1;
			j = j+1;			
		}while(j<=je);
		break;

		case 3:	//j-
		l1_max = abs(ie-is)+1;
		m1_max = abs(ke-ks)+1;

		l = 1;
		i = is;
		do
		{
			m = 1;
			k = ks;
			do
			{
				ShadowCell_A_cv[l][m] = A.FluidProp.Cell_cv[i][nd][k];
				ShadowCell_A_pv[l][m] = A.FluidProp.Cell_pv[i][nd][k];
				if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
					ShadowCell_A_ev[l][m] = A.FluidProp.Cell_ev[i][nd][k];
				m = m+1;
				k = k+1;		
			}while(k<=ke);
			l = l+1;
			i = i+1;		
		}while(i<=ie);
		break;

		case 4:	//j+
		l1_max = abs(ke-ks)+1;
		m1_max = abs(ie-is)+1;

		l = 1;
		k = ks;
		do
		{
			m = 1;
			i = is;
			do
			{
				ShadowCell_A_cv[l][m] = A.FluidProp.Cell_cv[i][A.Geometry.JM-nd][k];
				ShadowCell_A_pv[l][m] = A.FluidProp.Cell_pv[i][A.Geometry.JM-nd][k];
				if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
					ShadowCell_A_ev[l][m] = A.FluidProp.Cell_ev[i][A.Geometry.JM-nd][k];
				m = m+1;
				i = i+1;					
			}while(i<=ie);
			l = l+1;
			k = k+1;			
		}while(k<=ke);
		break;

		case 5:	//k-
		l1_max = abs(je-js)+1;
		m1_max = abs(ie-is)+1;

		l = 1;
		j = js;
		do
		{
			m = 1;
			i = is;
			do
			{
				ShadowCell_A_cv[l][m] = A.FluidProp.Cell_cv[i][j][nd];
				ShadowCell_A_pv[l][m] = A.FluidProp.Cell_pv[i][j][nd];
				if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
					ShadowCell_A_ev[l][m] = A.FluidProp.Cell_ev[i][j][nd];
				m = m+1;
				i = i+1;	
			}while(i<=ie);
			l = l+1;
			j = j+1;
		}while(j<=je);
		break;

		case 6:	//k+
		l1_max = abs(ie-is)+1;
		m1_max = abs(je-js)+1;

		l = 1;
		i = is;
		do
		{
			m = 1;
			j = js;
			do
			{
				ShadowCell_A_cv[l][m] = A.FluidProp.Cell_cv[i][j][A.Geometry.KM-nd];
				ShadowCell_A_pv[l][m] = A.FluidProp.Cell_pv[i][j][A.Geometry.KM-nd];
				if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
					ShadowCell_A_ev[l][m] = A.FluidProp.Cell_ev[i][j][A.Geometry.KM-nd];
				m = m+1;
				j = j+1;
			}while(j<=je);
			l = l+1;
			i = i+1;	
		}while(i<=ie);
		break;
	}
}


void CMesh::Trans(int ori, REAL dr)
{
	int l1, l2, m1, m2;

	switch(ori)
	{
	case 1:
		l2_max = l1_max;
		m2_max = m1_max;

		l1 = 1;
		l2 = 1;
		do
		{
			m1 = 1;
			m2 = m1_max;
			do
			{
				if(dr==0.0)
				{
					ShadowCell_B_cv[l2][m2] = ShadowCell_A_cv[l1][m1];
					ShadowCell_B_pv[l2][m2] = ShadowCell_A_pv[l1][m1];
				}
				else
				{
					ShadowCell_B_pv[l2][m2] = ShadowCell_A_pv[l1][m1];
					ShadowCell_B_cv[l2][m2] = ShadowCell_A_cv[l1][m1];
					ShadowCell_B_cv[l2][m2].xmom = ShadowCell_A_cv[l1][m1].xmom*cos(dr) - ShadowCell_A_cv[l1][m1].ymom*sin(dr);
					ShadowCell_B_cv[l2][m2].ymom = ShadowCell_A_cv[l1][m1].ymom*cos(dr) + ShadowCell_A_cv[l1][m1].xmom*sin(dr);
					ShadowCell_B_pv[l2][m2].uvel = ShadowCell_A_pv[l1][m1].uvel*cos(dr) - ShadowCell_A_pv[l1][m1].vvel*sin(dr);
					ShadowCell_B_pv[l2][m2].vvel = ShadowCell_A_pv[l1][m1].vvel*cos(dr) + ShadowCell_A_pv[l1][m1].uvel*sin(dr);
				}
				if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
					ShadowCell_B_ev[l2][m2] = ShadowCell_A_ev[l1][m1];

				m2 = m2-1;
				m1 = m1+1;
			}while(m1<=m1_max);
			l2 = l2+1;
			l1 = l1+1;			
		}while(l1<=l1_max);
		break;
		 
	case 2:
		l2_max = m1_max;
		m2_max = l1_max;

		l1 = 1;
		m2 = 1;
		do
		{
			m1 = 1;
			l2 = 1;
			do
			{
				if(dr==0.0)
				{					
					ShadowCell_B_cv[l2][m2] = ShadowCell_A_cv[l1][m1];
					ShadowCell_B_pv[l2][m2] = ShadowCell_A_pv[l1][m1];
				}
				else
				{
					ShadowCell_B_pv[l2][m2] = ShadowCell_A_pv[l1][m1];
					ShadowCell_B_cv[l2][m2] = ShadowCell_A_cv[l1][m1];
					ShadowCell_B_cv[l2][m2].xmom = ShadowCell_A_cv[l1][m1].xmom*cos(dr) - ShadowCell_A_cv[l1][m1].ymom*sin(dr);
					ShadowCell_B_cv[l2][m2].ymom = ShadowCell_A_cv[l1][m1].ymom*cos(dr) + ShadowCell_A_cv[l1][m1].xmom*sin(dr);
					ShadowCell_B_pv[l2][m2].uvel = ShadowCell_A_pv[l1][m1].uvel*cos(dr) - ShadowCell_A_pv[l1][m1].vvel*sin(dr);
					ShadowCell_B_pv[l2][m2].vvel = ShadowCell_A_pv[l1][m1].vvel*cos(dr) + ShadowCell_A_pv[l1][m1].uvel*sin(dr);
				}	
				if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
					ShadowCell_B_ev[l2][m2] = ShadowCell_A_ev[l1][m1];

				l2 = l2+1;
				m1 = m1+1;
			}while(m1<=m1_max);
			m2 = m2+1;
			l1 = l1+1;			
		}while(l1<=l1_max);
		break;

	case 3:
		l2_max = l1_max;
		m2_max = m1_max;

		l1 = 1;
		l2 = l1_max;
		do
		{
			m1 = 1;
			m2 = 1;
			do
			{
				if(dr==0.0)
				{
					ShadowCell_B_cv[l2][m2] = ShadowCell_A_cv[l1][m1];
					ShadowCell_B_pv[l2][m2] = ShadowCell_A_pv[l1][m1];
				}
				else
				{
					ShadowCell_B_pv[l2][m2] = ShadowCell_A_pv[l1][m1];
					ShadowCell_B_cv[l2][m2] = ShadowCell_A_cv[l1][m1];
					ShadowCell_B_cv[l2][m2].xmom = ShadowCell_A_cv[l1][m1].xmom*cos(dr) - ShadowCell_A_cv[l1][m1].ymom*sin(dr);
					ShadowCell_B_cv[l2][m2].ymom = ShadowCell_A_cv[l1][m1].ymom*cos(dr) + ShadowCell_A_cv[l1][m1].xmom*sin(dr);
					ShadowCell_B_pv[l2][m2].uvel = ShadowCell_A_pv[l1][m1].uvel*cos(dr) - ShadowCell_A_pv[l1][m1].vvel*sin(dr);
					ShadowCell_B_pv[l2][m2].vvel = ShadowCell_A_pv[l1][m1].vvel*cos(dr) + ShadowCell_A_pv[l1][m1].uvel*sin(dr);
				}
				if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
					ShadowCell_B_ev[l2][m2] = ShadowCell_A_ev[l1][m1];

				m2 = m2+1;
				m1 = m1+1;
			}while(m1<=m1_max);
			l2 = l2-1;
			l1 = l1+1;			
		}while(l1<=l1_max);
		break;

	case 4:
		l2_max = m1_max;
		m2_max = l1_max;

		l1 = 1;
		m2 = l1_max;
		do
		{
			m1 = 1;
			l2 = m1_max;
			do
			{
				if(dr==0.0)
				{
					ShadowCell_B_cv[l2][m2] = ShadowCell_A_cv[l1][m1];
					ShadowCell_B_pv[l2][m2] = ShadowCell_A_pv[l1][m1];
				}
				else
				{
					ShadowCell_B_pv[l2][m2] = ShadowCell_A_pv[l1][m1];
					ShadowCell_B_cv[l2][m2] = ShadowCell_A_cv[l1][m1];
					ShadowCell_B_cv[l2][m2].xmom = ShadowCell_A_cv[l1][m1].xmom*cos(dr) - ShadowCell_A_cv[l1][m1].ymom*sin(dr);
					ShadowCell_B_cv[l2][m2].ymom = ShadowCell_A_cv[l1][m1].ymom*cos(dr) + ShadowCell_A_cv[l1][m1].xmom*sin(dr);
					ShadowCell_B_pv[l2][m2].uvel = ShadowCell_A_pv[l1][m1].uvel*cos(dr) - ShadowCell_A_pv[l1][m1].vvel*sin(dr);
					ShadowCell_B_pv[l2][m2].vvel = ShadowCell_A_pv[l1][m1].vvel*cos(dr) + ShadowCell_A_pv[l1][m1].uvel*sin(dr);
				}
				if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
					ShadowCell_B_ev[l2][m2] = ShadowCell_A_ev[l1][m1];

				l2 = l2-1;
				m1 = m1+1;
			}while(m1<=m1_max);
			m2 = m2-1;
			l1 = l1+1;			
		}while(l1<=l1_max);
		break;
	}
}


void CMesh::ShadowIn(CBlocks &B, int fc, int nd, int is, int ie, int js, int je, int ks, int ke)
{
	int i, j, k, l, m;
	switch(fc)
	{
	case 1:	//i-

		l = 1;
		k = ks;
		do
		{
			m = 1;
			j = js;
			do
			{
				if(nd==1)
				{
					B.FluidProp.Cell_cv[0][j][k] = ShadowCell_B_cv[l][m];
					B.FluidProp.Cell_pv[0][j][k] = ShadowCell_B_pv[l][m];
					if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
						B.FluidProp.Cell_ev[0][j][k] = ShadowCell_B_ev[l][m];
				}								   
				B.FluidProp.DummyIS_cv[nd][j][k] = ShadowCell_B_cv[l][m];
				B.FluidProp.DummyIS_pv[nd][j][k] = ShadowCell_B_pv[l][m];

				m = m+1;
				j = j+1;
			}while(j<=je);
			l = l+1;
			k = k+1;
		}while(k<=ke);
		break;

	case 2:	//i+

		l = 1;
		j = js;
		do
		{
			m = 1;
			k = ks;
			do
			{
				if(nd==1) 
				{
					B.FluidProp.Cell_cv[B.Geometry.IM][j][k] = ShadowCell_B_cv[l][m];
					B.FluidProp.Cell_pv[B.Geometry.IM][j][k] = ShadowCell_B_pv[l][m];
					if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
						B.FluidProp.Cell_ev[B.Geometry.IM][j][k] = ShadowCell_B_ev[l][m];
				}
				B.FluidProp.DummyIP_cv[nd][j][k] = ShadowCell_B_cv[l][m];
				B.FluidProp.DummyIP_pv[nd][j][k] = ShadowCell_B_pv[l][m];

				m = m+1;
				k = k+1;			
			}while(k<=ke);
			l = l+1;
			j = j+1;	
		}while(j<=je);
		break;

		case 3:	//j-

		l = 1;
		i = is;
		do
		{
			m = 1;
			k = ks;
			do
			{
				if(nd==1) 
				{
					B.FluidProp.Cell_cv[i][0][k] = ShadowCell_B_cv[l][m];
					B.FluidProp.Cell_pv[i][0][k] = ShadowCell_B_pv[l][m];
					if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
						B.FluidProp.Cell_ev[i][0][k] = ShadowCell_B_ev[l][m];
				}
				B.FluidProp.DummyJS_cv[i][nd][k] = ShadowCell_B_cv[l][m];
				B.FluidProp.DummyJS_pv[i][nd][k] = ShadowCell_B_pv[l][m];

				m = m+1;
				k = k+1;		
			}while(k<=ke);
			l = l+1;
			i = i+1;		
		}while(i<=ie);
		break;

		case 4:	//j+

		l = 1;
		k = ks;
		do
		{
			m = 1;
			i = is;
			do
			{
				if(nd==1)
				{
					B.FluidProp.Cell_cv[i][B.Geometry.JM][k] = ShadowCell_B_cv[l][m];
					B.FluidProp.Cell_pv[i][B.Geometry.JM][k] = ShadowCell_B_pv[l][m];
					if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
						B.FluidProp.Cell_ev[i][B.Geometry.JM][k] = ShadowCell_B_ev[l][m];
				}
				B.FluidProp.DummyJP_cv[i][nd][k] = ShadowCell_B_cv[l][m];
				B.FluidProp.DummyJP_pv[i][nd][k] = ShadowCell_B_pv[l][m];

				m = m+1;
				i = i+1;	
			}while(i<=ie);
			l = l+1;
			k = k+1;
		}while(k<=ke);
		break;

		case 5:	//k-

		l = 1;
		j = js;
		do
		{
			m = 1;
			i = is;
			do
			{
				if(nd==1)
				{
					B.FluidProp.Cell_cv[i][j][0] = ShadowCell_B_cv[l][m];
					B.FluidProp.Cell_pv[i][j][0] = ShadowCell_B_pv[l][m];
					if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
						B.FluidProp.Cell_ev[i][j][0] = ShadowCell_B_ev[l][m];
				}
				B.FluidProp.DummyKS_cv[i][j][nd] = ShadowCell_B_cv[l][m];
				B.FluidProp.DummyKS_pv[i][j][nd] = ShadowCell_B_pv[l][m];

				m = m+1;
				i = i+1;
			}while(i<=ie);
			l = l+1;
			j = j+1;
				
		}while(j<=je);
		break;

		case 6:	//k+

		l = 1;
		i = is;
		do
		{
			m = 1;
			j = js;
			do
			{
				if(nd==1)
				{
					B.FluidProp.Cell_cv[i][j][B.Geometry.KM] = ShadowCell_B_cv[l][m];
					B.FluidProp.Cell_pv[i][j][B.Geometry.KM] = ShadowCell_B_pv[l][m];
					if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
						B.FluidProp.Cell_ev[i][j][B.Geometry.KM] = ShadowCell_B_ev[l][m];
				}
				B.FluidProp.DummyKP_cv[i][j][nd] = ShadowCell_B_cv[l][m];
				B.FluidProp.DummyKP_pv[i][j][nd] = ShadowCell_B_pv[l][m];

				m = m+1;
				j = j+1;
			}while(j<=je);
			l = l+1;
			i = i+1;			
		}while(i<=ie);
		break;
	}
}