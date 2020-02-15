


#include "defs.h"
#include "Geometry.h"
#include <iostream>

using namespace std;

void CGeometry::FaceVectorVolume()
{
	int i, j, k;
	TNode medA, medB;

#pragma omp parallel num_threads(NT) default(shared) private(i, j, k, medA, medB)
	{
#pragma omp for nowait 
		// FaceI Geometry : face vector, area, coordinate 
		for(i=1;i<=IM;i++)
		{
			for(j=1;j<=JM-1;j++)
			{
				for(k=1;k<=KM-1;k++)
				{
					medA = Node_coord[i][j+1][k] - Node_coord[i][j][k+1];
					medB = Node_coord[i][j+1][k+1] - Node_coord[i][j][k];
					FaceI[i][j][k].S.x = 0.5f*(medA.y*medB.z - medB.y*medA.z);
					FaceI[i][j][k].S.y = 0.5f*(medA.z*medB.x - medA.x*medB.z);
					FaceI[i][j][k].S.z = 0.5f*(medA.x*medB.y - medB.x*medA.y);
					FaceI[i][j][k].coord  = 0.25*(Node_coord[i][j][k] + Node_coord[i][j+1][k] + Node_coord[i][j+1][k+1] + Node_coord[i][j][k+1]);
					FaceI[i][j][k].r  = SQRT(FaceI[i][j][k].coord.x*FaceI[i][j][k].coord.x + FaceI[i][j][k].coord.y*FaceI[i][j][k].coord.y);
					FaceI[i][j][k].th = ATAN(FaceI[i][j][k].coord.y / FaceI[i][j][k].coord.x);
					FaceI[i][j][k].s  = SQRT(FaceI[i][j][k].S.x*FaceI[i][j][k].S.x + FaceI[i][j][k].S.y*FaceI[i][j][k].S.y + FaceI[i][j][k].S.z*FaceI[i][j][k].S.z);
				}
			}
		}
#pragma omp for nowait 
		// FaceJ Geometry : face vector, area, coordinate
		for(i=1;i<=IM-1;i++)
		{
			for(j=1;j<=JM;j++)
			{
				for(k=1;k<=KM-1;k++)
				{
					medA = Node_coord[i][j][k+1] - Node_coord[i+1][j][k];
					medB = Node_coord[i+1][j][k+1] - Node_coord[i][j][k];
					FaceJ[i][j][k].S.x = 0.5f*(medA.y*medB.z - medB.y*medA.z);
					FaceJ[i][j][k].S.y = 0.5f*(medA.z*medB.x - medA.x*medB.z);
					FaceJ[i][j][k].S.z = 0.5f*(medA.x*medB.y - medB.x*medA.y);
					FaceJ[i][j][k].coord  = 0.25*(Node_coord[i][j][k] + Node_coord[i+1][j][k] + Node_coord[i+1][j][k+1] + Node_coord[i][j][k+1]);
					FaceJ[i][j][k].r  = SQRT(FaceJ[i][j][k].coord.x*FaceJ[i][j][k].coord.x + FaceJ[i][j][k].coord.y*FaceJ[i][j][k].coord.y);
					FaceJ[i][j][k].th = ATAN(FaceJ[i][j][k].coord.y / FaceJ[i][j][k].coord.x);
					FaceJ[i][j][k].s  = SQRT(FaceJ[i][j][k].S.x*FaceJ[i][j][k].S.x + FaceJ[i][j][k].S.y*FaceJ[i][j][k].S.y + FaceJ[i][j][k].S.z*FaceJ[i][j][k].S.z);
				}
			}
		}
#pragma omp for
		// FaceK Geometry : face vector, area, coordinate
		for(i=1;i<=IM-1;i++)
		{
			for(j=1;j<=JM-1;j++)
			{
				for(k=1;k<=KM;k++)
				{
					medA = Node_coord[i+1][j+1][k] - Node_coord[i][j][k];
					medB = Node_coord[i][j+1][k] - Node_coord[i+1][j][k];
					FaceK[i][j][k].S.x = 0.5f*(medA.y*medB.z - medB.y*medA.z);
					FaceK[i][j][k].S.y = 0.5f*(medA.z*medB.x - medA.x*medB.z);
					FaceK[i][j][k].S.z = 0.5f*(medA.x*medB.y - medB.x*medA.y);
					FaceK[i][j][k].coord  = 0.25*(Node_coord[i][j][k] + Node_coord[i+1][j][k] + Node_coord[i+1][j+1][k] + Node_coord[i][j+1][k]);
					FaceK[i][j][k].r  = SQRT(FaceK[i][j][k].coord.x*FaceK[i][j][k].coord.x + FaceK[i][j][k].coord.y*FaceK[i][j][k].coord.y);
					FaceK[i][j][k].th = ATAN(FaceK[i][j][k].coord.y / FaceK[i][j][k].coord.x);
					FaceK[i][j][k].s  = SQRT(FaceK[i][j][k].S.x*FaceK[i][j][k].S.x + FaceK[i][j][k].S.y*FaceK[i][j][k].S.y + FaceK[i][j][k].S.z*FaceK[i][j][k].S.z);
				}
			}
		}
#pragma omp for 
		// Cell Geometry : cell volume, coordinate
		for(i=1;i<=IM-1;i++)
		{
			for(j=1;j<=JM-1;j++)
			{
				for(k=1;k<=KM-1;k++)
				{
					Cell[i][j][k].Vol = FaceI[i+1][j][k].coord*FaceI[i+1][j][k].S - FaceI[i][j][k].coord*FaceI[i][j][k].S + 
						FaceJ[i][j+1][k].coord*FaceJ[i][j+1][k].S - FaceJ[i][j][k].coord*FaceJ[i][j][k].S +
						FaceK[i][j][k+1].coord*FaceK[i][j][k+1].S - FaceK[i][j][k].coord*FaceK[i][j][k].S;
					Cell[i][j][k].Vol = Cell[i][j][k].Vol / 3.0f;


					Cell[i][j][k].coord = 0.125*(Node_coord[i][j][k]   + Node_coord[i+1][j][k]   + Node_coord[i][j+1][k]   + Node_coord[i+1][j+1][k]
					+ Node_coord[i][j][k+1] + Node_coord[i+1][j][k+1] + Node_coord[i][j+1][k+1] + Node_coord[i+1][j+1][k+1]);

					Cell[i][j][k].r		  = SQRT(Cell[i][j][k].coord.x*Cell[i][j][k].coord.x + Cell[i][j][k].coord.y*Cell[i][j][k].coord.y);
					Cell[i][j][k].th	  = ATAN(Cell[i][j][k].coord.y / Cell[i][j][k].coord.x);
					CellVol				  = CellVol + Cell[i][j][k].Vol;
					if(Cell[i][j][k].Vol<=0.0)
					{
						cout<<" Volume is minus: "<<endl;
					}
				}
			}
		}
	}
}

//  jocabian conversion coefficient
void CGeometry::JacobiConvertCoeff()
{
	int  i, j, k;

	REAL partA, partB, partC,
		x_kesi, y_kesi, z_kesi,
		x_yita, y_yita, z_yita,
		x_zita, y_zita, z_zita,
		Jocabian1, Jocabian, lim_zero = 1E-20f,		// less than lim_zero which can be as zero
		x_kesi_I, y_kesi_I, z_kesi_I,
		x_kesi_J, y_kesi_J, z_kesi_J,
		x_kesi_K, y_kesi_K, z_kesi_K,
		x_yita_I, y_yita_I, z_yita_I,
		x_yita_J, y_yita_J, z_yita_J,
		x_yita_K, y_yita_K, z_yita_K,
		x_zita_I, y_zita_I, z_zita_I,
		x_zita_J, y_zita_J, z_zita_J,
		x_zita_K, y_zita_K, z_zita_K,
		Jocabian_I, Jocabian_J, Jocabian_K;
#pragma omp parallel num_threads(NT) default(shared) private(i, j, k, partA, partB, partC, x_kesi, y_kesi, z_kesi, x_yita, y_yita, z_yita, x_zita, y_zita, z_zita, Jocabian1, Jocabian, \
	x_kesi_I, y_kesi_I, z_kesi_I, x_kesi_J, y_kesi_J, z_kesi_J, x_kesi_K, y_kesi_K, z_kesi_K, x_yita_I, y_yita_I, z_yita_I, x_yita_J, y_yita_J, z_yita_J, x_yita_K, y_yita_K, z_yita_K, \
	x_zita_I, y_zita_I, z_zita_I, x_zita_J, y_zita_J, z_zita_J, x_zita_K, y_zita_K, z_zita_K, Jocabian_I, Jocabian_J, Jocabian_K)
	{
		#pragma omp for nowait 
		// Cell I J K partial derivative
		for(i=1;i<=IM-1;i++)  
		{
			for(j=1;j<=JM-1;j++)
			{
				for(k=1;k<=KM-1;k++)
				{
					if (i==1)
					{
						x_kesi = Cell[i+1][j][k].coord.x - Cell[i][j][k].coord.x;
						y_kesi = Cell[i+1][j][k].coord.y - Cell[i][j][k].coord.y;
						z_kesi = Cell[i+1][j][k].coord.z - Cell[i][j][k].coord.z;
					}
					else if (i==IM-1)
					{
						x_kesi = Cell[i][j][k].coord.x - Cell[i-1][j][k].coord.x;
						y_kesi = Cell[i][j][k].coord.y - Cell[i-1][j][k].coord.y;
						z_kesi = Cell[i][j][k].coord.z - Cell[i-1][j][k].coord.z;
					}
					else
					{
						x_kesi = 0.5f*(Cell[i+1][j][k].coord.x - Cell[i-1][j][k].coord.x);
						y_kesi = 0.5f*(Cell[i+1][j][k].coord.y - Cell[i-1][j][k].coord.y);
						z_kesi = 0.5f*(Cell[i+1][j][k].coord.z - Cell[i-1][j][k].coord.z);
					}
					if (j==1)
					{
						x_yita = (Cell[i][j+1][k].coord.x - Cell[i][j][k].coord.x);
						y_yita = (Cell[i][j+1][k].coord.y - Cell[i][j][k].coord.y);
						z_yita = (Cell[i][j+1][k].coord.z - Cell[i][j][k].coord.z);
					}
					else if (j==JM-1)
					{
						x_yita = (Cell[i][j][k].coord.x - Cell[i][j-1][k].coord.x);
						y_yita = (Cell[i][j][k].coord.y - Cell[i][j-1][k].coord.y);
						z_yita = (Cell[i][j][k].coord.z - Cell[i][j-1][k].coord.z);
					}		   	 
					else	   	 
					{		   	 
						x_yita = 0.5f*(Cell[i][j+1][k].coord.x - Cell[i][j-1][k].coord.x);
						y_yita = 0.5f*(Cell[i][j+1][k].coord.y - Cell[i][j-1][k].coord.y);
						z_yita = 0.5f*(Cell[i][j+1][k].coord.z - Cell[i][j-1][k].coord.z);
					}		   	 
					if (k==1)   	 
					{		   	 
						x_zita = (Cell[i][j][k+1].coord.x - Cell[i][j][k].coord.x);
						y_zita = (Cell[i][j][k+1].coord.y - Cell[i][j][k].coord.y);
						z_zita = (Cell[i][j][k+1].coord.z - Cell[i][j][k].coord.z);
					}
					else if (k==KM-1)
					{
						x_zita = (Cell[i][j][k].coord.x - Cell[i][j][k-1].coord.x);
						y_zita = (Cell[i][j][k].coord.y - Cell[i][j][k-1].coord.y);
						z_zita = (Cell[i][j][k].coord.z - Cell[i][j][k-1].coord.z);
					}
					else
					{
						x_zita = 0.5f*(Cell[i][j][k+1].coord.x - Cell[i][j][k-1].coord.x);
						y_zita = 0.5f*(Cell[i][j][k+1].coord.y - Cell[i][j][k-1].coord.y);
						z_zita = 0.5f*(Cell[i][j][k+1].coord.z - Cell[i][j][k-1].coord.z);
					}

					partA    = (x_kesi * (y_yita*z_zita - y_zita*z_yita));
					partB    = (x_yita * (y_kesi*z_zita - y_zita*z_kesi));
					partC    = (x_zita * (y_kesi*z_yita - y_yita*z_kesi));
					Jocabian1 = (partA - partB + partC);
					if(ABS(Jocabian1) <= lim_zero)
						Jocabian1 = lim_zero;
					Jocabian = 1.0f / Jocabian1;
					Cell[i][j][k].gradx_kesi = Jocabian * (y_yita*z_zita - y_zita*z_yita);
					Cell[i][j][k].grady_kesi = Jocabian * (z_yita*x_zita - z_zita*x_yita);
					Cell[i][j][k].gradz_kesi = Jocabian * (y_yita*z_zita - y_zita*z_yita);
					Cell[i][j][k].gradx_yita = Jocabian * (y_zita*z_kesi - y_kesi*z_zita);
					Cell[i][j][k].grady_yita = Jocabian * (z_zita*x_kesi - z_kesi*x_zita);
					Cell[i][j][k].gradz_yita = Jocabian * (x_zita*y_kesi - x_kesi*y_zita);
					Cell[i][j][k].gradx_zita = Jocabian * (y_kesi*z_yita - y_yita*z_kesi);
					Cell[i][j][k].grady_zita = Jocabian * (z_kesi*x_yita - z_yita*x_kesi);
					Cell[i][j][k].gradz_zita = Jocabian * (x_kesi*y_yita - x_yita*y_kesi);
				}
			}
		}

		#pragma omp for nowait 
		//  Face I J K partial derivative
		for(i=1;i<=IM;i++) //I direction geometry var for viscous
		{
			for(j=1;j<=JM-1;j++)
			{
				for(k=1;k<=KM-1;k++)
				{
					//kesai direction
					if(i==1)
					{
						x_kesi_I = FaceI[i+1][j][k].coord.x - FaceI[i][j][k].coord.x;
						y_kesi_I = FaceI[i+1][j][k].coord.y - FaceI[i][j][k].coord.y;
						z_kesi_I = FaceI[i+1][j][k].coord.z - FaceI[i][j][k].coord.z;
					}
					else if(i==IM)
					{
						x_kesi_I = FaceI[i][j][k].coord.x - FaceI[i-1][j][k].coord.x;
						y_kesi_I = FaceI[i][j][k].coord.y - FaceI[i-1][j][k].coord.y;
						z_kesi_I = FaceI[i][j][k].coord.z - FaceI[i-1][j][k].coord.z;
					}
					else
					{
						x_kesi_I = 0.5f*(FaceI[i+1][j][k].coord.x - FaceI[i-1][j][k].coord.x);
						y_kesi_I = 0.5f*(FaceI[i+1][j][k].coord.y - FaceI[i-1][j][k].coord.y);
						z_kesi_I = 0.5f*(FaceI[i+1][j][k].coord.z - FaceI[i-1][j][k].coord.z);
					}

					// yita direction
					x_yita_I = 0.5f*(Node_coord[i][j+1][k].x + Node_coord[i][j+1][k+1].x) - 0.5f*(Node_coord[i][j][k].x + Node_coord[i][j][k+1].x);
					y_yita_I = 0.5f*(Node_coord[i][j+1][k].y + Node_coord[i][j+1][k+1].y) - 0.5f*(Node_coord[i][j][k].y + Node_coord[i][j][k+1].y);
					z_yita_I = 0.5f*(Node_coord[i][j+1][k].z + Node_coord[i][j+1][k+1].z) - 0.5f*(Node_coord[i][j][k].z + Node_coord[i][j][k+1].z);

					// zita  direction				    					   					      
					x_zita_I = 0.5f*(Node_coord[i][j][k+1].x + Node_coord[i][j+1][k+1].x) - 0.5f*(Node_coord[i][j][k].x + Node_coord[i][j][k].x);
					y_zita_I = 0.5f*(Node_coord[i][j][k+1].y + Node_coord[i][j+1][k+1].y) - 0.5f*(Node_coord[i][j][k].y + Node_coord[i][j][k].y);
					z_zita_I = 0.5f*(Node_coord[i][j][k+1].z + Node_coord[i][j+1][k+1].z) - 0.5f*(Node_coord[i][j][k].z + Node_coord[i][j][k].z);

					//FaceI coordinate partial derivative
					partA = (x_kesi_I * (y_yita_I*z_zita_I - y_zita_I*z_yita_I));
					partB = (x_yita_I * (y_kesi_I*z_zita_I - y_zita_I*z_kesi_I));
					partC = (x_zita_I * (y_kesi_I*z_yita_I - y_yita_I*z_kesi_I));
					Jocabian1 = (partA - partB + partC);
					if(ABS(Jocabian1) <= lim_zero)
						Jocabian1 = lim_zero;
					Jocabian_I = 1.0f / Jocabian1;
					FaceI[i][j][k].gradx_kesi = Jocabian_I * (y_yita_I*z_zita_I - y_zita_I*z_yita_I);
					FaceI[i][j][k].grady_kesi = Jocabian_I * (z_yita_I*x_zita_I - z_zita_I*x_yita_I);
					FaceI[i][j][k].gradz_kesi = Jocabian_I * (y_yita_I*z_zita_I - y_zita_I*z_yita_I);
					FaceI[i][j][k].gradx_yita = Jocabian_I * (y_zita_I*z_kesi_I - y_kesi_I*z_zita_I);
					FaceI[i][j][k].grady_yita = Jocabian_I * (z_zita_I*x_kesi_I - z_kesi_I*x_zita_I);
					FaceI[i][j][k].gradz_yita = Jocabian_I * (x_zita_I*y_kesi_I - x_kesi_I*y_zita_I);
					FaceI[i][j][k].gradx_zita = Jocabian_I * (y_kesi_I*z_yita_I - y_yita_I*z_kesi_I);
					FaceI[i][j][k].grady_zita = Jocabian_I * (z_kesi_I*x_yita_I - z_yita_I*x_kesi_I);
					FaceI[i][j][k].gradz_zita = Jocabian_I * (x_kesi_I*y_yita_I - x_yita_I*y_kesi_I);
				}
			}
		}

		#pragma omp for nowait 		
		for(i=1;i<=IM-1;i++) //J direction geometry var for viscous
		{
			for(j=1;j<=JM;j++)
			{
				for(k=1;k<=KM-1;k++)
				{
					// kesai direction
					x_kesi_J = 0.5f*(Node_coord[i+1][j][k].x + Node_coord[i+1][j][k+1].x) - 0.5f*(Node_coord[i][j][k].x + Node_coord[i][j][k+1].x);
					y_kesi_J = 0.5f*(Node_coord[i+1][j][k].y + Node_coord[i+1][j][k+1].y) - 0.5f*(Node_coord[i][j][k].y + Node_coord[i][j][k+1].y);
					z_kesi_J = 0.5f*(Node_coord[i+1][j][k].z + Node_coord[i+1][j][k+1].z) - 0.5f*(Node_coord[i][j][k].z + Node_coord[i][j][k+1].z);

					// yita direction
					if(j==1)
					{
						x_yita_J = (FaceJ[i][j+1][k].coord.x - FaceJ[i][j][k].coord.x);
						y_yita_J = (FaceJ[i][j+1][k].coord.y - FaceJ[i][j][k].coord.y);
						z_yita_J = (FaceJ[i][j+1][k].coord.z - FaceJ[i][j][k].coord.z);
					}
					else if(j==JM)
					{
						x_yita_J = (FaceJ[i][j][k].coord.x - FaceJ[i][j-1][k].coord.x);
						y_yita_J = (FaceJ[i][j][k].coord.y - FaceJ[i][j-1][k].coord.y);
						z_yita_J = (FaceJ[i][j][k].coord.z - FaceJ[i][j-1][k].coord.z);
					}
					else
					{
						x_yita_J = 0.5f*(FaceJ[i][j+1][k].coord.x - FaceJ[i][j-1][k].coord.x);
						y_yita_J = 0.5f*(FaceJ[i][j+1][k].coord.y - FaceJ[i][j-1][k].coord.y);
						z_yita_J = 0.5f*(FaceJ[i][j+1][k].coord.z - FaceJ[i][j-1][k].coord.z);
					}
					// zita direction
					x_zita_J = 0.5f*(Node_coord[i][j][k+1].x + Node_coord[i+1][j][k+1].x) - 0.5f*(Node_coord[i][j][k].x + Node_coord[i+1][j][k].x);
					y_zita_J = 0.5f*(Node_coord[i][j][k+1].y + Node_coord[i+1][j][k+1].y) - 0.5f*(Node_coord[i][j][k].y + Node_coord[i+1][j][k].y);
					z_zita_J = 0.5f*(Node_coord[i][j][k+1].z + Node_coord[i+1][j][k+1].z) - 0.5f*(Node_coord[i][j][k].z + Node_coord[i+1][j][k].z);

					//FaceJ coordinate partial derivative
					partA = (x_kesi_J * (y_yita_J*z_zita_J - y_zita_J*z_yita_J));
					partB = (x_yita_J * (y_kesi_J*z_zita_J - y_zita_J*z_kesi_J));
					partC = (x_zita_J * (y_kesi_J*z_yita_J - y_yita_J*z_kesi_J));

					Jocabian1 = (partA - partB + partC);
					if(ABS(Jocabian1) <= lim_zero)
						Jocabian1 = lim_zero;
					Jocabian_J = 1.0f / Jocabian1;
					FaceJ[i][j][k].gradx_kesi = Jocabian_J * (y_yita_J*z_zita_J - y_zita_J*z_yita_J);
					FaceJ[i][j][k].grady_kesi = Jocabian_J * (z_yita_J*x_zita_J - z_zita_J*x_yita_J);
					FaceJ[i][j][k].gradz_kesi = Jocabian_J * (y_yita_J*z_zita_J - y_zita_J*z_yita_J);
					FaceJ[i][j][k].gradx_yita = Jocabian_J * (y_zita_J*z_kesi_J - y_kesi_J*z_zita_J);
					FaceJ[i][j][k].grady_yita = Jocabian_J * (z_zita_J*x_kesi_J - z_kesi_J*x_zita_J);
					FaceJ[i][j][k].gradz_yita = Jocabian_J * (x_zita_J*y_kesi_J - x_kesi_J*y_zita_J);
					FaceJ[i][j][k].gradx_zita = Jocabian_J * (y_kesi_J*z_yita_J - y_yita_J*z_kesi_J);
					FaceJ[i][j][k].grady_zita = Jocabian_J * (z_kesi_J*x_yita_J - z_yita_J*x_kesi_J);
					FaceJ[i][j][k].gradz_zita = Jocabian_J * (x_kesi_J*y_yita_J - x_yita_J*y_kesi_J);

				}
			}
		}

		#pragma omp for nowait 	
		for(i=1;i<=IM-1;i++) //K direction geometry var for viscous
		{
			for(j=1;j<=JM-1;j++)
			{
				for(k=1;k<=KM;k++)
				{
					//kesai direction
					x_kesi_K = 0.5f*(Node_coord[i+1][j+1][k].x + Node_coord[i+1][j][k].x) - 0.5f*(Node_coord[i][j][k].x + Node_coord[i][j+1][k].x);
					y_kesi_K = 0.5f*(Node_coord[i+1][j+1][k].y + Node_coord[i+1][j][k].y) - 0.5f*(Node_coord[i][j][k].y + Node_coord[i][j+1][k].y);
					z_kesi_K = 0.5f*(Node_coord[i+1][j+1][k].z + Node_coord[i+1][j][k].z) - 0.5f*(Node_coord[i][j][k].z + Node_coord[i][j+1][k].z);

					//yita direction
					x_yita_K = 0.5f*(Node_coord[i+1][j+1][k].x + Node_coord[i][j+1][k].x) - 0.5f*(Node_coord[i+1][j][k].x + Node_coord[i][j][k].x);
					y_yita_K = 0.5f*(Node_coord[i+1][j+1][k].y + Node_coord[i][j+1][k].y) - 0.5f*(Node_coord[i+1][j][k].y + Node_coord[i][j][k].y);
					z_yita_K = 0.5f*(Node_coord[i+1][j+1][k].z + Node_coord[i][j+1][k].z) - 0.5f*(Node_coord[i+1][j][k].z + Node_coord[i][j][k].z);

					//zita direction
					if(k==1)
					{
						x_zita_K = (FaceK[i][j][k+1].coord.x - FaceK[i][j][k].coord.x);
						y_zita_K = (FaceK[i][j][k+1].coord.y - FaceK[i][j][k].coord.y);
						z_zita_K = (FaceK[i][j][k+1].coord.z - FaceK[i][j][k].coord.z);
					}
					else if(k==KM)
					{
						x_zita_K = (FaceK[i][j][k].coord.x - FaceK[i][j][k-1].coord.x);
						y_zita_K = (FaceK[i][j][k].coord.y - FaceK[i][j][k-1].coord.y);
						z_zita_K = (FaceK[i][j][k].coord.z - FaceK[i][j][k-1].coord.z);
					}
					else
					{
						x_zita_K = 0.5f*(FaceK[i][j][k+1].coord.x - FaceK[i][j][k-1].coord.x);
						y_zita_K = 0.5f*(FaceK[i][j][k+1].coord.y - FaceK[i][j][k-1].coord.y);
						z_zita_K = 0.5f*(FaceK[i][j][k+1].coord.z - FaceK[i][j][k-1].coord.z);
					}

					//FaceK coordinate partial derivative
					partA = (x_kesi_K * (y_yita_K*z_zita_K - y_zita_K*z_yita_K));
					partB = (x_yita_K * (y_kesi_K*z_zita_K - y_zita_K*z_kesi_K));
					partC = (x_zita_K * (y_kesi_K*z_yita_K - y_yita_K*z_kesi_K));

					Jocabian1 = (partA - partB + partC);
					if(ABS(Jocabian1) <= lim_zero)
						Jocabian1 = lim_zero;
					Jocabian_K = 1.0f / Jocabian1;
					FaceK[i][j][k].gradx_kesi = Jocabian_K * (y_yita_K*z_zita_K - y_zita_K*z_yita_K);
					FaceK[i][j][k].grady_kesi = Jocabian_K * (z_yita_K*x_zita_K - z_zita_K*x_yita_K);
					FaceK[i][j][k].gradz_kesi = Jocabian_K * (y_yita_K*z_zita_K - y_zita_K*z_yita_K);
					FaceK[i][j][k].gradx_yita = Jocabian_K * (y_zita_K*z_kesi_K - y_kesi_K*z_zita_K);
					FaceK[i][j][k].grady_yita = Jocabian_K * (z_zita_K*x_kesi_K - z_kesi_K*x_zita_K);
					FaceK[i][j][k].gradz_yita = Jocabian_K * (x_zita_K*y_kesi_K - x_kesi_K*y_zita_K);
					FaceK[i][j][k].gradx_zita = Jocabian_K * (y_kesi_K*z_yita_K - y_yita_K*z_kesi_K);
					FaceK[i][j][k].grady_zita = Jocabian_K * (z_kesi_K*x_yita_K - z_yita_K*x_kesi_K);
					FaceK[i][j][k].gradz_zita = Jocabian_K * (x_kesi_K*y_yita_K - x_yita_K*y_kesi_K);			
				}
			}
		}
	}
}

