
#include "Mesh.h"
#include <iostream>
#include <fstream>
#include "defs.h"

using namespace std;

void CMesh::WallDistanceSA()
{
	int i, j, k, nb, ib, jb, kb, iw;
	REAL dis_min, dis_temp, x_temp, y_temp, z_temp;
	
	
	for(nb=1;nb<=NB;nb++)
	{
		cout<<"	Block:"<<nb<<" wall distance Caculating..."<<endl;
#pragma omp parallel for num_threads(NT) default(shared) private(i,j,k, ib, jb, kb, iw, dis_min, dis_temp, x_temp, y_temp, z_temp)
		for(i=1;i<=Block[nb].Geometry.IM-1;i++)
		{
			for(j=1;j<=Block[nb].Geometry.JM-1;j++)
			{
				for(k=1;k<=Block[nb].Geometry.KM-1;k++)
				{
					dis_min = 10000.0;
					/*------------wall boundary cycle-----------------------------------*/
					for(iw=1;iw<=BndGeom.nWall;iw++)
					{
						if(BndGeom.Wall[iw].Face_type==1)
						{
							ib = 1;
							for(jb=BndGeom.Wall[iw].js;jb<=BndGeom.Wall[iw].je;jb++)
							{
								for(kb=BndGeom.Wall[iw].ks;kb<=BndGeom.Wall[iw].ke;kb++)
								{
									x_temp   = Block[nb].Geometry.Cell[i][j][k].coord.x - Block[BndGeom.Wall[iw].Block_ID].Geometry.FaceI[ib][jb][kb].coord.x;
									y_temp   = Block[nb].Geometry.Cell[i][j][k].coord.y - Block[BndGeom.Wall[iw].Block_ID].Geometry.FaceI[ib][jb][kb].coord.y;
									z_temp   = Block[nb].Geometry.Cell[i][j][k].coord.z - Block[BndGeom.Wall[iw].Block_ID].Geometry.FaceI[ib][jb][kb].coord.z;
									dis_temp = POW(x_temp*x_temp + y_temp*y_temp + z_temp*z_temp,0.5);
									if(dis_min>=dis_temp)
									{
										dis_min = dis_temp;
									}
								}
							}
						}
						else if(BndGeom.Wall[iw].Face_type==2)
						{
							ib = Block[BndGeom.Wall[iw].Block_ID].Geometry.IM;
							for(jb=BndGeom.Wall[iw].js;jb<=BndGeom.Wall[iw].je;jb++)
							{
								for(kb=BndGeom.Wall[iw].ks;kb<=BndGeom.Wall[iw].ke;kb++)
								{
									x_temp   = Block[nb].Geometry.Cell[i][j][k].coord.x - Block[BndGeom.Wall[iw].Block_ID].Geometry.FaceI[ib][jb][kb].coord.x;
									y_temp   = Block[nb].Geometry.Cell[i][j][k].coord.y - Block[BndGeom.Wall[iw].Block_ID].Geometry.FaceI[ib][jb][kb].coord.y;
									z_temp   = Block[nb].Geometry.Cell[i][j][k].coord.z - Block[BndGeom.Wall[iw].Block_ID].Geometry.FaceI[ib][jb][kb].coord.z;
									dis_temp = POW(x_temp*x_temp + y_temp*y_temp + z_temp*z_temp,0.5);
									if(dis_min>=dis_temp)
									{
										dis_min = dis_temp;
									}
								}
							}
						}
						else if(BndGeom.Wall[iw].Face_type==3)
						{
							jb = 1;
							for(ib=BndGeom.Wall[iw].is;ib<=BndGeom.Wall[iw].ie;ib++)
							{
								for(kb=BndGeom.Wall[iw].ks;kb<=BndGeom.Wall[iw].ke;kb++)
								{
									x_temp   = Block[nb].Geometry.Cell[i][j][k].coord.x - Block[BndGeom.Wall[iw].Block_ID].Geometry.FaceJ[ib][jb][kb].coord.x;
									y_temp   = Block[nb].Geometry.Cell[i][j][k].coord.y - Block[BndGeom.Wall[iw].Block_ID].Geometry.FaceJ[ib][jb][kb].coord.y;
									z_temp   = Block[nb].Geometry.Cell[i][j][k].coord.z - Block[BndGeom.Wall[iw].Block_ID].Geometry.FaceJ[ib][jb][kb].coord.z;
									dis_temp = POW(x_temp*x_temp + y_temp*y_temp + z_temp*z_temp,0.5);
									if(dis_min>=dis_temp)
									{
										dis_min = dis_temp;
									}
								}
							}
						}
						else if(BndGeom.Wall[iw].Face_type==4)
						{
							jb = Block[BndGeom.Wall[iw].Block_ID].Geometry.JM;
							for(ib=BndGeom.Wall[iw].is;ib<=BndGeom.Wall[iw].ie;ib++)
							{
								for(kb=BndGeom.Wall[iw].ks;kb<=BndGeom.Wall[iw].ke;kb++)
								{
									x_temp   = Block[nb].Geometry.Cell[i][j][k].coord.x - Block[BndGeom.Wall[iw].Block_ID].Geometry.FaceJ[ib][jb][kb].coord.x;
									y_temp   = Block[nb].Geometry.Cell[i][j][k].coord.y - Block[BndGeom.Wall[iw].Block_ID].Geometry.FaceJ[ib][jb][kb].coord.y;
									z_temp   = Block[nb].Geometry.Cell[i][j][k].coord.z - Block[BndGeom.Wall[iw].Block_ID].Geometry.FaceJ[ib][jb][kb].coord.z;
									dis_temp = POW(x_temp*x_temp + y_temp*y_temp + z_temp*z_temp,0.5);
									if(dis_min>=dis_temp)
									{
										dis_min = dis_temp;
									}
								}
							}
						}
						else if(BndGeom.Wall[iw].Face_type==5)
						{
							kb = 1;
							for(ib=BndGeom.Wall[iw].is;ib<=BndGeom.Wall[iw].ie;ib++)
							{
								for(jb=BndGeom.Wall[iw].js;jb<=BndGeom.Wall[iw].je;jb++)
								{
									x_temp   = Block[nb].Geometry.Cell[i][j][k].coord.x - Block[BndGeom.Wall[iw].Block_ID].Geometry.FaceK[ib][jb][kb].coord.x;
									y_temp   = Block[nb].Geometry.Cell[i][j][k].coord.y - Block[BndGeom.Wall[iw].Block_ID].Geometry.FaceK[ib][jb][kb].coord.y;
									z_temp   = Block[nb].Geometry.Cell[i][j][k].coord.z - Block[BndGeom.Wall[iw].Block_ID].Geometry.FaceK[ib][jb][kb].coord.z;
									dis_temp = POW(x_temp*x_temp + y_temp*y_temp + z_temp*z_temp,0.5);
									if(dis_min>=dis_temp)
									{
										dis_min = dis_temp;
									}
								}
							}
						}
						else if(BndGeom.Wall[iw].Face_type==6)
						{
							kb = Block[BndGeom.Wall[iw].Block_ID].Geometry.KM;
							for(ib=BndGeom.Wall[iw].is;ib<=BndGeom.Wall[iw].ie;ib++)
							{
								for(jb=BndGeom.Wall[iw].js;jb<=BndGeom.Wall[iw].je;jb++)
								{
									x_temp   = Block[nb].Geometry.Cell[i][j][k].coord.x - Block[BndGeom.Wall[iw].Block_ID].Geometry.FaceK[ib][jb][kb].coord.x;
									y_temp   = Block[nb].Geometry.Cell[i][j][k].coord.y - Block[BndGeom.Wall[iw].Block_ID].Geometry.FaceK[ib][jb][kb].coord.y;
									z_temp   = Block[nb].Geometry.Cell[i][j][k].coord.z - Block[BndGeom.Wall[iw].Block_ID].Geometry.FaceK[ib][jb][kb].coord.z;
									dis_temp = POW(x_temp*x_temp + y_temp*y_temp + z_temp*z_temp,0.5);
									if(dis_min>=dis_temp)
									{
										dis_min = dis_temp;
									}
								}
							}
						}
						Block[nb].Geometry.Cell[i][j][k].dis_Wall = dis_min;
					}
					/*------------wall boundary cycle-----------------------------------*/
				}
			}
		}
		cout<<"	Completed!"<<endl;
	}

	cout<<"	Output Wall Distance..."<<endl;
	ofstream stream;
	stream.open("WallDistance.plt");
	for(nb=1;nb<=NB;nb++)
	{
		stream<<"VARIABLES = x	y z walldis"<<endl;
		stream<<"zone i="<<Block[nb].Geometry.IM-1<<" j="<<Block[nb].Geometry.JM-1<<" k="<<Block[nb].Geometry.KM-1<<endl;
		for(k=1;k<=Block[nb].Geometry.KM-1;k++)
		{
			for(j=1;j<=Block[nb].Geometry.JM-1;j++)
			{
				for(i=1;i<=Block[nb].Geometry.IM-1;i++)
				{
					stream<<Block[nb].Geometry.Node_coord[i][j][k].x<<" "<<Block[nb].Geometry.Node_coord[i][j][k].y
						  <<" "<<Block[nb].Geometry.Node_coord[i][j][k].z<<" "<<Block[nb].Geometry.Cell[i][j][k].dis_Wall<<endl;
				}
			}
		}
	}
	stream.close();

	
	stream.open("WallDistance.dat");
	for(nb=1;nb<=NB;nb++)
	{
		stream<<nb<<endl;
		stream<<Block[nb].Geometry.IM-1<<" "<<Block[nb].Geometry.JM-1<<" "<<Block[nb].Geometry.KM-1<<endl;
		for(k=1;k<=Block[nb].Geometry.KM-1;k++)
		{
			for(j=1;j<=Block[nb].Geometry.JM-1;j++)
			{
				for(i=1;i<=Block[nb].Geometry.IM-1;i++)
				{
					stream<<i<<" "<<j<<" "<<k<<" "<<Block[nb].Geometry.Cell[i][j][k].dis_Wall<<endl;;
				}
			}
		}
	}
	stream.close();
	cout<<"	Completed!"<<endl;	
}