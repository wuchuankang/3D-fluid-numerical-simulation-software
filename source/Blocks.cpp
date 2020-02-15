

#include "Blocks.h"

CBlocks::CBlocks()
{
	Block_ID = 0;
	flag = nullptr;
}

CBlocks::~CBlocks()
{
	int i, j;
	if(flag != nullptr)
	{
		for(i=0;i<Geometry.IM;i++)
		{
			for(j=0;j<Geometry.JM;j++)
			{
				delete[] flag[i][j];
			}
			delete[] flag[i];
		}
		delete[] flag;
		flag = nullptr;
	}
}

void CBlocks::AllocateMemory()
{
	int i, j;
	flag = new int**[Geometry.IM];
	for(i=0;i<Geometry.IM;i++)
		flag[i] = new int*[Geometry.JM];
	for(i=0;i<Geometry.IM;i++)
	{
		for(j=0;j<Geometry.JM;j++)
		{
			flag[i][j] = new int[Geometry.KM];
		}
	}
}