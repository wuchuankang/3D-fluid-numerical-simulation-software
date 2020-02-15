

#include "FluidProps.h"
#include <iostream>
using namespace std;

void CFluidProps::PrimitiveVars()
{
	int i, j, k;
	REAL V2;

#pragma omp parallel for num_threads(NT) default(shared) private(i, j, k, V2)
	for(i=1;i<IM;i++)
	{
		for(j=1;j<JM;j++)
		{
			for(k=1;k<KM;k++)
			{
				Cell_pv[i][j][k].dens  = Cell_cv[i][j][k].dens;
				Cell_pv[i][j][k].uvel  = Cell_cv[i][j][k].xmom/Cell_cv[i][j][k].dens;
				Cell_pv[i][j][k].vvel  = Cell_cv[i][j][k].ymom/Cell_cv[i][j][k].dens;
				Cell_pv[i][j][k].wvel  = Cell_cv[i][j][k].zmom/Cell_cv[i][j][k].dens;
				V2 = Cell_pv[i][j][k].uvel*Cell_pv[i][j][k].uvel + Cell_pv[i][j][k].vvel*Cell_pv[i][j][k].vvel + Cell_pv[i][j][k].wvel*Cell_pv[i][j][k].wvel;
				Cell_pv[i][j][k].press = (Gamma-1.0f)*(Cell_cv[i][j][k].ener - Cell_cv[i][j][k].dens*V2/2.0f);
				Cell_pv[i][j][k].temp  = Cell_pv[i][j][k].press/Cell_pv[i][j][k].dens/Rgas;
			}
		}
	}
}