
#ifndef OUTPUT_H
#define OUTPUT_H

#include "defs.h"
#include "Blocks.h"
#include "Mesh.h"

class COutput
{
public:
	static REAL MassFlowCaculation(int K, CBlocks &Block);
	static REAL TotalPressureSecCal(int K, CBlocks &Block);
	static REAL TotalTemperatureSecCal(int K, CBlocks &Block);
	static REAL EfficiencyCal(REAL Pt_in, REAL Tt_in, REAL Pt_out, REAL Tt_out);
	static void OutputSolution(CMesh &Mesh);
	static void OutputRestartSolution(CMesh &Mesh);
	static void	OutputGrid(CMesh &Mesh);

	static void	OutputGeneralParameters(REAL &MassFlow_in, REAL &MassFlow_out, REAL &Efficiency, REAL &MassFlow_ratio, REAL &Press_ratio, CMesh &Mesh);
};
#endif

