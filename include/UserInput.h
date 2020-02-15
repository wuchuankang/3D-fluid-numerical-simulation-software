// Definition of the class providing output of plot files
// and of the convergence history.
#ifndef USERINPUT_H
#define USERINPUT_H

#include "Solver.h"
#include "UserInput.h"

class CUserInput
{
public:
	static void ReadConfig(CSolver &Solver);
	static void ReadGrids(CSolver &Solver);
	static void ReadBndGrid(CSolver &Solver);
	static void ReadWallDistance(CMesh &Mesh);
	static void ReadSolution(CMesh &Mesh);
	static void NumberOfMultigrid(CSolver &Solver);
	static void NonDimensionalization(CSolver &Solver);
};

#endif