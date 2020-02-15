

#ifndef SOLVER_H
#define SOLVER_H

#include "defs.h"
#include "Mesh.h"

class CSolver 
{
public:
	bool reUseSolution;       /**< use of previous solution for restart (true, false) */
	bool reUseWallDistance;   /**< use of previous wall distance for restart (true, false) */
	bool UseMultigrid;
    int  MaxIter,             /**< max. number of iterations */
		 OutStep,             /**< number of iterations between solution dumps */
		 iter,                /**< actual iteration number */
		 Coarse_iter,		  /**< 初始情况下，在粗网格上迭代步数 */
		 NG;				  /**< number of multigrid **/
	REAL ConvTol;             /**< convergence criterion (2-norm of density change for */
	REAL CFL_res;
	
	REAL MassFlow_in,
		 MassFlow_out,
		 Efficiency,
		 MassFlow_ratio,
		 Press_ratio;

	TConsVars Res;
	REAL Res_ev;

	CMesh  *Mesh;

	//质量流量出口所需要的量
	static bool UseFlowMassImposed;
	REAL MassFlow_Imposed,
		 relax_press,
		 resc_mass_imp;
	int  iter_mass_imp;
	void MassFlowImposed();

	// function
	CSolver();
	~CSolver();

	void Solver();

private:
		
	void Prolong(CMesh &Mesh_1, CMesh &Mesh_2);		//fine grid prolong to coarse grid
	void Interpolate(CMesh &Mesh_1, CMesh &Mesh_2);								//coarse grid prolong to fine grid
	void ComputeDeltConsVars(CMesh &Mesh_2);

	CSolver(const CSolver &);					// copy constructor 
	CSolver & operator = (const CSolver &);		// assignment operator
};

#endif