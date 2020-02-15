// Structured_FV3D.cpp : Defines the entry point for the console application.
//


#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include "Solver.h"
#include "Output.h"
#include "UserInput.h"
#include <omp.h>
#include <iomanip> 
#include "Error.h"

using namespace std;
void Info();

int main(int argc, char *argv[])
{
	CSolver  Solver;
	string   str;
	int i;

	cout << endl
		<< " *************************************************" << endl
		<< " *                                               *" << endl
		<< " *       3-D FLOW ON STRUCTURED  GRIDS           *" << endl
		<< " *                                               *" << endl
		<< " *       (c) Liu zhaowei AND Wu chuankang        *" << endl
		<< " *                                               *" << endl
		<< " *       Version 1.0 from 04/01/2016             *" << endl
		<< " *                                               *" << endl
		<< " *************************************************" << endl << endl;

	if (argc == 1)
	{
		cout<<"read user configure ... "<<endl;
		CUserInput::ReadConfig(Solver);
		cout<<"complete!"<<endl<<endl;
		cout<<"get number of multigrids ... "<<endl;
		CUserInput::NumberOfMultigrid(Solver);
		cout<<"complete!"<<endl<<endl;
		cout<<"read grids ..."<<endl;
		CUserInput::ReadGrids(Solver);
		cout<<"complete!"<<endl<<endl;
		cout<<"output grid..."<<endl;
		COutput::OutputGrid(Solver.Mesh[1]);
		cout<<"complete!"<<endl<<endl;
		cout<<"read geometry bondary..."<<endl;
		CUserInput::ReadBndGrid(Solver);	
		cout<<"complete!"<<endl<<endl;
		cout<<"NonDimensionalization..."<<endl;
		CUserInput::NonDimensionalization(Solver);
		cout<<"complete!"<<endl<<endl;
	}
	else Info();

	cout<<"allocate memory ..."<<endl;
	for(i=1;i<=Solver.NG;i++)
		Solver.Mesh[i].AllocateMemory();
	if(Solver.NG == 2)
		Solver.Mesh[2].AllocateMemory_Pro_Inter_Vars();		//分配多重网格插值量
	cout<<"complete!"<<endl<<endl;

	cout<<"compute geometry variables ..."<<endl;
	for(i=1;i<=Solver.NG;i++)
		Solver.Mesh[i].Geometry();
	cout<<"complete!"<<endl<<endl;

	if(CFluidProps::turbulenceType==TurbulenceTypes::SA)
	{
		if(Solver.reUseWallDistance==false)
		{
			cout<<"compute wall distance ..."<<endl;
			Solver.Mesh[1].WallDistanceSA();	//no use the turbulence model at coarse mesh
		}
		else
		{
			cout<<"read wall distance ..."<<endl;
			CUserInput::ReadWallDistance(Solver.Mesh[1]);
		}
		cout<<"complete!"<<endl<<endl;
	}

	if(Solver.reUseSolution==false)
	{
		cout<<"initialization ..."<<endl;
		if(Solver.NG==1)
			Solver.Mesh[1].Initial();
		else
		{
			Solver.Mesh[2].Initial();
			if(CFluidProps::turbulenceType==TurbulenceTypes::SA)
				Solver.Mesh[1].Initial();
		}	
	}
	else
	{
		cout<<"read solution ..."<<endl;
		CUserInput::ReadSolution(Solver.Mesh[1]);
	}
	cout<<"complete!"<<endl<<endl;

	// print some statisticsSolver.OutStep
	cout << "------------------------------------------------------------"<<endl;
	cout << " Mesh file name        : " << CGeometry::fnameGrid <<endl;
	cout << " use multigrid method  : " << Solver.UseMultigrid <<endl;
	cout << " use FlowMassImposed   : " << Solver.UseFlowMassImposed <<endl;
	cout << " imposed Flow Mass     : " << Solver.MassFlow_Imposed <<endl;
	cout << " No. of multigrids     : " << Solver.NG <<endl;
	cout << " No. of blocks         : " << CMesh::NB <<endl;
	cout << " No. of dummy layers   : " << CGeometry::NDummy <<endl;
	cout << " total press inlet     : " << CFluidProps::Press_ref <<endl;
	cout << " total temp  inlet     : " << CFluidProps::Temp_ref <<endl;
	cout << " static press outlet   : " << CBndConds::Pout_static <<endl;
	cout << " Equation type         : " << CFluidProps::equationName <<endl;
	cout << " Timestep type         : " << CTimeDiscr::timestepName <<endl;
	cout << " Time discrete type    : " << CTimeDiscr::timeDiscrName <<endl;
	cout << " Space discrete type   : " << CSpaceDiscr::spaceDiscrName <<endl;
	cout << " reconstruct variable  : " << CSpaceDiscr::reconstructVarName <<endl;
	cout << " reconstruct type      : " << CSpaceDiscr::reconstructTypeName <<endl;
	cout << " Turbulence type       : " << CFluidProps::turbulenceName <<endl;
	cout << " use variables limiter : " << CFluidProps::limiter_vars <<endl;
	cout << " use IRS               : " << CTimeDiscr::IRS <<endl;
	cout << " CFL number            : " << CTimeDiscr::CFL <<endl;
	cout << " Maximum iteration     : " << Solver.MaxIter <<endl;
	cout << " coarse iterations     : " << Solver.Coarse_iter <<endl;
	cout << " save every iteration  : " << Solver.OutStep <<endl;
	cout << " Convergence criterion : " << Solver.ConvTol <<endl;
	cout << "------------------------------------------------------------"<<endl<<endl;

	REAL start_time, end_time;
	ofstream out("history.dat");
	if(CFluidProps::turbulenceType==TurbulenceTypes::SA)
		cout<<"iter      res    Mass_in	   Mass_out     PR     Efficiency      res_ev"<<endl;
	else if (CFluidProps::turbulenceType==TurbulenceTypes::SST)
		cout<<"iter      res    Mass_in	   Mass_out     PR     Efficiency      res_epslion      res_k"<<endl;
	else
		cout<<"iter      res    Mass_in	   Mass_out     PR     Efficiency"<<endl;

	for(Solver.iter=1;Solver.iter<=Solver.MaxIter;Solver.iter++)
	{	
		start_time = float(omp_get_wtime());

		Solver.Solver();

//----------------------------------------------------------------------------------
		//compute the parameter
		if(Solver.NG==1)
			COutput::OutputGeneralParameters(Solver.MassFlow_in, Solver.MassFlow_out, Solver.Efficiency, Solver.MassFlow_ratio, Solver.Press_ratio, Solver.Mesh[1]);
		else
		{
			if(Solver.iter <= Solver.Coarse_iter)
				COutput::OutputGeneralParameters(Solver.MassFlow_in, Solver.MassFlow_out, Solver.Efficiency, Solver.MassFlow_ratio, Solver.Press_ratio, Solver.Mesh[2]);
			else
				COutput::OutputGeneralParameters(Solver.MassFlow_in, Solver.MassFlow_out, Solver.Efficiency, Solver.MassFlow_ratio, Solver.Press_ratio, Solver.Mesh[1]);
		}
		// output		
		if(Solver.iter%10==0)
		{
			cout<<" "<<endl;
			cout<<"CFL="<<CTimeDiscr::CFL<<endl;
			cout<<"MassFlow_in="<<Solver.MassFlow_in<<"  MassFlow_out="<<Solver.MassFlow_out<<" MassFlowRatio="<<Solver.MassFlow_ratio<<endl;
			if(CFluidProps::turbulenceType==TurbulenceTypes::SA)
				cout<<"iter	  res	   Mass_in	   Mass_out	      PR         Efficiency      res_ev"<<endl;
			else
				cout<<"iter      res    Mass_in	   Mass_out     PR     Efficiency"<<endl;
		}

		if(CFluidProps::turbulenceType==TurbulenceTypes::SA)
		{
			out <<Solver.iter<<"  "<<Solver.Res.dens<<"    "<<Solver.MassFlow_in<<"         "<<Solver.MassFlow_out<<"      "<<Solver.Press_ratio<<"      "<<Solver.Efficiency<<"    "<<Solver.Res_ev<<endl;
			cout<<setprecision(5)<<Solver.iter<<"  "<<Solver.Res.dens<<"    "<<Solver.MassFlow_in<<"         "<<Solver.MassFlow_out<<"	    "<<Solver.Press_ratio<<"     "<<Solver.Efficiency<<"	 "<<Solver.Res_ev<<endl;
		}
		else if (CFluidProps::turbulenceType==TurbulenceTypes::SST)
		{

		}
		else
		{
			out <<Solver.iter<<"  "<<Solver.Res.dens<<"    "<<Solver.MassFlow_in<<"         "<<Solver.MassFlow_out<<"      "<<Solver.Press_ratio<<"     "<<Solver.Efficiency<<endl;
			cout<<setprecision(5)<<Solver.iter<<"  "<<Solver.Res.dens<<"    "<<Solver.MassFlow_in<<"         "<<Solver.MassFlow_out<<"      "<<Solver.Press_ratio<<"     "<<Solver.Efficiency<<endl;
		}

		if(Solver.UseFlowMassImposed == 1)
		{
			if(Solver.Res.dens <= Solver.resc_mass_imp)
			{
				if(Solver.iter%Solver.iter_mass_imp == 0)
				{
					Solver.MassFlowImposed();
					cout<<"Modified outlet pressure:"<<CBndConds::Pout<<endl;
				}
			}
		}
		//---------------------------------------------------------------------------------------------------------------------------------				
		end_time = float(omp_get_wtime());
		cout<<"Calculation time cost: "<<end_time-start_time<<"s"<<endl;
//---------------------------------------------------------------------------------------------------------------------------------
	
		// send message::divergence, if the calculation is blow up	
		if(_finite(Solver.Res.dens) == 0 || _finite(Solver.MassFlow_in) == 0 || _finite(Solver.MassFlow_out) == 0)
		{
			cout<<"Calculation is divergence !"<<endl;
			system("pause");
		}

		if(Solver.iter%Solver.OutStep==0 && Solver.iter > Solver.Coarse_iter)
		{
			COutput::OutputSolution(Solver.Mesh[1]);
			COutput::OutputRestartSolution(Solver.Mesh[1]);			
		}		
		if((Solver.Res.dens<=Solver.ConvTol) && (Solver.MassFlow_ratio<=0.003))	
		{
			cout<<"Calculation is convergence !"<<endl;
			break;
		}
	}	
	out.close();
	COutput::OutputSolution(Solver.Mesh[1]);
	COutput::OutputRestartSolution(Solver.Mesh[1]);	
	system("pause");
	return 0;
}


void Info()
{
	printf("Usage:\n");
	printf("struct3d <input_file>\n\n");

	exit(EXIT_FAILURE);
}