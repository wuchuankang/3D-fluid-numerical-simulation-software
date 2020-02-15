/***********************
1.determine the no of dummy layers
*************************/


#include <stdexcept>
#include <string>
#include <sstream>
#include <iostream>
#include "UserInput.h"
#include "fstringIO.h"
#include "Solver.h"
#include "defs.h"

using namespace std;

void CUserInput::ReadConfig( CSolver &Solver)
{
	string str;

	ostringstream os;

	ifstream stream("config.dat");

	str = ReadLine( stream );
	CGeometry::fnameGrid = ReadLine( stream );

	str = ReadLine( stream );
	str = ReadLine( stream );
	CGeometry::fbndGrid  = ReadLine( stream );

	// physics general
	str = ReadLine( stream );
	str = ReadLine( stream );
	str = ReadLine( stream );

	str = ReadLine( stream );
	if (str=="a" || str=="A")
		CMesh::TurboType = TurboTypes::Axial;		
	else if(str=="c" || str=="C")
		CMesh::TurboType = TurboTypes::Centrifugal;	
	else
		CMesh::TurboType = TurboTypes::General;	
	

	str = ReadLine( stream );
	if (str=="e" || str=="E"){
		 CFluidProps::equationType = EquationTypes::Euler;
		 CFluidProps::equationName = "Euler";}
	else{
		 CFluidProps::equationType = EquationTypes::NavierStokes;
		 CFluidProps::equationName = "NavierStokes";}
		
	str = ReadLine( stream );
	if (str=="u" || str=="U")
		CMesh::InitialType = InitialTypes::Uniform;		
	else
		CMesh::InitialType = InitialTypes::Interpolate;	

	str = ReadLine( stream );	CFluidProps::Rgas     = stof( str );
	str = ReadLine( stream );	CFluidProps::Gamma    = stof( str );
	str = ReadLine( stream );	CFluidProps::Pr_lam   = stof( str );
	str = ReadLine( stream );	CFluidProps::Pr_tur   = stof( str );	

	// physics internal
	str = ReadLine( stream );
	str = ReadLine( stream );
	str = ReadLine( stream );

    str = ReadLine( stream );   CBndConds::Ptin		= stof( str );
    str = ReadLine( stream );   CBndConds::Ttin		= stof( str );			
    str = ReadLine( stream );   CBndConds::Pout     = stof( str );		 CBndConds::Pout_static = CBndConds::Pout;
    str = ReadLine( stream );   CBndConds::P21ratio = stof( str );
	str = ReadLine( stream );	CFluidProps::Ma0    = stof( str );
	str = ReadLine( stream );	CGeometry::Len_ref  = stof( str );
	str = ReadLine( stream );	CBndConds::ev_ratio = stof( str );
	str = ReadLine( stream );	stringstream(str)>>CBndConds::InletAirAngle.x>>CBndConds::InletAirAngle.y>>CBndConds::InletAirAngle.z;												
	str = ReadLine( stream );	stringstream(str)>>CBndConds::OutletAirAngle.x>>CBndConds::OutletAirAngle.y>>CBndConds::OutletAirAngle.z;											 							
	str = ReadLine( stream );	Solver.UseFlowMassImposed = (str=="y" || str=="Y");
	str = ReadLine( stream );	stringstream(str)>>Solver.MassFlow_Imposed>>Solver.relax_press>>Solver.iter_mass_imp>>Solver.resc_mass_imp;

	// iteration control
	str = ReadLine( stream );
	str = ReadLine( stream );
	str = ReadLine( stream );
	str = ReadLine( stream );	Solver.MaxIter     = stoi( str );
	str = ReadLine( stream );	Solver.Coarse_iter = stoi( str );
	str = ReadLine( stream );	Solver.OutStep	   = stoi( str );
	str = ReadLine( stream );	Solver.ConvTol	   = stof( str );
	str = ReadLine( stream );	NT				   = stoi( str );
	str = ReadLine( stream );	Solver.UseMultigrid      = (str=="y" || str=="Y");
	str = ReadLine( stream );	Solver.reUseSolution     = (str=="y" || str=="Y");	
	str = ReadLine( stream );	Solver.reUseWallDistance = (str=="y" || str=="Y");

	//Numerical parameters
	str = ReadLine( stream );
	str = ReadLine( stream );
	str = ReadLine( stream );
	str = ReadLine( stream );   CTimeDiscr::CFL           = stof( str );
	str = ReadLine( stream );   CTimeDiscr::C_timestep    = stof( str );
	str = ReadLine( stream );   CTimeDiscr::relax_factor  = stof( str );
	str = ReadLine( stream );   CSpaceDiscr::epsEntr      = stof( str );
	str = ReadLine( stream );   CSpaceDiscr::k2           = stof( str );
	str = ReadLine( stream );   CSpaceDiscr::k4		      = stof( str );		CSpaceDiscr::k4 = 1.0f/CSpaceDiscr::k4;
	str = ReadLine( stream );   CTimeDiscr::IRS			  = (str=="y" || str=="Y");
	str = ReadLine( stream );	CFluidProps::limiter_vars = (str=="y" || str=="Y");
	str = ReadLine( stream );
	if (str=="l" || str=="L"){
		CTimeDiscr::timestepType = TimestepTypes::Local;
		CTimeDiscr::timestepName = "Local";
	}
	else{
		CTimeDiscr::timestepType = TimestepTypes::Global;
		CTimeDiscr::timestepName = "Global";
	}

	str = ReadLine( stream );	stringstream(str)>>CSpaceDiscr::ReconstructVar;
	if (CSpaceDiscr::ReconstructVar==100){
		CSpaceDiscr::reconstructVar = ReconstructVars::ConsVars;
		CSpaceDiscr::reconstructVarName = "ConsVars";
	}
	else if (CSpaceDiscr::ReconstructVar==200){
		CSpaceDiscr::reconstructVar = ReconstructVars::PrimVars;
		CSpaceDiscr::reconstructVarName = "PrimVars";
	}
	else {
		CSpaceDiscr::reconstructVar = ReconstructVars::CharactVars;
		CSpaceDiscr::reconstructVarName = "CharactVars";
	}

	str = ReadLine( stream );	stringstream(str)>>CSpaceDiscr::ReconstructScheme;
	if (CSpaceDiscr::ReconstructScheme==0){
		CSpaceDiscr::reconstructType = ReconstructTypes::None;
		CSpaceDiscr::reconstructTypeName = "None";
		CGeometry::NDummy = 3;}
	else if (CSpaceDiscr::ReconstructScheme==100){
		CSpaceDiscr::reconstructType = ReconstructTypes::NND;
		CSpaceDiscr::reconstructTypeName = "NND";
		CGeometry::NDummy = 2;}
	else if (CSpaceDiscr::ReconstructScheme==200){                     
		CSpaceDiscr::reconstructType = ReconstructTypes::MUSCL2u;
		CSpaceDiscr::reconstructTypeName = "MUSCL2u";
		CGeometry::NDummy = 2;}
	else if (CSpaceDiscr::ReconstructScheme==300){                   
		CSpaceDiscr::reconstructType = ReconstructTypes::MUSCL2c;
		CSpaceDiscr::reconstructTypeName = "MUSCL2c";
		CGeometry::NDummy = 2;}
	else if (CSpaceDiscr::ReconstructScheme==400){                    
		CSpaceDiscr::reconstructType = ReconstructTypes::MUSCL2o;
		CSpaceDiscr::reconstructTypeName = "MUSCL2o";
		CGeometry::NDummy = 2;}
	else {                    
		CSpaceDiscr::reconstructType = ReconstructTypes::MUSCL3;
		CSpaceDiscr::reconstructTypeName = "MUSCL3";
		CGeometry::NDummy = 2;}
		
	
	str = ReadLine( stream );	stringstream(str)>>CSpaceDiscr::DiscrScheme;
	if (CSpaceDiscr::DiscrScheme==100){
		CSpaceDiscr::spaceDiscrType = SpaceDiscrTypes::JST;
		CSpaceDiscr::spaceDiscrName = "JST";
		CGeometry::NDummy = 3;
	}
	else if (CSpaceDiscr::DiscrScheme==200){
		CSpaceDiscr::spaceDiscrType = SpaceDiscrTypes::StegerWarming;
		CSpaceDiscr::spaceDiscrName = "StegerWarming";
	}
	else if (CSpaceDiscr::DiscrScheme==300){
		CSpaceDiscr::spaceDiscrType = SpaceDiscrTypes::Roe;
		CSpaceDiscr::spaceDiscrName = "Roe";
	}
	else if(CSpaceDiscr::DiscrScheme==400){
		CSpaceDiscr::spaceDiscrType = SpaceDiscrTypes::L_F;
		CSpaceDiscr::spaceDiscrName = "L_F";
	}
	else if(CSpaceDiscr::DiscrScheme==500){
		CSpaceDiscr::spaceDiscrType = SpaceDiscrTypes::Van_Leer;
		CSpaceDiscr::spaceDiscrName = "Van_Leer";
	}
	else if(CSpaceDiscr::DiscrScheme==600){
		CSpaceDiscr::spaceDiscrType = SpaceDiscrTypes::AUSM_Plus;
		CSpaceDiscr::spaceDiscrName = "AUSM+";
	}
	else if(CSpaceDiscr::DiscrScheme==700){
		CSpaceDiscr::spaceDiscrType = SpaceDiscrTypes::AUSM_PW;
		CSpaceDiscr::spaceDiscrName = "AUSM_PW";
	}
	else{
		CSpaceDiscr::spaceDiscrType = SpaceDiscrTypes::HLLC;
		CSpaceDiscr::spaceDiscrName = "HLLC";
	}

	str = ReadLine( stream );	stringstream(str)>>CTimeDiscr::TimeDiscrScheme;
	if (CTimeDiscr::TimeDiscrScheme==100){
		CTimeDiscr::timeDiscrType = TimeDiscrTypes::Runge_Kutta3;
		CTimeDiscr::timeDiscrName = "Runge_Kutta3";
		CTimeDiscr::timeDiscrGeneralName = "Runge_Kutta";
	}
	else if(CTimeDiscr::TimeDiscrScheme==200){
		CTimeDiscr::timeDiscrType = TimeDiscrTypes::Runge_Kutta4;
		CTimeDiscr::timeDiscrName = "Runge_Kutta4";
		CTimeDiscr::timeDiscrGeneralName = "Runge_Kutta";
	}
	else if(CTimeDiscr::TimeDiscrScheme==300){
		CTimeDiscr::timeDiscrType = TimeDiscrTypes::Runge_Kutta5;
		CTimeDiscr::timeDiscrName = "Runge_Kutta5";
		CTimeDiscr::timeDiscrGeneralName = "Runge_Kutta";
	}
	else if(CTimeDiscr::TimeDiscrScheme==400){
		CTimeDiscr::timeDiscrType = TimeDiscrTypes::Runge_Kutta5_Hybrid;
		CTimeDiscr::timeDiscrName = "Runge_Kutta5_Hybrid";
		CTimeDiscr::timeDiscrGeneralName = "Runge_Kutta";
	}
	else{
		CTimeDiscr::timeDiscrType = TimeDiscrTypes::LU_SGS;
		CTimeDiscr::timeDiscrName = "LU_SGS";
		CTimeDiscr::timeDiscrGeneralName = "LU_SGS";
	}

	str = ReadLine( stream );	stringstream(str)>>CFluidProps::turbulenceScheme;
	if(CFluidProps::turbulenceScheme==0){
		CFluidProps::turbulenceType = TurbulenceTypes::None;
		CFluidProps::turbulenceName = "None";
	}
	else if (CFluidProps::turbulenceScheme==100){
		CFluidProps::turbulenceType = TurbulenceTypes::BL;
		CFluidProps::turbulenceName = "BL";
	}
	else if (CFluidProps::turbulenceScheme==200){
		CFluidProps::turbulenceType = TurbulenceTypes::SA;
		CFluidProps::turbulenceName = "SA";
	}
	else{
		CFluidProps::turbulenceType = TurbulenceTypes::SST;
		CFluidProps::turbulenceName = "SST";
	}

	stream.close();
}

void CUserInput::ReadGrids(CSolver &Solver)
{
	int IM, JM, KM;
	int i, j, k, nb, Num=0;
	string str;

	ifstream stream(CGeometry::fnameGrid);

	str = ReadLine(stream);
	stringstream(str)>>CMesh::NB;
	Solver.Mesh[1].Block = new CBlocks [Solver.Mesh[1].NB+1];

	for(nb=1;nb<=CMesh::NB;nb++)
	{
		str = ReadLine(stream);
		stringstream(str) >> Solver.Mesh[1].Block[nb].Block_ID;

		str = ReadLine(stream);
		stringstream(str) >> Solver.Mesh[1].Block[nb].Geometry.IM>>Solver.Mesh[1].Block[nb].Geometry.JM>>Solver.Mesh[1].Block[nb].Geometry.KM;
		IM = Solver.Mesh[1].Block[nb].Geometry.IM;    JM = Solver.Mesh[1].Block[nb].Geometry.JM;    KM = Solver.Mesh[1].Block[nb].Geometry.KM;
		Num= Num + (IM-1)*(JM-1)*(KM-1);
		cout<<"Mesh[1].Block= "<<Solver.Mesh[1].Block[nb].Block_ID<<"   IM = "<<IM<<" JM = "<<JM<<" KM = "<<KM<<endl;
		Solver.Mesh[1].Block[nb].Geometry.Node_coord = new TNode **[IM+1];
		for (i=0;i<=IM;i++)	
		{
			Solver.Mesh[1].Block[nb].Geometry.Node_coord[i] = new TNode *[JM+1];
		}
		for (i=0;i<=IM;i++)
		{
			for (j=0;j<=JM;j++)
			{
				Solver.Mesh[1].Block[nb].Geometry.Node_coord[i][j] = new TNode [KM+1];
			}
		}
		for (k=1;k<=KM;k++)
		{
			for (j=1;j<=JM;j++)
			{
				for (i=1;i<=IM;i++)
				{
					str = ReadLine(stream);
					stringstream(str) >> Solver.Mesh[1].Block[nb].Geometry.Node_coord[i][j][k].x >> Solver.Mesh[1].Block[nb].Geometry.Node_coord[i][j][k].y >> Solver.Mesh[1].Block[nb].Geometry.Node_coord[i][j][k].z;
				}
			}
		}
	}
	cout<<"Total fine grids: "<<Num<<endl;
	stream.close();

	
	if(Solver.NG==2)				//确定第2层网格的参数
	{
		Solver.Mesh[2].Block = new CBlocks [CMesh::NB+1];

		Num = 0;
		for(nb=1;nb<=CMesh::NB;nb++)
		{
			Solver.Mesh[2].Block[nb].Block_ID = Solver.Mesh[1].Block[nb].Block_ID;
			Solver.Mesh[2].Block[nb].Geometry.IM = (Solver.Mesh[1].Block[nb].Geometry.IM+1)/2;
			Solver.Mesh[2].Block[nb].Geometry.JM = (Solver.Mesh[1].Block[nb].Geometry.JM+1)/2;
			Solver.Mesh[2].Block[nb].Geometry.KM = (Solver.Mesh[1].Block[nb].Geometry.KM+1)/2;
			IM = Solver.Mesh[2].Block[nb].Geometry.IM;    JM = Solver.Mesh[2].Block[nb].Geometry.JM;    KM = Solver.Mesh[2].Block[nb].Geometry.KM;
			Num= Num + (IM-1)*(JM-1)*(KM-1);
			cout<<"Mesh[2].Block= "<<Solver.Mesh[1].Block[nb].Block_ID<<"   IM = "<<IM<<" JM = "<<JM<<" KM = "<<KM<<endl;
			Solver.Mesh[2].Block[nb].Geometry.Node_coord = new TNode **[IM+1];
			for (i=0;i<=IM;i++)	
			{
				Solver.Mesh[2].Block[nb].Geometry.Node_coord[i] = new TNode *[JM+1];
			}
			for (i=0;i<=IM;i++)
			{
				for (j=0;j<=JM;j++)
				{
					Solver.Mesh[2].Block[nb].Geometry.Node_coord[i][j] = new TNode [KM+1];
				}
			}
		}

		for(nb=1;nb<=CMesh::NB;nb++)
		{	
			for(k=1;k<=Solver.Mesh[2].Block[nb].Geometry.KM;k++)
			{
				for(j=1;j<=Solver.Mesh[2].Block[nb].Geometry.JM;j++)
				{
					for(i=1;i<=Solver.Mesh[2].Block[nb].Geometry.IM;i++)
					{
						Solver.Mesh[2].Block[nb].Geometry.Node_coord[i][j][k] = Solver.Mesh[1].Block[nb].Geometry.Node_coord[2*i-1][2*j-1][2*k-1];					
					}			
				}
			}
		}
		cout<<"Total coarse grids: "<<Num<<endl;
	}
}

void CUserInput::ReadBndGrid(CSolver &Solver)
{
	int i, nm, temp;
	string str;

	ifstream stream(CGeometry::fbndGrid);

	// inlet geometry-------------------------------------
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);			Solver.Mesh[1].BndGeom.nInlet = stoi(str);
	nm = Solver.Mesh[1].BndGeom.nInlet;		TInOutlet *inlet;	inlet  = new TInOutlet[nm+1];
	str = ReadLine(stream);			
	for (i=1;i<=nm;i++) {
	str = ReadLine(stream);	
	stringstream(str)>>inlet[i].InOutlet_ID>>inlet[i].Block_ID>>inlet[i].Face_type
				     >>inlet[i].is>>inlet[i].ie>>inlet[i].js>>inlet[i].je>>inlet[i].ks>>inlet[i].ke;
	}
	Solver.Mesh[1].BndGeom.Inlet = inlet;

	// outlet geometry-------------------------------------
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);			Solver.Mesh[1].BndGeom.nOutlet = stoi(str);
	nm	= Solver.Mesh[1].BndGeom.nOutlet;	TInOutlet *outlet;	outlet = new TInOutlet[nm+1];
	str = ReadLine(stream);	
	for (i=1;i<=nm;i++) {
	str = ReadLine(stream);	
	stringstream(str)>>outlet[i].InOutlet_ID>>outlet[i].Block_ID>>outlet[i].Face_type
					 >>outlet[i].is>>outlet[i].ie>>outlet[i].js>>outlet[i].je>>outlet[i].ks>>outlet[i].ke;
	}
	Solver.Mesh[1].BndGeom.Outlet = outlet;

	// wall geometry-------------------------------------
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);			Solver.Mesh[1].BndGeom.nWall = stoi(str);
	nm	= Solver.Mesh[1].BndGeom.nWall;		TWall *wall;	wall = new TWall[nm+1];
	str = ReadLine(stream);
	for (i=1;i<=nm;i++) {
		str = ReadLine(stream);	
		stringstream(str)>>wall[i].Wall_ID>>wall[i].Block_ID>>wall[i].Face_type>>wall[i].rotate
						 >>wall[i].is>>wall[i].ie>>wall[i].js>>wall[i].je>>wall[i].ks>>wall[i].ke;
	}
	Solver.Mesh[1].BndGeom.Wall = wall;

	// periodic geometry--------------------------------------
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);				Solver.Mesh[1].BndGeom.nPeriodic = stoi(str);
	nm	= Solver.Mesh[1].BndGeom.nPeriodic;		TPeriodic	*periodic;	periodic = new TPeriodic[nm+1];
	str = ReadLine(stream);
	for (i=1;i<=nm;i++) {
		str = ReadLine(stream);	
		stringstream(str)>>periodic[i].Periodic_ID>>periodic[i].orient>>periodic[i].Block_ID_1>>periodic[i].Face_type_1
						 >>periodic[i].is_1>>periodic[i].ie_1>>periodic[i].js_1>>periodic[i].je_1>>periodic[i].ks_1>>periodic[i].ke_1
						 >>periodic[i].Block_ID_2>>periodic[i].Face_type_2
						 >>periodic[i].is_2>>periodic[i].ie_2>>periodic[i].js_2>>periodic[i].je_2>>periodic[i].ks_2>>periodic[i].ke_2;
	}
	Solver.Mesh[1].BndGeom.Periodic	= periodic;

	// interior geometry---------------------------------------
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);				Solver.Mesh[1].BndGeom.nInterior = stoi(str);
	nm	= Solver.Mesh[1].BndGeom.nInterior;		TInterior	*interior;	interior = new TInterior[nm+1];
	str = ReadLine(stream);
	for (i=1;i<=nm;i++) {
		str = ReadLine(stream);	
		stringstream(str)>>interior[i].Interior_ID>>interior[i].orient>>interior[i].Block_ID_1>>interior[i].Face_type_1
						 >>interior[i].is_1>>interior[i].ie_1>>interior[i].js_1>>interior[i].je_1>>interior[i].ks_1>>interior[i].ke_1
						 >>interior[i].Block_ID_2>>interior[i].Face_type_2
						 >>interior[i].is_2>>interior[i].ie_2>>interior[i].js_2>>interior[i].je_2>>interior[i].ks_2>>interior[i].ke_2;
	}
	Solver.Mesh[1].BndGeom.Interior	= interior;

	// Mixing plane------------------------------------------------
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);				Solver.Mesh[1].BndGeom.nMixingPlane = stoi(str);
	nm	= Solver.Mesh[1].BndGeom.nMixingPlane;	TMixingPlane	*mixingPlane;	mixingPlane = new TMixingPlane[nm+1];
	str = ReadLine(stream);
	for (i=1;i<=nm;i++) {
		str = ReadLine(stream);	
		stringstream(str)>>mixingPlane[i].MixingPlane_ID>>mixingPlane[i].orient>>mixingPlane[i].Block_ID_1>>mixingPlane[i].Face_type_1
						 >>mixingPlane[i].is_1>>mixingPlane[i].ie_1>>mixingPlane[i].js_1>>mixingPlane[i].je_1>>mixingPlane[i].ks_1>>mixingPlane[i].ke_1
						 >>mixingPlane[i].Block_ID_2>>mixingPlane[i].Face_type_2
						 >>mixingPlane[i].is_2>>mixingPlane[i].ie_2>>mixingPlane[i].js_2>>mixingPlane[i].je_2>>mixingPlane[i].ks_2>>mixingPlane[i].ke_2;
	}
	Solver.Mesh[1].BndGeom.MixingPlane	= mixingPlane;
	// Block------------------------------------------------
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);				nm = stoi(str);
	str = ReadLine(stream);
	for (i=1;i<=nm;i++) 
	{
		str = ReadLine(stream);	
		stringstream(str)>>temp>>temp>>temp>>temp>>Solver.Mesh[1].Block[i].Geometry.NBlade>>Solver.Mesh[1].Block[i].FluidProp.rpm;
		Solver.Mesh[1].Block[i].Geometry.delta_th = 2.0f*PI/Solver.Mesh[1].Block[i].Geometry.NBlade;
		Solver.Mesh[1].Block[i].FluidProp.w_rot = 2.0f*PI*Solver.Mesh[1].Block[i].FluidProp.rpm/60.0f;	// attention:没有无量纲化
	}
	stream.close();	

	
	if(Solver.NG == 2)			//确定第二层网格的参数
	{
		for (i=1;i<=CMesh::NB;i++) 
		{
			Solver.Mesh[2].Block[i].Geometry.NBlade = Solver.Mesh[1].Block[i].Geometry.NBlade;
			Solver.Mesh[2].Block[i].Geometry.delta_th = Solver.Mesh[1].Block[i].Geometry.delta_th;
			Solver.Mesh[2].Block[i].FluidProp.w_rot = Solver.Mesh[1].Block[i].FluidProp.w_rot;
		}

		Solver.Mesh[2].BndGeom.Inlet = new TInOutlet [CBndGeom::nInlet + 1];
		for(i=1;i<=CBndGeom::nInlet;i++)
		{
			Solver.Mesh[2].BndGeom.Inlet[i].InOutlet_ID = Solver.Mesh[1].BndGeom.Inlet[i].InOutlet_ID;
			Solver.Mesh[2].BndGeom.Inlet[i].Block_ID	= Solver.Mesh[1].BndGeom.Inlet[i].Block_ID;
			Solver.Mesh[2].BndGeom.Inlet[i].Face_type	= Solver.Mesh[1].BndGeom.Inlet[i].Face_type;
			if(Solver.Mesh[1].BndGeom.Inlet[i].is == Solver.Mesh[1].BndGeom.Inlet[i].ie)
			{
				Solver.Mesh[2].BndGeom.Inlet[i].is		  = (Solver.Mesh[1].BndGeom.Inlet[i].is+1)/2;
				Solver.Mesh[2].BndGeom.Inlet[i].ie		  = (Solver.Mesh[1].BndGeom.Inlet[i].ie+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.Inlet[i].is		  = (Solver.Mesh[1].BndGeom.Inlet[i].is+1)/2;
				Solver.Mesh[2].BndGeom.Inlet[i].ie		  =  Solver.Mesh[1].BndGeom.Inlet[i].ie/2;
			}
			if(Solver.Mesh[1].BndGeom.Inlet[i].js == Solver.Mesh[1].BndGeom.Inlet[i].je)
			{
				Solver.Mesh[2].BndGeom.Inlet[i].js		  = (Solver.Mesh[1].BndGeom.Inlet[i].js+1)/2;
				Solver.Mesh[2].BndGeom.Inlet[i].je		  = (Solver.Mesh[1].BndGeom.Inlet[i].je+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.Inlet[i].js		  = (Solver.Mesh[1].BndGeom.Inlet[i].js+1)/2;
				Solver.Mesh[2].BndGeom.Inlet[i].je		  =  Solver.Mesh[1].BndGeom.Inlet[i].je/2;
			}
			if(Solver.Mesh[1].BndGeom.Inlet[i].ks == Solver.Mesh[1].BndGeom.Inlet[i].ke)
			{
				Solver.Mesh[2].BndGeom.Inlet[i].ks		  = (Solver.Mesh[1].BndGeom.Inlet[i].ks+1)/2;
				Solver.Mesh[2].BndGeom.Inlet[i].ke		  = (Solver.Mesh[1].BndGeom.Inlet[i].ke+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.Inlet[i].ks		  = (Solver.Mesh[1].BndGeom.Inlet[i].ks+1)/2;
				Solver.Mesh[2].BndGeom.Inlet[i].ke		  =  Solver.Mesh[1].BndGeom.Inlet[i].ke/2;
			}
		}
				
		Solver.Mesh[2].BndGeom.Outlet = new TInOutlet [CBndGeom::nOutlet + 1];
		for(i=1;i<=CBndGeom::nOutlet;i++)
		{
			Solver.Mesh[2].BndGeom.Outlet[i].InOutlet_ID = Solver.Mesh[1].BndGeom.Outlet[i].InOutlet_ID;
			Solver.Mesh[2].BndGeom.Outlet[i].Block_ID	 = Solver.Mesh[1].BndGeom.Outlet[i].Block_ID;
			Solver.Mesh[2].BndGeom.Outlet[i].Face_type	 = Solver.Mesh[1].BndGeom.Outlet[i].Face_type;

			if(Solver.Mesh[1].BndGeom.Outlet[i].is == Solver.Mesh[1].BndGeom.Outlet[i].ie)
			{
				Solver.Mesh[2].BndGeom.Outlet[i].is		  = (Solver.Mesh[1].BndGeom.Outlet[i].is+1)/2;
				Solver.Mesh[2].BndGeom.Outlet[i].ie		  = (Solver.Mesh[1].BndGeom.Outlet[i].ie+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.Outlet[i].is		  = (Solver.Mesh[1].BndGeom.Outlet[i].is+1)/2;
				Solver.Mesh[2].BndGeom.Outlet[i].ie		  =  Solver.Mesh[1].BndGeom.Outlet[i].ie/2;
			}
			if(Solver.Mesh[1].BndGeom.Outlet[i].js == Solver.Mesh[1].BndGeom.Outlet[i].je)
			{
				Solver.Mesh[2].BndGeom.Outlet[i].js		  = (Solver.Mesh[1].BndGeom.Outlet[i].js+1)/2;
				Solver.Mesh[2].BndGeom.Outlet[i].je		  = (Solver.Mesh[1].BndGeom.Outlet[i].je+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.Outlet[i].js		  = (Solver.Mesh[1].BndGeom.Outlet[i].js+1)/2;
				Solver.Mesh[2].BndGeom.Outlet[i].je		  =  Solver.Mesh[1].BndGeom.Outlet[i].je/2;
			}
			if(Solver.Mesh[1].BndGeom.Outlet[i].ks == Solver.Mesh[1].BndGeom.Outlet[i].ke)
			{
				Solver.Mesh[2].BndGeom.Outlet[i].ks		  = (Solver.Mesh[1].BndGeom.Outlet[i].ks+1)/2;
				Solver.Mesh[2].BndGeom.Outlet[i].ke		  = (Solver.Mesh[1].BndGeom.Outlet[i].ke+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.Outlet[i].ks		  = (Solver.Mesh[1].BndGeom.Outlet[i].ks+1)/2;
				Solver.Mesh[2].BndGeom.Outlet[i].ke		  =  Solver.Mesh[1].BndGeom.Outlet[i].ke/2;
			}
		}
			
		Solver.Mesh[2].BndGeom.Wall = new TWall [CBndGeom::nWall + 1];
		for(i=1;i<=CBndGeom::nWall;i++)
		{
			Solver.Mesh[2].BndGeom.Wall[i].Wall_ID   = Solver.Mesh[1].BndGeom.Wall[i].Wall_ID;
			Solver.Mesh[2].BndGeom.Wall[i].Block_ID	 = Solver.Mesh[1].BndGeom.Wall[i].Block_ID;
			Solver.Mesh[2].BndGeom.Wall[i].Face_type = Solver.Mesh[1].BndGeom.Wall[i].Face_type;
			Solver.Mesh[2].BndGeom.Wall[i].rotate    = Solver.Mesh[1].BndGeom.Wall[i].rotate;

			if(Solver.Mesh[1].BndGeom.Wall[i].is == Solver.Mesh[1].BndGeom.Wall[i].ie)
			{
				Solver.Mesh[2].BndGeom.Wall[i].is		  = (Solver.Mesh[1].BndGeom.Wall[i].is+1)/2;
				Solver.Mesh[2].BndGeom.Wall[i].ie		  = (Solver.Mesh[1].BndGeom.Wall[i].ie+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.Wall[i].is		  = (Solver.Mesh[1].BndGeom.Wall[i].is+1)/2;
				Solver.Mesh[2].BndGeom.Wall[i].ie		  =  Solver.Mesh[1].BndGeom.Wall[i].ie/2;
			}
			if(Solver.Mesh[1].BndGeom.Wall[i].js == Solver.Mesh[1].BndGeom.Wall[i].je)
			{
				Solver.Mesh[2].BndGeom.Wall[i].js		  = (Solver.Mesh[1].BndGeom.Wall[i].js+1)/2;
				Solver.Mesh[2].BndGeom.Wall[i].je		  = (Solver.Mesh[1].BndGeom.Wall[i].je+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.Wall[i].js		  = (Solver.Mesh[1].BndGeom.Wall[i].js+1)/2;
				Solver.Mesh[2].BndGeom.Wall[i].je		  =  Solver.Mesh[1].BndGeom.Wall[i].je/2;
			}
			if(Solver.Mesh[1].BndGeom.Wall[i].ks == Solver.Mesh[1].BndGeom.Wall[i].ke)
			{
				Solver.Mesh[2].BndGeom.Wall[i].ks		  = (Solver.Mesh[1].BndGeom.Wall[i].ks+1)/2;
				Solver.Mesh[2].BndGeom.Wall[i].ke		  = (Solver.Mesh[1].BndGeom.Wall[i].ke+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.Wall[i].ks		  = (Solver.Mesh[1].BndGeom.Wall[i].ks+1)/2;
				Solver.Mesh[2].BndGeom.Wall[i].ke		  =  Solver.Mesh[1].BndGeom.Wall[i].ke/2;
			}
		}
				
		Solver.Mesh[2].BndGeom.Periodic = new TPeriodic[CBndGeom::nPeriodic + 1];
		for(i=1;i<=CBndGeom::nPeriodic;i++)
		{
			Solver.Mesh[2].BndGeom.Periodic[i].Periodic_ID  = Solver.Mesh[1].BndGeom.Periodic[i].Periodic_ID;
			Solver.Mesh[2].BndGeom.Periodic[i].Block_ID_1	= Solver.Mesh[1].BndGeom.Periodic[i].Block_ID_1;
			Solver.Mesh[2].BndGeom.Periodic[i].Face_type_1	= Solver.Mesh[1].BndGeom.Periodic[i].Face_type_1;
			Solver.Mesh[2].BndGeom.Periodic[i].orient		= Solver.Mesh[1].BndGeom.Periodic[i].orient;
			Solver.Mesh[2].BndGeom.Periodic[i].Block_ID_2	= Solver.Mesh[1].BndGeom.Periodic[i].Block_ID_2;
			Solver.Mesh[2].BndGeom.Periodic[i].Face_type_2	= Solver.Mesh[1].BndGeom.Periodic[i].Face_type_2;

			if(Solver.Mesh[1].BndGeom.Periodic[i].is_1 == Solver.Mesh[1].BndGeom.Periodic[i].ie_1)
			{
				Solver.Mesh[2].BndGeom.Periodic[i].is_1		  = (Solver.Mesh[1].BndGeom.Periodic[i].is_1+1)/2;
				Solver.Mesh[2].BndGeom.Periodic[i].ie_1		  = (Solver.Mesh[1].BndGeom.Periodic[i].ie_1+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.Periodic[i].is_1		  = (Solver.Mesh[1].BndGeom.Periodic[i].is_1+1)/2;
				Solver.Mesh[2].BndGeom.Periodic[i].ie_1		  =  Solver.Mesh[1].BndGeom.Periodic[i].ie_1/2;
			}
			if(Solver.Mesh[1].BndGeom.Periodic[i].js_1 == Solver.Mesh[1].BndGeom.Periodic[i].je_1)
			{
				Solver.Mesh[2].BndGeom.Periodic[i].js_1		  = (Solver.Mesh[1].BndGeom.Periodic[i].js_1+1)/2;
				Solver.Mesh[2].BndGeom.Periodic[i].je_1		  = (Solver.Mesh[1].BndGeom.Periodic[i].je_1+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.Periodic[i].js_1		  = (Solver.Mesh[1].BndGeom.Periodic[i].js_1+1)/2;
				Solver.Mesh[2].BndGeom.Periodic[i].je_1		  =  Solver.Mesh[1].BndGeom.Periodic[i].je_1/2;
			}
			if(Solver.Mesh[1].BndGeom.Periodic[i].ks_1 == Solver.Mesh[1].BndGeom.Periodic[i].ke_1)
			{
				Solver.Mesh[2].BndGeom.Periodic[i].ks_1		  = (Solver.Mesh[1].BndGeom.Periodic[i].ks_1+1)/2;
				Solver.Mesh[2].BndGeom.Periodic[i].ke_1		  = (Solver.Mesh[1].BndGeom.Periodic[i].ke_1+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.Periodic[i].ks_1		  = (Solver.Mesh[1].BndGeom.Periodic[i].ks_1+1)/2;
				Solver.Mesh[2].BndGeom.Periodic[i].ke_1		  =  Solver.Mesh[1].BndGeom.Periodic[i].ke_1/2;
			}


			if(Solver.Mesh[1].BndGeom.Periodic[i].is_2 == Solver.Mesh[1].BndGeom.Periodic[i].ie_2)
			{
				Solver.Mesh[2].BndGeom.Periodic[i].is_2		  = (Solver.Mesh[1].BndGeom.Periodic[i].is_2+1)/2;
				Solver.Mesh[2].BndGeom.Periodic[i].ie_2		  = (Solver.Mesh[1].BndGeom.Periodic[i].ie_2+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.Periodic[i].is_2		  = (Solver.Mesh[1].BndGeom.Periodic[i].is_2+1)/2;
				Solver.Mesh[2].BndGeom.Periodic[i].ie_2		  =  Solver.Mesh[1].BndGeom.Periodic[i].ie_2/2;
			}
			if(Solver.Mesh[1].BndGeom.Periodic[i].js_2 == Solver.Mesh[1].BndGeom.Periodic[i].je_2)
			{
				Solver.Mesh[2].BndGeom.Periodic[i].js_2		  = (Solver.Mesh[1].BndGeom.Periodic[i].js_2+1)/2;
				Solver.Mesh[2].BndGeom.Periodic[i].je_2		  = (Solver.Mesh[1].BndGeom.Periodic[i].je_2+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.Periodic[i].js_2		  = (Solver.Mesh[1].BndGeom.Periodic[i].js_2+1)/2;
				Solver.Mesh[2].BndGeom.Periodic[i].je_2		  =  Solver.Mesh[1].BndGeom.Periodic[i].je_2/2;
			}
			if(Solver.Mesh[1].BndGeom.Periodic[i].ks_2 == Solver.Mesh[1].BndGeom.Periodic[i].ke_2)
			{
				Solver.Mesh[2].BndGeom.Periodic[i].ks_2		  = (Solver.Mesh[1].BndGeom.Periodic[i].ks_2+1)/2;
				Solver.Mesh[2].BndGeom.Periodic[i].ke_2		  = (Solver.Mesh[1].BndGeom.Periodic[i].ke_2+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.Periodic[i].ks_2		  = (Solver.Mesh[1].BndGeom.Periodic[i].ks_2+1)/2;
				Solver.Mesh[2].BndGeom.Periodic[i].ke_2		  =  Solver.Mesh[1].BndGeom.Periodic[i].ke_2/2;
			}
		}

		Solver.Mesh[2].BndGeom.Interior = new TInterior[CBndGeom::nInterior + 1];
		for(i=1;i<=CBndGeom::nInterior;i++)
		{
			Solver.Mesh[2].BndGeom.Interior[i].Interior_ID  = Solver.Mesh[1].BndGeom.Interior[i].Interior_ID;
			Solver.Mesh[2].BndGeom.Interior[i].Block_ID_1	= Solver.Mesh[1].BndGeom.Interior[i].Block_ID_1;
			Solver.Mesh[2].BndGeom.Interior[i].Face_type_1	= Solver.Mesh[1].BndGeom.Interior[i].Face_type_1;
			Solver.Mesh[2].BndGeom.Interior[i].orient		= Solver.Mesh[1].BndGeom.Interior[i].orient;
			Solver.Mesh[2].BndGeom.Interior[i].Block_ID_2	= Solver.Mesh[1].BndGeom.Interior[i].Block_ID_2;
			Solver.Mesh[2].BndGeom.Interior[i].Face_type_2	= Solver.Mesh[1].BndGeom.Interior[i].Face_type_2;

			if(Solver.Mesh[1].BndGeom.Interior[i].is_1 == Solver.Mesh[1].BndGeom.Interior[i].ie_1)
			{
				Solver.Mesh[2].BndGeom.Interior[i].is_1		  = (Solver.Mesh[1].BndGeom.Interior[i].is_1+1)/2;
				Solver.Mesh[2].BndGeom.Interior[i].ie_1		  = (Solver.Mesh[1].BndGeom.Interior[i].ie_1+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.Interior[i].is_1		  = (Solver.Mesh[1].BndGeom.Interior[i].is_1+1)/2;
				Solver.Mesh[2].BndGeom.Interior[i].ie_1		  =  Solver.Mesh[1].BndGeom.Interior[i].ie_1/2;
			}
			if(Solver.Mesh[1].BndGeom.Interior[i].js_1 == Solver.Mesh[1].BndGeom.Interior[i].je_1)
			{
				Solver.Mesh[2].BndGeom.Interior[i].js_1		  = (Solver.Mesh[1].BndGeom.Interior[i].js_1+1)/2;
				Solver.Mesh[2].BndGeom.Interior[i].je_1		  = (Solver.Mesh[1].BndGeom.Interior[i].je_1+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.Interior[i].js_1		  = (Solver.Mesh[1].BndGeom.Interior[i].js_1+1)/2;
				Solver.Mesh[2].BndGeom.Interior[i].je_1		  =  Solver.Mesh[1].BndGeom.Interior[i].je_1/2;
			}
			if(Solver.Mesh[1].BndGeom.Interior[i].ks_1 == Solver.Mesh[1].BndGeom.Interior[i].ke_1)
			{
				Solver.Mesh[2].BndGeom.Interior[i].ks_1		  = (Solver.Mesh[1].BndGeom.Interior[i].ks_1+1)/2;
				Solver.Mesh[2].BndGeom.Interior[i].ke_1		  = (Solver.Mesh[1].BndGeom.Interior[i].ke_1+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.Interior[i].ks_1		  = (Solver.Mesh[1].BndGeom.Interior[i].ks_1+1)/2;
				Solver.Mesh[2].BndGeom.Interior[i].ke_1		  =  Solver.Mesh[1].BndGeom.Interior[i].ke_1/2;
			}


			if(Solver.Mesh[1].BndGeom.Interior[i].is_2 == Solver.Mesh[1].BndGeom.Interior[i].ie_2)
			{
				Solver.Mesh[2].BndGeom.Interior[i].is_2		  = (Solver.Mesh[1].BndGeom.Interior[i].is_2+1)/2;
				Solver.Mesh[2].BndGeom.Interior[i].ie_2		  = (Solver.Mesh[1].BndGeom.Interior[i].ie_2+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.Interior[i].is_2		  = (Solver.Mesh[1].BndGeom.Interior[i].is_2+1)/2;
				Solver.Mesh[2].BndGeom.Interior[i].ie_2		  =  Solver.Mesh[1].BndGeom.Interior[i].ie_2/2;
			}
			if(Solver.Mesh[1].BndGeom.Interior[i].js_2 == Solver.Mesh[1].BndGeom.Interior[i].je_2)
			{
				Solver.Mesh[2].BndGeom.Interior[i].js_2		  = (Solver.Mesh[1].BndGeom.Interior[i].js_2+1)/2;
				Solver.Mesh[2].BndGeom.Interior[i].je_2		  = (Solver.Mesh[1].BndGeom.Interior[i].je_2+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.Interior[i].js_2		  = (Solver.Mesh[1].BndGeom.Interior[i].js_2+1)/2;
				Solver.Mesh[2].BndGeom.Interior[i].je_2		  =  Solver.Mesh[1].BndGeom.Interior[i].je_2/2;
			}
			if(Solver.Mesh[1].BndGeom.Interior[i].ks_2 == Solver.Mesh[1].BndGeom.Interior[i].ke_2)
			{
				Solver.Mesh[2].BndGeom.Interior[i].ks_2		  = (Solver.Mesh[1].BndGeom.Interior[i].ks_2+1)/2;
				Solver.Mesh[2].BndGeom.Interior[i].ke_2		  = (Solver.Mesh[1].BndGeom.Interior[i].ke_2+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.Interior[i].ks_2		  = (Solver.Mesh[1].BndGeom.Interior[i].ks_2+1)/2;
				Solver.Mesh[2].BndGeom.Interior[i].ke_2		  =  Solver.Mesh[1].BndGeom.Interior[i].ke_2/2;
			}
		}
				
		Solver.Mesh[2].BndGeom.MixingPlane = new TMixingPlane[CBndGeom::nMixingPlane + 1];
		for(i=1;i<=CBndGeom::nMixingPlane;i++)
		{
			Solver.Mesh[2].BndGeom.MixingPlane[i].MixingPlane_ID   = Solver.Mesh[1].BndGeom.MixingPlane[i].MixingPlane_ID;
			Solver.Mesh[2].BndGeom.MixingPlane[i].Block_ID_1	= Solver.Mesh[1].BndGeom.MixingPlane[i].Block_ID_1;
			Solver.Mesh[2].BndGeom.MixingPlane[i].Face_type_1	= Solver.Mesh[1].BndGeom.MixingPlane[i].Face_type_1;
			Solver.Mesh[2].BndGeom.MixingPlane[i].orient		= Solver.Mesh[1].BndGeom.MixingPlane[i].orient;
			Solver.Mesh[2].BndGeom.MixingPlane[i].Block_ID_2	= Solver.Mesh[1].BndGeom.MixingPlane[i].Block_ID_2;
			Solver.Mesh[2].BndGeom.MixingPlane[i].Face_type_2	= Solver.Mesh[1].BndGeom.MixingPlane[i].Face_type_2;

			if(Solver.Mesh[1].BndGeom.MixingPlane[i].is_1 == Solver.Mesh[1].BndGeom.MixingPlane[i].ie_1)
			{
				Solver.Mesh[2].BndGeom.MixingPlane[i].is_1		  = (Solver.Mesh[1].BndGeom.MixingPlane[i].is_1+1)/2;
				Solver.Mesh[2].BndGeom.MixingPlane[i].ie_1		  = (Solver.Mesh[1].BndGeom.MixingPlane[i].ie_1+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.MixingPlane[i].is_1		  = (Solver.Mesh[1].BndGeom.MixingPlane[i].is_1+1)/2;
				Solver.Mesh[2].BndGeom.MixingPlane[i].ie_1		  =  Solver.Mesh[1].BndGeom.MixingPlane[i].ie_1/2;
			}
			if(Solver.Mesh[1].BndGeom.MixingPlane[i].js_1 == Solver.Mesh[1].BndGeom.MixingPlane[i].je_1)
			{
				Solver.Mesh[2].BndGeom.MixingPlane[i].js_1		  = (Solver.Mesh[1].BndGeom.MixingPlane[i].js_1+1)/2;
				Solver.Mesh[2].BndGeom.MixingPlane[i].je_1		  = (Solver.Mesh[1].BndGeom.MixingPlane[i].je_1+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.MixingPlane[i].js_1		  = (Solver.Mesh[1].BndGeom.MixingPlane[i].js_1+1)/2;
				Solver.Mesh[2].BndGeom.MixingPlane[i].je_1		  =  Solver.Mesh[1].BndGeom.MixingPlane[i].je_1/2;
			}
			if(Solver.Mesh[1].BndGeom.MixingPlane[i].ks_1 == Solver.Mesh[1].BndGeom.MixingPlane[i].ke_1)
			{
				Solver.Mesh[2].BndGeom.MixingPlane[i].ks_1		  = (Solver.Mesh[1].BndGeom.MixingPlane[i].ks_1+1)/2;
				Solver.Mesh[2].BndGeom.MixingPlane[i].ke_1		  = (Solver.Mesh[1].BndGeom.MixingPlane[i].ke_1+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.MixingPlane[i].ks_1		  = (Solver.Mesh[1].BndGeom.MixingPlane[i].ks_1+1)/2;
				Solver.Mesh[2].BndGeom.MixingPlane[i].ke_1		  =  Solver.Mesh[1].BndGeom.MixingPlane[i].ke_1/2;
			}


			if(Solver.Mesh[1].BndGeom.MixingPlane[i].is_2 == Solver.Mesh[1].BndGeom.MixingPlane[i].ie_2)
			{
				Solver.Mesh[2].BndGeom.MixingPlane[i].is_2		  = (Solver.Mesh[1].BndGeom.MixingPlane[i].is_2+1)/2;
				Solver.Mesh[2].BndGeom.MixingPlane[i].ie_2		  = (Solver.Mesh[1].BndGeom.MixingPlane[i].ie_2+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.MixingPlane[i].is_2		  = (Solver.Mesh[1].BndGeom.MixingPlane[i].is_2+1)/2;
				Solver.Mesh[2].BndGeom.MixingPlane[i].ie_2		  =  Solver.Mesh[1].BndGeom.MixingPlane[i].ie_2/2;
			}
			if(Solver.Mesh[1].BndGeom.MixingPlane[i].js_2 == Solver.Mesh[1].BndGeom.MixingPlane[i].je_2)
			{
				Solver.Mesh[2].BndGeom.MixingPlane[i].js_2		  = (Solver.Mesh[1].BndGeom.MixingPlane[i].js_2+1)/2;
				Solver.Mesh[2].BndGeom.MixingPlane[i].je_2		  = (Solver.Mesh[1].BndGeom.MixingPlane[i].je_2+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.MixingPlane[i].js_2		  = (Solver.Mesh[1].BndGeom.MixingPlane[i].js_2+1)/2;
				Solver.Mesh[2].BndGeom.MixingPlane[i].je_2		  =  Solver.Mesh[1].BndGeom.MixingPlane[i].je_2/2;
			}
			if(Solver.Mesh[1].BndGeom.MixingPlane[i].ks_2 == Solver.Mesh[1].BndGeom.MixingPlane[i].ke_2)
			{
				Solver.Mesh[2].BndGeom.MixingPlane[i].ks_2		  = (Solver.Mesh[1].BndGeom.MixingPlane[i].ks_2+1)/2;
				Solver.Mesh[2].BndGeom.MixingPlane[i].ke_2		  = (Solver.Mesh[1].BndGeom.MixingPlane[i].ke_2+1)/2;
			}
			else
			{
				Solver.Mesh[2].BndGeom.MixingPlane[i].ks_2		  = (Solver.Mesh[1].BndGeom.MixingPlane[i].ks_2+1)/2;
				Solver.Mesh[2].BndGeom.MixingPlane[i].ke_2		  =  Solver.Mesh[1].BndGeom.MixingPlane[i].ke_2/2;
			}
		}
	}
}

void CUserInput::NumberOfMultigrid(CSolver &Solver)
{
	int i, j, IM, JM, KM, nn, Num, temp;
	int is, ie, js, je, ks, ke;
	int nInlet, nOutlet, nWall, nPeriodic, nInterior, nMixingPlane, NB;
	int NG=4;			//number of multigrid;
	string str;

	ifstream stream(CGeometry::fbndGrid);

	// inlet geometry-------------------------------------
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);			nInlet = stoi(str);
	str = ReadLine(stream);
	for (i=1;i<=nInlet;i++) 
	{
		str = ReadLine(stream);	
		stringstream(str)>>temp>>temp>>temp>>is>>ie>>js>>je>>ks>>ke;

		if(is != ie)
		{
			is = is;
			ie = ie + 1;
		}
		if(js != je)
		{
			js = js;
			je = je + 1;
		}
		if(ks != ke)
		{
			ks = ks;
			ke = ke + 1;
		}

		nn = 2;
		Num= 1;
		if((is-1)%nn ==0 && (ie-1)%nn ==0 && (js-1)%nn ==0 && (je-1)%nn ==0 && (ks-1)%nn ==0 && (ke-1)%nn ==0)
		{
			nn = nn*2;
			Num= Num+1;
		}
		NG = MIN(Num, NG);			     
	}


	// outlet geometry-------------------------------------
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);			nOutlet = stoi(str);
	str = ReadLine(stream);
	for (i=1;i<=nOutlet;i++) 
	{
		str = ReadLine(stream);	
		stringstream(str)>>temp>>temp>>temp>>is>>ie>>js>>je>>ks>>ke;

		if(is != ie)
		{
			is = is;
			ie = ie + 1;
		}
		if(js != je)
		{
			js = js;
			je = je + 1;
		}
		if(ks != ke)
		{
			ks = ks;
			ke = ke + 1;
		}

		nn = 2;
		Num= 1;
		if((is-1)%nn ==0 && (ie-1)%nn ==0 && (js-1)%nn ==0 && (je-1)%nn ==0 && (ks-1)%nn ==0 && (ke-1)%nn ==0)
		{
			nn = nn*2;
			Num= Num+1;
		}
		NG = MIN(Num, NG);			     
	}

	// wall geometry-------------------------------------
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);			nWall = stoi(str);
	str = ReadLine(stream);
	for (i=1;i<=nWall;i++) 
	{
		str = ReadLine(stream);	
		stringstream(str)>>temp>>temp>>temp>>temp>>is>>ie>>js>>je>>ks>>ke;

		if(is != ie)
		{
			is = is;
			ie = ie + 1;
		}
		if(js != je)
		{
			js = js;
			je = je + 1;
		}
		if(ks != ke)
		{
			ks = ks;
			ke = ke + 1;
		}

		nn = 2;
		Num= 1;
		if((is-1)%nn ==0 && (ie-1)%nn ==0 && (js-1)%nn ==0 && (je-1)%nn ==0 && (ks-1)%nn ==0 && (ke-1)%nn ==0)
		{
			nn = nn*2;
			Num= Num+1;
		}
		NG = MIN(Num, NG);			     
	}

	// periodic geometry--------------------------------------
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);			nPeriodic = stoi(str);
	str = ReadLine(stream);
	for (i=1;i<=nPeriodic;i++) 
	{
		str = ReadLine(stream);	
		stringstream(str)>>temp>>temp>>temp>>temp>>is>>ie>>js>>je>>ks>>ke>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp;
					
		if(is != ie)
		{
			is = is;
			ie = ie + 1;
		}
		if(js != je)
		{
			js = js;
			je = je + 1;
		}
		if(ks != ke)
		{
			ks = ks;
			ke = ke + 1;
		}

		nn = 2;
		Num= 1;
		if((is-1)%nn ==0 && (ie-1)%nn ==0 && (js-1)%nn ==0 && (je-1)%nn ==0 && (ks-1)%nn ==0 && (ke-1)%nn ==0)
		{
			nn = nn*2;
			Num= Num+1;
		}
		NG = MIN(Num, NG);												 
	}

	// interior geometry---------------------------------------
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);				nInterior = stoi(str);
	str = ReadLine(stream);
	for (i=1;i<=nInterior;i++) 
	{
		str = ReadLine(stream);	
		stringstream(str)>>temp>>temp>>temp>>temp>>is>>ie>>js>>je>>ks>>ke>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp;
					
		if(is != ie)
		{
			is = is;
			ie = ie + 1;
		}
		if(js != je)
		{
			js = js;
			je = je + 1;
		}
		if(ks != ke)
		{
			ks = ks;
			ke = ke + 1;
		}

		nn = 2;
		Num= 1;
		if((is-1)%nn ==0 && (ie-1)%nn ==0 && (js-1)%nn ==0 && (je-1)%nn ==0 && (ks-1)%nn ==0 && (ke-1)%nn ==0)
		{
			nn = nn*2;
			Num= Num+1;
		}
		NG = MIN(Num, NG);												 
	}

	// Mixing plane------------------------------------------------
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);				nMixingPlane = stoi(str);
	str = ReadLine(stream);
	for (i=1;i<=nMixingPlane;i++) 
	{
		str = ReadLine(stream);	
		stringstream(str)>>temp>>temp>>temp>>temp>>is>>ie>>js>>je>>ks>>ke>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp;
		if(is != ie)
		{
			is = is;
			ie = ie + 1;
		}
		if(js != je)
		{
			js = js;
			je = je + 1;
		}
		if(ks != ke)
		{
			ks = ks;
			ke = ke + 1;
		}

		nn = 2;
		Num= 1;
		if((is-1)%nn ==0 && (ie-1)%nn ==0 && (js-1)%nn ==0 && (je-1)%nn ==0 && (ks-1)%nn ==0 && (ke-1)%nn ==0)
		{
			nn = nn*2;
			Num= Num+1;
		}
		NG = MIN(Num, NG);												 
	}

	// Block------------------------------------------------
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);
	str = ReadLine(stream);				NB = stoi(str);
	str = ReadLine(stream);
	for (i=1;i<=NB;i++) 
	{
		str = ReadLine(stream);	
		stringstream(str)>>temp>>IM>>JM>>KM>>temp>>temp;

		nn = 2;
		for(j=1;j<=MIN3(IM, JM, KM);j++)
		{
			nn = 2*j;
			if((IM-1)%nn ==0 && (JM-1)%nn ==0 && (KM-1)%nn ==0 && (IM-1)/nn>2 && (JM-1)/nn>2 && (KM-1)/nn>2)
				Num = j + 1;
			else
				break;
	
			NG = MIN(Num, NG);
		}					     
	}
	stream.close();

	if(Solver.UseMultigrid==1)
		Solver.NG = MIN(NG, 2);		//本程序现在只是用双重网格
	else
		Solver.NG = 1;			// 用单重网格时用

	cout<<"use number of multigrids: "<<Solver.NG<<endl;

	Solver.Mesh = new CMesh [Solver.NG+1]; 

	for(i=1;i<=Solver.NG;i++)	
		Solver.Mesh[i].Mesh_ID = i;
}

void CUserInput::NonDimensionalization(CSolver &Solver)
{
	// reference quantity
	CFluidProps::Temp_ref = CBndConds::Ttin;
	CFluidProps::Press_ref = CBndConds::Ptin;
	CFluidProps::Vel_ref = SQRT(CFluidProps::Gamma*CFluidProps::Rgas*CFluidProps::Temp_ref);	//以滞止声速作为参考速度
	CFluidProps::mue_ref = 1.716e-5f*POW(CFluidProps::Temp_ref/288.15f, 3.0f/2.0f) * (288.15f + 110.4f)/(CFluidProps::Temp_ref + 110.4f);	
	CFluidProps::dens_ref = CFluidProps::Press_ref/CFluidProps::Rgas/CFluidProps::Temp_ref;
	CFluidProps::Re_ref = CFluidProps::dens_ref*CFluidProps::Vel_ref*CGeometry::Len_ref/CFluidProps::mue_ref;
	CFluidProps::Ma_ref = 1.0f;              //Mach数    （以声速作为参考速度，因而参考Mach数为1）
	// non-dimensional quantity
	CBndConds::Pout = CBndConds::Pout/(CFluidProps::dens_ref*CFluidProps::Vel_ref*CFluidProps::Vel_ref);
	CBndConds::Ttin = 1.0f;			//CBndConds::Ttin/CFluidProps::Temp_ref
	CBndConds::Ptin = 1.0f/CFluidProps::Gamma;		//Ptin = dens*R*T
	CFluidProps::Rgas = 1.0f/(CFluidProps::Ma_ref*CFluidProps::Ma_ref*CFluidProps::Gamma);
	int i, j;
	for(i=1;i<=Solver.NG;i++)
		for(j=1;j<=CMesh::NB;j++)
			Solver.Mesh[i].Block[j].FluidProp.w_rot = Solver.Mesh[i].Block[j].FluidProp.w_rot/(CFluidProps::Vel_ref/CGeometry::Len_ref);

		//以总温下的层流动力粘性作为参考粘性，初始化ev_tur时，用的便是总温下层运动流粘性的3-5倍，所以无量纲下，初始的ev便是3.0-5.0
		//无量纲下的总密度dens=1，总温T=1，总mue=1, ev_lam=总mue/总密度dens=1，故ev_tur=(3-5)*ev_lam=3-5
}

void CUserInput::ReadWallDistance(CMesh &Mesh)
{
	int i, j, k, nb, med_a, med_b, med_c, med_d;
	ifstream read;
	read.open("WallDistance.dat");
	for(nb=1;nb<=Mesh.NB;nb++)
	{
		read>>med_d;
		read>>med_a>>med_b>>med_c;
		cout<<"	Reading Block:"<<nb<<" wall distance..."<<endl;
		for(k=1;k<=Mesh.Block[nb].Geometry.KM-1;k++)
		{
			for(j=1;j<=Mesh.Block[nb].Geometry.JM-1;j++)
			{
				for(i=1;i<=Mesh.Block[nb].Geometry.IM-1;i++)
				{
					read>>med_a>>med_b>>med_c>>Mesh.Block[nb].Geometry.Cell[i][j][k].dis_Wall;
				}
			}
		}
	}
	read.close();
}

void CUserInput::ReadSolution(CMesh &Mesh)
{
	int i, j, k, nb;
	int med_i, med_j, med_k, med_nb;
	REAL press, uvel, vvel, wvel, temp, dens, ener, ev;
	ifstream read;
	read.open("RestartSolution.dat");
	for(nb=1;nb<=Mesh.NB;nb++)
	{
		cout<<"Reading Mesh[1].Block "<<nb<<"  Resolution"<<endl;
		for(k=1;k<=Mesh.Block[nb].Geometry.KM-1;k++)
		{
			for(j=1;j<=Mesh.Block[nb].Geometry.JM-1;j++)
			{
				for(i=1;i<=Mesh.Block[nb].Geometry.IM-1;i++)
				{
					if (CFluidProps::turbulenceType==TurbulenceTypes::SA)
						read>>med_nb>>med_i>>med_j>>med_k>>press>>uvel>>vvel>>wvel>>temp>>dens>>ener>>ev;
					else
						read>>med_nb>>med_i>>med_j>>med_k>>press>>uvel>>vvel>>wvel>>temp>>dens>>ener;

					Mesh.Block[nb].FluidProp.Cell_cv[i][j][k].dens = dens;
					Mesh.Block[nb].FluidProp.Cell_cv[i][j][k].xmom = dens * uvel;
					Mesh.Block[nb].FluidProp.Cell_cv[i][j][k].ymom = dens * vvel;
					Mesh.Block[nb].FluidProp.Cell_cv[i][j][k].zmom = dens * wvel;
					Mesh.Block[nb].FluidProp.Cell_cv[i][j][k].ener = ener;

					Mesh.Block[nb].FluidProp.Cell_pv[i][j][k].dens = dens;
					Mesh.Block[nb].FluidProp.Cell_pv[i][j][k].uvel = uvel;
					Mesh.Block[nb].FluidProp.Cell_pv[i][j][k].vvel = vvel;
					Mesh.Block[nb].FluidProp.Cell_pv[i][j][k].wvel = wvel;
					Mesh.Block[nb].FluidProp.Cell_pv[i][j][k].temp = temp;
					Mesh.Block[nb].FluidProp.Cell_pv[i][j][k].press = press;

					if (Mesh.Block[nb].FluidProp.turbulenceType==TurbulenceTypes::SA)
						Mesh.Block[nb].FluidProp.Cell_ev[i][j][k].ev = ev;
				}
			}
		}
	}
	read.close();
}

