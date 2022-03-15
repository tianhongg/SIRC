//----------------------------------------------------------------------------------||
//-------------------                  domain.cpp                -------------------||
//----------------------------------------------------------------------------------||
//                               ____ ___ ____   ____                               ||
//                              / ___|_ _|  _ \ / ___|                              ||
//                              \___ \| || |_) | |                                  ||
//                               ___) | ||  _ <| |___                               ||
//                              |____/___|_| \_\\____|                              ||
//                                                                                  ||
//----------------------------------------------------------------------------------||
//--             (S)imple (I)ncoherent (R)adiation (C)alculation                   -||
//----------------------------------------------------------------------------------||
//---Author-----------           : Tianhong Wang                --------------------||
//---Starting---------           : Feb-16-2022                  --------------------||
//---Email------------           : tw474@cornell.edu            --------------------||
//---Group------------           : Dr. Gennady Shvets' Group    --------------------||
//---Copyright--------           : (C) 2022 by Tianhong Wang    --------------------||
//----------------------------------------------------------------------------------||
//----------------------------------------------------------------------------------||

#include "SIRC.h"

//Domain::p_D = NULL;
Domain* Domain::p_Domain= NULL;



Domain::Domain (char * infile, int rank) : NList("Domain") 
{

	Domain::p_Domain = this; 
	Rank = rank;
	Size_P=0;
	BunchXiMax = -1e10;
	BunchXiMin = +1e10;

	AddEntry((char*)"TStep", 			&dt,			1.0);
	AddEntry((char*)"MaxSteps",  		&MaxStep,		10);
	AddEntry((char*)"Refine",			&Refine, 		1);
	AddEntry((char*)"Wavelength", 		&lambda_L,		1.0);
	AddEntry((char*)"ReadType",   		&ReadType,		1);
	AddEntry((char*)"MovingFrame",  	&MovingFrame,	1);
	AddEntry((char*)"InputType",  		&InputType,		1);
	AddEntry((char*)"Normalization",  	&Normalization,	1);
	AddEntry((char*)"IntegrateOrder",  	&IntegrateOrder, 2);
	AddEntry((char*)"IncludePartI",  	&IncludePartI,   1);
	
	
	Log("Domain: Read Parameters From .ini File...");
	FILE *p_File = fopen(infile,"rt");
	if (p_File)
	{
		rewind(p_File);
		read(p_File);
	}
	fclose(p_File); 


	Refine = max(1,Refine);
	if(Refine>1)
	{
		Refine=ceil(Refine/2.0)*2;
		Log("Domain: Stepsize Refinement: [%d]; Cubic Interp Will be Used...",Refine);
	}
	else
	{
		Log("Domain: No Stepsize Refinement...");
	}

	IntegrateOrder=max(1,min(IntegrateOrder,2));


	IntegrateOrder = 1;
	Log("Domain: Integrate Order = [1] is Fixed in Current Version...");

	//normalize
	if(Normalization)
		dt *= 2*Constant::PI;

	//for the tick
	step = 0;
	n_tick=0;
	s_tick=0;

}


void Domain::Run()
{

	Log("Domain::Run--------------------");
	Log("Domain::Run: Read Trajectory...");
	if(this->ReadTrajectory())
	{
		Log("Domain::Run: Read Trajectory Failed");
		return;
	};

	this->ReduceBunchSize();

	Log("Domain::Run: Create Detector...");
	MyDetector = new Detector((char*)"SIRC.ini");

	Log("Domain::Run: TimeTagging Particles...");
	this->TagParticles();

	MPI_Barrier(MPI_COMM_WORLD);
	Log("Domain::Run: Start Calculation");
	this->OnCalculate();

	Log("Domain::Run: On Output");
	this->Output();

	Log("Domain::Run: Done!");
   	Log("Domain::Run: Total Runtime[%6.2f min]", RunTimeInSec()/60);

}


void Domain::ReduceBunchSize()
{
	double bmax;
	double bmin;
	MPI_Allreduce(&BunchXiMin, &bmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&BunchXiMax, &bmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	BunchXiMin = bmin;
	BunchXiMax = bmax;
	
	 Log("Domain::ReduceBunchSize(): Bunch Size (k^-1) in Moving Frame: [%.2f, %.2f].",bmin,bmax);
}


void Domain::TagParticles()
{

	if(N_Time==1) return;

	for(auto it = Particles.begin(); it!=Particles.end(); it++)
	{
		Particle * p = *it;		
		double dxi = (BunchXiMax-BunchXiMin)/(N_Time-1);
		p->t_bin = (ULONG)floor((p->xi - BunchXiMin)/dxi);
	}

}

//---------------------------- Domain::~Domain() -----------------------
Domain::~Domain()
{

	delete MyDetector;
	for(auto it : Particles)
	{
		delete it;
	}
	Particles.clear();
	if(GlobalVars::LogFile)
	{

		fclose(GlobalVars::LogFile);
	}
};