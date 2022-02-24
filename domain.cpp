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
	AddEntry((char*)"StepRefine",  		&Refine,		1);
	AddEntry((char*)"Wavelength", 		&lambda_L,		1.0);
	AddEntry((char*)"ReadType",   		&ReadType,		1);
	AddEntry((char*)"OutputInterval",  	&Out_dt,		1);
	AddEntry((char*)"MovingFrame",  	&MovingFrame,	1);
	AddEntry((char*)"InputType",  		&InputType,		1);
	AddEntry((char*)"Normalization",  	&Normalization,	1);
	AddEntry((char*)"IntegrateOrder",  	&IntegrateOrder, 2);
	
	Log("Domain: Read Parameters From .ini File...");
	FILE *p_File = fopen(infile,"rt");
	if (p_File)
	{
		rewind(p_File);
		read(p_File);
	}
	fclose(p_File); 

	Refine = max(1,Refine);

	if(IntegrateOrder<1) IntegrateOrder=1;
	if(IntegrateOrder>2) IntegrateOrder=2;

	//normalize
	if(Normalization)
		dt *= 2*Constant::PI;

	//for the tick
	step = 0;
	d_tick = MaxStep/200.01;
	n_out=n_tick=0;

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

	Log("Domain::Run: Done!");

}



void Domain::ReduceBunchSize()
{
	double bmax;
	double bmin;
	MPI_Allreduce(&BunchXiMin, &bmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&BunchXiMax, &bmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	BunchXiMin = bmin;
	BunchXiMax = bmax;
	
	 Log("Domain::ReduceBunchSize(): Bunch Size (kp) in Moving Frame: [%.2f, %.2f].",bmin,bmax);
	DLog("Domain::ReduceBunchSize(): Bunch Size (kp) in Moving Frame: [%.2f, %.2f].",bmin,bmax);

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