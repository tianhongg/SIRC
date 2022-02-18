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



	AddEntry((char*)"TStep", &dt,	1.0);
	AddEntry((char*)"TMax",  &Tmax,	1.0);
	AddEntry((char*)"Wavelength", 		&lambda_L,	1.0);
	AddEntry((char*)"ReadType",   		&ReadType,	10);
	AddEntry((char*)"OutputInterval",  	&Out_dt,	1);


	Log("Domain: Read Parameters From ini File...");
	FILE *p_File = fopen(infile,"rt");
	if (p_File)
	{
		rewind(p_File);
		read(p_File);
	}
	fclose(p_File); 

	//normalize
	dt     *= 2*Constant::PI;
	Tmax   *= 2*Constant::PI;
	Out_dt *= 2*Constant::PI;

	//for the tick
	time = 0.0;
	d_tick = Tmax/200.01;
	n_out=n_tick=0;

	Log("Domain: Create Detector...");

	MyDetector = new Detector((char*)"SIRC.ini");

}


void Domain::Run()
{

	Log("Domain::Run: Read Trajectory...");
	if(this->ReadTrajectory())
	{
		Log("Domain::Run: Read Trajectory Failed");
		return;
	};

	MPI_Barrier(MPI_COMM_WORLD);
	Log("Domain::Run: Start Calculation");
	this->OnCalculate();

	Log("Domain::Run: Done!");

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