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



	AddEntry((char*)"TStep", &dt,	1.0);
	AddEntry((char*)"Tmax",  &Tmax,	1.0);
	AddEntry((char*)"Wavelength", &lambda_L,	1.0);
	AddEntry((char*)"ReadType",   &ReadType,	1);


	Log("==== Domain: Read Parameters From ini File...");
	FILE *p_File = fopen(infile,"rt");
	if (p_File)
	{
		rewind(p_File);
		read(p_File);
	}
	fclose(p_File); 


	Log("==== Domain: Create Detector...");

	MyDetector = new Detector((char*)"SIRC.ini");

	
}


void Domain::Run()
{

	Log("==== Domain::Run: Read Trajectory...");
	if(this->ReadTrajectory())
	{
		Log("==== Domain::Run: Read Trajectory Failed");
		return;
	};

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