//----------------------------------------------------------------------------------||
//-------------------                  SIRC.cpp	                 -------------------||
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


int GlobalVars::Rank=0;
FILE* GlobalVars::LogFile=NULL;

int main(int argc, char** argv)
{

	int N_processor;
	int Rank;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &N_processor);
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

	GlobalVars::Rank=Rank;

	#ifdef _DEBUG
		char sFile[128];
		sprintf(sFile, "log_r_%d.dat", Rank);
		GlobalVars::LogFile = fopen(sFile, "w");
	#endif


	Log("=============================================");
	Log(" 	   ____ ___ ____   ____ 	");
	Log(" 	  / ___|_ _|  _ \\ / ___| 	");
	Log(" 	  \\___ \\| || |_) | |      ");
	Log(" 	   ___) | ||  _ <| |___ 	");
	Log(" 	  |____/___|_| \\_\\\\____|\n");
	Log("=============================================");
	Log("==== Starting Program: SIRC.");
	Log("==== Copyright: (C) 2022 by Tianhong Wang.");
	Log("==== Number of Total Ranks: [%d].",N_processor);
	Log("=============================================\n");
	DLog("Rank_%03d: Log file created.",Rank);



	Domain *domain = new Domain((char*)"SIRC.ini",Rank);
	domain-> Run();
	
	delete domain;

	MPI_Finalize();

	return 0;
}