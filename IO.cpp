//----------------------------------------------------------------------------------||
//-------------------                  IO.cpp	                 -------------------||
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


#include <hdf5.h>
#include "SIRC.h"

// general IO

void Log0(bool all, const char * fmt,...)
{

	if(!all && GlobalVars::Rank>0) return;

	const char *fmtn;
	string str(fmt);
	str+="\n";
	fmtn=str.c_str();
	va_list arglist;
    va_start( arglist, fmt );
    vprintf( fmtn, arglist);
    va_end(arglist);
}



void LogDebug(const char * fmt,...)
{

	#ifndef _DEBUG
		return;
	#endif

	const char *fmtn;
	string str(fmt);
	str+="\n";
	fmtn=str.c_str();
	va_list arglist;
    va_start( arglist, fmt );
    vfprintf(GlobalVars::LogFile,fmtn, arglist);
    va_end(arglist);

}



int Domain::ReadTrajectory()
{
	string str = FILENAME;
	char filename[128];

	if(ReadType==0)
	{
		sprintf(filename,"_%d_.h5", 0);
	}
	else
	{
		sprintf(filename,"_%d_.h5", Rank);
	}

	str+=string(filename);

	Log("==== Domain::ReadTrajectory: Read Trajectory Type [%d]", ReadType);
	DLog("Domain::ReadTrajectory: Reading Trajectories From File: [%s]",str.c_str());

	// open file 
	hid_t       file_id;  
	hsize_t    	size_p;
	herr_t      status; 

	file_id = H5Fopen (str.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	status  = H5Gget_num_objs(file_id, &size_p);

	// number of particles to read
	int N_Processor;
	MPI_Comm_size(MPI_COMM_WORLD, &N_Processor);

	ULONG N_read  = size_p; // number of particles
	ULONG Idx_read=0; // index to start reading

	if(ReadType==0)
	{
		N_read = size_p/N_Processor;
		if(N_read<1)
		{
			Log("==== Domain::ReadTrajectory: Too Many Cores, Not Enough Particles");
			return 1;
		}
		Idx_read = Rank*N_read;
		if(N_Processor-Rank<size_p%N_Processor+1)  
		{
			N_read++; 
			Idx_read += size_p%N_Processor-(N_Processor-Rank);
		}                        	
	}

	ALog("Domain::ReadTrajectory: Reading [%d] Trajectories [%d -> %d]", N_read, Idx_read, Idx_read+N_read-1);
	DLog("Domain::ReadTrajectory: Reading [%d] Trajectories [%d -> %d]", N_read, Idx_read, Idx_read+N_read-1);

	for(ULONG i=0; i<N_read; i++)
	{
		Particle* OneParticle = new Electron();
		Particles.push_back(OneParticle);

		if(OneParticle->Load(file_id, Idx_read, i))
			DLog("Domain::ReadTrajectory: Reading Trajectory [%d] Failed at Rank [%d]",Idx_read+i,Rank);
	}


	//close
	status = H5Fclose(file_id);


	return 0;

}



int Particle::Load(hid_t file_id,ULONG Idx_read, ULONG i)
{

	return 0;
}