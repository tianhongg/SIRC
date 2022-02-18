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

string return_current_time_and_date()
{
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);

    std::stringstream ss;
    ss<<"[";
    ss<<setfill('0')<<setw(4)<<GlobalVars::Rank;
    ss << std::put_time(std::localtime(&in_time_t), "][%H:%M:%S] ");
    return ss.str();
}
void Log0(bool all, const char * fmt,...)
{

	if(!all && GlobalVars::Rank>0) return;

	const char *fmtn;
	string str(fmt);
	str=return_current_time_and_date()+str+"\n";

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
	Log("Domain::ReadTrajectory--------------------");
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

	Log("Domain::ReadTrajectory: Read Trajectory Type [%d]", ReadType);
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
			Log("Domain::ReadTrajectory: Too Many Cores, Not Enough Particles");
			return 1;
		}
		Idx_read = Rank*N_read;
		if(N_Processor-Rank<size_p%N_Processor+1)  
		{
			N_read++; 
			Idx_read += size_p%N_Processor-(N_Processor-Rank);
		}                        	
	}

	ALog("Domain::ReadTrajectory: Reading [%d] Trajectories; Range [%d -> %d]", N_read, Idx_read, Idx_read+N_read-1);
	DLog("Domain::ReadTrajectory: Reading [%d] Trajectories; Range [%d -> %d]", N_read, Idx_read, Idx_read+N_read-1);

	for(ULONG i=0; i<N_read; i++)
	{
		Particle* OneParticle = new Electron();
		Particles.push_back(OneParticle);

		if(OneParticle->Load(file_id, Idx_read+i))
			DLog("Domain::ReadTrajectory: Reading Trajectory [%d] Failed ",Idx_read+i);
	}


	//close
	status = H5Fclose(file_id);

	ALog("Domain::ReadTrajectory: Size of Particle Data: [%d Kilobytes]", N_read*Size_P/1024);
	DLog("Domain::ReadTrajectory: Size of Particle Data: [%d Kilobytes]", N_read*Size_P/1024);

	

	return 0;

}



int Particle::Load(hid_t file_id, ULONG i)
{

	DLog("Particle::Load: Reading Trajectory [%d]",i);

	// h5
	hid_t	group_id, dataset_id; 
	hsize_t	num_obj;
	H5G_info_t group_info;

	// get name
	char *particle_name;
	ssize_t len = H5Gget_objname_by_idx(file_id, i, NULL, 0 );
	particle_name=new char[++len]; 
	H5Gget_objname_by_idx(file_id, i, particle_name, len);
	if(H5Gget_info_by_name( file_id, particle_name, &group_info, H5P_DEFAULT )<0)
	{
		ALog("Failed at Particle::Load::H5Gget_info_by_name()");
		DLog("Failed at Particle::Load::H5Gget_info_by_name()");
		return 1;
	}

	//get group
	group_id=H5Gopen( file_id, particle_name, H5P_DEFAULT);

	if(H5Gget_num_objs(group_id,&num_obj)<0)
	{
		ALog("Failed at Particle::Load::H5Gget_num_objs()");
		DLog("Failed at Particle::Load::H5Gget_num_objs()");
		return 1;
	}


	// read the position and get the size first
	dataset_id=H5Dopen1(group_id,"xx");
	hid_t dspace = H5Dget_space(dataset_id);
	const int ndims = H5Sget_simple_extent_ndims(dspace);
	hsize_t dims[ndims];
	H5Sget_simple_extent_dims(dspace, dims, NULL);
	ULONG nt = dims[0];
	//count the memory size
	p_domain()->Size_P = max(p_domain()->Size_P,sizeof(Particle)+sizeof(Vec3)*2*nt);


	DLog("Particle::Load: Particle [%d] Total Steps [%d]",i, nt);

	// make buff and allocate trajectory
	double databuff[nt];
	this->NStep = nt;
	Position = new Vec3[nt];
	Velocity = new Vec3[nt];

	//read x
	if(H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,databuff)<0)
	{
		DLog("Particle::Load: Failed at Particle [%i] xx",i);
		return 1;
	}
	for(int t=0; t<nt; t++)
		Position[t].x = databuff[t];

	//read y
	dataset_id=H5Dopen1(group_id,"yy");
	if(H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,databuff)<0)
	{
		DLog("Particle::Load: Failed at Particle [%i] yy",i);
		return 1;
	}
	for(int t=0; t<nt; t++)
		Position[t].y = databuff[t];

	//read z
	dataset_id=H5Dopen1(group_id,"zz");
	if(H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,databuff)<0)
	{
		DLog("Particle::Load: Failed at Particle [%i] zz",i);
		return 1;
	}
	for(int t=0; t<nt; t++)
		Position[t].z = databuff[t];


	//read px
	dataset_id=H5Dopen1(group_id,"px");
	if(H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,databuff)<0)
	{
		DLog("Particle::Load: Failed at Particle [%i] px",i);
		return 1;
	}
	for(int t=0; t<nt; t++)
		Velocity[t].x = databuff[t];


	//read py
	dataset_id=H5Dopen1(group_id,"py");
	if(H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,databuff)<0)
	{
		DLog("Particle::Load: Failed at Particle [%i] py",i);
		return 1;
	}
	for(int t=0; t<nt; t++)
		Velocity[t].y = databuff[t];


	//read pz
	dataset_id=H5Dopen1(group_id,"pz");
	if(H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,databuff)<0)
	{
		DLog("Particle::Load: Failed at Particle [%i] pz",i);
		return 1;
	}
	for(int t=0; t<nt; t++)
		Velocity[t].z = databuff[t];

	//weight
	dataset_id=H5Dopen1(group_id,"weight");
	if(H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,&(this->weight))<0)
	{
		DLog("Particle::Load: Failed at Particle [%i] weight",i);
		return 1;
	}


	//start
	dataset_id=H5Dopen1(group_id,"start");
	if(H5Dread(dataset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,&(this->start))<0)
	{
		DLog("Particle::Load: Failed at Particle [%i] start",i);
		return 1;
	}

	this->Normalize();

	// if(i==0)
	// {
	// 	printf("\n[%d,%d]\n",(int)this->weight,(int)this->start);
			
	// 	for(int t=0;t<2;t++)
	// 	{
	// 		printf("[%4.2f,%4.2f,%4.2f]\n",Velocity[t].x,Velocity[t].y,Position[t].x);
	// 	}
	// }

	return 0;
}


void Domain::Output(int n)
{



	
}