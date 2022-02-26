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
    fflush(GlobalVars::LogFile);
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

	ALog("Domain::ReadTrajectory: Estimated Size of Particle Data: [%d Kilobytes]", (sizeof(Particle)+sizeof(Vec3)*2*MaxStep)*N_read/1024);
	DLog("Domain::ReadTrajectory: Estimated Size of Particle Data: [%d Kilobytes]", (sizeof(Particle)+sizeof(Vec3)*2*MaxStep)*N_read/1024);

	for(ULONG i=0; i<N_read; i++)
	{
		Particle* OneParticle = new Electron();
		Particles.push_back(OneParticle);

		if(OneParticle->Load(file_id, Idx_read+i))
			DLog("Domain::ReadTrajectory: Reading Trajectory [%d] Failed ",Idx_read+i);
	}


	//close
	status = H5Fclose(file_id);

	ALog("Domain::ReadTrajectory: Size of Particle Data: [%d Kilobytes]", Size_P/1024);
	DLog("Domain::ReadTrajectory: Size of Particle Data: [%d Kilobytes]", Size_P/1024);

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
	p_domain()->Size_P += (sizeof(Particle)+sizeof(Vec3)*2*nt);


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

	// ALog("Particle starting step, %d",(this->start));

	this->Normalize();


	return 0;
}


int WriteHDF5RecordDOUBLE(hid_t file, const char* dataname, int nfields, double *source);

void Domain::Output()
{

	Log("Domain::Output-------------------------------------");

	// d^2I/domega/dOmega = e^2/(4*pi^2*c); 
	// here we output  d^2I/domega/dOmega when omega is normalized by laser wavelength 
	// so the prefactor is (mc^2) re/lambda_L/2/pi
	// output in the unit of mc^2 (let's convert it to eV)
	// and convert mc^2 to eV.

	// radius of electron, divided by laser wavelength, then convert to eV/
	double re = 2.8179403227e-9/(this->lambda_L)/2/Constant::PI*0.51099895e6;

	//output will be  re*|A|^2;

	hid_t   file, fdataset, fdataspace; 
	hsize_t  dimsfi[0], dimsf[4]; 
	hid_t rank4 = 4;
	herr_t     status;
	char fname[128];
	//create file
	if(Rank==0)
	{
		sprintf(fname,"Synchrotron_SIRC.h5");
		Log("Domain::Output: Create File: %s",fname);
		file = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

		//
		DLog("Domain::Output: Write E_Bin, Unit eV...");
		double* Axis = new double[MyDetector->N_Omega];

		for(ULONG i=0; i<MyDetector->N_Omega; i++)
			Axis[i] = (MyDetector->OmegaBin[i])*1.2398/lambda_L;

		dimsfi[0] = MyDetector->N_Omega;
		WriteHDF5RecordDOUBLE(file, "Energy[eV]", dimsfi[0], Axis);
		delete[] Axis;
		//
		DLog("Domain::Output: Write Theta_Y_Bin");
		dimsfi[0] = MyDetector->N_Theta_Y;
		WriteHDF5RecordDOUBLE(file, "ThetaY", dimsfi[0], MyDetector->ThetaYBin);
		//
		DLog("Domain::Output: Write Theta_Z_Bin");
		dimsfi[0] = MyDetector->N_Theta_Z;
		WriteHDF5RecordDOUBLE(file, "ThetaZ", dimsfi[0], MyDetector->ThetaZBin);

		//
		if(N_Time>1)
		{
			DLog("Domain::Output: Write Time_Bin Unit fs...");
			Axis = new double[MyDetector->N_Time];

			double dxi = (BunchXiMax-BunchXiMin)/(N_Time-1);
			for(ULONG i=0;  i<MyDetector->N_Time; i++)
				Axis[i] = (-dxi*i)*(lambda_L*10/3.0);

			dimsfi[0] = MyDetector->N_Time;
			WriteHDF5RecordDOUBLE(file, "Time[fs]", dimsfi[0], Axis);
			delete[] Axis;
		}


		DLog("Domain::Output: Create Dataspace");
		dimsf[0] = MyDetector->N_Omega;
		dimsf[1] = MyDetector->N_Theta_Y;
		dimsf[2] = MyDetector->N_Theta_Z;
		dimsf[3] = MyDetector->N_Time;
		fdataspace = H5Screate_simple(rank4, dimsf, NULL); 
		
	}
	//gather data;
	 Log("Domain::Output: Gather Data...");
	DLog("Domain::Output: Gather Data...");

	//declare data
	double* local_A; 
	double* global_A; 
	int N_Buffer = (MyDetector->N_Omega)*(MyDetector->N_Theta_Y)*(MyDetector->N_Theta_Z)*(MyDetector->N_Time);

	DLog("Domain::Output: Data Buffer Size Needed [%d Kilobytes]:", N_Buffer*sizeof(double)*2/1024);

	local_A  = new double[N_Buffer];
	global_A = new double[N_Buffer];


	// load data S1
	int k = 0;
	for(auto it = MyDetector->Pixels.begin(); it!=MyDetector->Pixels.end(); it++)
		for(int i=0; i<N_Time; i++)
			local_A[k++] = (*it)->S1[i]*re;
	MPI_Reduce(local_A, global_A, N_Buffer, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if(Rank==0)
	{
		 Log("Domain::Output: Write S1 = |Ay|^2+|Az|^2");
		DLog("Domain::Output: Write S1 = |Ay|^2+|Az|^2");
    	fdataset = H5Dcreate(file, "S1", H5T_IEEE_F64LE, fdataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    	status = H5Dwrite(fdataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_A);
    	status = H5Dclose(fdataset);
	}
	// load data S2
	k = 0;
	for(auto it = MyDetector->Pixels.begin(); it!=MyDetector->Pixels.end(); it++)
		for(int i=0; i<N_Time; i++)
			local_A[k++] = (*it)->S2[i]*re;
	MPI_Reduce(local_A, global_A, N_Buffer, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if(Rank==0)
	{
		 Log("Domain::Output: Write S2 = |Ay|^2-|Az|^2");
		DLog("Domain::Output: Write S2 = |Ay|^2-|Az|^2");
    	fdataset = H5Dcreate(file, "S2", H5T_IEEE_F64LE, fdataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    	status = H5Dwrite(fdataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_A);
    	status = H5Dclose(fdataset);
    }

    // load data S3
	k = 0;
	for(auto it = MyDetector->Pixels.begin(); it!=MyDetector->Pixels.end(); it++)
		for(int i=0; i<N_Time; i++)
			local_A[k++] = (*it)->S3[i]*re;
	MPI_Reduce(local_A, global_A, N_Buffer, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if(Rank==0)
	{	 
		 Log("Domain::Output: Write S3 = 2Re(AyAz*)");
		DLog("Domain::Output: Write S3 = 2Re(AyAz*)");
    	fdataset = H5Dcreate(file, "S3", H5T_IEEE_F64LE, fdataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    	status = H5Dwrite(fdataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_A);
    	status = H5Dclose(fdataset);
	}
	// load data S4
	k = 0;
	for(auto it = MyDetector->Pixels.begin(); it!=MyDetector->Pixels.end(); it++)
		for(int i=0; i<N_Time; i++)
			local_A[k++] = (*it)->S4[i]*re;
	MPI_Reduce(local_A, global_A, N_Buffer, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if(Rank==0)
	{
		 Log("Domain::Output: Write S4 =-2Im(AyAz*)");
		DLog("Domain::Output: Write S4 =-2Im(AyAz*)");
    	fdataset = H5Dcreate(file, "S4", H5T_IEEE_F64LE, fdataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    	status = H5Dwrite(fdataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_A);
    	status = H5Dclose(fdataset);
    }

    // load data I
	k = 0;
	for(auto it = MyDetector->Pixels.begin(); it!=MyDetector->Pixels.end(); it++)
		for(int i=0; i<N_Time; i++)
			local_A[k++] = (*it)->II[i]*re;
	MPI_Reduce(local_A, global_A, N_Buffer, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if(Rank==0)
	{
		 Log("Domain::Output: Write II = |Ax|^2+|Ay|^2+|Az|^2");
		DLog("Domain::Output: Write II = |Ax|^2+|Ay|^2+|Az|^2");
    	fdataset = H5Dcreate(file, "II", H5T_IEEE_F64LE, fdataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    	status = H5Dwrite(fdataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, global_A);
    	status = H5Dclose(fdataset);
	}
	

    //close
	if(Rank==0)
	{
		Log("Domain::Output: Close File: %s",fname);
    	status = H5Fclose(file);
    	status = H5Sclose(fdataspace);
	}    
    

	delete[] local_A;
	delete[] global_A; 

	Log("Domain::Output-------------------------------------");

}


int WriteHDF5RecordDOUBLE(hid_t file, const char* dataname, int nfields, double *source)
{
   hid_t rank1 = 1;
   hsize_t dimsfi[3];
   hid_t       dataset;         /* File and dataset            */
   hid_t       dataspace;   /* Dataspace handles           */
   herr_t      status;                /* Error checking              */

   assert (file >= 0);

   dimsfi[0] = nfields;
   dataspace = H5Screate_simple(rank1, dimsfi, NULL); 
   assert (dataspace >= 0);
   dataset = H5Dcreate(file, dataname, H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   assert (dataset >= 0);
   status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, source);
   assert (status >= 0);
   status = H5Dclose(dataset);
   assert (status >= 0);
   status = H5Sclose(dataspace);
   assert (status >= 0);
   int istatus = status;
   return istatus;
}