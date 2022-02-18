//----------------------------------------------------------------------------------||
//-------------------               radiation.cpp                -------------------||
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
// #include <unistd.h>


void Domain::OnCalculate()
{
	Log("Domain::OnCalculate--------------------");

	tic = chrono::system_clock::now();

	//time loop
	for(step = 0; step<MaxStep; step++)
	{
		this->Tick();

		//particle loop
		for(auto it = Particles.begin(); it!=Particles.end(); it++)
		{
			Particle* p = *it;

			//trajectory may have different starting time and length
			// so checck if the trajectory is in frame;
			if(p->IsInFrame(step)==false) continue;
			MyDetector->OnDeposit(p);
		}

	}

}


void Detector::OnDeposit(Particle* p)
{
	for(auto it = Pixels.begin(); it!=Pixels.end(); it++)
		(*it)->OnDeposit(p);
}

void Pixel::OnDeposit(Particle* p)
{	
	//calculate:   -nxnx\beta;
	Vec3 vm = p->Velocity[p->Current_Step-1]; 
	Vec3 vc = p->Velocity[p->Current_Step];
	Vec3 vp = p->Velocity[p->Current_Step+1];
	
}


//=====================================================================
void Domain::Tick()
{
	//tick
	// usleep(50000);
	if(step>=(n_tick)*d_tick)
	{
		n_tick++;
		chrono::time_point<std::chrono::system_clock> toc = chrono::system_clock::now();
		chrono::milliseconds diff = std::chrono::duration_cast<std::chrono::milliseconds>(toc - tic);
   		double delay = diff.count() / 100.0;


   		string fmt = "Domain::OnCalculate: [%5.1f%%] Time [%5.2f min "; 

   		int n = n_tick/10;
   		for(int i=0; i<n;	i++) fmt+="*";
   		for(int i=0; i<20-n;i++) fmt+="-";
   		fmt+=" %5.2f min]"; 

   		 Log(fmt.c_str(), n_tick/2.0, delay/60, 1.0*(200-n_tick)/n_tick*delay/60);
   		DLog(fmt.c_str(), n_tick/2.0, delay/60, 1.0*(200-n_tick)/n_tick*delay/60); 
	}

	//output;
	if((step+1)%Out_dt==0)
	{
		Log("Domain::OnCalculate: Output [%d].",++n_out);
		this->Output(n_out);
	}
}
