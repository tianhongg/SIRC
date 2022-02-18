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

#include <unistd.h>


void Domain::OnCalculate()
{
	tic = chrono::system_clock::now();

	//time loop
	for(time = 0; time<=Tmax+dt/2; time+=dt)
	{
		this->Tick();
		//particle loop
		for(auto p = Particles.begin(); p!=Particles.end(); p++)
		{

				// check if particle is in frame

				// get f1 f2 f3, g1, g2 ,g3

				// OnDeposit()


		}

	}

}



void Domain::Tick()
{
	//tick
	if(time>=(n_tick+1)*d_tick)
	{
		n_tick++;
		chrono::time_point<std::chrono::system_clock> toc = chrono::system_clock::now();
		chrono::milliseconds diff = std::chrono::duration_cast<std::chrono::milliseconds>(toc - tic);
   		double delay = diff.count() / 100.0;


   		string fmt = "Domain: [%5.1f%%] Time [%5.2f min "; 

   		int n = n_tick/10;
   		for(int i=0; i<n;	i++) fmt+=">";
   		for(int i=0; i<20-n;i++) fmt+="-";
   		fmt+=" %5.2f min]"; 

   		 Log(fmt.c_str(), n_tick/2.0, delay/60, 1.0*(200-n_tick)/n_tick*delay/60,(char)254u);
   		DLog(fmt.c_str(), n_tick/2.0, delay/60, 1.0*(200-n_tick)/n_tick*delay/60,(char)254u); 
	}

	//output;
	if(time>=(n_out+1)*Out_dt)
	{
		Log("Domain::OnCalculate: Output [%d].",++n_out);
		this->Output(n_out);
	}
}
