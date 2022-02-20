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
	double t  = p_domain()->GetTime();
	double dt = p_domain()->GetDt();
	dcom k = Constant::i*omega*dt*0.5;// i*omega*dt/2

	//[1] calculate:   f = -nxnx\beta;
	Vec3 vm = p->Velocity[p->Current_Step-1]; 
	Vec3 vc = p->Velocity[p->Current_Step];
	Vec3 vp = p->Velocity[p->Current_Step+1];

	Vec3 fm = -(n.Cross(n.Cross(vm)));
	Vec3 fc = -(n.Cross(n.Cross(vc)));
	Vec3 fp = -(n.Cross(n.Cross(vp)));


	//[2] calculate:   g =  t-n*r;
	Vec3 rm = p->Position[p->Current_Step-1]; 
	Vec3 rc = p->Position[p->Current_Step];
	Vec3 rp = p->Position[p->Current_Step+1];

	double gm = - (n.Dot(rm)) + t;
	double gc = - (n.Dot(rc)) + t;
	double gp = - (n.Dot(rp)) + t;

	//[3] calculate f1 f2 f3; g1, g2, g3
	Vec3 f1 = fc;
	Vec3 f2 = (fp-fm)/4.0;
	// Vec3 f3 = (fp-fc*2.0+fm)/4.0;

	double g1 = gc;
	double g2 = (gp-gm)/4.0;
	// double g3 = (gp-gc*2.0+gm)/4.0;


	//2nd-order?
	if(abs(g2)==0)
	{
		dcom tmp = (k*2.0*Constant::i)*exp(Constant::i*g1*omega)*omega;
		this->Ax[p->t_bin] += tmp*f1.x;
		this->Ay[p->t_bin] += tmp*f1.y;
		this->Az[p->t_bin] += tmp*f1.z;
	}
	else
	{	
		dcom tmp =(k*2.0)*exp(Constant::i*g1*omega)/(g2*g2*omega);
		this->Ax[p->t_bin] += tmp*(f2.x*g2*omega*cos(g2*omega)-(f2.x-Constant::i*f1.x*g2*omega)*sin(g2*omega));
		this->Ay[p->t_bin] += tmp*(f2.y*g2*omega*cos(g2*omega)-(f2.y-Constant::i*f1.y*g2*omega)*sin(g2*omega));
		this->Az[p->t_bin] += tmp*(f2.z*g2*omega*cos(g2*omega)-(f2.z-Constant::i*f1.z*g2*omega)*sin(g2*omega));
	}

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
   		double delay = diff.count() / 1000.0;

   		string fmt = "Domain::OnCalculate: [%5.1f%%] Time [%6.2f min "; 

   		int n = n_tick/10;
   		for(int i=0; i<n;	i++) fmt+="*";
   		for(int i=0; i<20-n;i++) fmt+="-";
   		fmt+=" %6.2f min]"; 

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
