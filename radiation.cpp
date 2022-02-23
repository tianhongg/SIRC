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
	int R = p_domain()->GetRefine();
	for(auto it = Pixels.begin(); it!=Pixels.end(); it++)
		for (int s=0; s<R; s++)
			(*it)->OnDeposit(p,s,R);
}


void Pixel::OnDeposit(Particle* p, int substep, int Refine)
{	
	double dt = p_domain()->GetDt();
	double shift = (2*substep+1-Refine)/2.0/Refine; //shift in the unit of dt;
	double t  = p_domain()->GetTime(); // time interpolation

	double R = (double)Refine;
	dcom k = dt*0.5*(p->weight); // weight*dt/2

	double hm = 1 + shift;
	double hp = 1 - shift;

	//[1] calculate:   f = -nxnx\beta;
	Vec3 vm = p->Velocity[p->Current_Step-1]; 
	Vec3 vc = p->Velocity[p->Current_Step  ];
	Vec3 vp = p->Velocity[p->Current_Step+1];

	Vec3 fm = -(n.Cross(n.Cross(vm)));
	Vec3 fc = -(n.Cross(n.Cross(vc)));
	Vec3 fp = -(n.Cross(n.Cross(vp)));

	//[2] calculate:   g =  t-n*r;
	Vec3 rm = p->Position[p->Current_Step-1]; 
	Vec3 rc = p->Position[p->Current_Step  ];
	Vec3 rp = p->Position[p->Current_Step+1];

	double gm = - (n.Dot(rm)) + t-dt;
	double gc = - (n.Dot(rc)) + t;
	double gp = - (n.Dot(rp)) + t+dt;

	//[3] calculate f1 f2 f3; g1, g2, g3
	Vec3   f1 = fc + (fp-fm)/2.0*shift + (fp-fc*2.0+fm)/2.0*shift*shift;  // with interpolation
	double g1 = gc + (gp-gm)/2.0*shift + (gp-gc*2.0+gm)/2.0*shift*shift;  // with interpolation

	Vec3   f2 = (fp-fc)/4.0*hm/hp + (fc-fm)/4.0*hp/hm; // with interpolation
	double g2 = (gp-gc)/4.0*hm/hp + (gc-gm)/4.0*hp/hm; // with interpolation
	
	// Vec3 f3 = (fp-fc*2.0+fm)/8.0;
	// double g3 = (gp-gc*2.0+gm)/8.0;

	//1nd-order is ok?
	if(abs(g2)==0)
	{
		dcom tmp = (k*2.0*Constant::i)*exp(Constant::i*g1*omega)*omega/R;
		this->Ax[p->t_bin] += tmp*f1.x;
		this->Ay[p->t_bin] += tmp*f1.y;
		this->Az[p->t_bin] += tmp*f1.z;
	}
	else
	{	
		dcom tmp =(k*2.0)*exp(Constant::i*g1*omega)/(g2*g2*omega*R);
		this->Ax[p->t_bin] += tmp*(f2.x*g2*omega*cos(g2*omega/R)-R*(f2.x-Constant::i*f1.x*g2*omega)*sin(g2*omega/R));
		this->Ay[p->t_bin] += tmp*(f2.y*g2*omega*cos(g2*omega/R)-R*(f2.y-Constant::i*f1.y*g2*omega)*sin(g2*omega/R));
		this->Az[p->t_bin] += tmp*(f2.z*g2*omega*cos(g2*omega/R)-R*(f2.z-Constant::i*f1.z*g2*omega)*sin(g2*omega/R));
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

   		string fmt = "Domain::OnCalculate: [%5.1f%%] Time[%6.2f min "; 

   		double prec = (step+1.0)/MaxStep*100;
   		int n = prec*20.0/100.0;
   		for(int i=0; i<n;	i++) fmt+="*";
   		for(int i=0; i<20-n;i++) fmt+="-";
   		fmt+=" %6.2f min]"; 
   		
   		 Log(fmt.c_str(), prec, delay/60, delay/60/prec*(100-prec));
   		DLog(fmt.c_str(), prec, delay/60, delay/60/prec*(100-prec)); 
	}

	//output;
	if((step+1)%Out_dt==0)
	{
		Log("Domain::OnCalculate: Output [%d].",++n_out);
		this->Output(n_out);
	}
}
