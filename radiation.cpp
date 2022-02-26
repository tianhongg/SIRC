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
#include "Faddeeva.h"

dcom ComputeA(double f1, double f2, double f3, double g1, double g2, double g3, double omega, int order);

void Domain::OnCalculate()
{
	Log("Domain::OnCalculate--------------------");

	tic = chrono::system_clock::now();

	//particle loop
	for(auto it = Particles.begin(); it!=Particles.end(); it++)
	{

		Particle* p = *it;

		//start and end point;
		if(MaxStep>p->start&&IncludePartI)
		{
			MyDetector->OnDeposit(p, Node::Start);
			MyDetector->OnDeposit(p, Node::End  );
		}
		
		//time loop
		for(step=0; step<MaxStep;step++)
		{
			this->Tick();
			//trajectory may have different starting time and length
			// so checck if the trajectory is in frame;
			if(p->IsInFrame(step)) MyDetector->OnDeposit(p);
		}

		MyDetector->OnGatherA();
	}

}

// start/end point deposit
void Detector::OnDeposit(Particle* p, Node which)
{
	for(auto it = Pixels.begin(); it!=Pixels.end(); it++)
		(*it)->OnDeposit(p,which);
}


//general integration
void Detector::OnDeposit(Particle* p)
{
	for(auto it = Pixels.begin(); it!=Pixels.end(); it++)
		(*it)->OnDeposit(p);

}


// start/end point deposit
void Pixel::OnDeposit(Particle* p, Node which)
{
	ULONG MaxStep = p_domain()->GetMaxStep();
	//get start and end step
	int step = which==Node::Start? 0: min(p->NStep-1,MaxStep-p->start);
	double t= p_domain()->GetDt()*(step+p->start);

	Vec3 vc = p->Velocity[step];
	Vec3 rc = p->Position[step];

	Vec3 fc = (n.Cross(n.Cross(vc)));
	double gc = t-(n.Dot(rc));
	double bc = 1-vc.Dot(n);

	dcom tmp = exp(Constant::i*omega*gc)/bc;
	CVec3 A(fc.x*tmp, fc.y*tmp, fc.z*tmp);

	if(which==Node::Start)
	{
		//overwrite to start a new particle
		this->Ax[p->t_bin] = -A.x;
		this->Ay[p->t_bin] = -A.y;
		this->Az[p->t_bin] = -A.z;
	}
	else 
	{	
		// //end point
		this->Ax[p->t_bin] += A.x;
		this->Ay[p->t_bin] += A.y;
		this->Az[p->t_bin] += A.z;
	}

}	

//general integration
void Pixel::OnDeposit(Particle* p)
{	

	// substep and refine is not used in this version, since they may not be useful
	// double D = MyDetector->Distance;

	double t  = p_domain()->GetTime();
	double dt = p_domain()->GetDt();
	int order = p_domain()->GetIntegrateOrder();
	dcom k = dt*0.5*sqrt(p->weight); // sqrt(weight)*dt/2

	// caclulate angle nn;

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
	Vec3 f1 = fc;
	Vec3 f2 = (fp-fm)/4.0;
	Vec3 f3 = (fp-fc*2.0+fm)/8.0;

	double g1 = gc;
	double g2 = (gp-gm)/4.0;
	double g3 = (gp-gc*2.0+gm)/8.0;

	//[4] calculate A
	dcom Ax = ComputeA(f1.x, f2.x, f3.x, g1, g2, g3, omega, order);
	dcom Ay = ComputeA(f1.y, f2.y, f3.y, g1, g2, g3, omega, order);
	dcom Az = ComputeA(f1.z, f2.z, f3.z, g1, g2, g3, omega, order);

	this->Ax[p->t_bin] += k*exp(Constant::i*g1*omega)*Ax;
	this->Ay[p->t_bin] += k*exp(Constant::i*g1*omega)*Ay;
	this->Az[p->t_bin] += k*exp(Constant::i*g1*omega)*Az;

}


dcom ComputeA(double f1, double f2, double f3, double g1, double g2, double g3, double omega, int order)
{

	dcom A;
	dcom I = Constant::i;

	dcom One(1.0,0.0); //careful for the sqrt

	double PI = Constant::PI;
	double g2w = g2*omega; 


	if(abs(g3)==0 || order==1)
	{
		if(abs(g2)==0) //g2=0 g3=0
		{
			dcom tmp = (2.0*I)*omega;
			A = tmp*(f1+f3/3.0);
		}
		else
		{	
			dcom tmp = 2.0/(g2*g2w*g2w);
			A = g2w*(2.0*I*f3+f2*g2w)*cos(g2w)+I*(-2.0*f3+I*f2*g2w+(f1+f3)*g2w*g2w)*sin(g2w);
			A = A*tmp;
		}
	}
	else //2nd order
	{
		dcom tmp = 0.125/sqrt(One*g3)/g3/g3*exp(-I*g2*g2*omega/4.0/g3);
		dcom e1 = I*sqrt(I*omega/g3)*(g2-2.0*g3)/2.0;
		dcom e2 = I*sqrt(I*omega/g3)*(g2+2.0*g3)/2.0;
		e1 = Faddeeva::erf(e1);
		e2 = Faddeeva::erf(e2);
			
		A = I*sqrt(I*PI/omega)*(2.0*I*f3*g3+f3*g2*g2w-2.0*f2*g3*g2w+4.0*f1*g3*g3*omega)*(e1-e2);
		A += 4.0*exp(I*(g2*g2+4.0*g3*g3)*omega/4.0/g3)*sqrt(One*g3)*(2.0*f3*g3*cos(g2w)-I*(f3*g2-2.0*f2*g3)*sin(g2w));
		A *= tmp;

	}

	return A;
}

//=====================================================================
void Domain::Tick()
{
	s_tick++;
	int N = MaxStep*Particles.size();
	double d_tick = N/200.01;
	//tick
	if(s_tick>=(n_tick)*d_tick)
	{
		n_tick++;
		chrono::time_point<std::chrono::system_clock> toc = chrono::system_clock::now();
		chrono::milliseconds diff = std::chrono::duration_cast<std::chrono::milliseconds>(toc - tic);
   		double delay = diff.count() / 1000.0;

   		string fmt = "Domain::OnCalculate: [%5.1f%%] Time[%6.2f min "; 

   		double prec = (s_tick+1.0)/N*100;
   		int n = prec*20.0/100.0;
   		for(int i=0; i<n;	i++) fmt+="*";
   		for(int i=0; i<20-n;i++) fmt+="-";
   		fmt+=" %6.2f min]"; 
   		
   		 Log(fmt.c_str(), prec, delay/60, delay/60/prec*(100-prec));
   		DLog(fmt.c_str(), prec, delay/60, delay/60/prec*(100-prec)); 
	}

}


