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

	//time loop
	for(step = 0; step<MaxStep; step++)
	{

		//should be infront of the p loop
		this->Tick();

		//particle loop
		for(auto it = Particles.begin(); it!=Particles.end(); it++)
		{
			Particle* p = *it;
			//trajectory may have different starting time and length
			// so checck if the trajectory is in frame;

			//starting point;
			if(p->start==step)
				MyDetector->OnDeposit(p, Node::Start);

			//integrate
			if(p->IsInFrame(step)==false) continue;
				MyDetector->OnDeposit(p);
		}

	}
}

// first part: end point
void Domain::OnCalculateEndPoint()
{
	Log("Domain::OnCalculateEndPoint--------------------");

	//Clean the End Point Contribution;
	for(auto it = MyDetector->Pixels.begin(); it!=MyDetector->Pixels.end(); it++)
		(*it)->CleanA1();
	
	for(auto it = Particles.begin(); it!=Particles.end(); it++)
	{
		Particle* p = *it;
		if(step<=p->start) continue;
		MyDetector->OnDeposit(p, Node::End);
	}

}


void Detector::OnDeposit(Particle* p, Node which)
{
	for(auto it = Pixels.begin(); it!=Pixels.end(); it++)
		(*it)->OnDeposit(p,which);
}


void Detector::OnDeposit(Particle* p)
{
	for(auto it = Pixels.begin(); it!=Pixels.end(); it++)
		(*it)->OnDeposit(p);

}

void Pixel::OnDeposit(Particle* p, Node which)
{
	return;
	int step = which==Node::Start? 0:min(p->NStep-1,p->Current_Step+1);

	double t  = p_domain()->GetTime();

	Vec3 vc = p->Velocity[step];
	Vec3 rc = p->Position[step];

	Vec3 fc = (n.Cross(n.Cross(vc)));
	double gc = - (n.Dot(rc)) + t;
	double bc = 1-vc.Dot(n);

	dcom tmp = exp(Constant::i*omega*gc)/bc;
	CVec3 A(fc.x*tmp, fc.y*tmp, fc.z*tmp);

	if(which==Node::Start)
	{
		this->Ax[p->t_bin]  += -A.x;
		this->Ay[p->t_bin]  += -A.y;
		this->Az[p->t_bin]  += -A.z;
	}
	else //end point
	{
		this->A1x[p->t_bin] += A.x;
		this->A1y[p->t_bin] += A.y;
		this->A1z[p->t_bin] += A.z;
	}

}	

void Pixel::OnDeposit(Particle* p,int substep, int Refine)
{	

	// substep and refine is not used in this version, since they may not be useful
	// double D = MyDetector->Distance;

	double t  = p_domain()->GetTime();
	double dt = p_domain()->GetDt();
	int order = p_domain()->GetIntegrateOrder();
	dcom k = dt*0.5*(p->weight); // weight*dt/2

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


// void Pixel::OnDeposit(Particle* p, int substep, int Refine)
// {	
// 	double dt = p_domain()->GetDt();
// 	double t  = p_domain()->GetTime();
// 	double shift = (2*substep+1-Refine)/2.0/Refine; //shift in the unit of dt; 

// 	double R = (double)Refine;
// 	dcom k = dt*0.5*(p->weight); // weight*dt/2

// 	double hm = 1 + shift;
// 	double hp = 1 - shift;

// 	//[1] calculate:   f = -nxnx\beta;
// 	Vec3 vm = p->Velocity[p->Current_Step-1]; 
// 	Vec3 vc = p->Velocity[p->Current_Step  ];
// 	Vec3 vp = p->Velocity[p->Current_Step+1];

// 	Vec3 fm = -(n.Cross(n.Cross(vm)));
// 	Vec3 fc = -(n.Cross(n.Cross(vc)));
// 	Vec3 fp = -(n.Cross(n.Cross(vp)));

// 	//[2] calculate:   g =  t-n*r;
// 	Vec3 rm = p->Position[p->Current_Step-1]; 
// 	Vec3 rc = p->Position[p->Current_Step  ];
// 	Vec3 rp = p->Position[p->Current_Step+1];

// 	double gm = - (n.Dot(rm)) + t-dt;
// 	double gc = - (n.Dot(rc)) + t;
// 	double gp = - (n.Dot(rp)) + t+dt;

// 	//[3] calculate f1 f2 f3; g1, g2, g3
// 	Vec3   f1 = fc + (fp-fm)/2.0*shift + (fp-fc*2.0+fm)/2.0*shift*shift;  // with interpolation
// 	double g1 = gc + (gp-gm)/2.0*shift + (gp-gc*2.0+gm)/2.0*shift*shift;  // with interpolation

// 	Vec3   f2 = (fp-f1)/4.0*hm/hp + (f1-fm)/4.0*hp/hm; // with interpolation
// 	double g2 = (gp-g1)/4.0*hm/hp + (g1-gm)/4.0*hp/hm; // with interpolation

// 	// Vec3 f3 = (fp-fc*2.0+fm)/8.0;
// 	// double g3 = (gp-gc*2.0+gm)/8.0;

// 	// //1nd-order is ok?
// 	if(abs(g2)==0)
// 	{
// 		dcom tmp = (k*2.0*Constant::i)*exp(Constant::i*g1*omega)*omega/R;
// 		this->Ax[p->t_bin] += tmp*f1.x;
// 		this->Ay[p->t_bin] += tmp*f1.y;
// 		this->Az[p->t_bin] += tmp*f1.z;
// 	}
// 	else
// 	{	
// 		dcom tmp =(k*2.0)*exp(Constant::i*g1*omega)/(g2*g2*omega*R);
// 		this->Ax[p->t_bin] += tmp*(f2.x*g2*omega*cos(g2*omega/R)-R*(f2.x-Constant::i*f1.x*g2*omega)*sin(g2*omega/R));
// 		this->Ay[p->t_bin] += tmp*(f2.y*g2*omega*cos(g2*omega/R)-R*(f2.y-Constant::i*f1.y*g2*omega)*sin(g2*omega/R));
// 		this->Az[p->t_bin] += tmp*(f2.z*g2*omega*cos(g2*omega/R)-R*(f2.z-Constant::i*f1.z*g2*omega)*sin(g2*omega/R));
// 	}

// }


//=====================================================================
void Domain::Tick()
{
	//tick
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
		this->OnCalculateEndPoint();
		this->Output(n_out);
	}
}


