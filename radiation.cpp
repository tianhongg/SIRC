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
		if(IncludePartI&&MaxStep>p->start)
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
	int step = which==Node::Start? 0: min(p->NStep-1,MaxStep-p->start)-1;
	double t= p_domain()->GetDt()*(step+0.5+p->start);

	Vec3 vc = (p->Velocity[step]+p->Velocity[step+1])/2.0;
	Vec3 rc = (p->Position[step]+p->Position[step+1])/2.0;

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
	int R = p_domain()->GetRefine();
	if(R>1) return OnDeposit(p,R);

	//-----------------integrate---------------
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
	Vec3 f1 =  fc;
	Vec3 f2 = (fp-fm)/4.0;
	Vec3 f3 = (fp-fc*2.0+fm)/8.0;

	double g1 =  gc;
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

//general integration with cubic interp
void Pixel::OnDeposit(Particle* p, int R)
{	
	double t  = p_domain()->GetTime();
	double dt = p_domain()->GetDt();
	int order = p_domain()->GetIntegrateOrder();
	dcom k = dt*0.5*sqrt(p->weight)/(R*1.0); // sqrt(weight)*dt/2

	// 
	Vec3 vmm, vpp, rmm, rpp, f1, f2, f3;
	double g1, g2, g3;

	//[1] calculate:   f = -nxnx\beta;
	Vec3 vm = p->Velocity[p->Current_Step-1]; 
	Vec3 vc = p->Velocity[p->Current_Step  ];
	Vec3 vp = p->Velocity[p->Current_Step+1];
	//[2] calculate:   g =  t-n*r;
	Vec3 rm = p->Position[p->Current_Step-1]; 
	Vec3 rc = p->Position[p->Current_Step  ];
	Vec3 rp = p->Position[p->Current_Step+1];

	if(p->Current_Step==1)
	{
		vmm = vm*2-vc; rmm = rm*2-rc;
	}
	else
	{
		vmm = p->Velocity[p->Current_Step-2]; rmm = p->Position[p->Current_Step-2];
	}
	if(p->Current_Step==p->NStep-2)
	{
		vpp = vp*2-vc; rpp = rp*2-rc;
	}
	else
	{
		vpp = p->Velocity[p->Current_Step+2]; rpp = p->Position[p->Current_Step+2];
	}

	Vec3 fmm = -(n.Cross(n.Cross(vmm)));
	Vec3 fm  = -(n.Cross(n.Cross(vm)));
	Vec3 fc  = -(n.Cross(n.Cross(vc)));
	Vec3 fp  = -(n.Cross(n.Cross(vp)));
	Vec3 fpp = -(n.Cross(n.Cross(vpp)));

	double gmm = - (n.Dot(rmm))+ t-dt*2;
	double gm  = - (n.Dot(rm)) + t-dt;
	double gc  = - (n.Dot(rc)) + t;
	double gp  = - (n.Dot(rp)) + t+dt;
	double gpp = - (n.Dot(rpp))+ t+dt*2;


	//==============first half
	double ag = (gc-gmm)/2.0-(gc-gm);
	double bg = (gm- gp)/2.0+(gc-gm);

	Vec3   af = (fc-fmm)/2.0-(fc-fm);
	Vec3   bf = (fm- fp)/2.0+(fc-fm);

	//cubic interp
	for(int sub=0;sub<R/2;sub++)
	{
		double t = 0.5+(2*sub+1)/2.0/R;

		g1 = gm*(1-t)+gc*t+( ag*(1-t)+bg*t )*t*(1-t);
		g2 = (gc-gm+bg*(2-3*t)*t)*R*R+(ag-bg+ag*R*R*(t-1)*(3*t-1)); g2/=2*R*R*R;
		g3 = bg*(1-3*t)+ag*(3*t-2); g3/=4*R*R;

		f1 = fm*(1-t)+fc*t+( af*(1-t)+bf*t )*t*(1-t);
		f2 = (fc-fm+bf*(2-3*t)*t)*R*R+(af-bf+af*R*R*(t-1)*(3*t-1)); f2/=2*R*R*R;
		f3 = bf*(1-3*t)+af*(3*t-2); f3/=4*R*R;

		//[4] calculate A
		dcom Ax = ComputeA(f1.x, f2.x, f3.x, g1, g2, g3, omega, order);
		dcom Ay = ComputeA(f1.y, f2.y, f3.y, g1, g2, g3, omega, order);
		dcom Az = ComputeA(f1.z, f2.z, f3.z, g1, g2, g3, omega, order);

		this->Ax[p->t_bin] += k*exp(Constant::i*g1*omega)*Ax;
		this->Ay[p->t_bin] += k*exp(Constant::i*g1*omega)*Ay;
		this->Az[p->t_bin] += k*exp(Constant::i*g1*omega)*Az;
	}

	//===============second half
	ag = (gp- gm)/2.0-(gp-gc);
	bg = (gc-gpp)/2.0+(gp-gc);
	af = (fp- fm)/2.0-(fp-fc);
	bf = (fc-fpp)/2.0+(fp-fc);

	for(int sub=0;sub<R/2;sub++)
	{
		double t = (2*sub+1)/2.0/R;

		g1 = gc*(1-t)+gp*t+( ag*(1-t)+bg*t )*t*(1-t);
		g2 = (gp-gc+bg*(2-3*t)*t)*R*R+(ag-bg+ag*R*R*(t-1)*(3*t-1)); g2/=2*R*R*R;
		g3 = bg*(1-3*t)+ag*(3*t-2); g3/=4*R*R;

		f1 = fc*(1-t)+fp*t+( af*(1-t)+bf*t )*t*(1-t);
		f2 = (fp-fc+bf*(2-3*t)*t)*R*R+(af-bf+af*R*R*(t-1)*(3*t-1)); f2/=2*R*R*R;
		f3 = bf*(1-3*t)+af*(3*t-2); f3/=4*R*R;

		//[4] calculate A
		dcom Ax = ComputeA(f1.x, f2.x, f3.x, g1, g2, g3, omega, order);
		dcom Ay = ComputeA(f1.y, f2.y, f3.y, g1, g2, g3, omega, order);
		dcom Az = ComputeA(f1.z, f2.z, f3.z, g1, g2, g3, omega, order);

		this->Ax[p->t_bin] += k*exp(Constant::i*g1*omega)*Ax;
		this->Ay[p->t_bin] += k*exp(Constant::i*g1*omega)*Ay;
		this->Az[p->t_bin] += k*exp(Constant::i*g1*omega)*Az;
	}

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

		dcom tmp = 0.25/sqrt(One*g3)/g3*exp(-I*g2*g2*omega/4.0/g3);

		dcom e1 = I*sqrt(I*omega/g3)*(g2-2.0*g3)/2.0;
		dcom e2 = I*sqrt(I*omega/g3)*(g2+2.0*g3)/2.0;

		e1 = Faddeeva::erf(e1);
		e2 = Faddeeva::erf(e2);
			
		A = I*sqrt(I*PI*omega)*(2.0*f1*g3-f2*g2)*(e1-e2);
		A += 4.0*I*exp(I*(g2*g2+4.0*g3*g3)*omega/4.0/g3)*sqrt(One*g3)*f2*sin(g2w);
		A *= tmp;

	}

	return A;
}




