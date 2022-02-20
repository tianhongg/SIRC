//----------------------------------------------------------------------------------||
//-------------------                  particle.cpp              -------------------||
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

void Particle::Normalize()
{

	double dt = p_domain()->GetDt();

	for(int i=0; i<NStep; i++)
	{
		Position[i] *= (2*Constant::PI);
		Velocity[i] /= (sqrt(1.0+Velocity[i].abs2()));

		if(p_domain()->IsMovingFrame())
		{
			xi = Position[i].x;
			p_domain()->BunchXiMax = max(p_domain()->BunchXiMax,xi);
			p_domain()->BunchXiMin = min(p_domain()->BunchXiMin,xi);
			//convert to lab frame
			Position[i].x = (start+i)*dt + xi;
		}
		else
		{
			xi =Position[i].x-(start+i)*dt;
			p_domain()->BunchXiMax = max(p_domain()->BunchXiMax,xi);
			p_domain()->BunchXiMin = min(p_domain()->BunchXiMin,xi);

		}
	}

}

bool Particle::IsInFrame(int step)
{
	this->Current_Step = step - start;

	// to ensure the derivative
	// Current_Step will be > 0 and < NStep-1;
	return (this->Current_Step>0&&this->Current_Step<NStep-1);

}
Particle::Particle()
{
	weight=1.0;
	start=0; //int: starting time step
	Position=NULL;
	Velocity=NULL;
	t_bin=0;
}


Particle::~Particle()
{
	delete[] Position;
	delete[] Velocity;

}

Electron::Electron() :Particle()
{
	q2m=1;
	type = ELECTRON;
}
Electron::~Electron()
{;};



