//----------------------------------------------------------------------------------||
//-------------------                  particle.h                -------------------||
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

#pragma once


#define ELECTRON 0
#define ION 1

class Particle
{
	friend class Domain;

public:

	int type;
	double q2m;
	double weight;  // how many real unit charge;
	double ts;		//trajectory start time;

private:



	Vec3* Position;
	Vec3* Velocity;



public:

	Particle(double weightp, double T_Start);

  	virtual ~Particle();


};

class Electron : public Particle {
	public:
	Electron(double weightp, double T_Start);
	~Electron();

};