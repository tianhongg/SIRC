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

#include <hdf5.h>

#define ELECTRON 0
#define ION 1

class Particle
{
	friend class Domain;
	friend class Detector;

public:


	Domain *p_domain() {return Domain::p_Domain;}; 

	int type;
	double q2m;
	double weight;     // how many real unit charge;
	double start;	   //trajectory start time;
	double dt;
	ULONG NStep; 	   //size of 

private:



	Vec3* Position;
	Vec3* Velocity;



public:

	int Load(hid_t file_id, ULONG i);
	void Normalize();

	Particle();

  	virtual ~Particle();


};

class Electron : public Particle 
{

public:

	Electron();
	~Electron();

};