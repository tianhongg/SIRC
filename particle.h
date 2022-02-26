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
	friend class Pixel;


private:

	Domain *p_domain() {return Domain::p_Domain;}; 

	
	
	double weight;     // how many real unit charges;

	int start;	   //trajectory start step;

	ULONG NStep; 	   //size of steps
	ULONG Current_Step;  //current integration step

	Vec3* Position;
	Vec3* Velocity;

	ULONG t_bin; //time bin it belongs to
	double xi;	//x-ct

public:

	double q2m;
	int type;

	int  Load(hid_t file_id, ULONG i);
	void Normalize();
	bool IsInFrame(int step);

	Particle();
  	virtual ~Particle();


};

class Electron : public Particle 
{

public:

	Electron();
	~Electron();

};