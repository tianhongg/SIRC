//----------------------------------------------------------------------------------||
//-------------------                  pixel.h 	                 -------------------||
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

class Detector;
class Domain;



class Pixel
{
	friend class Detector;
	friend class Domain;


public:


private:

	double theta_y;
	double theta_z;
	double omega;

	Detector* MyDetector;

	//Unit vector;
	Vec3 n;

	// number of time_bin
	int N_Time;

	
	dcom* Ax;
	dcom* Ay;
	dcom* Az;

	//for first term end points
	dcom* A1x;
	dcom* A1y;
	dcom* A1z;

	Domain *p_domain() {return Domain::p_Domain;}; 
	Pixel(double omg, double ty, double tz, int timebin, Detector* pd);
	void OnDeposit(Particle* p, int substep=0, int Refine=1);
	void OnDeposit(Particle* p, Node which);
	void CleanA1();

  	~Pixel();



};