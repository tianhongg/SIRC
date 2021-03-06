//----------------------------------------------------------------------------------||
//-------------------                  detector.h                -------------------||
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

class Domain;



class Detector : public NList
{
	friend class Domain;
	friend class Pixel;
private:

	int N_Pixel; //number of pixel

	int N_Omega; // number of omega bin
	int N_Theta_Y;
	int N_Theta_Z;
	int N_Time; //number of time bin

	int Omega_Scale;	//linear or log

	double Theta_Y_Max;
	double Theta_Z_Max;

	double Omega_Max;
	double Omega_Min;

	double* ThetaYBin;
	double* ThetaZBin;
	double* OmegaBin;

	double Distance;

	list<Pixel*> Pixels;

	Detector(char * infile);

  	~Detector();


public:
	Domain *p_domain() {return Domain::p_Domain;}; 
	void OnDeposit(Particle* p);
	void OnDeposit(Particle* p, Node which);
	void OnGatherA();

};