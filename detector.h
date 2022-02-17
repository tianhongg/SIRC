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
//-- 			 (S)imple (I)ncoherent (R)adiation (C)alculation 				  --||
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


public:


private:

	int N_Pixel;

	int N_Omega;
	int N_Theta_Y;
	int N_Theta_Z;

	int Omega_Scale;

	double Theta_Y_Max;
	double Theta_Z_Max;

	double Omega_Max;
	double Omega_Min;


	Detector(char * infile);

  	~Detector();



};