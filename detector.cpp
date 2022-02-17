//----------------------------------------------------------------------------------||
//-------------------                  detector.cpp              -------------------||
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

Detector::Detector(char * infile): NList("Detector") 
{


	AddEntry((char*)"ThetaY_Max", 	&Theta_Y_Max,	1.0);
	AddEntry((char*)"ThetaZ_Max", 	&Theta_Z_Max,	1.0);
	AddEntry((char*)"OmegaMax", 	&Omega_Max,		20.);
	AddEntry((char*)"OmegaMin", 	&Omega_Min,		0.1);
	AddEntry((char*)"ThetaY_Grid", 	&N_Theta_Y,		1);
	AddEntry((char*)"ThetaZ_Grid", 	&N_Theta_Z,		1);
	AddEntry((char*)"OmegaGrid", 	&N_Omega,		10);
	AddEntry((char*)"Energy_Scale", &Omega_Scale,   0);
	

}

Detector::~Detector()
{

	delete[] Pixels;
 
}