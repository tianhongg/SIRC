//----------------------------------------------------------------------------------||
//-------------------                  pixel.cpp                 -------------------||
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

Pixel::Pixel(double omg, double ty, double tz, int timebin, Detector* pd)
{

	omega   = omg;
	theta_y = ty;
	theta_z = tz;

	MyDetector = pd;
	double tmp =sqrt(1.0+sin(ty)*sin(ty)*sin(tz)*sin(tz));
	n.x = cos(ty)*cos(tz)/tmp;
	n.y = 		  sin(ty)/tmp;
	n.z =		  sin(tz)/tmp;

	N_Time = timebin;

	//potential
	Ax  = new dcom[N_Time];
	Ay  = new dcom[N_Time];
	Az  = new dcom[N_Time];

	//stock parameters and intensity
	S1 = new double[N_Time];  
	S2 = new double[N_Time];  
	S3 = new double[N_Time];  
	S4 = new double[N_Time];  
	II = new double[N_Time];

	// 
	for(int i=0;i<N_Time;i++)
	{
		S1[i]=S2[i]=S3[i]=S4[i]=II[i]=0.0;  
		Ax[i]=Ay[i]=Ay[i]=0.0;

	}
	
}

Pixel::~Pixel()
{
	delete[] Ax;
	delete[] Ay;
	delete[] Az;
	delete[] S1;
	delete[] S2;
	delete[] S3;
	delete[] S4;
	delete[] II;

}