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

	n.x = cos(ty)*cos(tz)/(1.0+sin(ty)*sin(ty)+sin(tz)*sin(tz));
	n.y = 		  sin(ty)/(1.0+sin(ty)*sin(ty)+sin(tz)*sin(tz));
	n.z =		  sin(tz)/(1.0+sin(ty)*sin(ty)+sin(tz)*sin(tz));

	N_Time = timebin;

	Ax = new dcom[N_Time];
	Ay = new dcom[N_Time];
	Az = new dcom[N_Time];

	A1x[0]=A1x[1]=0.0;
	A1y[0]=A1y[1]=0.0;
	A1z[0]=A1z[1]=0.0;
	
}


Pixel::~Pixel()
{
	delete[] Ax;
	delete[] Ay;
	delete[] Az;

}