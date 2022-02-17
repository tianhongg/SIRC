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

Pixel::Pixel(double omg, double ty, double tz)
{

	omega   = omg;
	theta_y = ty;
	theta_z = tz;

	nx = cos(ty)*cos(tz)/(1.0+sin(ty)*sin(ty)+sin(tz)*sin(tz));
	ny = 		 sin(ty)/(1.0+sin(ty)*sin(ty)+sin(tz)*sin(tz));
	nz =		 sin(tz)/(1.0+sin(ty)*sin(ty)+sin(tz)*sin(tz));
	
	
}


Pixel::~Pixel()
{
	
}