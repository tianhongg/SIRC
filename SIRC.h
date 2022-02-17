//----------------------------------------------------------------------------------||
//-------------------                   SIRC.h                   -------------------||
//----------------------------------------------------------------------------------||
//                               ____ ___ ____   ____                               ||
//                              / ___|_ _|  _ \ / ___|								||
//                              \___ \| || |_) | |    								||
//                               ___) | ||  _ <| |___ 								||
//                              |____/___|_| \_\\____|								||
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

#include <mpi.h>
#include <complex>
#include <cmath>
#include <time.h> 
#include "namelist.h"
#include "domain.h"
#include "detector.h"
#include "pixel.h"



#define dcom std::complex<double>
#define ci dcom(0.0,1.0)
#define PI 3.1415926536


enum class dim
{
	x=0,
	y,
	z
};

namespace GlobalVars 
{	
	extern int Rank;
	extern FILE *LogFile;

}

void Log0(const char * fmt,...);
void LogDebug(const char * fmt,...);

#define Log(fmt, ...)  Log0(fmt,##__VA_ARGS__)
#define DLog(fmt, ...) LogDebug(fmt,##__VA_ARGS__)



