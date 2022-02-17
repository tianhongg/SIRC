//----------------------------------------------------------------------------------||
//-------------------                   SIRC.h                   -------------------||
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

#include <mpi.h>
#include <complex>
#include <cmath>
#include <time.h> 
#include <string>
#include <list>
#include <array>

#include "struct.h"
#include "namelist.h"
#include "particle.h"
#include "domain.h"
#include "pixel.h"
#include "detector.h"


#define dcom std::complex<double>

namespace Constant
{
	const dcom ci(0.0,1.0);
	const double PI =  3.1415926536;
}

namespace GlobalVars 
{	
	extern int Rank;
	extern FILE *LogFile;

}

void Log0(bool all, const char * fmt,...);
void LogDebug(const char * fmt,...);

#define ALog(fmt, ...)  LogAll(true,fmt,##__VA_ARGS__)
#define Log(fmt, ...)  Log0(false,fmt,##__VA_ARGS__)
#define DLog(fmt, ...) LogDebug(fmt,##__VA_ARGS__)



