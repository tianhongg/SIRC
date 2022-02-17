//----------------------------------------------------------------------------------||
//-------------------                  IO.cpp	                 -------------------||
//----------------------------------------------------------------------------------||
//								 ____ ___ ____   ____ 								||
//								/ ___|_ _|  _ \ / ___|								||
//								\___ \| || |_) | |    								||
//								 ___) | ||  _ <| |___ 								||
//								|____/___|_| \_\\____|								||
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


#include <hdf5.h>
#include "SIRC.h"

// IO

void Log0(const char * fmt,...)
{

	if(GlobalVars::Rank>0) return;

	const char *fmtn;
	string str(fmt);
	str+="\n";
	fmtn=str.c_str();
	va_list arglist;
    va_start( arglist, fmt );
    vprintf( fmtn, arglist);
    va_end(arglist);
}



void LogDebug(const char * fmt,...)
{

	#ifndef _DEBUG
		return;
	#endif

	const char *fmtn;
	string str(fmt);
	str+="\n";
	fmtn=str.c_str();
	va_list arglist;
    va_start( arglist, fmt );
    vfprintf(GlobalVars::LogFile,fmtn, arglist);
    va_end(arglist);

}
