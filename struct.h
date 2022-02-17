//----------------------------------------------------------------------------------||
//-------------------                 struct.h                   -------------------||
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


class Vec3
{

public:

	double x;
	double y;
	double z;

	Vec3& operator/=(double d)
	{
		x /= d;
		y /= d;
		z /= d;
		return *this;
	}

	Vec3& operator*=(double d)
	{
		x *= d;
		y *= d;
		z *= d;
		return *this;
	}

	double abs()
	{
		return sqrt(x*x + y*y + z*z);
	}

	double abs2()
	{
		return (x*x + y*y + z*z);
	}

};
