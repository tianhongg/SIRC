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

	//--------------------------
	Vec3 operator-(void)
	{
		return Vec3(-x,-y,-z);
	}

	//--------------------------
	Vec3 operator-(Vec3 v)
	{
		Vec3 tmp(x-v.x, y-v.y, z-v.z);
		return tmp;
	}

	Vec3 operator+(Vec3 v)
	{
		Vec3 tmp(x+v.x, y+v.y, z+v.z);
		return tmp;
	}

	//--------------------------
	Vec3 operator*(double d)
	{
		Vec3 tmp(x*d, y*d, z*d);
		return tmp;
	}

	Vec3 operator/(double d)
	{
		Vec3 tmp(x/d, y/d, z/d);
		return tmp;
	}

	//--------------------------
	Vec3& operator-=(double d)
	{
		x -= d;
		y -= d;
		z -= d;
		return *this;
	}

	Vec3& operator+=(double d)
	{
		x += d;
		y += d;
		z += d;
		return *this;
	}

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

	//--------------------------
	double abs2()
	{
		return (x*x + y*y + z*z);
	}

	double Dot(Vec3 v)
	{
		return x*v.x+y*v.y+z*v.z;
	}

	Vec3 Cross(Vec3 v)
	{
		Vec3 tmp; 
		tmp.x = y*v.z - z*v.y;
		tmp.y = z*v.x - x*v.z;
		tmp.z = x*v.y - y*v.x;
		return tmp;
	}

	//--------------------------
	Vec3(){;};
	Vec3(double x0, double y0, double z0)
	{
		x=x0;
		y=y0;
		z=z0;
	}

};



class CVec3
{

public:

	dcom x;
	dcom y;
	dcom z;


	//--------------------------
	CVec3(){;};
	CVec3(dcom x0, dcom y0, dcom z0)
	{
		x=x0;
		y=y0;
		z=z0;
	}

};
