//----------------------------------------------------------------------------------||
//-------------------                   domian.h                 -------------------||
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


#include <stdio.h>
#include <iostream>
#include <fstream>
using namespace std;

class Detector;

class Domain : public NList
{

public:

	static Domain *p_Domain;
 
private:


	int Rank;
	int	ReadType;
  	double dt; 
  	double Tmax;
  	double lambda_L;

  	list<Particle*> Particles;
  	Detector* MyDetector;


  	
public: 

	void Run();
	int ReadTrajectory();

	Domain(char *infile, int rank);  
	~Domain();       


};