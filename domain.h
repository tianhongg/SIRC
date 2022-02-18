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
class Particle;

class Domain : public NList
{

public:

	static Domain *p_Domain;

	ULONG Size_P;
 
private:


	int Rank;
	int	ReadType;
  	double dt; 
  	double Tmax;
  	double Out_dt;
  	double lambda_L;
  	double time;

  	list<Particle*> Particles;
  	Detector* MyDetector;


  	//for tic function;
  	int n_out;
	int n_tick;
	double d_tick;
	chrono::time_point<std::chrono::system_clock> tic;


  	
public: 

	void Run();
	int ReadTrajectory();
	void OnCalculate();
	void Tick();
	void Output(int n);

	Domain(char *infile, int rank);  
	~Domain();       


};