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

	ULONG Size_P; //size of particle objs
	double BunchXiMax;
	double BunchXiMin;
	int N_Time;	//time bin sizes
 
private:

	int Rank;
	int	ReadType;
  	double dt; 
  	double lambda_L;

  	int MaxStep;
  	int step;
  	
  	int MovingFrame;
  	int InputType;
  	int Normalization;
  	int IntegrateOrder;
  	int IncludePartI;

  	list<Particle*> Particles;
  	Detector* MyDetector;

  	//misc: for tic function;
	int n_tick;
	int s_tick;
	chrono::time_point<std::chrono::system_clock> tic;

public: 
	void Run();
	int  ReadTrajectory();
	void OnCalculate();
	void Tick();
	void Output();
	void ReduceBunchSize();
	void TagParticles();

	double GetDt() {return dt;}
	double GetTime() {return dt*step;}
	int    GetStep()   {return step;}
	int    GetMaxStep()   {return MaxStep;}
	int    GetIntegrateOrder() {return IntegrateOrder;}
	bool   IsMovingFrame()   {return MovingFrame   == 1? true:false;}
	bool   IfNormalization() {return Normalization == 1? true:false;}
	bool   IfInputIsP() {return InputType == 1? true:false;}

	Domain(char *infile, int rank);  
	~Domain();       


};