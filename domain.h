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
  	int Out_dt;
  	int step;
  	int Refine; //step refinement; int
  	
  	int MovingFrame;
  	int InputType;
  	int Normalization;
  	int IntegrateOrder;

  	list<Particle*> Particles;
  	Detector* MyDetector;

  	//misc: for tic function;
  	int n_out;
	int n_tick;
	double d_tick;
	chrono::time_point<std::chrono::system_clock> tic;

public: 

	int tmpCounter=0;
	int tmpCounter2=0;
	
	void Run();
	int  ReadTrajectory();
	void OnCalculate();
	void Tick();
	void Output(int n);
	void ReduceBunchSize();
	void TagParticles();

	double GetDt() {return dt;}
	double GetTime() {return dt*step;}
	int    GetStep()   {return step;}
	int    GetRefine() {return Refine;}
	int    GetIntegrateOrder() {return IntegrateOrder;}
	bool   IsMovingFrame()   {return MovingFrame   == 1? true:false;}
	bool   IfNormalization() {return Normalization == 1? true:false;}
	bool   IfInputIsP() {return InputType == 1? true:false;}

	Domain(char *infile, int rank);  
	~Domain();       


};