//----------------------------------------------------------------------------------||
//-------------------                  detector.cpp              -------------------||
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

Detector::Detector(char * infile): NList("Detector") 
{

	Log("Detector--------------------");
	AddEntry((char*)"ThetaY_Max", 	&Theta_Y_Max,	1.0);
	AddEntry((char*)"ThetaZ_Max", 	&Theta_Z_Max,	1.0);
	AddEntry((char*)"OmegaMax", 	&Omega_Max,		20.);
	AddEntry((char*)"OmegaMin", 	&Omega_Min,		0.1);
	AddEntry((char*)"ThetaY_Grid", 	&N_Theta_Y,		1);
	AddEntry((char*)"ThetaZ_Grid", 	&N_Theta_Z,		1);
	AddEntry((char*)"OmegaGrid", 	&N_Omega,		10);
	AddEntry((char*)"Energy_Scale", &Omega_Scale,   0);
	AddEntry((char*)"Time_Grid",	&N_Time,		1);
	AddEntry((char*)"Distance",		&Distance,		1e5);
	
	
	Log("Detector: Read Parameters From .ini File...");
	FILE *p_File = fopen(infile,"rt");
	if (p_File)
	{
		rewind(p_File);
		read(p_File);
	}
	fclose(p_File); 


	// fix the parameters;
	N_Theta_Y = max(1,N_Theta_Y);
	N_Theta_Z = max(1,N_Theta_Z);
	N_Omega   = max(2,N_Omega);
	N_Time	  = max(1,N_Time);

	Theta_Y_Max = max(0.0,min(Constant::PI,Theta_Y_Max));
	Theta_Z_Max = max(0.0,min(Constant::PI,Theta_Z_Max));
	Omega_Min 	= max(1.0,Omega_Min);
	Omega_Max 	= max(Omega_Min+1,Omega_Max);

	//normalize
	if(p_domain()->IfNormalization())
		Distance *= 2*Constant::PI;

	//stick to odder number of grids for symmetry
	if(N_Theta_Y>1)
	{
		N_Theta_Y/=2; N_Theta_Y = N_Theta_Y*2+1;
	}
	if(N_Theta_Z>1)
	{
		N_Theta_Z/=2; N_Theta_Z = N_Theta_Z*2+1;
	}


	if(p_domain()->BunchXiMax==p_domain()->BunchXiMin) N_Time=1;
	p_domain()->N_Time = N_Time;


	// number of pixels
	N_Pixel = N_Omega*N_Theta_Y*N_Theta_Z;
	// size of pixel objects
	int Size_P = (N_Pixel*sizeof(Pixel) + sizeof(dcom)*8*N_Time)/1024;

	Log("Detector: N_Omega: [%d]; N_Theta_Y: [%d]; N_Theta_Z: [%d]",N_Omega, N_Theta_Y, N_Theta_Z);
	Log("		--- each pixel has [%d] time_bins",N_Time);
	Log("Detector: Create [%d] Pixels...",N_Pixel);
	Log("Detector: Size of All Pixel Objs: [%d Kilobytes]", Size_P);


	OmegaBin  = new double[N_Omega];
	ThetaYBin = new double[N_Theta_Y];
	ThetaZBin = new double[N_Theta_Z];

	//make the bin
	for(int i=0; i<N_Omega;i++)
	{
		if(Omega_Scale==0)
		{
			OmegaBin[i] = Omega_Min+(Omega_Max-Omega_Min)*i/(N_Omega-1);
		}
		else
		{
			OmegaBin[i] = Omega_Min*exp(log(Omega_Max/Omega_Min)/(N_Omega-1)*i);
		}
	}

	//make the bin
	if(N_Theta_Y>1)
		for(int i=0; i<N_Theta_Y;i++)
			ThetaYBin[i] = -Theta_Y_Max + 2*Theta_Y_Max/(N_Theta_Y-1)*i;
	else
		ThetaYBin[0]=0;

	//make the bin
	if(N_Theta_Z>1)
		for(int i=0; i<N_Theta_Z;i++)
			ThetaZBin[i] = -Theta_Z_Max + 2*Theta_Z_Max/(N_Theta_Z-1)*i;
	else
		ThetaZBin[0]=0;


	// create pixels
	for(int i=0; i<N_Omega;i++)
	{
		for(int j=0; j<N_Theta_Y;j++)
		{
			for(int k=0; k<N_Theta_Z;k++)
			{
				Pixels.push_back( new Pixel(OmegaBin[i], ThetaYBin[j], ThetaZBin[k], N_Time, this));
			}
		}

	}
}


void Detector::OnGatherA()
{
	for(auto it = Pixels.begin(); it!=Pixels.end(); it++)
	{
		Pixel* p = (*it);
		for(int t=0; t<N_Time; t++)
		{
			p->S1[t] += norm(p->Ay[t]) + norm(p->Az[t]);
			p->S2[t] += norm(p->Ay[t]) - norm(p->Az[t]);
			p->S3[t] +=  2.0*real((p->Ay[t])*conj(p->Az[t]));
			p->S4[t] += -2.0*imag((p->Ay[t])*conj(p->Az[t]));
			p->II[t] += norm(p->Ax[t]) + norm(p->Ay[t]) + norm(p->Az[t]);
			p->Ax[t]=p->Ay[t]=p->Az[t]=0.0;
		}
			
	}

}

Detector::~Detector()
{

	delete[] ThetaYBin;
	delete[] ThetaZBin;
	delete[] OmegaBin;

	for(auto it : Pixels)
		delete it;
	Pixels.clear();
 
}