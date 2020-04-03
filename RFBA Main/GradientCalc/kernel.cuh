#pragma once

struct SRRFConf {
	float RingRadius;	//0.5
	int	Magnification;	//5
	int	AxesinRing;		//6
	int MainRange;
	int XRange;
	int YRange;
	int ZRange;

	int Method;		//1
	int ITC;		//1
	int TRACOrder;		//2

	int RemoveCon;		//0
	int Renormalize;	//0
	int GradientSmooth;		//0
	int IntensityW;		//1

	int blockSizeX;	//16
	int blockSizeY;	//16
	int blockSizeZ;	//16
};

extern "C" int SRRFGradientCalcGPU(float * inImage, int * imageDim, float * outImage, SRRFConf * conf);



