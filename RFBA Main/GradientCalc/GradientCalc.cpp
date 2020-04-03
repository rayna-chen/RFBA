#include "GradientCalc.h"
#include "kernel.cuh"
#include "malloc.h"
#include "math.h"
#include "windows.h"
#include "stdio.h"




SRRFConf config;
int border = 2;


void loadConfig() {
	char cfgtemp[254];
	char * cfgFile = ".\\SRRF.ini";

	GetPrivateProfileString("SRRF", "RingRadius", "0.5", cfgtemp, 63, cfgFile);	
	config.RingRadius = (float)atof(cfgtemp);
	config.Magnification = GetPrivateProfileInt("SRRF", "Magnification", 5, cfgFile);	
	config.AxesinRing = GetPrivateProfileInt("SRRF", "AxesinRing", 6, cfgFile);		
	config.MainRange = GetPrivateProfileInt("SRRF", "MainRange", 1, cfgFile);		
	config.XRange = GetPrivateProfileInt("SRRF", "XRange", 0, cfgFile);		
	config.YRange = GetPrivateProfileInt("SRRF", "YRange", 0, cfgFile);		
	config.ZRange = GetPrivateProfileInt("SRRF", "ZRange", 0, cfgFile);		

	config.Method = GetPrivateProfileInt("Temproal", "Method", 1, cfgFile);		
	printf("SRRF Method: %d\n", config.Method);
	config.ITC = GetPrivateProfileInt("Temproal", "Integrate_Temporal_Correlations", 1, cfgFile);	
	config.TRACOrder = GetPrivateProfileInt("Temproal", "TRAC_ORDER", 2, cfgFile); 		

	config.RemoveCon = GetPrivateProfileInt("Radiality", "Remove_Positivity_Constraint", 0, cfgFile);		
	config.Renormalize = GetPrivateProfileInt("Radiality", "Renormalize", 0, cfgFile);	
	config.GradientSmooth = GetPrivateProfileInt("Radiality", "Gradient_Smoothing", 0, cfgFile);	
	config.IntensityW = GetPrivateProfileInt("Weighting", "Intensity_Weighting", 1, cfgFile);		

	//Split Block Parameters
	config.blockSizeX = GetPrivateProfileInt("Split", "BlockSizeX", 16, cfgFile);		
	config.blockSizeY = GetPrivateProfileInt("Split", "BlockSizeY", 16, cfgFile);		
	config.blockSizeZ = GetPrivateProfileInt("Split", "BlockSizeZ", 16, cfgFile);		
}




void getSplitCount(int * splitN, int * imageDim, int * maxDim) {
	for (int i = 0; i < 3; i++) {
		splitN[i] = (imageDim[i] - 2 * border) / (maxDim[i] - 2 * border);
		if ((imageDim[i] - 2 * border) % (maxDim[i] - 2 * border) != 0)
			splitN[i] = splitN[i] + 1;
	}
}

float * doSplit(float *inImage, int * imageDim, int * maxDim, int * splitN, int * posN, int * resDim) {
	int startPos[3];
	for (int i = 0; i < 3; i++) {
		if (posN[i] == splitN[i] - 1) {
			//last block
			resDim[i] = imageDim[i] - (splitN[i] - 1)*(maxDim[i] - 2 * border);
		}
		else {
			resDim[i] = maxDim[i];
		}

		//Calc Pos
		if (posN[i] == 0) {
			startPos[i] = 0;
		}
		else {
			startPos[i] = posN[i] * (maxDim[i] - 2 * border);
		}
	}
	//printf("resDim: %d\t%d\t%d\n", resDim[0], resDim[1], resDim[2]);
	//printf("startPos: %d\t%d\t%d\n", startPos[0], startPos[1], startPos[2]);

	resDim[3] = imageDim[3];
	float * result = (float *)malloc(resDim[0] * resDim[1] * resDim[2] * imageDim[3] * sizeof(float));
	//Do Copy
	for (int t = 0; t < imageDim[3]; t++)
		for (int z = 0; z < resDim[2]; z++)
			for (int y = 0; y < resDim[1]; y++)
				for (int x = 0; x < resDim[0]; x++) {
					int resIdx = t * resDim[0] * resDim[1] * resDim[2] + z* resDim[0] * resDim[1] + y* resDim[0] + x;
					int rawIdx = t * imageDim[0] * imageDim[1] * imageDim[2] + (z + startPos[2])*  imageDim[0] * imageDim[1] + (y + startPos[1])*imageDim[0] + x + startPos[0];
					result[resIdx] = inImage[rawIdx];
				}
	return(result);

}

void doPutback(float *outImage, int * imageDim, float * resImage, int * resDim, int * maxDim, int * posN) {
	int startPos[3];
	int borderM = border * config.Magnification;
	for (int i = 0; i < 3; i++) {
		//Calc Pos
		if (posN[i] == 0) {
			startPos[i] = 0;
		}
		else {
			startPos[i] = posN[i] * (maxDim[i] - 2 * border);
		}
		startPos[i] = startPos[i] * config.Magnification;
	}



	//Do Copy
	for (int z = borderM; z < resDim[2]; z++)
		for (int y = borderM; y < resDim[1]; y++)
			for (int x = borderM; x < resDim[0]; x++) {
				int resIdx = z* resDim[0] * resDim[1] + y* resDim[0] + x;
				int outIdx = (z + startPos[2])*  imageDim[0] * imageDim[1] + (y + startPos[1])*imageDim[0] + x + startPos[0];
				outImage[outIdx] = resImage[resIdx];
			}
}



int doSRRFGradientCalc(float * inImage, int * imageDim, float * outImage) {

	//Load Config
	loadConfig();

	int mag = config.Magnification;

	//Split
	int splitN[3];
	int maxDim[3];
	int posN[3];
	int resDim[4];
	int outDim[3];
	int outRes[3];
	maxDim[0] = config.blockSizeX;	maxDim[1] = config.blockSizeY;	maxDim[2] = config.blockSizeZ;
	outDim[0] = imageDim[0] * mag;
	outDim[1] = imageDim[1] * mag;
	outDim[2] = imageDim[2] * mag;


	getSplitCount(splitN, imageDim, maxDim);
	printf("BlockSize: %d\t%d\t%d\n", config.blockSizeX, config.blockSizeY, config.blockSizeZ);
	printf("Split: %d\t%d\t%d\n", splitN[0], splitN[1], splitN[2]);

	for (int k = 0; k < splitN[2]; k++)
		for (int j = 0; j < splitN[1]; j++)
			for (int i = 0; i < splitN[0]; i++) {
				posN[0] = i;
				posN[1] = j;
				posN[2] = k;
				float * stImage = doSplit(inImage, imageDim, maxDim, splitN, posN, resDim);
				outRes[0] = resDim[0] * mag;
				outRes[1] = resDim[1] * mag;
				outRes[2] = resDim[2] * mag;
				float * stOut = (float *)malloc(resDim[0] * resDim[1] * resDim[2] * mag*mag*mag * sizeof(float));
				SRRFGradientCalcGPU(stImage, resDim, stOut, &config);
				free(stImage);
				//Put back
				doPutback(outImage, outDim, stOut, outRes, maxDim, posN);
				free(stOut);
			}

	return(0);
}




//imageDim: x, y, z, t
__declspec(dllexport)  int SRRFGradientCalc(float * inImage, int * imageDim, float * outImage)
{

	return doSRRFGradientCalc(inImage, imageDim, outImage);
}


__declspec(dllexport) int SRRFGetMag() {
	return(GetPrivateProfileInt("SRRF", "Magnification", 5, "SRRF.ini"));
}