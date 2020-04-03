
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "kernel.cuh"
#include <stdio.h>
#include <malloc.h>
#include <math.h>


#define PI 3.141592654f





void buildRing(float * xRingCoordinates0, float * yRingCoordinates0, float * zRingCoordinates0, float spatialRadius, int nRingCoordinates) {
	float angleStep = (PI * 2.0f) / (float)nRingCoordinates;
	for (int angleIter1 = 0; angleIter1 < nRingCoordinates; angleIter1++) {
		for (int angleIter2 = 0; angleIter2 < nRingCoordinates; angleIter2++) {
			int idx = angleIter1 *nRingCoordinates + angleIter2;
			xRingCoordinates0[idx] = spatialRadius * sin(angleStep * angleIter1) * cos(angleStep * angleIter2);
			yRingCoordinates0[idx] = spatialRadius * sin(angleStep * angleIter1) * sin(angleStep * angleIter2);
			zRingCoordinates0[idx] = spatialRadius * cos(angleStep * angleIter1);
		}
	}
}


__device__ int getIdx3(int x, int y, int z, int t, int width, int height, int depth) {
	int pt = t * width * height *depth;
	int pf = z * width * height + y * width + x; // position within a frame
	return pt + pf;
}


//fixed blocksize 512
__global__ void gradKernel(float * pixels, float *GxArray, float *GyArray, float *GzArray, int width, int height, int depth, int nTimePoints, int doGradSmooth, int mR, int xR, int yR, int zR)
{
	int pixelIdx = blockIdx.x * blockDim.x + threadIdx.x;
	if (pixelIdx >= width*height*depth*nTimePoints)
		return;

	int x = pixelIdx % width;
	int y = (pixelIdx / width) % height;
	int z = (pixelIdx / (width*height)) % depth;
	int t = (pixelIdx / (width*height*depth)) % nTimePoints;

	int x_, y_, z_;
	float v;

	if (doGradSmooth == 0) {
		int x0 = fmaxf(x - 1, 0);
		int x1 = fminf(x + 1, width - 1);
		int y0 = fmaxf(y - 1, 0);
		int y1 = fminf(y + 1, height - 1);
		int z0 = fmaxf(z - 1, 0);
		int z1 = fminf(z + 1, depth - 1);
		int idx = getIdx3(x, y, z, t, width, height, depth);
		GxArray[idx] = -pixels[getIdx3(x0, y, z, t, width, height, depth)] + pixels[getIdx3(x1, y, z, t, width, height, depth)];
		GyArray[idx] = -pixels[getIdx3(x, y0, z, t, width, height, depth)] + pixels[getIdx3(x, y1, z, t, width, height, depth)];
		GzArray[idx] = -pixels[getIdx3(x, y, z0, t, width, height, depth)] + pixels[getIdx3(x, y, z1, t, width, height, depth)];
		return;
	}


	// Calculate Gx or Gy or Gz
	v = 0;
	for (int i = -mR; i <= mR; i++) {
		for (int j = -yR; j <= yR; j++) {
			for (int k = -zR; k <= zR; k++) {
				x_ = fminf(fmaxf((x + i), 0), width - 1);
				y_ = fminf(fmaxf((y + j), 0), height - 1);
				z_ = fminf(fmaxf((z + k), 0), depth - 1);
				if (i < 0) {
					v -= pixels[getIdx3(x_, y_, z_, t, width, height, depth)];
				}
				else if (i > 0) {
					v += pixels[getIdx3(x_, y_, z_, t, width, height, depth)];
				}
			}
		}
	}
	GxArray[getIdx3(x, y, z, t, width, height, depth)] = v;
	v = 0;
	for (int j = -mR; j <= mR; j++) {
		for (int i = -xR; i <= xR; i++) {
			for (int k = -zR; k <= zR; k++) {
				x_ = fminf(fmaxf((x + i), 0), width - 1);
				y_ = fminf(fmaxf((y + j), 0), height - 1);
				z_ = fminf(fmaxf((z + k), 0), depth - 1);
				if (j < 0) {
					v -= pixels[getIdx3(x_, y_, z_, t, width, height, depth)];
				}
				else if (j > 0) {
					v += pixels[getIdx3(x_, y_, z_, t, width, height, depth)];
				}
			}
		}
	}
	GyArray[getIdx3(x, y, z, t, width, height, depth)] = v;
	v = 0;
	for (int k = -mR; k <= mR; k++) {
		for (int i = -xR; i <= xR; i++) {
			for (int j = -yR; j <= yR; j++) {
				x_ = fminf(fmaxf((x + i), 0), width - 1);
				y_ = fminf(fmaxf((y + j), 0), height - 1);
				z_ = fminf(fmaxf((z + k), 0), depth - 1);
				if (k < 0) {
					v -= pixels[getIdx3(x_, y_, z_, t, width, height, depth)];
				}
				else if (k > 0) {
					v += pixels[getIdx3(x_, y_, z_, t, width, height, depth)];
				}
			}
		}
	}
	GzArray[getIdx3(x, y, z, t, width, height, depth)] = v;
}


__device__ float cubic(float xx) {
	float a = 0.5f; // Catmull-Rom interpolation
	float x = xx;
	if (x < 0.0f) x = -x;
	float z = 0.0f;
	if (x < 1.0f)
		z = x * x * (x * (-a + 2.0f) + (a - 3.0f)) + 1.0f;
	else if (x < 2.0f)
		z = -a * x * x * x + 5.0f * a * x * x - 8.0f * a * x + 4.0f * a;
	return z;
}


__device__ float interpolateG(float x, float y, float z, int t, int mag, int width, int height, int depth, float * Array) {

	x = x / mag;
	y = y / mag;
	z = z / mag;
	if (x < 1.5f || x > width - 1.5f || y < 1.5f || y > height - 1.5f || z < 1.5f || z > depth - 1.5f)
		return 0;

	int u0 = (int)floor(x - 0.5f);
	int v0 = (int)floor(y - 0.5f);
	int w0 = (int)floor(z - 0.5f);

	float s = 0.0f;
	for (int k = 0; k <= 3; k++) {
		int w = w0 - 1 + k;
		float q = 0.0f;
		for (int j = 0; j <= 3; j++) {
			int v = v0 - 1 + j;
			float p = 0.0f;
			for (int i = 0; i <= 3; i++) {
				int u = u0 - 1 + i;
				p = p + Array[getIdx3(u, v, w, t, width, height, depth)] * cubic(x - (u + 0.5f));
			}
			q = q + p * cubic(y - (v + 0.5f));
		}
		s = s + q * cubic(z - (w + 0.5f));
	}
	return s;
}




__device__ float _calculateDk(float x, float y, float z, float Xc, float Yc, float Zc, float Gx, float Gy, float Gz, float Gx2Gy2Gz2, float spatialRadius) {
	
	float tt = (Gx * (x - Xc) + Gy * (y - Yc) + Gz * (z - Zc)) / (Gx2Gy2Gz2 * Gx2Gy2Gz2);
	float xc = Gx * tt + Xc;
	float yc = Gy * tt + Yc;
	float zc = Gz * tt + Zc;
	float Dk = sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc) + (z - zc) * (z - zc));
	return Dk;

}


__global__ void calculateRadiality(float * xRingCoordinates0, float *yRingCoordinates0, float *zRingCoordinates0, int nRingCoordinates, float spatialRadius, float * radArray, int widthM, int heightM, int depthM, int nTimePoints,
	int mag, float * GxArray, float * GyArray, float * GzArray, int renormalize, int radialityPositivityConstraint) {

	// Position variables in Radiality magnified space
	int pixelIdx = blockIdx.x * blockDim.x + threadIdx.x;
	if (pixelIdx >= widthM*heightM*depthM)
		return;

	int X = pixelIdx % widthM;
	int Y = (pixelIdx / widthM) % heightM;
	int Z = (pixelIdx / (widthM*heightM)) % depthM;

	for (int t = 0; t < nTimePoints; t++) {
		float x0, y0, z0, Gx, Gy, Gz;


		float Xc = X + 0.5f;
		float Yc = Y + 0.5f;
		float Zc = Z + 0.5f;

		// Gradient Variables
		float GMag = 0;
		float GdotR = 0;

		// Radiality Variables
		float Dk, DivDFactor = 0;

		// Output
		float CGH = 0;
		for (int sampleIter1 = 0; sampleIter1 < nRingCoordinates; sampleIter1++) {
			for (int sampleIter2 = 0; sampleIter2 < nRingCoordinates; sampleIter2++) {

				int isDiverging = -1;
				// Sample (x, y, z) position
				x0 = Xc + xRingCoordinates0[sampleIter1 * nRingCoordinates + sampleIter2];
				y0 = Yc + yRingCoordinates0[sampleIter1 * nRingCoordinates + sampleIter2];
				z0 = Zc + zRingCoordinates0[sampleIter1 * nRingCoordinates + sampleIter2];

				Gx = interpolateG(x0, y0, z0, t, mag, widthM / mag, heightM / mag, depthM / mag, GxArray);
				Gy = interpolateG(x0, y0, z0, t, mag, widthM / mag, heightM / mag, depthM / mag, GyArray);
				Gz = interpolateG(x0, y0, z0, t, mag, widthM / mag, heightM / mag, depthM / mag, GzArray);
				GMag = sqrt(Gx * Gx + Gy * Gy + Gz * Gz);
				GdotR = (Gx * xRingCoordinates0[sampleIter1 * nRingCoordinates + sampleIter2] + Gy * yRingCoordinates0[sampleIter1 * nRingCoordinates + sampleIter2] + Gz * zRingCoordinates0[sampleIter1 * nRingCoordinates + sampleIter2]) / (GMag * spatialRadius);

				if (GdotR > 0 || GdotR != GdotR)
					isDiverging = 1;

				// Perpendicular distance from (Xc,Yc) to gradient line through (x,y)
				Dk = 1.0f - _calculateDk(x0, y0, z0, Xc, Yc, Zc, Gx, Gy, Gz, GMag, spatialRadius) / spatialRadius;
				Dk = Dk * Dk;
				// Accumulate Variables
				if (isDiverging == 1)
					DivDFactor -= Dk;
				else
					DivDFactor += Dk;

			}
		}
		DivDFactor = DivDFactor / nRingCoordinates;
		if (renormalize == 1)
			DivDFactor = 0.5f + (DivDFactor / 2);
		if (radialityPositivityConstraint == 1)
			CGH = fmaxf(DivDFactor, 0);
		else
			CGH = DivDFactor;
		radArray[getIdx3(X, Y, Z, t, widthM, heightM, depthM)] = CGH;
	}
}



__device__ float interpolateIntensity(float x, float y, float z, int t, int mag, int width, int height, int depth, float * pixels) {


	x = x / mag;
	y = y / mag;
	z = z / mag;

	if (x < 1.5f || x > width - 1.5f || y < 1.5f || y > height - 1.5f || z < 1.5f || z > depth - 1.5f)
		return 0;
	int u0 = (int)floor(x - 0.5f);
	int v0 = (int)floor(y - 0.5f);
	int w0 = (int)floor(z - 0.5f);


	float s = 0.0f;
	for (int k = 0; k <= 3; k++) {
		int w = w0 - 1 + k;
		float q = 0.0f;
		for (int j = 0; j <= 3; j++) {
			int v = v0 - 1 + j;
			float p = 0.0f;
			for (int i = 0; i <= 3; i++) {
				int u = u0 - 1 + i;
				p = p + pixels[getIdx3(u, v, w, t, width, height, depth)] * cubic(x - (u + 0.5f));
			}
			q = q + p * cubic(y - (v + 0.5f));
		}
		s = s + q * cubic(z - (w + 0.5f));
	}
	return fmaxf(s, 0.0f);
}



__global__ void calculateAVE(float * SRRFArray, float * radArray, float * imgd, int widthM, int heightM, int depthM, int nTimePoints, int mag, int doIntensityWeighting) {

	// Position variables in Radiality magnified space
	int pixelIdx = blockIdx.x * blockDim.x + threadIdx.x;
	if (pixelIdx >= widthM*heightM*depthM)
		return;

	int x = pixelIdx % widthM;
	int y = (pixelIdx / widthM) % heightM;
	int z = (pixelIdx / (widthM*heightM)) % depthM;


	float _x;
	float _y;
	float _z;

	float mean = 0;

	if (doIntensityWeighting == 1) {
		for (int t = 0; t < nTimePoints; t++) {
			_x = (float)x + 0.5f;
			_y = (float)y + 0.5f;
			_z = (float)z + 0.5f;
			mean += radArray[getIdx3(x, y, z, t, widthM, heightM, depthM)] * interpolateIntensity(_x, _y, _z, t, mag, widthM / mag, heightM / mag, depthM / mag, imgd);
		}
	}
	else {
		for (int t = 0; t < nTimePoints; t++) {
			mean += radArray[getIdx3(x, y, z, t, widthM, heightM, depthM)];
		}
	}

	mean = mean / nTimePoints;

	SRRFArray[z * widthM * heightM + y * widthM + x] = mean;

}




__global__ void calculateMIP(float * SRRFArray, float * radArray, float * imgd, int widthM, int heightM, int depthM, int nTimePoints, int mag, int doIntensityWeighting) {

	// Position variables in Radiality magnified space
	int pixelIdx = blockIdx.x * blockDim.x + threadIdx.x;
	if (pixelIdx >= widthM*heightM*depthM)
		return;

	int x = pixelIdx % widthM;
	int y = (pixelIdx / widthM) % heightM;
	int z = (pixelIdx / (widthM*heightM)) % depthM;


	float _x;
	float _y;
	float _z;

	float MIP = 0;

	if (doIntensityWeighting == 1) {
		for (int t = 0; t < nTimePoints; t++) {
			_x = (float)x + 0.5f;
			_y = (float)y + 0.5f;
			_z = (float)z + 0.5f;
			MIP = fmaxf(radArray[getIdx3(x, y, z, t, widthM, heightM, depthM)] * interpolateIntensity(_x, _y, _z, t, mag, widthM / mag, heightM / mag, depthM / mag, imgd), MIP);
		}
	}
	else {
		for (int t = 0; t < nTimePoints; t++) {
			MIP = fmaxf(radArray[getIdx3(x, y, z, t, widthM, heightM, depthM)], MIP);
		}
	}

	SRRFArray[z * widthM * heightM + y * widthM + x] = MIP;

}





__global__ void calculatePairwiseProductSum(float * SRRFArray, float * radArray, float * imgd, int widthM, int heightM, int depthM, int nTimePoints, int mag, int doIntensityWeighting) {

	// Position variables in Radiality magnified space
	int pixelIdx = blockIdx.x * blockDim.x + threadIdx.x;
	if (pixelIdx >= widthM*heightM*depthM)
		return;

	int x = pixelIdx % widthM;
	int y = (pixelIdx / widthM) % heightM;
	int z = (pixelIdx / (widthM*heightM)) % depthM;


	float _x;
	float _y;
	float _z;


	float pps = 0, r0 = 0, r1 = 0;
	float IW = 0;
	int counter = 0;


	for (int t0 = 0; t0 < nTimePoints; t0++) {
		r0 = fmaxf(radArray[getIdx3(x, y, z, t0, widthM, heightM, depthM)], 0);
		if (r0 > 0) {
			for (int t1 = t0; t1 < nTimePoints; t1++) {
				r1 = fmaxf(radArray[getIdx3(x, y, z, t1, widthM, heightM, depthM)], 0);
				pps += r0 * r1;
				counter++;
			}
		}
		else counter += nTimePoints - t0;
	}
	pps /= fmaxf(counter, 1);


	if (doIntensityWeighting == 1) {
		for (int t = 0; t < nTimePoints; t++) {
			_x = (float)x + 0.5f;
			_y = (float)y + 0.5f;
			_z = (float)z + 0.5f;
			IW += interpolateIntensity(_x, _y, _z, t, mag, widthM / mag, heightM / mag, depthM / mag, imgd) / nTimePoints;
		}
		pps *= IW;
	}
	SRRFArray[z * widthM * heightM + y * widthM + x] = pps;

}



__global__ void calculateACRF(float * SRRFArray, float * radArray, float * imgd, int widthM, int heightM, int depthM, int nTimePoints, int mag, int doIntensityWeighting, int doIntegrateLagTimes, int SRRForder) {

	// Position variables in Radiality magnified space
	int pixelIdx = blockIdx.x * blockDim.x + threadIdx.x;
	if (pixelIdx >= widthM*heightM*depthM)
		return;

	int x = pixelIdx % widthM;
	int y = (pixelIdx / widthM) % heightM;
	int z = (pixelIdx / (widthM*heightM)) % depthM;


	float _x;
	float _y;
	float _z;


	float mean = 0;
	float IW = 0;
	int t = 0;
	int tBin = 0;
	int nBinnedTimePoints = nTimePoints;
	float ABCD = 0;
	float ABC = 0;
	float AB = 0;
	float CD = 0;
	float AC = 0;
	float BD = 0;
	float AD = 0;
	float BC = 0;
	float A = 0;
	float B = 0;
	float C = 0;
	float D = 0;

	for (int ti = 0; ti < nTimePoints; ti++) {
		_x = (float)x + 0.5f;
		_y = (float)y + 0.5f;
		_z = (float)z + 0.5f;
		IW += radArray[getIdx3(x, y, z, ti, widthM, heightM, depthM)] * interpolateIntensity(_x, _y, _z, ti, mag, widthM / mag, heightM / mag, depthM / mag, imgd);
		mean += radArray[getIdx3(x, y, z, ti, widthM, heightM, depthM)];
	}
	mean = mean / nTimePoints;
	IW = IW / nTimePoints;


	if (doIntegrateLagTimes != 1) {
		t = 0;
		ABCD = 0;
		ABC = 0;
		AB = 0;
		CD = 0;
		AC = 0;
		BD = 0;
		AD = 0;
		BC = 0;
		while (t < (nTimePoints - SRRForder)) {
			AB += ((radArray[getIdx3(x, y, z, t, widthM, heightM, depthM)] - mean) * (radArray[getIdx3(x, y, z, t + 1, widthM, heightM, depthM)] - mean));
			if (SRRForder == 3) {
				ABC += ((radArray[getIdx3(x, y, z, t, widthM, heightM, depthM)] - mean) * (radArray[getIdx3(x, y, z, t + 1, widthM, heightM, depthM)] - mean) * (radArray[getIdx3(x, y, z, t + 2, widthM, heightM, depthM)] - mean));
			}
			if (SRRForder == 4) {
				A = radArray[getIdx3(x, y, z, t, widthM, heightM, depthM)] - mean;
				B = radArray[getIdx3(x, y, z, t + 1, widthM, heightM, depthM)] - mean;
				C = radArray[getIdx3(x, y, z, t + 2, widthM, heightM, depthM)] - mean;
				D = radArray[getIdx3(x, y, z, t + 3, widthM, heightM, depthM)] - mean;
				ABCD += A * B * C * D;
				CD += C * D;
				AC += A * C;
				BD += B * D;
				AD += A * D;
				BC += B * C;
			}
			t++;
		}
		if (SRRForder == 3) {
			SRRFArray[z * widthM * heightM + y * widthM + x] = abs(ABC) / (float)nTimePoints;
		}
		else if (SRRForder == 4) {
			SRRFArray[z *  widthM * heightM + y * widthM + x] = abs(ABCD - AB * CD - AC * BD - AD * BC) / (float)nTimePoints;
		}
		else {
			SRRFArray[z *  widthM * heightM + y * widthM + x] = abs(AB) / (float)nTimePoints;
		}
	}
	else {

		while (nBinnedTimePoints > SRRForder) {
			t = 0;
			AB = 0;

			while (t < (nBinnedTimePoints - SRRForder)) {
				tBin = t * SRRForder;
				AB += ((radArray[getIdx3(x, y, z, t, widthM, heightM, depthM)] - mean) * (radArray[getIdx3(x, y, z, t + 1, widthM, heightM, depthM)] - mean));
				if (SRRForder == 3) {
					ABC += ((radArray[getIdx3(x, y, z, t, widthM, heightM, depthM)] - mean) * (radArray[getIdx3(x, y, z, t + 1, widthM, heightM, depthM)] - mean) * (radArray[getIdx3(x, y, z, t + 2, widthM, heightM, depthM)] - mean));
				}
				if (SRRForder == 4) {
					A = radArray[getIdx3(x, y, z, t, widthM, heightM, depthM)] - mean;
					B = radArray[getIdx3(x, y, z, t + 1, widthM, heightM, depthM)] - mean;
					C = radArray[getIdx3(x, y, z, t + 2, widthM, heightM, depthM)] - mean;
					D = radArray[getIdx3(x, y, z, t + 3, widthM, heightM, depthM)] - mean;
					ABCD += A * B * C * D;
					CD += C * D;
					AC += A * C;
					BD += B * D;
					AD += A * D;
					BC += B * C;
				}
				radArray[getIdx3(x, y, z, t, widthM, heightM, depthM)] = 0;
				if (tBin < nBinnedTimePoints) {
					for (int _t = 0; _t < SRRForder; _t++) {
						radArray[getIdx3(x, y, z, t, widthM, heightM, depthM)] += radArray[getIdx3(x, y, z, tBin + _t, widthM, heightM, depthM)] / (float)SRRForder;
					}
				}
				t++;
			}
			if (SRRForder == 3) {
				SRRFArray[z * widthM * heightM + y * widthM + x] += abs(ABC) / (float)nBinnedTimePoints;
			}
			else if (SRRForder == 4) {
				SRRFArray[z * widthM * heightM + y * widthM + x] += abs((ABCD - AB * CD - AC * BD - AD * BC)) / (float)nBinnedTimePoints;
			}
			else {
				SRRFArray[z * widthM * heightM + y * widthM + x] += abs(AB) / (float)nBinnedTimePoints;
			}
			nBinnedTimePoints = nBinnedTimePoints / SRRForder;
		}
	}
	if (doIntensityWeighting == 1) {
		SRRFArray[z * widthM * heightM + y * widthM + x] = IW * SRRFArray[z * widthM * heightM + y * widthM + x];
	}
}



//imageDim: x, y, z, t
extern "C" int SRRFGradientCalcGPU(float * inImage, int * imageDim, float * outImage, SRRFConf * conf)
{
	int width = imageDim[0];
	int height = imageDim[1];
	int depth = imageDim[2];
	int nTimePoints = imageDim[3];
	int widthHeightDepth = width * height * depth;
	int widthHeightDepthTime = widthHeightDepth * nTimePoints;

	int nRingCoordinates = conf->AxesinRing * 2;
	float spatialRadius = conf->RingRadius;
	int mag = conf->Magnification;


	float * xRingCoordinates0, *yRingCoordinates0, *zRingCoordinates0;
	xRingCoordinates0 = (float *)malloc(nRingCoordinates * nRingCoordinates * sizeof(float));
	yRingCoordinates0 = (float *)malloc(nRingCoordinates * nRingCoordinates * sizeof(float));
	zRingCoordinates0 = (float *)malloc(nRingCoordinates * nRingCoordinates * sizeof(float));
	buildRing(xRingCoordinates0, yRingCoordinates0, zRingCoordinates0, spatialRadius, nRingCoordinates);
	//copy and free
	float * xRingCoordinates0d, *yRingCoordinates0d, *zRingCoordinates0d;
	cudaMalloc((void**)&xRingCoordinates0d, nRingCoordinates * nRingCoordinates * sizeof(float));
	cudaMalloc((void**)&yRingCoordinates0d, nRingCoordinates * nRingCoordinates * sizeof(float));
	cudaMalloc((void**)&zRingCoordinates0d, nRingCoordinates * nRingCoordinates * sizeof(float));
	cudaMemcpy(xRingCoordinates0d, xRingCoordinates0, nRingCoordinates * nRingCoordinates * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(yRingCoordinates0d, yRingCoordinates0, nRingCoordinates * nRingCoordinates * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(zRingCoordinates0d, zRingCoordinates0, nRingCoordinates * nRingCoordinates * sizeof(float), cudaMemcpyHostToDevice);
	free(xRingCoordinates0);
	free(yRingCoordinates0);
	free(zRingCoordinates0);


	//copy in image
	float * imgd, *GxArrayd, *GyArrayd, *GzArrayd;
	cudaMalloc((void**)&imgd, widthHeightDepthTime * sizeof(float));
	cudaMemcpy(imgd, inImage, widthHeightDepthTime * sizeof(float), cudaMemcpyHostToDevice);
	cudaMalloc((void**)&GxArrayd, widthHeightDepthTime * sizeof(float));
	cudaMalloc((void**)&GyArrayd, widthHeightDepthTime * sizeof(float));
	cudaMalloc((void**)&GzArrayd, widthHeightDepthTime * sizeof(float));
	cudaMemset(GxArrayd, 0, widthHeightDepthTime * sizeof(float));
	cudaMemset(GyArrayd, 0, widthHeightDepthTime * sizeof(float));
	cudaMemset(GzArrayd, 0, widthHeightDepthTime * sizeof(float));


	//Calc grad
	int blockSize = 512;
	int gridSize = width*height*depth*nTimePoints / blockSize + 1;
	gradKernel << <gridSize, blockSize >> > (imgd, GxArrayd, GyArrayd, GzArrayd, width, height, depth, nTimePoints, conf->GradientSmooth, conf->MainRange, conf->XRange, conf->YRange, conf->ZRange);


	float * radArrayd;
	cudaMalloc((void**)&radArrayd, widthHeightDepthTime * mag * mag * mag * sizeof(float));
	cudaMemset(radArrayd, 0, widthHeightDepthTime * mag * mag * mag * sizeof(float));
	//Calc Radiality
	gridSize = width*height*depth*mag*mag*mag / blockSize + 1;
	calculateRadiality << <gridSize, blockSize >> > (xRingCoordinates0d, yRingCoordinates0d, zRingCoordinates0d, nRingCoordinates, spatialRadius, radArrayd, width*mag, height*mag, depth*mag, nTimePoints, mag, GxArrayd, GyArrayd, GzArrayd, conf->Renormalize, conf->RemoveCon);

	float * SRRFArrayd;
	cudaMalloc((void**)&SRRFArrayd, widthHeightDepth * mag * mag * mag * sizeof(float));
	cudaMemset(SRRFArrayd, 0, widthHeightDepth * mag * mag * mag * sizeof(float));
	//Calc SRRF
	//0:TRM
	//1:TRA
	//2:TRPPM
	//3:TRAC
	if (conf->Method == 0) {
		printf("TRM\n");
		calculateMIP << <gridSize, blockSize >> > (SRRFArrayd, radArrayd, imgd, width*mag, height*mag, depth*mag, nTimePoints, mag, conf->IntensityW);
	}
	if (conf->Method == 1) {
		printf("TRA\n");
		calculateAVE << <gridSize, blockSize >> > (SRRFArrayd, radArrayd, imgd, width*mag, height*mag, depth*mag, nTimePoints, mag, conf->IntensityW);
	}
	if (conf->Method == 2) {
		printf("TRPPM\n");
		calculatePairwiseProductSum << <gridSize, blockSize >> > (SRRFArrayd, radArrayd, imgd, width*mag, height*mag, depth*mag, nTimePoints, mag, conf->IntensityW);
	}
	if (conf->Method == 3) {
		printf("TRAC\n");
		calculateACRF << <gridSize, blockSize >> > (SRRFArrayd, radArrayd, imgd, width*mag, height*mag, depth*mag, nTimePoints, mag, conf->IntensityW, conf->ITC, conf->TRACOrder);
	}


	cudaFree(xRingCoordinates0d);
	cudaFree(yRingCoordinates0d);
	cudaFree(zRingCoordinates0d);

	cudaFree(imgd);
	cudaFree(GxArrayd);
	cudaFree(GyArrayd);
	cudaFree(GzArrayd);
	cudaFree(radArrayd);

	cudaMemcpy(outImage, SRRFArrayd, widthHeightDepth * mag * mag * mag * sizeof(float), cudaMemcpyDeviceToHost);
	cudaFree(SRRFArrayd);

	return 0;
}
