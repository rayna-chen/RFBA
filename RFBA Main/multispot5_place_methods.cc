/*
   This file is part of B-cubed.

   Copyright (C) 2009, 2010, 2011, Edward Rosten and Susan Cox

   B-cubed is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 3.0 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>
*/

#include "tag/printf.h"
#undef make_tuple

#include <climits>
#include <algorithm>
#include <tuple>
#include "cvd/image.h"
#include "cvd/vector_image_ref.h"
#include "cvd/byte.h"
#include "gvars3/instances.h"
#include "TooN/TooN.h"

#include "multispot5.h"
#include "forward_algorithm.h"
#include "storm.h"
#include "debug.h"

using namespace std;
using namespace CVD;
using namespace GVars3;
using namespace TooN;
using namespace std::tr1;


extern int splitMode;

//RFBA: SRRF Using CUDA GPU Acceraltion, Added@20180306
#pragma comment(lib,"GradientCalc.lib")
_declspec(dllimport) int SRRFGradientCalc(float * inImage, int * imageDim, float * outImage);
_declspec(dllimport) int SRRFGetMag();

//RFBA: Update Parameters for 3D Image Computing, Modified@20180306
tuple<StateParameters, double, double> static get_defaults(const vector<Ref3D>& pixels)
{
	const double variance = 1; // it should be

							   //To scale the X axis of a log-normal distribution, only
							   //the mu parameter needs to be changed...
	const double intensity_mu = GV3::get<double>("intensity.rel_mu", 0., -1) + log(sqrt(variance));
	const double intensity_sigma = GV3::get<double>("intensity.rel_sigma", 0., -1);
	const double blur_mu = GV3::get<double>("blur.mu", 0., -1);
	const double blur_sigma = GV3::get<double>("blur.sigma", 0., -1);


	const double intensity = log_normal_mode(intensity_mu, intensity_sigma);
	const double blur = log_normal_mode(blur_mu, blur_sigma);

	StateParameters p;

	//Initialize the MT19937 RNG from a seed.
	p.rng = shared_ptr<MT19937>(new MT19937);
	p.rng->simple_seed(GV3::get<int>("seed", 0, 1));

	//Start at the beginning
	p.pass = 0;
	p.iteration = 0;

	//Remember which pixels are in use
	p.pixels = pixels;

	return make_tuple(p, intensity, blur);
}


//Function for accumulating an integer (counting)
static void acc(int & i, const Vector<2>&)
{
	i++;
}

//Function for accumulating a container (inserting)
static void acc(vector<Vector<2> >& i, const Vector<2> &r)
{
	i.push_back(r);
}

//Function for placing spots over a hexagonal grid 
template<class Ret>
Ret place_spots(double sp, Vector<2> centre, double radius, const Image<bool>& mask)
{
	double angle = M_PI / 180 * 6;
	Vector<2> a_axis = SO2<>(angle) * makeVector(1, 0);
	Vector<2> b_axis = SO2<>(M_PI / 3) * a_axis;

	Ret num = Ret();
	//The range is +- 2*r / sqrt(3) on each axis.
	//To prove:
	//Draw a curcle with a horizontal line through it,
	//brushing the top and brushing the bottom (the a axis).
	//
	//Now draw the same three lines, but rotated round by 60 degrees.
	//(the b axis).
	//
	//A bunch of 30/60/90 triangles with an opposite length of r
	//are formed. We need the hypotenuse which is 2r/sqrt(3).

	int n = (int)ceil(2 * radius / sqrt(3) / sp);
	for (int a = -n; a <= n; a++)
		for (int b = -n; b <= n; b++)
		{
			Vector<2> pos = centre + a*sp*a_axis + b*sp*b_axis;
			ImageRef p = ir(pos + makeVector(.5, .5));

			if (mask.in_image(p) && mask[p])
				acc(num, pos);
		}

	return num;
}

vector<Vector<2> > find_spacing(int target, const Image<bool>& mask)
{
	//First a bounding circle is required. The circle need not be tight.

	//Compute the approximate bounding circle
	//First compute the centroid
	Vector<2> centre = Zeros;
	int count = 0;
	for (int y = 0; y < mask.size().y; y++)
		for (int x = 0; x < mask.size().x; x++)
			if (mask[y][x])
			{
				centre += makeVector(x, y);
				count++;
			}
	centre /= count;

	double r2 = 0;
	//Now compute the radius
	for (int y = 0; y < mask.size().y; y++)
		for (int x = 0; x < mask.size().x; x++)
			if (mask[y][x])
				r2 = max(r2, norm_sq(makeVector(x, y) - centre));
	double radius = r2;

	//Perform a binary search to find an appropriate spacing. The function
	//is not monotonic, so the spacing is not guaranteed to be unique.
	double close = 0;
	int large_num = INT_MAX;

	double far = sqrt(mask.size().mag_squared());
	int small_num = place_spots<int>(far, centre, radius, mask);

	if (target > small_num)
		while (small_num != large_num && far - close > 1e-6)
		{
			double mid = (close + far) / 2;
			int mid_num = place_spots<int>(mid, centre, radius, mask);

			if (mid_num > target)
			{
				large_num = mid_num;
				close = mid;
			}
			else
			{
				small_num = mid_num;
				far = mid;
			}
		}

	//Pick the best, in case the algorithm terminated due to 
	//a too small disparity.
	double spacing;
	if (large_num - target < target - small_num)
		spacing = close;
	else
		spacing = far;

	//Use small_num or close as the spacing
	return place_spots<vector<Vector<2> > >(spacing, centre, radius, mask);
}


//RFBA: Removed@20180306
/*
StateParameters place_spots_uniform(int num_spots, const vector<ImageRef>& pixels, const ImageRef& size)
{
	Image<bool> pix(size);
	pix.fill(false);
	for (unsigned int i = 0; i < pixels.size(); i++)
		pix[pixels[i]] = true;


	vector<Vector<2> > spot_pos = find_spacing(num_spots, pix);


	StateParameters p;
	double intensity_mode, blur_mode;
	tie(p, intensity_mode, blur_mode) = get_defaults(pixels);

	//Create all the spots
	for (unsigned int i = 0; i < spot_pos.size(); i++)
		p.spots.push_back(makeVector(intensity_mode, blur_mode, spot_pos[i][0], spot_pos[i][1]));


	return p;
}
*/





//RFBA: Updated for 3D Image, Add SRRF Guided@20180306
StateParameters place_spots_intensity_sampled(int num_spots, const vector<Ref3D>& pixels, const vector<Image3D<float>> & ims, int dosrrfonly)
{
	//assert_same_size(ims);

	StateParameters p;
	double intensity_mode, blur_mode;
	tie(p, intensity_mode, blur_mode) = get_defaults(pixels);

	//Calc SRRF

	int imageDim[4];
	imageDim[0] = ims[0].size().x;
	imageDim[1] = ims[0].size().y;
	imageDim[2] = ims[0].size().z;
	imageDim[3] = ims.size();

	float *inImage = (float *)malloc(imageDim[0] * imageDim[1] * imageDim[2] * imageDim[3] * sizeof(float));
	for (int w = 0; w < imageDim[3]; w++)
		for (int k = 0; k < imageDim[2]; k++)
			for (int j = 0; j < imageDim[1]; j++)
				for (int i = 0; i < imageDim[0]; i++) {
					int idx = w * imageDim[0] * imageDim[1] * imageDim[2] + k * imageDim[0] * imageDim[1] + j* imageDim[0] + i;
					int idx2 = k * imageDim[0] * imageDim[1] + j* imageDim[0] + i;
					inImage[idx] = ims[w].data()[idx2];
				}

	//Parameters
	int mag = SRRFGetMag();

	float *outImage = (float *)malloc(imageDim[0] * imageDim[1] * imageDim[2] * mag * mag *mag * sizeof(float));

	if ((splitMode == 0) || (dosrrfonly == 1)) {

		SRRFGradientCalc(inImage, imageDim, outImage);

		float fmax = 1.0f;
		int totalpix = imageDim[0] * imageDim[1] * imageDim[2] * mag * mag *mag;
		for (int i = 0; i < totalpix; i++) {
			if (isnan(outImage[i]))
				outImage[i] = 0;
			if (outImage[i] < 0)
				outImage[i] = 0;
			if (outImage[i] > fmax)
				fmax = outImage[i];
		}


		FILE * fp = fopen("SRRF.raw", "wb");
		fwrite(outImage, sizeof(float), imageDim[0] * imageDim[1] * imageDim[2] * mag * mag *mag, fp);
		fclose(fp);

		free(inImage);

		if (dosrrfonly == 1) {
			free(outImage);
			return p;
		}
	}

	if (splitMode == 1) {
		//Load Data
		FILE * fp = fopen("SRRF.raw", "rb");
		fread(outImage, sizeof(float), imageDim[0] * imageDim[1] * imageDim[2] * mag * mag *mag, fp);
		fclose(fp);
	}



	//Build new pixels
	vector<Ref3D> pixelsH;
	pixelsH.reserve(pixels.size() * mag * mag *mag);
	for (int i = 0; i < pixels.size(); i++) {
		int x = pixels[i].x;
		int y = pixels[i].y;
		int z = pixels[i].z;
		for (int j1 = 0; j1 < mag; j1++) {
			for (int j2 = 0; j2 < mag; j2++) {
				for (int j3 = 0; j3 < mag; j3++) {
					Ref3D ref = Ref3D(x*mag + j1, y*mag + j2, z*mag + j3);
					pixelsH.push_back(ref);
				}
			}
		}
	}
	vector<float> intensitiesH(pixelsH.size(), 0);
	for (unsigned int i = 0; i < pixelsH.size(); i++) {
		int x = pixelsH[i].x;
		int y = pixelsH[i].y;
		int z = pixelsH[i].z;
		int idx = z*imageDim[0] * imageDim[1] * mag * mag + y*imageDim[0] * mag + x;
		intensitiesH[i] = outImage[idx];
	}
	free(outImage);

	double max_intensity = *max_element(intensitiesH.begin(), intensitiesH.end());
	printf("max Int: %lf\n", max_intensity);

	if (max_intensity < 0)
		return p;

	MT19937& rng = *(p.rng);

	while ((int)p.spots.size() < num_spots)
	{
		int element = floor(rng() * pixelsH.size());
		double y = rng() * max_intensity;

		if (y <= intensitiesH[element]) {
			float x = pixelsH[element].x + rng() - .5;
			float y = pixelsH[element].y + rng() - .5;
			float z = pixelsH[element].z + rng() - .5;
			x = x / mag;
			y = y / mag;
			z = z / mag;
			p.spots.push_back(makeVector(intensity_mode, blur_mode, x, y, z));
		}
	}
	return p;
}
