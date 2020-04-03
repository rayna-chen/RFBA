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

#include "gvars3/instances.h"
#include "cvd/image_io.h"
#include "cvd/convolution.h"
#include "TooN/wls.h"
#include <tuple>

#include "image3d.h"
#include <tiffio.h>
#include "storm_imagery.h"
#include "debug.h"
#include "utility.h"


using namespace TooN;
using namespace GVars3;
using namespace std;
using namespace std::tr1;


//RFBA: Removed@20180306
/**Load all images from disk and do the initial preprocessing.

@param names List of filenames to load.
@return preprocessed images.
@ingroup gStormImages
**/
//vector<CVD::Image<float> > load_and_preprocess_images2(const vector<string>& names)
//{
//	vector<CVD::Image<float> > ims;
//	//Load images
//	for (unsigned int i = 0; i < names.size(); i++)
//	{
//		CVD::Image<float> im = CVD::img_load(names[i]);
//		ims.push_back(im);
//
//		if (ims.back().size() != ims[0].size())
//		{
//			cerr << "Error with image " << names[i] << ":  all images must be the same size!\n";
//			exit(1);
//		}
//	}
//	double mean, variance;
//	tie(mean, variance) = mean_and_variance(ims);
//	{
//		for (unsigned int i = 0; i < ims.size(); i++)
//			transform(ims[i].begin(), ims[i].end(), ims[i].begin(), bind2nd(minus<double>(), mean));
//		for (unsigned int i = 0; i < ims.size(); i++)
//			transform(ims[i].begin(), ims[i].end(), ims[i].begin(), bind1st(multiplies<double>(), 1 / sqrt(variance)));
//	}
//
//	tie(mean, variance) = mean_and_variance(ims);
//
//	cerr << "Rescaled:\n";
//	cerr << "mean = " << mean << endl;
//	cerr << "std  = " << sqrt(variance) << endl;
//
//
//
//	//Normalize...
//
//	//Fit the background model
//	CVD::ImageRef size = ims[0].size();
//	Vector<10> p = Zeros;
//	p[6] = -3;
//	p[9] = -4;
//
//	CVD::Image<Vector<6> > monomials(size);
//	CVD::Image<double> polynomial(size);
//	for (int yy = 0; yy < size.y; yy++)
//		for (int xx = 0; xx < size.x; xx++)
//		{
//			double x = xx *2. / size.x - 1;
//			double x2 = x*x;
//			double y = yy *2. / size.y - 1;
//			double y2 = yy;
//			monomials[yy][xx] = makeVector(1, x, y, x2, x*y, y2);
//		}
//
//
//	for (int i = 0; i < 100; i++)
//	{
//		for (int yy = 0; yy < size.y; yy++)
//			for (int xx = 0; xx < size.x; xx++)
//				polynomial[yy][xx] = monomials[yy][xx] * p.slice<0, 6>();
//
//		WLS<10> wls;
//		for (unsigned int i = 0; i < ims.size(); i++)
//			for (int yy = 0; yy < size.y; yy++)
//				for (int xx = 0; xx < size.x; xx++)
//				{
//					double t = i *1. / ims.size();
//					double func = polynomial[yy][xx] * (exp(p[6] * t) + p[8] * exp(p[9] * t)) + p[7];
//
//					Vector<10> mJ;
//
//					mJ.slice<0, 6>() = exp(p[6] * t)* monomials[yy][xx];
//					//mJ.slice<3,3>() = Zeros;
//					mJ[6] = polynomial[yy][xx] * exp(p[6] * t) * t;
//					//mJ[6] = func  * t;
//					mJ[7] = 1;
//
//					mJ[8] = polynomial[yy][xx] * exp(p[9] * t);
//					mJ[9] = polynomial[yy][xx] * exp(p[9] * t) * t * p[8];
//
//					double err = ims[i][yy][xx] - func;
//
//					double w;
//
//
//					if (err > 0)
//						w = .01 / (abs(err) + .01);
//					else
//						w = 1;
//
//					wls.add_mJ(func - ims[i][yy][xx], -mJ, w);
//				}
//
//		wls.add_prior(10);
//		wls.compute();
//
//		p += wls.get_mu();
//
//		cout << p << endl << endl;
//	}
//
//	for (unsigned int i = 0; i < ims.size(); i++)
//		for (int yy = 0; yy < size.y; yy++)
//			for (int xx = 0; xx < size.x; xx++)
//			{
//				double x = xx *2. / size.x - 1;
//				double x2 = x*x;
//				double y = yy *2. / size.y - 1;
//				double y2 = yy;
//				double t = i *1. / ims.size();
//				Vector<6> f = makeVector(1, x, y, x2, x*y, y2);
//
//				double func = f * p.slice<0, 6>() * (exp(p[6] * t) + p[8] * exp(p[9] * t)) + p[7];
//				ims[i][yy][xx] -= func;
//			}
//
//	tie(mean, variance) = mean_and_variance(ims);
//
//	//A sanity check.
//	cerr << "The mean should be small compared to std:\n";
//	cerr << "mean = " << mean << endl;
//	cerr << "std  = " << sqrt(variance) << endl;
//
//	//Scale by the variance.
//	{
//		for (unsigned int i = 0; i < ims.size(); i++)
//			transform(ims[i].begin(), ims[i].end(), ims[i].begin(), bind1st(multiplies<double>(), 1 / sqrt(variance)));
//	}
//	tie(mean, variance) = mean_and_variance(ims);
//
//	cerr << "Rescaled:\n";
//	cerr << "mean = " << mean << endl;
//	cerr << "std  = " << sqrt(variance) << endl;
//
//	return ims;
//}




//RFBA: Load 3D Image, Added@20180306
Image3D<float> img_load3D(std::string name) {

	//Get Regions
	//TIFF
	CVD::Image<float> im = CVD::img_load(name);

	int dir = 0;
	int roix = im.size().x;
	int roiy = im.size().y;
	int roiz;
	short bit_sample;

	unsigned char * buffer;
	TIFF* tif = TIFFOpen(name.c_str(), "r");
	int dircount = 0;
	do {
		dircount++;
	} while (TIFFReadDirectory(tif));
	roiz = dircount;
	TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bit_sample);
	//printf("Image is %d bit\n", bit_sample);

	//Make im3D
	Image3D<float> im3D(Ref3D(roix, roiy, roiz));
	int stripSize = TIFFStripSize(tif);
	int stripMax = TIFFNumberOfStrips(tif);
	int imageOffset = 0;
	int result = 0;
	if (bit_sample == 8)
		buffer = (unsigned char *)malloc(roix * roiy * roiz * 1 * sizeof(unsigned char));
	else
		buffer = (unsigned char *)malloc(roix * roiy * roiz * 2 * sizeof(unsigned char));
	//do load
	for (dir = 0; dir < roiz; dir++) {
		TIFFSetDirectory(tif, dir);
		// Read in the possibly multiple strips
		for (int stripCount = 0; stripCount < stripMax; stripCount++) {
			if ((result = TIFFReadEncodedStrip(tif, stripCount, buffer + imageOffset, stripSize)) == -1) {
				printf("Read error on index %d input strip number %d\n", dir, stripCount);
			}
			imageOffset += result;
		}
	}
	TIFFClose(tif);

	if (bit_sample == 8) {
		for (int k = 0; k < roiz; k++) {
			for (int j = 0; j < roiy; j++) {
				for (int i = 0; i < roix; i++) {
					int idx = k*(roix*roiy) + j*roix + i;
					im3D.data()[idx] = buffer[idx];
				}
			}
		}
	}
	else {
		unsigned short * img = (unsigned short *)buffer;
		for (int k = 0; k < roiz; k++) {
			for (int j = 0; j < roiy; j++) {
				for (int i = 0; i < roix; i++) {
					int idx = k*(roix*roiy) + j*roix + i;
					im3D.data()[idx] = img[idx];
				}
			}
		}
	}
	free(buffer);
	return(im3D);

}


//RFBA: Update for 3D Image, Modified@20180306
/**Load all images from disk and do the initial preprocessing. Currently
this is a high pass filter to make the resultimg images zero mean.

The filter is controlled with the \c preprocess.lpf and \c preprocess.skip Gvars

See also load_and_preprocess_image()

@param names List of filenames to load.
@return preprocessed images.
@ingroup gStormImages
**/
vector<Image3D<float> > load_and_preprocess_images(const vector<string>& names)
{
	vector<Image3D<float> > ims;

	//float wide = GV3::get<float>("preprocess.lpf", 0., -1);
	//bool p = GV3::get<bool>("preprocess.skip", 0, -1);

	for (unsigned int i = 0; i < names.size(); i++)
	{
		Image3D<float> im = img_load3D(names[i]);

		ims.push_back(preprocess_image(im));

		if (ims.back().size() != ims[0].size())
		{
			cerr << "Error with image " << names[i] << ":  all images must be the same size!\n";
			exit(1);
		}
	}
	return ims;
}

//RFBA: Update for 3D Image, Modified@20180306
/**Compute the mean and variance of the (on average) darkest pixels, in order
to find the correct scaling, by examining hte background.
*/
pair<double, double> auto_fixed_scaling(const vector<Image3D<float> >& ims, double frac)
{
	assert_same_size(ims);

	//Compute the mean image (ish)
	Image3D<double> ave(ims[0].size());
	//RFBA: Update for 3D Image, Modified@20180306
	for (unsigned int i = 0; i < ims.size(); i++)
		for (int z = 0; z < ave.size().z; z++)
			for (int y = 0; y < ave.size().y; y++)
				for (int x = 0; x < ave.size().x; x++) {
					int idx = z*(ave.size().y*ave.size().x) + y*ave.size().x + x;
					ave.data()[idx] += ims[i].data()[idx];
				}


	//Find the smallest N% of the pixels...
	vector<pair<double, Ref3D> > pixels;
	for (int z = 0; z < ave.size().z; z++)
		for (int y = 0; y < ave.size().y; y++)
			for (int x = 0; x < ave.size().x; x++) {
				int idx = z*(ave.size().y*ave.size().x) + y*ave.size().x + x;
				pixels.push_back(make_pair(ave[idx], Ref3D(x, y, z)));
			}

	int npix = (int)floor(frac *pixels.size() + 0.5);
	npix = max(0, min(npix, (int)pixels.size()));

	nth_element(pixels.begin(), pixels.begin() + npix, pixels.end());

	pixels.resize(npix);

	//Now compute the mean and variance of those pixels.
	double sum = 0, sum2 = 0;

	for (unsigned int i = 0; i < ims.size(); i++)
	{
		for (unsigned int j = 0; j < pixels.size(); j++)
		{
			sum += ims[i][pixels[j].second];
			sum2 += sq(ims[i][pixels[j].second]);
		}
	}

	double num = 1.0 * pixels.size() * ims.size();
	double mean = sum / num;
	double std = sqrt(((sum2 / num) - mean*mean) * num / (num - 1));

	cout << "Automatic determination of fixed scaling:" << endl
		<< "mean       = " << mean << endl
		<< "std        = " << std << endl
		<< "sqrt(mean) = " << sqrt(mean * 255) / 255 << endl;

	return make_pair(mean, std);
}

//RFBA: Update for 3D Image, Modified@20180306
/**Wrapper for load_and_preprocess_images() to allow more flexible behaviour.

@param files List of filenames to load.
@return preprocessed images.
@ingroup gStormImages
**/
vector<Image3D<float> > load_and_normalize_images(const vector<string>& files)
{
	//Load the raw data, and then load the spot parameters.
	double mean, variance;

	vector<Image3D<float> > ims = load_and_preprocess_images(files);

	tie(mean, variance) = mean_and_variance(ims);

	if (GV3::get<bool>("preprocess.fixed_scaling", 0, FATAL_IF_NOT_DEFINED))
	{
		bool skip = GV3::get<bool>("preprocess.skip");
		if (!skip)
		{
			cerr << "WARNING WARNING WARNING WARNING!!!!!!!!!!!!!!!\n";
			cerr << "preprocessing and fixed scaling selected!!!\n";
			exit(1);
		}

		double sub, div;
		if (GV3::get<bool>("preprocess.fixed_scaling.auto", 0, FATAL_IF_NOT_DEFINED))
		{
			double frac = GV3::get<double>("preprocess.fixed_scaling.auto.proportion", 0, FATAL_IF_NOT_DEFINED);
			tie(sub, div) = auto_fixed_scaling(ims, frac);
		}
		else
		{
			sub = GV3::get<double>("preprocess.fixed_scaling.subtract", 0, FATAL_IF_NOT_DEFINED);
			div = GV3::get<double>("preprocess.fixed_scaling.divide", 0, FATAL_IF_NOT_DEFINED);
		}

		for (unsigned int i = 0; i < ims.size(); i++)
			for (int j = 0; j < ims[i].size().area(); j++)
				ims[i].data()[j] = (ims[i].data()[j] - sub) / div;
	}
	else
	{
		//A sanity check.
		cerr << "The mean should be small compared to std:\n";
		cerr << "mean = " << mean << endl;
		cerr << "std  = " << sqrt(variance) << endl;

		//Scale by the variance.
		{
			for (unsigned int i = 0; i < ims.size(); i++)
				for (int j = 0; j < ims[i].size().area(); j++)
					ims[i].data()[j] = ims[i].data()[j] / sqrt(variance);
		}
	}

	tie(mean, variance) = mean_and_variance(ims);

	//A sanity check.
	cerr << "Rescaled:\n";
	cerr << "mean = " << mean << endl;
	cerr << "std  = " << sqrt(variance) << endl;

	return ims;
}



