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

//The functiins from storm_imagery.h which don't make use of img_load
//Separated to allow JNI plugin to build more easily.
#include "debug.h"
#include "utility.h"
#include "storm_imagery.h"
#include <gvars3/instances.h>
#include <cvd/convolution.h>
#include "image3d.h"

using namespace CVD;
using namespace GVars3;
using namespace std;

/**Do the initial preprocessing on a loaded image. Currently
this is a high pass filter to make the resultimg images zero mean.

The filter is controlled with the \c preprocess.lpf and \c preprocess.skip Gvars

See also load_and_preprocess_images()

@param im  images
@return preprocessed images.
@ingroup gStormImages
*/
Image3D<float> preprocess_image(const Image3D<float>& im)
{
	
	float wide = GV3::get<float>("preprocess.lpf", 0., -1);
	bool p = GV3::get<bool>("preprocess.skip", 0, -1);

	//Highpass filter the images using blur and subtract
	if (!p)
	{
		//RFBA: Update for 3D Image, Modified@20180306
		Image3D<float> f3D = im.copy_from_me();
		for (int i = 0; i < im.size().z; i++) {
			ImageRef imgsz(im.size().x, im.size().y);
			Image<float> f(imgsz, 0), fwide(imgsz, 0), f2D(imgsz, 0);
			//Copy image2D
			for (int j = 0; j < imgsz.area(); j++) 
				f2D.data()[j] = im.data()[imgsz.area() * i + j];
			
			convolveGaussian_fir(f2D, fwide, wide);
			for (int r = 1; r < im.size().y - 1; r++)
				for (int c = 1; c < im.size().x - 1; c++)
					f[r][c] = f2D[r][c] - fwide[r][c];

			//Write back
			for (int j = 0; j < imgsz.area(); j++) 
			 f3D.data()[imgsz.area() * i + j] = f.data()[j];
		}	
		return f3D;
	}
	else {

		Image3D<float> haha = im.copy_from_me();
		return im.copy_from_me();
	}

}

//RFBA: Update for 3D Image, Modified@20180306
/**Find the mean and variance of a stack of images
@param images Image stack
@return (mean, variance)
@ingroup gStormImages
*/
pair<float, float> mean_and_variance(const vector<Image3D<float> >& images)
{
	assert_same_size(images);

	double sum = 0, sum2 = 0, area = 0;

	for (unsigned int i = 0; i < images.size(); i++)
	{
		area += images[i].size().area();
		for (int t = 0; t < images[i].size().z; t++)
			for (int r = 0; r < images[i].size().y; r++)
				for (int c = 0; c < images[i].size().x; c++)
				{
					int index = t*images[i].size().y*images[i].size().x + r*images[i].size().x + c;
					sum += images[i][index];
					sum2 += sq(images[i][index]);
				}
	}

	sum /= area;
	sum2 /= area;
	return make_pair(sum, sum2 - sq(sum));
}


