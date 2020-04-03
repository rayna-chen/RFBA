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

#ifndef STORM_INCLUDE_STORM_IMAGERY_H
#define STORM_INCLUDE_STORM_IMAGERY_H
#include <cvd/image.h>
#include <cvd/byte.h>
#include <utility>
#include <vector>
#include "image3d.h"


//! @cond Doxygen_Suppress

//RFBA: Removed@20180306
//std::vector<CVD::Image<float> > load_and_preprocess_images2(const std::vector<std::string>& names);


std::vector<Image3D<float> > load_and_preprocess_images(const std::vector<std::string>& names);
Image3D<float> preprocess_image(const Image3D<float>&);
std::pair<float, float> mean_and_variance(const std::vector<Image3D<float> >& images);
std::vector<Image3D<float> > load_and_normalize_images(const std::vector<std::string>& files);
//! @endcond
#endif



