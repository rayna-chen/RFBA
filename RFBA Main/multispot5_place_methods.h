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

#ifndef MULTISPOT5_PLACE_METHODS_H
#define MULTISPOT5_PLACE_METHODS_H
#include "multispot5.h"

StateParameters place_spots_uniform(int num_spots, const std::vector<CVD::ImageRef>& pixels, const CVD::ImageRef& size);

//RFBA: Update Parameters for 3D Image Computing, Modified@20180306
StateParameters place_spots_intensity_sampled(int num_spots, const std::vector<Ref3D>& pixels, const std::vector<Image3D<float> >&, int dosrrfonly);

#endif
