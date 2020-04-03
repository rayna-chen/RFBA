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

#ifndef MULTISPOT5_PLACE_CHOICE_H
#define MULTISPOT5_PLACE_CHOICE_H
#include "multispot5.h"

//RFBA: Update Parameters for 3D Image Computing, Modified@20180306
void place_and_fit_spots(const std::vector<Image3D<float> >& ims, const std::vector<Ref3D>& region, const Image3D<double>& log_ratios, std::string save_spots_file, FitSpotsGraphics& g, const std::string&s = "");
void place_and_fit_spots(const std::vector<Image3D<float> >& ims, const std::vector<Ref3D>& region, const Image3D<double>& log_ratios, std::ofstream& save_spots_file, FitSpotsGraphics& g, UserInterfaceCallback&, const std::string&s = "");


#endif
