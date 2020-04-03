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

#ifndef STORM_INCLUDE_UTILITY_H
#define STORM_INCLUDE_UTILITY_H
#include "cvd/image.h"
#include <vector>
#include <string>
#include <cstring>
#include <cerrno>
#include <cstdlib>
#include <utility>

/**computes the sign of x
@param x \e x
@return \f$ \begin{cases}
               1  & x \ge 0\\
			   -1 & x < 0
             \end{cases}\f$
@ingroup gUtility
*/
inline double sign(double x)
{
	return x>=0?1:-1;
}

/**The ubiquitous square function
@param f Number to square
@return square of the number
@ingroup gUtility
*/
inline float sq(float f) { return f*f; }

/**
	@overload
*/
inline double sq(double f) { return f*f; }


/**Cut sub images out of every member of a vector of images.
@param im Images to be cut
@param pos Top left corner
@param size Size of the patch
@returns subimages.
@ingroup gUtility
*/
const std::vector<CVD::SubImage<float> > sub_images(const std::vector<CVD::Image<float> >& im, CVD::ImageRef pos, CVD::ImageRef size);

/** Deffinition of a pixel aligned bounding box
@ingroup gUtility
*/
typedef std::pair<CVD::ImageRef, CVD::ImageRef> BBox;

/** Compute the bounding box of a set of points
@param all_spots List of points
@ingroup gUtility
*/
std::pair<CVD::ImageRef, CVD::ImageRef> boundingbox(const std::vector<CVD::ImageRef> & all_spots);


/**
@param save_spots Stream
@param save_spots_file File to open
@ingroup gUtility
*/
template<class Stream>
void open_or_die(Stream& save_spots, const std::string& save_spots_file)
{
	using std::cerr;
	using std::endl;
	using std::strerror;
	using std::exit;
	save_spots.open(save_spots_file.c_str());
	int err = errno;

	if(!save_spots.good())
	{
		cerr << "***********************************************************\n";
		cerr << "ERROR: failed to open " << save_spots_file << ": " <<strerror(err) << endl;
		cerr << "***********************************************************\n";
		exit(1);
	}
}

//RFBA: Get number of processors, Added@20180306
int GetNoOfProcessors();

#endif
