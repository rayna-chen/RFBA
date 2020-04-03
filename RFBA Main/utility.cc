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

#include "utility.h"
#include "debug.h"
#include <climits>
#include <algorithm>
#include <windows.h>

using namespace std;
using namespace CVD;

//! @cond Doxygen_Suppress
const std::vector<CVD::SubImage<float> > sub_images(const std::vector<CVD::Image<float> >& im, CVD::ImageRef pos, ImageRef size)
{
	assert_same_size(im);

	vector<SubImage<float> > subs;

	for(unsigned int i=0; i < im.size(); i++)
		subs.push_back(im[i].sub_image(pos, size));
	return subs;
}

pair<ImageRef, ImageRef> boundingbox(const vector<ImageRef> & all_spots)
{
	ImageRef lo(INT_MAX, INT_MAX), hi(INT_MIN, INT_MIN);
	for(unsigned int i=0; i < all_spots.size(); i++)
	{
		lo[0] = min(lo[0], all_spots[i][0]);
		lo[1] = min(lo[1], all_spots[i][1]);

		hi[0] = max(hi[0], all_spots[i][0]);
		hi[1] = max(hi[1], all_spots[i][1]);
	}

	return make_pair(lo, hi - lo + ImageRef(1,1));
}




//RFBA: Get number of processors, Added@20180306
int GetNoOfProcessors()
{
	SYSTEM_INFO si;
	GetSystemInfo(&si);
	return si.dwNumberOfProcessors;
}
//! @endcond 
