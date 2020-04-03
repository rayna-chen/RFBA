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

#ifndef STORM_INCLUDE_DEBUG_H
#define STORM_INCLUDE_DEBUG_H

#include <cassert>

/** Determines that all images in the incoming container are the 
same size, and that the container is not empty
@param images Container to check
@ingroup gDebug
*/
template<class C> void assert_same_size(const C& images)
{
	assert(!images.empty());
	for(typename C::const_iterator i=images.begin(); i != images.end(); i++)
		assert(i->size() == images.front().size());
}


#endif

