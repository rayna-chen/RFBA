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

#ifndef INC_DRIFT_H
#define INC_DRIFT_H


///Compute the integer step number from the frame number
///
///Drift is linear, with a piecewise constant approximation
///@ingroup gMultiSpotDrift
///@param frames Number of frames
///@param steps Number of steps
///@param frame Current frame number
int frame_to_step(int frames, int steps, int frame)
{
	return (frame * steps) / frames;
}


///Compute the approximate, real-valued frame number from the step number
///
///Drift is linear, with a piecewise constant approximation
///@ingroup gMultiSpotDrift
///@param frames Number of frames
///@param steps Number of steps
///@param frame Current frame number
double step_to_approximate_frame(int frames, int steps, int step)
{
	return (step + 0.5)/steps * frames;
}



//RFBA: Updated Vector for 3D Process, Modified@20180306

///Apply the drift model to the spot
///@ingroup gMultiSpotDrift
///@param spot Spot parameters (size/brightness/position)
///@param drift Drift vector
///@param frame_number Frame number, note that it is real-valued.
Vector<5> drift_spot(const Vector<5>& spot, const Vector<3>& drift, double frame_number)
{
	return spot + frame_number * makeVector(0, 0, drift[0], drift[1], drift[2]);
}


#endif
