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

#include "multispot5_place_methods.h"
#include "multispot5_place_choice.h"
#include "debug.h"
#include <gvars3/instances.h>
#include <algorithm>

using namespace std;
using namespace CVD;
using namespace GVars3;
using namespace TooN;


//RFBA: Flags for SRRF only Mode, Manual StartPoints Mode and Splitted Processing Mode, Added@20180306
extern int dosrrfonly;
extern int startPoints;
extern int splitMode;


void place_and_fit_spots(const vector<Image3D<float>>& ims, const vector<Ref3D>& region_, const Image3D<double>& log_ratios, string save_spots_file, FitSpotsGraphics& g, const string& extra)
{
	auto_ptr<UserInterfaceCallback> ui(null_ui());
	ofstream save_spots;
	open_or_die(save_spots, save_spots_file);
	place_and_fit_spots(ims, region_, log_ratios, save_spots, g, *ui, extra);
}


void place_and_fit_spots(const vector<Image3D<float>>& ims, const vector<Ref3D>& region_, const Image3D<double>& log_ratios, ofstream& save_spots, FitSpotsGraphics& g, UserInterfaceCallback& ui, const std::string& extra)
{
	assert_same_size(ims);
	//assert(ims[0].size() == log_ratios.size());

	vector<Ref3D> region = region_;
	sort(region.begin(), region.end());

	string mode = GV3::get<string>("mode", "new", 1);

	if (mode == "new")
	{
		string placement = GV3::get<string>("placement", "uniform", 1);
		save_spots << extra << endl;

		if (placement == "ye_olde")
		{
			//RFBA: Remove ye_olde Mode Support, Removed@20180306
			//StateParameters p(generate_state_parameters_ye_olde(log_ratios, ims, region));
			//fit_spots_new(ims, p, save_spots, g, ui);
			printf(" ye_olde NOT supported!\n");
		}
		else if (placement == "uniform")
		{
			//RFBA: Remove uniform Mode Support, Removed@20180306
			//int num = GV3::get<int>("placement.uniform.num_spots", 0, -1);
			//StateParameters p(place_spots_uniform(num, region, log_ratios.size()));
			//fit_spots_new(ims, p, save_spots, g, ui);
			printf("uniform NOT supported!\n");
		}
		else if (placement == "intensity_sampled")
		{
			int num = GV3::get<int>("placement.uniform.num_spots", 0, -1);
			if (splitMode == 1)
				num = startPoints;
			StateParameters p(place_spots_intensity_sampled(num, region, ims, dosrrfonly));
			if (dosrrfonly != 1)
				fit_spots_new(ims, p, save_spots, g, ui);
		}
		else
			cerr << "Mode must be intensity_sampled not `" + placement + "'.\n";
	}
	else
	{
		cerr << "Mode must be new , not `" + mode + "'.\n";
	}

}

