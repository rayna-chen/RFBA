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

#include "tag/printf.h"
#undef make_tuple

#include <tuple>
#include <algorithm>
#include <climits>
#include <iomanip>
#include <map>
#include <cvd/image_io.h>
#include <cvd/image_convert.h>
#include <cvd/morphology.h>
#include <cvd/connected_components.h>
#include <cvd/draw.h>
#include <cvd/vector_image_ref.h>
#include <iostream>
#include <io.h>
#include <string>
#include <gvars3/instances.h>

#include "storm_imagery.h"
#include "multispot5.h"
#include "multispot5_place_choice.h"
#include "utility.h"
#include "tiffio.h"


using namespace std;
using namespace std::tr1;
//using namespace CVD;
using namespace GVars3;
using namespace TooN;


//RFBA: Helper Function, Added@20180306
void printDir(const char* path, char ** start, int* num)
{
	struct _finddata_t data;
	*num = 0;
	//Get path:
	char * newpath = (char *)malloc(255 * sizeof(char));
	int i = 0;
	int lastchar = 0;
	for (i = strlen(path) - 1; i >= 0; i--) {
		if (path[i] == '\\') {
			lastchar = i;
			break;
		}
	}
	for (i = 00; i <= lastchar; i++) {
		newpath[i] = path[i];
	}
	newpath[i] = '\0';

	intptr_t hnd = _findfirst(path, &data);   
	if (hnd < 0)
	{
		perror(path);
	}
	int  nRet = (hnd < 0) ? -1 : 1;
	while (nRet >= 0)
	{
		if (data.attrib != _A_SUBDIR) {
			*(start + *num) = (char *)malloc(256 * sizeof(char));
			sprintf(*(start + *num), "%s\\%s", newpath, data.name);
			//printf("%s\n", *(start + *num));
			(*num)++;
		}
		nRet = _findnext(hnd, &data);
	}
	_findclose(hnd);    
}


vector<vector<CVD::ImageRef> > get_regions(const CVD::SubImage<double>& log_ratios)
{
	gvar3<double> radius("radius", 0, 1);

	//Set the liklihood ratio threshold/spot density prior
	//same thing.
	double threshold = GV3::get<double>("threshold", 0, -1);
	int edge = GV3::get<int>("edge", 0, -1);


	//Threshold image
	CVD::Image<CVD::byte> thresholded(log_ratios.size(), 0);
	for (int r = 0; r < thresholded.size().y; r++)
		for (int c = 0; c < min(thresholded.size().x, edge); c++)
			thresholded[r][c] = 255 * (log_ratios[r][c] > threshold);

	//Dilate
	CVD::Image<CVD::byte> dilated = morphology(thresholded, CVD::getDisc(*radius), CVD::Morphology::BinaryDilate<CVD::byte>());

	transform(dilated.begin(), dilated.end(), dilated.begin(), bind1st(multiplies<int>(), 255));

	//Connected components of dilated image
	vector<CVD::ImageRef> fg;
	for (int r = 0; r < thresholded.size().y; r++)
		for (int c = 0; c < min(thresholded.size().x, edge); c++)
			if (dilated[r][c])
				fg.push_back(CVD::ImageRef(c, r));

	vector<vector<CVD::ImageRef> > regions;
	connected_components(fg, regions);

	return regions;
}


//RFBA: Default Values, Added@20180306
int MAXTHREAD = 8;
int dosrrfonly = 0;
int startPoints = 15;
int maxIter = 1000;
int splitMode = 0;

int addGuide = 0;
double zFactor = 1;


void mmain(int argc, char** argv)
{

	//detect CPU num
	MAXTHREAD = GetNoOfProcessors();
	printf("Version: %d, CPU Number is %d\n", SWVER, MAXTHREAD);


	GUI.LoadFile("multispot5.cfg");
	int lastarg = GUI.parseArguments(argc, argv);
	if (lastarg >= argc)
	{
		cerr << "Specify the images to load\n";
		exit(1);
	}

	int filenum = 0;
	char ** filelist = (char **)malloc(8192 * sizeof(char *));
	printDir(*(argv + lastarg), filelist, &filenum);

	vector<string> files(filelist, filelist + filenum);


	//Save this now since the de-checkpointing code will kl0bber it 
	//when it reloads the gvars



	string save_spots_file = GV3::get<string>("save_spots", "", -1);
	string checkpoint_file = GV3::get<string>("load_checkpoint", "", 1);

	//RFBA: Parametes, Added@20180306
	string srrf_only = GV3::get<string>("srrf_only", "", 1);
	if (srrf_only == "1")
		dosrrfonly = 1;
	string split_mode = GV3::get<string>("split_mode", "", 1);
	if (split_mode == "1")
		splitMode = 1;
	string add_guide = GV3::get<string>("add_guide", "", 1);
	if (add_guide == "1")
		addGuide = 1;
	string z_factor = GV3::get<string>("z_factor", "1", 1);
	zFactor = atof(z_factor.c_str());
	printf("zFactor: %f\n", zFactor);



	//Load Points override and Max Iter
	string start_points = GV3::get<string>("start_points", "", 1);
	if (start_points != "")
		startPoints = atoi(start_points.c_str());

	string max_iter = GV3::get<string>("max_iter", "", 1);
	if (max_iter != "")
		maxIter = atoi(max_iter.c_str());


	if (checkpoint_file != "")
	{
		//Load and de-checkpointing
		ifstream chk;
		open_or_die(chk, checkpoint_file);

		StateParameters p;

		try {
			p = parse_log_file(chk);
		}
		catch (LogFileParseError e)
		{
			cerr << "SI TEH FUX0R11ONEone!oneleven: " << e.what << endl;
			exit(1);
		}

		vector<Image3D<float> > ims = load_and_normalize_images(files);

		//Restore kl0bbered variable
		GV3::get<string>("save_spots") = save_spots_file;

		ofstream save_spots;
		open_or_die(save_spots, save_spots_file);

		fit_spots_new(ims, p, save_spots, *null_graphics());

	}


	vector<Image3D<float> > ims = load_and_normalize_images(files);


	//Load the log_ratios image.
	//We will use this as a starting point for searching for spots.
	CVD::Image<double> log_ratios;
	try
	{
		log_ratios = CVD::img_load(GV3::get<string>("log_ratios", "", -1));
	}
	catch (CVD::Exceptions::All e)
	{
		cerr << "Error loading " << GV3::get<string>("log_ratios", "") << ": " << e.what << endl;
		exit(1);
	}



	gvar3<int> cluster_to_show("cluster_to_show", 0, -1);
	gvar3<int> use_largest("use_largest", 0, 1);

	vector<vector<CVD::ImageRef> > regions;

	regions = get_regions(log_ratios);


	if (*use_largest && !regions.empty())
	{
		*cluster_to_show = 0;
		for (unsigned int i = 1; i < regions.size(); i++)
			if (regions[i].size() > regions[*cluster_to_show].size())
				*cluster_to_show = i;

	}
	else
		*cluster_to_show = max(min(*cluster_to_show, (int)regions.size() - 1), 0);
	auto_ptr<FitSpotsGraphics> gr = null_graphics();

	//RFBA: Get Regions for 3D Image, Added@20180306
	int dir = 0;
	int roix = log_ratios.size().x;
	int roiy = log_ratios.size().y;
	int roiz;


	char * buffer;
	string str = GV3::get<string>("log_ratios", "", -1);
	TIFF* tif = TIFFOpen(str.c_str(), "r");

	int dircount = 0;
	do {
		dircount++;
	} while (TIFFReadDirectory(tif));
	roiz = dircount;
	//Make log_ratios3D
	Image3D<double> log_ratios3D(Ref3D(roix, roiy, roiz));

	printf("%d directories in %s\n", dircount, str.c_str());
	int stripSize = TIFFStripSize(tif);
	int stripMax = TIFFNumberOfStrips(tif);
	int imageOffset = 0;
	int result = 0;
	buffer = (char *)malloc(roix * roiy * roiz * 3 * sizeof(char));
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

	vector<Ref3D>  regions3D;
	for (int k = 0; k < roiz; k++) {
		for (int j = 0; j < roiy; j++) {
			for (int i = 0; i < roix; i++) {
				int idx = k*(roix*roiy) + j*roix + i;
				if ((buffer[idx * 3] != 0) || (buffer[idx * 3 + 1] != 0) || (buffer[idx * 3 + 2] != 0)) {
					Ref3D ref = Ref3D(i, j, k);
					log_ratios3D[ref] = 0.889f;
					regions3D.push_back(ref);
				}
			}
		}
	}
	free(buffer);

	place_and_fit_spots(ims, regions3D, log_ratios3D, save_spots_file, *gr);
}


int main(int argc, char** argv)
{
	try {
		mmain(argc, argv);
	}
	catch (CVD::Exceptions::All e)
	{
		cerr << "Fatal error: " << e.what << endl;
	}
}
