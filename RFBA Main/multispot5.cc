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

#include <cstdlib>
#include <cerrno>
#include <cstring>
#include <stack>
#include <algorithm>
#include <climits>
#include <iomanip>
#include <map>
#include <memory>
#include <cvd/image_io.h>
#include <cvd/image_convert.h>
#include <cvd/morphology.h>
#include <cvd/connected_components.h>
#include <cvd/draw.h>
#include <cvd/vector_image_ref.h>

#include <cvd/random.h>
#include <cvd/timer.h>
#include <gvars3/instances.h>

#include <TooN/functions/derivatives.h>
#include <TooN/determinant.h>
#include <TooN/SymEigen.h>
#include <TooN/optimization/conjugate_gradient.h>

#include "conjugate_gradient_only.h"
#include "forward_algorithm.h"
#include "numerical_derivatives.h"
#include "storm.h"
#include "storm_imagery.h"
#include "debug.h"
#include "sampled_multispot.h"
#include "mt19937.h"
#include "utility.h"
#include "multispot5.h"
#include "pthread.h"



//For benchmarking...
#define TIME(X) 
//#define TIME(X) X

using namespace std;
using namespace CVD;
using namespace GVars3;
using namespace TooN;
using namespace std::tr1;


extern int MAXTHREAD;
extern int maxIter;
extern int addGuide;
extern double zFactor;

//RFBA: Add Reference to SRRF GPU Lib, Added@20180306
#pragma comment(lib,"GradientCalc.lib")
_declspec(dllimport) int SRRFGetMag();


///Empty destructor
UserInterfaceCallback::~UserInterfaceCallback() {}

///User interface callback class which does nothing.
class NullUICallback : public UserInterfaceCallback
{
	void per_spot(int, int, int, int) {};
	void per_modification(int, int, int) {};

	void perhaps_stop() {};
};

///Factory function to generate an instance of NullGraphics
///@ingroup gStorm
auto_ptr<UserInterfaceCallback> null_ui()
{
	return auto_ptr<UserInterfaceCallback>(new NullUICallback);
}



//Declare the graphics classes.
//These provide all the debug drawing operations so that the code can run in GUI and
//headless mode easily and without macros.
FitSpotsGraphics::~FitSpotsGraphics() {}

///Graphics class which does absoloutely nothing
///@ingroup gStorm
class NullGraphics : public  FitSpotsGraphics
{
public:
	virtual void init(CVD::ImageRef) {}

	virtual void swap() {}
	virtual void draw_pixels(const std::vector<CVD::ImageRef>&, float, float, float, float) {}
	virtual void draw_bbox(const BBox&) {}
	virtual void glDrawCross(const TooN::Vector<2>&, int) {}
	virtual ~NullGraphics() {}
};

///Factory function to generate an instance of NullGraphics
///@ingroup gStorm
auto_ptr<FitSpotsGraphics> null_graphics()
{
	return auto_ptr<FitSpotsGraphics>(new NullGraphics);
}


///There are two sensible ways of storing the state vector of spot positions.
///This function converts between them. See also spots_to_vector.
///@param s list of spots to convert
///@ingroup gUtility
Vector<> spots_to_Vector(const vector<Vector<5> >& s)
{
	Vector<> r(s.size() * 5);
	for (unsigned int i = 0; i < s.size(); i++)
	{
		r[i * 5 + 0] = s[i][0];
		r[i * 5 + 1] = s[i][1];
		r[i * 5 + 2] = s[i][2];
		r[i * 5 + 3] = s[i][3];
		r[i * 5 + 4] = s[i][4];
	}
	return r;
}

///There are two sensible ways of storing the state vector of spot positions.
///This function converts between them. See also spots_to_Vector.
///@param s list of spots to convert
///@ingroup gUtility
vector<Vector<5> > spots_to_vector(const Vector<>& s)
{
	vector<Vector<5> > r(s.size() / 5);
	for (unsigned int i = 0; i < r.size(); i++)
		r[i] = s.slice<Dynamic, 5>(i * 5, 5);
	return r;
}

///Normalize an image for display purposes.
///@ingroup gUtility
Image<byte> scale_to_bytes(const Image<float>& im, float lo, float hi)
{
	Image<byte> out(im.size());
	for (int r = 0; r < out.size().y - 0; r++)
		for (int c = 0; c < out.size().x - 0; c++)
			out[r][c] = (int)floor((im[r][c] - lo) * 255 / (hi - lo));
	return out;
}

///Normalize an image for display purposes.
///@ingroup gUtility
Image<byte> scale_to_bytes(const Image<float>& im)
{
	float lo = *min_element(im.begin(), im.end());
	float hi = *max_element(im.begin(), im.end());
	Image<byte> out(im.size());
	for (int r = 0; r < out.size().y - 0; r++)
		for (int c = 0; c < out.size().x - 0; c++)
			out[r][c] = (int)floor((im[r][c] - lo) * 255 / (hi - lo));

	return out;
}

///Average the input image stack for display purposes
///@ingroup gUtility
Image3D<float> average_image(const vector<Image3D<float> >& ims)
{
	assert_same_size(ims);
	int x = ims[0].size().x;
	int y = ims[0].size().y;
	int z = ims[0].size().z;
	float * temp = (float *)malloc(x*y*z * sizeof(float));
	for (int i = 0; i < x*y*z; i++)
		temp[i] = 0;


	for (int j = 0; j < ims.size(); j++)
		for (int i = 0; i < x*y*z; i++)
			temp[i] += ims[j][i];

	for (int i = 0; i < x*y*z; i++)
		temp[i] = temp[i] / ims.size();

	Image3D<float> r(ims[0].size(), temp);
	free(temp);
	return r;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

///Closure hoding the data required do use GibbsSampler2
///See FitSpots for naming of variables.
///@ingroup gStorm
class DataForMCMC
{
protected:
	const vector<Ref3D>& pixels;
	const vector<vector<double> >& pixel_intensities;//[frame][pixel]
	const double mu_brightness, sigma_brightness, mu_blur, sigma_blur;
	const double variance;
	const int samples, sample_iterations;
	const Matrix<3> A;
	const Vector<3> pi;
	MT19937& rng;
public:

	MT19937& get_rng()const
	{
		return rng;
	}

public:
	DataForMCMC(const vector<Ref3D>& pixels_,
		const vector<vector<double> >& pixel_intensities_,
		double mu_brightness_,
		double sigma_brightness_,
		double mu_blur_,
		double sigma_blur_,
		double variance_,
		int samples_,
		int sample_iterations_,
		Matrix<3> A_,
		Vector<3> pi_,
		MT19937& rng_)
		:pixels(pixels_),
		pixel_intensities(pixel_intensities_),
		mu_brightness(mu_brightness_),
		sigma_brightness(sigma_brightness_),
		mu_blur(mu_blur_),
		sigma_blur(sigma_blur_),
		variance(variance_),
		samples(samples_),
		sample_iterations(sample_iterations_),
		A(A_),
		pi(pi_),
		rng(rng_)
	{}
};

///Class implementing the Kahan summation algorithm
///to allow accurate summation of very large numbers of
///doubles.
///@ingroup gUtility
class Kahan {
private:
	double y; ///< Input with the compensation removed. Temporary working space.
	double c; ///< Running compenstation for low-order bits.
	double t; ///< y + sum, which loses low order bits of y. Temporary working space.
public:
	double sum; ///< running sum

	Kahan()
		:c(0), sum(0)
	{}

	/// Add a number to the running sum
	/// @param i Number to add
	void add(double i)
	{
		//y = input -c
		y = i;
		y -= c;

		//t = sum + y
		t = sum;
		t += y;

		//c = (t - sum) - y
		//c = ((sum + y) - sum) - y)
		c = t;
		c -= sum;
		c -= y;
		sum = t;
	}
};


///Class for computing the negitve free energy using thermodynamic integration.
class NegativeFreeEnergy : public DataForMCMC
{
public:

	///@param d Necessary data
	NegativeFreeEnergy(const DataForMCMC& d)
		:DataForMCMC(d)
	{
	}

	///Give the noise variance given a sample number and growth parameters
	///@param sample Sample number
	///@param samples Total number of samples
	///@param base_sigma Starting value of sigme 
	///@param scale_pow Exponent scaling
	double variance_from_sample(double sample, double samples, double base_sigma, double scale_pow) const
	{
		double scale = pow(1.25, sample * 1. / samples * scale_pow);
		double sigma = base_sigma * scale;
		double new_variance = sq(sigma);

		return new_variance;
	}

	///Estimate free energy using the Slow Growth Thermodynamic Integration method given in
	///Probalistic Inference Using Markov Chain Monte Carlo Methods, Radford. M. Neal, 1993
	///Except using a 5 point stencil instead of forward differenceing in Eq 6.30
	///@param spots list of spots
	///@param spot_pixels Mask around each spot to use to save on computation
	double compute_with_mask(const Vector<>& spots, const vector<vector<int> >& spot_pixels) const
	{
		//Estimate free energy using the Slow Growth Thermodynamic Integration method given in
		//Probalistic Inference Using Markov Chain Monte Carlo Methods, Radford. M. Neal, 1993
		//Except using a 5 point stencil instead of forward differenceing in Eq 6.30
		double base_sigma = sqrt(variance); // should be 1
		double scale_pow = 100.0;

		const unsigned int nspots = spots.size() / 5;	//YZL
		const unsigned int nframes = pixel_intensities.size();
		const unsigned int npixels = pixels.size();
		assert(spots.size() % 5 == 0);
		assert(spot_pixels.size() == nspots);

		//Compute the intensities and derivatives for all spot
		vector<vector<double> > spot_intensity; //[spot][pixel]
		for (unsigned int i = 0; i < nspots; i++)
			spot_intensity.push_back(compute_spot_intensity(pixels, spots.slice<Dynamic, 5>(i * 5, 5), zFactor));	//YZL

		GibbsSampler2 sampler(pixel_intensities, spot_intensity, spots_to_vector(spots), spot_pixels, A, pi, variance, sample_iterations);

		double sum = 0;
		Kahan ksum;


		for (int sample = 0; sample < samples; sample++)
		{
			//Compute the positions of the surrounding steps
			double var1 = variance_from_sample(sample - 2, samples, base_sigma, scale_pow);
			double var2 = variance_from_sample(sample - 1, samples, base_sigma, scale_pow);
			double var3 = variance_from_sample(sample + 1, samples, base_sigma, scale_pow);
			double var4 = variance_from_sample(sample + 2, samples, base_sigma, scale_pow);



			//Take a sample
			sampler.set_variance(var2);
			sampler.next(DataForMCMC::get_rng());


			//Compute the SSD. This is fixed regardless of sigma.	
			double err_sum = 0;
			for (unsigned int frame = 0; frame < nframes; frame++)
				for (unsigned int pixel = 0; pixel < npixels; pixel++)
					err_sum -= sq(pixel_intensities[frame][pixel] - sampler.sample_intensities()[frame][pixel]);


			//Compute the derivative using a five point stencil
			//This could be done in a better but less clear way
			double e1 = err_sum / (2 * var1) - npixels*nframes*::log(2 * M_PI*var1) / 2;
			double e2 = err_sum / (2 * var2) - npixels*nframes*::log(2 * M_PI*var2) / 2;
			double e3 = err_sum / (2 * var3) - npixels*nframes*::log(2 * M_PI*var3) / 2;
			double e4 = err_sum / (2 * var4) - npixels*nframes*::log(2 * M_PI*var4) / 2;
			sum += (-e1 + 8 * e2 - 8 * e3 + e4) / 12;

			ksum.add(-e1 / 12);
			ksum.add(8 * e2 / 12);
			ksum.add(-8 * e3 / 12);
			ksum.add(e4 / 12);
		}



		double log_final = (log(variance_from_sample(samples, samples, base_sigma, scale_pow) * 2 * M_PI) / 2) * npixels * nframes;

		double priors = 0;
		for (unsigned int i = 0; i < nspots; i++)
		{
			priors += log_log_normal(spots[i * 5 + 0], mu_brightness, sigma_brightness);
			priors += log_log_normal(spots[i * 5 + 1], mu_blur, sigma_blur);
		}

		/*cout << "Thermo:\n";
		cout << "sum: " << sum -log_final << endl;
		cout << "kah: " << ksum.sum -log_final << endl;
		cout << "priors: " << priors + sum -log_final << endl;
		cout << "   kah: " << priors + ksum.sum -log_final << endl;
		*/
		//cout << log_final << endl;
		//cout << sum + log_final << endl;


		sampler.set_variance(variance);
		return -(sum + priors - log_final);
	}

	///Estimate free energy using the Slow Growth Thermodynamic Integration method given in
	///Probalistic Inference Using Markov Chain Monte Carlo Methods, Radford. M. Neal, 1993
	///Except using a 5 point stencil instead of forward differenceing in Eq 6.30
	///@param spots list of spots
	double operator()(const Vector<>& spots) const
	{
		double base_sigma = sqrt(variance); // should be 1
		double scale_pow = 100.0;

		const unsigned int nspots = spots.size() / 5;
		const unsigned int nframes = pixel_intensities.size();
		const unsigned int npixels = pixels.size();
		assert(spots.size() % 5 == 0);

		//Compute the intensities and derivatives for all spot
		vector<vector<double> > spot_intensity; //[spot][pixel]
		for (unsigned int i = 0; i < nspots; i++)
			spot_intensity.push_back(compute_spot_intensity(pixels, spots.slice<Dynamic, 5>(i * 5, 5), zFactor));

		GibbsSampler sampler(pixel_intensities, spot_intensity, spots_to_vector(spots), A, pi, variance, sample_iterations);

		double sum = 0;
		Kahan ksum;
		for (int sample = 0; sample < samples; sample++)
		{
			//Compute the positions of the surrounding steps
			double var1 = variance_from_sample(sample - 2, samples, base_sigma, scale_pow);
			double var2 = variance_from_sample(sample - 1, samples, base_sigma, scale_pow);
			double var3 = variance_from_sample(sample + 1, samples, base_sigma, scale_pow);
			double var4 = variance_from_sample(sample + 2, samples, base_sigma, scale_pow);

			//Take a sample
			sampler.set_variance(var2);
			sampler.next(DataForMCMC::get_rng());

			//Compute the SSD. This is fixed regardless of sigma.	
			double err_sum = 0;
			for (unsigned int frame = 0; frame < nframes; frame++)
				for (unsigned int pixel = 0; pixel < npixels; pixel++)
					err_sum -= sq(pixel_intensities[frame][pixel] - sampler.sample_intensities()[frame][pixel]);

			//Compute the derivative using a five point stencil
			//This could be done in a better but less clear way
			double e1 = err_sum / (2 * var1) - npixels*nframes*::log(2 * M_PI*var1) / 2;
			double e2 = err_sum / (2 * var2) - npixels*nframes*::log(2 * M_PI*var2) / 2;
			double e3 = err_sum / (2 * var3) - npixels*nframes*::log(2 * M_PI*var3) / 2;
			double e4 = err_sum / (2 * var4) - npixels*nframes*::log(2 * M_PI*var4) / 2;
			sum += (-e1 + 8 * e2 - 8 * e3 + e4) / 12;

			ksum.add(-e1 / 12);
			ksum.add(8 * e2 / 12);
			ksum.add(-8 * e3 / 12);
			ksum.add(e4 / 12);
		}

		double log_final = (log(variance_from_sample(samples, samples, base_sigma, scale_pow) * 2 * M_PI) / 2) * npixels * nframes;

		double priors = 0;
		for (unsigned int i = 0; i < nspots; i++)
		{
			priors += log_log_normal(spots[i * 5 + 0], mu_brightness, sigma_brightness);
			priors += log_log_normal(spots[i * 5 + 1], mu_blur, sigma_blur);
		}

		/*cout << "Thermo:\n";
		cout << "sum: " << sum -log_final << endl;
		cout << "kah: " << ksum.sum -log_final << endl;
		cout << "priors: " << priors + sum -log_final << endl;
		cout << "   kah: " << priors + ksum.sum -log_final << endl;
		*/
		//cout << log_final << endl;
		//cout << sum + log_final << endl;


		sampler.set_variance(variance);
		return -(sum + priors - log_final);
	}
};

/// Class for sorting a list of indexes to an array of spots lexicographically
/// according to the 2D positions of the spots.
///
///@param Cmp comparator function to specify less or greater
///@param First most significant position index
///@ingroup gMultiSpot
template<class Cmp, int First> struct IndexLexicographicPosition {
	const vector<Vector<5> >& spots;  ///< Keep around the array of spots for later comprison

	
	/// @param s Vector to sort indices of
	IndexLexicographicPosition(const vector<Vector<5> >& s)
		:spots(s)
	{}


	static const int Second = First == 2 ? 3 : 2; ///< Second most siginifcant position index for sorting

	/// Compare two indexes into the array of spots
	bool operator()(int a, int b)
	{
		Cmp cmp;

		if (cmp(spots[a][First], spots[b][First]))
			return true;
		else if (spots[a][First] == spots[b][First])
			return cmp(spots[a][Second], spots[b][Second]);
		else
			return false;
	}
};


///Closure holding image data generated using samples drawn from the model.
///NB this is used with one spot removed (i.e. set to dark). As a result, the data
///is treated as the background when considering that spot in isolation.
///See FitSpots for naming of variables.
///@ingroup gStorm
struct SampledBackgroundData
{
	const vector<vector<vector<double> > >& sample_intensities_without_spot;    //[sample][frame][pixel]
	const vector<vector<double> >& pixel_intensities;    //[frame][pixel]
	const vector<Ref3D> pixels;

	double mu_brightness, sigma_brightness, mu_blur, sigma_blur;
	const Matrix<3> A;
	const Vector<3> pi;
	double variance;

	const vector<int> O;

	SampledBackgroundData(
		const vector<vector<vector<double> > >& sample_intensities_without_spot_,
		const vector<vector<double> >& pixel_intensities_,
		const vector<Ref3D> pixels_,
		double mu_brightness_,
		double sigma_brightness_,
		double mu_blur_,
		double sigma_blur_,
		const Matrix<3> A_,
		const Vector<3> pi_,
		double variance_)
		:sample_intensities_without_spot(sample_intensities_without_spot_),
		pixel_intensities(pixel_intensities_),
		pixels(pixels_),
		mu_brightness(mu_brightness_),
		sigma_brightness(sigma_brightness_),
		mu_blur(mu_blur_),
		sigma_blur(sigma_blur_),
		A(A_),
		pi(pi_),
		variance(variance_),
		O(sequence(pixel_intensities.size()))
	{
	}
};



struct FWD_PARA {
	int idx;
	int step;
	const SampledBackgroundData * spb;
	vector<pair<double, Vector<5> > > * pspot_intensities;
	Vector<5> res;
};

//Simple Multithreading
static void  * forward_handler(void* para) {
	FWD_PARA * mypara = (FWD_PARA *)para;
	const SampledBackgroundData * spb = mypara->spb;
	int step = mypara->step;
	mypara->res = Zeros;

	for (int i = mypara->idx; i < spb->sample_intensities_without_spot.size(); i = i + step) {
		SpotWithBackground B(spb->sample_intensities_without_spot[i], *(mypara->pspot_intensities), spb->pixel_intensities, spb->variance);
		pair<double, Vector<5> > r = forward_algorithm_deriv(spb->A, spb->pi, B, spb->O);
		mypara->res += r.second;
	}
	return(0);
}


///Compute the derivative of the negative log probability with respect
///to the parameters of one spot, given some samples of the other spots.
///@ingroup gStorm
struct SpotNegProbabilityDiffWithSampledBackground : public SampledBackgroundData
{
	///@param d Necessary data for construction
	SpotNegProbabilityDiffWithSampledBackground(const SampledBackgroundData& d)
		:SampledBackgroundData(d)
	{
	}

	///Compute the probability of spot
	///@param spot Spot position
	Vector<5> operator()(const Vector<5>& spot) const
	{
		if (spot[0] <= 0 || spot[1] <= 0)
			return Ones * std::numeric_limits<double>::quiet_NaN();


		vector<pair<double, Vector<5> > > spot_intensities = compute_spot_intensity_derivatives(pixels, spot, zFactor);

		Vector<5> sum_diff_log = Zeros;

		//RFBA: Multithreading, Added@20180306
		int thread_num = (int)sample_intensities_without_spot.size();
		if (thread_num > MAXTHREAD)
			thread_num = MAXTHREAD;


		FWD_PARA * fp = (FWD_PARA *)malloc(thread_num * sizeof(FWD_PARA));
		for (int i = 0; i < thread_num; i++) {
			(fp + i)->idx = i;
			(fp + i)->pspot_intensities = &spot_intensities;
			(fp + i)->res = Zeros;
			(fp + i)->spb = this;
			(fp + i)->step = thread_num;
		}
		pthread_t *tid = (pthread_t *)malloc(thread_num * sizeof(pthread_t));

		for (int i = 0; i < thread_num; i++) {
			pthread_create((tid + i), NULL, forward_handler, (void *)(fp + i));
		}
		//Wait
		for (int i = 0; i < thread_num; i++) {
			pthread_join(*(tid + i), NULL);
		}

		//Sum
		for (int i = 0; i < thread_num; i++) {
			sum_diff_log += (fp + i)->res;
		}

		free(fp);
		free(tid);

		//cout << "22:   " << sum_diff_log << "\n" << endl;

		Vector<5> diff_log = sum_diff_log / sample_intensities_without_spot.size();

		//Compute the log probability of the prior
		Vector<5> logprior_deriv = makeVector(diff_log_log_normal(spot[0], mu_brightness, sigma_brightness),
			diff_log_log_normal(spot[1], mu_blur, sigma_blur), 0, 0, 0);


		return -(diff_log + logprior_deriv);
	}
};







///Class for computing the Hessian of the negative free energy.
///@ingroup gStorm
class FreeEnergyHessian : public DataForMCMC
{
public:

	///Constructor
	///@param d All data required
	FreeEnergyHessian(const DataForMCMC& d)
		:DataForMCMC(d)
	{
	}

	///Compute the Hessian
	///@param spots All spot positions
	///@param spot spot to compute hessian for
	Matrix<5> hessian(const vector<Vector<5> >& spots, int spot) const
	{		
		cout << "----Computing pure MCMC hessian\n";
		const unsigned int nspots = spots.size();
		const unsigned int nframes = pixel_intensities.size();
		const unsigned int npixels = pixels.size();
		cout << spot << " " << nspots << " " << nframes << " " << npixels << endl;

		vector<vector<double> > spot_intensity; //[spot][pixel]
		for (unsigned int i = 0; i < nspots; i++)
			spot_intensity.push_back(compute_spot_intensity(pixels, spots[i], zFactor));

		vector<tuple<double, Vector<5>, Matrix<5> > > spot_hess_etc = compute_spot_intensity_hessian(pixels, spots[spot], zFactor);

		GibbsSampler sampler(pixel_intensities, spot_intensity, spots, A, pi, variance, sample_iterations);

		//Compute derivative of log probability, summed (ie integrated)
		Matrix<5> sum_hess1 = Zeros(spots.size());
		Matrix<5> sum_hess2 = Zeros(spots.size());
		Vector<5> sum_deriv = Zeros(spots.size());

		for (int sample = 0; sample < samples; sample++)
		{
			sampler.next(rng);

			//Compute d log P(data | x, model) / d model, for a given sample
			//And the hessian
			Matrix<5> hess = Zeros(spots.size());
			Vector<5> deriv = Zeros(spots.size());
			for (unsigned int frame = 0; frame < nframes; frame++)
			{
				for (unsigned int pixel = 0; pixel < npixels; pixel++)
				{
					double e = pixel_intensities[frame][pixel] - sampler.sample_intensities()[frame][pixel];
					//Build up the derivative
					if (sampler.sample()[spot][frame] == 0)
					{
						hess += e * get<2>(spot_hess_etc[pixel]) - get<1>(spot_hess_etc[pixel]).as_col() * get<1>(spot_hess_etc[pixel]).as_row();
						deriv += e * get<1>(spot_hess_etc[pixel]);

					}
				}
			}

			hess[0][0] += hess_log_log_normal(spots[spot][0], mu_brightness, sigma_brightness);
			hess[1][1] += hess_log_log_normal(spots[spot][1], mu_blur, sigma_blur);
			sum_hess1 += hess;

			deriv[0] += diff_log_log_normal(spots[spot][0], mu_brightness, sigma_brightness);
			deriv[1] += diff_log_log_normal(spots[spot][1], mu_blur, sigma_blur);

			sum_hess2 += deriv.as_col() * deriv.as_row();

			sum_deriv = sum_deriv + deriv;
		}

		sum_hess1 /= (samples * variance);
		sum_hess2 /= (samples * variance);
		sum_deriv /= (samples * variance);


		cout << sum_hess1 << endl;
		cout << sum_hess2 << endl;
		cout << sum_deriv.as_col() * sum_deriv.as_row() << endl;

		cout << "......." << sum_deriv << endl;
		//Put in the prior
		//The derivative prior parts cancel out.
		//Rather sensibly this means that the second derivatives can be
		//computed without reference to the prior, and then the prior can
		//be added in later, i.e.: P(data, parameters) = P(data|parameter) P(parameters)

		//The second derivatives have been constructed to be diagonal
		DiagonalMatrix<5> hess_prior = Zeros(spots.size());

		cout << "sum of parts = \n" << sum_hess1 + sum_hess2 - sum_deriv.as_col() * sum_deriv.as_row() << endl;
		//TooN cannot currently add DiagonalMatrix to Matrix!!
		//sum_hess.diagonal_slice() += hess_prior.diagonal_slice();

		cout << "++++Done Computing pure MCMC hessian\n";
		return sum_hess1 + sum_hess2 - sum_deriv.as_col() * sum_deriv.as_row();
	}
};

///Comparator functor for the first element of a std::pair
struct LessSecond
{
	///Comparison function
	///@param a
	///@param b
	template<class A, class B> bool operator()(const pair<A, B>& a, const pair<A, B>& b) const
	{
		return a.second < b.second;
	}
};

//!@cond Doxygen_Suppress
struct NthDeriv {

	const SpotNegProbabilityDiffWithSampledBackground& compute_deriv;
	int i;

	NthDeriv(const SpotNegProbabilityDiffWithSampledBackground& c, int ii)
		:compute_deriv(c), i(ii)
	{}

	double  operator()(const Vector<5>& f) const
	{
		return compute_deriv(f)[i];
	}
};
//!@endcond




struct HESS_PARA {
	int idx;
	int step;
	int nframes;
	int npixels;
	int bs_iterations;
	const SampledBackgroundData * spb;
	vector<double>  * pspot_intensities;
	vector<tuple<double, Vector<5>, Matrix<5> > > * pspot_hess_etc;
	double * sumhess;
	double * sumderiv;
	MT19937 * rng;
};

//Simple Multithreading
static void  * hess_handler(void* para) {
	HESS_PARA * mypara = (HESS_PARA *)para;
	const SampledBackgroundData  d = *(mypara->spb);
	vector<double> spot_intensities = *(mypara->pspot_intensities);
	vector<tuple<double, Vector<5>, Matrix<5> > >  spot_hess_etc = *(mypara->pspot_hess_etc);
	int step = mypara->step;
	int nframes = mypara->nframes;
	int npixels = mypara->npixels;
	int bs_iterations = mypara->bs_iterations;

	Matrix<5> sum_hess = Zeros;
	Vector<5> sum_deriv = Zeros;


	vector<pair<Matrix<5>, Vector<5> > > hess_and_deriv_part(nframes);
	vector<State> current_sample;

	for (int i = mypara->idx; i < d.sample_intensities_without_spot.size(); i = i + step) {
		SpotWithBackground B(d.sample_intensities_without_spot[i], spot_intensities, d.pixel_intensities, d.variance);
		//Compute what the per-frame hess and deriv parts are
		//if the spot is on in a frame.

		for (int frame = 0; frame < nframes; frame++)
		{
			Matrix<5> hess = Zeros;
			Vector<5> deriv = Zeros;

			for (int pixel = 0; pixel < npixels; pixel++)
			{
				double e = d.pixel_intensities[frame][pixel] - (d.sample_intensities_without_spot[i][frame][pixel] + spot_intensities[pixel]);
				//Build up the derivative
				hess += e * get<2>(spot_hess_etc[pixel]) - get<1>(spot_hess_etc[pixel]).as_col() * get<1>(spot_hess_etc[pixel]).as_row();
				deriv += e * get<1>(spot_hess_etc[pixel]);
			}
			hess_and_deriv_part[frame] = make_pair(hess, deriv);
		}

		//Forward filtering
		std::vector<array<double, 3> > delta = forward_algorithm_delta(d.A, d.pi, B, d.O);

		for (int i = 0; i < bs_iterations; i++)
		{
			current_sample = backward_sampling<3, State>(d.A, delta, *(mypara->rng));

			Matrix<5> hess = Zeros;
			Vector<5> deriv = Zeros;
			for (int frame = 0; frame < nframes; frame++)
				if (current_sample[frame] == 0)
				{
					hess += hess_and_deriv_part[frame].first;
					deriv += hess_and_deriv_part[frame].second;
				}

			sum_hess += hess + deriv.as_col() * deriv.as_row();
			sum_deriv += deriv;
		}
	}

	for (int i = 0; i < 5; i++) {
		mypara->sumderiv[i] = sum_deriv[i];
		for (int j = 0; j < 5; j++) {
			mypara->sumhess[i * 5 + j] = sum_hess[i][j];
		}
	}


	return(0);
}



///Compute the Hessian of the log probability. The background is sampled rather sparsely, and the 
///spot in question is sampled much more densely using FFBS.
///@param spot Spot parameters
///@param d Background brightness (from other spots)
///@param bs_iterations Exter backward sampling iterations
///@param rng Random number generator
///@returns the Hessian of the log probability around the spot
Matrix<5> sampled_background_spot_hessian_ffbs(const Vector<5>& spot, const SampledBackgroundData& d, int bs_iterations, MT19937* rng)
{
	vector<tuple<double, Vector<5>, Matrix<5> > > spot_hess_etc = compute_spot_intensity_hessian(d.pixels, spot, zFactor);
	vector<double> spot_intensities = compute_spot_intensity(d.pixels, spot, zFactor);

	Matrix<5> sum_hess_log = Zeros;
	Matrix<5> sum_diff2_log = Zeros;

	vector<State> current_sample;

	const unsigned int nframes = d.pixel_intensities.size();
	const unsigned int npixels = d.pixels.size();

	Matrix<5> sum_hess = Zeros;
	Vector<5> sum_deriv = Zeros;

	//RFBA: Multithreading, Added@20180306
	int thread_num = (int)d.sample_intensities_without_spot.size();
	if (thread_num > MAXTHREAD)
		thread_num = MAXTHREAD;


	HESS_PARA * fp = (HESS_PARA *)malloc(thread_num * sizeof(HESS_PARA));
	for (int i = 0; i < thread_num; i++) {
		(fp + i)->idx = i;
		(fp + i)->step = thread_num;
		(fp + i)->nframes = nframes;
		(fp + i)->npixels = npixels;
		(fp + i)->bs_iterations = bs_iterations;
		(fp + i)->spb = &d;
		(fp + i)->pspot_intensities = &spot_intensities;
		(fp + i)->pspot_hess_etc = &spot_hess_etc;
		(fp + i)->sumhess = (double *)malloc(25 * sizeof(double));
		(fp + i)->sumderiv = (double *)malloc(5 * sizeof(double));
		(fp + i)->rng = rng + i;
	}
	pthread_t *tid = (pthread_t *)malloc(thread_num * sizeof(pthread_t));

	for (int i = 0; i < thread_num; i++) {
		pthread_create((tid + i), NULL, hess_handler, (void *)(fp + i));
	}
	//Wait
	for (int i = 0; i < thread_num; i++) {
		pthread_join(*(tid + i), NULL);
	}

	//Sum
	for (int k = 0; k < thread_num; k++) {

		for (int i = 0; i < 5; i++) {
			sum_deriv[i] += (fp + k)->sumderiv[i];
			for (int j = 0; j < 5; j++) {
				sum_hess[i][j] += (fp + k)->sumhess[i * 5 + j];
			}
		}

		free((fp + k)->sumderiv);
		free((fp + k)->sumhess);
	}

	free(fp);
	free(tid);


	//cout << "Turboderiv2:" << sum_deriv << endl;
	//cout << "Turbohess2:\n" << sum_hess << endl;

	sum_hess /= (bs_iterations * d.sample_intensities_without_spot.size() * d.variance);
	sum_deriv /= (bs_iterations * d.sample_intensities_without_spot.size() * d.variance);

	sum_hess -= sum_deriv.as_col() * sum_deriv.as_row();

	sum_hess[0][0] += hess_log_log_normal(spot[0], d.mu_brightness, d.sigma_brightness);
	sum_hess[1][1] += hess_log_log_normal(spot[1], d.mu_blur, d.sigma_blur);
	sum_deriv[0] += diff_log_log_normal(spot[0], d.mu_brightness, d.sigma_brightness);
	sum_deriv[1] += diff_log_log_normal(spot[1], d.mu_blur, d.sigma_blur);



	return sum_hess;
}



/// Convert a Vector into an image co-ordinate. Numbers are rounded
/// @param v The Vector to convert
/// @ingroup gImage
inline Ref3D ir_rounded3D(const TooN::Vector<3>& v)
{
	return Ref3D(
		static_cast<int>(v[0] > 0.0 ? v[0] + 0.5 : v[0] - 0.5),
		static_cast<int>(v[1] > 0.0 ? v[1] + 0.5 : v[1] - 0.5),
		static_cast<int>(v[2] > 0.0 ? v[2] + 0.5 : v[2] - 0.5));
}


///Which pixels belong to a given spot?
///Find the indices of those pixels
///@ingroup gStorm
void get_spot_pixels(const vector<Ref3D>& pixels, const Vector<5>& spot, vector<int>& out)
{
	//Go out to three sigma

	vector<Ref3D> pix = getBall(spot[1] * 6 + 1);
	out.resize(0);
	Ref3D offset = ir_rounded3D(spot.slice<2, 3>());		//YZL	
	for (unsigned int j = 0; j < pix.size(); j++)
	{
		int pos = lower_bound(pixels.begin(), pixels.end(), pix[j] + offset) - pixels.begin();
		if (pos != (int)pixels.size() && pixels[pos] == pix[j] + offset)
			out.push_back(pos);
	}

	if (out.size() == 0)
	{
		cout << "No pixels\n";
		cout << pix.size() << endl;
	}
}


///Tokenize a line
///@ingroup gUtility
vector<string> split(const string& line)
{
	vector<string> v;
	istringstream i(line);
	string s;

	while (!i.eof())
	{
		i >> s;
		if (i.fail())
			break;
		v.push_back(s);
	}
	return v;
}

///Generic version of itoa.
///How many times has this been reimplemented??
///@param x Value to convert to string
///@ingroup gUtility
template<class C> inline string xtoa(const C& x)
{
	ostringstream os;
	os << x;
	return os.str();
}

///Inverse of xtoa()
///How many times has this been reimplemented??
///@param s String to convert
///@param msg Mesage to print on error
///@ingroup gUtility
template<class C> inline C atox(const string& s, const string& msg)
{
	C c;
	istringstream i(s);
	i >> c;
	if (i.fail())
		throw LogFileParseError("Error parsing " + msg + ". Text is `" + s + "'.");
	return c;
}

/**Parser for multispot 5 log files.

Log files are mostly line oriented and contain various records

The main records are:

Iteraton: \#ITNUM
MAIN: \<list of spot parameters>

Pass: \#PASSNUM
MT19937 \<random number generator state>
PASS\#PASSNUM: \<list of spot parameters>
ENDCHECKPOINT

Note that MAIN is redundant since it will be the same as the following PASS 1 (or the first pass
computed if restoring from a checkpoint).

Data should only be considered valid after ENDCHECKPOINT has been read

Iteration is written once per iteration, not once per pass. (FIXME)

Which moron invented this file format???

Note that the file format hasn't beren fixed, so that the output can easily be compared to the
output of the historic version which is known to be good.

@param in Stream to parse file from
@ingroup gStorm
*/
StateParameters parse_log_file(istream& in)
{
	//A line read from the file
	string line;

	//State lines known to be OK
	string rngline, passline, iterationline;
	bool state_ok = 0;

	//State lines read in, with flags of goodness
	string new_rngline, new_passline, new_iterationline;
	bool new_rngline_ok = 0, new_passline_ok = 0, new_iterationline_ok = 0;

	unsigned int lineno = 0;
	bool doing_gvars = 0;

	vector<Ref3D> pixels;

	while (!in.eof())
	{
		getline(in, line);
		if (in.fail())
			break;

		lineno++;

		if (line == "ENDGVARLIST")
		{
			if (!doing_gvars)
				throw LogFileParseError("Spurious end of GVars");
			doing_gvars = 0;
		}
		else if (doing_gvars)
		{
			GUI.ParseLine(line);
		}
		else if (line == "BEGINGVARLIST")
		{
			doing_gvars = 1;
		}
		if (line.substr(0, 11) == "Iteration: ")
		{
			new_iterationline = line;
			new_iterationline_ok = true;
		}
		else if (line.substr(0, 4) == "PASS")
		{
			new_passline = line;
			if (new_passline_ok)
				throw LogFileParseError("Duplicate PASS on line " + xtoa(lineno));
			new_passline_ok = true;
		}
		else if (line.substr(0, 8) == "MT19937 ")
		{
			new_rngline = line;
			if (new_rngline_ok)
				throw LogFileParseError("Duplicate MT19937 on line " + xtoa(lineno));

			new_rngline_ok = true;

		}
		else if (line == "ENDCHECKPOINT")
		{
			if (new_passline_ok && new_rngline_ok && new_iterationline_ok)
			{
				iterationline = new_iterationline;
				rngline = new_rngline;
				passline = new_passline;
			}
			else
				throw LogFileParseError("Reached checkpoint with missing data: "
					"it=" + xtoa(new_iterationline_ok) +
					" pa=" + xtoa(new_passline_ok) +
					" rg=" + xtoa(new_rngline_ok) + " on line " + xtoa(lineno));

			//Don't reset iteration since it only appears once for each entire
			//set of passes. 
			new_rngline_ok = 0;
			new_passline_ok = 0;

			state_ok = true;
		}
		else if (line.substr(0, 7) == "PIXELS ")
		{
			vector<string> l = split(line);
			if ((l.size() - 1) % 3 == 0)
			{
				int n = (l.size() - 1) / 3;
				pixels.resize(n);
				for (int i = 0; i < n; i++)
				{
					pixels[i].x = atox<int>(l[i * 2 + 1 + 0], "pixels");
					pixels[i].y = atox<int>(l[i * 2 + 1 + 1], "pixels");
					pixels[i].z = atox<int>(l[i * 2 + 1 + 2], "pixels");
				}
			}
			else
				throw LogFileParseError("Bad PIXELS line");
		}
	}

	if (!state_ok)
		throw LogFileParseError("No state found");

	if (pixels.size() == 0)
		throw LogFileParseError("No pixels, or pixels is empty");

	//Now parse the lines
	StateParameters p;
	vector<string> l;

	//Parse the iterations
	l = split(iterationline);
	p.iteration = atox<int>(l[1], "iteration");

	//Parse the random number generator
	p.rng = shared_ptr<MT19937>(new MT19937);
	{
		istringstream rng_s(rngline);
		try {
			p.rng->read(rng_s);
		}
		catch (MT19937::ParseError p)
		{
			throw LogFileParseError("Error parsing MT19937");
		}
	}

	//Parse PASS and the listing of spots
	l = split(passline);
	if ((l.size() - 1) % 5 == 0)
	{
		p.pass = atox<int>(l[0].substr(5), "pass");	//YZL

		for (unsigned int i = 0; i < (l.size() - 1) / 5; i++)
		{
			cout << l[i * 5 + 1 + 0] << endl;
			cout << l[i * 5 + 1 + 1] << endl;
			cout << l[i * 5 + 1 + 2] << endl;
			cout << l[i * 5 + 1 + 3] << endl;
			cout << l[i * 5 + 1 + 4] << endl;
			p.spots.push_back(makeVector(
				atox<double>(l[i * 5 + 1 + 0], "spot"),
				atox<double>(l[i * 5 + 1 + 1], "spot"),
				atox<double>(l[i * 5 + 1 + 2], "spot"),
				atox<double>(l[i * 5 + 1 + 3], "spot"),
				atox<double>(l[i * 5 + 1 + 4], "spot")));
		}

	}
	else
		throw LogFileParseError("Wrong number of elements in PASS line");

	//Set up the pixels (oh god the pixels)
	p.pixels = pixels;

	return p;
}




//RFBA: Get Sphere for 3D Image Processing, Added@20180306
std::vector<Ref3D> getBall(float radius)
{
	std::vector<Ref3D> points;

	int r = (int)ceil(radius + 1);

	for (Ref3D p(0, 0, -r); p.z <= r; p.z++)
		for (p.y = -r; p.y <= r; p.y++)
			for (p.x = -r; p.x <= r; p.x++)
				if (p.mag_squared() <= radius*radius)
					points.push_back(p);

	return points;
}





///Very simple and inefficient dilation function.
///@ingroup gUtility
set<Ref3D>  dilate_mask(const vector<Ref3D>& v, double r)
{
	vector<Ref3D> m = getBall(r);
	set<Ref3D> ret;

	for (unsigned int i = 0; i < v.size(); i++)
		for (unsigned int j = 0; j < m.size(); j++)
			ret.insert(v[i] + m[j]);

	return ret;
}

///How far should steps in brightness be limited to?
///@ingroup gStorm
double brightness_motion_limit(double mu, double sigma, bool not_one)
{
	if (not_one)
		return log_normal_std(mu, sigma);
	else
		return 1;
}



struct OPT_PARA {
	void *fs;
	int id;
	vector<vector<double> > * spot_intensities; //[spot][pixel]
	vector<vector<int> > * spot_pixels; //[spot][pixel]
};

struct FE_PARA {
	void *fs;
	int id;
};

///Mega class which actually does the meat of the spot fitting.
///probably could be refactored a bit.
///@ingroup gStorm
class FitSpots
{
	const vector<Image3D<float> >& ims; ///< Input data
	FitSpotsGraphics& graphics;       ///< Graphics class.
	UserInterfaceCallback& ui;        ///< Callbacks to provide user interface
	const vector<Ref3D> pixels;    ///< Area in which to perform model fitting

	//Spot positions
	vector<Vector<5> > spots;         ///< State in terms of current spot positions

	//Clustering Result
	vector<vector<int>> clusters;		//cluseters
	vector<int> spots_to_erase;		//spots to erase
	pthread_mutex_t mutex;			//mutex


	//Starting point
	const int start_iteration;        ///< Starting iteration number (for restarting from checkpoint)
	int start_pass;                   ///< Starting pass (for restarting from checkpoint)

	MT19937& rng;                     ///< Random numbewr generator
	MT19937* rngMT;                     ///< Random numbewr generator


	const double variance;            ///< Variance of noise in data. Must be 1.
	const double intensity_mu;        ///< Prior for spot intensity
	const double intensity_sigma;     ///< Prior for spot intensity
	const double blur_mu;             ///< Prior for spot shape
	const double blur_sigma;          ///< Prior for spt shape

	//pixels is dilated slightly, by a disk of size area_extra_radius.
	//The result is stored in allowed_area. This is used to implement a uniform
	//prior over position, which is uniform within allowed_area, and zero
	//everywhere else.
	//
	//This is implemented as a set, because the mask may extend beyond the borders of the 
	//image slightly
	const double area_extra_radius;   ///< Extra size beyone marked region in which spots are allowed to exist
	set<Ref3D> allowed_area;       ///< Total allowed region, pixels dilated by area_extra_radius
	const int use_position_prior;     ///< Should a proper prior over position be uesd? A clue: yes.
	const double position_prior;      ///< Value for the posision prior, i.e. reciprocal of area

	//General optimizing and sampling parameters
	const double max_motion;          ///< Maximum motion on any axes for optimization. See ConjugateGradientOnly.
	const int sample_iterations;      ///< Number of mixing samples to use in Gibbs sampling

	//Task specific optimizing and sampling parameters:

	//Main optimization loop
	const int main_cg_max_iterations;  ///< Maximum iterations allowed for ConjugateGradientOnly in main optimization loop
	const int main_samples;            ///< Number of samples to use in main loop
	const int main_passes;             ///< Number of passes to perform per iteration in main loop
	const int outer_loop_iterations;   ///< Total number of iterations to perform

	//Spot selection loop
	const int add_remove_tries;        ///< Number of attemts to add/remove a spot
	const int add_remove_opt_samples;  ///< Number of samples to use in model modification phase
	const int add_remove_opt_retries;  ///< Number of attempts restarting the optimization to escape saddle points
	const int add_remove_opt_hess_inner_samples; ///< Number of extra FFBS samples to use for computing Hessian when testing for convergence to non-saddle point
	const int h_outer_samples;         ///< Number of samples used for computing Hessian as part of Laplace's approximation
	const int h_inner_samples;         ///< Number of additional FFBS samples to use for computing Hessian as part of Laplace's approximation
	const int tsamples;                ///< Number of samples to use in thermodynamic integration

	// Average of all inputs (used for debugging)
	const Image3D<float> ave;            ///< Average of input data: used for 


	//File to save output to
	ofstream& save_spots;               ///< Output stream for log file

	//Time accumulators for benchmarking
	double time_gibbs;                  ///< Benchmarking data
	double time_cg;                     ///< Benchmarking data

	///Motion limit for ConjugateGradientOnly
	///The x distance, y distance and spot size are all approximately the same scale
	///which is of order 1. The brigntness is not. By default, the brightness limit is also 1.
	///This flag controls whether is should be set to the standard deviation of the brightness prior
	///distribution. This setting will put the motion limit on the same scale as the other three
	///parameters.
	const bool scale_brightness_limit;
	const Vector<5> limit;              ///< Limit vector for ConjugateGradientOnly

	const Matrix<3> A;                  ///< Transition matrix 
	const Vector<3> pi;                 ///< Initial probabilities


	vector<vector<double> > pixel_intensities; ///<Pixel intensities for all images [frame][pixel]

	DataForMCMC data_for_t_mcmc;  ///<Aggergated data for thermodynamic integration
	DataForMCMC data_for_h_mcmc;  ///<Aggergated data for finding hessian

	vector<vector<int> > model2_spot_pixels_for_Calc;
	vector<Vector<5> >  model_2_for_Calc;
	vector<Vector<5> >  model_1_for_Calc;
	double model_1_energy;
	double model_2_energy;

	int iteration;                     ///< Iteration number
	//Guide
	int mag;
	vector<Ref3D> pixelsH;
	vector<float> intensitiesH;
	double max_intensity;

	TIME(cvd_timer timer;)
public:

	FitSpots(const vector<Image3D<float> >& ims_, FitSpotsGraphics& graphics_, UserInterfaceCallback& ui_, StateParameters& params, ofstream& save_spots_)
		:ims(ims_), graphics(graphics_), ui(ui_),

		//Set up the main parameters
		pixels(params.pixels),
		spots(params.spots),
		start_iteration(params.iteration),
		start_pass(params.pass),
		rng(*params.rng),



		//Set up all the system parameters
		variance(1), // it should be

		//To scale the X axis of a log-normal distribution, only
		//the mu parameter needs to be changed...
		intensity_mu(GV3::get<double>("intensity.rel_mu", 0., -1) + log(sqrt(variance))),
		intensity_sigma(GV3::get<double>("intensity.rel_sigma", 0., -1)),
		blur_mu(GV3::get<double>("blur.mu", 0., -1)),
		blur_sigma(GV3::get<double>("blur.sigma", 0., -1)),

		//Spot position prior
		area_extra_radius(GV3::get<double>("position.extra_radius", 0., -1)),
		allowed_area(dilate_mask(pixels, area_extra_radius)),
		use_position_prior(GV3::get<bool>("position.use_prior", true, -1)),
		position_prior(1.0 / allowed_area.size()),

		//General optimizing and sampling parameters
		max_motion(GV3::get<double>("cg.max_motion", 0., -1)),
		sample_iterations(GV3::get<int>("gibbs.mixing_iterations", 0, -1)),

		//Task specific optimizing and sampling parameters:

		//Main optimization loop
		main_cg_max_iterations(GV3::get<double>("main.cg.max_iterations", 0., -1)),
		main_samples(GV3::get<int>("main.gibbs.samples", 0, -1)),
		main_passes(GV3::get<int>("main.passes", 0, -1)),
		outer_loop_iterations(GV3::get<int>("main.total_iterations", 100000000, 1)),

		//Spot selection loop
		add_remove_tries(GV3::get<int>("add_remove.tries", 0, -1)),
		add_remove_opt_samples(GV3::get<int>("add_remove.optimizer.samples", 0, -1)),
		add_remove_opt_retries(GV3::get<int>("add_remove.optimizer.attempts", 0, -1)),
		add_remove_opt_hess_inner_samples(GV3::get<int>("add_remove.optimizer.hessian_inner_samples", 0, -1)),
		h_outer_samples(GV3::get<int>("add_remove.hessian.outer_samples", 0, -1)),
		h_inner_samples(GV3::get<int>("add_remove.hessian.inner_samples", 0, -1)),
		tsamples(GV3::get<int>("add_remove.thermo.samples", 0, -1)),

		ave(average_image(ims_)),

		save_spots(save_spots_),

		time_gibbs(0),
		time_cg(0),

		scale_brightness_limit(GV3::get<bool>("max_motion.use_brightness_std", 0, -1)),
		limit(makeVector(brightness_motion_limit(intensity_mu, intensity_sigma, scale_brightness_limit), 1, 1, 1, 1)*max_motion),

		A(GV3::get<Matrix<3> >("A", Zeros, 1)),
		pi(GV3::get<Vector<3> >("pi", Zeros, 1)),


		data_for_t_mcmc(pixels, pixel_intensities, intensity_mu, intensity_sigma, blur_mu, blur_sigma, variance, tsamples, sample_iterations, A, pi, rng),
		data_for_h_mcmc(pixels, pixel_intensities, intensity_mu, intensity_sigma, blur_mu, blur_sigma, variance, h_outer_samples, sample_iterations, A, pi, rng)

	{
		//assert_same_size(ims);

		//Pixel intensities for all images [frame][pixel]
		pixel_intensities.resize(ims.size(), vector<double>(pixels.size()));
		for (unsigned int frame = 0; frame < ims.size(); frame++)
			for (unsigned int p = 0; p < pixels.size(); p++)
				pixel_intensities[frame][p] = ims[frame][pixels[p]];

		//inti rngMT
		rngMT = (MT19937 *)malloc(MAXTHREAD * sizeof(MT19937));
		for (int i = 0; i < MAXTHREAD; i++) {
			(rngMT + i)->simple_seed(rng.rand_int());
		}

	}

	///Perform a complete iteration of the optimzation stage of the spot firrint algorithm
	void optimize_each_spot_in_turn_for_several_passes()
	{
		//Precompute the intensities for all spot pixels
		vector<vector<double> > spot_intensities; //[spot][pixel]
		for (unsigned int i = 0; i < spots.size(); i++)
			spot_intensities.push_back(compute_spot_intensity(pixels, spots[i], zFactor));

		//Which pixels does each spot have?
		vector<vector<int> > spot_pixels; //[spot][pixel]
		spot_pixels.resize(spots.size());
		for (unsigned int s = 0; s < spots.size(); s++)
			get_spot_pixels(pixels, spots[s], spot_pixels[s]);


		//Optimize the model, N spots at a time.
		//
		for (int pass = start_pass; pass < main_passes; pass++)
		{
			save_spots << "Pass: " << pass << endl;
			rng.write(save_spots);
			save_spots << endl;

			start_pass = 0; // This is only nonzero first time, since we might chekpoint midway through an iteration
			save_spots << "PASS" << pass << ": " << setprecision(20) << scientific << spots_to_Vector(spots) << endl;
			save_spots << "ENDCHECKPOINT" << endl << flush;

			//ui.per_pass(iteration, pass, spots);



			//Sort spots according to pass%4

			//Sort the spots so that the optimization runs in sweeps
			//This heiristic seems to increase the rate of propagation of information
			//about spot positions.

			//Create a list of indices
			vector<int> index = sequence(spots.size());

			int passs = pass + iteration;
			//Sort the indices according to the position of the spot that they index
			if (passs % 4 == 0)
				sort(index.begin(), index.end(), IndexLexicographicPosition<less<double>, 2>(spots));
			else if (passs % 4 == 1)
				sort(index.begin(), index.end(), IndexLexicographicPosition<greater<double>, 2>(spots));
			else if (passs % 4 == 1)
				sort(index.begin(), index.end(), IndexLexicographicPosition<less<double>, 3>(spots));
			else
				sort(index.begin(), index.end(), IndexLexicographicPosition<greater<double>, 3>(spots));

			//Reorder the spots and their intensities and their pixel lists
			{
				vector<Vector<5> > tmp_spot(index.size());
				vector<vector<double> > tmp_spot_intensities(index.size());
				vector<vector<int> > tmp_spot_pixels(index.size());
				for (unsigned int i = 0; i < index.size(); i++)
				{
					tmp_spot[i] = spots[index[i]];
					swap(tmp_spot_intensities[i], spot_intensities[index[i]]);
					swap(tmp_spot_pixels[i], spot_pixels[i]);
				}

				swap(tmp_spot, spots);
				swap(tmp_spot_intensities, spot_intensities);
				swap(tmp_spot_pixels, spot_pixels);
			}



			//Sweep through and optimize each spot in turn
			for (int s = 0; s < (int)spots.size(); s++)
			{
				ui.per_spot(iteration, pass, s, spots.size());
				ui.perhaps_stop();

				TIME(timer.reset();)
					//Get some samples with Gibbs sampling
					vector<vector<vector<State> > > sample_list; //[N][spot][frame]: list of samples drawn using Gibbs sampling
				vector<vector<vector<double> > > sample_intensities; //[sample][frame][pixel]

				GibbsSampler2 sampler(pixel_intensities, spot_intensities, spots, spot_pixels, A, pi, variance, sample_iterations);
				for (int i = 0; i < main_samples; i++)
				{
					sampler.next(rng);
					sample_list.push_back(sampler.sample());
					sample_intensities.push_back(sampler.sample_intensities());

					ui.perhaps_stop();
				}

				//First, remove the spot from all the samples.
				for (unsigned int i = 0; i < sample_list.size(); i++)
					remove_spot(sample_intensities[i], spot_intensities[s], sample_list[i][s]);

				//cout << timer.get_time() << endl;
				TIME(time_gibbs += timer.reset();)

					//Package up all the data
					SampledBackgroundData data(sample_intensities, pixel_intensities, pixels,
						intensity_mu, intensity_sigma, blur_mu, blur_sigma,
						A, pi, variance);

				//Derivative computer:
				SpotNegProbabilityDiffWithSampledBackground compute_deriv(data);


				//graphics.draw_pixels(pixels, 0, 0, 1, 1);
				//graphics.draw_krap(spots, scale_to_bytes(ave), boundingbox(pixels), s);
				//graphics.swap();

				//Optimize spot "s"
				//Optimize with the derivatives only since the actual probability
				//is much harder to compute
				ConjugateGradientOnly<5> cg(spots[s], compute_deriv, limit);


				cg.max_iterations = main_cg_max_iterations;


				//cout << "Starting opt... " << cg.x << endl;
				while (cg.iterate(compute_deriv))
				{
					//graphics.draw_krap(spots, scale_to_bytes(ave), boundingbox(pixels), s, cg.x);
					//graphics.draw_pixels(pixels, 0, 0, 1, .2);
					//graphics.swap();
					ui.perhaps_stop();
					//cout << cg.x << endl;
				}

				//Update to use the result of the optimization
				spots[s] = cg.x;
				//cout << "End: " << cg.x << endl;

				//graphics.draw_krap(spots, scale_to_bytes(ave), boundingbox(pixels), -1);
				//graphics.swap();

				//Recompute the new spot intensity, since the spot has changed
				spot_intensities[s] = compute_spot_intensity(pixels, spots[s], zFactor);

				//Recompute which are the useful pixels
				get_spot_pixels(pixels, spots[s], spot_pixels[s]);

				//Is the spot within the allowed area, i.e. is it's prior 0?
				//The prior is sero only if it we are using it and we're in an invalid area

				Ref3D quantized_spot_position = ir_rounded3D(spots[s].slice<2, 3>());
				bool zero_prior = use_position_prior && (allowed_area.count(quantized_spot_position) == 0);

				//Check to see if spot has been ejected. If spot_pixels is empty, then it has certainly been ejected.
				if (spot_pixels[s].empty() || zero_prior)
				{
					//Spot has been ejected. Erase it.
					cout << " Erasing ejected spot: " << spot_pixels[s].empty() << " " << zero_prior << endl;
					cout << spots[s] << endl;
					//GUI_Pause(1);

					spot_intensities.erase(spot_intensities.begin() + s);
					spot_pixels.erase(spot_pixels.begin() + s);
					spots.erase(spots.begin() + s);
					s--;
					//exit(0);
				}



				//cout << timer.get_time() << endl;
				TIME(time_cg += timer.reset();)


					//cout << "Times: " << time_gibbs << " " << time_cg << endl;
					//save_spots << "INTERMEDIATE: " << setprecision(15) << scientific << spots_to_Vector(spots) << endl;
			}
		}
	}


	//RFBA: Multithreading, Added@20180306
	void thread_handler_opt1(int id, vector<vector<double> > spot_intensitiesall, vector<vector<int> > spot_pixelsall)
	{
		vector<int> mycluster = clusters[id];
		//printf("Opt1 Cluster: %d, Points: %d\n", id, mycluster.size());

		vector<Vector<5> > myspots;         ///< State in terms of current spot positions
		for (int k = 0; k < (int)spots.size(); k++)
			myspots.push_back(spots[k]);

		//Precompute the intensities for all spot pixels
		vector<vector<double> > myspot_intensities; //[spot][pixel]
		for (int k = 0; k < (int)myspots.size(); k++)
			myspot_intensities.push_back(compute_spot_intensity(pixels, myspots[k], zFactor));

		//Which pixels does each spot have?
		vector<vector<int> > myspot_pixels; //[spot][pixel]
		myspot_pixels.resize(myspots.size());
		for (int k = 0; k < (int)myspots.size(); k++)
			get_spot_pixels(pixels, myspots[k], myspot_pixels[k]);


		//Sweep through and optimize each spot in cluster
		for (int id = 0; id < (int)mycluster.size(); id++)
		{
			int k = mycluster[id];
			//Get some samples with Gibbs sampling
			vector<vector<vector<State> > > sample_list; //[N][spot][frame]: list of samples drawn using Gibbs sampling
			vector<vector<vector<double> > > sample_intensities; //[sample][frame][pixel]

			GibbsSampler2 sampler(pixel_intensities, myspot_intensities, myspots, myspot_pixels, A, pi, variance, sample_iterations);
			for (int i = 0; i < main_samples; i++)
			{
				sampler.next(rng);
				sample_list.push_back(sampler.sample());
				sample_intensities.push_back(sampler.sample_intensities());

			}

			//First, remove the spot from all the samples.
			for (unsigned int i = 0; i < sample_list.size(); i++)
				remove_spot(sample_intensities[i], myspot_intensities[k], sample_list[i][k]);



			//Package up all the data
			SampledBackgroundData data(sample_intensities, pixel_intensities, pixels,
				intensity_mu, intensity_sigma, blur_mu, blur_sigma,
				A, pi, variance);

			//Derivative computer:
			SpotNegProbabilityDiffWithSampledBackground compute_deriv(data);


			//Optimize spot "s"
			//Optimize with the derivatives only since the actual probability
			//is much harder to compute
			ConjugateGradientOnly<5> cg(myspots[k], compute_deriv, limit);


			cg.max_iterations = main_cg_max_iterations;


			while (cg.iterate(compute_deriv))
			{

			}

			//Update to use the result of the optimization
			myspots[k] = cg.x;

			//Recompute the new spot intensity, since the spot has changed
			myspot_intensities[k] = compute_spot_intensity(pixels, myspots[k], zFactor);

			//Recompute which are the useful pixels
			get_spot_pixels(pixels, myspots[k], myspot_pixels[k]);

			//Is the spot within the allowed area, i.e. is it's prior 0?
			//The prior is sero only if it we are using it and we're in an invalid area

			Ref3D quantized_spot_position = ir_rounded3D(myspots[k].slice<2, 3>());
			bool zero_prior = use_position_prior && (allowed_area.count(quantized_spot_position) == 0);

			//Check to see if spot has been ejected. If spot_pixels is empty, then it has certainly been ejected.
			if (myspot_pixels[k].empty() || zero_prior)
			{

				pthread_mutex_lock(&mutex);
				//Spot has been ejected. Erase it.
				cout << " To Erasing ejected spot: " << myspot_pixels[k].empty() << " " << zero_prior << endl;
				cout << myspots[k] << endl;
				//spot_intensities.erase(spot_intensities.begin() + s);
				//spot_pixels.erase(spot_pixels.begin() + s);
				//spots.erase(spots.begin() + s);
				//s--;
				spots_to_erase.push_back(mycluster[k]);
				pthread_mutex_unlock(&mutex);
			}
		}
		//update
		pthread_mutex_lock(&mutex);
		for (int k = 0; k < (int)mycluster.size(); k++)
		{
			spot_intensitiesall[mycluster[k]] = myspot_intensities[k];
			spot_pixelsall[mycluster[k]] = myspot_pixels[k];
		}
		pthread_mutex_unlock(&mutex);
		//printf("Opt1 Cluster Done: %d\n", id);

	}
	//RFBA: Multithreading, Added@20180306
	static void * run_optimizer_step1(void* para)
	{
		OPT_PARA *p = (OPT_PARA *)para;
		int id = p->id;
		FitSpots * fs = (FitSpots *)p->fs;
		fs->thread_handler_opt1(id, *(p->spot_intensities), *(p->spot_pixels));
		return(0);

	}

	//RFBA: Multithreading, Added@20180306
	void thread_handler_fe(int id)
	{
		if (id == 0) {
			//Model 2
			model_2_energy = -NegativeFreeEnergy(data_for_t_mcmc).compute_with_mask(spots_to_Vector(model_2_for_Calc), model2_spot_pixels_for_Calc);
		}
		else {
			//Model1
			vector<vector<int> > model1_spot_pixels_for_Calc;
			model1_spot_pixels_for_Calc.reserve(model2_spot_pixels_for_Calc.size());
			for (int i = 0; i < model2_spot_pixels_for_Calc.size() - 1; i++)
				model1_spot_pixels_for_Calc.push_back(vector<int>(model2_spot_pixels_for_Calc[i]));
			model_1_energy = -NegativeFreeEnergy(data_for_t_mcmc).compute_with_mask(spots_to_Vector(model_1_for_Calc), model1_spot_pixels_for_Calc);
		}
	}

	//RFBA: Multithreading, Added@20180306
	static void * run_fe(void* para)
	{
		OPT_PARA *p = (OPT_PARA *)para;
		int id = p->id;
		FitSpots * fs = (FitSpots *)p->fs;
		fs->thread_handler_fe(id);
		return(0);
	}





	///Perform a complete iteration of the model size modification stage of the spot fitting algorithm
	void try_modifying_model()
	{

		for (int i = 0; i < add_remove_tries; i++)
		{
			ui.per_modification(iteration, i, add_remove_tries);
			ui.perhaps_stop();

			cout << endl << endl << "Modifying the model" << endl << "======================\n";
			//cout << "Hello\n";
			bool add_spot = (rng() > 0.5) || (spots.size() == 1);
			//cout << "World\n";

			vector<Vector<5> > model_1, model_2;
			int original; //What is the original model? Model 1 or Model 2?

			if (add_spot)
			{
				model_1 = spots;
				model_2 = model_1;

				//Pick a pixel within the thresholded ones as a starting point
				if (addGuide == 0) {
					int r;
					do
					{
						r = (int)(rng() * pixels.size());
						//cout << r << " " << log_ratios[pixels[r]] << " " << pixels[r] << " " << threshold << endl;
					} while (0);

					model_2.push_back(makeVector(log_normal_mode(intensity_mu, intensity_sigma),
						log_normal_mode(blur_mu, blur_sigma),
						pixels[r].x + rng() - .5, pixels[r].y + rng() - .5, pixels[r].z + rng() - .5));

				}
				else {
					//Add using srrf guide
					float x = 0;
					float y = 0;
					float z = 0;
					GetPoint(&x, &y, &z);
					printf("Point Guided by SRRF: %f, %f, %f!\n", x, y, z);

					model_2.push_back(makeVector(log_normal_mode(intensity_mu, intensity_sigma),
						log_normal_mode(blur_mu, blur_sigma),
						x, y, z));

				}

				//This version does not (yet?) suppotrt thresholding on log_ratios
				//for placing spots, since the purpose of log_ratios has diminished.
				//while(log_ratios[pixels[r]] < threshold);


				cout << "Adding a spot\n";

				original = 1;
			}
			else
			{
				//Pick a random point to optimize/remove
				int a_random_spot = static_cast<int>(rng() * spots.size());
				model_1 = spots;
				swap(model_1[model_1.size() - 1], model_1[a_random_spot]);
				model_2 = model_1;
				model_1.pop_back();
				cout << "Removing a spot\n";
				original = 2;
			}

			//The mobile spot is always the last spot of model_2
			const int spot = model_2.size() - 1;

			cout << "Original model: " << original << endl;

			//Precompute the intensities for all spot pixels
			//Model 2 is always one longer than model 1 and 
			//differs only on the extra element
			vector<vector<double> > model2_spot_intensities; //[spot][pixel]
			for (unsigned int i = 0; i < model_2.size(); i++)
				model2_spot_intensities.push_back(compute_spot_intensity(pixels, model_2[i], zFactor));

			//Which pixels does each spot have?
			vector<vector<int> > model2_spot_pixels(model_2.size()); //[spot][pixel]
			for (unsigned int s = 0; s < model_2.size(); s++)
				get_spot_pixels(pixels, model_2[s], model2_spot_pixels[s]);

			//Optimize spot:
			{
				cout << "Optimizing spot for model selection\n";


				//Get some samples with Gibbs sampling
				vector<vector<vector<State> > > sample_list; //[N][spot][frame]: list of samples drawn using Gibbs sampling
				vector<vector<vector<double> > > sample_intensities; //[sample][frame][pixel]

				GibbsSampler2 sampler(pixel_intensities, model2_spot_intensities, model_2, model2_spot_pixels, A, pi, variance, sample_iterations);
				for (int i = 0; i < add_remove_opt_samples; i++)
				{
					sampler.next(rng);
					sample_list.push_back(sampler.sample());
					sample_intensities.push_back(sampler.sample_intensities());
					ui.perhaps_stop();
				}

				//First, remove the spot from all the samples.
				for (unsigned int i = 0; i < sample_list.size(); i++)
					remove_spot(sample_intensities[i], model2_spot_intensities[spot], sample_list[i][spot]);

				//Package up all the data
				SampledBackgroundData data(sample_intensities, pixel_intensities, pixels,
					intensity_mu, intensity_sigma, blur_mu, blur_sigma,
					A, pi, variance);

				//Derivative computer:
				SpotNegProbabilityDiffWithSampledBackground compute_deriv(data);

				//graphics.draw_krap(model_2, scale_to_bytes(ave), boundingbox(pixels), spot);
				//graphics.swap();

				//Optimize spot "s"
				//Optimize with the derivatives only since the actual probability
				//is much harder to compute
				ConjugateGradientOnly<5> cg(model_2[spot], compute_deriv, limit);


				for (int attempt = 0; attempt < add_remove_opt_retries; attempt++)
				{
					cout << "Attempt " << attempt << " of " << add_remove_opt_retries << endl;
					ui.perhaps_stop();

					//Optimize with conjugate gradient
					while (cg.iterate(compute_deriv))
					{
						ui.perhaps_stop();
						//graphics.draw_krap(model_2, scale_to_bytes(ave), boundingbox(pixels), spot, cg.x);
						//graphics.swap();
					}

					if (attempt < add_remove_opt_retries - 1)
					{
						Matrix<5> hessian = sampled_background_spot_hessian_ffbs(cg.x, data, add_remove_opt_hess_inner_samples, rngMT);
						SymEigen<5> hess_decomp(hessian);


						cout << "Eigenvalues are: " << hess_decomp.get_evalues() << endl;

						if (hess_decomp.is_negdef())
							break;
						else
						{
							//Restart in the direction of the best uphill part
							cg.init(cg.x + 0.1 * hess_decomp.get_evectors()[3], (hess_decomp.get_evectors()[3]));


						}
					}
				}


				//Update to use the result of the optimization
				model_2[spot] = cg.x;

				//graphics.draw_krap(model_2, scale_to_bytes(ave), boundingbox(pixels), -1);
				//graphics.swap();


				//Update cached data based on the new spot position
				model2_spot_intensities[spot] = compute_spot_intensity(pixels, model_2[spot], zFactor);
				get_spot_pixels(pixels, model_2[spot], model2_spot_pixels[spot]);

				cout << "Done optimizing for model selection\n";
			}


			//Which model to keep?
			int keep = original;

			//Compute position prior (and we might be able to reject it really quickly here)
			bool zero_prior = use_position_prior && (allowed_area.count(ir_rounded3D(model_2[spot].slice<2, 3>())) == 0);

			if (zero_prior)
			{
				//Model 2 went bad, so we clearly keep model 1
				keep = 1;
			}
			else
			{

				//The position prior is independent
				//Compute the difference  model2 - model1
				//This is only valid if model2 is in the valid region
				double position_log_prior_model2_minus_model1;
				if (use_position_prior)
					position_log_prior_model2_minus_model1 = (model_2.size() - model_1.size()) * ln(position_prior);
				else
					position_log_prior_model2_minus_model1 = 0;


				//Integrate model_2
				//First compute the Hessian since this might go wrong.

				//FreeEnergyHessian hesscomputer(data_for_h_mcmc);
				Matrix<5> hess;// = hesscomputer.hessian(model_2, spot);



				//Use turbohess here since it is much faster, as the backwards sampling step is fast
				//We expect this hessian to be really quite different, if the spot has moved from
				//a long way away, since the sampling will have changed dramatically
				{
					//Get some samples with Gibbs sampling
					vector<vector<vector<State> > > sample_list; //[N][spot][frame]: list of samples drawn using Gibbs sampling
					vector<vector<vector<double> > > sample_intensities; //[sample][frame][pixel]

					GibbsSampler sampler(pixel_intensities, model2_spot_intensities, model_2, A, pi, variance, sample_iterations);
					for (int i = 0; i < h_outer_samples; i++)
					{
						ui.perhaps_stop();
						sampler.next(rng);
						sample_list.push_back(sampler.sample());
						sample_intensities.push_back(sampler.sample_intensities());
					}

					//First, remove the spot from all the samples.
					for (unsigned int i = 0; i < sample_list.size(); i++)
						remove_spot(sample_intensities[i], model2_spot_intensities[spot], sample_list[i][spot]);

					//Package up all the data
					SampledBackgroundData data(sample_intensities, pixel_intensities, pixels,
						intensity_mu, intensity_sigma, blur_mu, blur_sigma,
						A, pi, variance);


					hess = sampled_background_spot_hessian_ffbs(model_2[spot], data, h_inner_samples, rngMT);

				}



				double det = determinant(-hess / (M_PI * 2));
				SymEigen<5> hess_decomp(-hess);
				cout << "Hessien Eigenvalues are: " << hess_decomp.get_evalues() << endl;
				const double smallest_evalue = 1e-6;

				//If the determinant is negative, then we are still at a saddle point
				//despite the attempts above, so abandon the attempt.
				if (hess_decomp.get_evalues()[0] > smallest_evalue)
				{

					//Compute the free energy of model 2 at the MLE estimate
					cout << "Model 2:\n";

					model2_spot_pixels_for_Calc = model2_spot_pixels;
					model_2_for_Calc = model_2;
					model_1_for_Calc = model_1;

					pthread_t tid1, tid2;
					FE_PARA param[2];
					param[0].fs = this;
					param[0].id = 0;
					param[1].fs = this;
					param[1].id = 1;
					pthread_create(&tid1, NULL, run_fe, (void *)(param));
					pthread_create(&tid2, NULL, run_fe, (void *)(param + 1));

					pthread_join(tid1, NULL);
					pthread_join(tid2, NULL);

					cout << "Energy: " << model_2_energy << endl;
					//Combine the MLE energy and Hessian using Laplace's approximation
					double model_2_occam = -0.5 * log(det);
					double model_2_prob = model_2_energy + model_2_occam + position_log_prior_model2_minus_model1;

					cout << "Occam: " << model_2_occam << endl;
					cout << "Position: " << position_log_prior_model2_minus_model1 << endl;
					cout << "Prob: " << model_2_prob << endl;


					double model_1_prob = model_1_energy;

					cout << "Prob: " << model_1_prob << endl;


					model2_spot_pixels.pop_back();


					if (model_2_prob > model_1_prob)
						keep = 2;
					else
						keep = 1;

					cout << "Models evaluated\n";
				}
				else
				{
					cout << "Determinant has bad eigenvalues!\n";
					keep = original;
					cout << hess_decomp.get_evalues() << endl;
				}
			}

			if (keep == 2)
			{
				spots = model_2;
				cout << "Keeping model 2\n";

			}
			else
			{
				spots = model_1;
				cout << "Keeping model 1\n";
			}

			if (original != keep)
			{
				cout << "Model changed!\n";
				//break;
			}
		}
	}



	void InitGuide() {

		int imageDim[4];
		imageDim[0] = ims[0].size().x;
		imageDim[1] = ims[0].size().y;
		imageDim[2] = ims[0].size().z;
		imageDim[3] = ims.size();

		//Build new pixels
		mag = SRRFGetMag();
		pixelsH.reserve(pixels.size() * mag * mag *mag);
		for (int i = 0; i < pixels.size(); i++) {
			int x = pixels[i].x;
			int y = pixels[i].y;
			int z = pixels[i].z;
			for (int j1 = 0; j1 < mag; j1++) {
				for (int j2 = 0; j2 < mag; j2++) {
					for (int j3 = 0; j3 < mag; j3++) {
						Ref3D ref = Ref3D(x*mag + j1, y*mag + j2, z*mag + j3);
						pixelsH.push_back(ref);
					}
				}
			}
		}
		//Load
		float *outImage = (float *)malloc(imageDim[0] * imageDim[1] * imageDim[2] * mag * mag *mag * sizeof(float));
		FILE * fp = fopen("SRRF.raw", "rb");
		fread(outImage, sizeof(float), imageDim[0] * imageDim[1] * imageDim[2] * mag * mag *mag, fp);
		fclose(fp);

		intensitiesH.reserve(pixelsH.size());
		for (unsigned int i = 0; i < pixelsH.size(); i++) {
			int x = pixelsH[i].x;
			int y = pixelsH[i].y;
			int z = pixelsH[i].z;
			int idx = z*imageDim[0] * imageDim[1] * mag * mag + y*imageDim[0] * mag + x;
			intensitiesH[i] = outImage[idx];
		}
		free(outImage);

		max_intensity = *max_element(intensitiesH.begin(), intensitiesH.end());


	}


	void GetPoint(float * px, float *py, float *pz) {
		while (1)
		{
			int element = floor(rng() * pixelsH.size());
			double y = rng() * max_intensity;

			if (y <= intensitiesH[element]) {
				float x = pixelsH[element].x + rng() - .5;
				float y = pixelsH[element].y + rng() - .5;
				float z = pixelsH[element].z + rng() - .5;
				x = x / mag;
				y = y / mag;
				z = z / mag;
				*px = x;
				*py = y;
				*pz = z;
				break;
			}
		}

	}


	///Run the complete optimization algorithm.
	void run()
	{
		//graphics.init(ims[0].size());
		save_spots << "LOGVERSION 2" << endl;
		save_spots << "PIXELS";
		for (unsigned int i = 0; i < pixels.size(); i++)
			save_spots << " " << pixels[i].x << " " << pixels[i].y << " " << pixels[i].z;
		save_spots << endl;

		save_spots << "BEGINGVARLIST" << endl;
		GV3::print_var_list(save_spots, "", 1);
		save_spots << "ENDGVARLIST" << endl;

		//TODO all GVARS are set, so dump out gvars.

		cout << "Limit vector: " << limit << endl;

		//SRRF Guide
		if (addGuide == 1) {
			printf("Add: SRRF Guide!\n");
			InitGuide();
		}


		printf("MAX Iter: %d\n", maxIter);
		for (iteration = start_iteration; iteration < maxIter; iteration++)
		{
			save_spots << "Iteration: " << iteration << " (" << iteration *  main_passes << ")" << endl;
			save_spots << "MAIN: " << setprecision(20) << scientific << spots_to_Vector(spots) << endl;

			cout << endl << endl << "----------------------" << endl << "Optimizing:\n";
			cout << spots.size() << endl;

			optimize_each_spot_in_turn_for_several_passes();

			//spot_intensities is be correct here!
			try_modifying_model();


		}
		save_spots << "FINAL: " << setprecision(15) << scientific << spots_to_Vector(spots) << endl;
	}

};

///Wrapper function for using FitSpots
///@ingroup gStorm
void fit_spots_new(const vector<Image3D<float> >& ims, StateParameters& p, ofstream& save_spots, FitSpotsGraphics& gr)
{
	auto_ptr<UserInterfaceCallback> ui = null_ui();
	FitSpots fit(ims, gr, *ui, p, save_spots);
	fit.run();
}

///Wrapper function for using FitSpots
///@ingroup gStorm
void fit_spots_new(const vector<Image3D<float> >& ims, StateParameters& p, ofstream& save_spots, FitSpotsGraphics& gr, UserInterfaceCallback& ui)
{
	try {
		FitSpots fit(ims, gr, ui, p, save_spots);
		fit.run();
	}
	catch (UserInterfaceCallback::UserIssuedStop)
	{
	}
}




