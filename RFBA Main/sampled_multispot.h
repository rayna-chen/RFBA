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

#ifndef SPOT_WITH_BACKGROUND_H
#define SPOT_WITH_BACKGROUND_H

#include <vector>
#include <cvd/image_ref.h>
#include <tuple>
#include <TooN/TooN.h>
#include <pthread.h>
#include "mt19937.h"
#include "image3d.h"

#include "drift.h"

typedef char State;

namespace SampledMultispot
{

	using namespace std;
	using namespace CVD;
	using namespace TooN;
	using namespace std::tr1;



	//The changes to SpotWithBackground to operate with/without drift and masking are
	//irritatingly pervasive. I can't think of any other way of doing it.
#define SWBG_NAME SpotWithBackground
#include "spot_with_background.hh"


#define SWBG_NAME SpotWithBackgroundMasked
#define SWBG_HAVE_MASK
#include "spot_with_background.hh"

//#define SWBG_NAME SpotWithBackgroundDrift
//#define SWBG_HAVE_DRIFT
//#include "spot_with_background.hh"

//#undef SWBG_NAME
//#define SWBG_NAME SpotWithBackgroundDriftMasked
//#define SWBG_HAVE_DRIFT
//#define SWBG_HAVE_MASK
//#include "spot_with_background.hh"


	//RFBA: Multithreading Suppport Struct, Added@20180306
	struct GIBBS_PARA {
		void * pGibbsSampler2;
		int idx;
		int step;
		int * thread_finished;
		int * mainReady;
		pthread_mutex_t * mutex;
		pthread_cond_t * cond_c;
		pthread_cond_t * cond_p;
		MT19937 * rng;
	};




	inline double intensity(double i)
	{
		return i;
	}
	//RFBA: Update Parameter for 3D Image, Modified@20180306
	inline double intensity(const pair<double, Vector<5> >& i)
	{
		return i.first;
	}


	//Add and remove a spot over the entire region
	template<class T>
	void remove_spot(vector<vector<double> >& current_sample_intensities, const vector<T>& spot_intensities, const vector<State>& spot_sample)
	{
		for (unsigned int frame = 0; frame < current_sample_intensities.size(); frame++)
			if (spot_sample[frame] == 0) //Spot is on, so remove it
				for (unsigned int p = 0; p < spot_intensities.size(); p++)
					current_sample_intensities[frame][p] -= intensity(spot_intensities[p]);
	}

	template<class T>
	void add_spot(vector<vector<double> >& current_sample_intensities, const vector<T>& spot_intensities, const vector<State>& spot_sample)
	{
		for (unsigned int frame = 0; frame < current_sample_intensities.size(); frame++)
			if (spot_sample[frame] == 0) //Spot is on, so add it
				for (unsigned int p = 0; p < spot_intensities.size(); p++)
					current_sample_intensities[frame][p] += intensity(spot_intensities[p]);
	}


	//Add and remove a spot only over a mask. 
	template<class T>
	void remove_spot(vector<vector<double> >& current_sample_intensities, const vector<T>& spot_intensities, const vector<State>& spot_sample, const vector<int>& mask)
	{
		for (unsigned int frame = 0; frame < current_sample_intensities.size(); frame++)
			if (spot_sample[frame] == 0) //Spot is on, so remove it
				for (unsigned int p = 0; p < mask.size(); p++)
					current_sample_intensities[frame][mask[p]] -= intensity(spot_intensities[mask[p]]);
	}

	template<class T>
	void add_spot(vector<vector<double> >& current_sample_intensities, const vector<T>& spot_intensities, const vector<State>& spot_sample, const vector<int>& mask)
	{
		for (unsigned int frame = 0; frame < current_sample_intensities.size(); frame++)
			if (spot_sample[frame] == 0) //Spot is on, so add it
				for (unsigned int p = 0; p < mask.size(); p++)
					current_sample_intensities[frame][mask[p]] += intensity(spot_intensities[mask[p]]);
	}


	//Add and remove a drifty spot only over a mask.
	template<class T>
	void remove_spot(vector<vector<double> >& current_sample_intensities, const vector<vector<T> > & spot_intensities, const vector<State>& spot_sample, const vector<int>& mask)
	{
		const int steps = spot_intensities.size();
		const int frames = current_sample_intensities.size();

		for (int frame = 0; frame < frames; frame++)
		{
			int s = frame * steps / frames;

			if (spot_sample[frame] == 0) //Spot is on, so remove it
				for (unsigned int p = 0; p < mask.size(); p++)
					current_sample_intensities[frame][mask[p]] -= intensity(spot_intensities[s][mask[p]]);
		}
	}

	template<class T>
	void add_spot(vector<vector<double> >& current_sample_intensities, const vector<vector<T> >& spot_intensities, const vector<State>& spot_sample, const vector<int>& mask)
	{
		const int steps = spot_intensities.size();
		const int frames = current_sample_intensities.size();

		for (int frame = 0; frame < frames; frame++)
		{
			int s = frame * steps / frames;

			if (spot_sample[frame] == 0) //Spot is on, so add it
				for (unsigned int p = 0; p < mask.size(); p++)
					current_sample_intensities[frame][mask[p]] += intensity(spot_intensities[s][mask[p]]);
		}
	}

	//RFBA: Add for 3D Image, Added@20180306
	/// Convert an image co-ordinate into a Vector
	/// @param ir The ImageRef to convert
	/// @ingroup gImage
	inline TooN::Vector<3> vec3D(const Ref3D& ir)
	{
		TooN::Vector<3> r;
		r[0] = ir.x;
		r[1] = ir.y;
		r[2] = ir.z;
		return r;
	}


	//RFBA: Update for 3D Image, Modified@20180306
	//Compute the spot intensity for a given spot at each pixel
	//size brightness x y z 
	inline vector<double> compute_spot_intensity(const vector<Ref3D>& pixels, const Vector<5>& params, double zFactor)
	{
		vector<double> intensities(pixels.size());

		for (unsigned int i = 0; i < pixels.size(); i++)
			intensities[i] = spot_shape(vec3D(pixels[i]), params, zFactor);

		return intensities;
	}

	//RFBA: Update for 3D Image, Modified@20180306
	//Compute the spot intensity derivatives for a given spot at each pixel
	inline vector<pair<double, Vector<5> > > compute_spot_intensity_derivatives(const vector<Ref3D>& pixels, const Vector<5>& params, double zFactor)
	{
		vector<pair<double, Vector<5> > > derivatives(pixels.size());

		for (unsigned int i = 0; i < pixels.size(); i++)
			derivatives[i] = spot_shape_diff_position(vec3D(pixels[i]), params, zFactor);
		return derivatives;
	}

	//RFBA: Update for 3D Image, Modified@20180306
	inline vector<tuple<double, Vector<5>, Matrix<5> > > compute_spot_intensity_hessian(const vector<Ref3D>& pixels, const Vector<5>& params, double zFactor)
	{
		vector<tuple<double, Vector<5>, Matrix<5> > > hessian(pixels.size());

		for (unsigned int i = 0; i < pixels.size(); i++)
			hessian[i] = spot_shape_hess_position(vec3D(pixels[i]), params, zFactor);
		return hessian;
	}


	/**
	Create a sequence of integers. These can be used as observations
	in an observation class by forward_algorithm() and etc.
	@param n Length of sequence
	@ingroup gUtility
	*/
	inline vector<int> sequence(int n)
	{
		vector<int> v;
		for (int i = 0; i < n; i++)
			v.push_back(i);
		return v;
	}

	/*struct RndGrand48
	{
		double operator()()
		{
			return drand48();
		}
	};*/

	///Draw samples from the spot states given the spots positions and some data.
	///Variable naming matches that in FitSpots.
	///@ingroup gStorm
	class GibbsSampler
	{
		const vector<vector<double> >& pixel_intensities;
		const vector<vector<double> >& spot_intensities;
		const vector<Vector<5> > spots;
		const Matrix<3> A;
		const Vector<3> pi;
		const double base_variance;
		double variance;

		const int sample_iterations;
		const int num_frames, num_pixels;
		const vector<int> O;

		vector<vector<State> > current_sample;
		vector<vector<double> > current_sample_intensities;

	public:

		GibbsSampler(const vector<vector<double> >& pixel_intensities_,
			const vector<vector<double> >& spot_intensities_,
			const vector<Vector<5> >& spots_,
			const Matrix<3> A_,
			const Vector<3> pi_,
			double variance_,
			int sample_iterations_)
			:pixel_intensities(pixel_intensities_), //pixel_intensities: [frame][pixels]
			spot_intensities(spot_intensities_),   //spot_intensities: [spot][pixel]
			spots(spots_),
			A(A_),
			pi(pi_),
			base_variance(variance_),
			variance(variance_),
			sample_iterations(sample_iterations_),
			num_frames(pixel_intensities.size()),
			num_pixels(pixel_intensities[0].size()),
			//Observations vector. As usual for this application, the observations are just
			//numbered integers which refer to data held elsewhere.
			O(sequence(num_frames)),
			//Start all spots OFF, so the intensity is 0. OFF is 1 or 2, not 0!!!
			//sample_list: [sample][spot][frame]: list of samples drawn using Gibbs sampling
			current_sample(spots.size(), vector<State>(num_frames, 2)), //current sample [spot][frame]
			   //pixel intensities assosciated with the current sample [frame][pixel]
			current_sample_intensities(num_frames, vector<double>(num_pixels))
		{
			//Check a bunch of stuff
			assert_same_size(pixel_intensities);
			assert_same_size(spot_intensities);

		}

		///Update the noide variance. Used for adding thermal noise.
		///@param v noise variance.
		void set_variance(double v)
		{
			variance = v;
		}


		///Reset the gibbs sampler oro the initial state (all spots off)
		void reset()
		{
			vector<State> off(num_frames, 2);
			fill(current_sample.begin(), current_sample.end(), off);

			vector<double> black(num_pixels);
			fill(current_sample_intensities.begin(), current_sample_intensities.end(), black);
			variance = base_variance;
		}

		///Get the next sample
		///@param rng Random number generator
		template<class T> void next(T& rng)
		{

			for (int j = 0; j < sample_iterations; j++)
				for (int k = 0; k < (int)spots.size(); k++)
				{
					//Subtract off the spot we're interested in.
					remove_spot(current_sample_intensities, spot_intensities[k], current_sample[k]);

					//Now current_sample_intensities is the image value for every spot in every frame,
					//except the current spot, which is always set to off. This allows us to add it in 
					//easily.
					SpotWithBackground B(current_sample_intensities, spot_intensities[k], pixel_intensities, variance);
					vector<array<double, 3> > delta = forward_algorithm_delta(A, pi, B, O);
					current_sample[k] = backward_sampling<3, State, T>(A, delta, rng);

					//Put the newly sampled spot in
					add_spot(current_sample_intensities, spot_intensities[k], current_sample[k]);
				}
		}
		/*void next()
		{
			RngDrand48 rng;
			next(rng);
		}*/

		///Retrieve the current sample
		const vector<vector<State> >& sample() const
		{
			return current_sample;
		}
		///Retrieve the intensities for the current sample
		const vector<vector<double> >& sample_intensities() const
		{
			return current_sample_intensities;
		}

	};






	///Gibbs sampling class which masks spots to reduce computation.
	///
	///This draws samples from, the spot states given the spots positions and some data. It is
	///very similar to GibbsSampler, except that it only computes probabilities in a mask around each spot
	///to save on computation. Variable naming matches that in FitSpots.
	///@ingroup gStorm
	class GibbsSampler2
	{
		const vector<vector<double> >& pixel_intensities;
		const vector<vector<double> >& spot_intensities;
		const vector<Vector<5> > spots;
		const std::vector<std::vector<int> >& spot_pixels;
		const Matrix<3> A;
		const Vector<3> pi;
		const double base_variance;
		double variance;

		const int sample_iterations;
		const int num_frames, num_pixels;
		const vector<int> O;

		vector<vector<State> > current_sample;
		vector<vector<double> > current_sample_intensities;

		vector<double> cutout_spot_intensities;
		vector<vector<double> > cutout_pixel_intensities;
		vector<vector<double> > cutout_current_sample_intensities;

	public:

		GibbsSampler2(const vector<vector<double> >& pixel_intensities_,
			const vector<vector<double> >& spot_intensities_,
			const vector<Vector<5> >& spots_,
			const vector<vector<int> >& spot_pixels_,
			const Matrix<3> A_,
			const Vector<3> pi_,
			double variance_,
			int sample_iterations_)
			:pixel_intensities(pixel_intensities_), //pixel_intensities: [frame][pixels]
			spot_intensities(spot_intensities_),   //spot_intensities: [spot][pixel]
			spots(spots_),
			spot_pixels(spot_pixels_),
			A(A_),
			pi(pi_),
			base_variance(variance_),
			variance(variance_),
			sample_iterations(sample_iterations_),
			num_frames(pixel_intensities.size()),
			num_pixels(pixel_intensities[0].size()),
			//Observations vector. As usual for this application, the observations are just
			//numbered integers which refer to data held elsewhere.
			O(sequence(num_frames)),
			//Start all spots OFF, so the intensity is 0. OFF is 1 or 2, not 0!!!
			//sample_list: [sample][spot][frame]: list of samples drawn using Gibbs sampling
			current_sample(spots.size(), vector<State>(num_frames, 2)), //current sample [spot][frame]
			   //pixel intensities assosciated with the current sample [frame][pixel]
			current_sample_intensities(num_frames, vector<double>(num_pixels)),
			cutout_pixel_intensities(num_frames),
			cutout_current_sample_intensities(num_frames)
		{
			//Check a bunch of stuff
			assert_same_size(pixel_intensities);
			assert_same_size(spot_intensities);
		}

		///Update the noide variance. Used for adding thermal noise.
		///@param v noise variance.
		void set_variance(double v)
		{
			variance = v;
		}


		///Reset the gibbs sampler oro the initial state (all spots off)
		void reset()
		{
			vector<State> off(num_frames, 2);
			fill(current_sample.begin(), current_sample.end(), off);

			vector<double> black(num_pixels);
			fill(current_sample_intensities.begin(), current_sample_intensities.end(), black);
			variance = base_variance;
		}



		//RFBA: Simple Multithreading, Added@20180306
		static void  * gibss_handler(void* para) {
			GIBBS_PARA * mypara = (GIBBS_PARA *)para;
			int step = mypara->step;
			GibbsSampler2 * pgs = (GibbsSampler2 *)(mypara->pGibbsSampler2);
			pthread_mutex_t  mutex = *(mypara->mutex);
			pthread_cond_t cond_p = *(mypara->cond_p);
			pthread_cond_t cond_c = *(mypara->cond_c);
			int * thread_finished = mypara->thread_finished;
			int *mainReady = mypara->mainReady;

			std::vector<array<double, 3> > delta3;

			vector<vector<double> > current_sample_intensities;
			//Copy first
			current_sample_intensities.reserve(pgs->current_sample_intensities.size());


			for (int j = 0; j < pgs->sample_iterations; j++) {

				printf("tid: %d, start\n", mypara->idx);
				pthread_mutex_lock(&mutex);
				printf("tid: %d, wait main\n", mypara->idx);
				while ((*(mainReady + mypara->idx)) == 0) {
					pthread_cond_wait(&cond_c, &mutex);
				}
				printf("tid: %d, get main\n", mypara->idx);
				printf("tid: %d, finished %d\n", mypara->idx, *thread_finished);
				pthread_mutex_unlock(&mutex);
				printf("tid: %d, release go\n", mypara->idx);

				//update sample intensities
				for (unsigned int frame = 0; frame < pgs->current_sample_intensities.size(); frame++) {
					current_sample_intensities.push_back(vector<double>(pgs->current_sample_intensities[frame]));
				}

				for (int k = mypara->idx; k < pgs->spots.size(); k = k + step) {
					//Subtract off the spot we're interested in.
					remove_spot(current_sample_intensities, pgs->spot_intensities[k], pgs->current_sample[k], pgs->spot_pixels[k]);

					SpotWithBackgroundMasked B3(current_sample_intensities, pgs->spot_intensities[k], pgs->pixel_intensities, pgs->variance, pgs->spot_pixels[k]);
					forward_algorithm_delta2<3>(pgs->A, pgs->pi, B3, pgs->O, delta3);
					pgs->current_sample[k] = backward_sampling<3, State, MT19937>(pgs->A, delta3, *(mypara->rng));

					//Put the newly sampled spot in
					add_spot(current_sample_intensities, pgs->spot_intensities[k], pgs->current_sample[k], pgs->spot_pixels[k]);
				}
				pthread_mutex_lock(&mutex);
				(*thread_finished)++;
				*(mainReady + mypara->idx) = 0;
				//Inform for Main Thread 
				printf("tid: %d, finished %d\n", mypara->idx, *thread_finished);
				if (*thread_finished >= step) {
					pthread_cond_signal(&cond_p);
					printf("tid: %d, send\n", mypara->idx);
				}
				pthread_mutex_unlock(&mutex);
				printf("tid: %d, lopp %d\n", mypara->idx, j);
			}
			//Join
			printf("EXIT!\n");
			pthread_exit(NULL);
			return(0);
		}

		//RFBA: Simple Multithreading, Modified@20180306
		///Get the next sample
		///@param rng Random number generator
		template<class T> void next_new(T* rng)
		{
			int iter = 0;
			vector<vector<State> > current_sample_old;
			current_sample_old.reserve(current_sample.size());
			for (unsigned int i = 0; i < current_sample.size(); i++) {
				current_sample_old.push_back(vector<State>(current_sample[i]));
			}


			pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
			pthread_cond_t cond_p = PTHREAD_COND_INITIALIZER;
			pthread_cond_t cond_c = PTHREAD_COND_INITIALIZER;

			int finished_thread = 0;

			int finished_iter = 0;
			int thread_num = (int)spots.size();
			if (thread_num > MAXTHREAD)
				thread_num = MAXTHREAD;

			int *mainReady = (int *)malloc(thread_num * sizeof(int));

			GIBBS_PARA * fp = (GIBBS_PARA *)malloc(thread_num * sizeof(GIBBS_PARA));
			for (int i = 0; i < thread_num; i++) {
				(fp + i)->idx = i;
				(fp + i)->step = thread_num;
				(fp + i)->pGibbsSampler2 = (void *)this;
				(fp + i)->mutex = &mutex;
				(fp + i)->cond_c = &cond_c;
				(fp + i)->cond_p = &cond_p;
				(fp + i)->thread_finished = &finished_thread;
				(fp + i)->mainReady = mainReady;
				(fp + i)->rng = rng + i;
			}
			pthread_t *tid = (pthread_t *)malloc(thread_num * sizeof(pthread_t));

			for (unsigned int i = 0; i < thread_num; i++) {
				pthread_create((tid + i), NULL, gibss_handler, (void *)(fp + i));
			}


			pthread_cond_broadcast(&cond_c);

			while (finished_iter < sample_iterations) {
				printf("main: finished iter: %d\n", finished_iter);
				pthread_mutex_lock(&mutex);
				printf("main get!\n");

				pthread_cond_broadcast(&cond_c);
				for (int i = 0; i < thread_num; i++) {
					*(mainReady + i) = 1;
				}

				while (finished_thread < thread_num) {
					printf("main: wait: %d\n", finished_thread);
					pthread_cond_wait(&cond_p, &mutex);
					printf("main: wait2: %d\n", finished_thread);
				}

				printf("main go!\n");
				//update intensities
				for (int i = 0; i < thread_num; i++) {
					for (int k = i; k < spots.size(); k = k + thread_num) {
						remove_spot(current_sample_intensities, spot_intensities[k], current_sample_old[k], spot_pixels[k]);
						add_spot(current_sample_intensities, spot_intensities[k], current_sample[k], spot_pixels[k]);
					}
				}
				//keep old
				current_sample_old.clear();
				for (unsigned int i = 0; i < current_sample.size(); i++) {
					current_sample_old.push_back(vector<State>(current_sample[i]));
				}
				finished_thread = 0;
				finished_iter++;
				printf("main done!\n");
				pthread_mutex_unlock(&mutex);
			}

			free(mainReady);
			free(fp);
			free(tid);
		}



		//RFBA: Simple Multithreading, Added@20180306
		static void  * gibss_handler2(void* para) {
			GIBBS_PARA * mypara = (GIBBS_PARA *)para;
			int step = mypara->step;
			GibbsSampler2 * pgs = (GibbsSampler2 *)(mypara->pGibbsSampler2);
			std::vector<array<double, 3> > delta3;

			vector<vector<double> > current_sample_intensities;
			//Copy first
			current_sample_intensities.reserve(pgs->current_sample_intensities.size());


			for (int j = 0; j < pgs->sample_iterations; j++) {

				//update sample intensities
				for (unsigned int frame = 0; frame < pgs->current_sample_intensities.size(); frame++) {
					current_sample_intensities.push_back(vector<double>(pgs->current_sample_intensities[frame]));
				}

				for (int k = mypara->idx; k < pgs->spots.size(); k = k + step) {
					//Subtract off the spot we're interested in.
					remove_spot(current_sample_intensities, pgs->spot_intensities[k], pgs->current_sample[k], pgs->spot_pixels[k]);

					SpotWithBackgroundMasked B3(current_sample_intensities, pgs->spot_intensities[k], pgs->pixel_intensities, pgs->variance, pgs->spot_pixels[k]);
					forward_algorithm_delta2<3>(pgs->A, pgs->pi, B3, pgs->O, delta3);
					pgs->current_sample[k] = backward_sampling<3, State, MT19937>(pgs->A, delta3, *(mypara->rng));

					//Put the newly sampled spot in
					add_spot(current_sample_intensities, pgs->spot_intensities[k], pgs->current_sample[k], pgs->spot_pixels[k]);
				}
			}
			//Join		
			pthread_exit(NULL);
			return(0);
		}

		//RFBA: Simple Multithreading, Modified@20180306
		///Get the next sample
		///@param rng Random number generator
		template<class T> void next_2(T* rng)
		{
			int iter = 0;
			vector<vector<State> > current_sample_old;
			current_sample_old.reserve(current_sample.size());
			for (unsigned int i = 0; i < current_sample.size(); i++) {
				current_sample_old.push_back(vector<State>(current_sample[i]));
			}

			int thread_num = (int)spots.size();
			if (thread_num > MAXTHREAD)
				thread_num = MAXTHREAD;


			GIBBS_PARA * fp = (GIBBS_PARA *)malloc(thread_num * sizeof(GIBBS_PARA));
			for (int i = 0; i < thread_num; i++) {
				(fp + i)->idx = i;
				(fp + i)->step = thread_num;
				(fp + i)->pGibbsSampler2 = (void *)this;
				(fp + i)->rng = rng + i;
			}
			pthread_t *tid = (pthread_t *)malloc(thread_num * sizeof(pthread_t));

			for (unsigned int i = 0; i < thread_num; i++) {
				pthread_create((tid + i), NULL, gibss_handler2, (void *)(fp + i));
			}

			//Wait
			for (unsigned int i = 0; i < thread_num; i++) {
				pthread_join(*(tid + i), NULL);
			}

			//update intensities
			for (int k = 0; k < spots.size(); k = k + 1) {
				remove_spot(current_sample_intensities, spot_intensities[k], current_sample_old[k], spot_pixels[k]);
				add_spot(current_sample_intensities, spot_intensities[k], current_sample[k], spot_pixels[k]);
			}
			free(fp);
			free(tid);
		}
		



		///Get the next sample
		///@param rng Random number generator
		template<class T> void next(T& rng)
		{

			//double remove=0;
			//double cut=0;
			//double swb=0;
			//double ff_masked=0;
			//double bs=0;
			//double add=0;
			//cvd_timer t;
			std::vector<array<double, 3> > delta3;
			for (int j = 0; j < sample_iterations; j++)
				for (int k = 0; k < (int)spots.size(); k++)
				{					
					//Subtract off the spot we're interested in.
					remove_spot(current_sample_intensities, spot_intensities[k], current_sample[k], spot_pixels[k]);

					//Now current_sample_intensities is the image value for every spot in every frame,
					//except the current spot, which is always set to off. This allows us to add it in 
					//easily.
					SpotWithBackgroundMasked B3(current_sample_intensities, spot_intensities[k], pixel_intensities, variance, spot_pixels[k]);
					forward_algorithm_delta2<3>(A, pi, B3, O, delta3);	
					current_sample[k] = backward_sampling<3, State, T>(A, delta3, rng);	
					//Put the newly sampled spot in
					add_spot(current_sample_intensities, spot_intensities[k], current_sample[k], spot_pixels[k]);
				
				}
			//		cout << "remove=" <<remove << " cut=" << cut << " swb=" << swb<< " ff_mask=" << ff_masked << " bs=" <<bs << " add="<<add << endl;
		}



		/*	void next()
		{
			RngDrand48 rng;
			next(rng);
		}
	*/
	///Retrieve the current sample
		const vector<vector<State> >& sample() const
		{
			return current_sample;
		}
		///Retrieve the intensities for the current sample
		const vector<vector<double> >& sample_intensities() const
		{
			return current_sample_intensities;
		}

	};



}

using SampledMultispot::SpotWithBackground;
using SampledMultispot::remove_spot;
using SampledMultispot::add_spot;
using SampledMultispot::compute_spot_intensity;
using SampledMultispot::compute_spot_intensity_hessian;
using SampledMultispot::compute_spot_intensity_derivatives;
using SampledMultispot::sequence;
using SampledMultispot::GibbsSampler;
using SampledMultispot::GibbsSampler2;

#endif
