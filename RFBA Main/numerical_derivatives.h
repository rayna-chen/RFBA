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

#ifndef NUMERICAL_DERIVATIVES_H
#define NUMERICAL_DERIVATIVES_H
#include <TooN/TooN.h>

template<int S, class B, class Functor> TooN::Vector<S> debug_numerical_gradient(const Functor& f, const TooN::Vector<S,double, B>& x)
{
	using namespace TooN;
	Vector<S> grad(x.size());
	Vector<S> xh=x;
	const double h=1e-4;

	for(int i=0; i < x.size(); i++)
	{
		xh[i] += h;
		double fwd = f(xh);
		xh[i] = x[i] - h;
		double rev = f(xh);

		grad[i] = (fwd - rev) / (2*h);
		xh[i] = x[i];
	}

	return grad;
}

template<int S, class B, class Functor> TooN::Matrix<S> debug_numerical_hessian(const Functor& f, const TooN::Vector<S,double, B>& x)
{

	using namespace TooN;
	Matrix<S> hess(x.size(), x.size());
	Vector<S> xh=x;
	const double h=1e-3;

	for(int i=0; i < x.size(); i++)
		for(int j=i+1; j < x.size(); j++)
		{
			xh = x;
			xh[i] += h;
			xh[j] += h;
			double a = f(xh);

			xh = x;
			xh[i] -= h;
			xh[j] += h;
			double b = f(xh);

			xh = x;
			xh[i] += h;
			xh[j] -= h;
			double c = f(xh);

			xh = x;
			xh[i] -= h;
			xh[j] -= h;
			double d = f(xh);

			hess[i][j] = hess[j][i] = (a-b-c+d) / (4*h*h);
		}

	for(int i=0; i < x.size(); i++)
	{
		xh = x;
		double b = f(xh);

		xh[i] += h;
		double a = f(xh);
		
		xh[i] = x[i] - h;
		double c = f(xh);
		hess[i][i] = (a - 2*b + c) / (h*h);
	}

	return hess;
}


#endif
