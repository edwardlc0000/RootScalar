#pragma once

#ifndef ROOT_SCALAR_HPP
#define ROOT_SCALAR_HPP

#include <functional>
#include <cmath>
#include <limits>

namespace RootScalar
{
	using real_func = std::function<double(double)>;
	using real_func_derivative = std::function<double(real_func, double)>;

	bool almost_equal(double a, double b, int precision)
	{
		return std::abs(a - b) < std::pow(10, -precision);
	}

	double derivate(real_func f, double x)
	{
		double h = std::sqrt(std::numeric_limits<double>::epsilon());
		return (f(x + h) - f(x - h)) / (2.0 * h);
	}

	double newton(real_func f, double x0, int precision)
	{
		double x = x0;
		double f_prime;
		do
		{
			f_prime = derivate(f, x);
			x = x - f(x) / f_prime;
		} while (!almost_equal(f(x), 0.0, precision));
		return x;
	}
}

#endif // ROOT_SCALAR_HPP
