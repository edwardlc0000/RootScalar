#pragma once

#ifndef ROOT_SCALAR_HPP
#define ROOT_SCALAR_HPP

#include <functional>
#include <cmath>
#include <limits>
#include <stdexcept>

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

	double newton(real_func f, double x0, int precision, int nmax)
	{
		double x = x0;
		double f_prime;
		double pertubation = std::sqrt(std::numeric_limits<double>::epsilon());

		for(int i = 0; i < nmax; i++)
		{
			f_prime = derivate(f, x);

			if (f_prime == 0.0)
			{
				x += pertubation;
				f_prime = derivate(f, x);
				if (f_prime == 0.0)
				{
					throw std::runtime_error("Derivative is zero. No root found.");
				}
			}

			double xi = x - f(x) / f_prime;

			if (almost_equal(f(xi), 0.0, precision) || (std::abs(xi - x) < std::pow(10, -precision)))
			{
				return xi;
			}

			x = xi;
		}
		throw std::runtime_error("Maximum number of iterations reached without finding root.");
	}

	double bisection(real_func f, double a, double b, int precision, int nmax, int i = 0)
	{
		i++;

		if (std::signbit(f(a)) == std::signbit(f(b)))
		{
			throw std::invalid_argument("Function values at the interval endpoints must have opposite signs.");
		}

		double c = (a + b) / 2;
		if (almost_equal(f(c), 0.0, precision) || (std::abs(b - a) < std::pow(10, -precision)))
		{
			return c;
		}

		if (i < nmax)
		{
			if (std::signbit(f(a)) != std::signbit(f(c)))
			{
				return bisection(f, a, c, precision, nmax, i);
			}
			else
			{
				return bisection(f, b, c, precision, nmax, i);
			}
		}
		else
		{
			throw std::runtime_error("Maximum number of iterations reached without finding root.");
		}
	}

	double secant(real_func f, double x0, double x1, int precision, int nmax)
	{
		if (f(x1) == f(x0))
		{
			throw std::runtime_error("Division by zero encountered in secant method.");
		}

		for (int i = 0; i < nmax; i++)
		{
			double xi = x1 - f(x1) * ((x1 - x0) / (f(x1) - f(x0)));
			x0 = x1;
			x1 = xi;

			if (almost_equal(f(x1), 0.0, precision) || (std::abs(x1 - x0) < std::pow(10, -precision)))
			{
				return x1;
			}

		}
		throw std::runtime_error("Maximum number of iterations reached without finding root.");
	}

}

#endif // ROOT_SCALAR_HPP
