#pragma once

#ifndef ROOT_SCALAR_HPP
#define ROOT_SCALAR_HPP

#include <functional>
#include <cmath>
#include <limits>
#include <algorithm>
#include <stdexcept>

namespace RootScalar
{
	/**
	 * @brief Checks if two floating-point numbers are approximately equal.
	 *
	 * @tparam T The type of the numbers (e.g., float, double).
	 * @param a The first number.
	 * @param b The second number.
	 * @param precision The number of decimal places to compare.
	 * @return true if the numbers are approximately equal, false otherwise.
	 */
	template <typename T>
	bool almost_equal(T a, T b, unsigned int precision)
	{
		static_assert(std::is_floating_point<T>::value, "T must be a floating-point type.");
		return std::abs(a - b) < std::pow(10, -static_cast<int>(precision));
	}

	/**
	 * @brief Computes the numerical derivative of a function at a given point.
	 *
	 * @tparam T The type of the input and output values.
	 * @tparam Function The type of the function to differentiate.
	 * @param f The function to differentiate.
	 * @param x The point at which to compute the derivative.
	 * @return The derivative of the function at the given point.
	 */
	template <typename T, typename Function>
	T derivate(Function f, T x)
	{
		static_assert(std::is_floating_point<T>::value, "T must be a floating-point type.");
		static_assert(std::is_invocable<Function, T>::value, "Function must be invocable with T.");
		static_assert(std::is_floating_point<typename std::result_of<Function(T)>::type>::value,
			"Function must return a floating-point type.");
		T h = std::sqrt(std::numeric_limits<T>::epsilon());
		return (f(x + h) - f(x - h)) / (2.0 * h);
	}

	/**
	 * @brief Finds a root of a function using the Newton-Raphson method.
	 *
	 * @tparam T The type of the input and output values.
	 * @tparam Function The type of the function whose root is to be found.
	 * @param f The function whose root is to be found.
	 * @param x0 The initial guess for the root.
	 * @param precision The desired precision of the result.
	 * @param nmax The maximum number of iterations.
	 * @return The root of the function.
	 * @throws std::runtime_error If the derivative is zero or the maximum number of iterations is reached.
	 */
	template <typename T, typename Function>
	T newton(Function f, T x0, unsigned int precision, unsigned int nmax)
	{
		static_assert(std::is_floating_point<T>::value, "T must be a floating-point type.");
		static_assert(std::is_invocable<Function, T>::value, "Function must be invocable with T.");
		static_assert(std::is_floating_point<typename std::result_of<Function(T)>::type>::value,
			"Function must return a floating-point type.");

		T x;
		T f_prime;
		T pertubation = std::sqrt(std::numeric_limits<T>::epsilon());

		for (unsigned int i = 0; i < nmax; i++)
		{
			f_prime = derivate(f, x0);

			if (f_prime == 0.0)
			{
				x0 += pertubation;
				f_prime = derivate(f, x0);
				if (f_prime == 0.0)
				{
					throw std::runtime_error("Derivative is zero. No root found.");
				}
			}

			x = x0 - f(x0) / f_prime;

			if (almost_equal(f(x), static_cast<T>(0.0), precision) || almost_equal(x, x0, precision))
			{
				return x;
			}

			x0 = x;
		}
		throw std::runtime_error("Maximum number of iterations reached without finding root.");
	}

	/**
	 * @brief Finds a root of a function using the bisection method.
	 *
	 * @tparam T The type of the input and output values.
	 * @tparam Function The type of the function whose root is to be found.
	 * @param f The function whose root is to be found.
	 * @param a The lower bound of the interval.
	 * @param b The upper bound of the interval.
	 * @param precision The desired precision of the result.
	 * @param nmax The maximum number of iterations.
	 * @return The root of the function.
	 * @throws std::invalid_argument If the function values at the interval endpoints have the same sign.
	 * @throws std::runtime_error If the maximum number of iterations is reached.
	 */
	template <typename T, typename Function>
	T bisection(Function f, T a, T b, unsigned int precision, unsigned int nmax)
	{
		static_assert(std::is_floating_point<T>::value, "T must be a floating-point type.");
		static_assert(std::is_invocable<Function, T>::value, "Function must be invocable with T.");
		static_assert(std::is_floating_point<typename std::result_of<Function(T)>::type>::value,
			"Function must return a floating-point type.");

		if (std::signbit(f(a)) == std::signbit(f(b)))
		{
			throw std::invalid_argument("Function values at the interval endpoints must have opposite signs.");
		}
				
		T c;

		for (unsigned int i = 0; i < nmax; i++)
		{
			c = (a + b) / 2;
			if (almost_equal(f(c), static_cast<T>(0.0), precision) || almost_equal(a, b, precision))
			{
				return c;
			}

			if (std::signbit(f(a)) != std::signbit(f(c)))
			{
				b = c;
			}
			else
			{
				a = c;
			}
		}
		throw std::runtime_error("Maximum number of iterations reached without finding root.");
	}

	/**
	 * @brief Finds a root of a function using the secant method.
	 *
	 * @tparam T The type of the input and output values.
	 * @tparam Function The type of the function whose root is to be found.
	 * @param f The function whose root is to be found.
	 * @param x0 The first initial guess for the root.
	 * @param x1 The second initial guess for the root.
	 * @param precision The desired precision of the result.
	 * @param nmax The maximum number of iterations.
	 * @return The root of the function.
	 * @throws std::runtime_error If division by zero is encountered or the maximum number of iterations is reached.
	 */
	template <typename T, typename Function>
	T secant(Function f, T x0, T x1, unsigned int precision, unsigned int nmax)
	{
		static_assert(std::is_floating_point<T>::value, "T must be a floating-point type.");
		static_assert(std::is_invocable<Function, T>::value, "Function must be invocable with T.");
		static_assert(std::is_floating_point<typename std::result_of<Function(T)>::type>::value,
			"Function must return a floating-point type.");

		if (f(x1) == f(x0))
		{
			throw std::runtime_error("Division by zero encountered in secant method.");
		}
		
		for (unsigned int i = 0; i < nmax; i++)
		{
			T x = x1 - f(x1) * ((x1 - x0) / (f(x1) - f(x0)));
			x0 = x1;
			x1 = x;

			if (almost_equal(f(x1), static_cast<T>(0.0), precision))
			{
				return x1;
			}
		}
		throw std::runtime_error("Maximum number of iterations reached without finding root.");
	}

	/**
	 * @brief Finds a root of a function using Brent's method.
	 *
	 * @tparam T The type of the input and output values.
	 * @tparam Function The type of the function whose root is to be found.
	 * @param f The function whose root is to be found.
	 * @param a The lower bound of the interval.
	 * @param b The upper bound of the interval.
	 * @param precision The desired precision of the result.
	 * @param nmax The maximum number of iterations.
	 * @return The root of the function.
	 * @throws std::invalid_argument If the function values at the interval endpoints have the same sign.
	 * @throws std::runtime_error If the maximum number of iterations is reached without finding a root.
	 */
	template <typename T, typename Function>
	T brent(Function f, T a, T b, unsigned int precision, unsigned int nmax)
	{
		static_assert(std::is_floating_point<T>::value, "T must be a floating-point type.");
		static_assert(std::is_invocable<Function, T>::value, "Function must be invocable with T.");
		static_assert(std::is_floating_point<typename std::result_of<Function(T)>::type>::value,
			"Function must return a floating-point type.");

		if (std::signbit(f(a)) == std::signbit(f(b)))
		{
			throw std::invalid_argument("Function values at the interval endpoints must have opposite signs.");
		}
		
		if (f(a) < f(b))
		{
			std::swap(a, b);
		}

		T c = a;
		bool mflag = true;
		T s;
		T d;

		for (int i = 0; i < nmax; i++)
		{
			if ((f(a) != f(c)) && (f(b) != f(c)))
			{
				s = a * f(b) * f(c) / ((f(a) - f(b)) * (f(a) - f(c))) +
					b * f(a) * f(c) / ((f(b) - f(a)) * (f(b) - f(c))) +
					c * f(a) * f(b) / ((f(c) - f(a)) * (f(c) - f(b)));
			}
			else
			{
				s = b - f(b) * ((b - a) / (f(b) - f(a)));
			}

			if (!((((3 * a + b) / 4.0) < s) && (s < b)) ||
				(mflag && (std::abs(s - b) >= std::abs((b - c) / 2.0))) ||
				(!mflag && (std::abs(s - b) >= std::abs((c - d) / 2.0))) ||
				(mflag && almost_equal(b, c, precision)) ||
				(!mflag && almost_equal(c, d, precision)))
			{
				s = (a + b) / 2.0;
				mflag = true;
			}
			else
			{
				mflag = false;
			}

			d = c;
			c = b;

			if (std::signbit(f(a)) == std::signbit(f(s)))
			{
				a = s;
			}
			else
			{
				b = s;
			}

			if (std::abs(f(a)) < std::abs(f(b)))
			{
				std::swap(a, b);
			}

			if (almost_equal(f(b), static_cast<T>(0.0), precision) || almost_equal(a, b, precision))
			{
				return b;
			}
		}
		throw std::runtime_error("Maximum number of iterations reached without finding root.");
	}
}

#endif // ROOT_SCALAR_HPP