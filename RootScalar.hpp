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
    using real_func = std::function<double(double)>;
    using real_func_derivative = std::function<double(real_func, double)>;

    /**
     * @brief Checks if two floating-point numbers are almost equal.
     * 
     * @param a First number.
     * @param b Second number.
     * @param precision Number of decimal places to consider.
     * @return true if the numbers are almost equal, false otherwise.
     */
    bool almost_equal(double a, double b, unsigned int precision)
    {
        return std::abs(a - b) < std::pow(10, -static_cast<int>(precision));
    }

    /**
     * @brief Computes the numerical derivative of a function at a given point.
     * 
     * @param f The function to differentiate.
     * @param x The point at which to compute the derivative.
     * @return The numerical derivative of the function at the given point.
     */
    double derivate(real_func f, double x)
    {
        double h = std::sqrt(std::numeric_limits<double>::epsilon());
        return (f(x + h) - f(x - h)) / (2.0 * h);
    }

    /**
     * @brief Finds the root of a function using the Newton-Raphson method.
     * 
     * @param f The function for which to find the root.
     * @param x0 Initial guess for the root.
     * @param precision Number of decimal places to consider for convergence.
     * @param nmax Maximum number of iterations.
     * @return The root of the function.
     * @throws std::runtime_error if the root is not found within the maximum number of iterations.
     */
    double newton(real_func f, double x0, unsigned int precision, int nmax)
    {
        double x;
        double f_prime;
        double pertubation = std::sqrt(std::numeric_limits<double>::epsilon());
        double tolerance = std::pow(10, -static_cast<int>(precision));

        for(int i = 0; i < nmax; i++)
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

            if (almost_equal(f(x), 0.0, precision) || (std::abs(x - x0) < tolerance))
            {
                return x;
            }

            x0 = x;
        }
        throw std::runtime_error("Maximum number of iterations reached without finding root.");
    }

    /**
     * @brief Finds the root of a function using the bisection method.
     * 
     * @param f The function for which to find the root.
     * @param a The lower bound of the interval.
     * @param b The upper bound of the interval.
     * @param precision Number of decimal places to consider for convergence.
     * @param nmax Maximum number of iterations.
     * @return The root of the function.
     * @throws std::invalid_argument if the function values at the interval endpoints have the same sign.
     * @throws std::runtime_error if the root is not found within the maximum number of iterations.
     */
    double bisection(real_func f, double a, double b, unsigned int precision, int nmax, int i = 0)
    {
        if (std::signbit(f(a)) == std::signbit(f(b)))
        {
            throw std::invalid_argument("Function values at the interval endpoints must have opposite signs.");
        }

        double tolerance = std::pow(10, -static_cast<int>(precision));
        double c;

        for (int i = 0; i < nmax; i++)
        {
            c = (a + b) / 2;
            if (almost_equal(f(c), 0.0, precision) || (std::abs(b - a) < tolerance))
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
     * @brief Finds the root of a function using the secant method.
     * 
     * @param f The function for which to find the root.
     * @param x0 First initial guess for the root.
     * @param x1 Second initial guess for the root.
     * @param precision Number of decimal places to consider for convergence.
     * @param nmax Maximum number of iterations.
     * @return The root of the function.
     * @throws std::runtime_error if the root is not found within the maximum number of iterations or if division by zero is encountered.
     */
    double secant(real_func f, double x0, double x1, unsigned int precision, int nmax)
    {
        if (f(x1) == f(x0))
        {
            throw std::runtime_error("Division by zero encountered in secant method.");
        }

        double tolerance = std::pow(10, -static_cast<int>(precision));

        for (int i = 0; i < nmax; i++)
        {
            double x = x1 - f(x1) * ((x1 - x0) / (f(x1) - f(x0)));
            x0 = x1;
            x1 = x;

            if (almost_equal(f(x1), 0.0, precision) || (std::abs(x1 - x0) < tolerance))
            {
                return x1;
            }

        }
        throw std::runtime_error("Maximum number of iterations reached without finding root.");
    }

    /**
     * @brief Finds the root of a function using Brent's method.
     * 
     * @param f The function for which to find the root.
     * @param a The lower bound of the interval.
     * @param b The upper bound of the interval.
     * @param precision Number of decimal places to consider for convergence.
     * @param nmax Maximum number of iterations.
     * @return The root of the function.
     * @throws std::invalid_argument if the function values at the interval endpoints have the same sign.
     * @throws std::runtime_error if the root is not found within the maximum number of iterations.
     */
    double brent(real_func f, double a, double b, unsigned int precision, int nmax)
    {
        if (std::signbit(f(a)) == std::signbit(f(b)))
        {
            throw std::invalid_argument("Function values at the interval endpoints must have opposite signs.");
        }

        double tolerance = std::pow(10, -static_cast<int>(precision));

        if (f(a) < f(b))
        {
            std::swap(a, b);
        }

        double c = a;
        bool mflag = true;
        double s;
        double d;

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
                (mflag && (std::abs(b - c) < tolerance)) ||
                (!mflag && (std::abs(c - d) < tolerance))
                )
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

            if (almost_equal(f(b), 0.0, precision) || (std::abs(b - a) < tolerance))
            {
                return b;
            }

        }
        throw std::runtime_error("Maximum number of iterations reached without finding root.");

    }

}

#endif // ROOT_SCALAR_HPP
