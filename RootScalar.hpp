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
    template <typename T>
    bool almost_equal(T a, T b, unsigned int precision)
    {
        return std::abs(a - b) < std::pow(10, -static_cast<int>(precision));
    }

    template <typename T, typename Function>
    T derivate(Function f, T x)
    {
        T h = std::sqrt(std::numeric_limits<T>::epsilon());
        return (f(x + h) - f(x - h)) / (2.0 * h);
    }

    template <typename T, typename Function>
    T newton(Function f, T x0, unsigned int precision, int nmax)
    {
        T x;
        T f_prime;
        T pertubation = std::sqrt(std::numeric_limits<T>::epsilon());
        T tolerance = std::pow(10, -static_cast<int>(precision));

        for (int i = 0; i < nmax; i++)
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

            if (almost_equal(f(x), static_cast<T>(0.0), precision) || (std::abs(x - x0) < tolerance))
            {
                return x;
            }

            x0 = x;
        }
        throw std::runtime_error("Maximum number of iterations reached without finding root.");
    }

    template <typename T, typename Function>
    T bisection(Function f, T a, T b, unsigned int precision, int nmax)
    {
        if (std::signbit(f(a)) == std::signbit(f(b)))
        {
            throw std::invalid_argument("Function values at the interval endpoints must have opposite signs.");
        }

        T tolerance = std::pow(10, -static_cast<int>(precision));
        T c;

        for (int i = 0; i < nmax; i++)
        {
            c = (a + b) / 2;
            if (almost_equal(f(c), static_cast<T>(0.0), precision) || (std::abs(b - a) < tolerance))
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

    template <typename T, typename Function>
    T secant(Function f, T x0, T x1, unsigned int precision, int nmax)
    {
        if (f(x1) == f(x0))
        {
            throw std::runtime_error("Division by zero encountered in secant method.");
        }

        T tolerance = std::pow(10, -static_cast<int>(precision));

        for (int i = 0; i < nmax; i++)
        {
            T x = x1 - f(x1) * ((x1 - x0) / (f(x1) - f(x0)));
            x0 = x1;
            x1 = x;

            if (almost_equal(f(x1), static_cast<T>(0.0), precision) || (std::abs(x1 - x0) < tolerance))
            {
                return x1;
            }
        }
        throw std::runtime_error("Maximum number of iterations reached without finding root.");
    }

    template <typename T, typename Function>
    T brent(Function f, T a, T b, unsigned int precision, int nmax)
    {
        if (std::signbit(f(a)) == std::signbit(f(b)))
        {
            throw std::invalid_argument("Function values at the interval endpoints must have opposite signs.");
        }

        T tolerance = std::pow(10, -static_cast<int>(precision));

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
                (mflag && (std::abs(b - c) < tolerance)) ||
                (!mflag && (std::abs(c - d) < tolerance)))
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

            if (almost_equal(f(b), static_cast<T>(0.0), precision) || (std::abs(b - a) < tolerance))
            {
                return b;
            }
        }
        throw std::runtime_error("Maximum number of iterations reached without finding root.");
    }
}

#endif // ROOT_SCALAR_HPP