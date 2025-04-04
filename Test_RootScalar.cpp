#include "gtest/gtest.h"
#include "RootScalar.hpp"

//Test case for the Newton-Raphson method for valid function

TEST(RootScalarTest, NewtonValidFunction)
{
	//Function that converges easily: f(x) = x^2 -2
	auto func = [](double x) { return x * x - 2; };

	double initial_guess = 1.0;
	unsigned int precision = 6;
	int max_iterations = 100;

	double root = RootScalar::newton(func, initial_guess, precision, max_iterations);
	EXPECT_NEAR(root, std::sqrt(2), 1e-6);
}

// Test case for the Newton-Raphson method for a valid multivariate function with a fixed variable
TEST(RootScalarTest, NewtonMultivariateFunction)
{
	// Multivariate function: f(x, y) = x^2 + y^2 - 4
	// Fix y = 1 and find the root for x
	double y_fixed = 1.0;
	auto func = [y_fixed](double x) { return x * x + y_fixed * y_fixed - 4; };

	double initial_guess = 1.0;
	unsigned int precision = 6;
	int max_iterations = 100;

	double root = RootScalar::newton(func, initial_guess, precision, max_iterations);
	EXPECT_NEAR(root, std::sqrt(3), 1e-6);
}

// Test case for the Newton-Raphson for nonconvergence

TEST(RootScalarTest, NewtonMaxIterations)
{
	// Function that does not converge easily: f(x) = x^3 - 2x + 2
	auto func = [](double x) { return x * x * x - 2 * x + 2; };

	double initial_guess = 0.0;
	unsigned int precision = 6;
	int max_iterations = 10;

	EXPECT_THROW({
		RootScalar::newton(func, initial_guess, precision, max_iterations);
		}, std::runtime_error);
}

// Test case for the Bisection method for a valid function

TEST(RootScalarTest, BisectionValidFunction)
{
	// Function with a known root: f(x) = x^3 - x - 2
	auto func = [](double x) { return x * x * x - x - 2; };

	double a = 1.0;
	double b = 2.0;
	unsigned int precision = 6;
	int max_iterations = 100;

	double root = RootScalar::bisection(func, a, b, precision, max_iterations);
	EXPECT_NEAR(root, 1.5213797, 1e-6);
}

//Test case for the Bisection method for a valid multivariate function with a fixed variable

TEST(RootScalarTest, BisectionMultivariateFunction)
{
	// Multivariate function: f(x, y) = x^2 + y^2 - 4
	// Fix y = 1 and find the root for x
	double y_fixed = 1.0;
	auto func = [y_fixed](double x) { return x * x + y_fixed * y_fixed - 4; };

	double a = 1.0;
	double b = 2.0;
	unsigned int precision = 6;
	int max_iterations = 100;

	double root = RootScalar::bisection(func, a, b, precision, max_iterations);
	EXPECT_NEAR(root, std::sqrt(3), 1e-6);
}

// Test case for the Bisection method when function values do not have opposite signs

TEST(RootScalarTest, BisectionSameSignEndpoints)
{
	// Function with the same sign at both endpoints: f(x) = x^2 + 1
	auto func = [](double x) { return x * x + 1; };

	double a = -1.0;
	double b = 1.0;
	unsigned int precision = 6;
	int max_iterations = 100;

	EXPECT_THROW({
		RootScalar::bisection(func, a, b, precision, max_iterations);
		}, std::invalid_argument);
}

// Test case for the Bisection method when maximum number of iterations is reached without finding root

TEST(RootScalarTest, BisectionMaxIterations)
{
	// Function with a known root: f(x) = x^3 - x - 2
	auto func = [](double x) { return x * x * x - x - 2; };

	double a = 1.0;
	double b = 2.0;
	unsigned int precision = 6;
	int max_iterations = 2;

	EXPECT_THROW({
		RootScalar::bisection(func, a, b, precision, max_iterations);
		}, std::runtime_error);
}

// Test case for Brent's method for a valid function

TEST(RootScalarTest, BrentValidFunction)
{
	// Function with a known root: f(x) = x^3 - x - 2
	auto func = [](double x) { return x * x * x - x - 2; };

	double a = 1.0;
	double b = 2.0;
	unsigned int precision = 6;
	int max_iterations = 100;

	double root = RootScalar::brent(func, a, b, precision, max_iterations);
	EXPECT_NEAR(root, 1.5213797, 1e-6);
}

// Test case for Brent's method when function values at the interval endpoints do not have opposite signs

TEST(RootScalarTest, BrentSameSignEndpoints)
{
	// Function with the same sign at both endpoints: f(x) = x^2 + 1
	auto func = [](double x) { return x * x + 1; };

	double a = -1.0;
	double b = 1.0;
	unsigned int precision = 6;
	int max_iterations = 100;

	EXPECT_THROW({
		RootScalar::brent(func, a, b, precision, max_iterations);
		}, std::invalid_argument);
}

// Test case for Brent's method when maximum number of iterations is reached without finding root

TEST(RootScalarTest, BrentMaxIterations)
{
	// Function with a root: f(x) = x^3 - x - 2
	auto func = [](double x) { return x * x * x - x - 2; };

	double a = 1.0;
	double b = 2.0;
	unsigned int precision = 6;
	int max_iterations = 1; // Set to a low number to force max iterations

	EXPECT_THROW({
		RootScalar::brent(func, a, b, precision, max_iterations);
		}, std::runtime_error);
}

int main(int argc, char** argv) 
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}