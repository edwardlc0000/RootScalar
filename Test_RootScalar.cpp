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

int main(int argc, char** argv) 
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}