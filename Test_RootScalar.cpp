#include "RootScalar.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <stdexcept>

// Define a test fixture for typed tests
template <typename T>
class RootScalarErrorTest : public ::testing::Test {
protected:
    // Example function to test root-finding methods
    static T test_function(T x) {
        return x * x - T(4); // Root at x = ±2
    }

    // Example function with a flat derivative (to trigger errors in Newton's method)
    static T flat_function(T x) {
        return T(1); // Derivative is always 0
    }
};

// Define the types to test
using TestTypes = ::testing::Types<float, double, long double>;
TYPED_TEST_SUITE(RootScalarErrorTest, TestTypes);

// Test Newton-Raphson method for derivative being zero
TYPED_TEST(RootScalarErrorTest, NewtonDerivativeZero) {
    using T = TypeParam;
    EXPECT_THROW(
        RootScalar::newton(this->flat_function, T(1), 6, 100),
        std::runtime_error
    );
}

// Test Bisection method for invalid interval
TYPED_TEST(RootScalarErrorTest, BisectionInvalidInterval) {
    using T = TypeParam;
    EXPECT_THROW(
        RootScalar::bisection(this->test_function, T(3), T(4), 6, 100),
        std::invalid_argument
    );
}

// Test Secant method for division by zero
TYPED_TEST(RootScalarErrorTest, SecantDivisionByZero) {
    using T = TypeParam;
    auto constant_function = [](T) { return T(1); }; // Constant function
    EXPECT_THROW(
        RootScalar::secant(constant_function, T(1), T(1), 6, 100),
        std::runtime_error
    );
}

// Test Brent's method for invalid interval
TYPED_TEST(RootScalarErrorTest, BrentInvalidInterval) {
    using T = TypeParam;
    EXPECT_THROW(
        RootScalar::brent(this->test_function, T(3), T(4), 6, 100),
        std::invalid_argument
    );
}

// Test Newton-Raphson method for exceeding maximum iterations
TYPED_TEST(RootScalarErrorTest, NewtonMaxIterations) {
    using T = TypeParam;
    auto difficult_function = [](T x) { return std::exp(x) - T(100); }; // Root is far from initial guess
    EXPECT_THROW(
        RootScalar::newton(difficult_function, T(0), 6, 5), // Insufficient iterations
        std::runtime_error
    );
}

// Test Bisection method for exceeding maximum iterations
TYPED_TEST(RootScalarErrorTest, BisectionMaxIterations) {
    using T = TypeParam;
    EXPECT_THROW(
        RootScalar::bisection(this->test_function, T(0), T(3), 6, 1), // Insufficient iterations
        std::runtime_error
    );
}

// Test Secant method for exceeding maximum iterations
TYPED_TEST(RootScalarErrorTest, SecantMaxIterations) {
    using T = TypeParam;
    auto difficult_function = [](T x) { return std::exp(x) - T(100); }; // Root is far from initial guesses
    EXPECT_THROW(
        RootScalar::secant(difficult_function, T(0), T(1), 6, 2), // Insufficient iterations
        std::runtime_error
    );
}

// Test Brent's method for exceeding maximum iterations
TYPED_TEST(RootScalarErrorTest, BrentMaxIterations) {
    using T = TypeParam;
    EXPECT_THROW(
        RootScalar::brent(this->test_function, T(0), T(3), 6, 1), // Insufficient iterations
        std::runtime_error
    );
}
