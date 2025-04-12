#include "RootScalar.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <stdexcept>

// Define a test fixture for error tests
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

// Define a test fixture for valid tests
template <typename T>
class RootScalarValidTest : public ::testing::Test {
protected:
    // Example function to test root-finding methods
    static T test_function(T x) {
        return x * x - T(4); // Root at x = ±2
    }

    // Example function with a single root
    static T linear_function(T x) {
        return x - T(3); // Root at x = 3
    }

    // Example function with a more complex root
    static T cubic_function(T x) {
        return x * x * x - T(8); // Root at x = 2
    }
};

// Define the types to test
using TestTypes = ::testing::Types<float, double, long double>;
TYPED_TEST_SUITE(RootScalarErrorTest, TestTypes);
TYPED_TEST_SUITE(RootScalarValidTest, TestTypes);

// Error Tests
TYPED_TEST(RootScalarErrorTest, NewtonDerivativeZero) {
    using T = TypeParam;
    EXPECT_THROW(
        RootScalar::newton(this->flat_function, T(1), 6, 100),
        std::runtime_error
    );
}

TYPED_TEST(RootScalarErrorTest, BisectionInvalidInterval) {
    using T = TypeParam;
    EXPECT_THROW(
        RootScalar::bisection(this->test_function, T(3), T(4), 6, 100),
        std::invalid_argument
    );
}

TYPED_TEST(RootScalarErrorTest, SecantDivisionByZero) {
    using T = TypeParam;
    auto constant_function = [](T) { return T(1); }; // Constant function
    EXPECT_THROW(
        RootScalar::secant(constant_function, T(1), T(1), 6, 100),
        std::runtime_error
    );
}

TYPED_TEST(RootScalarErrorTest, BrentInvalidInterval) {
    using T = TypeParam;
    EXPECT_THROW(
        RootScalar::brent(this->test_function, T(3), T(4), 6, 100),
        std::invalid_argument
    );
}

TYPED_TEST(RootScalarErrorTest, NewtonMaxIterations) {
    using T = TypeParam;
    auto difficult_function = [](T x) { return std::exp(x) - T(100); }; // Root is far from initial guess
    EXPECT_THROW(
        RootScalar::newton(difficult_function, T(0), 6, 5), // Insufficient iterations
        std::runtime_error
    );
}

TYPED_TEST(RootScalarErrorTest, BisectionMaxIterations) {
    using T = TypeParam;
    EXPECT_THROW(
        RootScalar::bisection(this->test_function, T(0), T(3), 6, 1), // Insufficient iterations
        std::runtime_error
    );
}

TYPED_TEST(RootScalarErrorTest, SecantMaxIterations) {
    using T = TypeParam;
    auto difficult_function = [](T x) { return std::exp(x) - T(100); }; // Root is far from initial guesses
    EXPECT_THROW(
        RootScalar::secant(difficult_function, T(0), T(1), 6, 5), // Insufficient iterations
        std::runtime_error
    );
}

TYPED_TEST(RootScalarErrorTest, BrentMaxIterations) {
    using T = TypeParam;
    EXPECT_THROW(
        RootScalar::brent(this->test_function, T(0), T(3), 6, 1), // Insufficient iterations
        std::runtime_error
    );
}

// Valid Tests
TYPED_TEST(RootScalarValidTest, Newton) {
    using T = TypeParam;
    T root = RootScalar::newton(this->test_function, T(1), 6, 100);
    EXPECT_NEAR(root, T(2), T(1e-6));
}

TYPED_TEST(RootScalarValidTest, Bisection) {
    using T = TypeParam;
    T root = RootScalar::bisection(this->test_function, T(0), T(3), 6, 100);
    EXPECT_NEAR(root, T(2), T(1e-6));
}

TYPED_TEST(RootScalarValidTest, Secant) {
    using T = TypeParam;
    T root = RootScalar::secant(this->test_function, T(1), T(3), 6, 100);
    EXPECT_NEAR(root, T(2), T(1e-6));
}

TYPED_TEST(RootScalarValidTest, Brent) {
    using T = TypeParam;
    T root = RootScalar::brent(this->test_function, T(0), T(3), 6, 100);
    EXPECT_NEAR(root, T(2), T(1e-6));
}

TYPED_TEST(RootScalarValidTest, NewtonLinear) {
    using T = TypeParam;
    T root = RootScalar::newton(this->linear_function, T(0), 6, 100);
    EXPECT_NEAR(root, T(3), T(1e-6));
}

TYPED_TEST(RootScalarValidTest, BisectionLinear) {
    using T = TypeParam;
    T root = RootScalar::bisection(this->linear_function, T(0), T(5), 6, 100);
    EXPECT_NEAR(root, T(3), T(1e-6));
}

TYPED_TEST(RootScalarValidTest, SecantLinear) {
    using T = TypeParam;
    T root = RootScalar::secant(this->linear_function, T(0), T(5), 6, 100);
    EXPECT_NEAR(root, T(3), T(1e-6));
}

TYPED_TEST(RootScalarValidTest, BrentLinear) {
    using T = TypeParam;
    T root = RootScalar::brent(this->linear_function, T(0), T(5), 6, 100);
    EXPECT_NEAR(root, T(3), T(1e-6));
}

TYPED_TEST(RootScalarValidTest, NewtonCubic) {
    using T = TypeParam;
    T root = RootScalar::newton(this->cubic_function, T(1), 6, 100);
    EXPECT_NEAR(root, T(2), T(1e-6));
}

TYPED_TEST(RootScalarValidTest, BisectionCubic) {
    using T = TypeParam;
    T root = RootScalar::bisection(this->cubic_function, T(0), T(3), 6, 100);
    EXPECT_NEAR(root, T(2), T(1e-6));
}

TYPED_TEST(RootScalarValidTest, SecantCubic) {
    using T = TypeParam;
    T root = RootScalar::secant(this->cubic_function, T(1), T(3), 6, 100);
    EXPECT_NEAR(root, T(2), T(1e-6));
}

TYPED_TEST(RootScalarValidTest, BrentCubic) {
    using T = TypeParam;
    T root = RootScalar::brent(this->cubic_function, T(0), T(3), 6, 100);
    EXPECT_NEAR(root, T(2), T(1e-6));
}
