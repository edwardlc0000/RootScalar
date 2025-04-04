# RootScalar# RootScalar

RootScalar is a header-only C++ library for finding the roots of scalar functions using various numerical methods. It includes implementations of the Newton-Raphson method, Bisection method, Secant method, and Brent's method.

## Description

RootScalar provides a set of functions to find the roots of scalar functions. It is designed to be easy to use and integrate into your C++ projects. The library is header-only, which means you only need to include the header file in your project to start using it.

## Features

- Newton-Raphson method
- Bisection method
- Secant method
- Brent's method
- Numerical differentiation
- Precision control
- Error handling

## Installation

To use RootScalar in your project, simply include the `RootScalar.hpp` header file in your source code.

## Usage

### Newton-Raphson Method

```cpp

#include "RootScalar.hpp"
int main() 
{ 
    auto func = [](double x) { return x * x - 2; };
    
    double initial_guess = 1.0;
    unsigned int precision = 6;
    int max_iterations = 100;
    
    try
    {
        double root = RootScalar::newton(func, initial_guess, precision, max_iterations);
        std::cout << "Root: " << root << std::endl;
    } 
    catch (const std::runtime_error& e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
```

### Bisection Method

```cpp
#include "RootScalar.hpp"
int main()
{
    auto func = [](double x) { return x * x - 2; };

    double a = 0.0;
    double b = 2.0;
    unsigned int precision = 6;
    int max_iterations = 100;

    try 
    {
        double root = RootScalar::bisection(func, a, b, precision, max_iterations);
        std::cout << "Root: " << root << std::endl;
    }
    catch (const std::runtime_error& e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
```

### Secant Method

```cpp
int main() 
{
    auto func = [](double x) { return x * x - 2; };
    
    double x0 = 0.0;
    double x1 = 2.0;
    unsigned int precision = 6;
    int max_iterations = 100;
    
    try 
    {
        double root = RootScalar::secant(func, x0, x1, precision, max_iterations);
        std::cout << "Root: " << root << std::endl;
    } 
    catch (const std::runtime_error& e) 
    {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
```

### Brent's Method
```cpp
int main()
{
    auto func = [](double x) { return x * x - 2; };
    
    double a = 0.0; 
    double b = 2.0;
    unsigned int precision = 6;
    int max_iterations = 100;

    try 
    {
        double root = RootScalar::brent(func, a, b, precision, max_iterations);
        std::cout << "Root: " << root << std::endl;
    } 
    catch (const std::runtime_error& e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
    }

return 0;
}
```

## Contributing

Contributions are welcome! If you find a bug or have a feature request, please open an issue on GitHub. If you would like to contribute code, please fork the repository and submit a pull request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE.txt) file for details.

## Acknowledgements

- [Google Test](https://github.com/google/googletest) for the testing framework.
- [CMake](https://cmake.org/) for the build system.
- [Ninja](https://ninja-build.org/) for the build tool.
