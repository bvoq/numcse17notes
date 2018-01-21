#include <iomanip>
#include <iostream>
#include <cmath>

#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

#include "golubwelsh.hpp"

#define PI M_PI
#define PI_HALF M_PI_2

using namespace Eigen;

//! @brief Compute $\int_{-\infty}^\infty f(x) dx$ using transformation $x = \cot(t)$
//! @tparam Function template type for function handle f (e.g.\ lambda function)
//! @param[in] n number of Gauss points
//! @param[in] f integrand
//! @return Approximation of integral $\int_{-\infty}^\infty f(x) dx$


template <class Function>
double quadinf(const int n, Function&& f) {
    VectorXd w, x;
    
    //First create the substituted function g(u) = f(cot(u)) * 1/sin(u^2):
    auto g = [] (double u, Function&& f) { return f(cos(u)/sin(u)) * 1./(sin(u)*sin(u)); };

    // Compute nodes and weights of Gauss quadrature rule
    // using Golub-Welsh algorithm
    golubwelsh(n, w, x);

    double a = 0, b = PI;
    
    //Scale nodes and weights to interval 0,PI
    VectorXd xscaled(n), wscaled(n);
    for(int i = 0; i < n; ++i) {
        xscaled(i) = ((x(i) + 1)/2.)*(b-a)+a;
        wscaled(i) = (b-a)/2. * w(i);
    }
    
    double res = 0;
    //Evaluate the quadrature formula
    for(int i = 0; i < n; ++i) {
        res += g(xscaled(i),f) * wscaled(i);
    }
    
    return res;
}

int main() {
    // Number of max Gauss pts.
    const int N = 100;

    // Integrand and exact integral
    auto f = [] (double t) { return std::exp(-std::pow((t-1),2)); };

    // Exact value of integrand
    double I_ex = std::sqrt(PI);
    
    // NOTE: We observe exponential convergence
    int sep = 12;
    std::cout << std::setw(sep) << "Nodes"
              << std::setw(sep) << "Quadrature"
              << std::setw(sep) << "Exact"
              << std::setw(sep) << "Error"
              << std::endl;
    for(int n = 1; n <= N; ++n) {

        // Value of integrant approximated
        double QS = quadinf(n, f);

        std::cout << std::setw(sep) << n
                  << std::setw(sep) << QS
                  << std::setw(sep) << I_ex
                  << std::setw(sep) << std::abs(QS - I_ex)
                  << std::endl;
    }
    return 0;
}
