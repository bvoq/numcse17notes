#include <iomanip>
#include <iostream>
#include <cmath>

#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

#include "golubwelsh.hpp"

#define PI M_PI
#define PI_HALF M_PI_2

//! @brief Compute $\int_a^b f(x) dx \approx \sum w_i f(x_i)$ (with scaling of $w$ and $x$)
//! @tparam Function template type for function handle f (e.g.\ lambda function)
//! @param[in] f integrand
//! @param[in] w weights
//! @param[in] x nodes for interval $[-1,1]$
//! @param[in] a left boundary in $[a,b]$
//! @param[in] b right boundary in $[a,b]$
//! @return Approximation of integral $\int_a^b f(x) dx$
template <class Function>
double quad(const Function& f,
            const Eigen::VectorXd& w, const Eigen::VectorXd& x,
            double a, double b) {
    double I = 0;
    for(int i = 0; i < w.size(); ++i) {
        I += f( (x(i) + 1) * (b - a) / 2 + a ) * w(i);
    }
    return I * (b - a) / 2.;
}

//! @brief Compute $\int_{-\infty}^\infty f(x) dx$ using transformation $x = \cot(t)$
//! @tparam Function template type for function handle f (e.g.\ lambda function)
//! @param[in] n number of Gauss points
//! @param[in] f integrand
//! @return Approximation of integral $\int_{-\infty}^\infty f(x) dx$
template <class Function>
double quadinf(const int n, Function&& f) {
    Eigen::VectorXd w, x;

    // Compute nodes and weights of Gauss quadrature rule
    // using Golub-Welsh algorithm
    golubwelsh(n, w, x);

    //! NOTE: no function cot available in c++, need to resort to trigonometric identities
    //! Both lines below are valid, the first computes three trigonometric functions
    auto ftilde = [&f] (double x) { return f(std::cos(x)/std::sin(x)) / pow(std::sin(x),2); };
    /* auto ftilde = [&f] (double x) { double cot = std::tan(PI_HALF - x); return f(cot) * (1. + pow(cot,2)); }; */

    return quad(ftilde, w, x, 0, PI);
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
