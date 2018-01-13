#include <cmath>
#include <iostream>
#include <limits>

#include <Eigen/Dense>


using namespace Eigen;

/* @brief Newton's method to approximate $x^{(0)}$
 * @param[in] x0_ Initial guess
 * @param[out] x0 Final estimation of $x^{(0)}$, given convergence of Newton's method
 */
double newton_arctan(double x0_) {

    double x0 = x0_;
    double upd = 1;
    double eps = std::numeric_limits<double>::epsilon();

    while(std::abs(upd) > eps) {

        double x1 = x0 - (2*x0-(1+x0*x0)*std::atan(x0)) / (1-2*x0*std::atan(x0));
        upd = (x1 - x0) / x1;
        x0 =  x1;
        
        std::cout << "x0 = " << x1 << ", accuracy = " << upd << std::endl;
    }

    return x0;
}

int main() {
	// Initialization
    double x0_ = 2; // Initial guess

  // Netwon's method
    double x0 = newton_arctan(x0_);
}
