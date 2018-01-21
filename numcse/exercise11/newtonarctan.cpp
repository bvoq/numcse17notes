#include <cmath>
#include <iostream>
#include <limits>

#include <Eigen/Eigen>


using namespace Eigen;
using namespace std;

/* @brief Newton's method to approximate $x^{(0)}$
 * @param[in] x0_ Initial guess
 * @param[out] x0 Final estimation of $x^{(0)}$, given convergence of Newton's method
 */

double atandiff(double x) {
    return 1. / (1. + x*x);
}


//atan(x) has a root at x=0, therefore you can simply take an init close to 0.
//however the point of this exercise is to do it with damping.
template <class Function> double newton_arctannodamping(double x0_, const Function& f, const Function & fdiv) {
    //first use normal newton
    int steps = 100;
    double x = x0_;
    cout << "init: " << x << endl;
    
    for(unsigned i=1; i < steps+1; ++i) {
        x = x - f(x)/fdiv(x);
        cout << i << " = " << x << endl;
    }
    return x;
}

template <class Function> double newton_arctandamping(double x0_, const Function& f, const Function &fdiv) {
    //first use normal newton
    int steps = 100;
    double x = x0_;
    cout << "init: " << x << endl;

    double lambda = 1.; //initial damping factor (isn't halfed until needded)

    
    for(unsigned i=1; i < steps+1; ++i) {
        double xdelta = 0, xldelta = 0;
        lambda *= 2.;
        do {
            lambda /= 2.;
            xdelta  = 1./fdiv(x) * f(x); //factor according to
            xldelta = 1./fdiv(x) * f(x + lambda * xdelta);
        } while(!  (fabs(xldelta * lambda) <= (1. - lambda / 2.) * fabs(xdelta)));
        
        
        x = x - lambda * f(x)/fdiv(x);
        cout << i << " = " << x << "  l: " << lambda << endl;
    }
    return x;
}


int main() {
	// Initialization

    function<double(double)> f = [] (double x) { return atan(x) ; };
    function<double(double)> fdiv = [] (double x) { return 1. / (1. + x*x); };

    
    double xconv = 1.39174; // initial point for which the newton method without damping converge.
    double xnoconv = 1.39175; // initial point for which the newton method without damping does not converge.
    double x0_ = 2.; //also doesn't converge here.
    
    double res = newton_arctannodamping(xconv, f, fdiv);
    double resdamped = newton_arctandamping(xnoconv, f, fdiv);

    //Functions to figure out on which inputs of the newton iteration converges.
    function<double(double)> g = [] (double x) { return 2*x - (1.+x*x)*atan(x); };
    function<double(double)> gdiv = [] (double x) { return 2 - 2.*x*atan(x) - 1; };
    //which ironically use newton iteration to prove convergence:
    //double x1 = newton_arctandamping(x0_, g, gdiv);
}
