#include <utility>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

#include <Eigen/Dense>

using namespace Eigen;

double norm(double x) { return std::abs(x); }
double norm(const VectorXd & x) { return x.norm(); }


/*!
 *! \brief Implements a single step of the modified newton
 *! \tparam Scalar type of argument to function f: such as double, etc...
 *! \tparam Function type for the function f, likely a lambda function
 *! \tparam Jacobian type for the jacobian df, likely a lambda function
 *! \param[in] x previous value to use in step, also initial guess if needed
 *! \param[in] f function handle for f(x) = 0
 *! \param[in] df function handle for jacobian df of f
 *! \return next step x_{k+1} of modified Newton
 */
template <typename Scalar, class Function, class Jacobian>
Scalar mod_newt_step_scalar(const Scalar& x, 
                            const Function& f, 
                            const Jacobian& df) {
    Scalar y = x + f(x) / df(x);
    return y - f(y) / df(x);
}


/*!
 *! \brief Performs many steps of an iteration and terminate when convergence reached
 *! or maximum number of iterations has been reached.
 *! \tparam StepFunction type for the step function handle
 *! \tparam Vector argument type passed to the iteration function
 *! \tparam ErrorFunction type for the error function computing error of the method
 *! \param[in] step Function implementing the step x_{k+1} = step(x_{k}), signature Vector(const Vector&)
 *! \param[in,out] x initial data (as input) and final iteration (as output)
 *! \param[in] errf function implementing the norm of the error (errf(x)) for termination condition
 *! \param[in] eps tolerance to break iterations when res(x) < eps
 *! \param[in] max_itr maximal number of iterations
 */
template <class StepFunction, class Vector, class ErrorFunction>
bool sample_nonlinear_solver(const StepFunction& step,
                             Vector & x,
                             const ErrorFunction& errf,
                             double eps = 1e-8, int max_itr = 100) {
    // Temporary where to store new step
    Vector x_new = x;
    double r = 1;
    
    for(int itr = 0; itr < max_itr; ++itr) {
        
        r = errf(x);// Compute error (or residual)
        std::cout<<"[Step " <<itr<< "] Error: "<< r <<std::endl;
        
        x_new = step(x);// Advance to next step, $x_{new}$ becomes $x_{k+1}$
        
        // Termination criteria
        if (r < eps * norm(x)) {
            std::cout << "[CONVERGED] in " << itr << " it. due to err. err = " 
                      << r << " < " << eps << "." << std::endl;
            return true;
        }
        x = x_new;
    }
    
    // If max itr reached
    std::cout << "[NOT CONVERGED] due to MAX it. = " << max_itr 
              << " reached, err = " << r << "." << std::endl;
    return false;
}

/**
 *! \brief Solve a scalar non-linear eq. with the modified Newton.
 * Test the empirical order of convergence of the method.
 */
void mod_newt_ord() {
    // Setting up values, functions and jacobian
    const double a = 0.123;
    auto f = [&a] (double x) { return atan(x) - a; };
    auto df = [] (double x) { return 1. / (x*x+1.); };
    
    const double x_ex = tan(a); // Exact solution
    
    // Store solution and error at each time step
    std::vector<double> sol;
    std::vector<double> err;
    
    // Compute error and push back to err, used in 
    // general\_nonlinear\_solver as breaking condition errf(x) < eps
    auto errf = [&err, x_ex] (double & x) { 
        double e =  std::abs(x - x_ex); 
        err.push_back(e); 
        return e;
    };
    
    // Perform convergence study with Modified newton for scalar
    std::cout << std::endl 
              << "*** Modified Newton method (scalar) ***"
              << std::endl << std::endl;
    std::cout << "Exact: " << x_ex << std::endl;
    
    // Initial guess and compute initial error
    double x_scalar = 5;
    sol.push_back(x_scalar);
    errf(x_scalar);
    
    // Lambda performing the next step, used to define a proper 
    // function handle to be passed to general\_nonlinear\_solver
    auto newt_scalar_step = [&sol, &f, &df] (double x) -> double { 
        double x_new = mod_newt_step_scalar(x, f, df);
        sol.push_back(x_new);
        return x_new;
    };
    
    // Actually perform the solution
    sample_nonlinear_solver(newt_scalar_step, x_scalar, errf);
    
    // Print solution (final)
    std::cout << std::endl 
              << "x^*_scalar = " << x_scalar 
              << std::endl;

    // Print table of solutions, errors and EOOC
    auto space = std::setw(15);
    std::cout << space << "sol." 
              << space << "err." 
              << space << "order" 
              << std::endl;
    for(unsigned i = 0; i < sol.size(); ++ i) {
        std::cout << space << sol.at(i) 
                  << space << err.at(i);
        if(i >= 3) {
            std::cout << space << (log(err.at(i)) - log(err.at(i-1))) / (log(err.at(i-1)) - log(err.at(i-2)));
        }
        std::cout << std::endl;
    }
}

int main() {
    mod_newt_ord();
}
