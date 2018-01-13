# include <iostream>
# include <cmath>
# include <vector>
# include <Eigen/Dense>
# include "gaussquad.hpp"

/*!
 * \brief integrate Compute integral given quadrature rule.
 * Compute the integral $\int_a^b f(x) dx \approx \sum_{i=1}^n w_i f(c_i)$
 * for given quadrature rule $\{(w_i, x_i)\}_{i=1}^n$
 * \tparam Function A function object with operator()(double)
 * \param qr A QuadRule object passing nodes and weights
 * \param f A function handle passing the integrand
 * \return Value of integral of $f$ using quadrature rule Q
 */
template <class Function>
double integrate(const QuadRule& qr, const Function& f) {
    double I = 0;
    for (unsigned i = 0; i < qr.weights.size(); ++i) {
        I += qr.weights(i) * f(qr.nodes(i));
    }
    return I;
}

/*!
 * \brief nonSmoothIntegrand Approximte the integral $\int_{-1}^1 \arcsin(x) f(x) dx$
 *
 * \tparam Function A function object with operator()(double)
 * \param fh Will pass the integrand
 * \param I_ex Exact value of integral
 * \return Value of integral
 */
template <class Function>
void nonSmoothIntegrand(const Function& fh, const double I_ex) {

    std::vector<double> evals, // used to save no.\ of quad nodes
            error; // used to save the error
    // Build integrand
    auto f = [fh](double x) {
        return std::asin(x)*fh(x);
    };

    const unsigned N = 50; // Max. no. of nodes

    for (unsigned n = 1; n <= N; ++n) {
        QuadRule qr;
        gaussquad(n, qr); // Create quadrature rule

        double I = integrate(qr, f); // Compute integral

        evals.push_back(n); // Save no. of quadrature nodes
        error.push_back(std::abs(I - I_ex)); // Save error
    }
    
    Eigen::Map<Eigen::VectorXd> numNodes(evals.data(), evals.size());
    Eigen::Map<Eigen::VectorXd> quadErr(error.data(), error.size());
    
    for (int i=0; i<evals.size(); i++)
			std::cout << numNodes(i) << "\t" << quadErr(i) << std::endl;    

}

/*!
 * \brief smoothIntegrand Approximte the integral $\int_{-1}^1 \arcsin(x) f(x) dx$
 * Ensures that convergenge is expoenential using appropriate transformation,
 * provided $f$ is a smooth function.
 * \tparam Function A function object with operator()(double)
 * \param fh Will pass the integrand
 * \param I_ex Exact value of integral
 * \return Value of integral
 */
template <class Function>
void smoothIntegrand(const Function& f, const double I_ex) {
    std::vector<double> evals, // Used to save no. of quad nodes
                        error; // Used to save the error
    // Transform integrand
    auto g = [f] (double x) {
        return x*f(std::sin(x))*std::cos(x);
    };

    const unsigned N = 50; // Max. no. of nodes

    for (unsigned n = 1; n <= N; ++n) {
        QuadRule qr;
        gaussquad(n, qr); // Obtain quadrature rule

        // Transform nodes and weights to new interval
        Eigen::VectorXd w = qr.weights * M_PI/2;
        Eigen::VectorXd c = ( -M_PI/2 + M_PI/2*(qr.nodes.array() + 1) ).matrix();

        // Evaluate $g$ at quadrature nodes $c$
        Eigen::VectorXd gc = c.unaryExpr(g);

        // Same as $I = \sum_{i=1}^n w_i g(c_i)$
        double I = w.dot(gc);

        evals.push_back(n); // Save no. of quadrature nodes
        error.push_back(std::abs(I - I_ex)); // Save error
    }
    
    Eigen::Map<Eigen::VectorXd> numNodes(evals.data(), evals.size());
    Eigen::Map<Eigen::VectorXd> quadErr(error.data(), error.size());
    
    for (int i=0; i<evals.size(); i++)
			std::cout << numNodes(i) << "\t" << quadErr(i) << std::endl;
    
}

int main() {
    // "Exact" value of integral
    const double I_ex = 0.870267525725852642;

    // $f(x) = \sinh x$
    std::function<double(double)> f = [](double x) {
        return std::sinh(x);
    };

    // PART 1
    nonSmoothIntegrand(f, I_ex);

    // PART 2
    //smoothIntegrand(f, I_ex);

    return 0;
}
