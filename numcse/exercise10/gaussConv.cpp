# include <cmath>
# include <vector>
# include <Eigen/Dense>
# include "gaussquad.hpp"
# include <iomanip>
# include <mgl2/mgl.h>
using namespace Eigen;
using namespace std;
#include <iostream>
//#include <mgl2/wnd.h>
#include <mgl2/fltk.h>

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
    
    VectorXd plotx(50-1), ploty(50-1);
    
    for(int n = 1; n < 50; ++n) {
        QuadRule qr;
        gaussquad(n, qr); //nodes and weights of gauss quadrature on interval [-1,1
        double res=0;
        for(int j = 0; j < n; ++j) {
            res += std::asin(qr.nodes(j)) * fh(qr.nodes(j)) * qr.weights(j);
        }
        
        plotx(n-1) = n;
        ploty(n-1) = fabs(res - I_ex);
        cout << "n : " << n << " res: " << res  << " error: " << fabs(res - I_ex) << endl;
    }
    
    mglData datx, daty;
    datx.Link(plotx.data(), plotx.size());
    daty.Link(ploty.data(), ploty.size());
    mglGraph gr;
    gr.SetRanges(0,50,0,.0001);
    gr.Plot(datx, daty, "0");
    gr.Axis();
    //gr.Run();
    gr.WriteFrame("/Users/kdkdk/Desktop/ETH/numerics/numcse/numcse/exercise10/errorplotnonsmooth.png");
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
    VectorXd plotx(50-1), ploty(50-1);
    
    for(int n = 1; n < 50; ++n) {
        QuadRule qr;
        gaussquad(n, qr); //nodes and weights of gauss quadrature on interval [-1,1
        //scale nodes to [-pi/2, pi/2]
        for(int i = 0; i < qr.nodes.size(); ++i) qr.nodes(i) = qr.nodes(i) * M_PI/2.;
        for(int i = 0; i < qr.weights.size(); ++i) qr.weights(i) = qr.weights(i) * M_PI/2.;
        
        //cout << "n=" << n  << "   "; for(int i = 0; i < qr.nodes.size(); ++i) cout << qr.nodes(i) << " "; cout << endl;
        //n=8   -1.50842 -1.2514 -0.825504 -0.288138 0.288138 0.825504 1.2514 1.50842
        //n=8   -1.50842 -1.2514 -0.825504 -0.288138 0.288138 0.825504 1.2514 1.50842

        double res=0;
        for(int j = 0; j < qr.nodes.size(); ++j) {
            res += qr.nodes(j) * std::cos(qr.nodes(j)) * f(std::sin(qr.nodes(j))) * qr.weights(j);
        }
        //n=8   -1.50842 -1.2514 -0.825504 -0.288138 0.288138 0.825504 1.2514 1.50842
        plotx(n-1) = n;
        ploty(n-1) = fabs(res - I_ex);
        cout << "n : " << n << " res: " << res  << " error: " << fabs(res - I_ex) << endl;
    }
    
    mglData datx, daty;
    datx.Link(plotx.data(), plotx.size());
    daty.Link(ploty.data(), ploty.size());
    mglGraph gr;
    gr.SetRanges(0,50,0,.0001);
    gr.Plot(datx, daty, "0");
    gr.Axis();
    gr.WriteFrame("/Users/kdkdk/Desktop/ETH/numerics/numcse/numcse/exercise10/errorplotsmooth.png");

}

int main() {
    // "Exact" value of integral
    const double I_ex = 0.870267525725852642;

    // $f(x) = \sinh x$
    std::function<double(double)> f = [](double x) {
        return std::sinh(x);
    };
    
    std::cout << std::setprecision(16);

    // PART 1
    nonSmoothIntegrand(f, I_ex);

    // PART 2
    smoothIntegrand(f, I_ex);

    return 0;
}
