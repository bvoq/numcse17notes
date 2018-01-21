#include <iostream>
#include <set>
#include <cmath>
#include <mgl2/mgl.h>
#include "newtonIpol.hpp"

using namespace Eigen;
using namespace std;

// TODO: Implement the adaptive strategy using greedy algorithm
template <class Function>
void adaptivepolyintp(const Function& f, const double a, const double b,
                      const double tol, const unsigned N,
                      Eigen::VectorXd& adaptive_nodes,
                      Eigen::VectorXd& error_vs_step_no) {
    // Generate sampling points and evaluate $f$ there
    Eigen::VectorXd sampling_points = Eigen::VectorXd::LinSpaced(N, a, b);
    Eigen::VectorXd fvals_at_sampling_points = sampling_points.unaryExpr(f);
    
    set<int> usedSampleIndex = {(int)N/2};
    vector<double> tvec;
    vector<double> fvalvec;
    
    for(;;) {
        //interpolate T_n
        double max = -100000; int maxindex = 0;

        for(int i = 0; i < N; ++i) {
            //If not yet in the used samples:
            if(usedSampleIndex.count(i) == 0) {
                usedSampleIndex.insert(i); //insert the point
                tvec.push_back(sampling_points(i)); fvalvec.push_back(fvals_at_sampling_points(i));
                Eigen::Map<Eigen::VectorXd> teigen (tvec.data(), tvec.size());
                Eigen::Map<Eigen::VectorXd> fvaleigen (fvalvec.data(), fvalvec.size());

                VectorXd interpolcoeff = divDiff(teigen, fvaleigen);
                VectorXd evalvec (sampling_points(i));
                intPolyEval(teigen, interpolcoeff, evalvec, outevalvec);
                
                usedSampleIndex.remove(i); //exclude the point
            }
        }
    }
    

}


int main() {
  
  // Test the interpolant
  VectorXd t(3), y(3), coeffs;
  t << -1, 0, 1;
  y << 2, -4, 6;
  coeffs = divDiff(t, y);
  std::cout << coeffs.transpose() << std::endl; // Expected output: [2 -6 8]
    
    
  return 0;
}
