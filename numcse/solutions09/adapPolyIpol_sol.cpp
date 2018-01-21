#include <iostream>
#include <cmath>
#include <mgl2/mgl.h>
#include "newtonIpol_sol.hpp"

// Adaptive polynomial interpolation; uses greedy algorithm
template <class Function>
Eigen::VectorXd adapPolyIpol(const Function& f,
                      const double a, const double b,
                      const double tol, const unsigned N,
                      Eigen::VectorXd& adaptive_nodes) {
    
	// Generate sampling points and evaluate $f$ there
	Eigen::VectorXd sampling_points = Eigen::VectorXd::LinSpaced(N, a, b);
	Eigen::VectorXd fvals_at_sampling_points = sampling_points.unaryExpr(f);

	// Approximate $\max |f(x)|$
	double fmax = fvals_at_sampling_points.cwiseAbs().maxCoeff();

	// Adaptive mesh (initial node)
	std::vector<double> t { (a+b)/2. }; // Set of interpolation nodes
	std::vector<double> y { static_cast<double>(f((a+b)/2.)) }; // Values at nodes
	std::vector<double> errors; // Error at nodes

	for (int i = 0; i < N; ++i) {
		// *** Step 1: interpolate with current nodes
		//   need to convert std::vector to
		//   Eigen::VectorXd to use the function interpoyval
		Eigen::Map<Eigen::VectorXd> te(t.data(), t.size());
		Eigen::Map<Eigen::VectorXd> ye(y.data(), y.size());
        
		Eigen::VectorXd intpolyvals_at_sampling_points;
		intPolyEval(te, ye, sampling_points, intpolyvals_at_sampling_points);

		// *** Step 2: find node where error is the largest
		Eigen::VectorXd err = (fvals_at_sampling_points - intpolyvals_at_sampling_points).cwiseAbs();
		double max = 0; int idx = 0;
		max = err.maxCoeff(&idx); // see Eigen "Visitor"

		// Step 3: check termination criteria, return results
		if (max < tol * fmax) {
			adaptive_nodes = te;
			Eigen::Map<Eigen::VectorXd> errVsStep(errors.data(), errors.size());
			return errVsStep;
		}

		// Step 4: add this node to our set of nodes and save error
		errors.push_back(max);
		t.push_back(sampling_points(idx));
		y.push_back(fvals_at_sampling_points(idx));
	}
	std::cerr << "Desired accuracy could not be reached." << std::endl;
	
	adaptive_nodes = sampling_points; // return all sampling points
}


int main() {
  
  // Test the interpolant
  /*VectorXd t(3), y(3), coeffs;
  t << -1, 0, 1;
  y << 2, -4, 6;
  coeffs = divDiff(t, y);
  std::cout << coeffs.transpose() << std::endl; // Expected output: [2 -6 8]
  */

  // Declare test functions
	auto f1 = [](double t) { return std::sin(std::exp(2*t)); };
  auto f2 = [](double t) { return std::sqrt(t)/(1 + 16*t*t); };

  // Test interval
  const double a = 0, b = 1;

  // Get interpolation nodes and print runtimes
  const unsigned N = 1000; // no. of sampling points
  const double tol = 1e-6; // tolerance
  Eigen::VectorXd tf1, tf2, // nodes for f1 resp. f2
                  ef1, ef2; // errors for f1 resp. f2

  ef1 = adapPolyIpol(f1, a, b, tol, N, tf1);
	ef2 = adapPolyIpol(f2, a, b, tol, N, tf2);

	Eigen::VectorXd n1_ = Eigen::VectorXd::LinSpaced(ef1.size(), 0, ef1.size());
	Eigen::VectorXd n2_ = Eigen::VectorXd::LinSpaced(ef2.size(), 0, ef2.size());

  // Plot
  mglData data1, data2, n1, n2;
  data1.Link(ef1.data(), ef1.size());
  data2.Link(ef2.data(), ef2.size());
  n1.Link(n1_.data(), ef1.size());
  n2.Link(n2_.data(), ef2.size());
  	
  mglGraph *gr = new mglGraph;
  gr->Title("Error vs Step");
  gr->SetRanges(0,150,1e-6,1e+0);  gr->SetFunc("x","lg(y)");
  gr->Axis();
  gr->Plot(n1, data1,"b+"); gr->AddLegend("f_1(t) = sin(e^{2t})","b+");
  gr->Plot(n2, data2,"rs"); gr->AddLegend("f_2(t) = \\sqrt{t}/(1 + 16t^2)","rs");
  gr->Label('x',"No. of nodes",0);
  gr->Label('y', "Error",0);
  gr->Legend();
	gr->WriteFrame("adaptive.eps");

  return 0;
}
