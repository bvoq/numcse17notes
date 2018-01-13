// run with: g++ -std=gnu++11 -I /usr/include/eigen3/ conv_sol.cpp -lmgl
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <unsupported/Eigen/FFT>
#include <mgl2/mgl.h>

using namespace Eigen;

VectorXd PointA(const VectorXd &x, const VectorXd &y) {
	// returns discrete convolution of x and y (straightforward implementation)
	int m = x.size();
	int n = y.size();
	int L = m+n-1;
	VectorXd out = VectorXd::Zero(L);

	for (int i=0; i<L; i++) {
		for (int k=std::max(0,i+1-m); k<std::min(i+1,n); k++) {
			out(i) += x(i-k) * y(k);
		}
	}

	return out;
}

VectorXd PointB(const VectorXd &x, const VectorXd &y) {
	// returns an overlapping convolution of x and y
	int m = x.size();
	int n = y.size();
	VectorXd ym = VectorXd::Zero(m);
	ym.head(n) = y;

	FFT<double> fft;
	VectorXcd a = fft.fwd(x);
	VectorXcd b = fft.fwd(ym);
	return fft.inv(a.cwiseProduct(b).eval()).real();
}

VectorXd PointC(const VectorXd &x, const VectorXd &y) {
	// returns discrete convolution of x and y (fft implementation)
	int m = x.size();
	int n = y.size();
	int L = m+n-1;
	VectorXd yL = VectorXd::Zero(L);
	VectorXd xL = VectorXd::Zero(L);
	xL.head(m) = x;
	yL.head(n) = y;

	FFT<double> fft;
	VectorXcd a = fft.fwd(xL);
	VectorXcd b = fft.fwd(yL);
	return fft.inv(a.cwiseProduct(b).eval()).real();
}

int main() {
	srand(73);
	VectorXd input = VectorXd::Zero(100);
	VectorXd noise = VectorXd::Random(100); 
	input += (noise.array() - 0.5).matrix() / 5;

	// define truncated gaussian filter
	int N = 5;
	double a = .5;
	VectorXd gaussFilter(2*N+1);
	for (int j=0; j < 2*N+1; j++) {
		gaussFilter(j) = std::sqrt(a/3.1416) * std::exp(-a * ((j-N)) * (j-N));
	}

	// Plotting
	mglGraph gr, grA, grB, grC;
	mglData dat, datA, datB, datC;
	dat.Link(input.data(), input.size());
	gr.Plot(dat);
	gr.WriteFrame("/Users/kdkdk/Desktop/ETH/numerics/numcse/numcse/exercise06/solnoisySignal.png");
	auto convA = PointA(input, gaussFilter);
  	datA.Link(convA.data(), convA.size());
	gr.Plot(datA);
	grA.Plot(datA);
	grA.WriteFrame("/Users/kdkdk/Desktop/ETH/numerics/numcse/numcse/exercise06/solfilteredA.png");
	auto convB = PointB(input, gaussFilter);
  	datB.Link(convB.data(), convB.size());
	gr.Plot(datB);
	grB.Plot(datB);
	grB.WriteFrame("/Users/kdkdk/Desktop/ETH/numerics/numcse/numcse/exercise06/solfilteredB.png");
	auto convC = PointC(input, gaussFilter);
  	datC.Link(convC.data(), convC.size());
	gr.Plot(datC);
	grC.Plot(datC);
	grC.WriteFrame("/Users/kdkdk/Desktop/ETH/numerics/numcse/numcse/exercise06/solfilteredC.png");
	gr.WriteFrame("/Users/kdkdk/Desktop/ETH/numerics/numcse/numcse/exercise06/solnoisy+filters.png");
}
