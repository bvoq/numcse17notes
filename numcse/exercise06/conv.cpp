// run with: g++ -std=gnu++11 -I /usr/include/eigen3/ -g conv.cpp -lmgl 
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <unsupported/Eigen/FFT>
#include <mgl2/mgl.h>

using namespace Eigen;
using namespace std;
#define MOD(a,b) (((a%b)+b)%b)

VectorXd PointA(const VectorXd &x, const VectorXd &y) {
    VectorXd out = Vector::Zero(x.size()+y.size()-1);
    for(int k = 0; k < x.size()+y.size()-1; ++k) for(int j = 0; k-j >= 0 && j < x.size(); ++j) out(k) += x(j)* (k-j>=y.size() ? 0 : y(k-j));
    
	return out;
}

VectorXd PointB(const VectorXd &x, const VectorXd &y) {
    int appendYZeros = x.size() - y.size(); //pretty sure this is wrong in the exercise sheet.
    VectorXd yext = VectorXd::Zero(x.size());
    yext.head(y.size()) = y;
    
    Eigen::FFT<double> fft;
    VectorXcd tmp = fft.fwd((VectorXcd) x).cwiseProduct(fft.fwd((VectorXcd) yext));
    VectorXd out = fft.inv(tmp).real();
    return out;
}

VectorXd PointC(const VectorXd &x, const VectorXd &y) {
    VectorXcd xext (x.size() + y.size() - 1);
    VectorXcd yext (x.size() + y.size() - 1);
    xext << x, VectorXd::Zero(y.size()-1);
    yext << y, VectorXd::Zero(x.size()-1);
    Eigen::FFT<double> fft;
    VectorXcd tmp = fft.fwd((VectorXcd) xext).cwiseProduct(fft.fwd((VectorXcd) yext));
    VectorXd out = fft.inv(tmp).real();
    return out;
}

int main() {
	srand(73);
	//VectorXd input = VectorXd::Zero(100);
	VectorXd noise = VectorXd::Random(100); 
    //input += (noise.array() - 0.5).matrix() / 5;
    
    //Special sine input ascending
    VectorXd input (100);
    for(int i = 0; i < 100; ++i) input(i) = sin(i/50. * 3.1415926); //sine between [0,2*PI] sampled at 100 equidistant points.
    for(int i = 0; i < 100; ++i) input(i) += .25*rand()/RAND_MAX; //add noise
    
    //Goal: Get the derivative without smoothing.
    
    
	// define truncated gaussian filter
	/*
    int N = 5;
	double a = .5;
	VectorXd gaussFilter(2*N+1);
	
    for (int j=0; j < 2*N+1; j++) {
		gaussFilter(j) = std::sqrt(a/3.1416) * std::exp(-a * ((j-N)) * (j-N));
	}
    */
    
    //Other filters
    //Derivative filter dx / dy = (f(x+dt) - f(x-dt))/2*dt
    float dy = 2.*3.141592654/100.; //we sample at a rate of 2*pi/100 per unit
    VectorXd derivFilter(3);
    derivFilter(0) = -1. / (2.*dy);
    derivFilter(1) = 0;
    derivFilter(2) = 1. / (2.*dy);
    //Now we need to reverse the filter to do a convolution instead of a correlation.
    for(int i = 0; i < derivFilter.size()/2; ++i) swap(derivFilter(i),derivFilter(derivFilter.size()-i-1));
    VectorXd gaussFilter = derivFilter;
    
    

	// Plotting
	mglGraph gr, grA, grB, grC, grGauss;
	mglData dat, datA, datB, datC, datGauss;
	dat.Link(input.data(), input.size());
	gr.Plot(dat);
	gr.WriteFrame("/Users/kdkdk/Desktop/ETH/numerics/numcse/numcse/exercise06/noisySignal.png");
	auto convA = PointA(input, gaussFilter);
  	datA.Link(convA.data(), convA.size());
	gr.Plot(datA);
	grA.Plot(datA);
    grA.WriteFrame("/Users/kdkdk/Desktop/ETH/numerics/numcse/numcse/exercise06/filteredA.png");
    auto convB = PointB(input, gaussFilter);
  	datB.Link(convB.data(), convB.size());
	gr.Plot(datB);
	grB.Plot(datB);
	grB.WriteFrame("/Users/kdkdk/Desktop/ETH/numerics/numcse/numcse/exercise06/filteredB.png");
	auto convC = PointC(input, gaussFilter);
  	datC.Link(convC.data(), convC.size());
	gr.Plot(datC);
	grC.Plot(datC);
	grC.WriteFrame("/Users/kdkdk/Desktop/ETH/numerics/numcse/numcse/exercise06/filteredC.png");
	gr.WriteFrame("/Users/kdkdk/Desktop/ETH/numerics/numcse/numcse/exercise06/noisy+filters.png");
    //plot filter:
    datGauss.Link(gaussFilter.data(), gaussFilter.size());
    grGauss.Plot(datGauss);
    grGauss.WriteFrame("/Users/kdkdk/Desktop/ETH/numerics/numcse/numcse/exercise06/gaussFilter.png");
    

}
