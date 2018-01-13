#include <algorithm>
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <unsupported/Eigen/FFT>
#include <mgl2/mgl.h>
#include <mgl2/wnd.h>
#include <mgl2/fltk.h>

using namespace Eigen;
using namespace std;

VectorXd discreteconvolution(const VectorXd &x, const VectorXd &y) {
    VectorXd out = VectorXd::Zero(x.size()+y.size()-1);
    for(int k = 0; k < x.size()+y.size()-1; ++k) for(int j = 0; k-j >= 0 && j < x.size(); ++j) out(k) += x(j)* (k-j>=y.size() ? 0 : y(k-j));
    return out;
}


VectorXd lowpassfilter(const VectorXd & x, int freqCutoffRadius) {
    FFT<double> fft;
    VectorXcd y = fft.fwd(x);
    
    for(int i = y.size()/2 - freqCutoffRadius; i <= y.size()/2 + freqCutoffRadius; ++i) y(i) = 0;
    
    VectorXcd filtered = fft.inv(y);
    return filtered.real();
}

int main() {
    srand(73);
    //sine input
    VectorXd input (100);
    for(int i = 0; i < 100; ++i) input(i) = sin(i/50. * 3.1415926); //sine between [0,2*PI] sampled at 100 equidistant points.
    for(int i = 0; i < 100; ++i) input(i) += .25*rand()/RAND_MAX; //add noise
    
    //Exercise: Figure out a filter of size 2*N+1, which does the derivative and noise reduction on the data and apply it!
    int N = 9; //9 is just a random example, but it will be easier with 9 for sampling.
    VectorXd filter (2*N+1);
    
    //Solution
    
    //A gauss filter is not really the best noise reduction here, since the noise is not a continous function.
    double a = .25; //if a = 1, then the integral of the gauss function is 1 and therefore it would be appropriately scaled
    VectorXd gaussFilter(N);
    int mid = N/2;
    for (int j=0; j < N; j++) {
        gaussFilter(j) = std::sqrt(a/3.1416) * std::exp(-a * ((j-mid)) * (j-mid));
    }
    
    //Instead simply use a mean filter:
    VectorXd meanFilter(N);
    for(int i=0; i < N; ++i) meanFilter(i) = 1./N;
    
    
    VectorXd derivFilter = VectorXd::Zero(N);
    //Derivative filter dx / dy = (f(x+dt) - f(x-dt))/2*dt
    float dy = 2.*3.141592654/100.; //we sample at a rate of 2*pi/100 per unit
    derivFilter(mid-1) = -1. / (2.*dy);
    derivFilter(mid+1) = 1. / (2.*dy);
    
    //Now we need to reverse the filter to do a convolution instead of a correlation.
    for(int i = 0; i < derivFilter.size()/2; ++i) swap(derivFilter(i),derivFilter(derivFilter.size()-i-1));
    
    
    //Convolution (not correlation) is associative, so we can make one filter out of it before applying it on input:
    filter = discreteconvolution(meanFilter, derivFilter);
    //filter = discreteconvolution(gaussFilter, derivFilter);
    
    
    mglGraph gr, grA, grLow;
    mglData dat, datA, datLow;
    dat.Link(input.data(), input.size());
    gr.Plot(dat);
    gr.WriteFrame("/Users/kdkdk/Desktop/ETH/numerics/numcse/numcse/exercise06/inputSignal.png");
    VectorXd filteredsignal = discreteconvolution(input, filter);
    datA.Link(filteredsignal.data(), filteredsignal.size());
    gr.Plot(datA);
    grA.Plot(datA);
    grA.WriteFrame("/Users/kdkdk/Desktop/ETH/numerics/numcse/numcse/exercise06/filteredderivatedSignal.png");
    
    VectorXd lowpassfilteredsignal = lowpassfilter(input, 40);
    datLow.Link(lowpassfilteredsignal.data(), lowpassfilteredsignal.size() );
    grLow.Plot(datLow);
    grLow.WriteFrame("/Users/kdkdk/Desktop/ETH/numerics/numcse/numcse/exercise06/lowpassfilteredSignal.png");
    gr.Plot(datLow);
    
    gr.WriteFrame("/Users/kdkdk/Desktop/ETH/numerics/numcse/numcse/exercise06/inputSignal.png+filteredderivatedSignal+lowpass.png");

    
    mglFLTK grfltk;
    //1. set title
    //2. set range
    //3. set axis
    //4. plot data (IN THIS ORDER!)
    
    grfltk.Title("Filtered");
    grfltk.SetRanges(-2, 2, -2, 2); //Always set range BEFORE adding axis!
    grfltk.Axis();


    grfltk.Plot(datA);
    grfltk.Plot(dat);
    grfltk.Plot(datLow);
    grfltk.Run();
}


/*
int main() {
    srand(73);
    //sine input
    VectorXd input (100);
    for(int i = 0; i < 100; ++i) input(i) = sin(i/50. * 3.1415926); //sine between [0,2*PI] sampled at 100 equidistant points.
    for(int i = 0; i < 100; ++i) input(i) += .25*rand()/RAND_MAX; //add noise
    
    //Exercise: Figure out a filter of size 2*N+1, which does the derivative and noise reduction on the data and apply it!
    int N = 9; //9 is just a random example, but it will be easier with 9 for sampling.
    VectorXd filter (2*N+1);
    
    //Solution
    
    //A gauss filter is not really the best noise reduction here, since the noise is not a continous function.
    double a = .25; //if a = 1, then the integral of the gauss function is 1 and therefore it would be appropriately scaled
    VectorXd gaussFilter(N);
    int mid = N/2;
    for (int j=0; j < N; j++) {
        gaussFilter(j) = std::sqrt(a/3.1416) * std::exp(-a * ((j-mid)) * (j-mid));
    }
    
    //Instead simply use a mean filter:
    VectorXd meanFilter(N);
    for(int i=0; i < N; ++i) meanFilter(i) = 1./N;
    
    
    VectorXd derivFilter = VectorXd::Zero(N);
    //Derivative filter dx / dy = (f(x+dt) - f(x-dt))/2*dt
    float dy = 2.*3.141592654/100.; //we sample at a rate of 2*pi/100 per unit
    derivFilter(mid-1) = -1. / (2.*dy);
    derivFilter(mid+1) = 1. / (2.*dy);
    
    //Now we need to reverse the filter to do a convolution instead of a correlation.
    for(int i = 0; i < derivFilter.size()/2; ++i) swap(derivFilter(i),derivFilter(derivFilter.size()-i-1));
    
    
    //Convolution (not correlation) is associative, so we can make one filter out of it before applying it on input:
    filter = discreteconvolution(meanFilter, derivFilter);
    //filter = discreteconvolution(gaussFilter, derivFilter);
    
    
    mglGraph gr, grA;
    mglData dat, datA;
    dat.Link(input.data(), input.size());
    gr.Plot(dat);
    gr.WriteFrame("/Users/kdkdk/Desktop/ETH/numerics/numcse/numcse/exercise06/inputSignal.png");
    VectorXd filteredsignal = discreteconvolution(input, filter);
    datA.Link(filteredsignal.data(), filteredsignal.size());
    gr.Plot(datA);
    grA.Plot(datA);
    grA.WriteFrame("/Users/kdkdk/Desktop/ETH/numerics/numcse/numcse/exercise06/filteredderivatedSignal.png");
    gr.WriteFrame("/Users/kdkdk/Desktop/ETH/numerics/numcse/numcse/exercise06/inputSignal.png+filteredderivatedSignal.png");

}

*/
