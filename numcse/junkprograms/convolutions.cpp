#include <cassert>
#include <iostream>
#include <vector>

#include <Eigen/Eigen>
#include <unsupported/Eigen/FFT>
using namespace Eigen;
using namespace std;

vector<int> discreteconvolution (vector<int> x, vector<int> h) {
    vector<int> res(x.size()+h.size()-1, 0);
    //vector<int> res2(x.size()+h.size()-1, 0);

#define h(i) (i>=0&&i<h.size()?h[i]:0)
#define x(i) (i>=0&&i<x.size()?x[i]:0)
    for(int i = 0; i < res.size(); ++i) {
        for(int j = i; j >= 0; --j) {
            res[i] += x(i-j) * h(j);
        }
    }
    return res;
}


vector<int> periodicconvolution (vector<int> x, vector<int> h, int periodicity) {
    //1. convert h into p
    vector<int> nx (periodicity,0);
    for(int i = 0; i < x.size(); ++i) nx[i%periodicity] += x[i];
    
    vector<int> p (periodicity,0);
    for(int i = 0; i < h.size(); ++i) p[i%periodicity] += h[i];

    
    vector<int> res(periodicity,0);
    vector<int> res2(periodicity,0);

#define MOD(x,m) ((x%m)+m)%m

    for(int i = 0 ; i < periodicity; ++i) {
        for(int j = 0; j < periodicity; ++j) {
            res[i] += p[MOD(i-j,periodicity)] * nx[j]; //modded
            res2[i] += p[j] * nx[MOD(i-j,periodicity)];
        }
    }
    
    for(int i = 0; i < res.size(); ++i) assert(res[i] == res2[i]);
    return res2;
}



VectorXd periodicconvolutionFFT (vector<int> x, vector<int> h, int periodicity) {
    VectorXd nxvec = VectorXd::Zero(periodicity);
    for(int i = 0; i < x.size(); ++i) nxvec[i%periodicity] += x[i];
    
    VectorXd pvec = VectorXd::Zero(periodicity);
    for(int i = 0; i < h.size(); ++i) pvec[i%periodicity] += h[i];
    
    FFT<double > fft;
    VectorXcd tmp = (fft.fwd((VectorXcd)pvec)).cwiseProduct(fft.fwd((VectorXcd)
                                                                    nxvec));
    // doesn't work: VectorXcd y = fft.inv(tmp);
    // Alternative:
    VectorXcd y = fft.inv(tmp);
    VectorXcd y2 = 1./periodicity * fft.fwd(tmp.conjugate().eval()).conjugate();
    assert(y == y2);
    
    VectorXd realres = y2.real();
    
    return realres;
}


VectorXd discreteconvolutionFFT (vector<int> x, vector<int> h) {
    return periodicconvolutionFFT(x, h, x.size()+h.size()-1);
}

/*
float money=0;
float invest=1, previoust=1;
float exchange(float t) {
    float oldinvest = invest;
    invest = (t-previoust)/previoust;
    //sollte nicht invest = t/previoust sein?
    money += invest - 1;
    
    previoust = t;
}
*/
int main() {
    vector<int> x = {2,1,3,2}, h = {1,1,2};
    vector<int> res = discreteconvolution(x, h);

    cout << "discrete convolution: " << endl;
    for(int i = 0; i < res.size(); ++i) {
        cout << res[i] << ", ";
    }
    cout << endl;
    
    x = {2,1,3,2}, h = {1,1,2};
    cout << "periodic convolution of length 4" << endl;
    res = periodicconvolution(x, h, 4);
    for(int i = 0; i < res.size(); ++i) {
        cout << res[i] << ", ";
    }
    cout << endl;
    
    
    x = {2,1,3,2}, h = {1,1,2};
    cout << "6-periodic convolution of length 4" << endl;
    res = periodicconvolution(x, h, 6);
    for(int i = 0; i < res.size(); ++i) {
        cout << res[i] << ", ";
    }
    cout << endl;
    
    
    x = {2,1}, h = {1,1,2};
    cout << "2-periodic convolution of length 4" << endl;
    res = periodicconvolution(x, h, 2);
    for(int i = 0; i < res.size(); ++i) {
        cout << res[i] << ", ";
    }
    cout << endl;
    
    
    
    x = {2,1,3,2}, h = {1,1,2};
    cout << "periodic convolution using fft of length 4" << endl;
    VectorXd resvec = periodicconvolutionFFT(x, h, 4);
    for(int i = 0; i < resvec.size(); ++i) {
        cout << resvec[i] << ", ";
    }
    cout << endl;
    
    
    x = {2,1,3,2}, h = {1,1,2};
    cout << "periodic convolution using fft of length 6" << endl;
    resvec = periodicconvolutionFFT(x, h, 6);
    for(int i = 0; i < resvec.size(); ++i) {
        cout << resvec[i] << ", ";
    }
    cout << endl;
    
    
    x = {2,1}, h = {1,1,2};
    cout << "periodic convolution using fft of length 2" << endl;
    resvec = periodicconvolutionFFT(x, h, 2);
    for(int i = 0; i < resvec.size(); ++i) {
        cout << resvec[i] << ", ";
    }
    cout << endl;
    
    
    x = {2,1,3,2}, h = {1,1,2};
    cout << "discrete convolution using fft" << endl;
    resvec = discreteconvolutionFFT(x,h);
    for(int i = 0; i < resvec.size(); ++i) {
        cout << resvec[i] << ", ";
    }
    cout << endl;
    
    return 0;
}
