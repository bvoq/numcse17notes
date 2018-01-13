#include <Eigen/Eigen>
#include <unsupported/Eigen/FFT>

#include <iostream>
#include <vector>

using namespace std;
using namespace Eigen;

vector<double> multiplyPolynomials(vector<double> a, vector<double> b) {
    vector<double> c (a.size() + b.size() - 1);
    
    int paddedZeros = 11; //can be arbitrary since this is a convolution
    VectorXcd avec = VectorXcd::Zero(a.size() + b.size() - 1 + paddedZeros), bvec = VectorXcd::Zero(a.size() + b.size() - 1 + paddedZeros);
    
    for(int i = 0; i < a.size(); ++i) avec(i) = a[i];
    for(int i = 0; i < b.size(); ++i) bvec(i) = b[i];
    
    FFT<double> fft;
    VectorXcd outvec = fft.inv(fft.fwd(avec).cwiseProduct(fft.fwd(bvec)).eval());
    
    //this outvec should be real, since a convolution results in real numbers, so only look at real part of the result.
    for(int i = 0; i < c.size(); ++i) c[i] = outvec(i).real();
    
    return c;
}

int main() {
    vector<double> a = {3,0,4,5,0,7,8};
    vector<double> b = {2,5,1,2,3};
    vector<double> c = multiplyPolynomials(a,b);
    for(double k : c) cout << k << " "; cout << endl;
    return 0;
}
