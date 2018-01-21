//
//  chebbyapprox.cpp
//  numcse
//
//  Created by kdkdk on 15.01.18.
//  Copyright © 2018 a. All rights reserved.
//

#include <stdio.h>
#include <Eigen/Eigen>
#include <unsupported/Eigen/FFT>
#include <mgl2/mgl.h>
#include <string>

using namespace Eigen;
using namespace std;


double function1(double x) {
    return 1/(x*x + 5);
}

#include <cmath>
VectorXd chebyshevNodes(int n, double a, double b) {
    VectorXd nodes(n+1);
    for(int j = 0; j <= n; ++j) {
        double rootj = cos((2.*j+1)/(2.*(n+1)) * M_PI);
        nodes(j) = a+((rootj+1)*(b-a)/2.);
    }
    return nodes;
}

VectorXd newtonCoeffDivDiff(const VectorXd & t, const VectorXd& y) {
    VectorXd a = y;
    int n = t.size()-1;
    for(unsigned l=0; l<n; ++l) for(unsigned j=l; j<n; ++j) a(j+1) = (a(j+1)-a(l))/(t(j+1)-t(l));
    return a;
}

VectorXd evalNewton(const VectorXd & nodes, const VectorXd & coeff, const VectorXd & evalnodes) {
    int n = nodes.size()-1;
    int k = 0;
    VectorXd vec(evalnodes.size());
    for(int k = 0; k < evalnodes.size(); ++k) {
        double newtoncoeff = 1.;
        double res = 0;
        for(int i = 0; i <=n; ++i) {
            res += coeff(i) * newtoncoeff;
            newtoncoeff *= (evalnodes(k)-nodes(i));
        }
        vec(k) = res;
    }
    return vec;
}

VectorXd chebexp ( const VectorXd& y ) {
    const int n = y.size()-1;
    const std::complex<double> M_I(0, 1); //imaginaryunit
    // create vector z, see (6.1.110)
    VectorXcd b(2*(n + 1));
    const std::complex<double> om =-M_I*(M_PI*n)/((double)(n+1));
    for (int j = 0; j <= n; ++j) {
        b(j) = std::exp(om*double(j))*y(j);
        b(2*n+1-j) = std::exp(om*double(2*n+1-j))*y(j);
    }
    // Solve linear system (6.1.111) with effort O(n log n)
    Eigen::FFT<double> fft ; // EIGEN' helper class for DFT
    VectorXcd c = fft.inv(b) ; // -> c = ifft(z), inverse fourier transform
    // recover j, see (6.1.111)
    VectorXd beta(c.size());
    const std::complex<double> sc = M_PI_2 / ( n + 1)*M_I ;
    for (unsigned j = 0; j < c.size(); ++j)
        beta(j) = ( std::exp(sc*double(-n+j))*c[j] ).real(); // recover a, see (6.1.108)
    VectorXd alpha = 2*beta.tail(n); alpha(0) = beta(n);
    return alpha ;
}


double evalChebbyExp(const VectorXd & chebbycoeff, double t) {
    int n = chebbycoeff.size()-1;
    vector<double> chebyeval (n+1);
    chebyeval[0] = 1;
    chebyeval[1] = t;
    for(int i = 2; i <= n; ++i) chebyeval[i] = 2. * t *  chebyeval[i-1] - chebyeval[i-2];
    double res = 0;
    for(int i = 0; i < n; ++i) {
        double aprime = chebbycoeff(i);
        if(i == n-1) aprime = chebbycoeff(i) + 2 * t * chebbycoeff(i+1);
        if(i == n-2) aprime = chebbycoeff(i) - chebbycoeff(i+2);
        res += chebyeval[i] * aprime;
    }
    return -res;
}


double recclenshaw(const VectorXd& a, const double x) {
    const VectorXd::Index n = a.size() - 1;
    if (n == 0) return a(0) ; // Constant polynomial else if (n == 1) return (x∗a(1) + a(0));//Valueα1∗x+α0
    else {
        VectorXd new_a(n);
        new_a << a.head(n - 2) , a(n-2) - a(n) , a(n-1)+ 2*x*a(n) ;
        return recclenshaw(new_a,x) ; // recursion
    }
}

int main() {
    
    mglGraph grTot;

    for(int n = 1; n < 100; n++) {
        if(n>9) n+=5;
        int a = -1, b = 1;
        VectorXd nodes = chebyshevNodes(n, a, b);
        VectorXd f_nodes (n+1);
        for(int i = 0; i <= n; ++i) {
            double x = nodes(i);
            f_nodes(i) = 1 / (x*x+1);
        }
    
        VectorXd coeff = newtonCoeffDivDiff(nodes, f_nodes);
    
        VectorXd evalAt (301);
        for(int i = 0 ; i <= 300; ++i) evalAt(i) = a+(b-a)/300.*i;
        VectorXd res = evalNewton(nodes, coeff, evalAt);
    
        //doesn't work:
        //VectorXd chebbycoeff = chebexp(f_nodes);
        //VectorXd res2(evalAt.size());
        //for(int i = 0; i <= 300; ++i) res2(i) = recclenshaw(chebbycoeff,evalAt(i));
        
        mglGraph gr, gr2;
        mglData dat, dat2;
        dat.Link(res.data(), res.size());
        //dat2.Link(res2.data(), res2.size());
        gr.Plot(dat);
        grTot.Plot(dat);
        string out = "/Users/kdkdk/Desktop/ETH/numerics/numcse/numcse/junkprograms/approx"+to_string(n)+".png";
        gr.WriteFrame(out.c_str());
    }
    
    grTot.WriteFrame("/Users/kdkdk/Desktop/ETH/numerics/numcse/numcse/junkprograms/allthings.png");
    /*
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
    */
    return 0;
}
