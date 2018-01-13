#include <cmath>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "chebPolyEval.hpp"
#include <mgl2/mgl.h>

using namespace std;

// Compute the best approximation of the function $f$
// with Chebyshev polynomials.
// $\alpha$ is the output vector of coefficients.
template <typename Function>
void bestApproxCheb(const Function &f, Eigen::VectorXd &alpha) {
    int n=alpha.size()-1;
    Eigen::VectorXd fn(n+1);
    for (int k=0; k<n+1; k++) {
        double temp=cos(M_PI*(2*k+1)/2/(n+1));
        fn(k)=f(temp);
    }
    
    vector<double> V;
    Eigen::MatrixXd scal(n+1,n+1);
    for (int j=0; j<n+1; j++) {
        V=chebPolyEval(n,cos(M_PI*(2*j+1)/2/(n+1)));
        for (int k=0; k<n+1; k++) scal(j,k)=V[k];
    }
    
    for (int k=0; k<n+1; k++) {
        alpha(k)=0;
        for (int j=0; j<n+1; j++) {
            alpha(k)+=2*fn(j)*scal(j,k)/(n+1);
        }
    }
    alpha(0)=alpha(0)/2;
}


int main(){
	int n;
	
// Check the orthogonality of Chebyshev polynomials
    n=10;
    vector<double> V;
    Eigen::MatrixXd scal(n+1,n+1);
    for (int j=0; j<n+1; j++) {
        V=chebPolyEval(n,cos(M_PI*(2*j+1)/2/(n+1)));
        for (int k=0; k<n+1; k++) scal(j,k)=V[k];
    }
    
    double maxOrthErr = 1e-40;
    for (int k=0; k<n+1; k++)
        for (int l=k+1; l<n+1; l++)
            maxOrthErr = max( maxOrthErr, scal.col(k).dot(scal.col(l)) );
    cout<< "Maximum orthogonality error: " << maxOrthErr <<endl;

// Test the implementation
    auto f = [] (double & x) {return 1/(pow(5*x,2)+1);};
    n=20;
    Eigen::VectorXd alpha(n+1);
    bestApproxCheb(f, alpha);
    
    // Compute the error
    int nPts = 1e+6;
    Eigen::VectorXd X = Eigen::VectorXd::LinSpaced(nPts,-1,1);
    auto qn = [&alpha,&n] (double & x) {
        double val=0;
        vector<double> V=chebPolyEval(n,x);
        for (int k=0; k<n+1; k++) val+=alpha(k)*V[k];
        return val;
    };
    
    Eigen::VectorXd polyExact(nPts), polyApprox(nPts);
    for (int i=0; i<nPts; i++) {
    	polyExact(i) = f(X(i));
    	polyApprox(i) = qn(X(i));
    }
    
    double err_max=1e-20;
    for (int i=0; i<nPts; i++)
    	err_max=std::max(err_max,abs(polyExact(i)-polyApprox(i)));
    cout<<"Error: "<< err_max <<endl;
    
    // Plot
  	/*mglData data1, data2, coordX;
  	data1.Link(polyExact.data(), nPts);
  	data2.Link(polyApprox.data(), nPts);
  	coordX.Link(X.data(), nPts);
  	
  	mglGraph *gr = new mglGraph;
  	gr->Title("f(x) = \\1/(1 + (5x)^2)");
  	gr->SetRanges(-1,1,0,1.2);  gr->SetFunc("x","y");
  	gr->Axis();
  	gr->Plot(coordX, data1,"b"); gr->AddLegend("exact","b");
  	gr->Plot(coordX, data2,"r"); gr->AddLegend("bestApproxCheb","r");
  	gr->Label('x',"x",0);
  	gr->Label('y', "f(x)",0);
  	gr->Legend();
		gr->WriteFrame("bestApproxCheb.png");*/
    
}
