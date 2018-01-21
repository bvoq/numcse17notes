// to compile run: g++ -std=gnu++11 CubicSplines_sol.cpp -lmgl
#include <cassert>
#include <iostream>
#include <Eigen/Dense>
#include <mgl2/mgl.h>

using namespace Eigen;
using namespace std;

MatrixXd cubicSpline(const VectorXd &t, const VectorXd &y) {
	// returns the matrix representing the spline interpolating the data
	// with abscissae T and ordinatae Y. Each column represents the coefficients
	// of the cubic polynomial on a subinterval.
	// Assumes T is sorted, has no repeated elements and T.size() == Y.size().
    assert(t.size() == y.size() && t.size() > 0);
	int n = t.size() - 1; // T and Y have length n+1

	MatrixXd spline(4, n);
    
    //first compute h(i) = (t - t(j-1))
    VectorXd h(n+1);
    h(0) = 1./0.; //first one = NaN, since unused due to zero indexing.
    for(int i = 1; i <= n; ++i) h(i) = t(i) - t(i-1);
    
    MatrixXd M = MatrixXd::Zero(n-1,n-1);
    //first diagonal
    for(int i = 1; i < n-1; ++i) M(i-1,i) = h(i+1)/6.;
    //second diagonal
    for(int i = 0; i < n-1; ++i) M(i,i) = (h(i+1)+h(i+2))/3.;
    //third diagonal
    for(int i = 1; i < n-1; ++i) M(i,i-1) = h(i+1)/6.;
    
    cout << M << endl;
    VectorXd r (n-1);
    for(int i = 1; i <= n-1; ++i) r(i-1) = (y(i+1)-y(i))/h(i+1) - (y(i)-y(i-1))/h(i);
    
    /*
     0.735392
     -4.59058
     -2.16766
     12.3704
     2.0497
     -21.8324
     30.2767
     */
    cout << r << endl;
    
    VectorXd sigma (n+1);
    //natural/simple interpolation
    sigma(0) = 0;
    sigma(n) = 0;
    sigma.segment(1,n-1) = M.partialPivLu().solve(r);
    
    cout << sigma << endl;
    /*
     0
     6.88222
     -21.5234
     -30.0794
     87.4006
     57.4671
     -244.768
     305.569
     0
     */
    
    for(int i = 1; i <= n; ++i) {
        spline(0,i-1) = y(i-1); //a_i
        spline(1,i-1) = (y(i)-y(i-1))/h(i) - h(i) * (2*sigma(i-1)+sigma(i))/6.; //b_i
        spline(2,i-1) = sigma(i-1)/2; //c_i
        spline(3,i-1) = (sigma(i)-sigma(i-1)) / (6. * h(i)); //d_j
    }
    
    cout << spline.col(0) << endl;
    /*
     0
     0.153066
     0
     2.38866
     */
    return spline;
}

VectorXd evalCubicSpline(const MatrixXd &spline, const VectorXd &t, const VectorXd &evalt) {
	// Returns the values of the spline calculated in the points x.
	// Assumes t is sorted, with no repetetions.

	int n = evalt.size();
	VectorXd out(n);

    //if t is sorted simple linear search suffices.
    for(int i = 0, j=0; i < n; ++i) {
        while( !(t(j) <= evalt(i) && t(j+1) >= evalt(i)) ) j++; //should never overflow, if evalt in bound of t.
        double tj = t(j);
        double a=spline(0,j), b=spline(1,j),c=spline(2,j),d=spline(3,j);
        out(i) = a + b*tj + c*tj*tj + d*tj*tj*tj;
    }

	return out;
}

int main() {
	// tests
	VectorXd T(9);
	VectorXd Y(9);
	T << 0, 0.4802, 0.7634, 1, 1.232, 1.407, 1.585, 1.879, 2;
	Y << 0., 0.338, 0.7456, 0, -1.234, 0 , 1.62, -2.123, 0;

	int len = 1 << 9;
	VectorXd evalT = VectorXd::LinSpaced(len, T(0), T(T.size()-1));

	VectorXd evalSpline = evalCubicSpline(cubicSpline(T, Y), T, evalT);
    
 	mglData datx, daty;
	datx.Link(evalT.data(), len);
	daty.Link(evalSpline.data(), len);
	mglGraph gr;
	gr.SetRanges(0, 2, -3, 3);
	gr.Plot(datx, daty, "0");
	gr.WriteFrame("/Users/kdkdk/Desktop/ETH/numerics/numcse/numcse/exercise08/spline.eps");
}
	
