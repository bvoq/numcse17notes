#include <iostream>
#include <iomanip>
#include <vector>

#include <Eigen/Dense>
#include <mgl2/mgl.h>
#include "timer.h"

using namespace Eigen;

/* \brief compute $\mathbf{A}\mathbf{x}$
 * \mathbf{A} is defined by $(\mathbf{A})_{i,j} := \min {i,j}$
 * \param[in] x vector x for computation of A*x = y
 * \param[out] y = A*x
 */
void multAminSlow(const VectorXd & x, VectorXd & y) {
    unsigned int n = x.size();

    VectorXd one = VectorXd::Ones(n);
    VectorXd linsp = VectorXd::LinSpaced(n,1,n);
    y = ( ( one * linsp.transpose() )
          .cwiseMin( linsp * one.transpose()) ) * x;
}

/* \brief compute $\mathbf{A}\mathbf{x}$
 * \mathbf{A} is defined by $(\mathbf{A})_{i,j} := \min {i,j}$
 * Instead of a "Matlab style" construcion of the product,
 * we use simple loops.
 * \param[in] x vector x for computation of A*x = y
 * \param[out] y = A*x
 */
void multAminLoops(const VectorXd & x, VectorXd & y) {
    unsigned int n = x.size();

    MatrixXd A(n,n);

    for(unsigned int i = 0; i < n; ++i) {
        for(unsigned int j = 0; j < n; ++j) {
            A(i,j) = std::min(i+1,j+1);
        }
    }
    y = A * x;
}

/* \brief compute $\mathbf{A}\mathbf{x}$
 * This function has optimal complexity.
 * \mathbf{A} is defined by $(\mathbf{A})_{i,j} := \min {i,j}$
 * \param[in] x vector x for computation of A*x = y
 * \param[out] y = A*x
 */
void multAmin(const VectorXd & x, VectorXd & y) {
    unsigned int n = x.size();
    y = VectorXd::Zero(n);
    VectorXd v = VectorXd::Zero(n);
    VectorXd w = VectorXd::Zero(n);

    v(0) = x(n-1);
    w(0) = x(0);

    for(unsigned int j = 1; j < n; ++j) {
        v(j) = v(j-1) + x(n-j-1);
        w(j) = w(j-1) + (j+1)*x(j);
    }
    for(unsigned int j = 0; j < n-1; ++j) {
        y(j) = w(j) + v(n-j-2)*(j+1);
    }
    y(n-1) = w(n-1);
}

int main(void) {
    // Testing correctness of the code
    unsigned int M = 10;
    VectorXd xa = VectorXd::Random(M);
    VectorXd ys, yf;

    multAmin(xa, yf);
    multAminSlow(xa, ys);
    // Error should be small
    std::cout << "||ys-yf|| = " << (ys - yf).norm() << std::endl;


    unsigned int nLevels = 9;
		unsigned int *n = new unsigned int[nLevels];
		double *minTime = new double[nLevels];
		double *minTimeLoops = new double[nLevels];
		double *minTimeEff = new double[nLevels];
		
		n[0] = 4;
		for (unsigned int i=1; i<nLevels; i++)
			n[i] = 2*n[i-1];
    
    unsigned int nruns = 10;

    std::cout << "--> Timings:" << std::endl;
    // Header, see iomanip documentation
    std::cout << std::setw(15)
              << "N"
              << std::scientific << std::setprecision(3)
              << std::setw(15) << "multAminSlown"
              << std::setw(15) << "multAminLoops"
              << std::setw(15) << "multAmin"
              << std::endl;
    // From $2^4$ to $2^{13}$
    for(unsigned int i = 0; i<nLevels; i++) {
        Timer tm_slow, tm_slow_loops, tm_fast;
        // Compute runtime many times
        for(unsigned int r = 0; r < nruns; ++r) {
            VectorXd x = VectorXd::Random(n[i]);
            VectorXd y;

            // Runtime of slow method
            tm_slow.start();
            multAminSlow(x, y);
            tm_slow.stop();

            // Runtime of slow method with loops
            tm_slow_loops.start();
            multAminLoops(x, y);
            tm_slow_loops.stop();

            // Runtime of fast method
            tm_fast.start();
            multAmin(x, y);
            tm_fast.stop();
        }
        
        minTime[i] = tm_slow.min();
        minTimeLoops[i] = tm_slow_loops.min();
        minTimeEff[i] = tm_fast.min();

        std::cout << std::setw(15)
                  << n[i]
                  << std::scientific << std::setprecision(3)
                  << std::setw(15) << minTime[i]
                  << std::setw(15) << minTimeLoops[i]
                  << std::setw(15) << minTimeEff[i]
                  << std::endl;
    }
    
    // Plotting with MathGL
    double nMgl[nLevels];
    double ref1[nLevels], ref2[nLevels];
    for (int i=0; i<nLevels; i++) {
    	nMgl[i] = n[i];
    	ref1[i] = 1e-8*pow(n[i],2);
    	ref2[i] = 1e-7*n[i];
    }
    
    mglData matSize;
    matSize.Link(nMgl, nLevels);
    
    mglData data1, data2;
    mglData dataRef1, dataRef2;
  	data1.Link(minTime, nLevels);
  	data2.Link(minTimeEff, nLevels);
  	dataRef1.Link(ref1,nLevels);
  	dataRef2.Link(ref2,nLevels);
  	
  	mglGraph *gr = new mglGraph;
    gr->Title("Runtime of multAmin");
  	gr->SetRanges(n[0],n[0]*pow(2,nLevels-1),1e-6,1e+1);  gr->SetFunc("lg(x)","lg(y)");
  	gr->Axis();
  	gr->Plot(matSize,data1,"k +"); gr->AddLegend("slow","k +");
  	gr->Plot(matSize,data2,"r +"); gr->AddLegend("efficient","r +");
  	gr->Plot(matSize,dataRef1,"k"); gr->AddLegend("O(n^2)","k");
  	gr->Plot(matSize,dataRef2,"r"); gr->AddLegend("O(n)","r");
  	gr->Label('x',"Matrix size [n]",0);
  	gr->Label('y', "Runtime [s]",0);
    gr->Legend(2);
		gr->WriteFrame("multAmin_comparison.eps");


    // The following code is just for demonstration purposes.
    // Build Matrix B with dimension 10x10 such that B = inv(A)
    unsigned int nn = 10;
    MatrixXd B = MatrixXd::Zero(nn,nn);
    for(unsigned int i = 0; i < nn; ++i) {
        B(i,i) = 2;
        if(i < nn-1) B(i+1,i) = -1;
        if(i > 0) B(i-1,i) = -1;
    }
    B(nn-1,nn-1) = 1;

    std::cout << "B = " << std::endl
              << B << std::endl;

    // Check that B = inv(A) (up to machine precision)
    std::cout << "--> Test B = inv(A):" << std::endl;
    VectorXd x = VectorXd::Random(nn), y;
    multAmin(B*x, y);
    std::cout << "|y-x| = " << (y - x).norm() << std::endl;
    multAminSlow(B*x, y);
    std::cout << "|y-x| = " << (y - x).norm() << std::endl;
    multAminLoops(B*x, y);
    std::cout << "|y-x| = " << (y - x).norm() << std::endl;
}
