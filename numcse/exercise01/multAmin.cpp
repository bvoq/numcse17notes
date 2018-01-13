//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>

#include <Eigen/Dense>
#include <mgl2/mgl.h>
#include "timer.h"

using namespace std;
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
    // TODO: Implement an efifcient version of multAmin
    y = VectorXd(x.size());
    assert(x.size() == y.size());
    vector<double> sum = {x(0)};
    for(int i = 1; i < x.size(); ++i) sum[i] = sum[i-1]+x(i);
    double count = 0;
    for(int i = 0; i < x.size(); ++i) {
        count += sum[x.size()-1] - (i-1>=0?sum[i-1]:0);
        y(i) = count;
    }
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
    
    MatrixXd m(10,10);
    for(int i = 0; i < m.rows(); ++i) {
        for(int j = 0; j < m.cols(); ++j) {
            m(i,j) = min(i,j)+1;
        }
    }

    MatrixXd q = MatrixXd::Constant(10,10,0);
    for(int i = 0; i < q.rows(); ++i) {
        q(i,i)=2;
        if(i > 0) q(i-1,i)=-1;
        if(i > 0) q(i,i-1)=-1;
    }
    q(9,9)=1;
    MatrixXd lul = m*q;
}
