//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;
/* \brief Performs Gram-Schidt orthonormalization
 * Given a matrix $\mathbf{A}$ of linearly independent columns,
 * returns the result of a Gram-Schmidt orthonormalization.
 * Unstable GS algorithm: output is prone to cancellation issues.
 * \param[in] $\mathbf{A}$ Matrix of linearly independent columns
 * \return Matrix with ONB of $span(a_1, \cdots, a_n)$ as columns
 */
MatrixXd gram_schmidt(const MatrixXd & A) {
    // We create a matrix Q with the same size and data of A
    MatrixXd Q = MatrixXd::Constant(A.rows(), A.cols(), 0);
    for(int i = 0; i < A.cols(); ++i) {
        VectorXd vec = A.col(i);
        //3 normalization methods:
        for(int j = 0; j < i; ++j) vec -= Q.col(j).transpose() * A.col(i) * Q.col(j);
        if(vec.norm() > 0) vec.normalize();
        else {
            cerr << "Cannot normalize vector" << endl;
            return Q;
        }
        //vec = vec / vec.norm();
        //double accum = 0; for(int j = 0; j < vec.size(); ++j) accum += vec(j)*vec(j); vec = vec / sqrt(accum);
        Q.col(i) = vec;
    }

    return Q;
}

int main(void) {
    // Orthonormality test
    unsigned int n = 9;
    MatrixXd A = MatrixXd::Random(n,n);
    cout << A << endl << "-----" << endl;
    auto Q = gram_schmidt(A);
    cout << Q << endl;
    MatrixXd Q2 =  A.householderQr().householderQ();
    auto res = (Q - Q2);
    cout << "OTHER (should be equal with the exception of sign): " << endl;
    cout << Q2 << endl;
    // TODO: use gramschmidt to compute orthonormalization of
    // the matrix $\mathbf{A}$.
}
