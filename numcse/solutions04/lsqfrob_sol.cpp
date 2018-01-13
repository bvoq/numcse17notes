#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/unsupported/Eigen/KroneckerProduct>

using namespace Eigen;

MatrixXd min_frob(const VectorXd & z, const VectorXd & g) {

    int n = g.size();
    
    // Build temporary matrix C
    MatrixXd C = kroneckerProduct(MatrixXd::Identity(n,n),
                                  z.transpose());

    // Build system matrix and r.h.s.
    MatrixXd S(n*n+n,n*n+n);
    S << MatrixXd::Identity(n*n, n*n), C.transpose(),
         C                           , MatrixXd::Zero(n,n);
    VectorXd h(n*n+n);
    h << VectorXd::Zero(n*n), g;

    // Solve augmented system and return only head of solution
    // as a matrix (discard the rest).
    return Map<const MatrixXd>(S.fullPivLu().solve(h).eval().data(),
                               n, n).transpose();
    // It is possible to solve the system more efficiently
    // exploiting the special structure of the matrix.
}

int main() {

    int n = 10;

    VectorXd z = VectorXd::Random(n),
             g = VectorXd::Random(n);

    MatrixXd M_augmNormal = min_frob(z, g); 
    MatrixXd M_lagrange = g*z.transpose() / z.squaredNorm();

    std::cout << "Norm: "
              << (M_augmNormal - M_lagrange).norm()
              << std::endl;
}
