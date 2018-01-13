#include <iostream>
#include <limits>

#include <Eigen/Dense>
#include <Eigen/QR>

using namespace Eigen;

/* @brief Applies Givens rotation to a 2d vector
 * @param[in] x $\mathbb{R}^2$ vector
 * @param[out] G the Givens rotation matrix
 * @param[out] y = G*x is the new $\mathbb{R}^2$ vector
 */
void rotInPlane(const Vector2d& x, Matrix2d& G, Vector2d& y) {
		//TODO
}


/* @brief QR decomposition using Givens rotations
 * @param[in] A An $m \times n$ matrix
 * @param[out] R The upper triangular matrix from the QR decomposition of $A$
 * @param[out] Q The orthogonal matrix from the QR decomposition of $A$
 */
void givensQR(const MatrixXd & A, MatrixXd & Q, MatrixXd & R) {
		//TODO
}

int main() {
	size_t m = 3;
	size_t n = 2;
    
	MatrixXd A(m,n);
	double epsilon = std::numeric_limits<double>::denorm_min();
	A << 1, 2, 1, 1, sqrt(2), 1;
	std::cout << "A =" << std::endl << A << "\n" << std::endl;
    
	MatrixXd R, Q;

	givensQR(A, R, Q);
    
	HouseholderQR<MatrixXd> QR = A.householderQr();
	Q = QR.householderQ();
	R = QR.matrixQR().triangularView<Upper>();
	std::cout << "Direct QR: R =" << std::endl << R << "\n" << std::endl;
	std::cout << "Direct QR: Q =" << std::endl << Q << std::endl;

}
