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
	
	if (x(1) != 0.0) {
		
		double t, s, c;
		// to avoid cancellation
		if (std::abs(x(1)) > std::abs(x(0))) {
			t = x(0)/x(1);
			s = 1./std::sqrt(1 + t*t);
			c = s*t;
		} else {
			t = x(1)/x(0);
			c = 1./std::sqrt(1 + t*t);
			s = c*t;
		}

		G << c,s,-s,c; // $2 \times 2$ Givens rotation matrix
		
	} else G.setIdentity();
	y << x.norm(), 0; // $y = Gx$

}


/* @brief QR decomposition using Givens rotations
 * @param[in] A An $m \times n$ matrix
 * @param[out] R The upper triangular matrix from the QR decomposition of $A$
 * @param[out] Q The orthogonal matrix from the QR decomposition of $A$
 */
void givensQR(const MatrixXd & A, MatrixXd & Q, MatrixXd & R) {
		
	unsigned int m = A.rows();
	unsigned int n = A.cols();
	Q.setIdentity();
	R = A;
		
	Vector2d x, y;
	Matrix2d G;
		
	for (int j=0; j<n; j++)
		for (int i=m-1; i>j; i--) {
			x(0) = R(i-1,j);
			x(1) = R(i,j);
			rotInPlane(x, G, y);
			R.block(i-1,j,2,n-j) = G*R.block(i-1,j,2,n-j);
			Q.block(0,i-1,m,2) = Q.block(0,i-1,m,2)*G.transpose();
		}
	
}

int main() {
	size_t m = 3;
	size_t n = 2;
    
	MatrixXd A(m,n);
	double epsilon = std::numeric_limits<double>::denorm_min();
	A << 1, 2, 1, 1, sqrt(2), 1;
	std::cout << "A =" << std::endl << A << "\n" << std::endl;
    
	MatrixXd R(m,n), Q(m,m);

	givensQR(A, Q, R);
	std::cout << "Givens QR: R =" << std::endl << R << "\n" << std::endl;
	std::cout << "Givens QR: Q =" << std::endl << Q << std::endl;
    
	HouseholderQR<MatrixXd> QR = A.householderQr();
	Q = QR.householderQ();
	R = QR.matrixQR().triangularView<Upper>();
	std::cout << "Householder QR: R =" << std::endl << R << "\n" << std::endl;
	std::cout << "Householder QR: Q =" << std::endl << Q << std::endl;

}
