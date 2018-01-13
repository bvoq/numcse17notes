#include <iostream>

#include <Eigen/Dense>

using namespace Eigen;

/* @brief 
 * @param[in] z An $n$-dimensional vector containing one side of input data
 * @param[in] c An $n$-dimensional vector containing the other side of input data
 * @param[out] x The vector of parameters $(\alpha,\beta)$, intercept and slope of the line fitted
 */
VectorXd lsqEst(const VectorXd &z, const VectorXd &c)
{
    // Initialization
	int n = z.size();
    assert( z.size() == c.size() && "z and c must have same size");

	VectorXd x(2);

	MatrixXd A(n,2);
	A.col(0) = z;
	A(0,1) = z(1);
	for(size_t i=1; i<n-1; ++i) {
		A(i,1) = z(i-1) + z(i+1);
	}
	A(n-1,1) = z(n-2);
	
	// Normal equations
	MatrixXd lhs = A.transpose() * A; // Left-hand side
	VectorXd rhs = A.transpose() * c; // Right-hand side
	x = lhs.fullPivLu().solve(rhs);
    
	return x;
}

int main() {
    // Initialization
    unsigned int n = 10;
    VectorXd z(n), c(n);
    for(size_t i=0; i<n; ++i) {
		z(i) = i+1;
		c(i) = n-i;
	}

	VectorXd x = lsqEst(z, c);

	std::cout << "alpha = " << x(0) << std::endl;
	std::cout << "beta = "  << x(1) << std::endl;
}
