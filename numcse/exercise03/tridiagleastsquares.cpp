#include <iostream>

#include <Eigen/Dense>

using namespace Eigen;
using namespace std;
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
    
    MatrixXd M(z.size(),2);
    for(int i = 0; i < z.size(); ++i) {
        M(i,0) = z(i);
        M(i,1) = (i>0?z(i-1):0) + (i+1<z.size()?z(i+1):0);
    }
    
    JacobiSVD<MatrixXd> svd(M,ComputeFullU | ComputeFullV);
    VectorXd x1 = svd.solve(c);
    
    //Alternative solution
    MatrixXd LHS = M.transpose() * M;
    assert(LHS.fullPivLu().rank() == 2);
    VectorXd RHS = M.transpose() * c;
    cout << "RHS Size: " << RHS.size() << endl;
    cout << "LHS size: " << LHS.rows() << " " << LHS.cols() << endl;
    VectorXd x2 = LHS.fullPivLu().solve(RHS);
    
    //Third solution:
    auto QR = M.householderQr();
    VectorXd x3 = QR.solve(c);
    
    cout << "X1: " << endl << x1 << endl << "X2: " << endl << x2 << endl << "X3: " << endl << x3 << endl;
    assert(fabs(x1.norm() - x2.norm()) < numeric_limits<float>::epsilon()); //float not double!
    assert(fabs(x1.norm() - x3.norm()) < numeric_limits<float>::epsilon());
    assert(fabs(x2.norm() - x3.norm()) < numeric_limits<float>::epsilon());
    
    x=x3;
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
