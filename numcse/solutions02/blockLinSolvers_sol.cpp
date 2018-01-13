#include <iostream>
int main() {}
/*
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

// Point (c)

VectorXd solve_backSub(const MatrixXd& R, const VectorXd& c)
{
	int n = R.rows();
	VectorXd y(n);

	// Since R is upper triangular, we can solve by backwards substitution
	for (int i = n-1; i >= 0; --i)
	{
		y(i) = c(i);
		for (int j = n-1; j > i ; --j)
		{
			y(i) -= R(i,j) * y(j);
		}
		y(i) /= R(i,i);
	}
	return y;
}

void solve_blockLU(const MatrixXd &R,
				   const VectorXd &u,
				   const VectorXd &v,
				   const VectorXd &b,
				   VectorXd &x) {

	int n = R.rows();
	VectorXd y(n+1);
	//
	// Solve Ly = b by forward substitution.
	y.head(n) = b.head(n);
	y(n) = b(n) - u.transpose() * solve_backSub(R, y.head(n));
	
	// Solve Ux = y by backward substitution.
	MatrixXd U(n+1,n+1);
	U << R,							 	v,
		 VectorXd::Zero(n).transpose(), -u.transpose()*solve_backSub(R,v);

	x = solve_backSub(U,y);
}


// Point (f)

void solve_blockGauss(const MatrixXd &R,
				  	  const VectorXd &u,
					  const VectorXd &v,
					  const VectorXd &b,
					  VectorXd &x) {

	int n = R.rows();

	const TriangularView<const MatrixXd, Upper> & triR = R.triangularView<Upper>();

	double sinv = - 1. / u.dot(triR.solve(v));
	double bs = (b(n) - u.dot(triR.solve(b.head(n))));
	double sinvbs = sinv*bs;

	x << triR.solve(b.head(n) - v*sinvbs),
		 sinvbs;
}

// some tests on random matrices
int main() {
	int N = 5;
	MatrixXd R = MatrixXd::Random(N,N);
	R = R.triangularView<Upper>();
	VectorXd u = VectorXd::Random(N);
	VectorXd v = VectorXd::Random(N);
	VectorXd b = VectorXd::Random(N+1);
	VectorXd x1 = VectorXd::Random(N+1);
	VectorXd x2 = VectorXd::Random(N+1);

	solve_blockLU(R, u, v, b, x1);
	solve_blockGauss(R, u, v, b, x2);

	MatrixXd A(N+1,N+1);
	A << R,				v,
		 u.transpose(), 0;

	cout << A.fullPivLu().solve(b) << "\n\n";
	cout << x1 << "\n\n";
	cout << x2 << "\n\n";

}
*/
