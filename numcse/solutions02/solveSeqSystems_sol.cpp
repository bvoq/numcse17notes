#include <iostream>
#include <eigen3/Eigen/Dense>
#include <ctime>

using namespace Eigen;
using namespace std;

// Point (a)

void solve_naive(const MatrixXd &A, const MatrixXd &b, MatrixXd &X) {

	for(int i = 0; i < b.cols(); i++) {
		X.col(i) = A.fullPivLu().solve(b.col(i));
	}
}

// Point (b)

void solve_LU(const MatrixXd &A, const MatrixXd &b, MatrixXd &X) {

	FullPivLU<MatrixXd> LU = A.fullPivLu();

	for(int i = 0; i < b.cols(); i++) {
		X.col(i) = LU.solve(b.col(i));
	}
}

// some tests with random matrices

int main() {
	int n = 3;
	int m = 5;

	srand(time(0));
	MatrixXd A = MatrixXd::Random(n,n);
	MatrixXd b = MatrixXd::Random(n,m);
	MatrixXd X1(n,m);
	MatrixXd X2(n,m);

	solve_naive(A,b,X1);
	cout << A*X1 << "\n\n" << b << "\n\n";

	solve_LU(A,b,X2);
	cout << A*X2 << "\n\n" << b << endl;
}
