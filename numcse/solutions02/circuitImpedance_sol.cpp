#include <iostream>
#include <eigen3/Eigen/Dense>

using namespace Eigen;

void changeRx(MatrixXd &ARx, double newRx) {
// modifies ARx so that Rx = newRx
	ARx(13,14) = -1/newRx;
	ARx(14,13) = -1/newRx;
	ARx(13,13) = 3 + 1/newRx;
	ARx(14,14) = 4 + 1/newRx;
}


// Point (c)
double calc_impedance(const MatrixXd &A1inv, const VectorXd &b, double Rx, double V) {

	// create u and v
	double f = 1/Rx;
	VectorXd u = VectorXd::Zero(b.size());
	VectorXd v = VectorXd::Zero(b.size());
	u(13) = 1;
	u(14) = -1;
	v(13) = f-1;
	v(14) = 1-f;

	// Use SMW formula
	double alpha = 1 + v.transpose() * A1inv * u;
	VectorXd y = A1inv * b;
	VectorXd x = y - A1inv * u * v.transpose() * y / alpha;

	return V / (V - x(5));
}

int main() {
	MatrixXd ARx(15,15), A1inv(15,15);
	VectorXd b = VectorXd::Zero(15);
	VectorXd x = VectorXd::Zero(15);

	int V = 1;
	b(5) = V; // source term from node 16
	
	// initializes ARx with Rx = 1
	ARx <<  2,-1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	   	   -1, 4,-1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0,
 			0,-1, 3,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,
 			0, 0,-1, 3, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0,-1,
       	   -1,-1, 0, 0, 4, 0,-1, 0, 0, 0, 0, 0, 0,-1, 0,
 			0, 0, 0,-1, 0, 4, 0, 0,-1, 0, 0, 0, 0, 0,-1,
 			0, 0, 0, 0,-1, 0, 4, 0, 0,-1,-1, 0, 0, 0, 0,
 			0, 0, 0, 0, 0, 0, 0, 4,-1, 0, 0,-1,-1, 0,-1,
 			0, 0, 0, 0, 0,-1, 0,-1, 3, 0, 0, 0,-1, 0, 0,
 			0, 0, 0, 0, 0, 0,-1, 0, 0, 2,-1, 0, 0, 0, 0,
 			0, 0, 0, 0, 0, 0,-1, 0, 0,-1, 4,-1, 0, 0, 0,
 			0, 0, 0, 0, 0, 0, 0,-1, 0, 0,-1, 3,-1, 0, 0,
 			0, 0, 0, 0, 0, 0, 0,-1,-1, 0, 0,-1, 3, 0, 0,
 			0,-1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 4,-1,
 			0, 0,-1,-1, 0,-1, 0,-1, 0, 0, 0, 0, 0,-1, 5;

	A1inv = ARx.inverse();
	
	// Point (d)
	for (int i = 0; i <= 10; i++) {
		std::cout << "For Rx = 2^" << i << " the impedance is "
				  << calc_impedance(A1inv, b, 1 << i, V) << std::endl;
	}
}
