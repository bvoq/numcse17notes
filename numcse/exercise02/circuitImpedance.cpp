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
double calc_impedance(const MatrixXd &A1inv, const VectorXd &b, double Rx) {
	// TODO: calculate impedance
	// Use Sherman-Morrison-Woodbury
	return 0;
}


int main() {
	MatrixXd ARx(15,15), A1inv(15,15);
	VectorXd b = VectorXd::Zero(15);
	VectorXd x = VectorXd::Zero(15);

	b(5) = 1; // source term from node 16
	
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

	// TODO: calculate inverse of A1 
	
	// Point (d)
	// TODO: run some test
}
