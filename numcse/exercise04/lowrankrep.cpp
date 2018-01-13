//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati, aldabrow
//// This file is part of the NumCSE repository.
////
#include <iostream>
#include <limits>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/QR>
#include <eigen3/Eigen/SVD>

using namespace Eigen;

void factorize(const MatrixXd & X, int k,
					 MatrixXd & A, MatrixXd & B) {
	// TODO					 	
}

void svd_AB(const MatrixXd & A, const MatrixXd & B,
				  MatrixXd & U, MatrixXd & S, MatrixXd & V) {
	
	// TODO
}

void rank_k_approx(const MatrixXd &Ax, const MatrixXd &Ay,
				   const MatrixXd &Bx, const MatrixXd &By,
				   MatrixXd &Az, MatrixXd &Bz) {
				   
	// TODO
}

int main() {
	int m = 3;
	int n = 2;
	int k = 2;
	MatrixXd X(m,n);
	X << 5, 0, 2, 1, 7, 4;
	
	MatrixXd A, B;
	factorize(X, k, A, B);
	
	std::cout << "A =" << std::endl << A << std::endl;
	std::cout << "B =" << std::endl << B << std::endl;
	
	A.resize(m,k); B.resize(n,k);
	A << 2, 1, 2, 3, 6, 1;
	B << 4, 4, 5, 0;
	MatrixXd U, S, V;
	
	svd_AB(A, B, U, S, V);
	
	std::cout << "U =" << std::endl << U << std::endl;
	std::cout << "S =" << std::endl << S << std::endl;
	std::cout << "V =" << std::endl << V << std::endl;
	
	MatrixXd Ax(m,k), Ay(m,k), Bx(n,k), By(n,k), Az, Bz;
	Ax << 1,  0, 9, 2, 6, 3;
	Ay << 8, -2, 3, 4, 5, 8;
	Bx << 2, 1, 2, 3;
	By << 4, 4, 5, 0;
	
	rank_k_approx(Ax, Ay, Bx, By, Az, Bz);
	
	std::cout << "Az =" << std::endl << Az << std::endl;
	std::cout << "Bz =" << std::endl << Bz << std::endl;
}
