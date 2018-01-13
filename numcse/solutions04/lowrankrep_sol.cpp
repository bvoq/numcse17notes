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
using std::min;

void factorize(const MatrixXd & X, int k,
					 MatrixXd & A, MatrixXd & B) {
					 	
	int m = X.rows();
	int n = X.cols();
	
	// You can ask for thin U or V to be computed
	JacobiSVD<MatrixXd> svd(X, ComputeThinU | ComputeThinV);
	
	VectorXd s = svd.singularValues().head(k);
	MatrixXd U = svd.matrixU();
	MatrixXd V = svd.matrixV();

	int ns = s.size();
	int safeMax = min(k, ns); // in case k > number sing. values
	
	A = U.leftCols(safeMax).array().rowwise() * 
		s.head(safeMax).transpose().array();
	B = V.leftCols(safeMax);
}

void svd_AB(const MatrixXd & A, const MatrixXd & B,
				  MatrixXd & U, MatrixXd & S, MatrixXd & V) {
	
	int m = A.rows();
	int n = B.rows();
	int k = A.cols();
	
	// QA: m x k; QB: n x k; RA,RB k x k
	HouseholderQR<MatrixXd> QRA = A.householderQr();
	MatrixXd QA = QRA.householderQ() * MatrixXd::Identity(m, min(m, k));
	MatrixXd RA = MatrixXd::Identity(min(m, k), m) *
		QRA.matrixQR().triangularView<Upper>();
	HouseholderQR<MatrixXd> QRB = B.householderQr();
	MatrixXd QB = QRB.householderQ() * MatrixXd::Identity(n, min(n, k));
	MatrixXd RB = MatrixXd::Identity(min(n, k), n) *
		QRB.matrixQR().triangularView<Upper>();
	
	// U,V: k x k
	JacobiSVD<MatrixXd> svd(RA*RB.transpose(), ComputeFullU | ComputeFullV);
	// Thin matrices are unnecessary here as RA*RB' is
	// a square k times k matrix
		
	S = svd.singularValues().asDiagonal();
	U = svd.matrixU();
	V = svd.matrixV();
	
	// U: m x k; V: n x k
	U = QA*U;
	V = QB*V;
}

void rank_k_approx(const MatrixXd &Ax, const MatrixXd &Ay,
				   const MatrixXd &Bx, const MatrixXd &By,
				   MatrixXd &Az, MatrixXd &Bz) {
	
	MatrixXd A(Ax.rows(), Ax.cols()+Ay.cols());
	A << Ax, Ay;
	MatrixXd B(Bx.rows(), Bx.cols()+By.cols());
	B << Bx, By;

	// U: m x 2k; S: 2k x 2k; V: n x 2k
	MatrixXd U, S, V;
	svd_AB(A, B, U, S, V);
	
	int k = Ax.cols();
	Az = U.leftCols(k) * S;
	Bz = V.leftCols(k);
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
