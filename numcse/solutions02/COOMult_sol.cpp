#include <algorithm>
#include <random>
#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <chrono>

using namespace Eigen;
using std::cout;
using std::endl;
using std::vector;

class MatrixCOO {
public:
	vector<Triplet<double>> list;

	MatrixCOO() {
	}

	MatrixCOO(vector<Triplet<double>> &list_) {
		// constructor from raw list (performs deep copy)
		list = list_;
	}

	MatrixCOO(MatrixXd &A) {
		// constructor from Eigen dense matrix
		for (int i=0; i < A.rows(); i++) {
			for (int j=0; j < A.cols(); j++) {
				if (A(i,j) != 0) {
					Triplet<double> t(i, j, A(i,j));
					list.push_back(t);
				}
			}
		}
	}

	void sort_byRow() {
		std::sort(list.begin(), list.end(),
			[] (auto t1, auto t2) {
				return t1.row() < t2.row();
			});
	}

	void sort_byCol() {
		std::sort(list.begin(), list.end(),
			[] (auto t1, auto t2) {
				return t1.col() < t2.col();
			});
	}

	int cols() {
		// returns maximum column which appears among the triplets
		int maxCol = -1;
		for (auto t : list) {
			if (maxCol < t.col()) {
				maxCol = t.col();
			}
		}
		return maxCol+1;
	}

	int rows() {
		// returns maximum row which appears among the triplets
		int maxRow = -1;
		for (auto t : list) {
			if (maxRow < t.row()) {
				maxRow = t.row();
			}
		}
		return maxRow+1;
	}

	MatrixXd toDense() {
		// returns an Eigen dense version
		MatrixXd dense = MatrixXd::Zero(rows(), cols());
		for(auto t : list) {
			dense(t.row(), t.col()) += t.value();
		}
		return dense;
	}

	
	static MatrixCOO mult_naive(MatrixCOO &A, MatrixCOO &B) {
		// naive implementation of COO matrices multiplication
		MatrixCOO result;
		for (auto t1 : A.list) {
			for (auto t2: B.list) {
				if (t1.col() == t2.row()) {
					Triplet<double> t(t1.row(), t2.col(), t1.value() * t2.value());  
					result.list.push_back(t);
				}
			}
		}
		return result;
	}

	static MatrixCOO mult_efficient(MatrixCOO &A1, MatrixCOO &A2) {
		// efficient implementation of COO matrices multiplication
		
		if (&A1 == &A2) {// to avoid complications when A1 and A2 are the
						 // same object, (deep) copy by value A2.
			MatrixCOO copyA2(A2.list);
			return mult_efficient(A1, copyA2);
		}
			
		MatrixCOO result;

		// convenience renamings
		vector<Triplet<double>> &l1 = A1.list;
		vector<Triplet<double>> &l2 = A2.list;
		vector<Triplet<double>> &lr = result.list;

		// useful way of sorting the triplets
		A1.sort_byCol();
		A2.sort_byRow();

		// build vectors of indices b1 and b2, which contain the indices
		// of elements which begin a new column for A1 and which begin a
		// new row for A2.
		vector<int> b1;
		vector<int> b2;
		b1.push_back(0);
		b2.push_back(0);

		for (int i=0; i < l1.size()-1; i++) {
			if (l1[i].col() != l1[i+1].col()) {
				b1.push_back(i+1);
			}
		}
		b1.push_back(l1.size());
				
		for (int j=0; j < l2.size()-1; j++) {
			if (l2[j].row() != l2[j+1].row()) {
				b2.push_back(j+1);
			}
		}
		b2.push_back(l2.size());

		// exploiting the special sorting of l1, l2 we are able to pass
		// in a single loop all and only the elements of the matrices 
		// A1, A2 which will be paired when one computes A1*A2.
		int i = 0;
		int j = 0;
		
		while (i < b1.size()-1 && j < b2.size()-1) {
			if (l1[b1[i]].col() == l2[b2[j]].row()) {
				// in this case, all the associated couples should be
				// multiplied and the resulting triplets
				// have to be added to result.
				int c1, c2;
				for (c1 = b1[i]; c1 < b1[i+1]; c1++) {
					for (c2 = b2[j]; c2 < b2[j+1]; c2++) {

						Triplet<double> t(l1[c1].row(), l2[c2].col(),
							l1[c1].value() * l2[c2].value());  

						lr.push_back(t);
					}
				}
				i++;
				j++;
			} else {
				if (l1[b1[i]].col() < l2[b2[j]].row()) {
					// due to sorting, if the column of l1 is smaller than
					// the row of l2, eventual couples which can contribute
					// to the product can be found only by increasing i.
					i++;
				}
				if (l1[b1[i]].col() > l2[b2[j]].row()) {
					// same as previous, only the roles are reversed.
					j++;
				}
			}
		}
		return result;
	}
};


int main() {
	srand(time(0));

	int N = 100; // size of matrix
	double sparseCoeff = 1./std::log(N); // how much sparse the matrix
										 // should be (btw 0 and 1)

	// generate full N*N matrix with random elements
	MatrixXd A = MatrixXd::Random(N,N);
	MatrixXd B = MatrixXd::Random(N,N);

	// convert to COO format
	MatrixCOO Acoo(A);
	MatrixCOO Bcoo(B);

	// define some random number generator
    std::random_device rd;
    std::mt19937 g(rd());

	// shuffle randomly the elements of the COO lists
	std::shuffle(Acoo.list.begin(), Acoo.list.end(), g);
	std::shuffle(Bcoo.list.begin(), Bcoo.list.end(), g);

	// keep only part of COO lists, proportionally to sparseCoeff
	int n = sparseCoeff*N*N;
	Acoo.list.resize(n);
	Bcoo.list.resize(n);
	cout << "Multiplication of " << N << "x" << N << " matrices with "  
		<< n << " non-zero elements. ";

	// time naive multiplication
	auto start_naive = std::chrono::system_clock::now();
	MatrixCOO::mult_naive(Acoo,Bcoo);
	auto end_naive = std::chrono::system_clock::now();


	// time efficient multiplication
	auto start_eff = std::chrono::system_clock::now();
	MatrixCOO::mult_efficient(Acoo,Bcoo);
	auto end_eff = std::chrono::system_clock::now();

	cout << "Time ratio naive/efficient implementation: " 
		<< (double) (end_naive-start_naive).count() / (end_eff-start_eff).count();
}

