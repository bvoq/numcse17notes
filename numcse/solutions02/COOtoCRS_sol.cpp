#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

using namespace Eigen;

class MatrixCOO {
public:
	std::vector<Triplet<double>> list;

	MatrixCOO() {
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
		// std::sort can sort by applying a specific comparison function,
		// which is what is passed hereafter as third parameter using
		// the lambda expressions syntax.
		std::sort(list.begin(), list.end(),
			[] (auto t1, auto t2) {
				return t1.row() < t2.row();
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
		MatrixXd dense(this->rows(), this->cols());
		for (auto t : list) {
			dense(t.row(), t.col()) = t.value();
		}
		return dense;
	}

};

class MatrixCRS {
public:
	std::vector<double> val;
	std::vector<int> col_ind;
	std::vector<int> row_ptr;

	MatrixCRS() {
	}

	MatrixCRS(MatrixXd &A) {
		// constructor from Eigen dense matrix
		int nnz = 0; // number of non-zero elements
		int prevNotEmptyRow = -1;

		// fill up val and colind by looping through triplets
		// some additional care must be taken for empty rows for rowptr
		for (int i=0; i < A.rows(); i++) {
			bool emptyRow = true;

			for (int j=0; j < A.cols(); j++) {
				if (A(i,j) != 0) {
					nnz++;
					val.push_back(A(i,j));
					col_ind.push_back(j);
					emptyRow = false;
				}
			}

			if (!emptyRow) {
				for (int c = prevNotEmptyRow; c < i; c++) {
					row_ptr.push_back(nnz);
				}
				prevNotEmptyRow = i;
			}
		}
	}


	MatrixCRS(MatrixCOO &A) {
		// constructor from COO matrix
		A.sort_byRow(); // bring A to a more maneagable form

		int prevNotEmptyRow = -1;
		std::vector<Triplet<double>> &l = A.list;

		for (int i=0; i < l.size(); i++) {

			while (prevNotEmptyRow < l[i].row()) {
				// if there are some empty rows, fill up rowptr
				// following the conventions.
				prevNotEmptyRow++;
				this->row_ptr.push_back(i);
			}

			double currVal = l[i].value();

			while (l[i].col() == l[i+1].col() &&
				   l[i].row() == l[i+1].row()) {
				// add up the values corresponding to the same indices
				currVal += l[i].value();
				i++;
			}

			this->val.push_back(currVal);
			this->col_ind.push_back(l[i].col());
		}
	}

	int cols() {
		// returns maximum column which appears 
		int maxCol = -1;
		for (auto c : col_ind) {
			if (maxCol < c) {
				maxCol = c;
			}
		}
		return maxCol+1;
	}

	int rows() {
		// returns maximum row which appears 
		return row_ptr.size();
	}

	MatrixXd toDense() {
		// conversion to Eigen dense matrix
		int nRows = this->rows();
		MatrixXd dense(nRows, this->cols());
		int row = 0;
		int i = 0;

		while (row < nRows) {
			do { // some extra care for empty rows
				row++;
			} while (row_ptr[row-1] == row_ptr[row]);

			while (i < row_ptr[row]) {
				// insert all elements of a certain row
				dense(row-1, col_ind[i]) = val[i];
				i++;
			}
		}

		while (i < val.size()) {
			// insert last row; this step would not be necessary if
			// we had directly appended to rowptr the size of val.
			dense(row-1, col_ind[i]) = val[i];
			i++;
		}

		return dense;
	}
};

// correctness tests
int main() {
	MatrixXd A(18,17);
	A<< 2,-1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	   -1, 4,-1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0,
 		0,-1, 3,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,
 		0, 0,-1, 3, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0,-1,
 		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       -1,-1, 0, 0, 4, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0,-1, 0,
 		0, 0, 0,-1, 0, 0, 0, 4, 0, 0,-1, 0, 0, 0, 0, 0,-1,
 		0, 0, 0, 0,-1, 0, 0, 0, 4, 0, 0,-1,-1, 0, 0, 0, 0,
 		0, 0, 0, 0, 0, 0, 0, 0, 0, 4,-1, 0, 0,-1,-1, 0,-1,
 		0, 0, 0, 0, 0, 0, 0,-1, 0,-1, 3, 0, 0, 0,-1, 0, 0,
 		0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 2,-1, 0, 0, 0, 0,
 		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 		0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0,-1, 4,-1, 0, 0, 0,
 		0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0,-1, 3,-1, 0, 0,
 		0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1, 0, 0,-1, 3, 0, 0,
 		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 		0,-1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4,-1,
 		0, 0,-1,-1, 0, 0, 0,-1, 0,-1, 0, 0, 0, 0, 0,-1, 5;
	
	MatrixCOO B(A);
	MatrixCRS C(B);
	MatrixXd D = C.toDense();
	std::cout << D;
}
