#include <iostream>
#include <vector>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

class MatrixCOO {
	// TODO: implement COO data structures, conversion to CRS and to/from dense Eigen matrix.
	// The level of abstraction is left to your choice:
	// you can delete this and use raw vectors;
	// (recommended) you can implement this as a wrapper and consider
	// only matrices with entries of the double type;
	// you can implement this as a template class.
	// The same for conversions:
	// you can work directly with the raw vectors;
	// (recommended) you can implement functions toCRS and toDense for your wrappers;
	// you can overload operator= for implicit/explicit conversion.
public:
    size_t rows, cols;
    vector<size_t> id, jd;
    vector<double> val;
    MatrixCOO(MatrixXd & m) { //from MatrixXd
        rows = m.rows();
        cols = m.cols();
        for(size_t i = 0; i < m.rows(); ++i) {
            for(size_t j = 0; j < m.cols(); ++j) {
                if(m(i,j) != 0) {
                    val.push_back(m(i,j));
                    id.push_back(i);
                    jd.push_back(j);
                }
            }
        }
    }
    MatrixXd toDenseMatrix() const {
        MatrixXd out = MatrixXd::Zero(rows, cols);
        for(int i = 0; i < val.size(); ++i) out(id[i],jd[i]) = val[i];
        return out;
    }
    friend ostream& operator<<(ostream&,const MatrixCOO &);
};
ostream& operator<<(ostream& os, const MatrixCOO & m) {
    os << m.toDenseMatrix() << endl;
    return os;
}

class MatrixCRS {
public:
	// TODO: implement CRS data structures, conversion to COO and to/from dense Eigen matrix.
    size_t rows, cols;
    vector<double> val, col_ind, row_compressed;
    MatrixCRS(MatrixXd & A) {
        size_t pcount = 0;
        rows = A.rows(); cols = A.cols();
        //optional:
        row_compressed.push_back(pcount);
        for(size_t i = 0; i < A.rows(); ++i) {
            for(size_t j = 0; j < A.cols(); ++j) {
                if(A(i,j) != 0) {
                    pcount++;
                    val.push_back(A(i,j));
                    col_ind.push_back(j);
                }
            }
            row_compressed.push_back(pcount);
        }
    }
    
    
    MatrixXd fromRowCompressed() {
        MatrixXd res  = MatrixXd::Zero(rows, cols);
        for(size_t i = 0; i < row_compressed.size()-1; ++i) {
            for(size_t j = row_compressed[i]; j < row_compressed[i+1]; ++j) {
                res(i,col_ind[j]) = val[j];
            }
        }
        return res;
    }
};





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
	
    MatrixCOO cooo (A);
    
    MatrixCRS rowform (A);
    cout << rowform.fromRowCompressed() << endl;
    Triplet<
	// TODO: test correctness of implementend conversion methods on A
}
