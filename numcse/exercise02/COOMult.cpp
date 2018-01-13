#include <iostream>
#include <Eigen/Dense>
#include <chrono>
#include <vector>

#define sep "-------------" << endl

using namespace Eigen;
using namespace std;

struct Tripplet {
    size_t row, col;
    double val;
};

class MatrixCOO { // re-use the class you implemented in exercise 2.4
				  // or delete this and work with raw vectors
public:
    size_t rows, cols;
    vector<Tripplet> tripplets;
    MatrixCOO() {}
    MatrixCOO(MatrixXd & m) { //from MatrixXd
        rows = m.rows();
        cols = m.cols();
        for(size_t i = 0; i < m.rows(); ++i) {
            for(size_t j = 0; j < m.cols(); ++j) {
                if(m(i,j) != 0) {
                    tripplets.push_back({i,j,m(i,j)});
                }
            }
        }
    }
    MatrixXd toDenseMatrix() const {
        MatrixXd out = MatrixXd::Zero(rows, cols);
        for(int i = 0; i < tripplets.size(); ++i) out(tripplets[i].row,tripplets[i].col) += tripplets[i].val;
        return out;
    }
	// Point (b)
	static MatrixCOO mult_naive (MatrixCOO &A, MatrixCOO &B) {
		assert(A.cols == B.rows);
        
        MatrixCOO result;
        result.cols = A.cols;
        result.rows = B.rows;
        int rescount = 0;
        for(int i = 0; i < A.tripplets.size(); ++i) {
            for(int j = 0; j < B.tripplets.size(); ++j) {
                if(A.tripplets[i].col == B.tripplets[j].row) {
                    result.tripplets.push_back({A.tripplets[i].row, B.tripplets[j].col, A.tripplets[i].val * B.tripplets[j].val});
                }
            }
        }
		return result;
	}
    
    //VOILÀ
	static MatrixCOO mult_efficient (MatrixCOO &A, MatrixCOO &B) {
		MatrixCOO result;
        result.cols = A.cols; result.rows = B.rows;
        vector<Tripplet> a = A.tripplets, b = B.tripplets;
        
        sort(a.begin(), a.end(), [](const Tripplet & a, const Tripplet & b) { return a.col < b.col; });
        sort(b.begin(), b.end(), [](const Tripplet & a, const Tripplet & b) { return a.row < b.row; });
        for(int i = 0; i < a.size(); ++i) {
            static int previousBPos = 0, j = 0;
            if(i-1>=0 && a[i-1].col == a[i].col) j = previousBPos; previousBPos = j;
            for(; j < b.size() && a[i].col < b[j].row; ++j);
            for(; j < b.size() && a[i].col == b[j].row; ++j) result.tripplets.push_back({a[i].row, b[j].col, a[i].val * b[j].val});
        }
        
		return result;
	}
};

int main() {
	// Point (f)
	// TODO
    MatrixXd Adense (4,4);
    Adense <<   1,2,3,4,
                5,6,7,8,
                9,8,7,6,
                5,4,3,2;
    MatrixCOO A (Adense);
    MatrixXd Bdense (4,4);
    Bdense <<   4,4,4,4,
                4,4,4,4,
                4,4,4,4,
                4,4,4,4;
    MatrixCOO B (Bdense);
    
    MatrixCOO C = MatrixCOO::mult_naive(A, B);
    MatrixCOO D = MatrixCOO::mult_efficient(A,B);

    cout << C.toDenseMatrix() << endl;
    cout << sep << Adense*Bdense << endl;
    cout << sep << D.toDenseMatrix() << endl;
}
