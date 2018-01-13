#include <iostream>
#include <limits>
#include <vector>
#include <Eigen/Eigen>
#include <Eigen/SVD>
using namespace Eigen;
using namespace std;

#define EQU(a,b) (fabs(a-b) < numeric_limits<double>::epsilon())
bool EQUM(const MatrixXd & A, const MatrixXd & B) {
    if(A.cols() != B.cols() || A.rows() != B.rows()) return false;
    for(int i = 0; i < A.size(); ++i) {
        if(!EQU(A(i),B(i))) return false;
    }
    return true;
}
int main(void) {
    MatrixXd A (4,3);
    A << 1,2,3,   2,4,6,    3,3,3,    1,2,3;
    
    cout << A.fullPivLu().rank() << endl;
    cout << A.jacobiSvd().rank() << endl;
    cout << A.colPivHouseholderQr().rank() << endl;
    
    cout << A.fullPivLu().matrixLU() << endl;
    
    JacobiSVD<MatrixXd> SVD(A, ComputeFullU | ComputeFullV);
    cout << "U: " << endl << SVD.matrixU() << endl;
    cout << "V: " << endl << SVD.matrixV() << endl;
    cout << "lcols R(A): " << endl << SVD.matrixU().leftCols(SVD.rank()) << endl << endl;
    cout << "rcols N(A^T): " << endl << SVD.matrixU().rightCols(A.rows()-SVD.rank()) << endl;

    //assert(A.jacobiSvd().nonzeroSingularValues().size() == A.jacobiSvd().rank());
    
    /*
    MatrixXd U (3,3);
    U << 3,3,3,0,2,1,0,0,-0.5;
    MatrixXd L(4,3);
    L << 1,0,0, 0.333333333,1,0,   0.6666666666,0.5,1,  0.33333333,1,0;
    cout << "Should be A: " << endl << A.fullPivLu().permutationP().inverse() * L * U * A.fullPivLu().permutationQ() << endl;
    */
    cout << "TEST FOR EQUALITY:\n" << EQUM(A,A.transpose()) << endl;
    return 0;
}
