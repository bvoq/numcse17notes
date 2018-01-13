#include <Eigen/QR>
#include <iostream>
using namespace Eigen;
using namespace std;

pair<MatrixXd, MatrixXd> qr_full(const MatrixXd & A) {
    
    HouseholderQR<MatrixXd> qr(A);
    auto Q = qr.householderQ();
    auto R = qr.matrixQR().template triangularView<Upper>();
    return {Q,R};
}


pair<MatrixXd, MatrixXd> qr_economic(const MatrixXd & A) {
    HouseholderQR<MatrixXd> qr(A);
    auto Q = qr.householderQ() * MatrixXd::Identity(A.rows(), A.cols());
    auto R = qr.matrixQR().block(0,0,A.cols(),A.cols()).template triangularView<Upper>();
    return {Q,R};
}
int main() {
    MatrixXd A (3,3);
    A << 1,1,2,    2,-3,0,    2,4,-4;
    
    cout << qr_full(A).second << endl;
    cout << qr_full(A).first << endl;
    
    
    cout << qr_economic(A).first << endl;
    cout << qr_economic(A).second << endl;
    return 0;
}
