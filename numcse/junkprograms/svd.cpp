#include <iostream>
#include <tuple>
#include <Eigen/Dense>
#include <Eigen/SVD>
using namespace Eigen;
using namespace std;

tuple<MatrixXd,MatrixXd,MatrixXd> svd_full(const MatrixXd & A) {
    JacobiSVD<MatrixXd> svd(A, ComputeFullU | ComputeFullV);
    MatrixXd U = svd.matrixU(); //unitary square
    MatrixXd V = svd.matrixV(); //unitary square
    VectorXd sv = svd.singularValues(); //including zero singular values
    cout << "Size is: " << sv.size() << endl;
    MatrixXd Sigma = MatrixXd::Zero(A.rows(), A.cols());
    for(int i = 0; i < sv.size(); ++i) Sigma(i,i) = sv(i);
    
    assert(((U*Sigma*V.transpose()).norm() - A.norm()) < numeric_limits<double>::epsilon());
    return {U,Sigma,V};
}

tuple<MatrixXd,MatrixXd,MatrixXd> svd_partial(const MatrixXd & A) {
    JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);
    MatrixXd U = svd.matrixU(); //not-necessairly square but column-wise orthonormal.
    MatrixXd V = svd.matrixV(); //not-necessairly square but column-wise orthonormal.
    MatrixXd Sigma = svd.singularValues().asDiagonal(); //including zero singular values
    return {U,Sigma,V};
}

int main() {
    MatrixXd A (4,3);
    A << 1,2,4,    2,4,8,    3,6,12,      3,4,1;
    auto tp = svd_full(A);
    cout << "U: "<<endl << get<0>(tp) << endl; //
    cout << "Sigma:" << endl << get<1>(tp) << endl; //
    cout << "V: " << endl << get<2>(tp) << endl; //
    
    cout << A << endl;
    cout << endl << "almost equal to : " << endl << get<0>(tp)*get<1>(tp)*get<2>(tp).transpose() << endl;
    
    tp = svd_partial(A);
    cout << "Up: " << endl << get<0>(tp) << endl; //
    cout << "Sigmap: " << endl << get<1>(tp) << endl; //
    cout << "Vp: " << endl << get<2>(tp) << endl;
    
    //U and V are unitary
    assert(get<0>(tp).col(0).transpose() * get<0>(tp).col(1)  < numeric_limits<double>::epsilon());
    assert(get<0>(tp).col(0).transpose() * get<0>(tp).col(2)  < numeric_limits<double>::epsilon());
    assert(get<0>(tp).col(1).transpose() * get<0>(tp).col(2)  < numeric_limits<double>::epsilon());
    
    cout << A << endl;
    cout << endl << "almost equal to : " << endl << get<0>(tp)*get<1>(tp)*get<2>(tp).transpose() << endl;
    
   
}
