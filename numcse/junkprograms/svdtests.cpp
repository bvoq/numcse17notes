#include <iostream>
#include <Eigen/Eigen>

using namespace Eigen;
using namespace std;


double frobeniusNorm(const MatrixXd & M) {
    double accum = 0;
    for(int i = 0; i < M.size(); ++i) accum += M(i) * M(i);
    return sqrt(accum);
    
}
int main() {
    MatrixXd A = MatrixXd::Random(5, 5);
    cout << "A: " << endl << A << endl << endl;
    auto svd = A.jacobiSvd(ComputeFullU | ComputeFullV);
    auto U = svd.matrixU();
    auto V = svd.matrixV();
    
    auto S = svd.singularValues();
    MatrixXd Sigma = MatrixXd::Zero(S.size(), S.size());
    Sigma.diagonal() << S;
    cout << "SIGMA: " << endl << Sigma << endl;
    
    cout << "NORM A; 2-norm: " << A.squaredNorm() << " Frobenius norm: " << frobeniusNorm(A) <<endl;
    cout << "NORM Sig; 2-norm: " << Sigma.squaredNorm() << " Frobenius norm: "  << frobeniusNorm(Sigma) << endl;
    cout << "Computed 2-norm: " << S(0) / S(S.size()-1) << endl;
    cout << "errr: " << S(0) << endl;

    MatrixXd m(2,2); m << 1,-2,  -3,4;
    cout << "feige: " << m.squaredNorm() << endl;
    cout << "sing values: " << m.jacobiSvd(ComputeFullU|ComputeFullV).singularValues() << endl;
    
    MatrixXd K (4,4); K << 1,-2,2,4,  2,-4,5,9,  3,-6,8,14,  5,-10,12,22;
    auto svd2 = K.jacobiSvd(ComputeFullU | ComputeFullV);
    cout << "yo: " << svd2.rank() << endl;
    int r = svd2.rank();
    MatrixXd U2 = svd2.matrixU();
    MatrixXd V2 = svd2.matrixV();
    MatrixXd nSp =  V2.block(0,r,4,4-r);
    cout << nSp << endl;
    
    return 0;
}
