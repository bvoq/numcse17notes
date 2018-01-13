#include <iostream>
#include <vector>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;
std::string sep = "\n---------------\n";

int signOfPermMatrix(MatrixXd P) {
    vector<int> pos(P.rows());
    for(int i = 0; i < P.rows(); ++i) {
        for(int j = 0; j < P.cols(); ++j) {
            if(P(i,j) == 1) pos[i] = j;
        }
    }

    int swapCount = 0;
    for(int c = 0; c < pos.size(); ++c) {
        for(int i = 0; i < pos.size(); ++i) {
            if(pos[i] == c) {
                int j = i;
                swapCount += abs(c-i);
                while(j < c) swap(pos[j],pos[++j]);
                while(j > c) swap(pos[j],pos[--j]);
            }
        }
    }
    return swapCount % 2 == 0 ? 1 : -1;
}

double determinant(MatrixXd A) {
    long long n = A.rows();
    assert(A.rows() == A.cols());
    auto PLU = A.partialPivLu();
    auto LU = PLU.matrixLU();
    auto P = PLU.permutationP(); //inverse
    
    double prod = 1;
    for(long long i = 0; i < n; ++i) prod *= LU(i,i);
    prod *= signOfPermMatrix(P);
    //cout << "sign " << signOfPermMatrix(P) << " " << P.determinant() << endl;
    assert(signOfPermMatrix(P) == P.determinant());
    assert(fabs(prod - A.determinant()) <0.00001);
    return prod;
}

int main() {
    MatrixXd N (4,4); N << 1,2,6,4,2,3,4,5,6,4,5,6,4,6,6,7;
    
    //auto PLU = N.partialPivLu();
    //auto LU = PLU.matrixLU();
    auto PiLUQi = N.fullPivLu();
    auto Pi = (MatrixXd) PiLUQi.permutationP();
    auto Qi = (MatrixXd) PiLUQi.permutationQ();
    auto LU = PiLUQi.matrixLU();
    
    cout << "Sign: " << signOfPermMatrix(Pi) << endl;

    
    int n = N.rows();
    MatrixXd seperatedL = MatrixXd::Zero(n, n);
    MatrixXd seperatedU = MatrixXd::Zero(n, n);
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < i; ++j) seperatedL(i,j) = LU(i,j);
        for(int j = i; j < n; ++j) seperatedU(i,j) = LU(i,j);
        seperatedL(i,i) = 1;
    }
    assert(PiLUQi.isInvertible());
    assert(PiLUQi.inverse().norm() == N.inverse().norm());
    assert(determinant(N) == N.determinant());
    assert(N.determinant() != 0);
    assert((Pi.inverse() * seperatedL * seperatedU * Qi.inverse()).norm() == N.norm());
    VectorXd b (n); b << 0, 1, 0, 1;
    VectorXd x = PiLUQi.solve(b);
    double relative_error = (N * x - b).norm() / b.norm();
    cout << "Relative error: " << relative_error << endl;
    
    /*
    for(int i = 0; i < 1000000; ++i) {
        MatrixXd T = MatrixXd::Random(9, 9);
        if(i % 1000 == 0) cout << i << " DET: " << endl << determinant(T) << endl;
    }
    */
    
    return 0;
    //add transpose on .eval()
}
