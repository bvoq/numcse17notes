#include <iostream>
#include <Eigen/Eigen>
using namespace Eigen;
using namespace std;

/*
 void randomMatrixAsComparison() {
 int n = 3;
 MatrixXd A = MatrixXd::Random(n,n);
 
 auto svd = A.jacobiSvd(ComputeThinU | ComputeThinV);
 cout << "The rank of A (computed using jacobi SVD is): " << svd.rank() << endl;
 cout << "The matrix A is: " << endl << A << endl << endl;
 
 auto LU = ((MatrixXd)A.fullPivLu()).eval();
 cout << "LU: " << endl;
 cout << (MatrixXd)LU << endl;
 }*/


void LUdecLeadsToTinyNumbers() {
    for(int n = 5; n < 15; ++n) {
        MatrixXd A = MatrixXd::Zero(n,n);
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < i; ++j) {
                A(i,j) = -1;
            }
            A(i,i) = 1;
            A(i,n-1) = 1;
        }
        
        auto svd = A.jacobiSvd(ComputeThinU | ComputeThinV);
        auto singularvalues = svd.singularValues();
        double norm = singularvalues(0) / singularvalues(singularvalues.size()-1);
        cout << "The rank of A (computed using jacobi SVD is): " << svd.rank() << endl;
        cout << "Condition number of A is " << norm << " ~ " << (float)pow(2,n+1) << endl;
        //cout << "The matrix A is: " << endl << A << endl << endl;
        
        auto LU = ((MatrixXd)A.partialPivLu().matrixLU()).eval();
        cout << "LU: " << endl;
        cout << (MatrixXd)LU << endl;
    }
}



int main() {
    LUdecLeadsToTinyNumbers();
    
    //cout << "A solution is: " << endl << endl;
    
    
    return 0;
}
