#include <Eigen/Eigen>
#include <iostream>
using namespace Eigen;
using namespace std;

int main() {
    int n = 10;
    MatrixXd M = MatrixXd::Zero(n,n);
    for(int i = 0; i < n; ++i) {
        M(i,i) = 2;
        for(int j = 0; j < i; ++j) M(i,j) = 1;
    }
    
    auto LU = M.lu();
    
    
    return 0;
}
