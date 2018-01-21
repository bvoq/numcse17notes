#include <iostream>
#include <Eigen/Eigen>
using namespace Eigen; using namespace std;

int main() {
    MatrixXd M(3,3); M << 1,2,3, 4,5,6, 7,8,9;
    VectorXd q (3); q << 5, 9, 7;
    VectorXd x (3); x << -0.666667, 1.33333, 0;
    cout << M.fullPivLu().rank() << endl;
    
    cout << M.lu().solve(q) << endl;
    cout << M * x << endl;
    cout << M.lu().matrixLU().eval() << endl;
    
    auto LU = M.lu();
    cout << LU.solve(q) << endl;
    cout << LU.solve(x) << endl;
    auto QR = M.colPivHouseholderQr();
    
    
    return 0;
}
