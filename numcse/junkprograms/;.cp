#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
using namespace Eigen;
using namespace std;

int main() {
    Triplet<double> a (4,4,2.5);
    SparseMatrix<double> spm(6,4);
    spm.setIdentity();
    cout << spm << endl;
    
    return 0;
}
