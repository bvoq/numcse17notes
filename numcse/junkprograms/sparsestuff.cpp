#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
using namespace Eigen;
using namespace std;

int main() {
    Triplet<double> a (4,4,2.5);
    SparseMatrix<double> spm(6,6);
    spm.reserve(2);
    spm.setIdentity();
    spm.coeffRef(2, 4) = 3;
    spm.makeCompressed();
    vector<Triplet<double> > yo;
    yo.push_back(Triplet<double>(4,3,2.5));
    spm.setFromTriplets(yo.begin(), yo.end()); //overwrites previous spm
    
    
    MatrixXd denseM (4,3);
    denseM <<
    1,2,3,
    0,0,0,
    0,0,0,
    0,4,0;
    
    
    SparseMatrix<double> spm2 = denseM.sparseView();
    cout << spm2 << endl;
    
    
    MatrixXd q(3,3);
    q << 1,-1,-1,   -1,1,0,    1,0,1;
    cout << q.inverse() << endl;
    
    MatrixXd u(3,3);
    u << 1,2,3,    0,4,5,    0,0,6;
    cout << u.transpose() << endl;
    
    return 0;
}
