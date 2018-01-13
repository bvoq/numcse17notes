#include <iostream>
#include <Eigen/Eigen>

using namespace Eigen;
using namespace std;
int main() {
    MatrixXd A (3,3);
    A << 1,2,3,   4,4,4,    7,7,5;
    auto svd = A.jacobiSvd(ComputeFullU | ComputeFullV);
    auto sv = svd.singularValues();
    cout << "Rank A: " << svd.rank() << endl;
    
    cout << "Norm of A: " << sv(0) / sv(sv.size()-1) << endl;
    MatrixXd AI = A.inverse();
    cout << "More than the norm: " << A.norm() << " " << AI.norm() << " |A|*|A^-1|" << A.norm() * AI.norm() << endl;
    
    
    return 0;
}
