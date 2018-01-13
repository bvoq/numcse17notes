#include <iostream>
#include <vector>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

// Point (c)

VectorXd solve_general_blockLU(MatrixXd & A, MatrixXd & B, MatrixXd & C, MatrixXd & D, VectorXd & b) {
    size_t n = A.rows(), m = C.rows();
    //    [ I         0 ]      [ R              v ]
    //L = [ C*A^-1    I ]  U = [ 0    -u^T R^-1 v ]
    //Solve L c = b
    VectorXd c (n+m);
    c.head(n) = b.head(n);
    //Fancy idea: C*A^-1*c(0..n-1) + 1 * c(n..n+m-1) = b(n..n+m-1)
    VectorXd z = A.fullPivLu().solve(c.head(n)); //A^-1 * c(0..n-1)
    c.tail(m) = b.tail(m)-C*z;
    
    MatrixXd S = D-C*A.fullPivLu().solve(B);
    MatrixXd U (n+m,n+m);
    MatrixXd Z = MatrixXd::Zero(m,n);
    U << A, B,
    Z, S;
    VectorXd x = U.fullPivLu().solve(c);
    //IF A TRIANGULAR: U.triangularView<Upper>().solve(c);
    return x;
}

void solve_blockLU(const MatrixXd &R,
                   const VectorXd &u,
                   const VectorXd &v,
                   const VectorXd &b,
                   VectorXd &x)
{
    //The matrix looks like:
    //    [ I         0 ]      [ R              v ]
    //L = [ u^T*R^-1  I ]  U = [ 0    -u^T R^-1 v ]
    
    size_t n = R.rows(), m = 1;
    
    //Solve L c = b
    VectorXd c (n+m);
    c.head(n) = b.head(n);
    // Equation to solve: u.transpose() * R^-1 * c(0..n-1) + 1 * x = b(n..n+m-1);
    // R^-1 * c(0..n-1) = z <=> c(0..n-1) = R * z
    // u.transpose() * z + 1 * x(n..n+m-1) = b(n..n+m-1)
    VectorXd z = R.triangularView<Upper>().solve(c.head(n));
    c.tail(m) = b.tail(m)-u.transpose()*z;
    
    //Solve U x = c
    MatrixXd U (n+m,n+m);
    U <<    R,                              v,
    VectorXd::Zero(n).transpose(),  -u.transpose()*R.triangularView<Upper>().solve(v);
    
    //Since U is triangular this should go faster, but hey:
    x = U.triangularView<Upper>().solve(c);
    //c(n+1) = b(n+1) / (u.transpose() * c.head(n));
    
}

// Point (f)

void solve_blockGauss(const MatrixXd &R,
                      const VectorXd &u,
                      const VectorXd &v,
                      const VectorXd &b,
                      VectorXd &x) {
    
    int n = R.rows();
    
    const TriangularView<const MatrixXd, Upper> & triR = R.triangularView<Upper>();
    
    double sinv = - 1. / u.dot(triR.solve(v));
    double bs = (b(n) - u.dot(triR.solve(b.head(n))));
    double sinvbs = sinv*bs;
    
    x << triR.solve(b.head(n) - v*sinvbs),
    sinvbs;
}

int main() {
    int N = 5;
    MatrixXd R = MatrixXd::Random(N,N);
    R = R.triangularView<Upper>();
    VectorXd u = VectorXd::Random(N);
    VectorXd v = VectorXd::Random(N);
    VectorXd b = VectorXd::Random(N+1);
    VectorXd x1 = VectorXd::Random(N+1);
    VectorXd x2 = VectorXd::Random(N+1);
    
    solve_blockLU(R, u, v, b, x1);
    solve_blockGauss(R, u, v, b, x2);
    
    MatrixXd A(N+1,N+1);
    A << R,                v,
    u.transpose(), 0;
    
    cout << A.fullPivLu().solve(b) << "\n\n";
    cout << x1 << "\n\n";
    cout << x2 << "\n\n";
    
    
    
    MatrixXd kunen = MatrixXd::Random(15,15);
    MatrixXd bunen = MatrixXd::Random(15,15);
    MatrixXd yo = kunen.inverse()*bunen;
    MatrixXd yo2 = kunen.fullPivLu().solve(bunen);
    cout << yo << "\n---------\n" << yo2 << endl;
    cout << "Norm: " << yo.norm() << " norm2: " << yo2.norm() << endl;
    
    
    /*
    for(int i = 0; i < 10000; ++i) {
        MatrixXd A = MatrixXd::Random(15,15);
        MatrixXd B = MatrixXd::Random(15,10);
        MatrixXd C = MatrixXd::Random(10,15);
        MatrixXd D = MatrixXd::Random(10,10);
        VectorXd b = VectorXd::Random(15+10);
        MatrixXd M(15+10,15+10); M << A,B,C,D;
        if(M.fullPivLu().isInvertible()) {
            VectorXd x = solve_general_blockLU(A, B, C, D, b);
            VectorXd x2 = M.fullPivLu().solve(b);
            //cout << x << "\n-------\n" << x2 << endl;
            assert(fabs(x.norm() - x2.norm() < 0.000001));
        }
    }*/
}
