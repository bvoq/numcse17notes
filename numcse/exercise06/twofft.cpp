#include <iostream>
#include <Eigen/Eigen>
#include <unsupported/Eigen/FFT>
#include <mgl2/mgl.h>
#include <mgl2/wnd.h>
#include <mgl2/fltk.h>

using namespace Eigen;
using namespace std;

MatrixXcd fft2(MatrixXcd in) {
    in.eval(); //make sure the representation of in is actually written.
    MatrixXcd temp(in.rows(), in.cols()), out(in.rows(), in.cols());
    FFT<double> fft;
    for(size_t i = 0; i < in.rows(); ++i) temp.row(i) = (VectorXcd)fft.fwd((VectorXcd)in.row(i));
    for(size_t i = 0; i < in.cols(); ++i) out.col(i) = (VectorXcd) fft.fwd((VectorXcd)temp.col(i));
    return out;
}

MatrixXcd ifft2(MatrixXcd in) {
    return 1./(in.rows() * in.cols()) * fft2(in.conjugate()).conjugate();
}


MatrixXd periodicconvol2(MatrixXd X, MatrixXd Y) {
    assert(X.rows() == Y.rows() && X.cols() == Y.cols());
    size_t m = X.rows(), n = X.cols();
    MatrixXd Z(m,n);
    for(size_t k = 0; k < m; ++k) for(size_t l = 0; l < n; ++l) for(size_t i = 0; i < m; ++i) for(size_t j = 0; j < n; ++j) {
        Z(k,l) = X(i,j) * Y((k-i)%m, (l-j)%n);
    }
    return Z;
}

MatrixXd discreteconvol2(MatrixXd A, MatrixXd B) {
    size_t Lr = A.rows() + B.rows() - 1;
    size_t Lc = A.cols() + B.cols() - 1;
    MatrixXd Az = MatrixXd::Zero(Lr, Lc);
    MatrixXd Bz = MatrixXd::Zero(Lr, Lc);
    Az.block(0,0,A.rows(),A.cols()) = A;
    Bz.block(0,0,B.rows(),B.cols()) = B;
    return periodicconvol2(Az, Bz);
}


int main() {
    
    MatrixXcd m (3,2);
    m << 2,3,4,5,-1,complex<double>(6,4);
    
    cout << m << endl << "-------" << endl << endl;
    
    MatrixXcd res = ifft2(m);
    cout << endl<< endl <<res << endl;
    return 0;
}


