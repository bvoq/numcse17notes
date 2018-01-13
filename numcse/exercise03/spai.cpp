#include <iostream>
#include <Eigen/Eigen>
using namespace Eigen;
using namespace std;
SparseMatrix<double> spai(SparseMatrix<double> & A) {
    
    SparseMatrix<double> X (A.rows(), A.cols());
    vector<Triplet<double> > Xtrip;
    
    //Assumes A is in colmajor form
    
    double* vp = A.valuePtr();
    int * iip = A.innerIndexPtr();
    int * oip = A.outerIndexPtr();
    for(int j = 0; j < A.nonZeros(); ++j) cout << iip[j] << " ";
    cout << endl;
    int nnz = A.nonZeros();
    
    cout << "this is A: " << endl << A << endl;
    
    //argmin ||e_i - A x_i||
    cout << endl;
    

    for(int i = 0; i < A.cols(); ++i) { //for each column in A do:
        int nnzincoli = oip[i+1]-oip[i];
        SparseMatrix<double> D(A.rows(), nnzincoli); // the D matrix, should contain one more A though
        vector<Triplet<double > > Dtrip;
        cout << "i=" << i << " from oip[i]="<<oip[i] << " to oip[i+1]="<<oip[i+1] << endl;
        for(int j = oip[i]; j < oip[i+1]; ++j) { //loop over nnz column entries of i-th column in A
            int col = i;
            int row = iip[j];
            //now insert all nnz of row into C:
            for(int k = oip[row]; k < oip[row+1]; ++k) { //add all nnz elements of the k-th column of A to D
                cout << "inserting: " << iip[k] << " " << j-oip[i] << " with val " << -vp[k]  << "     j: " << j << endl;
                Dtrip.push_back(Triplet<double>(iip[k], j-oip[i], vp[k]));
            }
            //row = iip[j];
            //col = i;
        }
        D.setFromTriplets(Dtrip.begin(), Dtrip.end());
        D.makeCompressed();
        
        /*
        SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > spq(D);
        VectorXd ei = VectorXd::Zero(A.rows());
        ei(i) = 1;
        VectorXd x = spq.solve(ei);
        */
        
       
        SparseMatrix<double> S = D.transpose() * D;
        VectorXd xt = D.row(i).transpose();
        SparseLU<SparseMatrix<double>> spLU(S);
        VectorXd b = spLU.solve(xt);
        VectorXd x = b;
        //reconstruct actual x by adding zero elements to the column vector at the appropriate places
        //VectorXd xn = VectorXd::Zero(A.cols());
        //for(int j = oip[i]; j < oip[i+1]; ++j) {
        //    xn(iip[j]) = x(j-oip[i]);
        //}
        for(int k = oip[i]; k < oip[i+1]; ++k) {
            Xtrip.emplace_back(Triplet<double>(iip[k],i,b(k-oip[i])));
        }
        //directly add solution tripplets to X
        /*for(int j = oip[i]; j < oip[i+1]; ++j) {
            cout << i << ", " << iip[j] << ", ";
            cout << "x(" << j-oip[i] << ") = " << x(j-oip[i]) << endl;
            Xtrip.push_back(Triplet<double>(iip[j],i,x(j-oip[i])));
        }*/
        cout << "nnzincoli(" << i << ") = " << nnzincoli << endl;
    }

    X.setFromTriplets(Xtrip.begin(), Xtrip.end());
    X.makeCompressed();
    cout << endl << "X:" << endl << X << endl;
    //Least squares for sparse matrices:
    return X;
}

bool small_test = true, big_test = false;
int main(int argc, char **argv) {
    
    signed int n = 60;
    if(argc > 1) {
        n = std::stoi(argv[1]);
    }
    srand(time(NULL));
    
    // First test: compute SPAI with very small matrix
    // Check correctness
    if(small_test)
    {
        n = 5;
        SparseMatrix<double> M(5,5);
        M.coeffRef(3,4) = 1;
        M.coeffRef(4,3) = 2;
        M.coeffRef(1,4) = 3;
        M.coeffRef(3,3) = 4;
        M.coeffRef(3,2) = 4;
        M.coeffRef(2,3) = 4;
        M.coeffRef(2,2) = 5;
        M.coeffRef(3,1) = 6;
        M.coeffRef(0,0) = 9;
        M.makeCompressed();
        SparseMatrix<double> N = spai(M);
        SparseMatrix<double>I(n,n);
        I.setIdentity();
        
        std::cout << "M = " << std::endl
        << M
        << std::endl;
        std::cout << "N = " << std::endl
        << N
        << std::endl;
        
        //std::cout << "Error:"
        //          << (I - M*N).norm()
        //          << std::endl;
        
    }
    // Big test: test with large, sparse matrix
    /*if(big_test)
    {
        SparseMatrix<double> M(n*n,n*n);
        
        SparseMatrix<double>I(n,n);
        I.setIdentity();
        MatrixXd R = MatrixXd::Random(n,n);
        M = kroneckerProduct(R, I);
        
        Timer tm;
        tm.start();
        SparseMatrix<double> N = spai(M);
        tm.stop();
        
        SparseMatrix<double>Ibig(n*n,n*n);
        Ibig.setIdentity();
        
        std::cout << "Error (n = " << n*n << "):"
        << (Ibig - M*N).norm()
        << std::endl
        << "Elapsed:"
        << tm.duration() << " s"
        << std::endl;
    }*/
    
    return 1;
}
/*
int main() {
    SparseMatrix<double> spm (4,3);
    vector<Triplet<double> > trip;
    trip.push_back(Triplet<double>(0,0,1));
    trip.push_back(Triplet<double>(0,1,2));
    trip.push_back(Triplet<double>(0,2,3));
    trip.push_back(Triplet<double>(3,1,4));
    spm.setFromTriplets(trip.begin(), trip.end());
    spai(spm);
}
*/
