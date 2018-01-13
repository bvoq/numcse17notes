#include <iostream>
#include <cmath>
#include <ctime>
#include <random>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <unsupported/Eigen/KroneckerProduct>

#include "timer.h"

using namespace Eigen;

using index_t = int;

/* @brief Compute $B = \argmin_{X \in P(A)} |I-AX|_F$
 * @param[in] A An $n \times n$ matrix
 * @param[out] B The $n \times n$ matrix $= \argmin_{X \in P(A)} |I-AX|_F$
 */
SparseMatrix<double> spai(SparseMatrix<double> & A) {
    // Size check
    assert(A.rows() == A.cols() &&
    "Matrix must be square!");
    unsigned int N = A.rows();

    A.makeCompressed();

    // Obtain pointers to data of A
    double* valPtr = A.valuePtr();
    index_t* innPtr = A.innerIndexPtr();
    index_t* outPtr = A.outerIndexPtr();
    
    
    // Create vector for triplets of B and reserve enough space
    std::vector<Triplet<double>> B_triplets;
    B_triplets.reserve(A.nonZeros());

    // Project $\mathcal{P}(A)$ onto $\mathcal{P}(a_i)$ and compute $b_i$
    for(unsigned int i = 0; i < N; ++i) {
        // Number of non-zeros in $a_i$
        index_t nnz_i = outPtr[i+1] - outPtr[i];
        if(nnz_i == 0) continue; // skip column, if empty

        // Smaller and denser matrix to store
        // non-zero elements relevant for computing $b_i$ 
        SparseMatrix<double> C(N, nnz_i);
        std::vector<Triplet<double>> C_triplets;
        C_triplets.reserve(nnz_i*nnz_i);

        // Build matrix $C$.
        for(unsigned int k = outPtr[i]; k < outPtr[i+1]; ++k) {
            // Row index of non-zero element in column $b_i$
            index_t row_k = innPtr[k];
            // Number of non-zero entries in ($row_k$)-th column
            index_t nnz_k = outPtr[row_k+1] - outPtr[row_k];
            // Loop over all non-zeros of $row_k$th-column
            // Store the triplet for the non-zero element 
            for(unsigned int l = 0; l < nnz_k; ++l) {
                unsigned int innIdx = outPtr[row_k] + l;
                C_triplets.emplace_back(Triplet<double>(innPtr[innIdx], k - outPtr[i], valPtr[innIdx]));
            }
        }
        C.setFromTriplets(C_triplets.begin(), C_triplets.end());
        C.makeCompressed();
        
        // Normal equation method: $b_i = (C^{\top} C)^{-1} C^{\top} e_i$ 
        SparseMatrix<double> S = C.transpose() * C;
        VectorXd xt = C.row(i).transpose();
        SparseLU<SparseMatrix<double>> spLU(S);
        VectorXd b = spLU.solve(xt);
        
        // store the triplets for elements $b_i$
        for(unsigned int k = 0; k < b.size(); ++k) {
            B_triplets.emplace_back(Triplet<double>(innPtr[outPtr[i] + k], i, b(k)));
        }
    }

    // Build and return $B$
    SparseMatrix<double> B = SparseMatrix<double>(N,N);
    B.setFromTriplets(B_triplets.begin(), B_triplets.end());
    B.makeCompressed();
    return B;
}


// Run conditionally only one test
const bool small_test = true;
const bool big_test = false;

int main(int argc, char **argv) {
    unsigned int n = 60;
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
    if(big_test)
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
    }
    
    return 1;
}
