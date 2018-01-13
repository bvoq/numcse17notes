//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>
#include <iomanip>

#include <Eigen/Dense>
#include <mgl2/mgl.h>
#include "timer.h"

using namespace std;
using namespace Eigen;

/* @brief Build an "arrow matrix" and compute A*A*y
 * Given vectors $a$ and $d$, returns A*A*x in $y$, where A is built from a, d
 * @param[in] d An n-dimensional vector
 * @param[in] a An n-dimensional vector
 * @param[in] x An n-dimensional vector
 * @param[out] y The vector y = A*A*x
 */
void arrow_matrix_2_times_x(const VectorXd &d, const VectorXd &a,
                            const VectorXd &x, VectorXd &y) {
    assert(d.size() == a.size() && a.size() == x.size() &&
           "Vector size must be the same!");
    int n = d.size();

    VectorXd d_head = d.head(n-1);
    VectorXd a_head = a.head(n-1);
    MatrixXd d_diag = d_head.asDiagonal();

    MatrixXd A(n,n);

    A << d_diag,             a_head,
         a_head.transpose(), d(n-1);

    y = A*A*x;
}

/* A SKETCHd:
d1 0  0  a1
0  d2 0  a2
:  0  \\ :
a1 a2 .. dn
*/

/* @brief Build an "arrow matrix"
 * Given vectors $a$ and $b$, returns A*A*x in $y$, where A is build from a,d
 * @param[in] d An n-dimensional vector
 * @param[in] a An n-dimensional vector
 * @param[in] x An n-dimensional vector
 * @param[out] y The vector y = A*A*x
 */
void efficient_arrow_matrix_2_times_x(const VectorXd &d,
                                      const VectorXd &a,
                                      const VectorXd &x,
                                      VectorXd &y) {
    assert(d.size() == a.size() && a.size() == x.size() &&
         "Vector size must be the same!");
    int n = d.size();
    VectorXd t = VectorXd(n); //=A*x
    for(int i = 0; i < n-1; ++i) t(i) = d(i)*x(i) + a(i)*x(n-1);
    t(n-1)=0;
    for(int i = 0; i < n-1; ++i) t(n-1) += a(i)*x(i);
    t(n-1) += d(n-1)*x(n-1);
    y = VectorXd(n); //will the vector by default be zero initialised?
    for(int i = 0; i < n-1; ++i) y(i) = d(i)*t(i) + a(i)*t(n-1);
    y(n-1)=0;
    for(int i = 0; i < n-1; ++i) y(n-1) += a(i)*t(i);
    y(n-1) += d(n-1)*t(n-1);
}

/* \brief Compute the runtime of arrow matrix multiplication.
 * Repeat tests 10 times, and output the minimal runtime
 * amongst all times. Test both the inefficient and the efficient
 * versions.
*/
#include <chrono>
#include "timer.h"
void runtime_arrow_matrix() {
    // TODO: your code here, time the codes
    srand(time(0));
    std::chrono::time_point<std::chrono::system_clock> start, end;
    Timer t;
    long long startms, endms;
    for(int i = 1; i <= 10; ++i) {
        int n = pow(2,i);
        cout << "Computing for n=" << n << "="<<"2^"<<i<<"." << endl;
        VectorXd d(n), a(n), x(n);
        for(int i = 0; i < n; ++i) d(i) = 100.*rand()/(RAND_MAX-1);
        for(int i = 0; i < n; ++i) a(i) = 100.*rand()/(RAND_MAX-1);
        for(int i = 0; i < n; ++i) x(i) = 100.*rand()/(RAND_MAX-1);
        VectorXd y1, y2;
        t.reset();
        t.start();
        start = std::chrono::system_clock::now();
        startms = time(0);
        arrow_matrix_2_times_x(d, a, x, y1);
        end = std::chrono::system_clock::now();
        endms = time(0);
        t.stop();

        cout << "Inefficient took " << (end-start).count() << " ticks and ~ " << endms - startms << "s and " << t.duration() << "s" << endl;
        start = std::chrono::system_clock::now();
        t.reset();
        t.start();
        startms = time(0);
        efficient_arrow_matrix_2_times_x(d, a, x, y2);
        end = std::chrono::system_clock::now();
        endms = time(0);
        t.stop();
        cout << "Inefficient took " << (end-start).count() << " ticks and ~ " << endms - startms << "s and " << t.duration() << "s" << endl;

    }
}


int main(void) {
    // Test vectors
    VectorXd a(5);
    a << 1., 2., 3., 4., 5.;
    VectorXd d(5);
    d <<1., 3., 4., 5., 6.;
    VectorXd x(5);
    x << -5., 4., 6., -8., 5.;
    VectorXd yi;
    // Run both functions
    arrow_matrix_2_times_x(a,d,x,yi);
    VectorXd ye(yi.size());
    efficient_arrow_matrix_2_times_x(a,d,x,ye);

    // Compute error
    double err = (yi - ye).norm();

    // Output error
    std::cout << "--> Correctness test." << std::endl;
    std::cout << "Error: " << err << std::endl;

    // Print out runtime
    std::cout << "--> Runtime test." << std::endl;
    runtime_arrow_matrix();

    // Final test: exit with error if error is too big
    double eps = std::numeric_limits<double>::denorm_min();
    exit(err < eps);
}
