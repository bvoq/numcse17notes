#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "timer.h"

using namespace Eigen;

/*
 * @brief Evaluate a polynomial and its derivative using Horner scheme
 * @param[in] vector c of size $n$, coefficients of the polynomial p
 * @param[in] double x, where the polynomial has to be evaluated
 * @param[out] vector containing p(x),p'(x)
 */
template <typename CoeffVec>
Eigen::Vector2d dpEvalHorner (const CoeffVec& c, const double x) {
    int n = c.size()-1;
    double val = c[n];
    for(int i = n-1; i >= 0; --i) {
        val *= x;
        val += c[i];
    }
    double dval = c[n] * n;
    for(int i = n-1; i >= 1; --i) {
        dval *= x;
        dval += c[i] * i;
    }
    Vector2d out; out << val, dval;
    return out;
}

/*
 * @brief Evaluate a polynomial and its derivative using a naive implementation
 * @param[in] vector c of size $n$, coefficients of the polynomial p
 * @param[in] double x, where the polynomial has to be evaluated
 * @param[out] vector containing p(x),p'(x)
 */
template <typename CoeffVec>
Eigen::Vector2d dpEvalNaive(const CoeffVec& c, const double x) {
  
    double val = 0, dval = 0;
    for(int i = 0; i < c.size(); ++i) val += c[i] * pow(x,i);
    for(int i = 1; i < c.size(); ++i) dval += i * c[i] * pow(x,i-1);
    Vector2d out; out << val, dval;
    return out;
}

int main() {
  std::vector<double> c {3, 1, 5, 7, 9};
  double x = .123;
    
  // Check implementations
  Eigen::Vector2d p, p_naive;
  p = dpEvalHorner(c,x);
  std::cout << "Using horner scheme:\n"
            << "p(x) = " << p(0) << 
            ", dp(x) = " << p(1) << "\n";
    
  p_naive = dpEvalNaive(c,x);
  std::cout << "Using monomial approach:\n"
            << "p(x) = " << p_naive(0) 
            << ", dp(x) = " << p_naive(1) << "\n";

 
}
