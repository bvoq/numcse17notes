#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "timer.h"


/*
 * @brief Evaluate a polynomial and its derivative using Horner scheme
 * @param[in] vector c of size $n$, coefficients of the polynomial p
 * @param[in] double x, where the polynomial has to be evaluated
 * @param[out] vector containing p(x),p'(x)
 */
template <typename CoeffVec>
Eigen::Vector2d dpEvalHorner (const CoeffVec& c, const double x) {
  Eigen::Vector2d val;
  double px, dpx;
  int s = c.size();

  px = c[0];
  for (int i = 1; i < s; ++i) { px = x*px+c[i]; }
    
  dpx = (s-1)*c[0];
  for (int i = 1; i < s-1; ++i) { dpx = x*dpx+(s-i-1)*c[i]; }
    
  val(0) = px;
  val(1) = dpx;

  return val;
}

/*
 * @brief Evaluate a polynomial and its derivative using a naive implementation
 * @param[in] vector c of size $n$, coefficients of the polynomial p
 * @param[in] double x, where the polynomial has to be evaluated
 * @param[out] vector containing p(x),p'(x)
 */
template <typename CoeffVec>
Eigen::Vector2d dpEvalNaive(const CoeffVec& c, const double x) {
  Eigen::Vector2d val;
  double px,dpx;
  int n=c.size();
    
  px = c[0]*std::pow(x, n-1);
  for (int i = 1; i < n; ++i) { px = px + c[i]*std::pow(x, n-i-1); }
    
  dpx = (n-1)*c[0]*std::pow(x, n-2);
  for (int i = 1; i < n-1; ++i) {dpx = dpx + (n-i-1)*c[i]*std::pow(x, n-i-2); }
    
  val(0) = px;
  val(1) = dpx;

  return val;
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
    
  // Compare runtimes
  const unsigned repeats = 10;
   
  std::cout << std::setw(10) << "n" << std::setw(25) << "Horner scheme:" 
            << std::setw(25) << "Monomial approach:" << "\n"
            << " ================================================================\n";

  std::vector<double> e, h, m;
  for (unsigned k = 2; k <= 20; ++k) {
    Timer tm_slow, tm_fast;
    std::vector<double> c;
        
    const int n = std::pow(2, k);
    for (int i = 0; i < n; ++i) { c.push_back(i+1); }

    for (unsigned r = 0; r < repeats; ++r) {
      tm_slow.start();
      p_naive = dpEvalNaive(c,x);
      tm_slow.stop();
                
      tm_fast.start();
      p = dpEvalHorner(c,x);
      tm_fast.stop();
    }
    
    std::cout << std::setw(10) << n << std::setw(25) << tm_fast.min() 
              << std::setw(25) << tm_slow.min() << "\n";
  }
}
