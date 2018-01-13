//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>

#include <Eigen/Eigen>
#include <unsupported/Eigen/FFT>
#include <mgl2/mgl.h>

using namespace Eigen;

/*!
 * \brief triPolyFit Find best trigonometric polynomial
 * passing trough distances $\mathbf{d}$.
 * Using Least Squares formulation and assuming best fitting
 * polynomial has degree $m$.
 * \param d Vector of size $n$ distances at angle $2*\pi*i/n$.
 * \param m Degree of the trigonometric polynomial $p$.
 * \return The coefficients of the trigonometric polinomial.
 */
VectorXd triPolyFit(const VectorXd & d, unsigned int m) {
/*
    // We will use a real to complex, discrete Fourier transform.
    FFT<double> fft;
    
    VectorXcd fft_d = fft.fwd(d);
    VectorXcd ifft_d = fft.inv(d);
    VectorXd lul;
    return lul;
    */
    unsigned int n = d.size();
    
    // We will use a real to complex, discrete Fourier transform.
    FFT<double> fft;
    // Computing Fourier coefficients
    
    VectorXcd fft_d = fft .fwd(d) ;
    //fft_d = fft_d.head(m).eval();
    
    
    
    VectorXcd dconj = d.conjugate().eval();
    VectorXcd ifft_d = 1./n * fft.fwd(dconj).real().eval();

    
    VectorXcd coeffs = .5*(n*ifft_d+fft_d);
    VectorXcd res = coeffs.head(m);
    res.head(1) *= 1./n;
    res.tail(m-1) *= 2./n;
    std::cout << "computing " << std::endl;
    
    //Alternative method as described in the solution.
    VectorXd res2 = fft_d .head(m) . real () / n;
    res2.tail(m-1)*=2;
    
    assert(fabs(res.real().norm() - res2.norm()) < 0.000000001); //comparing to alternative solution
    return res.real();

    //VectorXcd coeffs = .5*(n*ifft_d+fft_d);
   // return (VectorXd) coeffs;
}


/*!
 * \brief eval_p Given polynomial coefficients, return value of polynomial
 * at $n$ equidistant points.
 *
 * \param p Coefficient vector of trigonometrix polynomial.
 * \param n Number of equidistant points at which to evaluate.
 * \return Value of polynomial $p$ at $2\pi i / n$.
 */
VectorXd eval_p(VectorXd c, unsigned int n) {

    // Degree of polynomial
    unsigned int m = c.size();

    VectorXd ret(n);
    // Loop over all points
    for (unsigned int i = 0; i < n; ++i) {
        double r = 0;
        // Loop over all coefficients
        for (unsigned int j = 0; j < m; ++j) {
            r += c(j) * std::cos(2 * M_PI * i * j / n);
        }
        ret(i) += r;
    }
    return ret;
}

int main(int argc, char **argv) {

    // Degree of trigonometric polynomial
    unsigned int m = 3;
    if(argc > 1) {
        m = std::stoi(argv[1]);
    }

    // Test points
    unsigned int npoints = 10;
    VectorXd d(npoints);
    d << 0.987214,
         1.03579,
        0.997689,
        0.917471,
         1.00474,
         0.92209,
         1.03517,
         1.08863,
        0.904992,
        0.956089;

    VectorXd g;
    // TODO: Find coefficients that best fit the data in d
    // using trig. poly of degree $m$, store the result in $g$.
    g=triPolyFit(d, m);
    std::cout << g
              << std::endl;

    // Find coordinates of best poly coeff.
    unsigned int neval = 100;
    VectorXd e = eval_p(g, neval);
    VectorXd x, y;
    x.resizeLike(e);
    y.resizeLike(e);
    for(unsigned int i = 0; i < neval; ++i) {
        x(i) = std::sin(2. * M_PI * i / neval) * e(i);
        y(i) = std::cos(2. * M_PI * i / neval) * e(i);
    }

    // Find coordinates of given data points
    VectorXd x_p, y_p;
    x_p.resizeLike(d);
    y_p.resizeLike(d);
    for(unsigned int i = 0; i < npoints; ++i) {
        x_p(i) = std::sin(2 * M_PI * i / npoints) * d(i);
        y_p(i) = std::cos(2 * M_PI * i / npoints) * d(i);
    }

    // Plot points and poly
    mglData datx, daty;
    mglData datxp, datyp;
  	datx.Link(x.data(), neval);
  	daty.Link(y.data(), neval);
  	datxp.Link(x_p.data(), npoints);
  	datyp.Link(y_p.data(), npoints);
		mglGraph gr;
		gr.Title("Orbit of the planet");
		gr.SetRanges(-1.5,1.5,-1.5,1.5);
		gr.Axis();
		gr.Plot(datx, daty, "k");
		gr.Plot(datxp, datyp, "r +");
		gr.Label('x',"X",0);
  	gr.Label('y',"Y",0);
		gr.WriteFrame("/Users/kdkdk/Desktop/ETH/numerics/numcse/numcse/exercise05/FFTLSQ/orbit.png");
		
}
