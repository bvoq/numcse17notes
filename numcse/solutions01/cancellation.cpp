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
#include <mgl2/mgl.h>
#include <mgl2/wnd.h>
#include <mgl2/fltk.h>

using namespace Eigen;

int main() {
    // Array of values of $h$
    ArrayXd h = ArrayXd::LinSpaced(21, -20, 0.)
        .unaryExpr([] (double i) {
            return std::pow(10., i);
        });
    // Dummy array where to evaluate the derivative (1.2)
    ArrayXd x = ArrayXd::Constant(h.size(), 1.2);

    // Derivative
    ArrayXd g1 = (sin(x +h) - sin(x)) / h; // naive
    ArrayXd g2 = 2 * cos(x+0.5*h) * sin(0.5 * h) / h; // better
    ArrayXd ex = cos(x); // exact

    //// Print error

    // Table header
    std::cout << std::setw(15) << "h"
              << std::setw(15) << "exact"
              << std::setw(15) << "cancellation"
              << std::setw(15) << "error"
              << std::setw(15) << "improved"
              << std::setw(15) << "error" << std::endl;
    for(unsigned int i = 0; i < h.size(); ++i) {
        // Table entries
        std::cout << std::setw(15) << h(i)
                  << std::setw(15) << ex(i)
                  << std::setw(15) << g1(i)
                  << std::setw(15) << std::abs(g1(i) - ex(i))
                  << std::setw(15) << g2(i)
                  << std::setw(15) << std::abs(g2(i) - ex(i)) << std::endl;
    }

    // Plotting with MathGL
    double ref[h.size()];
    for (int i=0; i<h.size(); i++)
    	ref[i] = h[i];
    	
    VectorXd error1 = (g1-ex).abs().matrix();
    VectorXd error2 = (g2-ex).abs().matrix();
    
    mglData stepSize, dataRef;
    stepSize.Link(h.matrix().data(),h.size());
    dataRef.Link(ref,h.size());
    
    mglData data1, data2;
  	data1.Link(error1.data(), h.size());
  	data2.Link(error2.data(), h.size());
  	
  	 mglFLTK *gr = new  mglFLTK;
    gr->Title("Error of approximation of f'(x_0)");
  	gr->SetRanges(h(0),1,1e-20,1e+0);
  	gr->Label('x',"h",0);
  	gr->Label('y', "Error",0);
  	gr->SetFunc("lg(x)","lg(y)");
  	
  	double xTicks[] = {1e-20,1e-15,1e-10,1e-5,1e+0};
  	double yTicks[] = {1e-15,1e-10,1e-5,1e+0};
  	gr->SetTicksVal('x', mglData(5,xTicks), "10^{-20}\n10^{-15}\n10^{-10}\n\\10^{-5}\n\\10^{0}");
  	gr->SetTicksVal('y', mglData(4,yTicks), "10^{-15}\n10^{-10}\n\\10^{-5}\n\\10^{0}");
		gr->Axis();
		
  	gr->Plot(stepSize,data1,"k+"); gr->AddLegend("g_1","k+");
  	gr->Plot(stepSize,data2,"r +"); gr->AddLegend("g_2","r +");
  	gr->Plot(stepSize,dataRef,"b|"); gr->AddLegend("O(h)","b|");
    gr->Legend(1);
  	gr->WriteFrame("error_cancellation.eps");
    gr->Run();
}
