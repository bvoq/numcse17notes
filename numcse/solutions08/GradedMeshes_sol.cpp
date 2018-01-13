#include <iostream>
#include <Eigen/Dense>
#include <mgl2/mgl.h>

using namespace Eigen;

VectorXd evalPiecewiseInterp(const VectorXd &T, const VectorXd &Y, const VectorXd &evalT) {
	// returns the values of the piecewise linear interpolant in evalT.

	int n = evalT.size();
	VectorXd out(n);

	for (int i=0; i < n; i++) {
		for (int j=0; j < T.size()-1; j++) {
			if (evalT(i) < T(j+1) || j==T.size()-2) {
				double slope = (Y(j+1) - Y(j)) / (T(j+1) - T(j));
				out(i) = Y(j) + slope * (evalT(i) - T(j));
				break;
			}
		}
	}

	return out;
}

double maxInterpError(double a, VectorXd T, VectorXd evalT) {

	int nInterpPts = T.size();
	VectorXd Y(nInterpPts);

	for (int i=0; i<nInterpPts; i++) {
		Y(i) = std::pow(T(i), a);
	}

	VectorXd evalInterp = evalPiecewiseInterp(T, Y, evalT);

	double maxError = 0;
	for (int i=0; i<evalT.size(); i++) {
		double error = std::abs(evalInterp(i) - std::pow(evalT(i), a));
		if (error > maxError)
			maxError = error;
	}
	
	return maxError;
}


int main() {
	// Numerical experiments
	
	{// vary alpha, keep fixed n, uniform meshes
	int nInterpNodes = 8;
	int nEvalPts = 1 << 9;
	VectorXd T = VectorXd::LinSpaced(nInterpNodes, 0, 1);
	VectorXd evalT = VectorXd::LinSpaced(nEvalPts, 0, 1);

	int nParamVals = 100;
	VectorXd paramVals = VectorXd::LinSpaced(nParamVals, 0, 2); 
	VectorXd maxErrors(nParamVals);
	
	for (int i=0; i<nParamVals; i++) 
		maxErrors(i) = maxInterpError(paramVals(i), T, evalT);
	
 	mglData datx, daty;
	datx.Link(paramVals.data(), nParamVals);
	daty.Link(maxErrors.data(), nParamVals);
	mglGraph gr;
	gr.SetRanges(0, 2, 0, 2);
	gr.Axis();
	gr.Plot(datx, daty);
	gr.WriteFrame("uniformMesh_interpMaxError_varA.eps");
	}
	
	{// vary n, keep fixed alpha, uniform meshes
	int nTests = 10;
	int nEvalPts = 1 << 12;
	VectorXd evalT = VectorXd::LinSpaced(nEvalPts, 0, 1);
	VectorXd maxErrors(nTests);

	for (int n=0; n<nTests; n++) {
		int nInterpNodes = 1 << n;
		VectorXd T = VectorXd::LinSpaced(nInterpNodes, 0, 1);
		maxErrors(n) = std::log(maxInterpError(0.531, T, evalT));
	}
	
 	mglData daty;
	daty.Link(maxErrors.tail(nTests).data(), nTests);
	mglGraph gr;
	gr.SetRanges(0, nTests, -6, 0);
	gr.Axis();
	gr.Plot(daty);
	gr.WriteFrame("uniformMesh_interpMaxErrorLog_varN.eps");
	}

	double a = 0.5; // varying this manually
	double b = 2/a; // varied this manually to find best value;
					// better idea: vary it automatically, using
					// linear regression to estimate convergence rate.

	{// vary n, keep fixed alpha = 1/2, graded mesh
	int nTests = 10;
	int nEvalPts = 1 << 12;
	VectorXd evalT = VectorXd::LinSpaced(nEvalPts, 0, 1);
	VectorXd maxErrors(nTests);

	for (int n=0; n<nTests; n++) {
		int nInterpNodes = 1 << n;
		VectorXd T = VectorXd::LinSpaced(nInterpNodes, 0, 1);
		for (int i=0; i<nInterpNodes; i++) {
			T(i) = std::pow(T(i), b);
		}

		maxErrors(n) = std::log(maxInterpError(0.5, T, evalT));
	}
	
 	mglData daty;
	daty.Link(maxErrors.tail(nTests).data(), nTests);
	mglGraph gr;
	gr.SetRanges(0, nTests, -16, 0);
	gr.Axis();
	gr.Plot(daty);
	gr.WriteFrame("gradedMesh_interpMaxErrorLog_varN.eps");
	}

}
