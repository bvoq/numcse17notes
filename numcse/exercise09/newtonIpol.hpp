# include <Eigen/Dense>

using Eigen::VectorXd;

// IN: t = node set (mutually different)
// y = nodal values
// OUT: coeffs = coefficients of polynomial in Newton basis
VectorXd divDiff(const VectorXd& t, const VectorXd& y) {
    VectorXd a = y;
    int n = t.size()-1;
    for(unsigned l=0; l<n; ++l) for(unsigned j=l; j<n; ++j) a(j+1) = (a(j+1)-a(l))/(t(j+1)-t(l));
    return a;
}


// Evaluation of polynomial in Newton basis (divided differences)
// IN: t = nodes (mutually different)
// y = values in t
// x = evaluation points (as Eigen::Vector)
// OUT: p = values in x */
void intPolyEval (const VectorXd& t, const VectorXd& y, const VectorXd& x, VectorXd& p) {

	const unsigned n = y.size()-1 ;
	// get Newton coefficients of polynomial (non in-situ implementation!)
	VectorXd coeffs; 
	coeffs = divDiff(t, y);

	// evaluate
	VectorXd ones = VectorXd::Ones(x.size());
	p = coeffs(n)*ones;
	for (int j=n-1; j>= 0; --j)
		p = (x - t(j)*ones).cwiseProduct(p) + coeffs(j)*ones;

}


