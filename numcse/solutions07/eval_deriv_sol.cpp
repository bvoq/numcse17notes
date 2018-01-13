#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

/*!
 * \brief dipoleval
 * \param t
 * \param y
 * \param x
 * \return
 */
double derivIpolEvalAN(const VectorXd & t,
                 const VectorXd & y,
                 const double x) {
    
	assert(t.size() == y.size());
		
	VectorXd p(y);
	VectorXd dP = VectorXd::Zero(y.size());

	for(int i = 1; i < y.size(); ++i) {
		for(int k = i-1; k >= 0; --k) {

			dP(k) = (p(k+1) + (x-t(k))*dP(k+1) - p(k) - (x-t(i))*dP(k))
     					/ (t(i) - t(k));
                  
			p(k) = ((x-t(k))*p(k+1) - (x-t(i))*p(k)) / (t(i) - t(k));
		}
	}

	return dP(0);
}


int main() {

	VectorXd t(3), y(3);
	t << -1, 0, 1;
	y << 2, -4, 6;
    
	double x = 0.5;
	std::cout<< "Polynomial derivative @ x = " << x << " is: "<< derivIpolEvalAN(t,y,x) << std::endl;

}
