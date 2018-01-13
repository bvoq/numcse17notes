#include <complex>
#include <mgl2/mgl.h>
#include <mgl2/wnd.h>
#include <mgl2/fltk.h>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

int main(int ,char **)
{
    Eigen::VectorXd v1 (6);
    v1 << 1,2,3,4,5,6;
    cout << v1.head(5) << endl;
    MatrixXd diagonal = v1.asDiagonal();
    cout << diagonal << endl;
    /*
  	mglGraph gr;
  	gr.Alpha(true);   gr.Light(true);
	c.plot(&gr);
    
  	gr.WritePNG("test.png");  // Don't forget to save the result!
    return 0;
     */
}
