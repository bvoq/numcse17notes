#include <complex>
#include <mgl2/mgl.h>
#include <mgl2/wnd.h>
#include <mgl2/fltk.h>
#include <iostream>
#include <vector>

#include <Eigen/Dense>

using namespace std;


class Point {
private: void plotGraph(mglData & x, mglData & y, int i) {
    cout << "x: " << this->x << " y: " << this->y << endl;
		 	x.a[i] =(float)(this->x);
    y.a[i] =(float)(this->y);
    if(next != nullptr) next->plotGraph(x,y,i+1);
}
public:
    double x, y;
    Point * next;
    Point() { next = nullptr; }
    Point(double _x, double _y) : x(_x), y(_y) { next = nullptr;}
    Point(double _x, double _y, Point & _next) : x(_x), y(_y) { next = &_next;}
    ~Point() {}
    
    double length() {
        return next ? (next->length() + (next->x-x)*(next->x-x)+(next->y-y)*(next->y-y) ) : 0;
    }
    
    size_t size() {
        return 1+ (next!=nullptr ? next->size() : 0);
    }
    
    void plot(mglGraph * gr) {
        cout << "plotting"  << endl;
        mglData x(size()), y(size());
        cout << "init with size: " << size() << endl;
        plotGraph(x,y,0);
        gr->Plot(x,y);
    }
};


int sample(mglGraph *gr)
{
    Point a(0.3,0.2);
    Point b(0.4,0.1,a);
    Point c(0.8,0.9,b);
    c.plot(gr);
    
    return 0;
}


using namespace Eigen;
int main()
{
    MatrixXd m(3,3);
    m << 1,2,3,4,5,6,7,8,9;
    cout << "m =\n" << endl << m << endl;
    for(int i = 0; i < m.rows(); ++i) {
        for(int j = 0; j < m.cols(); ++j) {
            cout << m(i,j) << ",";
        }
    }cout << endl;
    
    for(int i = 0; i < m.size(); ++i) {
        cout << m(i) << ",";
    }cout << endl;
    
    //FAST
    for(int i = 0; i < m.cols(); ++i) {
        cout << i << "-th col: " << endl << m.col(i) << endl << "-----" << endl;
    }
    
    //SLOW, will give matrix 1xn
    for(int i = 0; i < m.rows(); ++i) {
        cout << i << "-th row: " << endl << m.row(i) << endl << "-----" << endl;
    }
    
    //EVEN SLOWER, will give vector
    for(int i = 0; i < m.rows(); ++i) {
        cout << i << "-th row: " << endl << (VectorXd) m.row(i) << endl << "-----" << endl;
    }
    VectorXd veccol = m.col(2);
    VectorXd vecrow = m.row(2);
    
    cout << "OK STORED COL VECTOR:" << endl << veccol << endl;
    cout << "OK STORED ROW VECTOR:" << endl << vecrow << endl;
    
    //m << {1,2,3},{2,3,4},{3,4,5};
    
    
    mglFLTK gr(sample,"MathGL examples");
    return gr.Run();
    /* ALTERNATIVE WITHOUT USING FLTK
     mglGraph gr;
     gr.Alpha(true);   gr.Light(true);
     c.plot(&gr);
     
     gr.WritePNG("test.png");  // Don't forget to save the result!
     return 0;
     */
}
