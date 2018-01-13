/******************************************************************************

Welcome to GDB Online.
GDB online is an online compiler and debugger tool for C, C++, Python.
Code, Compile, Run and Debug online from anywhere in world.

*******************************************************************************/
#include <cstdlib>
#include <ctime>
#include <vector>
#include <iostream>
#include <math.h>       /* sqrt */
class Point{
    Point * _next = nullptr;
    double _x, _y;
public:
    Point () {} // default Constructor
    ~Point () {} // default Destructor
    
    Point(double x, double y) :
    _x(x), _y(y) {}
    
    void set(double x, double y){
        _x = x;
        _y = y;
    }
    void setNext(Point * p){
        _next = p;
    }
    double getX(){
        return _x;
    }
    double getY(){
        return _y;
    }
    Point * getNext(){
        return _next;
    }
    double curveLength(){
        if(_next){
            return sqrt( pow(_x-_next->getX(),2) + pow(_y-_next->getY(),2) ) + _next->curveLength();
            
        }
        return 0;
    }
    void print(){
        std::cout << _x << " " << _y << std::endl;
    }
    void plot(){
        this->print();
        if(_next){
            _next->plot();
        }
    }
    void fractalizeSegment(){
        if(_next){
            Point *p1 = new Point((_x+_next->getX())/3,(_y+_next->getY())/3);
            Point *p2 = new Point();
            Point *p3 = new Point(2*(_x+_next->getX())/3,2*(_y+_next->getY())/3);
            double len = 1/4 * (sqrt( pow(_x-_next->getX(),2) + pow(_y-_next->getY(),2))); //height of the triangle
            double x2,y2;
            if( _x <= _next->getX()){       //we are lower, go to the left
                x2 = (_x+_next->getX())/2 - len * ((_y+_next->getY()));
                y2 = (_y+_next->getY())/2 + len * ((_x+_next->getX()));
            }else{
                x2 = (_x+_next->getX())/2 + len * ((_y+_next->getY()));
                y2 = (_y+_next->getY())/2 - len * ((_x+_next->getX()));
            }
            p2->set(x2,y2);
            
            p1->setNext(p2);
            p2->setNext(p3);
            p3->setNext(_next);
            _next = p1;
        }
    }
};

int main(int argc, char const* argv[]) {
	const int n = 2;
	std::vector<Point> p(n);
	Point p1(0,0);
	Point p2(1,1);
	p1.setNext(&p2);
	p1.fractalizeSegment();
	p1.plot();
	/*for (int i = 0; i < n; i++) {
		p[i].print();
	}*/
	// std::cout << p[0].curveLength() << std::endl;
	return 0;
}
