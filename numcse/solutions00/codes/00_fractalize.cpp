// to compile run: g++ -std=gnu++11 01_fractalize.cpp -lmgl

#include <iostream>
#include <cmath>
#include <string>
#include <mgl2/mgl.h>
#include <mgl2/wnd.h>
#include <mgl2/fltk.h>


using namespace std;

class Point {
    public:
    double x,y;
    Point *next;    
    
    Point (double x, double y) {
        this->x = x;
        this->y = y;
        this->next = NULL;
    }
    
    static double dist(Point &a, Point &b) {
		// returns euclidean distance between two points
		double dx = a.x - b.x;
		double dy = a.y - b.y;
        return sqrt(dx*dx + dy*dy);
    }
    
    double length() {
		// returns euclidean length of polygonal curve 
		Point *iter = this;
		double len = 0;
        while (iter->next != NULL) {
			len += dist(*iter, *iter->next);
        	iter=iter->next;
		}
		return len;
    }
    
    void fractalizeSegment() {
		// constructs a triangular bump on the first segment

		Point *n = this->next;
		// calculation of coordinates of the points used in the construction
        double xonethird = (2*this->x + n->x)/3;
        double yonethird = (2*this->y + n->y)/3;
        double xmid = 1./2*(this->x + n->x);
        double ymid = 1./2*(this->y + n->y);
        double xtwothird = (this->x + 2*n->x)/3;
        double ytwothird = (this->y + 2*n->y)/3;
		double yperp = xtwothird - xonethird;
		double xperp = - ytwothird + yonethird;
        
		// creation of points and pointer referencing
		Point *p1 = new Point(xonethird, yonethird);
        Point *p2 = new Point(xmid + xperp*sqrt(3)/2, ymid + yperp*sqrt(3)/2);
        Point *p3 = new Point(xtwothird, ytwothird);
        p3->next = n;
        this->next = p1;
        p1->next = p2;
        p2->next = p3;
    }    
        
	void fractalize() {
		// repeats fractalizeSegment on every subsegment of a iteratively
		Point *iter1 = this;
		Point *iter2 =this->next;
		while (iter2 != NULL) {
			iter1->fractalizeSegment();
			iter1 = iter2;
			iter2 = iter1->next;
		}
	}

	string toString() {
		// returns a description of the curve as sequence of coordinates
		string s = string();
		Point *iter = this;
		s += "[("+to_string(iter->x)+","+to_string(iter->y)+")";
		while (iter->next != NULL) {
			iter = iter->next;
			s += ", ("+to_string(iter->x)+","+to_string(iter->y)+")";
		}
		s += "]";
		return s;
	}

	void toArrays(double* xcoords, double* ycoords, int len) {
		// writes in *xcoords and *ycoords the coordinates of
		// the first len points of the curve
		Point *iter = this;
		for (int i=0; i<len; i++) {
			xcoords[i] = iter->x;
			ycoords[i] = iter->y;
			iter = iter->next;
		}
	}

	int countPoints() {
		// returns number of points until end of curve
		// (current point included in the count)
		int counter = 1;		
		Point *iter = this->next;
        while (iter != NULL) {
			counter++;
        	iter=iter->next;
		}
		return counter;
	}

	void plot(const char *name) {
		// saves plot of the curve in the file *name
		int len = this->countPoints();
		double* xcoords = new double[len];
		double* ycoords = new double[len];
		this->toArrays(xcoords, ycoords, len);

  		mglData datx, daty;
  		datx.Link(xcoords, len);
  		daty.Link(ycoords, len);
		mglGraph gr;
		gr.SetRanges(0,1,0,1);
		gr.Plot(datx, daty, "0");
		mglPoint pleftdown(0,0), prightup(1,1);
		gr.SetCutBox(pleftdown, prightup);
		gr.WriteFrame(name);
	}
};

int nStages = 5;
	double *lengths = new double[nStages];
int sample(mglGraph * gr) {
    mglData data;
    data.Link(lengths, nStages);
    gr->SetRanges(0,nStages,0,nStages*3);
    gr->Plot(data);
    gr->Axis();
    gr->Label('x',"Stage");
    gr->Label('y',"Length");
    gr->WriteFrame("lengthsPlot.png");
    return 0;
}

int main() {
	// Initialization
    Point p1(0,0);
    Point p2(1,0);
    p1.next = &p2;

	// fractalizes the segment and calculates the lengths
	for (int k=0; k<nStages; k++) {
		p1.fractalize();
		lengths[k] = p1.length();
	}

	p1.plot("sample1.eps");
	// plots length of curve as a function of stage
    mglFLTK gr(sample,"MathGL examples");
    return gr.Run();
}
