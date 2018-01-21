#include <iostream>
#include <cmath>
#include <cassert>
using namespace std;

double func1 (double x) {
    return 10./ (x*x*x-1);
}
double func2 (double x) {
    return x*x*x*x - 10;
}

double func3 (double x) {
    return pow(x+10,0.25);
}

template <class Function> void printevalsquared(const Function& f, double x) {
    cout << f(x*x)*f(x*x) << endl;
}

template <class Function> double bisect(double a, double b, const Function & f, int depth) {
    if(depth > 10) return a;
    assert(f(a) * f(b) < 0);
    double mid = (a+b)/2.;
    cout << depth+1 << ": " << "["<<a<<","<<b<<"]  size: "<<b-a<<" f(" << mid << ")="<<f(mid) << "diff: " << endl;

    if(f(mid) == 0) return mid;
    if(f(a)*f(mid) < 0) return bisect(a,mid,f, depth+1);
    else return bisect(mid,b,f,depth+1);
}

template <class Function> double fixedPoint(double e, const Function & g, int depth) {
    
    cout << depth << ": " << e << endl;
    if(depth > 20) return e;
    return fixedPoint(g(e), g, depth+1);
}

template <class Function> double rootFindWithNewtonRaphson(double e, const Function & f, const Function & fderiv) {
    auto g = [f,fderiv] (double x) { return x - f(x) / fderiv(x); };
    return fixedPoint(e,g,0);
}



#include <iomanip>
int main() {
    double r = -3.5;
    std::function<double(double)> f = [] (double x) { return 3*x+sin(x)-exp(x);};
    std::function<double(double)> g = [] (double x) { return x*x*x-5; };
    std::function<double(double)> h = [] (double x) { return (x-2)*(x-2)*(x-2)+(x-2)*(x-2)-1; };
    std::function<double(double)> i = [] (double x) { return pow(2,x)+pow(2,-x)-3; };

    std::cout << std::setprecision(25);
    
    std::function<double(double)> n1 = [] (double x) { return 4*x*x-8; };
    std::function<double(double)> n1deriv = [] (double x) { return 8*x; };

    
    std::function<double(double)> fix = [] (double x) { return exp(-x); };
    fixedPoint(100000000000,fix,0);
    //rootFindWithNewtonRaphson(12, n1, n1deriv);
    
    //double res = bisect(1,2,i,-1);
    
}
