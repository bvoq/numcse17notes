#include <Eigen/Eigen>
#include <cassert>
#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <iomanip>

using namespace Eigen;
using namespace std;

VectorXd vandermondeCoeff(const VectorXd & t, const VectorXd & y) {
	int n = y.size()-1;
	MatrixXd m(n+1,n+1);
	for(int i = 0; i <= n; ++i) {
		m(i,0) = 1;
		double value = t(i);
		for(int j = 1; j <= n; ++j) {
			m(i,j) = value;
			value *= t(i);
		}
	}
    cout << m << endl;

	assert(m.fullPivLu().rank() == n+1); //solvable
	//VectorXd a = m.fullPivLu().solve(y);
    VectorXd atemp = m.jacobiSvd(ComputeThinU | ComputeThinV).solve(y);
	return atemp;
	//return a;
}



VectorXd lagrangeCoeff(const VectorXd & t, const VectorXd & y, double teval) {
    int n = t.size()-1;
    VectorXd coeff(n+1);
    for(int i = 0; i <= n; ++i) {
        double res = 1;
        for(int j = 0; j <= n; ++j) {
            if(i != j) {
                res *= (teval - t(j)) / (t(i) - t(j));
            }
        }
        coeff(i) = res;
    }
    return coeff;
}

double lagrangeEval(const VectorXd & t, const VectorXd & y, double teval) {
    VectorXd L = lagrangeCoeff(t,y,teval);
    cout << "L(" << teval << ") = "; for(int i = 0; i < L.size(); ++i) cout << L(i) <<" "; cout << endl;

    int n = t.size()-1;
    VectorXd out (n+1);
    double res = 0;
    for(int i = 0; i < y.size(); ++i) {
        res += y(i)*L(i);
    }
    return res;
}

VectorXd lambdas(const VectorXd & t, const VectorXd & y) {
    int n = t.size() - 1; VectorXd out = VectorXd::Ones(n+1);
    for(int i = 0; i <= n; ++i) for(int j = 0; j <= n; ++j) if(i != j) out(i) *= 1./(t(i) - t(j));
    return out;
}

VectorXd evalBarycentric(const VectorXd & t, const VectorXd & y, const VectorXd & evalPoints) {
    VectorXd lam = lambdas(t, y);
    int n = t.size() - 1;
    VectorXd res(evalPoints.size());
    for(int tit = 0; tit < evalPoints.size(); ++tit) {
        double teval = evalPoints(tit);
        double top = 0, bot = 0;
        for(int i = 0; i <= n; ++i) top += y(i) * lam(i) / (teval - t(i));
        for(int i = 0; i <= n; ++i) bot += lam(i) / (teval-t(i));
        res(tit) = top/bot;
    }
    return res;
}

VectorXd newtonCoeffDivDiff(const VectorXd & t, const VectorXd& y) {
    VectorXd a = y;
    int n = t.size()-1;
    for(unsigned l=0; l<n; ++l) for(unsigned j=l; j<n; ++j) a(j+1) = (a(j+1)-a(l))/(t(j+1)-t(l));
    return a;
}

/*
VectorXd dividedDiffCoeffExt(const VectorXd & aold, double t, double y) {
    int n = (prevComputedCoeff.size()+1)-1;
    VectorXd anew (n); aold << prevComputedCoeff, y;
    int l=n-1;
    for(unsigned j=l;j<n;++j) {}
    for(unsigned l=0; l<n; ++l) for(unsigned j=l; j<n; ++j) outy(j+1) = (outy(j+1)-outy(l))/(t(j+1)-t(l));
    return outy;
}*/

map<pair<int,int>, double> dp;
static double divcoeff(int a, int b, const VectorXd & t, const VectorXd & y) {
    if(dp.count({a,b}) == 0) {
        if(a == b) return y[a];
        else {
            double res = (divcoeff(a+1,b, t,y) - divcoeff(a,b-1, t,y)) / (t(b) - t(a));
            dp[{a,b}] = res;
            return res;
        }
    } else return dp[{a,b}];
}
VectorXd dividedDiffCoeffStart(const VectorXd & t, const VectorXd & y) {
    dp.clear(); //assume new values
    VectorXd a (t.size());
    for(int i = 0; i < t.size(); ++i) a(i) = divcoeff(0,i, t,y);
    return a;
}

VectorXd dividedDiffExt(const VectorXd & oldt, const VectorXd & oldy, const VectorXd & newt, const VectorXd & newy) {
    //fuse the two vectors
    VectorXd t(oldt.size() + newt.size()); t << oldt, newt;
    VectorXd y(oldy.size() + newy.size()); y << oldy, newy;
    
    VectorXd a(t.size());
    for(int i = 0; i < t.size(); ++i) a(i) = divcoeff(0,i, t,y);
    return a;
}
VectorXd newtonCoeff(const VectorXd & t, const VectorXd & y) {
    int n = t.size() - 1;
    MatrixXd M = MatrixXd::Zero(n+1,n+1);
    
    /*for(int j = 0; j < n+1; ++j) {
        for(int i = j; i < n+1; ++i) {
            //Slow way of computing Ni
            double Ni = 1;
            for(int k = 0; k < j; ++k) {
                Ni *= t(i) - t(k);
            }
            
            M(i,j) = Ni;
        }
    }*/
    
    for(int i = 0; i < n+1; ++i) {
        double Ni = 1;
        for(int j = 0; j <= i; ++j) {
            M(i,j) = Ni;
            Ni *= t(i) - t(j);
        }
    }
    VectorXd acoeff = (M.triangularView<Lower>()).solve(y);
    return acoeff;
}


double newtonEval(const VectorXd & acoeff, const VectorXd & t, const VectorXd & y, double teval) {
    int n = t.size() - 1;
    double res = 0;
    double Ni = 1;
    for(int i = 0; i < n+1; ++i) {
        res += acoeff(i) * Ni;
        Ni *= (teval - t(i));
    }
    return res;
}

double newtonEval(const VectorXd & t, const VectorXd & y, double teval) {
    VectorXd acoeff = newtonCoeff(t,y);
    return newtonEval(acoeff, t, y, teval);
}


map<pair<int,int>, double> dpAN;
double aitkenNeville(int a, int b, const VectorXd & t, const VectorXd & y, double teval) {
    if(a == b) return y(a);
    else if(dpAN.count({a,b}) == 0) {
        double res = (aitkenNeville(a+1,b, t,y,teval) * (teval - t(a)) - aitkenNeville(a,b-1, t,y,teval) * (teval - t(b))) / (t(b) - t(a));
        dpAN[{a,b}] = res;
        return res;
    } else return dpAN[{a,b}];
}

double aitkenNevilleStart(const VectorXd & t, const VectorXd & y, double teval) {
    dpAN.clear();
    return aitkenNeville(0, t.size()-1, t,y,teval);
}

double aitkenNevilleExt(const VectorXd & oldt, const VectorXd & oldy, const VectorXd & newt, const VectorXd & newy, double teval) {
    VectorXd t (oldt.size() + newt.size()); t << oldt, newt;
    VectorXd y (oldy.size() + newy.size()); y << oldy, newy;
    return aitkenNeville(0, t.size()-1, t,y,teval);
}

int main() {
	VectorXd a(4); a << 1,2,4,3;
	VectorXd b(4); b << -6,2,12,-10;
	VectorXd c = vandermondeCoeff(a, b);
    for(int i = 0; i < c.size(); ++i) cout << "+"<<c(i) << "x^"<<i<<" "; cout << endl;

    double d = lagrangeEval(a, b, 1.5);
    cout << "lagrangeeval: " << d << endl;
    
    VectorXd e = lambdas(a, b);
    for(int i = 0; i < e.size(); ++i) cout << e(i) <<" "; cout << endl;

    VectorXd multipleEval(2); multipleEval << 2.5, 1.5;
    VectorXd f = evalBarycentric(a, b, multipleEval);
    for(int i = 0; i < f.size(); ++i) cout << f(i) <<" "; cout << endl;

    VectorXd g = newtonCoeff(a, b);
    for(int i = 0; i < g.size(); ++i) cout << g(i) <<" "; cout << endl;
    
    double gd = newtonEval(a,b,2.5);
    cout << "newtoneval: " <<  gd << endl;
    
    VectorXd h = newtonCoeffDivDiff(a,b);
    cout << "ethdivdiff: "; for(int i = 0; i < h.size(); ++i) cout << h(i) <<" "; cout << endl;

    VectorXd k = dividedDiffCoeffStart(a, b);
    cout << "div-diff rec: "; for(int i = 0; i < k.size(); ++i) cout << k(i) <<" "; cout << endl;

    VectorXd aext(2); aext <<-1,-2;
    VectorXd bext(2); bext <<15,2;
    VectorXd atot(a.size() + aext.size()); atot << a, aext;
    VectorXd btot(b.size() + bext.size()); btot << b, bext;
    VectorXd l = dividedDiffExt(a,b,aext,bext);
    cout << "div-diff-ext: "; for(int i = 0; i < l.size(); ++i) cout << l(i) <<" "; cout << endl;
    cout << "div-diff-ext eval(1.5): " << newtonEval(l,atot,btot,2.5) << endl;
    VectorXd lcheck = newtonCoeff(atot, btot);
    cout << "div-diff-ext double-check: "; for(int i = 0; i < lcheck.size(); ++i) cout << lcheck(i) <<" "; cout << endl;
    cout << "div-diff-ext double-check eval(1.5): "; cout << newtonEval(lcheck,atot,btot,2.5) << endl;
    
    cout << "aitken neville eval 1.5: " << aitkenNevilleStart(a,b,1.5) << endl;
    cout << "aitken neville eval 2.5: " << aitkenNevilleStart(a,b,2.5) << endl;
    cout << "aitken neville eval 2.5 ext: " << aitkenNevilleExt(a,b,aext,bext,2.5) << endl;
    
    
    
    //For gauss integration, check weights
    
    //Midpoint rule:
    VectorXd azero(1); azero << 2;
    VectorXd bzero(1); bzero << 1;
    VectorXd czero = vandermondeCoeff(azero, bzero);
    std::cout << std::setprecision(20);
    for(int i = 0; i < czero.size(); ++i) cout << "+"<<czero(i) << "x^"<<i<<" "; cout << endl;
    
    return 0;
}
