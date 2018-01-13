#include <cmath>
#include <iostream>
#include <iomanip>
#include <Eigen/Eigen>

using namespace Eigen;
using namespace std;

struct Level {
    double width;
    vector<double> seglengths; //cm
    double N, R, totA; //a*N+b*R;
};
Level L1() {
    Level L;
    L.width = 0.45/100.; //use m instead of cm
    L.seglengths = {1,2,3,4,5,6,7,8};
    for(int i = 0; i < L.seglengths.size(); ++i) L.seglengths[i] *= L.width;
    L.N = L.seglengths.size()-.5; //= number of corners during the trip + .5 // The plus .5 comes from start and end of acceleration.
    L.R = 0; L.totA = 0;
    for(int i = 0; i < L.seglengths.size(); ++i) L.totA += L.seglengths[i];
    for(int i = 0; i < L.seglengths.size(); ++i) {
        L.R += log2((double)L.seglengths[i]/(L.totA*L.width)+1);
    }
    return L;
}

Level L2() {
    Level L;
    L.width = 0.325/100.;
    L.seglengths = {1,2,3,4,5,6,7,8,9,10,11,12,13,13};
    for(int i = 0; i < L.seglengths.size(); ++i) L.seglengths[i] *= L.width;
    L.N = L.seglengths.size()-.5; //= number of corners during the trip + .5 // The plus .5 comes from start and end of acceleration.
    L.R = 0; L.totA = 0;
    for(int i = 0; i < L.seglengths.size(); ++i) L.totA += L.seglengths[i];
    for(int i = 0; i < L.seglengths.size(); ++i) {
        L.R += log2((double)L.seglengths[i]/(L.totA*L.width)+1);
    }
    return L;
}

Level L3() {
    Level L;
    L.width = 0.45/100.;
    L.seglengths = {9,2,9,2,9,2,9};
    for(int i = 0; i < L.seglengths.size(); ++i) L.seglengths[i] *= L.width;
    L.N = L.seglengths.size()-.5; //= number of corners during the trip + .5 // The plus .5 comes from start and end of acceleration.
    L.R = 0; L.totA = 0;
    for(int i = 0; i < L.seglengths.size(); ++i) L.totA += L.seglengths[i];
    for(int i = 0; i < L.seglengths.size(); ++i) {
        L.R += log2((double)L.seglengths[i]/(L.totA*L.width)+1);
    }
    return L;
}

Level L4() {
    Level L;
    L.width = 0.235/100.;
    L.seglengths = {19,2,19,2,19,2,19};
    for(int i = 0; i < L.seglengths.size(); ++i) L.seglengths[i] *= L.width;
    L.N = L.seglengths.size()-.5; //= number of corners during the trip + .5 // The plus .5 comes from start and end of acceleration.
    L.R = 0; L.totA = 0;
    for(int i = 0; i < L.seglengths.size(); ++i) L.totA += L.seglengths[i];
    for(int i = 0; i < L.seglengths.size(); ++i) {
        L.R += log2((double)L.seglengths[i]/(L.totA*L.width)+1);
    }
    return L;
}



pair<double,double> computeAandB(vector<vector<double> > data) {
    assert(data.size() > 0);
    assert(data[0].size() == 4);
    
    //Least squares mit allen Daten
    //for each data set compute a&b from 4 samples:
    MatrixXd M(4*data.size(),2); //first col = N for each level, second col = R = sum log2...
    VectorXd B(4*data.size()); //measured times
    for(int i = 0; i < data.size(); i++) {
        M(4*i+0,0) = L1().N;
        M(4*i+1,0) = L2().N;
        M(4*i+2,0) = L3().N;
        M(4*i+3,0) = L4().N;
    }
    
    for(int i = 0; i < data.size(); i++) {
        M(4*i+0,1) = L1().R;
        M(4*i+1,1) = L2().R;
        M(4*i+2,1) = L3().R;
        M(4*i+3,1) = L4().R;
    }
    
    for(int i = 0; i < data.size(); ++i) {
        B(4*i+0) = data[i][0];
        B(4*i+1) = data[i][1];
        B(4*i+2) = data[i][2];
        B(4*i+3) = data[i][3];
    }
    
    //make sure the matrix has a low condition number
    assert(M.jacobiSvd().rank() == 2); // M has cond number 59, which seems alright considering the values inside it are also ~ 59.
    
    VectorXd sv = M.jacobiSvd().singularValues();
    //cout << "Condition number: " << sv(0) / sv(sv.size()-1) << endl;
    
    //Solve for b using 3 methods and check whether the results are equal.
    VectorXd res = M.jacobiSvd(ComputeThinU | ComputeThinV).solve(B);
    VectorXd res2 = M.colPivHouseholderQr().solve(B); //just to be certain
    VectorXd res3 = (M.transpose().eval()*M).inverse().eval() * M.transpose().eval() * B;
    assert(fabs(res.norm() - res2.norm()) <.0000001);
    assert(fabs(res.norm() - res3.norm()) <.0000001);
    assert(fabs(res2.norm() - res3.norm()) <.0000001);
    return {(double)res(0),(double)res(1)};
}


void computeAandBForEveryOne(vector<vector<double> > data) {
    for(int i = 0; i < data.size(); ++i) {
        //for each data set compute a&b from 4 samples:
        MatrixXd M(4,2); //first col = N for each level, second col = R = sum log2...
        M(0,0) = L1().N;
        M(1,0) = L2().N;
        M(2,0) = L3().N;
        M(3,0) = L4().N;
        
        M(0,1) = L1().R;
        M(1,1) = L2().R;
        M(2,1) = L3().R;
        M(3,1) = L4().R;
        
        VectorXd B(4); //measured times
        B(0) = data[i][0];
        B(1) = data[i][1];
        B(2) = data[i][2];
        B(3) = data[i][3];
        
        //make sure the matrix has a low condition number
        assert(M.jacobiSvd().rank() == 2); // M has cond number 59, which seems alright considering the values inside it are also ~ 59.
        
        //VectorXd sv = M.jacobiSvd().singularValues();
        //cout << "Condition number: " << sv(0) / sv(sv.size()-1) << endl;
        
        //Solve for b using 3 methods and check whether the results are equal.
        VectorXd res = M.jacobiSvd(ComputeThinU | ComputeThinV).solve(B);
        VectorXd res2 = M.colPivHouseholderQr().solve(B); //just to be certain
        VectorXd res3 = (M.transpose().eval()*M).inverse().eval() * M.transpose().eval() * B;
        assert(fabs(res.norm() - res2.norm()) <.00000001);
        assert(fabs(res.norm() - res3.norm()) <.00000001);
        assert(fabs(res2.norm() - res3.norm()) <.00000001);
        
        cout << "a["<<i<<"]: " << res(0) << "\t\tb[" << i << "] " << res(1) << endl;
        //VectorXd res = pqrpT.solve(b);
    }
}

int main() {
    cout << std::setprecision(9) << endl;
    vector<vector<double> > twofingerdata =
    {   {11519856457,17433822480,8083502388,16368155692},
        {12281081029,31750333183,17665569143,32648219627},
        {21018899936,33797036841,16885252424,19000336099},
        {13965432259,32548659159,13000436128,30649841148},
        {13020136343,28466825425,13429342065,18433927301},
        {7400984317,21599795126,10865010691,12417041657},
        {13669093668,22167036945,10297827519,15884154446},
    };
    

    
    vector<vector<double> > onefingerdata = {
        {7265906735,25914588504,8982671165,18666016992},
        {13669093668,22167036945,10297827519,15884154446},
        {4718789472,13799121623,4969440105,9568281468}
    };
    
    //convert from ns to s
    for(auto & a : twofingerdata) for(auto & b : a) b /= 1000000000.; // /= 10^9
    for(auto & a : onefingerdata) for(auto & b : a) b /= 1000000000.;

    vector<vector<double> > data;
    for(const auto & a : onefingerdata) data.push_back(a);
    for(const auto & a : twofingerdata) data.push_back(a);
    
    
    //MatrixXd L1(twofingerdata.size(),2);
    cout << "L1: " << L1().R << endl;
    cout << "L2: " << L2().R << endl;
    cout << "L3: " << L3().R << endl;
    cout << "L4: " << L4().R << endl;
    
    
    
    cout << "a: " << computeAandB(onefingerdata).first  << " b: " << computeAandB(onefingerdata).second << endl;
    computeAandBForEveryOne(onefingerdata);

    
    
        //VectorXd res = pqrpT.solve(b);


    
    return 0;
}
