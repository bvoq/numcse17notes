#include <iostream> 
#include <fstream>
#include <sstream>
#include <string>
#include <Eigen/Eigen>

using namespace Eigen;
using namespace std;

VectorXd load_pgm(const std::string &filename) {
	// returns a picture as a flattened vector

	int row = 0, col = 0, rows = 0, cols = 0;

	std::ifstream infile(filename);
	std::stringstream ss;
	std::string inputLine = "";

	// First line : version
	std::getline(infile,inputLine);

	// Second line : comment
	std::getline(infile,inputLine);

	// Continue with a stringstream
	ss << infile.rdbuf();
	// Third line : size
	ss >> cols >> rows;

	VectorXd picture(rows*cols);

	// Following lines : data
	for(row = 0; row < rows; ++row) {
		for (col = 0; col < cols; ++col) {
			int val;
			ss >> val;
			picture(col*rows + row) = val;
		}
	}

	infile.close();
	return picture;
}

int main() {
	
	int h = 231;
	int w = 195;
	int M = 15;

	MatrixXd faces(h*w, M);
	VectorXd meanFace(h*w);
	
	// loads pictures as flattened vectors into faces
	for (int i=0; i<M; i++) {
		std::string filename = "/Users/kdkdk/Desktop/ETH/numerics/numcse/numcse/exercise04/basePictures/subject"+
			std::to_string(i+1) + ".pgm";
		VectorXd flatPic = load_pgm(filename);
		faces.col(i) = flatPic;
	}
    cout << "BLA " << endl;

    //1. Find the mean face
    for(int i = 0; i < faces.rows(); ++i) {
        double mean = 0;
        for(int j = 0; j < faces.cols(); ++j) {
            mean += faces(i,j);
        }
        mean /= (double)faces.cols();
        meanFace(i) = mean;
    }
    cout << "done meaning " << endl;

    //2. Add vectors pointing away from the mean face.
    MatrixXd A(h*w,M);
    for(int i = 0; i < faces.cols(); ++i) {
        A.col(i) = faces.col(i) - meanFace;
        cout << "err " << i << endl;
    }
    
    // cannot be computed, to large MatrixXd Acov = A * A.transpose().eval();
    
    //Ex 4b-d)
    //A*A^T has size h*w,h*w and can at most have min(M,h*w) non-zero eigenvalues, since rank(A*A^T) = rank(A) = # of non-zero eigenvalues of A / A*A^T
    
    //Eigenvectors of A*A^T = U[1:m]  A^T*A = V[1:n]
    auto svd = A.jacobiSvd(ComputeThinU | ComputeThinV);
    MatrixXd U = svd.matrixU(); //left singular values = eigenvectors of A*A^T
    
	string testPicName = "/Users/kdkdk/Desktop/ETH/numerics/numcse/numcse/exercise04/testPictures/subject03.sleepy.pgm";
	VectorXd newFace = load_pgm(testPicName);

    
    //4e)
    // Given a new face y, implement C++ code which computes its projection on the space spanned by the eigenfaces, that is finds x such that Ux = (y− mean face), where U is the matrix of singular vectors of A.
    assert(newFace.size() == meanFace.size());
    VectorXd newFaceMeaned = newFace - meanFace;
    //U*x=b=(newFace-meanFace)
    VectorXd x = U.transpose().eval()*newFaceMeaned;

    
    //4f) Implement C++ code which computes the distance between the projection of the new face and the
    // projection of a column k for a generic k, and print the k which minimizes the distance.
    //assert(x.size() == (U.transpose().eval()*A.col(0)).size());
    double mindist = (x-U.transpose().eval()*A.col(0)).norm();
    int minK = 0;
    for(int i = 0; i < A.cols(); ++i) {
        double dist = (x-U.transpose().eval()*A.col(i)).norm();
        if(dist < mindist) {
            mindist = dist;
            minK = i;
        }
    }
    cout << "k: " << minK+1 << endl;
    
    
}
