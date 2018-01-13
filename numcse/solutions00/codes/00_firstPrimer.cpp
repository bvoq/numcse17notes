#include <iostream>
#include <climits>
#include <float.h>

using namespace std;

int recurFact(int n) {
	// recursive implementation of the factorial
	if (n==0) {
		return 1;
	}
	return n * recurFact(n-1);
}

int iterFact(int n) {
	// iterative implementation of the factorial
	int out = 1;
	while (n>0) {
		out *= n;
		n--;
	}
	return out;
}


long iterFact(long n) {
	// iterative implementation of the factorial
	long out = 1;
	while (n>0) {
		out *= n;
		n--;
	}
	return out;
}

double iterFact(double n) {
	// iterative implementation of the factorial
	double out = 1;
	while (n>0) {
		out *= n;
		n--;
	}
	return out;
}

int main() {
	// To find the maximum admissible value for the factorial, consider the
	// maximum value of the type considered and divide it for
	// 2, 3, 4, ... until you get something smaller than 1. 
	int i = 1;
	double x = INT_MAX;
	while (x >= 1) {
		x /= ++i;
	}	
	cout << "Factorial of " << i-1 << " is " << iterFact(i-1)
		 << ", but factorial of " << i << " is not "
		 << iterFact(i) << endl;

		
	x = LONG_MAX, i = 1;
	while (x >= 1) {
		x /= ++i;
	}	
	cout << "Factorial of " << i-1 << " is " << iterFact((long) i-1)
		 << ", but factorial of " << i << " is not "
		 << iterFact((long) i) << endl;
	
	x = DBL_MAX, i = 1;
	while (x >= 1) {
		x /= ++i;
	}	
	cout << "Factorial of " << i-1 << " is " << iterFact((double) i-1)
		 << ", but factorial of " << i << " is not "
		 << iterFact((double) i) << endl;
}
