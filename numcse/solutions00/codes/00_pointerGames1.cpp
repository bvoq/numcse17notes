#include <iostream>

int refAssignment() {
	int *p = new int;
	int r = *p;
	r = 2;
	return *p;
}

int main() {
	std::cout << refAssignment();	
}
