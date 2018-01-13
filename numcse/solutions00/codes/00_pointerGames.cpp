#include <iostream>

int refAssignment() {
	int *p = new int;
	*p = 1;
	int &r = *p;
	r = 2;
	return *p;
}

int pointerShift() {
    int *p = new int[5];
    
    for (int i=0; i<5; i++) 
        p[i] = i;
    
    p++;
    
    return p[2];
}

int main() {
	std::cout << pointerShift() << std::endl;
	std::cout << refAssignment();
}
