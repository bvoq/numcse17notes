#include <iostream>

void copyAddOne(int n) {
	n++;
}

void refAddOne(int &n) {
	n++;
}

int main ()
{
	int n = 1;

	copyAddOne(n);
	std::cout << n << std::endl;

	refAddOne(n);
	std::cout << n;

}
