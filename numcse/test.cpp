#include <iostream>
using namespace std;
int main() {
    int currentlyloadedin = 0;
	int secondloadedin = 0;
	cout << &currentlyloadedin << " " << &secondloadedin << endl;
    int q; scanf("%i",q);
    int * p =  (int*) q;
    int * a  = (int*) q + 0x2;
    *a = 0;
	*p = -1;

	cout << a << " " << p << endl;
	return 0;
}
