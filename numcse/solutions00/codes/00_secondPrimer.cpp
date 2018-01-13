#include <iostream>
#include <chrono>

using namespace std;

template<typename T>
T recFact(T n) {
	if (n<=0) return 1;
	return n * recFact(n-1);
}

template<typename T>
T iterFact(T n) {
	T out = 1;
	while (n>0) {
		out *= n;
		n--;
	}
	return out;
}

template<typename T>
double funcTimer (T (f)(T), int* input, int inputL) {
	// returns the time the function f took to run on
	// all elements of input.
	auto start = std::chrono::system_clock::now();
	for (int i=0; i<inputL; i++) {
		f(input[i]);
	}
	auto end = std::chrono::system_clock::now();
	auto elapsed = end - start;
	return elapsed.count();
}

int* createRandArray(int N, int k1, int k2) {
	// returns N integers chosen randomly from [k1,k2]
	int* randNums = new int[N];
	srand(time(NULL));
	for (int i=0; i<N; i++) {
		randNums[i] = k1 + rand() % (1 + k2 - k1);
	}
	return randNums;
}

int main ()
{
	int N = 10000000, k1 = 5, k2 = 25;
	int *randNums = createRandArray(N, k1, k2);
	double timeRecFact = funcTimer(recFact<long>, randNums, N);
	double timeIterFact = funcTimer(iterFact<long>, randNums, N);
	cout << timeRecFact / timeIterFact;
}
