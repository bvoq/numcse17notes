//
//  main.cpp
//  numcse
//
//  Created by kdkdk on 26.09.17.
//  Copyright © 2017 a. All rights reserved.
//

#include <complex>
#include <mgl2/mgl.h>
#include <mgl2/wnd.h>
#include <mgl2/fltk.h>

#include <cassert>
#include <cmath>
#include <iostream>
using namespace std;

/*1a*/int recursiveFact(int n) {return !n?!n:n*recursiveFact(n-1);}
/*1b*/int iterfact(int n, int fact=1) {for(int i=1;i<=n;++i)fact*=i;return fact;}
/*1c - See SysProg*/
/*1d - When two functions share the same name, but not the same arguments in a namespace.*/
/*1e - Both functions return a valid result up to n = 16 with an int.*/
/*2a - Pass by reference: void p(int&a) {} - Pass by value: void p(int a) {}*/
/*2b - 1\\n2 */
/*3b - 2*/
/*3c - memory at location x[0] (the size of int)*/
/*3d - cool stuff*/
/*3e - 3*/

void blub(int a, int b) {
  cout << "yo " << a*a+b*b << endl;
}

//doesn't work without pointers.
void vfunc(void (*func)(int,int)) {
  func(3,4);
}

void simplecall(char * (*ptr)(char*,const char*)) {
  char arr [] = "damn";
  char dstarr [strlen(arr)];
  ptr(dstarr,arr);
  cout << "str: " << arr << " " << dstarr << endl;
}

void pointertests() {
  int n = 73;
  int * p = &n; // (&n) will return the adress of n.
  cout << *p << endl; //dereference operator (*p)
  //a dereference operator is an l-expression, so you can write to it.
  (*p)=3; //(*p) dereferences the pointer

  //MISCONCEPTION 1
  int *q, r;
  /* is equivalent to:
	 int *q;
	 int r;
	 */
  cout << n << endl;

  int ham [] = {5,6,7}; //c - arrays
  cout << *ham << " == " << ham[0] << endl; //should dereference ham

  //pointers to the first element of the array are called decayed pointers.
  int * decaying1 = &ham[0];
  int * decaying2 = &(*ham);
  cout << decaying1 << " == " << decaying2 << endl;

  //difference between array pointer and normal pointer.
  //re-assigning an array pointer is illegal
  //sizeof(arraypointer) = actual size of array in bytes, so sizeof(ham)/sizeof(int) returns size of array.
  //sizeof(decaying1) = size of pointer (usually 8 bytes on a 64 bit machine)
  //sizeof(*decaying1) = size of the integer which is likely 4.
  cout << sizeof(ham) << " " << sizeof(decaying1) << " " << sizeof(*decaying1) << endl;
  //ham = {2,3};

  int    pa =  3;
  int   *pb = &pa;
  int  **pc = &pb;
  int ***pd = &pc;
  assert(***pd == **pc && **pc == *pb && *pb == pa);// == *pb == pa == 3);

  int aa = 44, ba = 45, ca = 46;
  int const *ptr_a = &aa; //only disallows assignment
  int *const ptr_b = &ba; //works as expected
  //Impossible ptr_a = 32;
  ptr_a++; //possible


  char src[] = "hello";
  cout << "strlen: " << strlen(src) << " " << src << endl;
  char dst[strlen(src)];// = nullptr;
  strcpy(dst,src);
  printf("%s and %s\n", src, dst);
  cout << *dst << endl;


  char *(*strcpy_ptr)(char *dst, const char *src); // Pointer to strcpy-like function
  char *(*strcpy_ptr2)(char *, const char *); // variables don't have to be named.
  //char *strcpy(char *dst, const char *src); // An ordinary function declaration, for reference

  strcpy_ptr = &strcpy;

  char dst2[strlen(src)];
  strcpy_ptr(dst2, src);
  printf("%s and %s\n", src, dst2);

  simplecall(strcpy_ptr);
  vfunc(&blub);
  cout << "ptrs: " << " " << strcpy_ptr << endl;

  struct pointerteststruct {
	int damnboi = 4;
  };

  pointerteststruct ptstr;
  pointerteststruct * pptstr = &ptstr;
  (*pptstr).damnboi  = 4; // <=>:  pptstr->damnboi = 4;
  cout << "ELPitch: " << pptstr->damnboi << endl;

}


void heapinitstuff() {
  cout << "\n\nHEAPSTUFF\n-------\n";
  int * p = new int;
  *p = 938;
  cout << "heap stuff" << *p << endl;
  int & rp = *p; //reftyped
  rp = 2;
  cout << "you guessed it: " << *p << endl;
  delete p;
  //now p and rp are wild pointers / dangling pointers!
  cout << "ERROR WHEN ACCESSING POINTER" << endl;
  char *dp = NULL;
  {
	char c = 'o';
	dp = &c;
  }
  cout << "yo " << dp << endl;
  cout << rp << endl;



  //not in c++11char *read_only_stringptr = "Hello World"; //THIS CREATES A SUPPOSEDLY READ-ONLY POINTER.
  // some_memory[0] = 'G'; illegal
  char normal_stringptr[] = "Hello World";
  normal_stringptr[0] = 'G';
  cout << normal_stringptr << endl;
}


void scanstuff() {
  int a,b;
  scanf("%i",&a); //allows decimal, hexadecimal as 0x123 and octal using 0123.
  scanf("%d",&b); //only allows decimal
  //no difference when printing using %i vs. %d.
  printf("\n%i\n",a);
  printf("%d\n",a);
  //recommended to use %i
}

int pointerShift() {
  int *p = new int[5];
  for (int i=0; i<5; i++)
	p[i] = i;
  p++;
  return p[2];
}


class Point {
  public:
	double x, y;
	Point * next;

	double length() {
	  return sqrt(x*x+y*y)+next->length();
	}
};

void plot() {

}

int sample(mglGraph *gr)
{
    gr->Rotate(60,40);
	  gr->Box();
    return 0;
}

int main(int ,char **)
{
  pointertests();
  heapinitstuff();
  cout << "pshift " << pointerShift() << endl;
  mglFLTK gr(sample,"MathGL examples");
  return gr.Run(); /*
  gr.Alpha(true);   gr.Light(true);
  sample(&gr);              // The same drawing function.
  gr.WritePNG("test.png");  // Don't forget to save the result!
                      */
  return 0;
}
