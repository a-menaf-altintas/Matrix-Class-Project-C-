// This is a test file to test our matrix class

#include "matrix_real.h"
#include <iostream>
#include <complex>
#include <cmath>
using namespace std::complex_literals;



int main(int argc, char **argv) {
  int n=3;
  int m=4;
  RMatrix<double> A(m, n, 2);
  RMatrix<double> B(n, n, 2);
  std::cout << "A and B have been initilized" <<std::endl;


  RMatrix<double> C(m,n,0.0);
  C = A * B;
  std::cout << "A and B have been multiplied and C has been reassigned" <<std::endl;
  RMatrix<double> D = C;

  A *= B;
  C(m-1,n-1) = 3.0;


  RMatrix<double> E = A;
  double a = 4.0;
  E = A / a; 

  int length = n;
  //length++;
  //length--;
  std::vector<double> v(length, 2.0);

  std::vector<double> new_vec = A * v; // Matrix vector multiplication

  std::cout << std::endl;
  for (unsigned i=0; i<new_vec.size(); i++) {
      std::cout << new_vec[i] << " ";
    }
   std::cout << std::endl;
   std::cout << std::endl;


  for (unsigned i=0; i<C.get_nrows(); i++) {
    for (unsigned j=0; j<C.get_ncols(); j++) {
      std::cout << C(i,j) << " ";
    }
    std::cout << std::endl;
  }

  std::cout <<std::endl;

  for (unsigned i=0; i<A.get_nrows(); i++) {
    for (unsigned j=0; j<A.get_ncols(); j++) {
      std::cout << A(i,j) << " ";
    }
    std::cout << std::endl;
  }

  return 0;
}

