// This is a test file to test our matrix class

#include "matrix_real.h"
#include <iostream>
#include <complex>
#include <cmath>
using namespace std::complex_literals;



int main(int argc, char **argv) {
  int n=3;
  int m=3;
  RMatrix<double> A(n, n, 5.0);
  RMatrix<double> B(m, m, 3.1);


  RMatrix<double> C = A + B;
  RMatrix<double> D = C;

  A += B;
  C(0,0) = 3.0;
  C -= A;

  RMatrix<double> E = A;
  double a = 4.0;
  E = A / a; 

  int length = E.get_ncols();
  //length++;
  //length--;
  std::vector<double> v(length, 2.0);

  std::vector<double> new_vec = A * v; // Matrix vector multiplication

  std::cout << std::endl;
  for (unsigned i=0; i<C.get_nrows(); i++) {
      std::cout << v[i] << " ";
    }
   std::cout << std::endl;
   std::cout << std::endl;


  for (unsigned i=0; i<C.get_nrows(); i++) {
    for (unsigned j=0; j<C.get_ncols(); j++) {
      std::cout << E(i,j) << " ";
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

