// This is a test file to test our matrix class

#include "matrix_real.hpp"
#include <iostream>
#include <complex>
#include <cmath>
using namespace std::complex_literals;



int main(int argc, char **argv) {
  int n=3;
  RMatrix<double> A(n, n, 5.0);
  RMatrix<double> B(n, n, 2.0);


  RMatrix<double> C = A + B;
  RMatrix<double> D = C;
  C(0,0) = 3.0;
  C -= A;

  RMatrix<double> E = A;
  double a = 4.0;
  E = A / a; 


  for (unsigned i=0; i<C.get_nrows(); i++) {
    for (unsigned j=0; j<C.get_ncols(); j++) {
      std::cout << E(i,j) << " ";
    }
    std::cout << std::endl;
  }

  std::cout <<std::endl;

  for (unsigned i=0; i<D.get_nrows(); i++) {
    for (unsigned j=0; j<D.get_ncols(); j++) {
      std::cout << D(i,j) << " ";
    }
    std::cout << std::endl;
  }

  return 0;
}

