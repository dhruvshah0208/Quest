#include <iostream>
using namespace std;
#include "structures.cpp"

double f(double &x,Matrix_3 &B,Vec &z){
  return (((x*x - (B.trace())*(B.trace()) + ((B + B.T()).adjoint().trace()))*(x*x - (B.trace())*(B.trace()) - z.dot(z))) - ((B + B.T()).det()+(z.dot((B + B.T())*z)))*(x - (B.trace()))  - (z.dot((B + B.T())*((B + B.T())*z))));
}

double f_bar(double &x,Matrix_3 &B,Vec &z){
  return (2*x*(2*x*x - 2*(B.trace())*(B.trace()) + ((B + B.T()).adjoint()).trace() - z.dot(z)) - z.dot((B + B.T())*z) - (B + B.T()).det());
}
double lambda_max(double guess,Matrix_3 &B,Vec &z){

  double x = guess; // sum of weights
  for (size_t i = 0; i < 5; i++) {
    x -= f(x,B,z)/f_bar(x,B,z);
  }
  return x;
}
