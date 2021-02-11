#include<stdio.h>
typedef struct Vec Vec;
struct Vec{
  double x,y,z;
};
void Vec_construct(Vec *this, double x,double y,double z)
{
  this->x = x;
  this->y = y;
  this->z = z;
}
int main(){
  Vec a;
  int x = 2;
  Vec_construct(&a,2,3,4);
  int test(int y){
    return x;
  }
  printf("Hello World");
}

/*
double f(double x,Matrix_3 *B_ptr,Matrix_3 *S_ptr,Vec *z_ptr){
  return (((x*x - (B.trace())*(B.trace()) + ((B + B.T()).adjoint().trace()))*(x*x - (B.trace())*(B.trace()) - z.dot(z))) -
          ((B + B.T()).det()+(z.dot((B + B.T())*z)))*(x - (B.trace()))  - (z.dot((B + B.T())*((B + B.T())*z))));
  double trace_B;
  Matrix_3 adj;
  adjoiont(&adj,S_ptr);
  double K = trace(&adj);
  trace_B = trace(&B_ptr);

  double first = (x*x - trace_B*trace_B + K)*(x*x - trace_B*trace_B - ((z_ptr->x)*(z_ptr->x) + (z_ptr->y)*(z_ptr->y) + (z_ptr->z)*(z_ptr->z)));
  double second = (x - trace_B)*();
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
*/
