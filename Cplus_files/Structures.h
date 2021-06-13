#include <iostream>
using namespace std;

class Vec{
public:
  double x,y,z;
  Vec(double p=0,double q = 0,double r = 0);
  Vec cross(Vec &b);
  double dot(Vec b);
  Vec operator* (double &t);
  Vec operator+ (Vec b);
  Vec operator- (Vec b);
//  Matrix_3 outer_product(Vec b); // Define this as a global function
};

class Matrix_3{
public:
  double elements[3][3];
  Matrix_3(){};
  Matrix_3(int x);
  double &operator ()(int i,int j);
  void assign(int i,int j,double z);
  Matrix_3 operator+ (Matrix_3 b);
  Matrix_3 operator- (Matrix_3 b);
  double det();
  Matrix_3 T();
  Matrix_3 adjoint();
  double trace();
  Matrix_3 operator* (double x);
  Vec operator* (Vec v);
};