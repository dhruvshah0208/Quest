#include <iostream>
#include "Structures.h"
using namespace std;

Vec::Vec(double p,double q,double r){
    x = p;y = q;z = r;
  }
Vec Vec::cross(Vec &b){
    Vec v;
    v.x = y*b.z - z*b.y;
    v.y = z*b.x - x*b.z;
    v.z = x*b.y - y*b.x;
    return v;
  }
double Vec::dot(Vec b){         // Inner Product
    return (x*b.x + y*b.y + z*b.z);
  }

Vec Vec::operator* (double &t){
  Vec c;
  c.x = x*t;
  c.y = y*t;
  c.z = z*t;
  return c;
}
Vec Vec::operator+ (Vec b){
  Vec c;
  c.x = x + b.x;
  c.y = y + b.y;
  c.z = z + b.z;
  return c;
}
Vec Vec::operator- (Vec b){
  Vec c;
  c.x = x - b.x;
  c.y = y - b.y;
  c.z = z - b.z;
  return c;
}


Matrix_3::Matrix_3(int x){
    if (x == 0){ // Null Matrix
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          elements[i][j] = 0;
        }
      }
    }
    if (x == 1){  // Identity Matrix
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          if (i == j) elements[i][j] = 1;
          else  elements[i][j] = 0;
        }
      }
    }
  }

double& Matrix_3:: operator ()(int i,int j){
    return elements[i][j];
  }                       // Accessing the elements
void Matrix_3::assign(int i,int j,double z){      // Assigning the elements
    elements[i][j] = z;
  }                       // Accessing the elements

Matrix_3 Matrix_3::operator*(double x){
  Matrix_3 c;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        c.assign(i,j,elements[i][j]*x);
    }
  }
  return c;
}                       // Accessing the elements


Vec Matrix_3::operator* (Vec v){
    Vec b;
    b.x = elements[0][0]*v.x + elements[0][1]*v.y + elements[0][2]*v.z;
    b.y = elements[1][0]*v.x + elements[1][1]*v.y + elements[1][2]*v.z;
    b.z = elements[2][0]*v.x + elements[2][1]*v.y + elements[2][2]*v.z;

    return b;
  }
Matrix_3 Matrix_3::operator+ (Matrix_3 b){
    Matrix_3 c;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        c.assign(i,j,elements[i][j] + b(i,j));
      }
    }
    return c;
  }

Matrix_3 Matrix_3::operator- (Matrix_3 b){
    Matrix_3 c;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        c.assign(i,j,elements[i][j] - b(i,j));
      }
    }
    return c;
  }

double Matrix_3::det(){
    double x;
    x = elements[0][0] * (elements[1][1]*elements[2][2] - elements[1][2]*elements[2][1])
      - elements[0][1] * (elements[1][0]*elements[2][2] - elements[1][2]*elements[2][0])
      + elements[0][2] * (elements[1][0]*elements[2][1] - elements[1][1]*elements[2][0]);
    return x;
  }

double Matrix_3::trace(){
  return (elements[0][0] + elements[1][1] + elements[2][2]);
}

Matrix_3 Matrix_3::T(){                                                                  // Transpose Operation
    Matrix_3 b;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        b.assign(i,j,elements[j][i]);
      }
    }
    return b;        // #REVISIT pass value by reference
  }

Matrix_3 Matrix_3::adjoint(){
    Matrix_3 b; // Transpose has been taken care of
    b.assign(0,0,elements[1][1]*elements[2][2] - elements[1][2]*elements[2][1]);
    b.assign(1,0,elements[1][2]*elements[2][0] - elements[1][0]*elements[2][2]);
    b.assign(2,0,elements[1][0]*elements[2][1] - elements[2][0]*elements[1][1]);
    b.assign(0,1,elements[2][1]*elements[0][2] - elements[0][1]*elements[2][2]);
    b.assign(1,1,elements[0][0]*elements[2][2] - elements[2][0]*elements[0][2]);
    b.assign(2,1,elements[2][0]*elements[0][1] - elements[0][0]*elements[2][1]);
    b.assign(0,2,elements[0][1]*elements[1][2] - elements[1][1]*elements[0][2]);
    b.assign(1,2,elements[1][0]*elements[0][2] - elements[0][0]*elements[1][2]);
    b.assign(2,2,elements[0][0]*elements[1][1] - elements[1][0]*elements[0][1]);

    return b;
  }

Matrix_3 outer_product(const Vec &a,const Vec &b){ // REVISIT Pass by reference
    Matrix_3 z;
    z.assign(0,0,a.x*b.x);
    z.assign(0,1,a.x*b.y);
    z.assign(0,2,a.x*b.z);
    z.assign(1,0,a.y*b.x);
    z.assign(1,1,a.y*b.y);
    z.assign(1,2,a.y*b.z);
    z.assign(2,0,a.z*b.x);
    z.assign(2,1,a.z*b.y);
    z.assign(2,2,a.z*b.z);
    return z;
  } // Giving Transposed Result
/*
int main(){
  Vec r1 = Vec(1,2,3);
  Vec b1 = Vec(4,1,2);
  Vec z;
  double x = 2;
  for (int i = 0; i < 3; i++) z = z + (b1.cross(r1))*x; // Optimization Done

}
*/
