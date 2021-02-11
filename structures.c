#include<stdio.h>
typedef struct Vec Vec;
typedef struct Matrix_3 Matrix_3;

void Vec_construct(Vec *this, double x,double y,double z);
void cross(Vec* this,Vec* vec1,Vec* vec2);
double dot(Vec* vec1,Vec* vec2);
void scale_vec(Vec *this,double k);
void add_vec(Vec *this,Vec* vec1,Vec* vec2);
void matrix_construct(Matrix_3 *this,unsigned char x);
void add_matrix(Matrix_3 *this,Matrix_3 *m1,Matrix_3 *m2);
void scale_matrix(Matrix_3 *this,Matrix_3 *M,double k); // This modifies the input matrix
double det(Matrix_3 *this);
double trace(Matrix_3 *this);
void T(Matrix_3 *this,Matrix_3 *m1);
void adjoint(Matrix_3 *this,Matrix_3* m);
void outer_product(Matrix_3* this,Vec* v1,Vec* v2);
void matmul(Vec* this,Matrix_3* M,Vec *v);

struct Vec{
  double x,y,z;
};
void Vec_construct(Vec *this, double x,double y,double z)
{
  this->x = x;
  this->y = y;
  this->z = z;
}
void cross(Vec* this,Vec* vec1,Vec* vec2){
  this->x = vec1->y*vec2->z - vec1->z*vec2->y;
  this->y = vec1->z*vec2->x - vec1->x*vec2->z;
  this->z = vec1->x*vec2->y - vec1->y*vec2->x;
}
double dot(Vec* vec1,Vec* vec2){
  return (vec1->x*vec2->x + vec1->y*vec2->y + vec1->z*vec2->z);
}
void scale_vec(Vec *this,double k){  // This modifies the input matrix
  this->x *= k;this->y *= k;this->z *= k;
}
void add_vec(Vec *this,Vec* vec1,Vec* vec2){
  this->x = vec1->x + vec1->x;
  this->y = vec1->y + vec1->y;
  this->z = vec1->z + vec1->z;
}
struct Matrix_3{
  double elements[3][3];
// Vec operator* (Vec v);
};

void matrix_construct(Matrix_3 *this,unsigned char x)
{
  if (x == 0){ // Null Matrix
    for (unsigned char i = 0; i < 3; i++) {
      for (unsigned char j = 0; j < 3; j++) {
        (this->elements)[i][j] = 0;
      }
    }
  }
  if (x == 1){  // Identity Matrix
    for (unsigned char i = 0; i < 3; i++) {
      for (unsigned char j = 0; j < 3; j++) {
        if (i == j) (this->elements)[i][j] = 1;
        else (this->elements)[i][j] = 0;
      }
    }
  }
}
void add_matrix(Matrix_3 *this,Matrix_3 *m1,Matrix_3 *m2){
  for (unsigned char i = 0; i < 3; i++) {
    for (unsigned char j = 0; j < 3; j++) {
      (this->elements)[i][j] = (m1->elements)[i][j] + (m2->elements)[i][j];
    }
  }
}

void scale_matrix(Matrix_3 *this,Matrix_3 *M,double k){
  for (unsigned char i = 0; i < 3; i++) {
    for (unsigned char j = 0; j < 3; j++) {
        this->elements[i][j] = (M->elements[i][j])*k;
    }
  }
}
double det(Matrix_3 *this){
  double x;
  x = this->elements[0][0] * (this->elements[1][1]*this->elements[2][2] - this->elements[1][2]*this->elements[2][1])
    - this->elements[0][1] * (this->elements[1][0]*this->elements[2][2] - this->elements[1][2]*this->elements[2][0])
    + this->elements[0][2] * (this->elements[1][0]*this->elements[2][1] - this->elements[1][1]*this->elements[2][0]);
  return x;
}
double trace(Matrix_3 *this){
  return (this->elements[0][0] + this->elements[1][1] + this->elements[2][2]);
}
void T(Matrix_3 *this,Matrix_3 *m1){                                                                  // Transpose Operation
    for (unsigned char i = 0; i < 3; i++) {
      for (unsigned char j = 0; j < 3; j++) {
        this->elements[i][j] = this->elements[j][i];
      }
    }        // #REVISIT pass value by reference
  }
void adjoint(Matrix_3 *this,Matrix_3* m){
  this->elements[0][0] = m->elements[1][1]*m->elements[2][2] - m->elements[1][2]*m->elements[2][1];
  this->elements[1][0] = m->elements[1][2]*m->elements[2][0] - m->elements[1][0]*m->elements[2][2];
  this->elements[2][0] = m->elements[1][0]*m->elements[2][1] - m->elements[2][0]*m->elements[1][1];
  this->elements[0][1] = m->elements[2][1]*m->elements[0][2] - m->elements[0][1]*m->elements[2][2];
  this->elements[1][1] = m->elements[0][0]*m->elements[2][2] - m->elements[2][0]*m->elements[0][2];
  this->elements[2][1] = m->elements[2][0]*m->elements[0][1] - m->elements[0][0]*m->elements[2][1];
  this->elements[0][2] = m->elements[0][1]*m->elements[1][2] - m->elements[1][1]*m->elements[0][2];
  this->elements[1][2] = m->elements[1][0]*m->elements[0][2] - m->elements[0][0]*m->elements[1][2];
  this->elements[2][2] = m->elements[0][0]*m->elements[1][1] - m->elements[1][0]*m->elements[0][1];
}
void outer_product(Matrix_3* this,Vec* v1,Vec* v2){
  this->elements[0][0] = v1->x*v2->x;
  this->elements[0][1] = v1->x*v2->y;
  this->elements[0][2] = v1->x*v2->z;
  this->elements[1][0] = v1->y*v2->x;
  this->elements[1][1] = v1->y*v2->y;
  this->elements[1][2] = v1->y*v2->z;
  this->elements[2][0] = v1->z*v2->x;
  this->elements[2][1] = v1->z*v2->y;
  this->elements[2][2] = v1->z*v2->z;
}
void matmul(Vec* this,Matrix_3* M,Vec *v){
  this->x = M->elements[0][0]*v->x + M->elements[0][1]*v->y + M->elements[0][2]*v->z;
  this->y = M->elements[1][0]*v->x + M->elements[1][1]*v->y + M->elements[1][2]*v->z;
  this->z = M->elements[2][0]*v->x + M->elements[2][1]*v->y + M->elements[2][2]*v->z;
}
