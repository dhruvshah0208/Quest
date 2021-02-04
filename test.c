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
  Vec_construct(&a,2,3,4);
  printf("%lf\n",v1->x );
  printf("Hello World");
}
  this->x = M->elements[0][0]*v->x + M->elements[0][1]*v->y + M->elements[0][2]*v->z;
  this->y = M->elements[1][0]*v->x + M->elements[1][1]*v->y + M->elements[1][2]*v->z;
  this->z = M->elements[2][0]*v->x + M->elements[2][1]*v->y + M->elements[2][2]*v->z;
