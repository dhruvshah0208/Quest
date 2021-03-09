#include<stdio.h>
#include<math.h>
typedef struct Vec Vec;
typedef struct Matrix_3 Matrix_3;

void Vec_construct(Vec *temp, double x,double y,double z);
void cross(Vec* temp,Vec* vec1,Vec* vec2);
double dot(Vec* vec1,Vec* vec2);
void scale_vec(Vec *temp,double k);
void add_vec(Vec *temp,Vec* vec1,Vec* vec2);
void matrix_construct(Matrix_3 *temp,int x);
void add_matrix(Matrix_3 *temp,Matrix_3 *m1,Matrix_3 *m2);
void scale_matrix(Matrix_3 *temp,Matrix_3 *M,double k); // This modifies the input matrix
double det(Matrix_3 *temp);
double trace(Matrix_3 *temp);
void T(Matrix_3 *temp,Matrix_3 *m1);
void adjoint(Matrix_3 *temp,Matrix_3* m);
void outer_product(Matrix_3* temp,Vec* v1,Vec* v2);
void matmul(Vec* temp,Matrix_3* M,Vec *v);

struct Vec{
  double x,y,z;
};
void Vec_construct(Vec *temp, double x,double y,double z)
{
  temp->x = x;
  temp->y = y;
  temp->z = z;
}
void cross(Vec* temp,Vec* vec1,Vec* vec2){
  temp->x = vec1->y*vec2->z - vec1->z*vec2->y;
  temp->y = vec1->z*vec2->x - vec1->x*vec2->z;
  temp->z = vec1->x*vec2->y - vec1->y*vec2->x;
}
double dot(Vec* vec1,Vec* vec2){
  return (vec1->x*vec2->x + vec1->y*vec2->y + vec1->z*vec2->z);
}
void scale_vec(Vec *temp,double k){  // This modifies the input matrix
  temp->x *= k;temp->y *= k;temp->z *= k;
}
void add_vec(Vec *temp,Vec* vec1,Vec* vec2){
  temp->x = vec1->x + vec2->x;
  temp->y = vec1->y + vec2->y;
  temp->z = vec1->z + vec2->z;
}
struct Matrix_3{
  double elements[3][3];
// Vec operator* (Vec v);
};

void matrix_construct(Matrix_3 *temp,int x)
{
  if (x == 0){ // Null Matrix
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        (temp->elements)[i][j] = 0;
      }
    }
  }
  if (x == 1){  // Identity Matrix
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        if (i == j) (temp->elements)[i][j] = 1;
        else (temp->elements)[i][j] = 0;
      }
    }
  }
}
void add_matrix(Matrix_3 *temp,Matrix_3 *m1,Matrix_3 *m2){
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (temp->elements)[i][j] = (m1->elements)[i][j] + (m2->elements)[i][j];
    }
  }
}

void scale_matrix(Matrix_3 *temp,Matrix_3 *M,double k){
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        temp->elements[i][j] = (M->elements[i][j])*k;
    }
  }
}
double det(Matrix_3 *temp){
  double x;
  x = temp->elements[0][0] * (temp->elements[1][1]*temp->elements[2][2] - temp->elements[1][2]*temp->elements[2][1])
    - temp->elements[0][1] * (temp->elements[1][0]*temp->elements[2][2] - temp->elements[1][2]*temp->elements[2][0])
    + temp->elements[0][2] * (temp->elements[1][0]*temp->elements[2][1] - temp->elements[1][1]*temp->elements[2][0]);
  return x;
}
double trace(Matrix_3 *temp){
  return (temp->elements[0][0] + temp->elements[1][1] + temp->elements[2][2]);
}
void T(Matrix_3 *temp,Matrix_3 *m1){                                                                  // Transpose Operation
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        temp->elements[i][j] = temp->elements[j][i];
      }
    }        // #REVISIT pass value by reference
  }
void adjoint(Matrix_3 *temp,Matrix_3* m){
  temp->elements[0][0] = m->elements[1][1]*m->elements[2][2] - m->elements[1][2]*m->elements[2][1];
  temp->elements[1][0] = m->elements[1][2]*m->elements[2][0] - m->elements[1][0]*m->elements[2][2];
  temp->elements[2][0] = m->elements[1][0]*m->elements[2][1] - m->elements[2][0]*m->elements[1][1];
  temp->elements[0][1] = m->elements[2][1]*m->elements[0][2] - m->elements[0][1]*m->elements[2][2];
  temp->elements[1][1] = m->elements[0][0]*m->elements[2][2] - m->elements[2][0]*m->elements[0][2];
  temp->elements[2][1] = m->elements[2][0]*m->elements[0][1] - m->elements[0][0]*m->elements[2][1];
  temp->elements[0][2] = m->elements[0][1]*m->elements[1][2] - m->elements[1][1]*m->elements[0][2];
  temp->elements[1][2] = m->elements[1][0]*m->elements[0][2] - m->elements[0][0]*m->elements[1][2];
  temp->elements[2][2] = m->elements[0][0]*m->elements[1][1] - m->elements[1][0]*m->elements[0][1];
}
void outer_product(Matrix_3* temp,Vec* v1,Vec* v2){
  temp->elements[0][0] = v1->x*v2->x;
  temp->elements[0][1] = v1->x*v2->y;
  temp->elements[0][2] = v1->x*v2->z;
  temp->elements[1][0] = v1->y*v2->x;
  temp->elements[1][1] = v1->y*v2->y;
  temp->elements[1][2] = v1->y*v2->z;
  temp->elements[2][0] = v1->z*v2->x;
  temp->elements[2][1] = v1->z*v2->y;
  temp->elements[2][2] = v1->z*v2->z;
}
void matmul(Vec* temp,Matrix_3* M,Vec *v){
  temp->x = M->elements[0][0]*v->x + M->elements[0][1]*v->y + M->elements[0][2]*v->z;
  temp->y = M->elements[1][0]*v->x + M->elements[1][1]*v->y + M->elements[1][2]*v->z;
  temp->z = M->elements[2][0]*v->x + M->elements[2][1]*v->y + M->elements[2][2]*v->z;
}

int main()
{
  int N = 5;
  Vec b[N];Vec r[N];

  double epsilon = 0.15;

  double a[N];
  a[0] = 1.0/N;a[1] = 1.0/N;a[2] = 1.0/N;a[3] = 1.0/N;a[4] = 1.0/N;
  double guess = 1;
  Matrix_3 B,S;
  matrix_construct(&B,0);
  matrix_construct(&S,0);  
  Vec z;
  Vec_construct(&z,0,0,0);
  int counter = 0;
  double K;
  double norm_2_z;
  double det_S;
  double trace_B;
  Matrix_3 adj_S;
  Matrix_3 check;
  Vec_construct(&(b[0]),0.4242640687,
  0.8434513391,
  0.3295297233);
  Vec_construct(&b[1],0.2672612419,
  -0.06201818702,
  0.9616263167);
  Vec_construct(&b[2],-0.32929278,
  0.8594680978,
  0.3909998131);
  Vec_construct(&b[3],-0.4364357805,
  0.247453544,
  0.8650378911);
  Vec_construct(&b[4],-0.1632993162,
  0.09551972887,
  0.9819416046);

  Vec_construct(&r[0],0.4242640687,
  0.5656854249,
  0.7071067812);
  Vec_construct(&r[1],0.2672612419,
  -0.5345224838,
  0.8017837257);
  Vec_construct(&r[2],-0.32929278,
  0.5488212999,
  0.7683498199);
  Vec_construct(&r[3],-0.4364357805,
  -0.2182178902,
  0.8728715609);
  Vec_construct(&r[4],-0.1632993162,
  -0.4082482905,
  0.898146239);
  // Perform Quest
  do {
    if(counter == 1){ // Perform rotation abt x --- {x,-y,-z}
      for (int i = 0; i < N; i++){
        r[i].y = -1*r[i].y;
        r[i].z = -1*r[i].z;
      }
    }
    if(counter == 2){ // Perform rotation abt y --- {-x,y,-z}
      for (int i = 0; i < N; i++){
        r[i].x = -1*r[i].x;
        r[i].z = -1*r[i].z;
      }
    }
    if(counter == 3){ // Perform rotation abt z --- {-x,-y,z}
      for (int i = 0; i < N; i++){
        r[i].y = -1*r[i].y;
        r[i].x = -1*r[i].x;
      }
    }

    for (int i = 0; i < N; i++){
       Matrix_3 temp;
       outer_product(&temp,&(b[i]),&(r[i]));
       scale_matrix(&temp,&temp,a[i]);
       add_matrix(&B,&B,&temp);
    }
    for (int i = 0; i < N; i++)
    {
      Vec temp;
      cross(&temp,&(b[i]),&(r[i]));
      scale_vec(&temp,a[i]);
      add_vec(&z,&z,&temp);  
    }
    //rho = lambda_max(guess,B,z) + B.trace();

    trace_B = trace(&B);
    T(&S,&B);
    add_matrix(&S,&S,&B); // construct S
    adjoint(&adj_S,&S);
    K = trace(&adj_S);
    norm_2_z = z.x*z.x + z.y*z.y + z.z*z.z;
    det_S = det(&S);
    Vec temp;
    matmul(&temp,&S,&z);
    double zTSz = dot(&z,&temp);
    matmul(&temp,&S,&temp);
    double zTSSz = dot(&z,&temp);

    double rho = guess; // sum of weights
    double f(){
      double first = (rho*rho - trace_B*trace_B + K)*(rho*rho - trace_B*trace_B - norm_2_z);
      double second = (rho - trace_B)*(zTSz + det_S);
      return (first - second - zTSSz);
    }
    double f_bar(){
      double temp = 2*rho*(2*rho*rho - 2*trace_B*trace_B + K - norm_2_z);
      return (temp - zTSz - det_S);
    }
    for (int i = 0; i < 5; i++) {
      rho -= f()/f_bar();
    }
    rho += trace_B;
    Matrix_3 I;
    matrix_construct(&I,1);
    scale_matrix(&I,&I,rho);
    scale_matrix(&check,&S,-1);
    add_matrix(&check,&check,&I);

    counter++;
  } while((det(&check)) < epsilon);

  Vec q_3;
  double q_4 = det(&check);
  adjoint(&check,&check);
  matmul(&q_3,&check,&z); //  (I3*rho - (B + B.T())).adjoint() * z

  if(counter == 2){
    Vec e;
    Vec_construct(&e,1,0,0);
    double new_4 = -dot(&e,&q_3);
    cross(&q_3,&e,&q_3);
    scale_vec(&e,q_4);
    add_vec(&q_3,&q_3,&e);
    q_4 = new_4;
  }
  else if(counter == 3){
    Vec e;
    Vec_construct(&e,0,1,0);
    double new_4 = -dot(&e,&q_3);
    cross(&q_3,&e,&q_3);
    scale_vec(&e,q_4);
    add_vec(&q_3,&q_3,&e);
    q_4 = new_4;
  }
  else if(counter == 4){
    Vec e;
    Vec_construct(&e,0,0,1);
    double new_4 = -dot(&e,&q_3);
    cross(&q_3,&e,&q_3);
    scale_vec(&e,q_4);
    add_vec(&q_3,&q_3,&e);
    q_4 = new_4;
  }
  double alpha; 
  alpha = sqrt((q_3.x)*(q_3.x) + (q_3.y)*(q_3.y) + (q_3.z)*(q_3.z) + (q_4)*(q_4));
  
  q_3.x /= alpha;
  q_3.y /= alpha;
  q_3.z /= alpha;
  q_4 /= alpha;
  printf("q_3.x = %lf \n", q_3.x);    
  printf("q_3.y = %lf \n", q_3.y);    
  printf("q_3.z = %lf \n", q_3.z);    
  printf("q_4 = %lf \n", q_4);    
}