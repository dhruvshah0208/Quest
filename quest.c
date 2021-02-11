#include <stdio.h>
#include "structures.c"
#include <math.h>
int main()
{
  unsigned char N = 5;
  Vec b[N];Vec r[N];

  double epsilon = 0.15;

  double a[N];
  a[0] = 1.0/N;a[1] = 1.0/N;a[2] = 1.0/N;a[3] = 1.0/N;a[4] = 1.0/N;
  double guess = 1;
  Matrix_3 B,S;
  Vec z;
  unsigned char counter;
  double K;
  double norm_2_z;
  double det_S;
  double trace_B;
  Matrix_3 adj_S;
  Matrix_3 check;
  // Perform Quest
  do {
    if(counter == 1){ // Perform rotation abt x --- {x,-y,-z}
      for (unsigned char i = 0; i < N; i++){
        r[i].y = -1*r[i].y;
        r[i].z = -1*r[i].z;
      }
    }
    if(counter == 2){ // Perform rotation abt y --- {-x,y,-z}
      for (unsigned char i = 0; i < N; i++){
        r[i].x = -1*r[i].x;
        r[i].z = -1*r[i].z;
      }
    }
    if(counter == 3){ // Perform rotation abt z --- {-x,-y,z}
      for (unsigned char i = 0; i < N; i++){
        r[i].y = -1*r[i].y;
        r[i].x = -1*r[i].x;
      }
    }

    for (unsigned char i = 0; i < N; i++){
       Matrix_3 temp;
       outer_product(&temp,&(b[i]),&(r[i]));
       scale_matrix(&temp,&temp,a[i]);
       add_matrix(&B,&B,&temp);
    }
    for (unsigned char i = 0; i < N; i++)
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
    for (unsigned char i = 0; i < 5; i++) {
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


}
