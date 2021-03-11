#include <iostream>
#include "newton_raphson.cpp"
using namespace std;
#include <cmath>

int main(){
  // INPUT FORMAT is a_i,b_i,r_i
  int N ;  // #REVISIT - No of stars
  cin >> N;
  double a[N];

  double epsilon = 0.15;

  a[0] = 1.0/N;
  a[1] = 1.0/N;
  a[2] = 1.0/N;
  a[3] = 1.0/N;
  a[4] = 1.0/N;

  double guess = 1;
  Vec b[5];Vec r[5];
  b[0] = Vec(0.4242640687,
  0.8434513391,
  0.3295297233);
  b[1] = Vec(0.2672612419,
  -0.06201818702,
  0.9616263167);
  b[2] = Vec(-0.32929278,
  0.8594680978,
  0.3909998131);
  b[3] = Vec(-0.4364357805,
  0.247453544,
  0.8650378911);
  b[4] = Vec(-0.1632993162,
  0.09551972887,
  0.9819416046);

  r[0] = Vec(0.4242640687,
  0.5656854249,
  0.7071067812);
  r[1] = Vec(0.2672612419,
  -0.5345224838,
  0.8017837257);
  r[2] = Vec(-0.32929278,
  0.5488212999,
  0.7683498199);
  r[3] = Vec(-0.4364357805,
  -0.2182178902,
  0.8728715609);
  r[4] = Vec(-0.1632993162,
  -0.4082482905,
  0.898146239);

  Matrix_3 B(0);
  Vec z;
  Matrix_3 I3(1);
  double rho = 0;
  int counter = 0;
  // Perform Quest
  do {
    if(counter == 1){ // Perform rotation abt x --- {x,-y,-z}
      B = Matrix_3(0);
      z = Vec();
      for (int i = 0; i < N; i++){
        r[i].y = -1*r[i].y;
        r[i].z = -1*r[i].z;
      }
    }
    if(counter == 2){ // Perform rotation abt y --- {-x,y,-z}
      B = Matrix_3(0);
      z = Vec();
      for (int i = 0; i < N; i++){
        r[i].x = -1*r[i].x;
        r[i].z = -1*r[i].z;
      }
    }
    if(counter == 3){ // Perform rotation abt z --- {-x,-y,z}
      B = Matrix_3(0);
      z = Vec();
      for (int i = 0; i < N; i++){
        r[i].y = -1*r[i].y;
        r[i].x = -1*r[i].x;
      }
    }

    for (int i = 0; i < N; i++) B = B + (outer_product(b[i],r[i]))*a[i]; // Optimization Done
    for (int i = 0; i < N; i++) z = z + (b[i].cross(r[i]))*a[i]; // Optimization Done
    rho = lambda_max(guess,B,z) + B.trace();
    counter++;
  } while((I3*rho - (B + B.T())).det() < epsilon);

  Vec q_3 = (I3*rho - (B + B.T())).adjoint() * z;  // Test #REVISIT
  double q_4 = (I3*rho - (B + B.T())).det();
  double alpha = sqrt((q_3.x)*(q_3.x)+ (q_3.y)*(q_3.y) + (q_3.z)*(q_3.z) + (q_4)*(q_4));

  if(counter == 2){
    Vec e1(1,0,0);
    double new_4 = -e1.dot(q_3);
    q_3 = e1*q_4 + e1.cross(q_3);
    q_4 = new_4;

  }
  else if(counter == 3){
    Vec e2(0,1,0);
    double new_4 = -e2.dot(q_3);
    q_3 = e2*q_4 + e2.cross(q_3);
    q_4 = new_4;

  }
  else if(counter == 4){
    Vec e3(0,0,1);
    double new_4 = -e3.dot(q_3);
    q_3 = e3*q_4 + e3.cross(q_3);
    q_4 = new_4;

  }
  cout << q_3.x/alpha << endl;
  cout << q_3.y/alpha << endl;
  cout << q_3.z/alpha << endl;
  cout << q_4/alpha << endl;
}
