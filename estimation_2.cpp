#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <boost/algorithm/string.hpp>
using namespace std;
#include "newton_raphson.cpp"
#include <cmath>
#include <stdlib.h>
/*
long convert(string s){

  long x = 0;
  int i;
  // Traverse the string
  // Assumming the numbers are always less than 1
  if(s[0] == '-'){
    for (i = 3;i < 13; i++){
      if (i < s.length())
      {
        int t = (int)(s[i]) - 48;
        x = x * 10 + t;
      }
      else
        x = x * 10;
    }
    x = -1*x;
  }
  else{
    for (i = 2;i < 12; i++){
      if (i < s.length())
      {
        int t = (int)(s[i]) - 48;
        x = x * 10 + t;
      }
      else
        x = x * 10;
    }
  }
  return x;

}
*/
///*
double convert(string s){

  double x = 0,j;
  int i;
  // Traverse the string
  // Assumming the numbers are always less than 1
  if(s[0] == '-'){
    for (i = 3,j = 0.1;i < s.length(); i++,j/=10){
      int t = (int)(s[i]) - 48;
      x = x + t*j;
    }
    x = -1*x;
  }
  else{
    for (i = 2,j = 0.1;i < s.length(); i++,j/=10){
      int t = (int)(s[i]) - 48;
      x = x + t*j;
    }
  }
  return x;
}
//*/

int main()
{
  string fileName;
  cin >> fileName;
  ifstream file(fileName);
  string line = "";
  int i = 0;
  int N = 5;
  Vec b[N];Vec r[N];

  double epsilon = 0.15;

// Iterate through each line and split the content using delimeter
while (getline(file, line))
{

    vector<string> vec;
    boost::algorithm::split(vec, line, boost::is_any_of(","));
    //
    for (int j = 0; j < 2*N; j++) {
      if(i == 0){
        if (j < N){
          r[j].x = convert(vec[j]);
        }
        else
          b[j-N].x = convert(vec[j]);
      }
      else if(i == 1){
        if (j < N)
          r[j].y = convert(vec[j]);
        else
          b[j-N].y = convert(vec[j]);
      }
      if(i == 2){
        if (j < N)
          r[j].z = convert(vec[j]);
        else
          b[j-N].z = convert(vec[j]);
      }
    }
    i++;
}
  // Close the File
  file.close();

  double a[N];
  a[0] = 1.0/N;a[1] = 1.0/N;a[2] = 1.0/N;a[3] = 1.0/N;a[4] = 1.0/N;
  double guess = 1;
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
