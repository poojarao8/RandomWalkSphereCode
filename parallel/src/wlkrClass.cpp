#include <iostream>
#include <random>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <sstream>
#include <string>
#include <fstream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <iterator>
#include "main.h"

using namespace std;

void diffusion_dis(double* a)
{
   double del_xd[3]={0,0,0};

   double r1 = normalRand();
   double r2 = normalRand();
   double r3 = normalRand();

   del_xd[0] = sqrt (2.0*alpha_s*del_tw)*r1;
   del_xd[1] = sqrt (2.0*alpha_s*del_tw)*r2;
   del_xd[2] = 0.0; // sqrt (2.0*alpha_s*del_tw)*r3;
   double ds = vectorNorm(del_xd);

   // convert to the unit vector
   double norm_a = vectorNorm(a);
   a[0] = a[0]/norm_a;
   a[1] = a[1]/norm_a;
   a[2] = a[2]/norm_a;

   // define north pole 
   double b[3]={0, 0, 1};  //length 3 and initialize as 0
   double v[3]={};

   vectorCrossProd(a, b, v);
   double s = vectorNorm(v);
   double cp = vectorDotProd(a,b,3);

   double V[3][3]={}; //Initialize to 0
   double I[3][3]={1,0,0,0,1,0,0,0,1}; //Identity matrix
   double R[3][3]={}; //Rotation matrix

   V[0][0] = 0.0;
   V[0][1] = -v[2];
   V[0][2] =  v[1];

   V[1][0] = -V[0][1];
   V[1][1] = 0.0;
   V[1][2] = -v[0];

   V[2][0] = -V[0][2];
   V[2][1] = -V[1][2];
   V[2][2] = 0.0;

   if (approximatelyEqual(cp, -1.0, 1e-4))
    {
      for (int ii=0; ii<3; ++ii)
      for (int jj=0; jj<3; ++jj)
       {
         R[ii][jj] = I[ii][jj];
       }

       R[2][2] = -1.0;
    }
   else
    {
      double V2[3][3]={};
      matrixMul(V, V, V2, 3, 3, 3, 3);
      for (int iii=0; iii<3; ++iii)
      for (int jjj=0; jjj<3; ++jjj)
       {
         R[iii][jjj] = I[iii][jjj] + V[iii][jjj] + (1.0/(1.0+ cp))*V2[iii][jjj];
       }

    }

   //calculate the arc length
   double phi_p = ds/rad;
   //double ru =0.1177474;
   double theta_p = 2.0*PI*unifRand();
   //double theta_p = 2.0*PI*ru;
   //print phi_p, theta_p 
   // new point from the pole
   double xn_p[3]={0,0,0};

   xn_p[0] = cos(theta_p)*sin(phi_p);
   xn_p[1] = sin(theta_p)*sin(phi_p);
   xn_p[2] = cos(phi_p);

   // new point moved back using the inverse rotation
   double invR[3][3];
   MatrixInversion(R, 3, invR);
   matVectorMul(invR, xn_p, a, 3, 1);

   a[0] = rad*a[0];
   a[1] = rad*a[1];
   a[2] = rad*a[2];

  //cout << a[0] << endl;
  //cout << a[1] << endl;
  //cout << a[2] << endl;

}


void Walker::moveWalker(double* a)
{
  diffusion_dis(a);
}


