#include <iostream>
#include <vector>
#include <math.h>
#include <random>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <sstream>
#include <string>
#include "stl.h"

#define PI 3.14159265

using namespace std;


//sphC = np.asarray([0.000125, 0.000125, 0.0004]) # center of the sphere
double r = (0.000080)/3.0; // radius of the sphere
double del_tw = 0.000005; // time step for the walker 
double alpha_s = 0.279/(1200.*3200.); // thermal diffusivity

// Generate a random number between 0 and 1
// return a uniform number in [0,1].
double unifRand()
{
    return rand() / double(RAND_MAX);
}

// Reset the random number generator with the system clock.
void seed()
{
    srand(time(0));
}


double normalRand()
{
  std::random_device rd;
  float sample;
  // Mersenne twister PRNG, initialized with seed from previous random device instance
  std::mt19937 gen(rd());
  std::normal_distribution<float> d(0.0, 1.0);
  sample = d(gen);
  return sample;
  //std::cout << sample << std::endl;
}

// passing the vector by reference

void randPtPick(vector<double> vt1, vector<double> vt2, vector<double> vt3, vector<double>& d)
{
  float a;
  float b;

  seed();
  a = unifRand();
  b = unifRand();

  d[0] = (1.0-sqrt(a))*vt1[0] + sqrt(a)*(1.0-b)*vt2[0] + b*sqrt(a)*vt3[0];
  d[1] = (1.0-sqrt(a))*vt1[1] + sqrt(a)*(1.0-b)*vt2[1] + b*sqrt(a)*vt3[1];
  d[2] = (1.0-sqrt(a))*vt1[2] + sqrt(a)*(1.0-b)*vt2[2] + b*sqrt(a)*vt3[2];
  //cout << d[0] << endl;
  //cout << d[1] << endl;
  //cout << d[2] << endl;

}

/* Source: https://stackoverflow.com/questions/17333/what-is-the-most-effective-way-for-float-and-double-comparison */

bool approximatelyEqual(float a, float b, float epsilon)
{
    return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

void vectorSub(vector<double> v1, vector<double> v2, vector<double>& vs)
{
    vs[0] = v1[0]-v2[0];
    vs[0] = v1[1]-v2[1];
    vs[0] = v1[2]-v2[2];  
}

void vectorCrossProd(vector<double> v1, vector<double> v2, vector<double>& vs)
{
  vs[0] = v1[1] * v2[2] - v1[2] * v2[1];
  vs[1] = v1[2] * v2[0] - v1[2] * v2[2];
  vs[2] = v1[0] * v2[1] - v1[1] * v2[0];
 
}

double vectorNorm(vector<double> v1)
{
  double vnorm;
  vnorm = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]) ; 
  return vnorm;
}


double areaTri(vector<double> vt1, vector<double> vt2, vector<double> vt3)
{
  vector<double> vt2_vt1;
  vector<double> vt3_vt1;
  vector<double> vcp;

  vectorSub(vt2, vt1, vt2_vt1);
  vectorSub(vt3, vt1, vt3_vt1);
  vectorCrossProd(vt2_vt1, vt3_vt1, vcp);

  return vectorNorm(vcp);
 
}

bool isPtInTri (vector<double> Pt, vector<double> vt1, vector<double> vt2, vector<double> vt3)
{
    double area2, alpha, beta, gamma, sum_weights;

    area2 = areaTri(vt1, vt2, vt3);
    alpha = areaTri(Pt, vt2, vt3) / area2 ;
    beta = areaTri(Pt, vt2, vt3) / area2 ;
    gamma = areaTri(Pt, vt1, vt2) / area2 ;
    
    sum_weights = alpha + beta + gamma; //this should be sent to the calling function
    // TODO: make sum_weights a pointer
    if (alpha>=0.0 && beta>=0.0 && gamma>=0.0 && approximatelyEqual(sum_weights, 1.0, 0.1))
      return true;
    else
      return false;
     
}

void vectorUnit(vector<double>& v1)
{
  v1[0] = v1[0]/vectorNorm(v1);
  v1[1] = v1[1]/vectorNorm(v1);
  v1[2] = v1[2]/vectorNorm(v1);
}

double vectorDotProd(vector<double> v1, vector<double> v2)
{
   double dp = 0.0;
   dp = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
   return dp;
}

void matrixMul(double firstMatrix[][3], double secondMatrix[][3], double mult[][3], int rowFirst, int columnFirst, int rowSecond, int columnSecond)
{
  int i, j, k;

  // Initializing elements of matrix mult to 0.
  for(i = 0; i < rowFirst; ++i)
   for(j = 0; j < columnSecond; ++j)
    {
      mult[i][j] = 0;
    }

  // Multiplying matrix firstMatrix and secondMatrix and storing in array mult.
  for(i = 0; i < rowFirst; ++i)
   for(j = 0; j < columnSecond; ++j)
    for(k=0; k<columnFirst; ++k)
     {
       mult[i][j] += firstMatrix[i][k] * secondMatrix[k][j];
     }
}

double vecDotProduct(const double *x, const double *y, int n)
{
    double res = 0.0;
    int i; 
    for (i = 0; i < n; i++)
    {   
        res += x[i] * y[i];
    }
    return res;
}

void matVectorMul(const double mat[][3], const double vec[], vector<double>& result, int rows, int cols)
{ // in matrix form: result = mat * vec;
    int i;
    for (i = 0; i < rows; i++)
    {
        result[i] = vecDotProduct(mat[i], vec, 3);
    }
}

void diffusion_dis(vector<double> a, vector<double>& xn)
{
   vector<double> del_xd={0,0,0};

   double r1 = normalRand();
   double r2 = normalRand();
   double r3 = normalRand();
   
   del_xd[0] = sqrt (2.0*alpha_s*del_tw)*r1;
   del_xd[1] = sqrt (2.0*alpha_s*del_tw)*r2;
   del_xd[2] = sqrt (2.0*alpha_s*del_tw)*r3;
   cout << del_xd[0] << endl;
   cout << del_xd[1] << endl;
   cout << del_xd[2] << endl;
   double ds = vectorNorm(del_xd);

   /* convert to the unit vector */
   a[0] = a[0]/vectorNorm(a);
   a[1] = a[1]/vectorNorm(a);
   a[2] = a[2]/vectorNorm(a);

   /* define north pole */
   vector<double> b(3,0);  //length 3 and initialize as 0
   vector<double> v(3,0);

   vectorCrossProd(a, b, v); 
   double s = vectorNorm(v);
   double cp = vectorDotProd(a,b);

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
   double phi_p = ds/r;
   double theta_p = 2.0*PI*unifRand();
   //print phi_p, theta_p 
   // new point from the pole
   double xn_p[3]={0,0,0};
 
   xn_p[0] = 2.0;//cos(theta_p)*sin(phi_p);
   xn_p[1] = sin(theta_p)*sin(phi_p);
   xn_p[2] = cos(phi_p);

   // new point moved back using the inverse rotation
   // This is what gets "returned"
   //xn = r*np.matrixMul(np.linalg.inv(R), xn_p)
   
  matVectorMul(R, xn_p, xn, 3, 1);
 
}

int main(int argc, char** argv) 
{
  vector<double> p1={1,2,3};
  vector<double> p2={2,3,4};
  vector<double> p3={3,3,4};
  vector<double> randPt(3);
  
  randPtPick(p1, p2, p3, randPt);
  //cout << approximatelyEqual(0.0345, 0.03467, 0.01)  << endl;
  diffusion_dis(p1, randPt);

  c_Tri_Element *Geom_Panels;
  Geom_Panels = new c_Tri_Element[1];

  Allocate_geom(Geom_Panels); 
  read_geometry_STL(Geom_Panels);
  cout << Geom_Panels[0].vert1[0].x[0] << endl;
  //vectorNorm(te);
  
  return 0;
}
