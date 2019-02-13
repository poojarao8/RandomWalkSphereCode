#include<fstream>
#include<math.h>
#include<iostream>
#include<sstream>
#include<string>
#include<vector>
#include<map>

using namespace std;

#define PI 3.14159265

extern const double rad; // center of the sphere
extern const double rad; // radius of the sphere
extern const double del_tw; // time step for the walker 
extern const double alpha_s; //thermal diffusivity

extern int glo_num_primitives;
extern int glo_panels_per_primitive;
extern int glo_total_panels;
 
struct dim
{
  double x[3];
};

class Walker
{
  public:
  double initPos[3]; // may not be needed anymore
  double curPos[3];
  int posInd;
  int tag;

  void moveWalker(double* a);
};

class c_Tri_Element
{
  public :
     dim *vert1;
     dim *vert2;
     dim *vert3;
     dim *edge1;
     dim *edge2;
     dim *norm;
     //double *temp;
     void tempDistribution(double rTemp, int worldSize, double* Temp);
     void ComputeEdge1(int i_pan);
     void ComputeEdge2(int i_pan);
     void ComputeNormal(int i_pan);
     void RandPtPick(int i_pan, double*);
     bool isPtInTri (int i_pan, double* Pt, double& sum_weights);
     int indexLocate(double* position, int i_pan, std::map<int, 
                     vector<int>>& ntris);
     void calVolumeSphere(double* dVol);
     void nearTris(map<int, vector<int>>& ntris);

     void writeInSTL(string filename,  int glo_num_subparts, 
                     int*&glo_panels_per_subpart, 
                     int*& primitiveSubpartRemove, 
                     string*& subpartStringPrim, double* temp);

     int adjNumWlkrs(double TMaxCell, int i_maxPanel, double curTemp,                               double wTemp, vector<Walker> &listofwalkers, 
                     int numProcs, int myRank);

};


// functions in randNumGen.cpp
void seed();
double unifRand();
double normalRand();

// functions in vecMatOps.cpp
bool approximatelyEqual(float a, float b, float epsilon);
bool vectorEqual(double* v1, double* v2);
void vectorSubtract(double* v1, double* v2, double* vs);
void vectorScalarMul(double* v, double s, double* sv);
double vectorNorm(double* v1);
double vectorDotProd(const double *x, const double *y, int n);
void vectorCrossProd(double* v1, double* v2, double* vs);

void matVectorMul(const double mat[][3], const double vec[], 
                  double* result, int rows, int cols);

void matrixMul(double firstMatrix[][3], double secondMatrix[][3], 
               double mult[][3], int rowFirst, int columnFirst, 
               int rowSecond, int columnSecond);

void MatrixInversion(double A[][3], int order, double Y[][3]);

// functions in geometry.cpp
void sphereTrans(double* v1, double* vt);
double areaTri(double* vt1, double* vt2, double* vt3);
void Allocate_geom(c_Tri_Element *&Geom_Panels);
void read_geometry_STL(c_Tri_Element *&Geom_Panels);

// function in wlkrClass.cpp
void diffusion_dis(double* a);

