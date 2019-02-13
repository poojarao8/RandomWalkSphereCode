#include <iostream>
#include <vector>
#include <math.h>
#include <random>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <sstream>
#include <string>
#include <fstream>
#include <algorithm>
#include <iterator>
#include "stl_serial.h"

using namespace std;

const double sphC[3] = {0.000125, 0.000125, 0.0004}; // center of the sphere
const double rad = 0.000080/3; // radius of the sphere
const double del_tw = 0.000005; // time step for the walker 
const double alpha_s = 0.279/(1200.*3200.); // thermal diffusivity

dim *prim_minpos;
dim *prim_maxpos;
dim glo_PanelMinPos;
dim glo_PanelMaxPos;

// Generate a random number between 0 and 1/
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
  
}


/* Source: https://stackoverflow.com/questions/17333/what-is-the-most-effective-way-for-float-and-double-comparison */

bool approximatelyEqual(float a, float b, float epsilon)
{
    return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

void vectorSubtract(double* v1, double* v2, double* vs)
{
  for(int i=0; i < 3; i++)
   {
     vs[i] = v1[i] - v2[i]; 
   }
}

void vectorScalarMul(double* v, double s, double* sv)
{
  for(int i=0; i < 3; i++)
   {
     sv[i] = s*v[i];
   }
}

double vectorNorm(double* v1)
{
  double vnorm;
  vnorm = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]) ;
  return vnorm;
}

double vectorDotProd(const double *x, const double *y, int n)
{
    double res = 0.0;
    int i;
    for (i = 0; i < n; i++)
    {
        res += x[i] * y[i];
    }
    return res;
}

void vectorCrossProd(double* v1, double* v2, double* vs)
{
  vs[0] = v1[1] * v2[2] - v1[2] * v2[1];
  vs[1] = v1[2] * v2[0] - v1[0] * v2[2];
  vs[2] = v1[0] * v2[1] - v1[1] * v2[0];

}

void matVectorMul(const double mat[][3], const double vec[], double* result, int rows, int cols)
{ // in matrix form: result = mat * vec;
    int i;
    for (i = 0; i < rows; i++)
    {
        result[i] = vectorDotProd(mat[i], vec, 3);
    }
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


double areaTri(double* vt1, double* vt2, double* vt3)
{
  double vt2_vt1[3]={vt2[0]-vt1[0], vt2[1]-vt1[1], vt2[2]-vt1[2]};
  double vt3_vt1[3]={vt3[0]-vt1[0], vt3[1]-vt1[1], vt3[2]-vt1[2]};
  double vcp[3];
  
  vcp[0] = vt2_vt1[1] * vt3_vt1[2] - vt2_vt1[2] * vt3_vt1[1];
  vcp[1] = vt2_vt1[2] * vt3_vt1[0] - vt2_vt1[0] * vt3_vt1[2];
  vcp[2] = vt2_vt1[0] * vt3_vt1[1] - vt2_vt1[1] * vt3_vt1[0];
  //cout << "cross prod" << vcp[0] << " "<< vcp[1] <<"" <<vcp[2]<<endl; 
  return vectorNorm(vcp);

}


void diffusion_dis(double* a)
{
   double del_xd[3]={0,0,0};

   double r1 = normalRand();
   double r2 = normalRand();
   double r3 = normalRand();

   del_xd[0] = sqrt (2.0*alpha_s*del_tw)*r1;
   del_xd[1] = sqrt (2.0*alpha_s*del_tw)*r2;
   //del_xd[2] = sqrt (2.0*alpha_s*del_tw)*r3;
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

bool vectorEqual(double* v1, double* v2)
{
  bool eql=false;
  if (approximatelyEqual(v1[0], v2[0], 1e-4) and approximatelyEqual(v1[1], v2[1], 1e-4) and approximatelyEqual(v1[2], v2[2], 1e-4))
    eql = true;

  return eql;
}

void c_Tri_Element::nearTris(map<int, vector<int>>& ntris)
{
  for (int k=0; k < glo_total_panels;  k++)
  { int count = 0;
   for (int kk=0; kk < glo_total_panels;  kk++)
    {
       if (vectorEqual(vert1[k].x, vert1[kk].x) or vectorEqual(vert1[k].x, vert2[kk].x) or vectorEqual(vert1[k].x, vert3[kk].x) or vectorEqual(vert2[k].x, vert1[kk].x) or vectorEqual(vert2[k].x, vert2[kk].x) or vectorEqual(vert2[k].x, vert3[kk].x) or vectorEqual(vert3[k].x, vert1[kk].x) or vectorEqual(vert3[k].x, vert2[kk].x) or vectorEqual(vert3[k].x, vert3[kk].x) )
      {  
         ntris[k].push_back(kk);
         //cout << ntris[k][count] << endl;
         //count++;

      }
  
    }

   }
     
}


void c_Tri_Element::calVolumeSphere(double* dVol)
{
 for (int k=0; k < glo_total_panels;  k++)
 {
 
  double e1[3];
  double e2[3];
  double e3[3];

  double se1[3];
  double se2[3];
  double se3[3];

  vectorSubtract(vert1[k].x, vert2[k].x, e1);
  vectorSubtract(vert2[k].x, vert3[k].x, e2);
  vectorSubtract(vert3[k].x, vert1[k].x, e3);

  vectorScalarMul(e1,1./3., se1);
  vectorScalarMul(e2,1./3., se2);
  vectorScalarMul(e3,1./3., se3);
 
  double a = vectorNorm(se1);
  double b = vectorNorm(se2);
  double c = vectorNorm(se3);

  // calculate the semi-perimeter
  double s = (a + b + c) / 2.0;

  // calculate the area
  double area = pow (s*(s-a)*(s-b)*(s-c),  0.5);
  dVol[k] = area*2e-6; // vol = area*thickness
  
 }
}


void c_Tri_Element::ComputeEdge1(int i_pan)
{
  /* Edge 1 is vector from vert1 to vert2 */
  for(int i_dim=0;i_dim<3;i_dim++)
  {
    edge1[i_pan].x[i_dim] = vert2[i_pan].x[i_dim] - vert1[i_pan].x[i_dim];
  }
}


void c_Tri_Element::ComputeEdge2(int i_pan)
{
  /* Edge 2 is vector from vert1 to vert3 */
  for(int i_dim=0;i_dim<3;i_dim++)
  {
    edge2[i_pan].x[i_dim] = vert3[i_pan].x[i_dim] - vert1[i_pan].x[i_dim];
  }
}

void sphereTrans(double* v1, double* vt)
{
  vt[0] = (v1[0] - sphC[0])/3.0;
  vt[1] = (v1[1] - sphC[1])/3.0;
  vt[2] = (v1[2] - sphC[2])/3.0; 

}

// added by PRAO
void c_Tri_Element::RandPtPick(int i_pan, double* d)
{
  float a;
  float b;
  
  a = unifRand();
  b = unifRand();
 
  //cout << "a, b " << a <<","<< b << endl;
  double v1_adj[3];
  double v2_adj[3];
  double v3_adj[3];

  sphereTrans(vert1[i_pan].x, v1_adj);
  sphereTrans(vert2[i_pan].x, v2_adj);
  sphereTrans(vert3[i_pan].x, v3_adj);
  
  sphereTrans(vert1[i_pan].x, v1_adj);
   
  d[0] = (1.0-sqrt(a))*v1_adj[0] + sqrt(a)*(1.0-b)*v2_adj[0] + b*sqrt(a)*v3_adj[0];
  d[1] = (1.0-sqrt(a))*v1_adj[1] + sqrt(a)*(1.0-b)*v2_adj[1] + b*sqrt(a)*v3_adj[1];
  d[2] = (1.0-sqrt(a))*v1_adj[2] + sqrt(a)*(1.0-b)*v2_adj[2] + b*sqrt(a)*v3_adj[2];

}


bool c_Tri_Element::isPtInTri (int i_pan, double* Pt, double& sum_weights)
{
    double area2, alpha, beta, gamma;//, sum_weights;

    double v1_adj[3];
    double v2_adj[3];
    double v3_adj[3];

    sphereTrans(vert1[i_pan].x, v1_adj);
    sphereTrans(vert2[i_pan].x, v2_adj);
    sphereTrans(vert3[i_pan].x, v3_adj);

    area2 = areaTri(v1_adj, v2_adj, v3_adj);
    alpha = areaTri(Pt, v2_adj, v3_adj) / area2 ;
    beta = areaTri(Pt, v3_adj, v1_adj) / area2 ;
    gamma = areaTri(Pt, v1_adj, v2_adj) / area2 ;

    //cout << Pt[0] << " " << Pt[1] << " " << Pt[2] << endl;

    sum_weights = alpha + beta + gamma; //this should be sent to the calling function
    //cout << "sum_weights " << sum_weights << " i=" << i_pan <<endl;


    if (alpha>=0.0 && beta>=0.0 && gamma>=0.0 && approximatelyEqual(sum_weights, 1.0, 0.2))
    {
      return true;
    }
    else
      return false;
    

}
//TODO: fix the indexLocate function
int c_Tri_Element::indexLocate(double* position, int i_pan, map<int, vector<int>>& ntris) 
{
  // check what panel does the point lie on
  bool ptLoc=false;
  double sum_wts=0.0;

  assert(ntris.find(i_pan) != ntris.end());
  vector<double> list_wts={}; //create a list of potential
  vector<int> list_potNbrs={};//create a list of potential nbrs

  for (auto iter = ntris[i_pan].begin(); iter != ntris[i_pan].end(); ++iter) 
   {
     //cout << *iter << endl;
     int i_nbr = *iter; // pass in neighbor's index

     ptLoc = isPtInTri(i_nbr, position, sum_wts);

     if (ptLoc==true){
      list_wts.push_back(sum_wts);
      list_potNbrs.push_back(i_nbr);
     }

   }

 if (list_wts.empty())
  {
   cout << "No Match Found for " << i_pan <<"!!!!"<<endl;
   cout << "Point " << position[0] << " " <<position[1] << " " << position[2] << endl; 
   cout << "Neighbors : " << endl;
  for (auto iter = ntris[i_pan].begin(); iter != ntris[i_pan].end(); ++iter)
     cout << *iter << endl;

    exit(0);
  }
 else
  {
   vector<double>::iterator result = min_element(std::begin(list_wts), std::end(list_wts));
   int i_fin = list_potNbrs[distance(std::begin(list_wts), result)];
   //cout << "J RESULT: " << i_fin << " " << endl;
   return i_fin;
  }


}

void c_Tri_Element::ComputeNormal(int i_pan)
{
  double x_e1 = edge1[i_pan].x[0];
  double y_e1 = edge1[i_pan].x[1];
  double z_e1 = edge1[i_pan].x[2];

  double x_e2 = edge2[i_pan].x[0];
  double y_e2 = edge2[i_pan].x[1];
  double z_e2 = edge2[i_pan].x[2];

  double nx = y_e1*z_e2 - z_e1*y_e2;
  double ny = -(x_e1*z_e2 - z_e1*x_e2);
  double nz = x_e1*y_e2 - y_e1*x_e2;

  double abs_n = sqrt(nx*nx + ny*ny + nz*nz);

  norm[i_pan].x[0] = nx/abs_n;
  norm[i_pan].x[1] = ny/abs_n;
  norm[i_pan].x[2] = nz/abs_n;

}

void c_Tri_Element::adjNumWlkrs(double TMaxCell, int i_maxPanel, double curTemp, double wTemp, vector<Walker> &listofwalkers)
{
  double dTemp = TMaxCell - curTemp; // curTemp is current temp of the max cell
  // TMaxCell is the max temperature that is maintained at fixed temperature
  //cout <<  "TMaxCell, wTemp" << TMaxCell << " " <<wTemp << endl;
  int nWlkrs = abs(dTemp/wTemp);
  int sWlkrs =  listofwalkers.size(); // no. of walkers currently
  // create walkers
  if (dTemp > 0.0)
  {
   //cout << "# walkers to be created = " << nWlkrs << endl;
   int tag_val = 1;
   for (int i=0; i <nWlkrs; ++i) // check!!!!
    {
     // initiate a new walker
     listofwalkers.push_back(Walker());
     // tag the walker
     listofwalkers[sWlkrs+i].tag = 1; //tag for max cell 
     // put in the current index
     listofwalkers[sWlkrs+i].posInd = i_maxPanel;
     // randomly pick a starting position in the panel to assign as a starting position for the walker
     //cout << "before pos " << listofwalkers[sWlkrs+i].curPos[0] << endl;
     RandPtPick(i_maxPanel, listofwalkers[sWlkrs+i].curPos);
     //cout << "pos " << listofwalkers[sWlkrs+i].curPos[0] << endl; 
    }
  }

  //kill extra walkers
  if (dTemp < 0.0 && nWlkrs > 0)
  {
    //cout << "# walkers to be killed in cell = " << nWlkrs << endl;
    int count = 0;
    for (auto it1 = begin (listofwalkers); it1 != end (listofwalkers); ++it1) 
     { 
        if (it1->posInd==i_maxPanel)
         {
            listofwalkers.erase( it1 );
            count++;
            if (count >= nWlkrs)
             break;
         }
     }
      
  }
 
  
}


void Walker::moveWalker(double* a)
{
  diffusion_dis(a);
}

void Allocate_geom(c_Tri_Element *&Geom_Panels)
{
   Geom_Panels = new c_Tri_Element[1];
   prim_minpos = new dim[glo_num_primitives];
   prim_maxpos = new dim[glo_num_primitives];

   for(int i = 0;i<1;i++)
   {
      Geom_Panels[i].vert1 = new dim[glo_total_panels];
      Geom_Panels[i].vert2 = new dim[glo_total_panels];
      Geom_Panels[i].vert3 = new dim[glo_total_panels];
      Geom_Panels[i].edge1 = new dim[glo_total_panels];
      Geom_Panels[i].edge2 = new dim[glo_total_panels];
      Geom_Panels[i].norm  = new dim[glo_total_panels];
   }

}

void read_geometry_STL(c_Tri_Element *&Geom_Panels)
{
  //cout << "welcome to read in geom " << endl;
  string *mystr = new string;
  ifstream stl_file;
  stl_file.open("AvcoatSingleSphere.stl");

  glo_PanelMinPos.x[0] = 500; glo_PanelMinPos.x[1] = 500; glo_PanelMinPos.x[2] = 500;
  glo_PanelMaxPos.x[0] = -1; glo_PanelMaxPos.x[1] = -1; glo_PanelMaxPos.x[2] = -1;

  int pan_id = 0;
  for(int i_prim=0;i_prim<glo_num_primitives;i_prim++)
  {
    getline(stl_file,*mystr); // SOLID_GEOM
    prim_minpos[i_prim].x[0] = 500; prim_minpos[i_prim].x[1] = 500 ; prim_minpos[i_prim].x[2] = 500;
    prim_maxpos[i_prim].x[0] = -1;  prim_maxpos[i_prim].x[1] = -1;   prim_maxpos[i_prim].x[2] = -1;

    for(int i_pan=0;i_pan<glo_panels_per_primitive;i_pan++)
    {
      char dummy1[8], dummy2[6];
      getline(stl_file,*mystr);
      stringstream ss_norm(*mystr);
      ss_norm >> dummy1 >> dummy2 >> Geom_Panels[0].norm[pan_id].x[0] >> Geom_Panels[0].norm[pan_id].x[1]
                                  >> Geom_Panels[0].norm[pan_id].x[2];
      getline(stl_file,*mystr); // outer loop

      getline(stl_file,*mystr); stringstream ss_v1(*mystr);
      ss_v1 >> dummy2 >> Geom_Panels[0].vert1[pan_id].x[0] >> Geom_Panels[0].vert1[pan_id].x[1] >> Geom_Panels[0].vert1[pan_id].x[2];

      getline(stl_file,*mystr); stringstream ss_v2(*mystr);
      ss_v2 >> dummy2 >> Geom_Panels[0].vert2[pan_id].x[0] >> Geom_Panels[0].vert2[pan_id].x[1] >> Geom_Panels[0].vert2[pan_id].x[2];

      getline(stl_file,*mystr); stringstream ss_v3(*mystr);
      ss_v3 >> dummy2 >> Geom_Panels[0].vert3[pan_id].x[0] >> Geom_Panels[0].vert3[pan_id].x[1] >> Geom_Panels[0].vert3[pan_id].x[2];

     // cout << Geom_Panels[0].vert3[pan_id].x[0] << endl;
    for(int i_dim=0;i_dim<3;i_dim++)
    {
      if(Geom_Panels[0].vert1[pan_id].x[i_dim] < glo_PanelMinPos.x[i_dim])
      {
        glo_PanelMinPos.x[i_dim] = Geom_Panels[0].vert1[pan_id].x[i_dim];
      }
      if(Geom_Panels[0].vert2[pan_id].x[i_dim] < glo_PanelMinPos.x[i_dim])
      {
        glo_PanelMinPos.x[i_dim] = Geom_Panels[0].vert2[pan_id].x[i_dim];
      }
      if(Geom_Panels[0].vert3[pan_id].x[i_dim] < glo_PanelMinPos.x[i_dim])
      {
        glo_PanelMinPos.x[i_dim] = Geom_Panels[0].vert3[pan_id].x[i_dim];
      }


      if(Geom_Panels[0].vert1[pan_id].x[i_dim] > glo_PanelMaxPos.x[i_dim])
      {
        glo_PanelMaxPos.x[i_dim] = Geom_Panels[0].vert1[pan_id].x[i_dim];
      }
      if(Geom_Panels[0].vert2[pan_id].x[i_dim] > glo_PanelMaxPos.x[i_dim])
      {
        glo_PanelMaxPos.x[i_dim] = Geom_Panels[0].vert2[pan_id].x[i_dim];
      }
      if(Geom_Panels[0].vert3[pan_id].x[i_dim] > glo_PanelMaxPos.x[i_dim])
      {
        glo_PanelMaxPos.x[i_dim] = Geom_Panels[0].vert3[pan_id].x[i_dim];
      }

    }

      getline(stl_file,*mystr); // end loop
      getline(stl_file,*mystr); // end facet
      pan_id++ ;
    }
    getline(stl_file,*mystr); // end GEOM
  }

  stl_file.close();
}

void c_Tri_Element::writeInSTL(string filename,  int glo_num_subparts, int*&glo_panels_per_subpart, int*& primitiveSubpartRemove, string*& subpartStringPrim, double *temp){
   ofstream stl_file;

   string filenameSTL = filename+".stl";
   stl_file.open(filenameSTL.c_str());
   cout << "Writing stl geometry to file: " << filenameSTL.c_str() << endl;
   int startPanID = 0;
   for(int i_sub=0; i_sub<glo_num_subparts; i_sub++)
   {
    if( primitiveSubpartRemove[i_sub] == 0){
    stl_file << subpartStringPrim[i_sub] << endl;

    for(int i_pan=startPanID; i_pan<(startPanID+glo_panels_per_subpart[i_sub]);i_pan++)
    {

      stl_file << "  facet normal " << norm[i_pan].x[0] <<" "<< norm[i_pan].x[1] <<" "<< norm[i_pan].x[2] << endl;
      stl_file << "    outer loop"<< endl;
      stl_file << "     vertex "<< vert1[i_pan].x[0] <<" "<<vert1[i_pan].x[1] <<" "<< vert1[i_pan].x[2] << endl;
      stl_file << "     vertex "<< vert2[i_pan].x[0] <<" "<<vert2[i_pan].x[1] <<" "<< vert2[i_pan].x[2] << endl;
      stl_file << "     vertex "<< vert3[i_pan].x[0] <<" "<<vert3[i_pan].x[1] <<" "<< vert3[i_pan].x[2] << endl;
      stl_file << "    endloop"<< endl;
      stl_file << "  endfacet"<< endl;

    }
    stl_file << "end"<<  subpartStringPrim[i_sub] << endl;
    }
    startPanID += glo_panels_per_subpart[i_sub];
   } //subpart
  stl_file.close();

   ofstream tempFile;
   string filenameTemp = filename+"_Temp.dat";
   tempFile.open(filenameTemp.c_str());
   cout << "Writing stl geometry to file: " << filenameTemp.c_str() << endl;
   startPanID = 0;
   for(int i_sub=0; i_sub<glo_num_subparts; i_sub++)
   {
    if( primitiveSubpartRemove[i_sub] == 0){
     for(int i_pan=startPanID; i_pan<(startPanID+glo_panels_per_subpart[i_sub]);i_pan++)
     {
        tempFile << temp[i_pan] << endl;
//        cout << "temp writing for pan: " << i_pan << setw(15) << temp[i_pan] << endl;
     }
    }

   startPanID += glo_panels_per_subpart[i_sub];
   } //subpart
   tempFile.close();
}

/*
void WriteGeomGrid(c_Tri_Element *&Geom_Panels, double** timeToDegrade, int** remove )
{
   int num_vert_points = glo_total_panels*3;
   char char_panels[256] ; sprintf(char_panels,"%08d",glo_total_panels);
   char char_vert_points[256] ; sprintf(char_vert_points,"%08d",num_vert_points);
   string str_numvert(char_vert_points);
   string str_numpanels(char_panels);
   ofstream IB_grid;
   string IB_grid_filename="Output/geometryOutput.tec";
   IB_grid.open(IB_grid_filename.c_str());
   IB_grid << "title = \" Immbersed Body Grid \" " << endl;
   IB_grid << "variables = \"x\", \"y\", \"z\", \"heatFlux [W/cm^2]\", \"timeToDegrade [microSeconds]\", \"remove\",\"temperature\",\"nearbyGasTemp\"" << endl;
//   IB_grid << "variables = \"x\", \"y\", \"z\",\"heatFlux\"" << endl;
   IB_grid << "zone n="<<str_numvert<<", e="<<str_numpanels;
   IB_grid << ", DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE" << endl;
   IB_grid << " VARLOCATION = ([4-8] = CELLCENTERED)" << '\n';
   IB_grid << " "<< endl;

   for(int i_prim=0;i_prim<glo_num_primitives;i_prim++)   {
     for(int i_pan=0;i_pan<glo_panels_per_primitive[i_prim];i_pan++)   {

     IB_grid << setprecision(16);
     IB_grid << scientific ;
     IB_grid << setw(20) \
             << Geom_Panels[i_prim].vert1[i_pan].x[0] << "\t" \
             << Geom_Panels[i_prim].vert2[i_pan].x[0] << "\t" \
             << Geom_Panels[i_prim].vert3[i_pan].x[0] << endl;
     }
   }
   IB_grid << endl;
   for(int i_prim=0;i_prim<glo_num_primitives;i_prim++)   {
     for(int i_pan=0;i_pan<glo_panels_per_primitive[i_prim];i_pan++)   {

     IB_grid << setprecision(16);
     IB_grid << scientific ;
     IB_grid << setw(20) \
             << Geom_Panels[i_prim].vert1[i_pan].x[1] << "\t" \
             << Geom_Panels[i_prim].vert2[i_pan].x[1] << "\t" \
             << Geom_Panels[i_prim].vert3[i_pan].x[1] << endl;
     }
   }
   IB_grid << endl;
   for(int i_prim=0;i_prim<glo_num_primitives;i_prim++)   {
     for(int i_pan=0;i_pan<glo_panels_per_primitive[i_prim];i_pan++)   {
        IB_grid << setprecision(16);
        IB_grid << scientific ;
        IB_grid << setw(20) \
                << Geom_Panels[i_prim].vert1[i_pan].x[2] << "\t" \
                << Geom_Panels[i_prim].vert2[i_pan].x[2] << "\t" \
                << Geom_Panels[i_prim].vert3[i_pan].x[2] << endl;
     }
   }
  IB_grid << endl;
   for(int i_prim=0;i_prim<glo_num_primitives;i_prim++)   {
    for(int i_pan=0;i_pan<glo_panels_per_primitive[i_prim];i_pan++)   {
        IB_grid << setprecision(15);
        IB_grid << scientific ;
        IB_grid << setw(20) \
                << (Geom_Panels[i_prim].heatFlux[i_pan]/1e+04) << "\t" << endl;
    }
   }
   IB_grid << endl;
   for(int i_prim=0;i_prim<glo_num_primitives;i_prim++)   {
    for(int i_pan=0;i_pan<glo_panels_per_primitive[i_prim];i_pan++)   {
        IB_grid << setprecision(15);
        IB_grid << scientific ;
        IB_grid << setw(20) \
                << timeToDegrade[i_prim][i_pan] << "\t" << endl;
    }
   }
   IB_grid << endl;
   for(int i_prim=0;i_prim<glo_num_primitives;i_prim++)   {
    for(int i_pan=0;i_pan<glo_panels_per_primitive[i_prim];i_pan++)   {
        IB_grid << setprecision(15);
        IB_grid << scientific ;
        IB_grid << setw(20) \
                << remove[i_prim][i_pan] << "\t" << endl;
    }
   }
   IB_grid << endl;
   for(int i_prim=0;i_prim<glo_num_primitives;i_prim++)   {
    for(int i_pan=0;i_pan<glo_panels_per_primitive[i_prim];i_pan++)   {
        IB_grid << setprecision(15);
        IB_grid << scientific ;
        IB_grid << setw(20) \
                << Geom_Panels[i_prim].temp[i_pan] << "\t" << endl;
    }
   }
   IB_grid << endl;
   for(int i_prim=0;i_prim<glo_num_primitives;i_prim++)   {
    for(int i_pan=0;i_pan<glo_panels_per_primitive[i_prim];i_pan++)   {
        IB_grid << setprecision(15);
        IB_grid << scientific ;
        IB_grid << setw(20) \
                << Geom_Panels[i_prim].nearbyGasTemp[i_pan] << "\t" << endl;
    }
   }
   int node_sum=0;
   for(int i_prim=0;i_prim<glo_num_primitives;i_prim++)   {
    for(int i_pan=0;i_pan<glo_panels_per_primitive[i_prim];i_pan++)   {
     for(int i_vert=0;i_vert<3;i_vert++)
     {
       node_sum += 1 ;
       IB_grid << node_sum << "\t";
     }
     IB_grid << endl;
    }
   }
   IB_grid.close();
   cout <<"Geometry output is written! " << endl;
}
*/


int fast_mod(const int input, const int ceil) {
    // apply the modulo operator only when needed
    // (i.e. when the input is greater than the ceiling)
    return input >= ceil ? input % ceil : input;
    // NB: the assumption here is that the numbers are positive
}
