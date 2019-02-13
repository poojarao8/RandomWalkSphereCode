#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <math.h>
#include <vector>
#include <iterator>
#include <algorithm>
#include "main.h"

using namespace std;

int glo_num_primitives=1;
int glo_panels_per_primitive=140;
int glo_total_panels=140;

double dir_z[3] ={0,0,-1}; //needed for initial temperature distribution

const double sphC[3] = {0.000125, 0.000125, 0.0004}; // center of the sphere
const double rad = 0.000080/3; // radius of the sphere
const double del_tw = 0.000005; // time step for the walker 
const double alpha_s = 0.279/(1200.*3200.); // thermal diffusivity

dim *prim_minpos;
dim *prim_maxpos;
dim glo_PanelMinPos;
dim glo_PanelMaxPos;

void sphereTrans(double* v1, double* vt)
{
  vt[0] = (v1[0] - sphC[0])/3.0;
  vt[1] = (v1[1] - sphC[1])/3.0;
  vt[2] = (v1[2] - sphC[2])/3.0;
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

void c_Tri_Element::tempDistribution(double rTemp, int worldSize, double* Temp)
{
  for (int i=0; i < glo_total_panels ; ++i)
   {
    //Temp[i] = 500.*vectorDotProd(dir_z, Geom_Panels[0].norm[i].x, 3.0);
    
    Temp[i] = 500.*vectorDotProd(dir_z, norm[i].x, 3.0);

    if (Temp[i] > 0.0)
      Temp[i] = rTemp + Temp[i];
    else if (Temp[i] < 0.0)
      Temp[i] = rTemp;
    //cout << temp[i] << endl; 

    //PARALLEL: Instead of maintaining total temp, i will maintain the perturbation temp per processor
    Temp[i] = (Temp[i] - rTemp)/worldSize;
   }

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

int c_Tri_Element::adjNumWlkrs(double TMaxCell, int i_maxPanel, double curTemp, double wTemp, vector<Walker> &listofwalkers, int numProcs, int myRank)
{
  double dTemp = TMaxCell - curTemp; // curTemp is current temp of the max cell
  // TMaxCell is the max temperature that is maintained at fixed temperature
  //cout <<  "TMaxCell, wTemp" << TMaxCell << " " <<wTemp << endl;
  int nWlkrs = abs(dTemp/wTemp);
  int sgn;

  //PARALLEL: We need to adjust the number of walkers created/destroyed here to account for parallel
  //          To do this we take the floor (accomplished through integer division)
  //          of the total numWalkers/numProcs
  int nWlkrsFloor = nWlkrs/numProcs;
 
  //PARALLEL: Here we are taking the remainder from the integer division and adding those walkers on
  //          the lowest processors, so we get the correct number of created or removed walkers
  if (myRank < (nWlkrs % numProcs))
     nWlkrsFloor = nWlkrsFloor + 1;

  nWlkrs = nWlkrsFloor;

  int sWlkrs =  listofwalkers.size(); // no. of walkers currently
  // create walkers
  if (dTemp > 0.0)
  {
   sgn = 1;
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
    sgn = -1;
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
 
  return nWlkrs*sgn;
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
  stl_file.open("/home/prao/rw_par++/parallel/input/AvcoatSingleSphere.stl");

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

