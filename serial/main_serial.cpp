#include <iostream>
#include <vector>
#include <math.h>
#include <random>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <sstream>
#include <string>
#include "stl_serial.h"

using namespace std;

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const string currentDateTime() 
{
  time_t     now = time(0);
  struct tm  tstruct;
  char       buf[80];
  tstruct = *localtime(&now);
  // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
  // for more information about date/time format
  strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

  return buf;
}

// global variables

int glo_num_primitives=1;
int glo_panels_per_primitive=140;
int glo_total_panels=140;

//double temp[glo_total_panels];
double temp[140];
double dir_z[3] ={0,0,-1}; 

// simulation parameters
double RTemp = 300.; // reference temp
double nWlkrs = 3000.; // number of walkers
double rho_s = 1200.0; // density of the solid
double Cp_s = 3200.0; // specific heat of the solid
double tot_delH = 8.25746164447e-06; //TODO: add full calculation
double dH_wlkr = tot_delH/nWlkrs; // estimate of enthalpy per walker
int sizeOfWlkrs;

int main(int argc, char** argv) 
{
  // Load the stl file
  c_Tri_Element *Geom_Panels;
  Geom_Panels = new c_Tri_Element[1];
  Allocate_geom(Geom_Panels); 
  read_geometry_STL(Geom_Panels);

  // Create the temperature distribution
  for (int i=0; i < glo_total_panels ; ++i)
   {
    temp[i] = 500.*vectorDotProd(dir_z, Geom_Panels[0].norm[i].x, 3.0);
    
    if (temp[i] > 0.0)
      temp[i] = 300. + temp[i];
    else if (temp[i] < 0.0)
      temp[i] = 300;
    //cout << temp[i] << endl; 
   } 

  double dV[140]={};
  double WTemp[140]={};
  Geom_Panels[0].calVolumeSphere(dV);

  //cout <<"dH_wlkr " <<  dH_wlkr << endl;
  // Temperature of walker in each panel
  for (int j=0; j < glo_total_panels;  j++)
  {  
    WTemp[j] = dH_wlkr/(rho_s*Cp_s*dV[j]);
    //cout << "walker temp " << WTemp[j] << endl;
  }
 
  // Store the neighbouring triangles of each tri in a map
   map<int, vector<int>> nbrTris;
   Geom_Panels[0].nearTris(nbrTris);
   /*for (int jjj=0; jjj<14; jjj++)
   {
    cout << nbrTris[38][jjj] << endl;
   }*/
   
  // Create walkers
  vector<Walker> listofwalkers;
  seed(); //seed for generating random numbers

  int cnt_wlkrs = 0;
  for (int i = 0; i < 140; i++)
  {
     double dTemp = temp[i] - RTemp;
     int nW = dTemp/WTemp[i];
     //cout << "at panel=" << i << "# walkers" << nW<<endl;
          
     if (nW > 0)
     {
       for (int ii=0; ii<nW; ii++)
       {
         listofwalkers.push_back(Walker());
         if (i==41)//index for the hottest cell
           listofwalkers[cnt_wlkrs].tag = 1;
         else listofwalkers[cnt_wlkrs].tag = 0;

         listofwalkers[cnt_wlkrs].posInd = i;
         Geom_Panels[0].RandPtPick(i, listofwalkers[cnt_wlkrs].curPos);
         //cout << listofwalkers[cnt_wlkrs].curPos[0] << " " << listofwalkers[cnt_wlkrs].curPos[1] << " " << listofwalkers[cnt_wlkrs].curPos[2]<<endl; 
         cnt_wlkrs++;
         
       }
      }
    }  
  
  cout << "total walkers = " << cnt_wlkrs << endl;

  vector<Walker>::iterator it;

   // main time loop
   //Move walkers and update the temperature
   for (int JJ=0; JJ < 5000; JJ++)
   {
     int count = 0;
     // Loop over all the walkers and move them
     //cout << "timestep : " << JJ << endl;
    // cout << "temp at start of timestep  " << temp[41] << endl;
    //if (JJ%100 == 0)
    if(fast_mod(JJ, 100)==0) 
    {
     //cout << JJ << endl;
     cout << JJ*del_tw<<" " << temp[20]<< " # walkers " << listofwalkers.size() <<  " time: " << currentDateTime() << endl;
    } 
     for (it = listofwalkers.begin(); it != listofwalkers.end(); ++it)      
     {
        //cout << "initial position" << it->curPos[0] << " "<<it->curPos[1] << " "<<it->curPos[2] << endl;
        //cout << "initial index" << it->posInd << endl;
        it->moveWalker(it->curPos);
        //cout << "final position" << it->curPos[0] <<" " << it->curPos[1] << " "<<it->curPos[2] << endl;
        //cout << it->curPos[0] << endl;
        int inPos = it->posInd;
        int finPos = Geom_Panels[0].indexLocate(it->curPos, it->posInd, nbrTris);
        it->posInd = finPos;
        //cout << "final index" << finPos << endl;
        // Update tempearture in each panel
        //cout << "finPos " << finPos << endl; 
        temp[inPos] = temp[inPos] - WTemp[inPos];
        temp[finPos] = temp[finPos] + WTemp[finPos];
        if (finPos==41)
         {
           it->tag = 1;
         }
        else it->tag = 0;

  
      } 
      Geom_Panels[0].adjNumWlkrs(785.6, 41, temp[41], WTemp[41], listofwalkers);    
      temp[41] = 785.6; 
      sizeOfWlkrs = listofwalkers.size();
      //cout << " size of listofwalkers" << sizeOfWlkrs<< endl; 
      //cout << "end of timestep " << JJ << endl;
      //cout << "temp at end of timestep " << temp[20] << endl;
      
    }

  return 0;

}
