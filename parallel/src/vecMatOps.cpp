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

// function prototyping for MatrixInversion()
int GetMinor(double src[][3], double dest[][3], int row, int col, int order);
double CalcDeterminant( double mat[][3], int order);


bool approximatelyEqual(float a, float b, float epsilon)
{
    return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool vectorEqual(double* v1, double* v2)
{
  bool eql=false;
  if (approximatelyEqual(v1[0], v2[0], 1e-4) and approximatelyEqual(v1[1], v2[1], 1e-4) and approximatelyEqual(v1[2], v2[2], 1e-4))
    eql = true;

  return eql;
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


void MatrixInversion(double A[][3], int order, double Y[][3])
{
    // get the determinant of a
    double det = 1.0/CalcDeterminant(A,order);

    // memory allocation
    double temp[2][2]={};
    double minor[3][3]={};

    for(int j=0;j<order;j++)
    {
        for(int i=0;i<order;i++)
        {
            // get the co-factor (matrix) of A(j,i)
            GetMinor(A,minor,j,i,order);
            Y[i][j] = det*CalcDeterminant(minor,order-1);
            if( (i+j)%2 == 1)
                Y[i][j] = -Y[i][j];
        }
    }

    // release memory
    //delete [] minor[0];
}


// calculate the cofactor of element (row,col)
int GetMinor(double src[][3], double dest[][3], int row, int col, int order)
{
    // indicate which col and row is being copied to dest
    int colCount=0,rowCount=0;

    for(int i = 0; i < order; i++ )
    {
        if( i != row )
        {
            colCount = 0;
            for(int j = 0; j < order; j++ )
            {
                // when j is not the element
                if( j != col )
                {
                    dest[rowCount][colCount] = src[i][j];
                    colCount++;
                }
            }
            rowCount++;
        }
    }

    return 1;
}

// Calculate the determinant recursively.
double CalcDeterminant( double mat[][3], int order)
{
    // order must be >= 0
    // stop the recursion when matrix is a single element
    if( order == 1 )
        return mat[0][0];

    // the determinant value
    double det = 0;

    // allocate the cofactor matrix
    double minor[3][3];

    for(int i = 0; i < order; i++ )
    {
        // get minor of element (0,i)
        GetMinor( mat, minor, 0, i , order);
        // the recusion is here!

        det += (i%2==1?-1.0:1.0) * mat[0][i] * CalcDeterminant(minor,order-1);
        //det += pow( -1.0, i ) * mat[0][i] * CalcDeterminant( minor,order-1 );
    }

    // release memory

    return det;

}


