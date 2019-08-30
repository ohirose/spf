// Copyright (c) 2014 Osamu Hirose
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>

#define SQ(x) ((x)*(x))
#define CB(x) ((x)*(x)*(x))
#define SS(x,y,z) (((x)*(x))+((y)*(y))+((z)*(z)))

double dist2(const double x1[3], const double x2[3], const int D){
  assert(D==2||D==3);
  return D==3 ? SQ(x1[0]-x2[0])+SQ(x1[1]-x2[1])+SQ(x1[2]-x2[2])
              : SQ(x1[0]-x2[0])+SQ(x1[1]-x2[1]);
}

double quad(const double u[3], const double A[3][3]){
  int i,j;double val=0,v[3]={0,0,0};
  for(i=0;i<3;i++)for(j=0;j<3;j++) v[i]+=A[i][j]*u[j];
  for(i=0;i<3;i++)                 val +=u[i]*v[i];
  return val;
}

double det3(const double A[3][3]){ return
   A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])
  -A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
  +A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
}

void eig3(double lmd[3], const double A[3][3]){
  double p,q,v,B[3][3]; double c=(2.0*M_PI)/3.0;
  q=(A[0][0]+A[1][1]+A[2][2])/3.0;

  B[0][0]=A[0][0]-q; B[0][1]=B[1][0]=A[0][1];
  B[1][1]=A[1][1]-q; B[0][2]=B[2][0]=A[0][2];
  B[2][2]=A[2][2]-q; B[1][2]=B[2][1]=A[1][2];

  p=sqrt((SS(B[0][0],B[1][1],B[2][2])+2*SS(B[0][1],B[0][2],B[1][2]))/6.0); 
  v=acos(det3((const double (*)[3])B)/(2.0*CB(p)))/3.0;

  lmd[0]=q+2*cos(v    )*p;
  lmd[1]=q+2*cos(v+2*c)*p;
  lmd[2]=q+2*cos(v+  c)*p;
} 

int inv3(double iA[3][3], const double A[3][3]){ 
  double det=det3(A); if(fabs(det)<1e-10) goto error;

  iA[0][0]=(A[1][1]*A[2][2]-A[1][2]*A[2][1])/det;
  iA[0][1]=(A[0][2]*A[2][1]-A[0][1]*A[2][2])/det;
  iA[0][2]=(A[0][1]*A[1][2]-A[0][2]*A[1][1])/det;

  iA[1][0]=(A[1][2]*A[2][0]-A[1][0]*A[2][2])/det;
  iA[1][1]=(A[0][0]*A[2][2]-A[0][2]*A[2][0])/det;
  iA[1][2]=(A[0][2]*A[1][0]-A[0][0]*A[1][2])/det;

  iA[2][0]=(A[1][0]*A[2][1]-A[1][1]*A[2][0])/det;
  iA[2][1]=(A[0][1]*A[2][0]-A[0][0]*A[2][1])/det;
  iA[2][2]=(A[0][0]*A[1][1]-A[0][1]*A[1][0])/det;
  return 0;

  error:
    printf("Inverting a singular matrix. Abort.\n");
    exit  (EXIT_FAILURE);
}

