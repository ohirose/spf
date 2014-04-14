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
#include<assert.h>
#include<math.h>
#include"util.h"
#include"dpmeans.h"

int dpmeans(
      int           *  z,    /* OUTPUT: cluster id | size N           */ 
      double        ** m,    /* OUTPUT: means      | size K(N) x P    */
      double        ** d,    /* OUTPUT: distance^2 | size N    x K(N) */
      int           *  l,    /* OUTPUT: #members   | size K(N)        */
      int           *  K,    /* OUTPUT: #clusters  | size 1           */ 
      const double  ** X,    /* INPUT : data       | size N   x P     */ 
      const int        N,    /* INPUT : #samples   | size 1           */ 
      const int        P,    /* INPUT : #variables | size 1           */ 
      const double     lmd,  /* INPUT : parameter  | size 1           */ 
      const double     dz,   /* INPUT : zscale     | size 1           */
      const int        nlp   /* INPUT : #loops     | size 1           */
  ){

  int    n,p,k,i;
  int    kmin,nmin,Kt=1;
  double val,min=1.0E200;
  //const double lmd2 = lmd*lmd;

  assert(P==2||P==3);

  /* Initialization */
  for(n=0;n<N;n++){z[n]=0;l[n]=0;}
  for(k=0;k<N;k++){for(n=0;n<N;n++)d[n][k]=0;for(p=0;p<P;p++)m[k][p]=0;}
  for(p=0;p<P;p++){for(n=0;n<N;n++)m[0][p]+=X[n][p];m[0][p]/=(double)N;} 
  for(n=0;n<N;n++){val=wdist(X[n],m[0],P,dz);if(val<min){min=val;nmin=n;}}    /* Finding the nearest center */ 
  for(p=0;p<P;p++)m[0][p]=X[nmin][p];
  
  /* Main computation */
  for(i=0;i<nlp;i++){
    for(n=0;n<N;n++){min=1.0E200;kmin=-1;
      for(k=0;k<Kt;k++)d[n][k]=wdist(X[n],m[k],P,dz);                         /* Calcuating distance        */
      for(k=0;k<Kt;k++)if(d[n][k]<min){min=d[n][k];kmin=k;}                   /* Finding the nearest center */ 
      if(min>lmd){z[n]=Kt;for(p=0;p<P;p++)m[Kt][p]=X[n][p];Kt++;} else z[n]=kmin;  
    }
    for(n=0;n<N ;n++)l[n]=0;for(n=0;n<N ;n++)l[z[n]]++;                       /* Counting cluster members   */

    /* Update means */
    for(k=0;k<Kt;k++)for(p=0;p<P;p++)m[k][p]=0;
    for(n=0;n<N ;n++)for(p=0;p<P;p++)m[z[n]][p]+=X[n][p];
    for(k=0;k<Kt;k++)for(p=0;p<P;p++){
      assert(l[k]>0);m[k][p]/=(double)l[k];
    }
  } 
  *K=Kt;
  
  return 0;
} 

  
double wdist(const double * x1, const double * x2, const int P, const double dz){ 
  int p; double v,d=0.0; 
  for(p=0;p<P;p++){v=(x1[p]-x2[p])*(x1[p]-x2[p]); d+=(p<P-1)?v:dz*dz*v;}
  return sqrt(d);
} 

int cutoff  (byte *y,const int *imsize, const byte cut){ 
  int l,L=imsize[0]*imsize[1]*imsize[2]; for(l=0;l<L;l++)if(y[l]<cut)y[l]=0;
  return 0;
}

int localmax(double **vx, int *nvx, const byte *y, const int *imsize, const int *width){
  int  i,j,k,l,l0,I,J,K,L,a,b,c,A,B,C,n=0;

  I=imsize[0];J=imsize[1];K=imsize[2];L=I*J*K; A=width[0];B=width[1];C=width[2]; *nvx=0;

  for(i=0;i<I;i++)for(j=0;j<J;j++)for(k=0;k<K;k++){l0=i+j*I+k*I*J;if(y[l0]==0) goto tag;
    for(a=-A;a<=A;a++)for(b=-B;b<=B;b++)for(c=-C;c<=C;c++){l=i+a+(j+b)*I+(k+c)*I*J;
         if((i+a>=0)&&(i+a<I)&&(j+b>=0)&&(j+b<J)&&(k+c>=0)&&(k+c<K)&&l!=l0)if(y[l0]<y[l])goto tag;
    }
    vx[n][0]=(double)i;vx[n][1]=(double)j;vx[n][2]=(double)k;n++; tag:;
  }*nvx=n;
  
  return 0;
} 

