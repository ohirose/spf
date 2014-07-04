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

#include"util.h"
#include"tree.h"
#include"dpmeans.h"
#include"io.h"

#define NVXLIMIT 2048

int    mappsf       (char *y, const int *imsize, const double **X, const int V, const int *objsize, const double sgm);
int    normal       (double *v, int nv);

/* NOTE --------------------------------------------------------------*/               
/* Use 'Mersenne Twister (mt19937ar.c)' as a rundom number generator. */
/* Source code is available in the author's webpage:                  */
/* http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html           */
/* The following two functions are defined in the 'mt19937ar.c' file. */
/* -------------------------------------------------------------------*/               
void          init_genrand   (unsigned long s);
double        genrand_real1  (void);


int main (int argc, char** argv){

  FILE   *fp; const int P=3;

  // Parameters
  double sgmt[3],sgms[3];
  int    objsize[3];
  double alpha,beta,zscale;
  int    seed;
  char   init[256],data[256],traj[256];

  int    imsize[4];
  int    i,t,k,p,T,K,k1,k2,u,L;
  int    root,rad;
  int    cut=0;

  char   *y;
  double *eta,*eta0,*mov;
  double **D,***noise; 
  double ***X,*buf;
  int    **G,*U,*V,*Ks; 

  fp=fopen("conf-sim.txt", "r");if(!fp){printf("File: \'param.txt\' Not Found.\n"); exit(1);}
  fscanf(fp,"imsize:%d,%d,%d,%d\n", imsize,imsize+1,imsize+2,imsize+3);
  fscanf(fp,"K:%d\n",               &K);
  fscanf(fp,"root:%d\n",            &root);
  fscanf(fp,"seed:%d\n",            &seed);
  fscanf(fp,"alpha:%lf\n",          &alpha);
  fscanf(fp,"beta:%lf\n",           &beta);
  fscanf(fp,"sgmt:%lf,%lf,%lf\n",   sgmt,sgmt+1,sgmt+2);
  fscanf(fp,"sgms:%lf,%lf,%lf\n",   sgms,sgms+1,sgms+2);
  fscanf(fp,"objsize:%d,%d,%d\n",   objsize,objsize+1,objsize+2);
  fscanf(fp,"init:%s\n",            init);
  fscanf(fp,"data:%s\n",            data);
  fscanf(fp,"traj:%s\n",            traj);
  fclose(fp);fp=NULL; init_genrand(seed);

  rad=objsize[0]/2; zscale=objsize[0]/objsize[2];
  L=imsize[0]*imsize[1]*imsize[2];T=imsize[3];
  fclose(fp);fp=NULL;

  X     = calloc3d (T,K,P+1);
  noise = calloc3d (T,K,P);
  D     = calloc2d (K,K);
  G     = calloc2i (K,K);

  buf   = calloc   (T*K*P,sizeof(double));   
  y     = calloc   (L,    sizeof(char));   
  eta0  = calloc   (P,    sizeof(double));
  eta   = calloc   (P,    sizeof(double));
  mov   = calloc   (P,    sizeof(double));
  V     = calloc   (K,    sizeof(int));     
  U     = calloc   (K,    sizeof(int));    
  Ks    = calloc   (T,    sizeof(int)); for(t=0;t<T;t++)Ks[t]=K;  

  fp = fopen(init,"r");if(!fp){printf("File: \'%s\' Not Found.\n",init);exit(1);}
  for(k=0;k<K;k++) fscanf(fp,"%lf\t%lf\t%lf\t%lf\n",X[0][k],X[0][k]+1,X[0][k]+2,X[0][k]+3);
  fclose(fp);fp=NULL;

  for(k1=0;k1<K;k1++)for(k2=0;k2<K;k2++) D[k1][k2]=wdist(X[0][k1],X[0][k2],P,zscale);
  mstree(G,(const double**)D,K); 
  makeiter(V,U,(const int**)G,K,root);

  normal(buf,T*K*P);
  for(t=0;t<T;t++)for(i=0;i<K;i++)for(p=0;p<P;p++)
    noise[t][V[i]][p]=buf[p+k*P+t*P*K]*((V[i]==root)?sgmt[p]:sgms[p]); 

  fp = fopen(data,"wb"); fwrite(imsize,sizeof(int),4,fp);
  for(t=0;t<T;t++){
    if(t)for(i=0;i<K;i++){k=V[i];u=U[k];assert(k>=0&&k<K&&u>=-1&&u<K); 
      if(i)for(p=0;p<P;p++){
        eta [p]=X[t-1][k][p]-X[t-1][u][p];
        eta0[p]=X  [0][k][p]-X  [0][u][p];
        mov [p]=alpha*eta[p]+(1-alpha)*eta0[p];
        X   [t][k][p]=X[t][u][p]+mov[p]+noise[t][k][p];
      }
      else for(p=0;p<P;p++) X[t][k][p]=0.97*X[t-1][k][p]+0.03*X[0][k][p]+noise[t][k][p]; 
    }
    mappsf(y,imsize,(const double**)X[t],K,objsize,2.0);
    fwrite(y,1,L,fp);
  } 
  fclose(fp); 
  write(traj,(const double**)X,Ks,imsize,objsize,cut);  

  return 0;
}


int normal(double *v, int nv){
  int i,n=(nv/2)-1,odd=(nv%2==1)?1:0; double u1,u2;
  if(n)for(i=0;i<n;i++){u1=genrand_real1();u2=genrand_real1();
    v[2*i  ]=sqrt(-2*log(u1))*cos(2*M_PI*u2);
    v[2*i+1]=sqrt(-2*log(u1))*sin(2*M_PI*u2);
  } 
  if(odd){u1=genrand_real1();u2=genrand_real1();v[nv-1]=sqrt(-2*log(u1))*cos(2*M_PI*u2);}
  return 0;
}

int mappsf(char *y, const int *imsize, const double **X, const int V, const int *objsize, const double sgm){
  int a,b,c,l,v,A,B,C,I,J,K,L,P=3;
  double x[3]; 
  double zscale = objsize[0]/(double)objsize[2];

  I=imsize[0];J=imsize[1];K=imsize[2];L=I*J*K; A=objsize[0]/2;B=objsize[1]/2;C=objsize[2]/2;
  for(l=0;l<L;l++)y[l]=0;

  for(v=0;v<V;v++){
    for(a=-A;a<=A;a++)for(b=-B;b<=B;b++)for(c=-C;c<=C;c++){
      x[0]=floor(X[v][0])+a;x[1]=floor(X[v][1])+b;x[2]=floor(X[v][2])+c;l=x[0]+x[1]*I+x[2]*I*J; 
      if(x[0]>=0&&x[0]<I&&x[1]>=0&&x[1]<J&&x[2]>=0&&x[2]<K)
        y[l]=255*exp(-pow(wdist(X[v],x,P,zscale),2)/(sgm*sgm));
    }
  }

  return 0;
}


