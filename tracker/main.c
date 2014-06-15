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
#include"dpmeans.h"
#include"tree.h"
#include"io.h"

#define NVXLIMIT 2048

int    findroot     (int *root, const byte *y, 
                     const double **mx, const int K, const int *imsize, const double zscale);
int    gcenter      (double *g, const byte *y, const int *imsize);
double gety         (const byte *y, const double *x, const int *imsize);
int    normal       (double *v, int nv);
int    swaponoff    (double **X, double *W, int *Non, const int N);
int    revert       (double **X, const double *xp, const int Non, const int N);
int    filtering    (double *W, const int N, const int nf, 
                     const byte *y, const double **X, const int* wlik,const int *imsize);
int    resampling   (int *nums, const double *W, const int N);
int    objinimage   (const double *x, const int *imsize, const int *objsize, const int *margin);
int    smoothing    (double **W2, int *buf, 
                     const int **H, const int **Nf, const int *V, const int *U, const int K, const int N);
int    expectation_i(double *e, const int   *Nf, const double **X, const int N);
int    expectation  (double *e, const double *W, const double **X, const int N);

/* NOTE --------------------------------------------------------------*/               
/* Use 'Mersenne Twister (mt19937ar.c)' as a rundom number generator. */
/* Source code is available in the author's webpage:                  */
/* http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html           */
/* The following two functions are defined in the 'mt19937ar.c' file. */
/* -------------------------------------------------------------------*/               
void          init_genrand   (unsigned long s);
double        genrand_real1  (void);


int main (int argc, char** argv){

  FILE   *fp,*fpp; int P=3;

  // Parameters
  double sgmt[3],sgms[3];
  int    wmax[3],wlik[3],objsize[3],margin[3];
  double alpha,beta,zscale;
  int    seed,cut,dploop,N,tune;
  char   in[256],out[256];

  int    imsize[4];
  int    i,j,t,k,n,n1,p,T,K,k1,k2,u,nf,L,ofs;
  int    nvx,root,rad;

  byte   *y;
  double **vx,**dist;
  double *eta,*eta0,*mov,**g;
  double **x0,***X,**W,**W2,**D,***noise; 
  double ***xyf,***xys,*xf,*xs,*xp;
  double *buf;
  int    **G,*U,*V,*Ks; 
  int    *nmem,*lbl,**Nf;
  int    **H,*tmp;

  fpp=fopen("param.txt", "r");if(!fpp){printf("File: \'param.txt\' Not Found.\n"); exit(1);}

  fscanf(fpp,"seed:%d\n",            &seed);
  fscanf(fpp,"cutoff:%d\n",          &cut);
  fscanf(fpp,"dploop:%d\n",          &dploop);
  fscanf(fpp,"root:%d\n",            &root); root--;
  fscanf(fpp,"N:%d\n",               &N);
  fscanf(fpp,"alpha:%lf\n",          &alpha);
  fscanf(fpp,"beta:%lf\n",           &beta);
  fscanf(fpp,"sgmt:%lf,%lf,%lf\n",   &sgmt[0],&sgmt[1],&sgmt[2]);
  fscanf(fpp,"sgms:%lf,%lf,%lf\n",   &sgms[0],&sgms[1],&sgms[2]);
  fscanf(fpp,"wmax:%d,%d,%d\n",      &wmax[0],&wmax[1],&wmax[2]);
  fscanf(fpp,"wlik:%d,%d,%d\n",      &wlik[0],&wlik[1],&wlik[2]);
  fscanf(fpp,"objsize:%d,%d,%d\n",   &objsize[0],&objsize[1],&objsize[2]);
  fscanf(fpp,"margin:%d,%d,%d\n",    &margin [0],&margin [1],&margin [2]);
  fscanf(fpp,"tune:%d\n",            &tune);
  fscanf(fpp,"input:%s\n",           in);
  fscanf(fpp,"output:%s\n",          out);
  fclose(fpp);fpp=NULL; init_genrand(seed);
  rad=objsize[0]/2; zscale=objsize[0]/objsize[2];

  fp =fopen(in,"rb");if(!fp){printf("File: \'%s\' Not Found.\n",in);exit(1);}

  fread(imsize,sizeof(int),4,fp); 
  L=imsize[0]*imsize[1]*imsize[2];T=imsize[3];
  printf("T=%d\n",T);

  y     = malloc   (L*sizeof(byte));
  vx    = calloc2d (NVXLIMIT,P);
  g     = calloc2d (T,P);
  lbl   = calloc   (NVXLIMIT,sizeof(int));
  nmem  = calloc   (NVXLIMIT,sizeof(int));
  x0    = calloc2d (NVXLIMIT,P);
  dist  = calloc2d (NVXLIMIT,NVXLIMIT);
  
  /* Determination of Initial positions */
  fread(y,L,sizeof(byte),fp);
  gcenter(g[0],y,imsize);cutoff(y,imsize,cut);localmax(vx,&nvx,y,imsize,wmax); 
  dpmeans (lbl,x0,dist,nmem,&K,(const double**)vx,nvx,P,objsize[1],zscale,dploop);

  X     = calloc3d (K,N,P);               // Particles
  noise = calloc3d (K,N,P);               // 
  W     = calloc2d (K,N);                 // Particle weights
  W2    = calloc2d (K,N);                 // Particle weights (smoothing)
  Nf    = calloc2i (K,N);                 // Filter ensemble
  H     = calloc2i (K,N);                 // Particle history 
  D     = calloc2d (K,K);                 // Distance matrix 
  G     = calloc2i (K,K);                 // Graph (MST)
  xyf   = calloc3d (T,K,P+1);             // Filter means
  xys   = calloc3d (T,K,P+1);             // Smoother means

  xp    = calloc   (P,    sizeof(double));  // Prediction mean (buffer)
  xf    = calloc   (P,    sizeof(double));  // Filter     mean (buffer)
  xs    = calloc   (P,    sizeof(double));  // Smoother   mean (buffer)
  eta0  = calloc   (P,    sizeof(double));
  eta   = calloc   (P,    sizeof(double));
  mov   = calloc   (P,    sizeof(double));
  V     = calloc   (K,    sizeof(int));     
  U     = calloc   (K,    sizeof(int));    
  Ks    = calloc   (T,    sizeof(int));   
  buf   = calloc   (K*N*P,sizeof(double));  
  tmp   = calloc   (N,    sizeof(double));

  /* Construction of MRF tree */ 
  for(k1=0;k1<K;k1++)for(k2=0;k2<K;k2++) D[k1][k2]=wdist(x0[k1],x0[k2],P,zscale);
  mstree(G,(const double**)D,K); 

  if(root<0)findroot(&root,y,(const double**)x0,K,imsize,zscale); makeiter(V,U,(const int**)G,K,root);
  for(i=0;i<K;i++)for(p=0;p<P;p++){k=V[i];xyf[0][k][p]=xys[0][k][p]=x0[k][p];}
  for(i=0;i<K;i++){k=V[i];xyf[0][k][3]=xys[0][k][3]=gety(y,(const double*)x0[k],imsize);}Ks[0]=K;
  printf("root=%d\n",root+1);
   
  /* Spatial particle filter */
  for(t=1;t<T;t++){
    if(t%10==9) fprintf(stderr,"t=%d\n",t+1);
    fread(y,L,sizeof(byte),fp); gcenter(g[t],y,imsize); 

    normal(buf,K*N*P);Ks[t]=K;
    for(i=0;i<K;i++)for(n=0;n<N;n++)for(p=0;p<P;p++)
      noise[V[i]][n][p]=buf[i+n*K+p*K*N]*((V[i]==root)?sgmt[p]:sgms[p]); 

    for(i=0;i<K;i++){k=V[i];u=U[k];assert(k>=0&&k<K&&u>=-1&&u<K); 
      /* Prediction */
      if(i){
        for(p=0;p<P;p++){
          eta [p]=xyf[t-1][k][p]-xyf[t-1][u][p];
          eta0[p]=xyf  [0][k][p]-xyf  [0][u][p];
          mov [p]=alpha*eta[p]+(1-alpha)*eta0[p];
          xp  [p]=xyf[t][u][p]+mov[p];
        }
        ofs=0;
        for(n=0;n<N;n++){n1=Nf[u][n];
          if(n1)for(j=0;j<n1;j++){assert(ofs>=0&&j+ofs<N); H[k][j+ofs]=n;
            for(p=0;p<P;p++) X[k][j+ofs][p]=X[u][n][p]+mov[p]+noise[k][n][p];
            if(wdist(X[k][j+ofs],X[u][n],P,zscale)<0.8*rad) /* Avoiding collisions */
              for(p=0;p<P;p++) X[k][j+ofs][p]=X[u][n][p]+mov[p]-0.5*noise[k][n][p];
          }
          ofs+=n1;
        }
      } 
      else {
        for(p=0;p<P;p++)xp[p]=xyf[t-1][k][p]+g[t][p]-g[t-1][p];
        for(n=0;n<N;n++)for(p=0;p<P;p++) X[k][n][p]=xyf[t-1][k][p]+(g[t][p]-g[t-1][p])+noise[k][n][p]; 
      }
    
      nf=N; 
      /* Filtering */
      if(objinimage(xp,imsize,objsize,margin)){
        filtering  (W[k],N,nf,y,(const double**)X[k],wlik,imsize);
        resampling (Nf[k],W[k],nf); 
      }
      else for(n=0;n<N;n++)Nf[k][n]=1;             

      expectation_i(xf,(const int*)Nf[k],(const double**)X[k],N);
      for(p=0;p<P;p++)xyf[t][k][p]=xf[p];xyf[t][k][3]=gety(y,(const double*)xf,imsize);

    }

    /* Fine-tuning by k-means */ 
    if(tune){
      for(i=0;i<K;i++)for(p=0;p<P;p++)x0[i][p]=(double)xyf[t][i][p];
      cutoff(y,imsize,cut);localmax(vx,&nvx,y,imsize,wmax); 
      kmeans(lbl,x0,dist,nmem,(const double**)vx,nvx,K,P,zscale,dploop,tune);
      for(i=0;i<K;i++){for(p=0;p<P;p++)xyf[t][i][p]=x0[i][p];xyf[t][i][3]=gety(y,(const double*)x0[i],imsize);}
    }
  } 

  write(out,(const double**)xyf,Ks,imsize,objsize,cut);

  return 0;
}


int expectation(double *e, const double *W, const double **X, const int N){
  int n,p,P=3; 
  for(p=0;p<P;p++)e[p]=0;
  for(p=0;p<P;p++)for(n=0;n<N;n++)e[p]+=W[n]*X[n][p];
  return 0;
}

int expectation_i(double *e, const int *Nf, const double **X, const int N){
  int n,p,P=3;
  for(p=0;p<P;p++)e[p]=0;
  for(p=0;p<P;p++)for(n=0;n<N;n++)e[p]+=Nf[n]*X[n][p];
  for(p=0;p<P;p++)e[p]/=(double)N;
  return 0;
}

inline int allzero(const double *W, const int N){
  int n,flag=1; for(n=0;n<N;n++)if(W[n]>=1.0e-10){flag=0;break;} return flag;
} 

int smoothing(double **W2, int *buf, 
              const int **H, const int **Nf, const int *V, const int *U, const int K, const int N){
  int i,n,k=-1; double val=0;

  for(i=0;i<K;i++)for(n=0;n<N;n++) W2[i][n]=Nf[i][n];
  for(i=1;i<K;i++){k=V[K-i];
    for(n=0;n<N;n++)buf[n]=0;
    for(n=0;n<N;n++)buf[H[k][n]]+=Nf[k][n];
    for(n=0;n<N;n++)W2 [U[k]][n]*=buf[n];
    if(allzero(W2[k],N))for(n=0;n<N;n++) W2[k][n]=Nf[k][n];
  }
  for(k=0;k<K;k++){val=0;
    for(n=0;n<N;n++)val+=W2[k][n];
    for(n=0;n<N;n++)W2[k][n]/=val;
  }

  return 0;
}
 
int findroot(int *root, const byte *y, const double **mx, const int K, const int *imsize, const double zscale){
  int k,kmax=-1; double val,max =0;
  for(k=0;k<K;k++){val=gety(y,mx[k],imsize);if(val>max){max=val;kmax=k;}}
  *root=kmax; assert(kmax!=-1);
  return 0;
}

int gcenter (double *g, const byte *y, const int *imsize){
  int  p,i,j,k,l,I,J,K,L; I=imsize[0];J=imsize[1];K=imsize[2];L=I*J*K;
  long s=0; double w,P=3; byte v;

  for(l=0;l<L;l++)s+=y[l]; assert(s>0);
  for(p=0;p<P;p++)g[p]=0;
  for(i=0;i<I;i++)for(j=0;j<J;j++)for(k=0;k<K;k++){v=y[i+j*I+k*I*J];
    if(v){w=v/(double)s;g[0]+=w*i;g[1]+=w*j;g[2]+=w*k;}
  }

  return 0;
}

int objinimage(const double *x, const int *imsize, const int *objsize, const int *margin){
  return    x[0]-objsize[0]-margin[0]>=0 &&  x[0]+objsize[0]+margin[0]<imsize[0]
         && x[1]-objsize[1]-margin[1]>=0 &&  x[1]+objsize[1]+margin[1]<imsize[1]
         && x[2]-objsize[2]-margin[2]>=0 &&  x[2]+objsize[2]+margin[2]<imsize[2];
}

int swaponoff(double **X, double *W, int *Non, const int N){
  int i=0,j=N-1,ct=0; double *tmp,val;

  while(1){ 
    while(i<N && W[i]> 1.0e-9) i++; 
    while(j>=0&& W[j]<=1.0e-9) j--; 
    if(i<j){ct++;
      tmp=X[j];X[j]=X[i];X[i]=tmp;
      val=W[j];W[j]=W[i];W[i]=val;
    } 
    else break;
  }*Non=N-ct;

  return 0;
}

int revert(double **X, const double *xp, const int Non, const int N){
  int n,p,P=3; for(n=Non-1;n<N;n++)for(p=0;p<P;p++) X[n][p]=xp[p]; return 0;
}

int filtering(double *W, const int N, const int nf, 
              const byte *y, const double **X, const int* wlik,const int *imsize){
  int n,i,j,k,a,b,c,A,B,C,I,J,K; double s,val=0,max=0;

  A=wlik  [0]; B=wlik  [1]; C=wlik  [2];
  I=imsize[0]; J=imsize[1]; K=imsize[2]; 

  for(n=0;n<nf;n++){s=0;i=floor(X[n][0]);j=floor(X[n][1]);k=floor(X[n][2]);
    for(a=-A;a<=A;a++)for(b=-B;b<=B;b++)for(c=-C;c<=C;c++)
      if(i+a>=0&&i+a<I&&j+b>=0&&j+b<J&&k+c>=0&&k+c<K)s+=y[i+a+(j+b)*I+(k+c)*I*J];
    W[n]=s;
  }

  for(n=0;n<nf;n++)if(W[n]>max) max=W[n];
  if(max>0){
    for(n=0;n<nf;n++)W[n]-=max;
    for(n=0;n<nf;n++)W[n]=exp(-W[n]*W[n]);
    for(n=0;n<nf;n++)val+=W[n];
    for(n=0;n<nf;n++)W[n]/=val;
    for(n=nf;n<N;n++)W[n]=0;
  }
  else {
    for(n=0;n<nf;n++)W[n]=1/(double)nf;
    for(n=nf;n<N;n++)W[n]=0;
  }

  return 0;
}

int resampling(int *nums, const double *probs, const int N){
  int i; double u =genrand_real1()/(double)N; 
  for(i=0;i<N;i++){
    nums[i]=floor((probs[i]-u)*N)+1; 
    u+=(nums[i]/(double)N)-probs[i];
  }
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

double gety(const byte *y, const double *x, const int *imsize){
  int p,P=3,l[3],out=0;  
  for(p=0;p<P;p++)if(x[p]<0||x[p]>=imsize[p]){out=1;break;} else l[p]=floor(x[p]);
  return out? 0:(double) y[l[0]+l[1]*imsize[0]+l[2]*imsize[0]*imsize[1]];
}

