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
#include<string.h>
#include<math.h>
#include<assert.h>
#include<time.h>
#include<unistd.h>
#include<ctype.h>
#include<omp.h>

#include"../common/util.h"
#include"../common/io.h"
#include"tree.h"

#define HEADER_SIZE     16
#define CONFIG_LINE_MAX 2048

int    findroot     (int *root, const byte *y, const double **mx, const double **D, const int K, const int *imsize);
int    gcenter      (double *g, const byte *y, const int *imsize);
double gety         (const byte *y, const double *x, const int *imsize);
int    normal       (double *v, int nv);
int    resampling   (int *nums, const double *W, const int N);
int    objinimage   (const double *x, const int *imsize, const int *objsize, const int *margin);
int    smoothing    (double **W2, int *buf, 
                     const int **H, const int **Nf, const int *V, const int *U, const int K, const int N);
int    expectation_i(double *e, const int   *Nf, const double **X, const int N);
int    expectation  (double *e, const double *W, const double **X, const int N);

int    subimage     (short *sy, const byte *y, const double *x, const int *wdw, const int *imsize);
int    filtering    (double *W, short *yt1, short **Y, const byte *y, const double **X, const short *yt0,
                     const int N, const int nf, const int *wdw, const double gamma, const double sgml, const int *imsize);
int    progress     (const int n, const int N, const int times, const int w, const double realtime, const double cputime);
int    writetraj    (const char *file, const double ***xy, const int T, const int K);
int    categorical  (const double *probs, const int N);
double wdist        (const double * x1, const double * x2, const int P, const double dz);
void   readPrms     (char *file);
void   printPrms    (void);

/* NOTE --------------------------------------------------------------*/               
/* Use 'Mersenne Twister (mt19937ar.c)' as a rundom number generator. */
/* Source code is available in the author's webpage:                  */
/* http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html           */
/* The following two functions are defined in the 'mt19937ar.c' file. */
/* -------------------------------------------------------------------*/               
void          init_genrand   (unsigned long s);
double        genrand_res53  (void);

// Parameters (with default values)
int    seed=1,N=1000,root=-1,order=1;
int    wlik[3]={6,6,2},margin[3]={2,2,0};
double alpha=0.8,beta=0.0,gma=0.1,rho=0.8,zscale=3.5,sgml=10.0;
double sgmt[3]={3.0,3.0,0.1},sgms[3]={1.5,1.0,0.1};
double sst=1.0,logit_c=30.0,logit_k=0.3;

double logit(double x, double k, double c){ return 1.0/(1.0+exp(-k*(x-c))); }

int main (int argc, char** argv){

  FILE   *fp; int P=3;
  char   *fcn,*fin,fdet[256];

  ImageHeader h; 
  struct _tracker det;

  int    imsize[4];
  int    i,j,t,k,l,n,n1,p,T,K,k1,k2,u,nf,L,Ls,ofs;
  int    tau,q,Q;

  byte   *y,*ytmp;
  double lambda;
  double **eta,**g,*probs;
  double **x0,***X,**W,**W2,**D,***noise; 
  double ***xyf,***xys,*xf,*xs,*xp;
  double *buf,val;
  int    **G,*U,*V,*Ks; 
  int    **Nf;
  int    **H,*tmp;
  short  **Y,**Yt1,**Yt0;
  clock_t cputime=clock();
  time_t  realtime=time(NULL);

  if(argc!=3){printf("./track <in:conf> <in:data>\n");exit(EXIT_SUCCESS);}
  fcn=argv[1];fin=argv[2]; strcpy(fdet,fin); strcat(fdet,".det");

  readPrms(fcn);lambda=wlik[0];  
  printPrms();

  init_genrand(seed); nf=(int)N*(1-beta);
  Ls=(2*wlik[0]+1)*(2*wlik[1]+1)*(2*wlik[2]+1);

  loadImageInfo(&h,fin);
  loadTracker(&det,fdet); K=det.nums[0];
  fp=fopen(h.fnRawImage,"rb");if(!fp){printf("File: \'%s\' Not Found.\n",h.fnRawImage);exit(1);}
  fseek(fp,HEADER_SIZE,SEEK_SET);
  imsize[0]=h.resolution.X;
  imsize[1]=h.resolution.Y;
  imsize[2]=h.resolution.Z;
  imsize[3]=h.nFrames; 
  zscale=h.voxelAspectRatio.z;
  L=imsize[0]*imsize[1]*imsize[2]; T=imsize[3];

  y     = malloc   (L*sizeof(byte));
  ytmp  = malloc   (L*sizeof(byte));
  g     = calloc2d (T,P);
  x0    = calloc2d (K,P);

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
  Yt1   = calloc2s (K,Ls);
  Yt0   = calloc2s (K,Ls);
  Y     = calloc2s (N,Ls);

  xp    = calloc   (P,    sizeof(double));  // Prediction mean (buffer)
  xf    = calloc   (P,    sizeof(double));  // Filter     mean (buffer)
  xs    = calloc   (P,    sizeof(double));  // Smoother   mean (buffer)
  eta   = calloc2d (order+1,P);
  probs = calloc   (order+1,sizeof(double));
  V     = calloc   (K,    sizeof(int));     
  U     = calloc   (K,    sizeof(int));    
  Ks    = calloc   (T,    sizeof(int));   
  buf   = calloc   (K*N*P,sizeof(double));  
  tmp   = calloc   (N,    sizeof(double));

  fread(y,L,sizeof(byte),fp);for(l=0;l<L;l++)ytmp[l]=y[l];gcenter(g[0],ytmp,imsize);
  for(k=0;k<K;k++)for(p=0;p<P;p++) x0[k][p]=det.data[p+k*(P+1)];

  /* Construction of MRF tree */ 
  for(k1=0;k1<K;k1++)for(k2=0;k2<K;k2++) D[k1][k2]=wdist(x0[k1],x0[k2],P,zscale);
  mstree(G,(const double**)D,K); 

  if(root<0){findroot(&root,y,(const double**)x0,(const double **)D,K,imsize);} makeiter(V,U,(const int**)G,K,root);
  for(i=0;i<K;i++)for(p=0;p<P;p++){k=V[i];xyf[0][k][p]=xys[0][k][p]=x0[k][p];}
  for(i=0;i<K;i++){k=V[i];xyf[0][k][3]=xys[0][k][3]=gety(y,(const double*)x0[k],imsize);}Ks[0]=K;

  for(i=0;i<K;i++){k=V[i];
    subimage(Yt0[k],y,x0[k],wlik,imsize);
    subimage(Yt1[k],y,x0[k],wlik,imsize);
  }

  printf("\n");
  printf("  o------------------o------o\n");
  printf("  | #Time points     | %4d |\n",T);
  printf("  | #Detected cells  | %4d |\n",K);
  printf("  | Root cell ID     | %4d |\n",root+1);
  printf("  o------------------o------o\n");
  printf("\n");

  /* Spatial particle filter */
  for(t=1;t<T;t++){Q=(t-1<=order)?t-1:order;assert(Q>=0);assert(Q<t);
    progress(t,T,50,50,(double)time(NULL)-realtime,(double)(clock()-cputime)/CLOCKS_PER_SEC);
    fread(y,L,sizeof(byte),fp); gcenter(g[t],y,imsize); 

    probs[0]=Q?1-alpha:1;for(q=1;q<=Q;q++)probs[q]=probs[q-1]*alpha;
    val=0;
    for(q=1;q<=Q;q++){val+=probs[q];} if(t!=1){assert(val>0);}
    for(q=1;q<=Q;q++){probs[q]/=val;} for(q=1;q<=Q;q++)probs[q]*=alpha;

    normal(buf,K*N*P);Ks[t]=K;
    for(i=0;i<K;i++)for(n=0;n<N;n++)for(p=0;p<P;p++)
      noise[V[i]][n][p]=buf[i+n*K+p*K*N]*((V[i]==root)?sgmt[p]:sgms[p]); 

    for(i=0;i<K;i++){k=V[i];u=U[k];assert(k>=0&&k<K&&u>=-1&&u<K); 
      /* Prediction */
      if(i){
        for(p=0;p< P;p++) xp[p]=xyf[t][u][p];
        for(q=0;q<=Q;q++){tau=q?t-q:0;assert(tau>=0);assert(tau<t);
          for(p=0;p<P;p++){eta[q][p] = xyf[tau][k][p]-xyf[tau][u][p];}
          for(p=0;p<P;p++){xp [p]   +=(xyf[tau][k][p]-xyf[tau][u][p])*probs[q];}
        }
        ofs=0;
        for(n=0;n<N;n++){n1=Nf[u][n];
          if(n1)for(j=0;j<n1;j++){assert(ofs>=0&&j+ofs<N); H[k][j+ofs]=n;
            for(p=0;p<P;p++) X[k][j+ofs][p]=X[u][n][p]+noise[k][n][p];
            q=categorical(probs,Q+1);for(p=0;p<P;p++)X[k][j+ofs][p]+=eta[q][p];
            if(wdist(X[k][j+ofs],X[u][n],P,zscale)<rho*lambda) /* Avoiding collisions */
              for(p=0;p<P;p++) X[k][j+ofs][p]-=0.5*noise[k][n][p];
          }
          ofs+=n1;
        }
      } 
      else {
        for(p=0;p<P;p++)xp[p]=xyf[t-1][k][p]+g[t][p]-g[t-1][p];
        for(n=0;n<N;n++)for(p=0;p<P;p++) X[k][n][p]=xp[p]+noise[k][n][p]; 
      }
    
      /* Filtering */
      if(objinimage(xp,imsize,wlik,margin)){
        filtering(W[k],Yt1[k],Y,y,(const double**)X[k],Yt0[k],N,nf,wlik,gma,sgml,imsize);
        resampling (Nf[k],W[k],nf); 
        for(n=nf;n<N;n++)Nf[k][n]=1;
      }
      else for(n=0;n<N;n++)Nf[k][n]=1;             

      expectation_i(xf,(const int*)Nf[k],(const double**)X[k],N);
      for(p=0;p<P;p++){xyf[t][k][p]=xf[p];} xyf[t][k][3]=gety(y,(const double*)xf,imsize);

    }

    /* Backward refinement */
    smoothing(W2,tmp,(const int**)H,(const int**)Nf,(const int*)V,(const int*)U,K,N);
    for(i=0;i<K;i++){k=V[i];
      expectation (xs,(const double*)W2[k],(const double**)X[k],N);
      for(p=0;p<P;p++)xys[t][k][p]=xs[p];
      xys[t][k][3]=gety(y,(const double*)xs,imsize);
    }

  } 
  writetraj(fin,(const double***)xyf,T,K);

  return 0;
}

void readPrms(char *file){
  int j,l,m,L,M=13; FILE *fp; char found,*s,*b, line[CONFIG_LINE_MAX];
  int  nargs[13]={1,1,1,1, 1,1,1,1,1,1, 3,3,3};
  char *form []={"%d","%lf","%lf,%lf,%lf","%d,%d,%d"};
  char *names[]={"seed","order","np","root","alpha","beta","gamma","rho","zscale","sgml","sgmt","sgms","wlik"};

  fp=fopen(file,"r");if(!fp){printf("File not found: %s\n",file);exit(EXIT_FAILURE);}
  while(fgets(line,CONFIG_LINE_MAX,fp)){found=0;
    if(strchr("#|+\n",line[0])){continue;} L=strlen(line);
    for(l=0;l<L;l++) line[l]=tolower(line[l]);
    for(m=0;m<M;m++){s=names[m];if(!strstr(line,s))continue; b=line+strlen(s);found=1;
      switch(m){
        case  0: j=sscanf(b,form[0],&seed  );            break;
        case  1: j=sscanf(b,form[0],&order );            break;
        case  2: j=sscanf(b,form[0],&N     );            break;
        case  3: j=sscanf(b,form[0],&root  ); root--;    break;
        case  4: j=sscanf(b,form[1],&alpha );            break;
        case  5: j=sscanf(b,form[1],&beta  );            break;
        case  6: j=sscanf(b,form[1],&gma   );            break;
        case  7: j=sscanf(b,form[1],&rho   );            break;
        case  8: j=sscanf(b,form[1],&zscale);            break;
        case  9: j=sscanf(b,form[1],&sgml  );            break;
        case 10: j=sscanf(b,form[2],sgmt,sgmt+1,sgmt+2); break;
        case 11: j=sscanf(b,form[2],sgms,sgms+1,sgms+2); break;
        case 12: j=sscanf(b,form[3],wlik,wlik+1,wlik+2); break;
      } if(j!=nargs[m]) printf("%s: Incorrect argument(s).\n",names[m]);  
    } if(!found) printf("Can't read the line: %s",line);
  } fclose(fp);

  return;
}

double wdist(const double * x1, const double * x2, const int P, const double dz){ 
  int p; double v,d=0.0; 
  for(p=0;p<P;p++){v=(x1[p]-x2[p])*(x1[p]-x2[p]); d+=(p<P-1)?v:dz*dz*v;}
  return sqrt(d);
} 

int writetraj(const char *file, const double ***xy, const int T, const int K){
  int t,k,num; char fni[256],fnb[256],fnt[256]; FILE *fpi,*fpb;
  strcpy(fni,file);strcat(fni,".trk.info");
  strcpy(fnb,file);strcat(fnb,".trk.bin" );
  strcpy(fnt,file);strcat(fnt,".cellnames.txt" );

  fpi=fopen(fni,"w");assert(fpi);
  fprintf(fpi,"HeaderSize           %d\n",HEADER_SIZE);
  fprintf(fpi,"NumberOfFrames       %d\n",T);
  fprintf(fpi,"MaxNumberOfTrackers  %d\n",K);
  fprintf(fpi,"TrackerNames         %s\n",fnt);
  fclose(fpi);

  fpb=fopen(fnb,"wb");assert(fpb);fseek(fpb,HEADER_SIZE,SEEK_SET);
  for(t=0;t<T;t++)for(k=0;k<K;k++){num=fwrite(xy[t][k],sizeof(double),4,fpb);assert(num==4);}
  fclose(fpb);

  return 0;
}

void printPrms(void){
  printf("\n");
  printf("*** Parameters ********************\n");
  printf("  np:     %d\n",               N);
  printf("  alpha:  %.2lf\n",            alpha);
  printf("  sgmt:   %.2lf,%.1lf,%.1lf\n",sgmt[0],sgmt[1],sgmt[2]);
  printf("  sgms:   %.2lf,%.1lf,%.1lf\n",sgms[0],sgms[1],sgms[2]);
  printf("  sgml:   %.2lf\n",            sgml);
  printf("  wlik:   %d,%d,%d\n",         wlik[0],wlik[1],wlik[2]);
  printf("***********************************\n");
  //printf("  zscale: %.2lf\n",            zscale);
  //printf("  seed:   %d\n",               seed);
  //printf("  rho:    %.2lf\n",            rho);
  //printf("  order:  %d\n",               order);
  //printf("  beta:   %.2lf\n",            beta);
  //printf("  gma:    %.2lf\n",            gma);
  //printf("  margin: %d,%d,%d\n",         margin[0],margin[1],margin[2]);
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
  int n,flag=1; for(n=0;n<N;n++)if(W[n]>=1.0e-250){flag=0;break;} return flag;
} 

int smoothing(double **W2, int *buf, 
              const int **H, const int **Nf, const int *V, const int *U, const int K, const int N){
  int i,n,k=-1; double val=0;

  for(i=0;i<K;i++)for(n=0;n<N;n++) W2[i][n]=Nf[i][n];
  for(i=1;i<K;i++){k=V[K-i];
    for(n=0;n<N;n++)buf[n]=0;
    for(n=0;n<N;n++)buf[H[k][n]]+=Nf[k][n];
    for(n=0;n<N;n++)W2 [U[k]][n]*=buf[n];
    if(allzero(W2[U[k]],N))for(n=0;n<N;n++) W2[U[k]][n]=Nf[k][n];
  }
  for(k=0;k<K;k++){val=0;
    for(n=0;n<N;n++)val+=W2[k][n];
    assert(fabs(val)>1e-10);
    for(n=0;n<N;n++)W2[k][n]/=val;
  }

  return 0;
}
 
int findroot(int *root, const byte *y, const double **mx, const double **D, const int K, const int *imsize){
  sortbox *sb; int n,nn=8;
  int k,kmax=-1; double v,*val,max=0,th1=40,th2=0.20;

  sb  =(sortbox*)calloc(K,sizeof(sortbox));
  val =(double*) calloc(K,sizeof(double ));

  /* Computing average distance from nearest neighbors */
  for(k=0;k<K;k++){v=0;
    prepare_sortbox(sb,D[k],K);qsort(sb,K,sizeof(sortbox),cmp_sortbox); assert(sb[0].val==0);
    v=0;for(n=1;n<=nn;n++)v+=sb[n].val; val[k]=v/nn;
  }
  for(k=0;k<K;k++)val[k]*=1/(1+exp(th1-(double)gety(y,mx[k],imsize)));
  for(k=0;k<K;k++){v=mx[k][0]/(double)imsize[0];val[k]*=(v>th2&&v<1-th2)?1:0;}
  max=0;for(k=0;k<K;k++)if(val[k]>max){max=val[k];kmax=k;}
  *root=kmax; assert(kmax!=-1);

  free(sb);free(val);
  return 0;
}

int gcenter (double *g, const byte *y, const int *imsize){
  int  p,i,j,k,l,I,J,K,L,P=3; I=imsize[0];J=imsize[1];K=imsize[2];L=I*J*K;
  long s=0; double w; byte v;

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


int categorical(const double *probs, const int N){
  int i; double u =genrand_res53();
  for(i=0;i<N;i++){u-=probs[i];if(u<=0)break;}
  return i;
}

int resampling(int *nums, const double *probs, const int N){
  int i; double u =genrand_res53()/(double)N;
  for(i=0;i<N;i++){
    nums[i]=floor((probs[i]-u)*N)+1; 
    u+=(nums[i]/(double)N)-probs[i];
  }
  return 0;
}

int normal(double *v, int nv){
  int i,n=(nv/2),odd=(nv%2==1)?1:0; double u1,u2;
  if(n)for(i=0;i<n;i++){u1=genrand_res53();u2=genrand_res53();
    v[2*i  ]=sqrt(-2*log(u1))*cos(2*M_PI*u2);
    v[2*i+1]=sqrt(-2*log(u1))*sin(2*M_PI*u2);
  } 
  if(odd){u1=genrand_res53();u2=genrand_res53();v[nv-1]=sqrt(-2*log(u1))*cos(2*M_PI*u2);}
  return 0;
}

double gety(const byte *y, const double *x, const int *imsize){
  int p,P=3,l[3],out=0;  
  for(p=0;p<P;p++)if(x[p]<=0.0||x[p]>=imsize[p]){out=1;break;} else l[p]=floor(x[p]);
  return out? 0:(double) y[l[0]+l[1]*imsize[0]+l[2]*imsize[0]*imsize[1]];
}

int subimage(short *sy, const byte *y, const double *x, const int *wdw, const int *imsize){
  int i,j,k,a,b,c,A,B,C,I,J,K,l=-1;

  A=wdw   [0]; B=wdw   [1]; C=wdw   [2];
  I=imsize[0]; J=imsize[1]; K=imsize[2];

  i=floor(x[0]);j=floor(x[1]);k=floor(x[2]);
  for(a=-A;a<=A;a++)for(b=-B;b<=B;b++)for(c=-C;c<=C;c++){l++;
    sy[l]=(i+a>=0&&i+a<I&&j+b>=0&&j+b<J&&k+c>=0&&k+c<K)? y[i+a+(j+b)*I+(k+c)*I*J]:-1;
  }

  return 0;
}

inline int wmeanimage(short *wy, const short **Y, const double *W, const int N, const int *wdw){
  int ct,n,l,L=(2*wdw[0]+1)*(2*wdw[1]+1)*(2*wdw[2]+1); double val,sum;
  for(l=0;l<L;l++){ct=val=sum=0;
    for(n=0;n<N;n++)if(Y[n][l]){ct++;sum+=W[n];val+=W[n]*Y[n][l];}
    wy[l]=ct?(short)(val/sum):-1;
  }
  return 0;
}

inline int mergeimage(short *y, const short *y0, const int *wdw, const double gamma){
  int l,L=(2*wdw[0]+1)*(2*wdw[1]+1)*(2*wdw[2]+1);
  for(l=0;l<L;l++) y[l]=(y[l]>=0&&y0[l]>=0)?gamma*y[l]+(1-gamma)*y0[l]:-1;
  return 0;
}

inline double likelihood(const short *y, const short *y0, const int *wdw, const double sgml){
  int l,ct=0,L=(2*wdw[0]+1)*(2*wdw[1]+1)*(2*wdw[2]+1), lim=1; double tmp,val=0;
  for(l=0;l<L;l++){if((y[l]<0||y0[l]<0)||(y[l]==0&&y0[l]==0))continue; ct++;
    tmp=y[l]-(double)y0[l]; val+=tmp*tmp;
  } val = ct>lim?exp(-val/(ct*sgml*sgml)):0.0;
  return val;
}

int filtering (double *W/*O*/, short *yt1/*IO*/, short **Y/*B*/, const byte *y, const double **X,  const short *yt0,
               const int N, const int nf, const int *wdw, const double gamma, const double sgml, const int *imsize){
  int n; double val=0;

  #pragma omp parallel for
  for(n=0;n<nf;n++) subimage(Y[n],y,X[n],wdw,imsize);
  #pragma omp parallel for
  for(n=0;n<nf;n++) W[n]=likelihood(Y[n],yt1,wdw,sgml);

  for(n=0;n<nf;n++) val+=W[n];
  if(val>0)for(n=0;n<N;n++)W[n]/=val;
  else{assert(nf);
    for(n=0;n<nf;n++)W[n]=1.0/nf;
    for(n=nf;n<N;n++)W[n]=0;
  }

  //wmeanimage(yt1,(const short**)Y,W,N,wdw);
  mergeimage(yt1,yt0,wdw,gamma);

  return 0;
}
 
int progress(const int n, const int N, const int times, const int w, const double realtime, const double cputime) {
  int i,c; double r;

  if(n%(1+N/times)!=0){r=(double)(n+1)/N;c=r*w;
    printf("  %3d%% [", (int)(r*100));
    for(i=0;i<c;i++) printf("=");
    for(i=c;i<w;i++) printf(" ");
    printf("]\n  Time: %d (sec)   CPU Time: %.2f (sec)",(int)realtime,cputime);
    if(n<N-1)printf("\n\033[F\033[J\033[F\033[J");
    else     printf("\n\n");
  }

  return 0;
}

