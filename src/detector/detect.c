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
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<ctype.h>
#include"../common/util.h"
#include"../common/io.h"
#include"../common/mat3.h"

#define NVXLIMIT        4096
#define CONFIG_LINE_MAX 2048
#define MAX_WINDOW_SIZE 256
#define PRT(n) printf("%d\n",n)

void LoGfilter(double *y1, const double *y, const int imsize[3], const double sd, const int wsize[3], const double dz){
  int i,j,J,l,m,n,L,M,N,LM,w,s[3],jx,jy,jz,Jx,Jy,Jz,Wx,Wy,Wz; double val,G[3][MAX_WINDOW_SIZE],H[3][MAX_WINDOW_SIZE];
  /* check arguments */
  for(i=0;i<3;i++) assert(MAX_WINDOW_SIZE>=wsize[i]); assert(sd>0);

  L=imsize[0]; M=imsize[1]; N=imsize[2]; LM=L*M;J=LM*N;
  Wx=wsize[0]; s[0]=1;  jx=Wx*s[0]; Jx=J-jx;
  Wy=wsize[1]; s[1]=L;  jy=Wy*s[1]; Jy=J-jy;
  Wz=wsize[2]; s[2]=LM; jz=Wz*s[2]; Jz=J-jz;

  /* initialization */
  for(j=0;j<J;j++) y1[j]=0;

  /* build Gauss window */
  for(w=-Wx;w<Wx;w++) {val=   w/sd;val*=val; G[0][w+Wx]=exp(-0.5*val)/sd;}
  for(w=-Wy;w<Wy;w++) {val=   w/sd;val*=val; G[1][w+Wy]=exp(-0.5*val)/sd;}
  for(w=-Wz;w<Wz;w++) {val=dz*w/sd;val*=val; G[2][w+Wz]=exp(-0.5*val)/sd;}
  /* build window (negative 2nd derivative) */
  for(w=-Wx;w<Wx;w++) {val=   w/sd;val*=val; H[0][w+Wx]=exp(-0.5*val)*(1.0-val)/(sd*sd);}
  for(w=-Wy;w<Wy;w++) {val=   w/sd;val*=val; H[1][w+Wy]=exp(-0.5*val)*(1.0-val)/(sd*sd);}
  for(w=-Wz;w<Wz;w++) {val=dz*w/sd;val*=val; H[2][w+Wz]=exp(-0.5*val)*(1.0-val)/(sd*sd);}

  /* convolution */
  for(j=0;j<J;j++){l=j%L;     for(w=-Wx;w<Wx;w++)if(l+w>=0&&l+w<L)y1[j]+=H[0][w+Wx]*y[j+w*s[0]];}
  for(j=0;j<J;j++){m=(j%LM)/L;for(w=-Wy;w<Wy;w++)if(m+w>=0&&m+w<M)y1[j]+=G[1][w+Wy]*y[j+w*s[1]];}
  for(j=0;j<J;j++){n=j/LM;    for(w=-Wz;w<Wz;w++)if(n+w>=0&&n+w<N)y1[j]+=G[2][w+Wz]*y[j+w*s[2]];}

  for(j=0;j<J;j++){l=j%L;     for(w=-Wx;w<Wx;w++)if(l+w>=0&&l+w<L)y1[j]+=G[0][w+Wx]*y[j+w*s[0]];}
  for(j=0;j<J;j++){m=(j%LM)/L;for(w=-Wy;w<Wy;w++)if(m+w>=0&&m+w<M)y1[j]+=H[1][w+Wy]*y[j+w*s[1]];}
  for(j=0;j<J;j++){n=j/LM;    for(w=-Wz;w<Wz;w++)if(n+w>=0&&n+w<N)y1[j]+=G[2][w+Wz]*y[j+w*s[2]];}

  for(j=0;j<J;j++){l=j%L;     for(w=-Wx;w<Wx;w++)if(l+w>=0&&l+w<L)y1[j]+=G[0][w+Wx]*y[j+w*s[0]];}
  for(j=0;j<J;j++){m=(j%LM)/L;for(w=-Wy;w<Wy;w++)if(m+w>=0&&m+w<M)y1[j]+=G[1][w+Wy]*y[j+w*s[1]];}
  for(j=0;j<J;j++){n=j/LM;    for(w=-Wz;w<Wz;w++)if(n+w>=0&&n+w<N)y1[j]+=H[2][w+Wz]*y[j+w*s[2]];}

  /* scale multiplication */
  for(j=0;j<J;j++) y1[j]*=sd;

  return;
}

double wdist(const double * x1, const double * x2, const int P, const double dz){
  int p; double v,d=0.0;
  for(p=0;p<P;p++){v=(x1[p]-x2[p])*(x1[p]-x2[p]); d+=(p<P-1)?v:dz*dz*v;}
  return sqrt(d);
}

int dpmeans(
      int           *  z,    /* OUTPUT: cluster id | size N           */
      double        ** m,    /* OUTPUT: means      | size K(N) x P    */
      double        ** d,    /* OUTPUT: distance   | size N    x K(N) */
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
  int    kmin=-1,Kt=1;
  double *tmp,min=1.0E200;

  assert(P==2||P==3);

  /* Initialization */
  for(n=0;n<N;n++){z[n]=0;l[n]=0;}
  for(k=0;k<N;k++){for(n=0;n<N;n++)d[n][k]=0;for(p=0;p<P;p++)m[k][p]=0;}
  for(p=0;p<P;p++)m[0][p]=X[0][p];

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
      if(!l[k]) {Kt--;tmp=m[k];m[k]=m[Kt];m[Kt]=tmp;l[k]=l[Kt];}
      m[k][p]/=(double)l[k];
    }
  }
  *K=Kt;

  return 0;
}

int cutoff  (byte *y,const int *imsize, const byte cut){
  int l,L=imsize[0]*imsize[1]*imsize[2]; for(l=0;l<L;l++)if(y[l]<cut)y[l]=0;
  return 0;
}

int localmax(double **vx, int *nvx, const byte *y, const int *imsize, const int *width){
  int  i,j,k,l,l0,I,J,K,a,b,c,A,B,C,n=0;
  I=imsize[0];J=imsize[1];K=imsize[2];A=width[0];B=width[1];C=width[2]; *nvx=0;

  for(i=0;i<I;i++)for(j=0;j<J;j++)for(k=0;k<K;k++){l0=i+j*I+k*I*J;if(y[l0]==0) goto tag;
    for(a=-A;a<=A;a++)for(b=-B;b<=B;b++)for(c=-C;c<=C;c++){l=i+a+(j+b)*I+(k+c)*I*J;
      if((i+a>=0)&&(i+a<I)&&(j+b>=0)&&(j+b<J)&&(k+c>=0)&&(k+c<K)&&l!=l0)if(y[l0]<y[l])goto tag;
    } if(n>=NVXLIMIT){printf("Increase the cutoff value and retry (Abort).\n");exit(EXIT_FAILURE);}
    vx[n][0]=(double)i;vx[n][1]=(double)j;vx[n][2]=(double)k;n++; tag:;
  }*nvx=n;

  return 0;
}

int lcov(double S[3][3], const double g[3], const byte *y, const int *imsize, const int *width, const double dz){
  int  i,j,k,l,I,J,K,a,b,c,A,B,C,p,q,ct=0,flag=0; double v,N;

  I=imsize[0];J=imsize[1];K=imsize[2];A=width[0];B=width[1];C=width[2];
  i=(int)g[0];j=(int)g[1];k=(int)g[2];

  for(p=0;p<3;p++)for(q=0;q<3;q++)S[p][q]=0;

  for(a=-A;a<=A;a++)for(b=-B;b<=B;b++)for(c=-C;c<=C;c++){
    if(i+a<0||i+a>=I||j+b<0||j+b>=J||k+c<0||k+c>=K) continue;
    l=i+a+(j+b)*I+(k+c)*I*J;ct+=y[l]; v=dz*c;

    S[0][0]+=y[l]*a*a; S[0][1]+=y[l]*a*b;
    S[1][1]+=y[l]*b*b; S[1][2]+=y[l]*b*v;
    S[2][2]+=y[l]*v*v; S[0][2]+=y[l]*a*v;
  } N=(double)ct;
  if(det3(S)<1e-5) {flag=1;goto tag;}

  S[0][0]/=N; S[0][1]/=N; S[1][0]=S[0][1];
  S[1][1]/=N; S[1][2]/=N; S[2][1]=S[1][2];
  S[2][2]/=N; S[0][2]/=N; S[2][0]=S[0][2]; tag:

  return flag;
}

int segment(byte *y, const int id, const double g[3], const double S[3][3],
            const double sd, const int *imsize, const int *width, const double dz){
  int  i,j,k,l,I,J,K,a,b,c,A,B,C,flag=0; double x[3],iS[3][3];

  I=imsize[0];J=imsize[1];K=imsize[2];A=width[0];B=width[1];C=width[2];
  i=(int)g[0];j=(int)g[1];k=(int)g[2];

  if(det3(S)<1e-5){flag=1;goto tag;} inv3(iS,S);

  for(a=-A;a<=A;a++)for(b=-B;b<=B;b++)for(c=-C;c<=C;c++){
    if(i+a<0||i+a>=I||j+b<0||j+b>=J||k+c<0||k+c>=K) continue;
    l=i+a+(j+b)*I+(k+c)*I*J; x[0]=a;x[1]=b;x[2]=dz*c;
    y[l]=quad(x,iS)<sd?y[l]=id:y[l];
  }
  tag:

  return flag;
}

/* Scaling & Reduction to 8bit */
void scaling(byte *buf, const double *y, const int size){
    int p; double min=1e250,max=-1e250,scale;
    for(p=0;p<size;p++) max =(y[p]>max)?y[p]:max;
    for(p=0;p<size;p++) min =(y[p]<min)?y[p]:min;
    scale=255.0/(max-min);
    for(p=0;p<size;p++) buf[p]=(byte)(scale*(y[p]-min));
}

double gety(const byte *y, const double *x, const int *imsize){
  int p,P=3,l[3],out=0;
  for(p=0;p<P;p++)if(x[p]<=0.0||x[p]>=imsize[p]){out=1;break;} else l[p]=floor(x[p]);
  return out? 0:(double) y[l[0]+l[1]*imsize[0]+l[2]*imsize[0]*imsize[1]];
}

struct {
  int    cutoff,dploop,maxK;
  int    wpeak[3],wcell[3];
  double sd,lmd,dz;
} prms;


void readPrms (const char *file){
  int j,n,l,L,N=8; FILE *fp; char *s,*b, line[CONFIG_LINE_MAX];
  int   nargs[]={1,1,1,1,3,3,1,1};
  char *form []={"%d","%d,%d,%d","%lf"};
  char *names[]={"intensity cutoff","number of loops in dpmeans","max number of cells",
                 "lambda in dpmeans","peak bounding box","cell bounding box","voxel aspect ratio z",
                 "standard dev. in log"};

  fp=fopen(file,"r");if(!fp){printf("File not found: %s\n",file);exit(EXIT_FAILURE);}
  while(fgets(line,CONFIG_LINE_MAX,fp)){
    if(strchr("#|+\n",line[0])) continue; L=strlen(line);
    for(l=0;l<L;l++) line[l]=tolower(line[l]);
    for(n=0;n<N;n++){s=names[n];if(!strstr(line,s))continue; b=line+strlen(s);
      switch(n){
        case  0: j=sscanf(b,form[0],&(prms.cutoff));                                      break;
        case  1: j=sscanf(b,form[0],&(prms.dploop));                                      break;
        case  2: j=sscanf(b,form[0],&(prms.maxK));                                        break;
        case  3: j=sscanf(b,form[2],&(prms.lmd));                                         break;
        case  4: j=sscanf(b,form[1],&(prms.wpeak[0]),&(prms.wpeak[1]),&(prms.wpeak[2]));  break;
        case  5: j=sscanf(b,form[1],&(prms.wcell[0]),&(prms.wcell[1]),&(prms.wcell[2]));  break;
        case  6: j=sscanf(b,form[2],&(prms.dz));                                          break;
        case  7: j=sscanf(b,form[2],&(prms.sd));                                          break;
      } if(j!=nargs[n]) printf("Can't read the line: %s",line);
    }
  } fclose(fp);
  return;
}

void printPrms(){
  printf("\nParameters:\n");
  printf("  intensity cutoff            %d\n",       prms.cutoff);
  printf("  number of loops in dpmeans  %d\n",       prms.dploop);
  printf("  max number of cells         %d\n",       prms.maxK);
  printf("  lambda in dpmeans           %f\n",       prms.lmd);
  printf("  peak bounding box           %d,%d,%d\n", prms.wpeak[0],prms.wpeak[1],prms.wpeak[2]);
  printf("  cell bounding box           %d,%d,%d\n", prms.wcell[0],prms.wcell[1],prms.wcell[2]);
  printf("  standard dev. in LoG        %f\n\n",     prms.sd);
}

int main(int argc, char **argv) {
  byte   *y,*yseg; ImageHeader hI;
  double *buf,*buf1,dz;
  int    l,t,k,I,J,L,uT,T,K,hsize,nvx,num,imsize[4];
  int    *lbl,*nmem;
  double **vx,**vx1,**x0,**dist;
  FILE   *fpd,*fpt,*fpti,*fptn;
  char   fnt[MAX_FILENAME_LENGTH],fnti[MAX_FILENAME_LENGTH],fntb[MAX_FILENAME_LENGTH],fntn[MAX_FILENAME_LENGTH];

  if(argc<3){printf("./detect <conf> <data>\n");exit(EXIT_SUCCESS);}

  readPrms(argv[1]);
  printPrms();

  loadImageInfo(&hI,argv[2]);
  hsize=hI.rawHeaderSize;

  I=imsize[0]=hI.resolution.X;
  J=imsize[1]=hI.resolution.Y;
  K=imsize[2]=hI.resolution.Z;
  T=imsize[3]=hI.nFrames; L=I*J*K;
  dz=hI.voxelAspectRatio.z;

  uT=1;
  strcpy(fnt, argv[2]); strcat(fnt, ".det");
  strcpy(fntn,argv[2]); strcat(fntn,".cellnames.txt");
  strcpy(fnti,argv[2]); strcat(fnti,".det.info");
  strcpy(fntb,argv[2]); strcat(fntb,".det.bin" );

  printf("4D image:\n");
  printf("  Filename:   %s\n",           argv[2]);
  printf("  Resolution: %d x %d x %d\n", I,J,K);
  printf("  Z aspect:   %.2lf\n\n",      dz);

  printf("Output:\n");
  printf("  output file (centroid):     %s\n\n",  fnt);

  fpti=fopen_e(fnti,"w");
  fprintf(fpti,"HeaderSize           %d\n",hsize);
  fprintf(fpti,"NumberOfFrames       %d\n",uT);
  fprintf(fpti,"MaxNumberOfTrackers  %d\n",prms.maxK);
  fprintf(fpti,"TrackerNames         %s\n",fntn);
  fclose(fpti);

  fpd=fopen (hI.fnRawImage,"rb");assert(fpd);
  fpt=fopen (fntb,         "wb");assert(fpt);

  y     = (byte*)  malloc (L*sizeof(byte));
  yseg  = (byte*)  malloc (L*sizeof(byte));
  buf   = (double*)malloc (L*sizeof(double));
  buf1  = (double*)malloc (L*sizeof(double));
  lbl   = (int*)   calloc (NVXLIMIT,sizeof(int));
  nmem  = (int*)   calloc (NVXLIMIT,sizeof(int));
  vx    = calloc2d (NVXLIMIT,3);
  vx1   = calloc2d (NVXLIMIT,3);
  x0    = calloc2d (NVXLIMIT,4);
  dist  = calloc2d (NVXLIMIT,NVXLIMIT);

  fseek(fpd,hsize,SEEK_SET);
  fseek(fpt,hsize,SEEK_SET);

  for(t=0;t<uT;t++){ num=fread(y,sizeof(byte),L,fpd); assert(num==L);
    for(l=0;l<L;l++) {yseg[l]=0;buf1[l]=0;buf[l]=(double)y[l];}
    LoGfilter (buf1,buf,imsize,prms.sd,prms.wcell,dz);
    scaling   (y,buf1,L);
    cutoff    (y,imsize,prms.cutoff);
    localmax  (vx,&nvx,y,imsize,prms.wpeak);
    dpmeans   (lbl,x0, dist,nmem,&K,(const double**)vx,nvx,3,prms.lmd,prms.dz,prms.dploop);

    for(k=0;k<prms.maxK;k++) {
      x0[k][3]=k<K?gety(y,x0[k],imsize):nan("NaN");
      num=fwrite(x0[k],sizeof(double),4,fpt); assert(num==4);
    }
  }
  if(t==1) printf("Result:\n  %d cells found (among %d LoG peaks).\n\n", K,nvx);
  fptn=fopen(fntn,"w"); for(k=0;k<K;k++) fprintf(fptn,"%d\n",k+1);
  fclose(fptn);
  fclose(fpd);
  fclose(fpt);

  return 0;
}
