#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"../util.h"
#include"../dpmeans.h"
#include"../io.h"

#define NVXLIMIT 10240 
#define KLIMIT   1024 

int main (int argc, char** argv){

  const int P=3;
  int    nvx,t,k,T,K,l1,l2,l3,L1,L2,L3,L;
  double **x,***Xy,**vx,**dist;
  int    *Ks,*nmem,*lbl;
  byte   *y;
  char   in[256],out[256];
  FILE   *fp;

  // Parameters
  int    dploop,cut; double zscale;  
  int    imsize[4],wmax[3],objsize[3];

  fp=fopen("conf-detect.txt", "r");if(!fp){printf("File: \'conf-detect.txt\' Not Found.\n"); exit(1);}
  fscanf(fp,"cutoff:%d\n",          &cut);
  fscanf(fp,"dploop:%d\n",          &dploop);
  fscanf(fp,"wmax:%d,%d,%d\n",      &wmax[0],&wmax[1],&wmax[2]);
  fscanf(fp,"objsize:%d,%d,%d\n",   &objsize[0],&objsize[1],&objsize[2]);
  fscanf(fp,"input:%s\n",           in);
  fscanf(fp,"output:%s\n",          out);
  fclose(fp);fp=NULL; 

  zscale = objsize[0]/(double)objsize[2];

  /* Input, Allocation, & Initialization */
  fp=fopen(in,"rb");if(!fp){printf("File: \'%s\' Not Found.\n",in); exit(1);}
  fread(imsize,sizeof(int),4,fp);
  L1=imsize[0];L2=imsize[1];L3=imsize[2];T=imsize[3];
  L=L1*L2*L3;
 
  for(t=0;t<4;t++)printf("imsize[%d]=%d\n",t,imsize[t]);

  y     = calloc   (L, sizeof(byte));
  lbl   = calloc   (L, sizeof(int));
  Ks    = calloc   (T, sizeof(int));
  nmem  = calloc   (NVXLIMIT,sizeof(int));
  vx    = calloc2d (NVXLIMIT,   P  );
  Xy    = calloc3d (T,KLIMIT,   P+1);
  x     = calloc2d (NVXLIMIT,   P  );
  dist  = calloc2d (NVXLIMIT,  NVXLIMIT);

  for(t=0;t<T;t++){fread(y,L,sizeof(byte),fp); 
    if(t%100==99) printf("t=%d\n",t+1);
    cutoff   (y,imsize,(unsigned char)cut);
    localmax (vx,&nvx,y,imsize,wmax);
    dpmeans  (lbl,x,dist,nmem,&K,(const double**)vx,nvx,P,objsize[1],zscale,dploop);

    Ks[t]=K;
    for(k=0;k<K;k++){
//      if(t>846)printf("(A)k,k1,y %d,%d,%d\n",k,k1,y[l1+l2*L1+l3*L1*L2]);
//      if(t>846)printf("(A)l: %d,%d,%d\n",l1,l2,l3);
//      if(t>846)printf("x: %f\n",x[k][0]);
      l1=floor(x[k][0]); Xy[t][k][0]=l1;
      l2=floor(x[k][1]); Xy[t][k][1]=l2;
      l3=floor(x[k][2]); Xy[t][k][2]=l3; 
      Xy[t][k][3]=y[l1+l2*L1+l3*L1*L2];
    } 
  }

  write("result.bin",(const double**)Xy,Ks,imsize,objsize,cut);

  return 0;
}




















   
