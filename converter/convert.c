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
#include<opencv/cv.h>
#include<opencv/highgui.h>

double median (double *a, double *w, const int N);

int main(int argc, char **argv) {
  int      ct,l,N,M,a,b,c,i,j,k,t,A,B,C,I,J,K,T;
  double   *y0,*y1,*w0,*w1,max,mean,val,scale;
  int      ofs,size[4],wind[3];
  char     fn[1024],pfx[1024],typ[64];
  FILE     *fpo,*fpp;
  unsigned char *buf;
  IplImage *img;

  if(!(fpp=fopen("info.txt", "r"))){printf("File: \'info.txt\' Not Found.\n");exit(1);}
  fscanf(fpp,"prefix:%s\n",         pfx);
  fscanf(fpp,"imtype:%s\n",         typ);
  fscanf(fpp,"offset:%d\n",         &ofs); ofs++;
  fscanf(fpp,"imsize:%d,%d,%d,%d\n",&size[0],&size[1],&size[2],&size[3]);
  fscanf(fpp,"window:%d,%d,%d\n",   &wind[0],&wind[1],&wind[2]);
  fclose(fpp);fpp=NULL;

  I=size[0];J=size[1];K=size[2];T=size[3];N=I*J*K;
  A=wind[0];B=wind[1];C=wind[2];M=(2*A+1)*(2*B+1)*(2*C+1);

  y0 =malloc(N*sizeof(double)); w0=malloc(M*sizeof(double));
  y1 =malloc(N*sizeof(double)); w1=malloc(M*sizeof(double));
  buf=malloc(N*sizeof(unsigned char));
  
  fpo=fopen("data.bin","wb"); fwrite(size,sizeof(int),4,fpo);
  for(t=0;t<T;t++){printf("t=%d\n",t+1);
    /* Load images */
    for(k=0;k<K;k++){l=k+t*K+ofs;
      sprintf(fn,"%s%d.%s",pfx,l,typ); 
      img=cvLoadImage(fn,CV_LOAD_IMAGE_ANYCOLOR|CV_LOAD_IMAGE_ANYDEPTH);
      for(i=0;i<I;i++)for(j=0;j<J;j++) y0[i+j*I+(K-k-1)*I*J]=(double)cvGet2D(img,J-j-1,i).val[0]; 
      cvReleaseImage(&img);
    }
    /* Average filter */
    for(k=0;k<K;k++){val=0;
      for(i=0;i<I;i++)for(j=0;j<J;j++){val+=y0[i+j*I+k*I*J];} mean=val/(I*J);
      for(i=0;i<I;i++)for(j=0;j<J;j++){val =y0[i+j*I+k*I*J]-mean; y0[i+j*I+k*I*J]=(val>0)?val:0;}
    }
    /* Median filter */
    for(i=0;i<I;i++)for(j=0;j<J;j++)for(k=0;k<K;k++){ct=0;
      for(a=-A;a<=A;a++)for(b=-B;b<=B;b++)for(c=-C;c<=C;c++)
        if(i+a>=0&&i+a<I&&j+b>=0&&j+b<J&&k+c>=0&&k+c<K){w0[ct]=y0[i+a+(j+b)*I+(k+c)*I*J];ct++;}
      y1[i+j*I+k*I*J]=median(w0,w1,ct);
    } 
    /* Scaling & Reduction to 8bit */
    max=0;
    for(i=0;i<I;i++)for(j=0;j<J;j++)for(k=0;k<K;k++){val=y1[i+j*I+k*I*J];if(val>max)max=val;} scale=255.0/max;
    for(i=0;i<I;i++)for(j=0;j<J;j++)for(k=0;k<K;k++) y1[i+j*I+k*I*J]*=scale;
    for(i=0;i<I;i++)for(j=0;j<J;j++)for(k=0;k<K;k++) buf[i+j*I+k*I*J]=(unsigned char)y1[i+j*I+k*I*J]; 
    /* Output */ 
    fwrite(buf,sizeof(unsigned char),N,fpo); 
  }
  fclose(fpo);fpo=NULL;
  	
  return 0;  	
}


