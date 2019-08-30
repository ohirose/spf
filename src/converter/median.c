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

double median(double *a, double *w, const int N){
  const int c = N/2;
  int     i,j,k,l,u,e,ofs,size;
  double  *tmp, p;/*pivot*/

  int t=0;
  size=N;ofs=0;
  while(1){i=j=k=0;p=a[0];e=1;t++;
    for(i=1;i<size;i++){
      if(a[i]< p) {a[j]=a[i];j++;}
      if(a[i]> p) {w[k]=a[i];k++;}
      if(a[i]==p) e++; 
    } l=ofs+j;u=l+e-1;

    #ifdef DEBUG
    printf("[t=%d] BEFORE: ofs=%d,c=%d,l=%d,u=%d,size=%d\n",t,ofs,c,l,u,size);
    printf("[t=%d] L:%c",t,j==0?'\n':' ');for(i=0;i<j;i++)printf("%f%c",a[i],i==j-1?'\n':' '); 
    printf("[t=%d] P:%c",t,e==0?'\n':' ');for(i=0;i<e;i++)printf("%f%c",p,   i==e-1?'\n':' '); 
    printf("[t=%d] U:%c",t,k==0?'\n':' ');for(i=0;i<k;i++)printf("%f%c",w[i],i==k-1?'\n':' '); 
    #endif

    if      (c<l) {size=j;} 
    else if (c>u) {tmp=a;a=w;w=tmp;ofs=u+1;size=k;}
    else break;
  }

  return p;
}

//void          init_genrand   (unsigned long s);
//double        genrand_real1  (void);
//
//int main(int argc, char **argv){
//
//  int i;
//  int     N = atoi(argv[1]);
//  double *a = malloc(N*sizeof(double));
//  double *w = malloc(N*sizeof(double));
//
//  for(i=0;i<N;i++) a[i]=genrand_real1();
//  printf("ANSWER: %f\n",median(a,w,N));
//
//  return 0;
//}
  
