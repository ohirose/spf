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
#include"../common/util.h"
#define LIMIT 1024


int mstree(int          ** G,  /* OUTPUT: K x K (MST; Array of neighbor indeces) */ 
           const double ** D,  /* INPUT:  K x K (Distance matrix)                */
           const int       K   /* INPUT:        (#nodes)                         */){

  int     u,v,a,l,e,ft,ct=0;
  int     nu,nv,ne;
  int     *U,*V,*f; 
  double  *W;
  sortbox *sb;

  ne = (K*(K-1))/2;
  f  = malloc(K  * sizeof(int    ));
  U  = malloc(ne * sizeof(int    ));
  V  = malloc(ne * sizeof(int    ));
  W  = malloc(ne * sizeof(double ));
  sb = malloc(ne * sizeof(sortbox));

  /* Initialization */
  for(u=0;u<K;u++){G[u][K-1]=0;f[u]=u;}

  /* Sorting edges by weights */ 
  for(v=0;v<K;v++)for(u=0;u<v;u++){e=u+v*(v-1)/2; W[e]=D[u][v]; U[e]=u; V[e]=v; assert(e<ne);}
  prepare_sortbox(sb,W,ne); qsort(sb,ne,sizeof(sortbox),cmp_sortbox);

  /* Selecting edges */
  for(l=0;l<ne;l++){e=sb[l].idx;u=U[e];v=V[e];
    if(f[u]!=f[v]){nu=G[u][K-1];nv=G[v][K-1];ct++;
      G[u][nu ]=v; G[v][nv ]=u;
      G[u][K-1]++; G[v][K-1]++;
      ft=f[v]; for(a=0;a<K;a++)if(f[a]==ft)f[a]=f[u];
    }
    if(ct==K-1)break;
  }

  free(f);free(U);free(V);free(W);free(sb);

  return 0;
}  

int makeiter(int *iter, int *u, const int ** G, const int K, const int root){
  int k,c,n,nn,S[LIMIT];
  int i=0,top=-1; 

  for(k=0;k<K;k++){u[k]=-1;} top++;S[top]=root;
  while(top>=0){k=S[top];top--;iter[i]=k;nn=G[k][K-1];i++;
    for(n=0;n<nn;n++){c=G[k][n];
      if(c!=u[k]){u[c]=k;top++;S[top]=c;}
    }
  }

  return 0;
} 



