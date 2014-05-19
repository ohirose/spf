#include<stdio.h>
#include<stdlib.h>
#include<assert.h>


int main(int argc, char** argv){

  unsigned char *buf; 
  FILE *fp  = fopen(argv[1],"rb");
  int  cut  = (argc<3)? 0:atoi(argv[2]);
  int  t1   = (argc<5)? 0:atoi(argv[3]);
  int  t2   = (argc<5)? 0:atoi(argv[4]);
  int  s[4],i,j,k,l,t,I,J,K,L,T;

  if(argc==1){printf("read_d <data> <cutoff (opt)> <frame begin (opt)> <frame end (opt)>\n");exit(1);}

  fread(s,sizeof(int),4,fp); I=s[0];J=s[1];K=s[2];T=s[3]; L=I*J*K;
  assert(t1<=t2 && t2<T); if(t2==0) t2=T-1;

  printf("Resolution:  %dx%3dx%2d\n",I,J,K);
  printf("Frames:      %d\n",        T    );
  printf("Frame Begin: %d\n",        t1   );
  printf("Frame End:   %d\n",        t2   );
  printf("Cutoff:      %d\n\n",      cut  );

  fseek (fp,4+t1*L,SEEK_SET);

  buf=malloc(L);
  for(t=t1;t<=t2;t++){fread(buf,L,1,fp); printf("Frame:%.4d\n",t);
    for(i=0;i<I;i++)for(j=0;j<J;j++)for(k=0;k<K;k++){l=i+j*I+k*I*J;
      if(buf[l]>cut)printf("signal: %3d\tvoxel: %d %d %d\n",buf[l],i,j,k);
    }
  }

  return 0;
}
  

  
   
  



