#include<stdio.h>
#include<stdlib.h>

int main(int argc, char** argv){

  int s,k,K,t,end,*buf,l,m,n,L,M,N;
  FILE *fp  = fopen(argv[1],"rb");

  fseek(fp,0,SEEK_END);end=ftell(fp);fseek(fp,0,SEEK_SET); 
  buf = malloc((1+end)*sizeof(int)); 

  fread(buf,sizeof(int),8,fp); 
  printf("Resolution:   %3dx%3dx%2d\n",buf[0],buf[1],buf[2]);
  printf("Tracker size: %3dx%3dx%2d\n",buf[4],buf[5],buf[6]);
  printf("#Frames:      %3d\n",        buf[3]);
  printf("Cutoff:       %3d\n\n",      buf[7]);

  L=buf[0];M=buf[1];N=buf[2];

  while(ftell(fp)<end){
    fread(buf,sizeof(int),2,fp); t=buf[0];K=buf[1];
    printf("Frame: %.3d (#targets=%d)\n",t,K);

    fread(buf,sizeof(int),4*K,fp); 
    for(k=0;k<K;k++){l=buf[4*k];m=buf[4*k+1];n=buf[4*k+2];s=buf[4*k+3];
      if(argc==2) printf("%4d%4d%4d%4d\n",  l,  m,  n,s);
      if(argc> 2) printf("%4d%4d%4d%4d\n",M-m,l+1,N-n,s);
    } printf("\n");
  }

  return 0;
}
  

  
   
  



