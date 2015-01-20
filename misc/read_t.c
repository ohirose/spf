#include<stdio.h>
#include<stdlib.h>

int main(int argc, char** argv){

  int k,K,t,end,buf[8],L,M,N,flag; float s,l,m,n,*xy;
  FILE *fp  = fopen(argv[1],"rb");

  flag= (argc>2)?1:0;
  fseek(fp,0,SEEK_END);end=ftell(fp);fseek(fp,0,SEEK_SET); 
  xy = malloc((1+end)*sizeof(float));

  fread(buf,sizeof(int),8,fp); 
  if(!flag) printf("Resolution:   %3dx%3dx%2d\n",buf[0],buf[1],buf[2]);
  if( flag) printf("Resolution:   %3dx%3dx%2d\n",buf[1],buf[0],buf[2]);
  printf("Tracker size: %3dx%3dx%2d\n",buf[4],buf[5],buf[6]);
  printf("#Frames:      %3d\n",        buf[3]);
  printf("Cutoff:       %3d\n\n",      buf[7]);

  L=buf[0];M=buf[1];N=buf[2];

  while(ftell(fp)<end){
    fread(buf,sizeof(int),2,fp); t=buf[0];K=buf[1];
    printf("Frame: %.3d (#targets=%d)\n",t,K);

    fread(xy,sizeof(float),4*K,fp);
    for(k=0;k<K;k++){l=xy[4*k];m=xy[4*k+1];n=xy[4*k+2];s=xy[4*k+3];
      if(!flag) printf("%.3f\t%.3f\t%.3f\t%.3f\n",  l,  m,  n,s);
      if( flag) printf("%.3f\t%.3f\t%.3f\t%.3f\n",M-m,l+1,N-n,s);
    } printf("\n");
  }

  return 0;
}
  

  
   
  



