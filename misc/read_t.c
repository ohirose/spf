#include<stdio.h>
#include<stdlib.h>

int main(int argc, char** argv){

  int i,k,K,t,end,*buf;
  FILE *fp  = fopen(argv[1],"rb");

  fseek(fp,0,SEEK_END);end=ftell(fp);fseek(fp,0,SEEK_SET); 
  buf = malloc((1+end)*sizeof(int)); 

  fread(buf,sizeof(int),8,fp); 
  printf("Resolution:   %3dx%3dx%2d\n",buf[0],buf[1],buf[2]);
  printf("Tracker size: %3dx%3dx%2d\n",buf[4],buf[5],buf[6]);
  printf("#Frames:      %3d\n",        buf[3]);
  printf("Cutoff:       %3d\n\n",      buf[7]);

  while(ftell(fp)<end){
    fread(buf,sizeof(int),2,fp); t=buf[0];K=buf[1];
    printf("Frame: %.3d (#targets=%d)\n",t,K);

    fread(buf,sizeof(int),4*K,fp); 
    for(k=0;k<K;k++)for(i=0;i<4;i++) printf("%3d%c",buf[4*k+i],i==3?'\n':' ');
    printf("\n");
  }

  return 0;
}
  

  
   
  



