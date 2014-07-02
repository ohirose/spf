#include<stdio.h>
#include<stdlib.h>
#include<assert.h>


int main(int argc, char** argv){

  unsigned char *buf; 
  long int L;
  int  s[4],t,T;
  FILE *fpi  = fopen(argv[1],"rb");
  FILE *fpo  = fopen(argv[2],"wb");
  int    tu  = atoi (argv[3]);

  if(argc<=3){printf("./uturn <input> <output> <until>\n");exit(1);}

  fread(s,sizeof(int),4,fpi); L=s[0]*s[1]*s[2];T=s[3];

  printf("Resolution:  %dx%3dx%2d\n",s[0],s[1],s[2]);
  printf("Frames:      %d\n",        T    );

  buf=malloc(L); s[3]=2*tu-1;

  fwrite(s,sizeof(int),4,fpo);
  for(t=0;t<tu;t++){
    if((t+1)%10==0) printf("t=%d\n",t+1);
    fread (buf,1,L,fpi); 
    fwrite(buf,1,L,fpo);
  }

  for(t=tu-2;t>=0;t--){
    if((t+1)%10==0) printf("t=%d\n",t+1);
    fseek (fpi,t*L+(4*sizeof(int)),SEEK_SET); 
    fread (buf,1,L,fpi); 
    fwrite(buf,1,L,fpo);
  }
  
  fclose(fpi);
  fclose(fpo);

  return 0;
}
  

  
   
  



