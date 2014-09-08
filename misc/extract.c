#include<stdio.h>
#include<stdlib.h>
#include<assert.h>


int main(int argc, char** argv){

  unsigned char *buf; 
  FILE *fp  = fopen(argv[1],"rb");
  FILE *fpo = fopen(argv[4],"wb");
  int  t1   = atoi(argv[2]);
  int  t2   = atoi(argv[3]);
  int  s[4]; assert(t2>=t1);
  int  t,L;

  if(argc==1){printf("./extract <input> <frame begin> <frame end> <output>\n");exit(1);}
  fread (s,sizeof(int),4,fp ); L=s[0]*s[1]*s[2];s[3]=t2-t1+1;
  fwrite(s,sizeof(int),4,fpo); 
  
  fseek (fp,sizeof(int)*4+t1*L,SEEK_SET);buf=malloc(L);
  for(t=t1;t<=t2;t++){
    fread (buf,L,1,fp ); 
    fwrite(buf,L,1,fpo); 
  }

  fclose(fp );
  fclose(fpo);

  return 0;
}
  

  
   
  



