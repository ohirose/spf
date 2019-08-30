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
#include<string.h>
#include<ctype.h>
#include<math.h>
#include"util.h"
#include"io.h"

#define CONFIG_LINE_MAX 2048

void loadImageInfo (ImageHeader* h, const char *file){
  int  j,n,l,L,N=5; FILE *fpi; char *s,*b, line[CONFIG_LINE_MAX];
  char fni[MAX_FILENAME_LENGTH], fnb[MAX_FILENAME_LENGTH];
  int  nargs []={1,1,3,1,1};
  char *names[]={"rawimage","headersize","resolution","numberofframes","voxelaspectratio.z"};
  char *form []={"%s","%d","%d,%d,%d","%lf"};
  strcpy(fni,file);strcat(fni,".info");
  strcpy(fnb,file);strcat(fnb,".bin" ); strcpy(h->fnRawImage,fnb);

  fpi=fopen_e(fni,"r");
  while(fgets(line,CONFIG_LINE_MAX,fpi)){
    if(strchr("#|+\n",line[0])){continue;} L=strlen(line);
    for(l=0;l<L;l++) line[l]=tolower(line[l]);
    for(n=0;n<N;n++){s=names[n];if(!strstr(line,s))continue; b=line+strlen(s);
      switch(n){
        case  0: j=sscanf(b,form[0], h->fnRawImage);                                                  break;
        case  1: j=sscanf(b,form[1], &(h->rawHeaderSize));                                            break;
        case  2: j=sscanf(b,form[2], &(h->resolution.X),&(h->resolution.Y),&(h->resolution.Z));       break;
        case  3: j=sscanf(b,form[1], &(h->nFrames));                                                  break;
        case  4: j=sscanf(b,form[3], &(h->voxelAspectRatio.z));                                       break;
      } if(j!=nargs[n]) printf("Can't read the line: %s",line);
    }
  } fclose(fpi);
  return;
}

void countTracker (struct _tracker *tracker) {
  int T    = tracker->header.nFrames;
  int maxK = tracker->header.maxnTrackers;
  int t,k,D=3;
  for(t=0;t<T;t++) tracker->nums[t]=0;
  for(t=0;t<T;t++)for(k=0;k<maxK;k++)if(!isnan(tracker->data[3+k*(D+1)+t*maxK*(D+1)])) tracker->nums[t]++;
  return;
}

void calculateTrackerMeans (struct _tracker *tracker, const int nFrames, const int maxnTrackers) {
    int T    = nFrames;
    int maxK = maxnTrackers;
    int D = 3;
    int t,k;
    int cnt;
    double x,y,z;
    for ( k=0; k<maxK; k++ ) {
        cnt = 0;
        x = y = z = 0.0;
        for ( t=0; t<T; t++ ) {
            if (!isnan(tracker->data[3+k*(D+1)+t*maxK*(D+1)])) {
                x += tracker->data[0+k*(D+1)+t*maxK*(D+1)];
                y += tracker->data[1+k*(D+1)+t*maxK*(D+1)];
                z += tracker->data[2+k*(D+1)+t*maxK*(D+1)];
                cnt++;
            }
        }
        tracker->means[0+k*D] = x/cnt;
        tracker->means[1+k*D] = y/cnt;
        tracker->means[2+k*D] = z/cnt;
    }
}

void loadObjectNames( struct _tracker *tracker, const int maxnTrackers, const char *filename ){
    int  k,K,lim=256;
    FILE *fp;

    if ((fp=fopen( filename, "r" ))==NULL) {
        printf("File: \"%s\" Not found.\n",filename);
        tracker->names = NULL;
        return;
    }

    K = maxnTrackers;
    tracker->names = (char**)malloc(K*sizeof(char*));
    for ( k=0; k<K; k++ ) {
        tracker->names[k] = (char*)malloc(lim*sizeof(char));
        fscanf(fp,"%s\n",tracker->names[k]);
    }

    fclose(fp);

    return;
}

void loadTrackerInfo(TrackerHeader* h, const char *file){
  int j,n,l,L,N=5; FILE *fp; char *s,*b, line[CONFIG_LINE_MAX];
  int  nargs []={1,1,1,1,1};
  char *names[]={"trackerdata","trackernames","headersize","maxnumberoftrackers","numberofframes" };
  char *form []={"%s","%d"};

  fp=fopen_e((char*)file,"r"); h->fnTrackerRawData[0]='\0';
  while(fgets(line,CONFIG_LINE_MAX,fp)){
    if(strchr("#|+\n",line[0])){continue;} L=strlen(line);
    for(l=0;l<L;l++) line[l]=tolower(line[l]);
    for(n=0;n<N;n++){s=names[n];if(!strstr(line,s))continue; b=line+strlen(s);
      switch(n){
        case  0: j=sscanf(b,form[0], h->fnTrackerRawData);  break;
        case  1: j=sscanf(b,form[0], h->fnTrackerNames);    break;
        case  2: j=sscanf(b,form[1], &(h->rawHeaderSize));  break;
        case  3: j=sscanf(b,form[1], &(h->maxnTrackers));   break;
        case  4: j=sscanf(b,form[1], &(h->nFrames));        break;
      } if(j!=nargs[n]) fprintf(stderr,"Can't read the line: %s",line);
    }
  } fclose(fp);

  return;
}

void loadTracker (struct _tracker *tracker, const char *file){
  int T,maxK,num; FILE *fp;
  char fni[MAX_FILENAME_LENGTH], fnb[MAX_FILENAME_LENGTH];

  strcpy(fni,file);strcat(fni,".info");
  strcpy(fnb,file);strcat(fnb,".bin" );

  loadTrackerInfo(&(tracker->header),fni);
  if(!strlen(tracker->header.fnTrackerRawData))
    strcpy(tracker->header.fnTrackerRawData,fnb);

  fp=fopen_e(tracker->header.fnTrackerRawData,"rb");
  fseek(fp,tracker->header.rawHeaderSize,SEEK_SET);

  T    = tracker->header.nFrames;
  maxK = tracker->header.maxnTrackers;

  /* allocate tracker */
  tracker->means = (double*) malloc(  maxK*(sizeof(double))*3);
  tracker->data  = (double*) malloc(T*maxK*(sizeof(double))*4);
  tracker->nums  = (int*)    malloc(T*(sizeof(int)));

  num=fread(tracker->data,sizeof(double),T*maxK*4,fp);
  if(num!=T*maxK*4){fprintf(stderr,"Failed to load %s.(%d,%s)\n",fnb,__LINE__,__FILE__);
    printf("num=%d,T=%d,maxK=%d,T*maxK=%d,T*maxK*4=%d,nbyte=%d\n",num,T,maxK,T*maxK,T*maxK*4,num*8);exit(EXIT_FAILURE);}
  fclose(fp);

  calculateTrackerMeans(tracker,T,maxK);
  countTracker    (tracker);
  loadObjectNames ( tracker, maxK, tracker->header.fnTrackerNames );

  return;
}

