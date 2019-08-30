// Copyright (c) 2014-2019 Shotaro Kagaguchi, Masayuki Tamura, and Osamu Hirose
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
#include<time.h>
#include<GL/freeglut.h>
#include"../common/util.h"
#include"../common/io.h"

#define CONFIG_LINE_MAX 128
#define INTENSITY_X_MAX_SIZE 2048
#define CONFIGURE_VERSION 0
#define SEGMENTATION_COLOR_MAX 255

#define MOUSE_BUTTON_CNT 12

#if !defined(__HEADER_MOUSE_BUTTON_ADDITION)
#define __HEADER_MOUSE_BUTTON_ADDITION
#define OS_CSEKU_LEFT_BUTTON GLUT_LEFT_BUTTON
#define OS_CSEKU_RIGHT_BUTTON GLUT_RIGHT_BUTTON
#define OS_CSEKU_MIDDLE_BUTTON GLUT_MIDDLE_BUTTON
#define OS_CSEKU_WHEEL_UP 3
#define OS_CSEKU_WHEEL_DOWN 4
#define OS_CSEKU_BACK 7
#define OS_CSEKU_NEXT 8
#endif


#if !defined(__HEADER_BOOL)
#define __HEADER_BOOL
typedef int bool;
const bool true = 1;
const bool false = 0;
#endif

const int success = 0;
const int failure = 1;

typedef struct RGB {
    double r,g,b;
} RGB;

const RGB Red     = { 1.0, 0.0, 0.0 }; const RGB Fuchsia = { 0.5, 0.0, 0.5 };
const RGB Green   = { 0.0, 0.5, 0.0 }; const RGB Yellow  = { 1.0, 1.0, 0.0 };
const RGB Blue    = { 0.0, 0.0, 1.0 }; const RGB Lime    = { 0.0, 1.0, 0.0 };
const RGB Aqua    = { 0.0, 1.0, 1.0 }; const RGB Olive   = { 0.5, 0.5, 0.0 };
const RGB Teal    = { 0.0, 0.5, 0.5 }; const RGB White   = { 1.0, 1.0, 1.0 };
const RGB Navy    = { 0.0, 0.0, 0.5 }; const RGB Gray    = { 0.4, 0.4, 0.4 };
const RGB Purple  = { 0.5, 0.0, 0.5 }; const RGB Black   = { 0.0, 0.0, 0.0 };
const RGB Maroon  = { 0.5, 0.0, 0.0 };

typedef struct Quaternion {
    double w,x,y,z;
} Quaternion;

struct {
    int frame;
    unsigned char cutoff;
    int trackerIdx;
    struct { int x; int y; int z; } rotation;
    int zoomPercentage;
    struct {
        struct {
            bool intensity;
            struct { bool box; bool point; } tracker;
            bool segmentation;
            bool names;
        } plot;
        bool move;
        bool resize;
        bool focus;
        bool parts;
        struct { bool x; bool y; bool z; } rotation;
    } flag;
    double rotate[16];
    Quaternion current;
    struct {
        struct {
            bool button[MOUSE_BUTTON_CNT];
            struct { int x,y; } position;
        } mouse;
    } device;
    struct { int x,y; } translation;
    struct {
        bool enable;
        struct { int x,y,z; } position;
    } cursor;
} state;

struct _parameters{
    struct {
        struct { int w,h; } size;
        struct { int x,y; } position;
    } window;
    int plotPointSize;
    RGB backgroundColor;
    RGB foregroundColor;
    RGB cursorColor;
    struct { double x,y,z; } viewpoint;
    struct { int colorCnt; RGB color[SEGMENTATION_COLOR_MAX]; } segmentation;
    struct { double x,y; } mouseSensitivity;
    struct { double x,y,z; } trackerBoxSize;
    int intervalTime_ms;
    int trackerColorType;
    char scrShotFileName[CONFIG_LINE_MAX];
    int  defaultCutoff;
} parameters;

struct {
    struct {
        ImageHeader header;
        FILE *fp;
    } intensity;
    struct _tracker tracker;
    struct {
        ImageHeader header;
        FILE *fp;
    } segmentation;
} data;


void Q2R ( double *r, Quaternion q ) {
    double xx,yy,zz,xy,yz,zx,xw,yw,zw;
    xx = q.x * q.x * 2.0; yy = q.y * q.y * 2.0; zz = q.z * q.z * 2.0;
    xy = q.x * q.y * 2.0; yz = q.y * q.z * 2.0; zx = q.z * q.x * 2.0;
    xw = q.x * q.w * 2.0; yw = q.y * q.w * 2.0; zw = q.z * q.w * 2.0;

    r[ 0]=1.0-yy-zz; r[ 1]=xy+zw;     r[ 2]=zx-yw;     r[ 3]=0.0;
    r[ 4]=xy-zw;     r[ 5]=1.0-zz-xx; r[ 6]=yz+xw;     r[ 7]=0.0;
    r[ 8]=zx+yw;     r[ 9]=yz-xw;     r[10]=1.0-xx-yy; r[11]=0.0;
    r[12]=0.0;       r[13]=0.0;       r[14]=0.0;       r[15]=1.0;
}

void QuaternionMult ( Quaternion *q1, Quaternion q2 ) {
    Quaternion q0 = {
        q1->w * q2.w - q1->x * q2.x - q1->y * q2.y - q1->z * q2.z,
        q1->w * q2.x + q1->x * q2.w + q1->y * q2.z - q1->z * q2.y,
        q1->w * q2.y - q1->x * q2.z + q1->y * q2.w + q1->z * q2.x,
        q1->w * q2.z + q1->x * q2.y - q1->y * q2.x + q1->z * q2.w
    };
    q1->w=q0.w;
    q1->x=q0.x;
    q1->y=q0.y;
    q1->z=q0.z;
}

void resetRotate ( void ) {
    state.current.w = 1.0;
    state.current.x = 0.0;
    state.current.y = 0.0;
    state.current.z = 0.0;
    Q2R(state.rotate,state.current);
}

void resetTranslation ( void ) {
    state.translation.x = 0;
    state.translation.y = 0;
}

void initializeState ( void ) {
    int i;
    state.frame = 0;
    state.zoomPercentage = 100;
    state.flag.plot.intensity = true;
    state.flag.plot.tracker.box  = false;
    state.flag.plot.tracker.point= false;
    state.flag.plot.segmentation = false;
    state.flag.plot.names = true;
    state.flag.move = true;
    state.flag.resize = false;
    state.flag.focus = false;
    state.cutoff = 60;
    state.trackerIdx = 0;
    resetRotate();
    resetTranslation();
    for ( i=0; i<MOUSE_BUTTON_CNT; i++ ) {
        state.device.mouse.button[i] = false;
    }
    state.cursor.enable = false;
    state.cursor.position.x = 0;
    state.cursor.position.y = 0;
    state.cursor.position.z = 0;
}

void initializeParameters ( void ) {
    parameters.window.size.w = 1300;
    parameters.window.size.h = 800;
    parameters.window.position.x = 0;
    parameters.window.position.y = 0;
    parameters.plotPointSize = 2;
    parameters.backgroundColor = White;
    parameters.foregroundColor = Black;
    parameters.cursorColor = Blue;
    parameters.mouseSensitivity.x = 1.0;
    parameters.mouseSensitivity.y = 1.0;
    parameters.trackerBoxSize.x = 8.0;
    parameters.trackerBoxSize.y = 8.0;
    parameters.trackerBoxSize.z = 8.0;
    parameters.segmentation.colorCnt = 6;
    parameters.segmentation.color[ 0] = Red;
    parameters.segmentation.color[ 1] = Green;
    parameters.segmentation.color[ 2] = Blue;
    parameters.segmentation.color[ 3] = Aqua;
    parameters.segmentation.color[ 4] = Fuchsia;
    parameters.segmentation.color[ 5] = Lime;
    parameters.segmentation.color[ 6] = Olive;
    parameters.segmentation.color[ 7] = Teal;
    parameters.segmentation.color[ 8] = Navy;
    parameters.segmentation.color[ 9] = Purple;
    parameters.segmentation.color[10] = Maroon;
    parameters.segmentation.color[11] = Yellow;
    parameters.intervalTime_ms = 0;
    parameters.trackerColorType = 0;
    strcpy(parameters.scrShotFileName,"screenshot.rgb");
}

void string2lower ( char *str ) {
    int i=0;
    while ( str[i] != '\0' ) {
        str[i] = tolower(str[i]);
        i++;
    }
}

int loadConfig(char *file){
  int j,n,l,L,N=14; FILE *fp; char *s,*b, line[CONFIG_LINE_MAX]; struct _parameters *p=&parameters;
  int  nargs[14]={2,2,1,3, 3,3,3,2, 2,1,1,1, 1,1};
  char *form []={"%s","%d","%d,%d","%d,%d,%d","%lf","%lf,%lf","%lf,%lf,%lf"};
  char *names[]={
     "window.size","window.position","plotpointsize","backgroundcolor","foregroundcolor",
     "cursorcolor","viewpoint","mousesensitivity","trackerboxsize","trackercolortype",
     "segmentation.colorcnt","intervaltime_ms","scrshotfilename","defaultcutoff"
   };

  fp=fopen(file,"r");if(!fp){printf("File not found: %s\n",file);exit(EXIT_FAILURE);}
  while(fgets(line,CONFIG_LINE_MAX,fp)){
    if(strchr("#|+\n",line[0])) continue; L=strlen(line);
    for(l=0;l<L;l++) line[l]=tolower(line[l]);
    for(n=0;n<N;n++){s=names[n];if(!strstr(line,s))continue; b=line+strlen(s);
      switch(n){
        case  0: j=sscanf(b,form[2], &(p->window.size.w),&(p->window.size.h));                                  break;
        case  1: j=sscanf(b,form[2], &(p->window.position.x),&(p->window.position.y));                          break;
        case  2: j=sscanf(b,form[1], &(p->plotPointSize));                                                      break;
        case  3: j=sscanf(b,form[6], &(p->backgroundColor.r),&(p->backgroundColor.g),&(p->backgroundColor.b));  break;
        case  4: j=sscanf(b,form[6], &(p->foregroundColor.r),&(p->foregroundColor.g),&(p->foregroundColor.b));  break;
        case  5: j=sscanf(b,form[6], &(p->cursorColor.r),&(p->cursorColor.g),&(p->cursorColor.b));              break;
        case  6: j=sscanf(b,form[6], &(p->viewpoint.x),&(p->viewpoint.y),&(p->viewpoint.z));                    break;
        case  7: j=sscanf(b,form[5], &(p->mouseSensitivity.x), &(p->mouseSensitivity.y));                       break;
        case  8: j=sscanf(b,form[5], &(p->trackerBoxSize.x),&(p->trackerBoxSize.y));                            break;
        case  9: j=sscanf(b,form[1], &(p->trackerColorType));                                                   break;
        case 10: j=sscanf(b,form[1], &(p->segmentation.colorCnt));                                              break;
        case 11: j=sscanf(b,form[1], &(p->intervalTime_ms));                                                    break;
        case 12: j=sscanf(b,form[0], &(p->scrShotFileName));                                                    break;
        case 13: j=sscanf(b,form[1], &(p->defaultCutoff));                                                      break;
      } if(j!=nargs[n]) printf("Can't read the line: %s",line);
    }
  } fclose(fp);

  return success;
}

int checkArguments ( int argc, char *argv[] ) {
    if ( argc < 2 ) {
        fprintf(stderr,"Too few arguments.(%d,%s)\n",__LINE__,__FILE__);
        fputs("1st argument : configure file\n",stderr);
        fputs("2nd argument : intensity file\n",stderr);
        fputs("3rd argument : tracker file\n",stderr);
        return failure;
    }
    if ( loadConfig ( argv[1] ) == failure ) {
        return failure;
    }

    loadImageInfo ( &(data.intensity.header),argv[2]);
    data.intensity.fp = fopen_e(data.intensity.header.fnRawImage,"rb");
    state.cutoff = parameters.defaultCutoff;
    data.tracker.data = NULL;

    if ( argc > 3 ) {
        loadTracker ( &(data.tracker), argv[3]);
        state.flag.plot.tracker.box   = false;
        state.flag.plot.tracker.point = true;
    }
    data.segmentation.fp = NULL;
    if ( argc > 4 ) {
        loadImageInfo (  &(data.segmentation.header),argv[4]);
        state.flag.plot.segmentation = true;
        data.segmentation.fp=fopen_e(data.segmentation.header.fnRawImage,"rb");
    }
    return success;
}

void setRenderingColor ( RGB rgb ) {
    glColor3f ( rgb.r, rgb.g, rgb.b );
}

void clearBackgroundColor ( RGB rgb ) {
    glClearColor(rgb.r, rgb.g, rgb.b, 0);
    glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
}

void renderString ( double x, double y, const char *str ) {
    double z = -1.0;
    int i;
    glRasterPos3f ( x, y, z );
    for ( i=0; str[i] != '\0'; i++ ) {
        glutBitmapCharacter ( GLUT_BITMAP_HELVETICA_18, str[i] );
    }
}

void renderString3 ( double x, double y, double z, const char *str ) {
    int i;
    glRasterPos3f ( x, y, z );
    for ( i=0; str[i] != '\0'; i++ ) {
        glutBitmapCharacter ( GLUT_BITMAP_HELVETICA_12, str[i] );
    }
}

void renderFrameIndex ( RGB rgb) {
    char str[256];
    setRenderingColor(rgb);
    sprintf(str,"frame : %04d / %04d" , state.frame+1,
            data.intensity.header.nFrames);
    renderString ( -0.25 , -0.19, str );
    renderString ( -0.0, -0.19, "s : " );
    if ( state.flag.move ) {
        setRenderingColor(Red);
    } else {
        setRenderingColor(rgb);
    }
    renderString ( 0.016, -0.19, "start" );
    setRenderingColor(rgb);
    renderString ( 0.040, -0.19, " / " );
    if ( !state.flag.move ) {
        setRenderingColor(Red);
    } else {
        setRenderingColor(rgb);
    }
    renderString ( 0.050, -0.19, "stop" );
}

void renderFocus ( RGB rgb) {
    setRenderingColor(rgb);
    renderString( -0.0, -0.21, "f : " );
    setRenderingColor(Red);
    if ( state.flag.focus ) {
        setRenderingColor(Red);
    } else {
        setRenderingColor(rgb);
    }
    renderString( 0.016, -0.21, "focus mode" );
}

void renderCutoff ( RGB rgb) {
    char str[256];
    sprintf(str, "cutoff : %03d / 255", state.cutoff );
    setRenderingColor(rgb);
    renderString( -0.25, -0.21, str );
}

void renderParts( RGB rgb) {
    char str[256];
    setRenderingColor(rgb);
    if ( state.flag.parts ) {
        sprintf(str,"plot : %02d / %02d" , state.trackerIdx + 1,
                data.tracker.header.maxnTrackers);
        renderString( -0.25, -0.23, str );
    } else {
        renderString( -0.25, -0.23, "plot : all" );
    }
}

void renderCursor ( RGB rgb ) {
    char str[256];
    setRenderingColor(parameters.foregroundColor);
    if ( state.cursor.enable ) {
        sprintf(str,"Cursor : %d,%d,%d",state.cursor.position.x,
            state.cursor.position.y,state.cursor.position.z);
        renderString( -0.0, -0.23, str );
    } else {
        renderString( -0.0, -0.23, "cursor : disable" );
    }
}

void drawCursor ( int x, int y, int z ) {
    int cx,cy,cz;
    cx = state.cursor.position.x;
    cy = state.cursor.position.y;
    cz = state.cursor.position.z;
    setRenderingColor(parameters.cursorColor);
    glBegin(GL_LINES);
    glVertex3d(cx,cy,0<cz?0:cz);
    glVertex3d(cx,cy,z>cz?z:cz);
    glVertex3d(cx,0<cy?0:cy,cz);
    glVertex3d(cx,y>cy?y:cy,cz);
    glVertex3d(0<cx?0:cx,cy,cz);
    glVertex3d(x>cx?x:cx,cy,cz);
    glEnd();
}

void drawVoxel ( int mode ) {
  int x,y,z,n,idx,mark; double scale; FILE *fp;
  unsigned char  value[INTENSITY_X_MAX_SIZE];
  unsigned long int X = data.intensity.header.resolution.X;
  unsigned long int Y = data.intensity.header.resolution.Y;
  unsigned long int Z = data.intensity.header.resolution.Z;
  unsigned long int hs= mode?data.segmentation.header.rawHeaderSize:data.intensity.header.rawHeaderSize;
  unsigned long int offset = hs + X*Y*Z*(state.frame);

  fp  =mode?(data.segmentation.fp):(data.intensity.fp); fseek(fp,offset,SEEK_SET);
  mark=mode?0:state.cutoff;

  for(z=0;z<Z;z++)for(y=0;y<Y;y++){n=fread(value,1,X,fp);
    if(n!=X){fprintf(stderr,"Failed to load.(%d,%s)\n",__LINE__,__FILE__);}
    for(x=0;x<X;x++){if(value[x]<=mark) continue;
      switch(mode){
        case 0: /* intensity */
          scale=value[x]>=150.0?1.0:value[x]/150.0;
          glColor3d ( 0.0 , scale, 0.0 ); break;
        case 1: /* segmentation */
          idx=value[x]%parameters.segmentation.colorCnt;
          setRenderingColor(parameters.segmentation.color[idx]); break;
      }
      if    (state.flag.resize){
        glPointSize(parameters.plotPointSize);
        glBegin(GL_POINTS);glVertex3d( x, y, z*data.intensity.header.voxelAspectRatio.z ); glEnd();
      }
      else/*!state.flag.resize*/ {
        glPointSize(parameters.plotPointSize);
        glBegin(GL_POINTS);glVertex3d( x, y, z ); glEnd();
      }
    }
  }
}

void drawBox ( int x, int y, int z, int X, int Y, int Z ) {
    X+=x;Y+=y;Z+=z;
    glBegin(GL_LINES );
    glVertex3d(x,y,z); glVertex3d(X,y,z); glVertex3d(x,y,z); glVertex3d(x,Y,z);
    glVertex3d(x,y,z); glVertex3d(x,y,Z); glVertex3d(X,Y,Z); glVertex3d(X,Y,z);
    glVertex3d(X,Y,Z); glVertex3d(X,y,Z); glVertex3d(X,Y,Z); glVertex3d(x,Y,Z);
    glVertex3d(X,Y,z); glVertex3d(x,Y,z); glVertex3d(X,Y,z); glVertex3d(X,y,z);
    glVertex3d(x,Y,Z); glVertex3d(x,Y,z); glVertex3d(x,Y,Z); glVertex3d(x,y,Z);
    glVertex3d(X,y,Z); glVertex3d(x,y,Z); glVertex3d(X,y,Z); glVertex3d(X,y,z);
    glEnd();
}

void drawTracker ( void ) {
    int i,D=3;
    int maxK = data.tracker.header.maxnTrackers;
    int beginIdx,endIdx;
    if ( state.flag.parts ||state.flag.focus) {
        beginIdx = state.trackerIdx;
        endIdx = state.trackerIdx + 1;
    } else {
        beginIdx = 0;
        endIdx = maxK;
    }
    for ( i=beginIdx; i<endIdx; i++ ) {
        float x,y,z,X,Y,Z;
        int t = state.frame;
        int idx;
        if ( !(data.tracker.data[3+i*(D+1)+t*(D+1)*maxK]<= 0) &&
             !(data.tracker.data[3+i*(D+1)+t*(D+1)*maxK]>= 0) ) {
            continue;
        }
        x = data.tracker.data[0+i*(D+1)+t*(D+1)*maxK] - parameters.trackerBoxSize.x/2.0;
        y = data.tracker.data[1+i*(D+1)+t*(D+1)*maxK] - parameters.trackerBoxSize.y/2.0;
        z = data.tracker.data[2+i*(D+1)+t*(D+1)*maxK] - parameters.trackerBoxSize.z/2.0;
        X = parameters.trackerBoxSize.x;
        Y = parameters.trackerBoxSize.y;
        Z = parameters.trackerBoxSize.z;

        if ( state.flag.resize ) {
            z *= data.intensity.header.voxelAspectRatio.z;
            Z *= data.intensity.header.voxelAspectRatio.z;
        }
        if ( parameters.trackerColorType == 0 ) {
            idx=i%parameters.segmentation.colorCnt;
            setRenderingColor ( parameters.segmentation.color[idx] );
        } else if ( parameters.trackerColorType == 1 ) {
            if ( data.tracker.data[3+i*(D+1)+t*(D+1)*maxK] < state.cutoff ) {
                setRenderingColor(Gray);
            } else {
                double scale = (255-data.tracker.data[3+i*(D+1)+t*(D+1)*maxK]) /255.0;
                glColor3d ( 1.0, scale, 0.0 );
            }
        } else {
            idx=i%parameters.segmentation.colorCnt;
            setRenderingColor ( parameters.segmentation.color[idx] );
        }
        if ( state.flag.plot.tracker.box ) {
            drawBox ( (int)x, (int)y, (int)z, (int)X, (int)Y, (int)Z );
        }
        if ( state.flag.plot.tracker.point ) {
            glPointSize (parameters.plotPointSize+3.0 );
            glDisable   ( GL_DEPTH_TEST );
            glBegin     ( GL_POINTS );
            glVertex3d  ( x+X/2.0, y+Y/2.0, z+Z/2.0 );
            glEnd();
            glEnable    ( GL_DEPTH_TEST );
        }

        if(state.cursor.enable && (state.flag.parts||state.flag.focus)) {
            state.cursor.position.x = x+X/2.0;
            state.cursor.position.y = y+Y/2.0;
            state.cursor.position.z = z+Z/2.0;
        }

        if( state.flag.plot.names ){
            setRenderingColor(parameters.foregroundColor);
            glDisable(GL_DEPTH_TEST);
            if(data.tracker.names!=NULL) renderString3 ( x+X/2.0, y+Y/2.0, z+Z/2.0, data.tracker.names[i] );
            glEnable (GL_DEPTH_TEST);
        }
    }
}

void gotoNextFrame ( void ) {
    state.frame = (state.frame+1) % data.intensity.header.nFrames;
}
void gotoPreviousFrame ( void ) {
    state.frame--;
    if ( state.frame < 0 ) {
        state.frame += data.intensity.header.nFrames;
    }
}

void DISPLAY ( void ) {
    double x,y,z; int D=3;
    struct { double x,y,z; } basepoint;
    clearBackgroundColor(parameters.backgroundColor);
    glViewport ( 0, 0, parameters.window.size.w, parameters.window.size.h );
    glMatrixMode( GL_MODELVIEW ); glLoadIdentity();
    if ( state.flag.focus && state.flag.parts ) {
        basepoint.x = data.tracker.means[0+state.trackerIdx*D];
        basepoint.y = data.tracker.means[1+state.trackerIdx*D];
        basepoint.z = data.tracker.means[2+state.trackerIdx*D];
    } else {
        basepoint.x = basepoint.y = basepoint.z = 0.0;
    }
    if ( state.flag.resize ) {
        basepoint.z *= data.intensity.header.voxelAspectRatio.z;
    }
    gluLookAt ( 0, 0, parameters.viewpoint.x * state.zoomPercentage/100.0 , 0,0,0, 0,1,0);
    glMultMatrixd ( state.rotate );
    if ( state.flag.rotation.x ) { state.rotation.x = (state.rotation.x+1)%360; }
    if ( state.flag.rotation.y ) { state.rotation.y = (state.rotation.y+1)%360; }
    if ( state.flag.rotation.z ) { state.rotation.z = (state.rotation.z+1)%360; }
    glRotated ( (double) state.rotation.x , 1, 0, 0 );
    glRotated ( (double) state.rotation.y , 0, 1, 0 );
    glRotated ( (double) state.rotation.z , 0, 0, 1 );
    x = data.intensity.header.resolution.X;
    y = data.intensity.header.resolution.Y;
    z = data.intensity.header.resolution.Z;
    if ( state.flag.resize ) {
        z *= data.intensity.header.voxelAspectRatio.z;
    }
    if ( state.flag.focus && state.flag.parts ) {
        glTranslated( -basepoint.x, -basepoint.y, -basepoint.z );
    } else {
        glTranslated ( -x/2.0, -y/2.0, -z/2.0 );
    }
    glTranslated(state.translation.x,state.translation.y,0);
    glMatrixMode ( GL_PROJECTION ); glLoadIdentity();
    gluPerspective ( 30.0, (double)parameters.window.size.w/(double) parameters.window.size.h, 0.01, 10000.0 );
    setRenderingColor(Black);
    drawBox(0,0,0,x,y,z);

    if(state.cursor.enable)          drawCursor (x,y,z);
    if(state.flag.plot.intensity)    drawVoxel  (0);
    if(state.flag.plot.segmentation) drawVoxel  (1);
    if(data.tracker.data)            drawTracker();
    if(state.flag.move)              gotoNextFrame();

    glMatrixMode( GL_MODELVIEW ); glLoadIdentity();
    renderFrameIndex (parameters.foregroundColor);
    renderFocus      (parameters.foregroundColor);
    renderCutoff     (parameters.foregroundColor);
    renderParts      (parameters.foregroundColor);
    renderCursor     (parameters.foregroundColor);

    glutSwapBuffers();
}

void RESHAPE ( int w, int h ) {
    parameters.window.size.w = w;
    parameters.window.size.h = h;
}

void zoomIn ( void ) {
    state.zoomPercentage -= 2;
    if ( state.zoomPercentage < 0 ) {
        state.zoomPercentage = 0;
    }
}

void zoomOut ( void ) {
    state.zoomPercentage += 2;
    if ( state.zoomPercentage > 200 ) {
        state.zoomPercentage = 200;
    }
}

void MOUSE ( int button, int state_l, int x, int y ) {
    if ( button < MOUSE_BUTTON_CNT ) {
        if ( state_l == GLUT_DOWN ) {
            state.device.mouse.button[button] = true;
        } else if ( state_l == GLUT_UP ) {
            state.device.mouse.button[button] = false;
        }
    }
    switch ( button ) {
        case OS_CSEKU_LEFT_BUTTON  : break;
        case OS_CSEKU_RIGHT_BUTTON : break;
        case OS_CSEKU_MIDDLE_BUTTON :
            if ( state_l == GLUT_DOWN ) {
                resetRotate();
                resetTranslation();
            }
            break;
        case OS_CSEKU_WHEEL_UP :
            if ( state_l == GLUT_DOWN ) {
                zoomIn();
            }
            break;
        case OS_CSEKU_WHEEL_DOWN :
            if ( state_l == GLUT_DOWN ) {
                zoomOut();
            }
            break;
        case OS_CSEKU_BACK :
            if ( state_l == GLUT_DOWN ) {
                gotoPreviousFrame();
            }
            break;
        case OS_CSEKU_NEXT :
            if ( state_l == GLUT_DOWN ) {
                gotoNextFrame();
            }
            break;
        default : break;
    }
    state.device.mouse.position.x = x;
    state.device.mouse.position.y = y;
}

void MOTION ( int x, int y ) {
    int X,Y;
    double dx,dy;
    double length;
    X = state.device.mouse.position.x;
    Y = state.device.mouse.position.y;
    if ( state.device.mouse.button[OS_CSEKU_RIGHT_BUTTON] == true ) {
        dx = (x-X) * parameters.mouseSensitivity.x;
        dy = (y-Y) * parameters.mouseSensitivity.y;
        state.translation.x += dx;
        state.translation.y -= dy;
    }
    if ( state.device.mouse.button[OS_CSEKU_LEFT_BUTTON] == true ) {
        dx = (x-X) * parameters.mouseSensitivity.x/parameters.window.size.w;
        dy = (y-Y) * parameters.mouseSensitivity.y/parameters.window.size.h;
        length = sqrt ( dx*dx + dy*dy );
        if ( length != 0.0 ) {
            double pi = atan(1.0) * 4;
            double radian = length * pi;
            double theta = sin ( radian ) / length;
            Quaternion after = { cos(radian) , dy*theta, dx*theta, 0.0 };
            QuaternionMult (&(state.current),after);
            Q2R ( state.rotate,state.current);
        }
    }
    state.device.mouse.position.x = x;
    state.device.mouse.position.y = y;
}

void PASSIVE_MOTION ( int x, int y ) {
    state.device.mouse.position.x = x;
    state.device.mouse.position.y = y;
}

void toggleFlag ( bool *flag ) {
    *flag = !(*flag);
}

void trackerDisplayMode ( void ){
    bool b = state.flag.plot.tracker.box;
    bool p = state.flag.plot.tracker.point;
    int  m = 2*b+p;
    switch (m) {
        case 0: b=false; p=true;  break;
        case 1: b=true;  p=false; break;
        case 2: b=true;  p=true;  break;
        case 3: b=false; p=false; break;
    }
    state.flag.plot.tracker.box   = b;
    state.flag.plot.tracker.point = p;
}

void incrementTrackerIdx ( void ) {
    state.trackerIdx ++;
    if ( state.trackerIdx >= data.tracker.header.maxnTrackers ) {
        state.trackerIdx -= data.tracker.header.maxnTrackers;
    }
}

void decrementTrackerIdx ( void ) {
    state.trackerIdx --;
    if ( state.trackerIdx < 0 ) {
        state.trackerIdx += data.tracker.header.maxnTrackers;
    }
}

void incrementTrackerIdx50 ( void ) {
    state.trackerIdx += 50;
    if ( state.trackerIdx > data.tracker.header.maxnTrackers ) {
        state.trackerIdx %= data.tracker.header.maxnTrackers;
    }
}

void decrementTrackerIdx50 ( void ) {
    state.trackerIdx -= 50;
    while ( state.trackerIdx < 0 ) {
        state.trackerIdx += data.tracker.header.maxnTrackers;
    }
}

void resetRotation ( void ) {
    state.flag.rotation.x = false; state.rotation.x = 0;
    state.flag.rotation.y = false; state.rotation.y = 0;
    state.flag.rotation.z = false; state.rotation.z = 0;
}

void finish ( void ) {
    glutLeaveMainLoop();
    fclose ( data.intensity.fp );
    if ( data.segmentation.fp != NULL ) {
        fclose(data.intensity.fp);
    }
    exit ( 0 );
}

/* Viewr operation
    a : null
    b : Toggle display mode regarding trackers: all/one-by-one
    c : Increase the cutoff of voxel intensities
    d : Decrease the cutoff of voxel intensities
    e : Enable/Disable the voxel aspect ratio along z-axis
    f : Enable/Disable focus mode
    h : null
    i : null
    j : null
    k : null
    l : null
    m : Next tracker in forcus mode
    n : Previous tracker in forcus mode
    o : Draw/Delete raw intensity data
    q : Exit viewer
    r : Reset zoom
    s : Start/Stop
    t : Draw/Delete tracker positions
    u : null
    v : null
    w : Reset rotation originated from keyboard operation
    x : rotate around x-axis
    y : rotate around y-axis
    z : rotate around z-axis
  ESC : Exit viewer
    . : Reset rotation originated from mouse operation
    , : Reset translation originated from mouse operation
*/
void KEYBOARD ( unsigned char key, int x, int y ) {
    switch ( key ) {
        case 's' : toggleFlag(&state.flag.move); break;
        case 'o' : toggleFlag(&state.flag.plot.intensity); break;
        case 't' : trackerDisplayMode(); break;
        case 'a' : toggleFlag(&state.flag.plot.names); break;
        case 'b' : toggleFlag(&state.flag.parts); break;
        case 'm' : incrementTrackerIdx(); break;
        case 'n' : decrementTrackerIdx(); break;
        case '+' : incrementTrackerIdx50(); break;
        case '-' : decrementTrackerIdx50(); break;
        case 'f' : toggleFlag(&state.flag.focus); break;
        case 'c' : if (++state.cutoff>255){state.cutoff=255;}; break;
        case 'd' : if (--state.cutoff<  0){state.cutoff=  0;}; break;
        case 'x' : toggleFlag(&state.flag.rotation.x); break;
        case 'y' : toggleFlag(&state.flag.rotation.y); break;
        case 'z' : toggleFlag(&state.flag.rotation.z); break;
        case '.' : resetRotate(); break;
        case ',' : resetTranslation(); break;
        case 'w' : resetRotate(); resetRotation();break;
        case 'e' : toggleFlag(&state.flag.resize);break;
        case 'r' : state.zoomPercentage=100; break;
        case 'q' : finish(); break;
        case 'p' : toggleFlag(&state.cursor.enable);break;
        case '[' : state.cursor.position.x--; break;
        case ']' : state.cursor.position.x++; break;
        case '<' : state.cursor.position.y--; break;
        case '>' : state.cursor.position.y++; break;
        case '{' : state.cursor.position.z--; break;
        case '}' : state.cursor.position.z++; break;
        case '\033' : finish(); break; // ESC
        default : break;
    }
}

void scrShot ( void ) {
    int t,T;
    int tmp_frameIdx;
    FILE *fp;
    int w,h;
    struct { unsigned char r,g,b; } *buffer,*pixels;
    if ( (fp=fopen(parameters.scrShotFileName,"wb")) == NULL ) {
        fprintf(stderr,"%s is unwritable.(%d,%s)\n"
                ,parameters.scrShotFileName,__LINE__,__FILE__);
        return;
    }
    w = parameters.window.size.w;
    h = parameters.window.size.h;
    buffer = malloc ( 3*h*w );
    pixels = malloc ( 3*h*w );
    T = data.intensity.header.nFrames;
    tmp_frameIdx = state.frame;
    for ( t=0; t<T; t++ ) {
        int i,j;
        state.frame = t;
        DISPLAY();
        glReadPixels(0,0,w,h,GL_BGR,GL_UNSIGNED_BYTE,buffer);
        for ( i=0; i<h; i++ ) {
            for ( j=0; j<w; j++ ) {
                pixels[j+w*i] = buffer[j+(h-i-1)*w];
            }
        }
        fwrite(pixels,3,h*w,fp);
    }
    state.frame = tmp_frameIdx;
    free(buffer);
    free(pixels);
    fclose(fp);
}

void SPECIAL ( int key, int x, int y ) {
    switch ( key ) {
        case GLUT_KEY_UP :    zoomIn();            break;
        case GLUT_KEY_DOWN :  zoomOut();           break;
        case GLUT_KEY_RIGHT : gotoNextFrame();     break;
        case GLUT_KEY_LEFT :  gotoPreviousFrame(); break;
        case GLUT_KEY_F2 :    scrShot();           break;
    }
}

void IDLE ( void ) {
    static clock_t current,previous=(clock_t)-1;
    int t;
    do {
        if ( (current = clock())==(clock_t)-1 ) {
            fprintf(stderr,"Clock error.(%d,%s)\n",__LINE__,__FILE__);
            finish();
        }
        if ( previous == (clock_t)-1 ) {
            break;
        }
        t = (current-previous)*1000/CLOCKS_PER_SEC;
    } while ( t < parameters.intervalTime_ms );
    previous = current;
    glutPostRedisplay();
}

void Init ( void ) {
    glClearColor( 0, 0, 0, 0 );
    glEnable( GL_DEPTH_TEST );
    glEnable( GL_POINT_SMOOTH );
    glEnable( GL_CULL_FACE );
    glEnable( GL_BACK );
    Q2R ( state.rotate, state.current );
    parameters.viewpoint.x = data.intensity.header.resolution.X * 1.5;
    parameters.viewpoint.y = data.intensity.header.resolution.Y * 1.5;
    parameters.viewpoint.z = data.intensity.header.resolution.Z * 1.5;
}

int main ( int argc, char *argv[] ) {
    char title[2048];

    initializeState();
    initializeParameters();
    if ( checkArguments( argc, argv ) == failure ) {
        return failure;
    }

    sprintf(title,"[viewer] %s %s %s %s\n",
        argv[1],argv[2],argc>3?argv[3]:"",argc>4?argv[4]:"");
    glutInitWindowPosition ( parameters.window.position.x ,
            parameters.window.position.y);
    glutInitWindowSize ( parameters.window.size.w ,
            parameters.window.size.h);
    glutInit ( &argc, argv );
    glutCreateWindow(title);
    glutInitDisplayMode ( GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH );
    glutDisplayFunc ( DISPLAY );
    glutReshapeFunc ( RESHAPE );
    glutMouseFunc ( MOUSE );
    glutMotionFunc ( MOTION );
    glutPassiveMotionFunc ( PASSIVE_MOTION );
    glutKeyboardFunc ( KEYBOARD );
    glutSpecialFunc ( SPECIAL );
    glutIdleFunc ( IDLE );
    Init();
    glutMainLoop();
    return 0;
}

