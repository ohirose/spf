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

#define MAX_FILENAME_LENGTH 256

typedef struct ImageHeader {
    char   fnRawImage[MAX_FILENAME_LENGTH];
    int    rawHeaderSize;
    int    nFrames;
    struct { int X; int Y; int Z;} resolution;
    struct { double y; double z; } voxelAspectRatio;
} ImageHeader;

typedef struct TrackerHeader {
    char fnTrackerRawData [MAX_FILENAME_LENGTH];
    char fnTrackerNames   [MAX_FILENAME_LENGTH];
    int  rawHeaderSize;          
    int  nFrames;            
    int  maxnTrackers;          
} TrackerHeader;

struct _tracker{
    TrackerHeader header;
    int    *nums;
    double *data;
    double *means;
    char   **names;
};

void loadImageInfo         (ImageHeader   *h, const char *file);
void loadTrackerInfo       (TrackerHeader *h, const char *file);
void loadTracker           (struct _tracker *tracker, const char *file);

