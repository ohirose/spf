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

double wdist    (const double * x1, const double * x2, const int P, const double dz);
int    localmax (double **vx, int *nvx, const byte *y, const int *imsize, const int *width);
int    cutoff   (byte *y, const int *imsize, const byte cut); 

int dpmeans(
      int           *  z,    /* OUTPUT: cluster id | size N           */ 
      double        ** m,    /* OUTPUT: means      | size K(N) x P    */
      double        ** d,    /* OUTPUT: distance^2 | size N    x K(N) */
      int           *  l,    /* OUTPUT: #members   | size K(N)        */
      int           *  K,    /* OUTPUT: #clusters  | size 1           */ 
      const double  ** X,    /* INPUT : data       | size N   x P     */ 
      const int        N,    /* INPUT : #samples   | size 1           */ 
      const int        P,    /* INPUT : #variables | size 1           */ 
      const double     lmd,  /* INPUT : parameter  | size 1           */ 
      const double     dz,   /* INPUT : zscale     | size 1           */
      const int        nlp   /* INPUT : #loops     | size 1           */
  );

