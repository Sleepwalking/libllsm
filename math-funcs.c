/*
libllsm - Low Level Speech Model
===
Copyright (c) 2015 Kanru Hua. All rights reserved.

Developed by: Kanru Hua

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal with
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

    Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimers.
    Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimers in the documentation
and/or other materials provided with the distribution.
    Neither the names of development group, institution, nor the names of its
contributors may be used to endorse or promote products derived from this
Software without specific prior written permission.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE CONTRIBUTORS
OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE SOFTWARE.
*/

#include "math-funcs.h"
#include <math.h>
#include <string.h>

void llsm_idft(FP_TYPE* xr, FP_TYPE* xi, FP_TYPE* yr, FP_TYPE* yi, int n) {
  FP_TYPE tpon = 2.0 * M_PI / n;
  for(int k = 0; k < n; k ++) {
    FP_TYPE sumr = 0;
    FP_TYPE sumi = 0;
    FP_TYPE tponk = tpon * k;
    if(xr != NULL && xi != NULL)
      for(int i = 0; i < n; i ++) {
        FP_TYPE re = cos(tponk * i);
        FP_TYPE im = sin(tponk * i);
        sumr += xr[i] * re - xi[i] * im;
        sumi += xr[i] * im + xi[i] * re;
      }
    else if(xr != NULL && xi == NULL)
      for(int i = 0; i < n; i ++) {
        FP_TYPE re = cos(tponk * i);
        FP_TYPE im = sin(tponk * i);
        sumr += xr[i] * re;
        sumi += xr[i] * im;
      }
    else if(xr == NULL && xi != NULL)
      for(int i = 0; i < n; i ++) {
        FP_TYPE re = cos(tponk * i);
        FP_TYPE im = sin(tponk * i);
        sumr -= xi[i] * im;
        sumi += xi[i] * re;
      }
    if(yr != NULL) yr[k] = sumr / n;
    if(yi != NULL) yi[k] = sumi / n;
  }
}

FP_TYPE* llsm_winfir(int order, FP_TYPE cutoff, char* type, char* window) {
  FP_TYPE cutk = cutoff * order;
  FP_TYPE* freqrsp = calloc(order, sizeof(FP_TYPE));
  FP_TYPE* timersp = calloc(order, sizeof(FP_TYPE));

  FP_TYPE* w = NULL;
  if(! strcmp(window, "hanning"))
    w = hanning(order);
  else if(! strcmp(window, "hamming"))
    w = hamming(order);
  else if(! strcmp(window, "blackman_harris"))
    w = blackman_harris(order);
  
  if(! strcmp(type, "lowpass")) {
    for(int i = 0; i <= floor(cutk); i ++)
      freqrsp[i] = 1;
    freqrsp[(int)ceil(cutk) + 1] = fmod(cutk, 1.0);
  } else
  if(! strcmp(type, "highpass")) {
    for(int i = order / 2; i > floor(cutk); i --)
      freqrsp[i] = 1;
    freqrsp[(int)floor(cutk)] = 1.0 - fmod(cutk, 1.0);
  }
  
  complete_symm(freqrsp, order);
  idft(freqrsp, NULL, timersp, NULL, order);
  FP_TYPE* h = fftshift(timersp, order);
  for(int i = 0; i < order; i ++)
    h[i] *= w[i];
  
  free(w);
  free(freqrsp);
  free(timersp);
  return h;
}

FP_TYPE* llsm_convolution(FP_TYPE* x, FP_TYPE* h, int nx, int nh) {
  FP_TYPE* xpad = calloc(nx + nh * 2 - 1, sizeof(FP_TYPE));
  FP_TYPE* y = calloc(nx + nh - 1, sizeof(FP_TYPE));
  for(int i = 0; i < nx; i ++)
    xpad[i + nh - 1] = x[i];
  for(int i = 0; i < nx + nh - 1; i ++) {
    FP_TYPE sum = 0;
    for(int k = 0; k < nh; k ++)
      sum += h[k] * xpad[i + nh - 1 - k];
    y[i] = sum;
  }
  free(xpad);
  return y;
}

FP_TYPE* llsm_interp(FP_TYPE* xi, FP_TYPE* yi, int ni, FP_TYPE* x, int nx) {
  FP_TYPE* y = calloc(nx, sizeof(FP_TYPE));
  int srcidx = 0;
  for(int i = 0; i < nx; i ++) {
    FP_TYPE dstx = x[i];
    while(srcidx + 1 < ni && xi[srcidx + 1] < dstx) srcidx ++;
    int i0 = srcidx == 0 ? 0 : srcidx;
    int i1 = srcidx == ni - 1 ? srcidx : srcidx + 1;
    if(i0 != i1 && dstx > xi[0])
      y[i] = (yi[i1] - yi[i0]) * (dstx - xi[i0]) / (xi[i1] - xi[i0]) + yi[i0];
    else
      y[i] = yi[i0];
  }
  return y;
}

