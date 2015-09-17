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
#include "filter-coef.h"

void llsm_idft(FP_TYPE* xr, FP_TYPE* xi, FP_TYPE* yr, FP_TYPE* yi, int n) {
  FP_TYPE tpon = 2.0 * M_PI / n;
  for(int k = 0; k < n; k ++) {
    FP_TYPE sumr = 0;
    FP_TYPE sumi = 0;
    FP_TYPE tponk = tpon * k;
    if(xr != NULL && xi != NULL)
      for(int i = 0; i < n; i ++) {
        FP_TYPE re = cos_3(tponk * i);
        FP_TYPE im = sin_3(tponk * i);
        sumr += xr[i] * re - xi[i] * im;
        sumi += xr[i] * im + xi[i] * re;
      }
    else if(xr != NULL && xi == NULL)
      for(int i = 0; i < n; i ++) {
        FP_TYPE re = cos_3(tponk * i);
        FP_TYPE im = sin_3(tponk * i);
        sumr += xr[i] * re;
        sumi += xr[i] * im;
      }
    else if(xr == NULL && xi != NULL)
      for(int i = 0; i < n; i ++) {
        FP_TYPE re = cos_3(tponk * i);
        FP_TYPE im = sin_3(tponk * i);
        sumr -= xi[i] * im;
        sumi += xi[i] * re;
      }
    if(yr != NULL) yr[k] = sumr / n;
    if(yi != NULL) yi[k] = sumi / n;
  }
}

FP_TYPE* llsm_winfir(int order, FP_TYPE cutoff, FP_TYPE cutoff2, char* type, char* window) {
  FP_TYPE cutk  = cutoff  * order;
  FP_TYPE cutk2 = cutoff2 * order;
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
    freqrsp[(int)ceil(cutk)] = fmod(cutk, 1.0);
  } else
  if(! strcmp(type, "highpass")) {
    for(int i = order / 2; i > floor(cutk); i --)
      freqrsp[i] = 1;
    freqrsp[(int)floor(cutk)] = 1.0 - fmod(cutk, 1.0);
  }
  if(! strcmp(type, "bandpass")) {
    for(int i = floor(cutk2); i > floor(cutk); i --)
      freqrsp[i] = 1;
    freqrsp[(int)floor(cutk)] = 1.0 - fmod(cutk, 1.0);    
    freqrsp[(int)ceil(cutk2)] = fmod(cutk2, 1.0);
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

int llsm_get_iir_filter(FP_TYPE cutoff, char* type, FP_TYPE** a, FP_TYPE** b) {
  int n = max(0, round(cutoff * 2.0 / step_freq - 1));
  if(n >= filter_number) n = filter_number - 1;
  int order = coef_size;
  *a = calloc(order, sizeof(FP_TYPE));
  *b = calloc(order, sizeof(FP_TYPE));
  const FP_TYPE* a_line, *b_line;
  if(! strcmp(type, "lowpass")) {
    a_line = cheby_l_a + n * coef_size;
    b_line = cheby_l_b + n * coef_size;
  } else {
    a_line = cheby_h_a + n * coef_size;
    b_line = cheby_h_b + n * coef_size;
  }
  for(int i = 0; i < order; i ++) {
    (*a)[i] = a_line[i];
    (*b)[i] = b_line[i];
  }
  return order;
}

FP_TYPE* llsm_convolution(FP_TYPE* x, FP_TYPE* h, int nx, int nh) {
  FP_TYPE* xpad = calloc(nx + nh * 2 - 1, sizeof(FP_TYPE));
  FP_TYPE* y = calloc(nx + nh - 1, sizeof(FP_TYPE));
  for(int i = 0; i < nx; i ++)
    xpad[i + nh - 1] = x[i];
  for(int i = 0; i < nx + nh - 1; i ++)
    for(int k = 0; k < nh; k ++)
      y[i] += h[k] * xpad[i + nh - 1 - k];
  free(xpad);
  return y;
}

static FP_TYPE* llsm_filter_order6(FP_TYPE* b, FP_TYPE* a, FP_TYPE* x, int nx) {
  FP_TYPE* y = calloc(nx + 5, sizeof(FP_TYPE));
  for(int i = 5; i < nx; i ++) {
      y[i] -= a[1] * y[i - 1] + a[2] * y[i - 2] + a[3] * y[i - 3] + a[4] * y[i - 4] + a[5] * y[i - 5];
      y[i] += b[0] * x[i - 0] + b[1] * x[i - 1] + b[2] * x[i - 2] + b[3] * x[i - 3] + b[4] * x[i - 4] + b[5] * x[i - 5];
  }
  return y;
}

FP_TYPE* llsm_filter(FP_TYPE* b, int nb, FP_TYPE* a, int na, FP_TYPE* x, int nx) {
  if(na == 6 && nb == 6)
    return llsm_filter_order6(b, a, x, nx);

  int nh = max(na, nb);
  FP_TYPE* y = calloc(nx + nh - 1, sizeof(FP_TYPE));
  for(int i = na - 1; i < nx; i ++) {
    for(int k = 1; k < na; k ++)
      y[i] -= a[k] * y[i - k];
    for(int k = 0; k < nb; k ++)
      y[i] += b[k] * x[i - k];
  }
  return y;
}

FP_TYPE* llsm_chebyfilt(FP_TYPE* x, int nx, FP_TYPE cutoff1, FP_TYPE cutoff2, char* type) {
  if(! strcmp(type, "bandpass")) {
    FP_TYPE* x1 = llsm_chebyfilt(x , nx, cutoff1, 0, "highpass");
    FP_TYPE* y  = llsm_chebyfilt(x1, nx, cutoff2, 0, "lowpass");
    free(x1);
    return y;
  }
  FP_TYPE* a, *b;
  int order = llsm_get_iir_filter(cutoff1, type, &a, &b);
  FP_TYPE* y = llsm_filter(b, order, a, order, x, nx);
  free(a);
  free(b);
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

