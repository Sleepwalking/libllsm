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

#ifndef LLSM_MFUNCS
#define LLSM_MFUNCS

#include "common.h"
#include <stdlib.h>

#define def_singlepass(name, op, init) \
inline FP_TYPE name(FP_TYPE* src, int n) { \
  FP_TYPE ret = init; \
  for(int i = 0; i < n; i ++) \
    ret = op(ret, src[i]); \
  return ret; \
}

#define def_add(a, b) ((a) + (b))
#define def_max(a, b) ((a) > (b) ? (a) : (b))
#define def_min(a, b) ((a) < (b) ? (a) : (b))

def_singlepass(sumfp, def_add, 0)
def_singlepass(maxfp, def_max, src[0])
def_singlepass(minfp, def_min, src[0])

inline FP_TYPE* boxcar(int n) {
  FP_TYPE* ret = calloc(n, sizeof(FP_TYPE));
  for(int i = 0; i < n; i ++)
    ret[i] = 1.0;
  return ret;
}

inline FP_TYPE* hanning(int n) {
  FP_TYPE* ret = calloc(n, sizeof(FP_TYPE));
  for(int i = 0; i < n; i ++)
    ret[i] = 0.5 * (1 - cos(2 * M_PI * i / (n - 1)));
  return ret;
}

inline FP_TYPE* hamming(int n) {
  FP_TYPE* ret = calloc(n, sizeof(FP_TYPE));
  for(int i = 0; i < n; i ++)
    ret[i] = 0.54 - 0.46 * cos(2 * M_PI * i / (n - 1));
  return ret;
}

inline FP_TYPE* mltsine(int n) {
  FP_TYPE* ret = calloc(n, sizeof(FP_TYPE));
  for(int i = 0; i < n; i ++)
    ret[i] = sin(M_PI / n * (i + 0.5));
  return ret;
}

inline FP_TYPE* blackman_harris(int n) {
  FP_TYPE* ret = calloc(n, sizeof(FP_TYPE));
  const FP_TYPE a0 = 0.35875;
  const FP_TYPE a1 = 0.48829;
  const FP_TYPE a2 = 0.14128;
  const FP_TYPE a3 = 0.01168;
  for(int i = 0; i < n; i ++)
    ret[i] = a0 - a1 * cos(2.0 * M_PI * i / n) +
                  a2 * cos(4.0 * M_PI * i / n) -
                  a3 * cos(6.0 * M_PI * i / n);
  return ret;
}

void cdft(int n, int isgn, FP_TYPE* a);
void rdft(int n, int isgn, FP_TYPE* a);
void llsm_idft(FP_TYPE* xr, FP_TYPE* xi, FP_TYPE* yr, FP_TYPE* yi, int n);
FP_TYPE* llsm_winfir(int order, FP_TYPE cutoff, char* type, char* window);
FP_TYPE* llsm_convolution(FP_TYPE* x, FP_TYPE* h, int nx, int nh);
FP_TYPE* llsm_interp(FP_TYPE* xi, FP_TYPE* yi, int ni, FP_TYPE* x, int nx);

inline double fastatan2(double y, double x) {
  double coeff_1 = M_PI / 4.0;
  double coeff_2 = 3.0 * coeff_1;
  double abs_y = fabs(y) + 1e-10; // kludge to prevent 0/0 condition
  double angle = 0;
  if(x >= 0) {
    double r = (x - abs_y) / (x + abs_y);
    angle = coeff_1 - coeff_1 * r;
  } else {
    double r = (x + abs_y) / (abs_y - x);
    angle = coeff_2 - coeff_1 * r;
  }
  if(y < 0)
    return -angle; // negate if in quad III or IV
  else
   return angle;
}

inline void fft_core(FP_TYPE* xr, FP_TYPE* xi, FP_TYPE* yr, FP_TYPE* yi, int n, FP_TYPE* buffer, FP_TYPE mode) {
  for(int i = 0; i < n; i ++) {
    buffer[i * 2] = xr == NULL ? 0 : xr[i];
    buffer[i * 2 + 1] = xi == NULL ? 0 : xi[i];
  }
  cdft(2 * n, mode, buffer);
  if(mode < 0)
    for(int i = 0; i < n; i ++) {
      if(yr != NULL) yr[i] = buffer[i * 2];
      if(yi != NULL) yi[i] = buffer[i * 2 + 1];
    }
  else
    for(int i = 0; i < n; i ++) {
      if(yr != NULL) yr[i] = buffer[i * 2] / n;
      if(yi != NULL) yi[i] = buffer[i * 2 + 1] / n;
    }
}

inline void fft(FP_TYPE* xr, FP_TYPE* xi, FP_TYPE* yr, FP_TYPE* yi, int n, FP_TYPE* buffer) {
  fft_core(xr, xi, yr, yi, n, buffer, -1.0);
}

inline void ifft(FP_TYPE* xr, FP_TYPE* xi, FP_TYPE* yr, FP_TYPE* yi, int n, FP_TYPE* buffer) {
  fft_core(xr, xi, yr, yi, n, buffer, 1.0);
}

inline void idft(FP_TYPE* xr, FP_TYPE* xi, FP_TYPE* yr, FP_TYPE* yi, int n) {
  llsm_idft(xr, xi, yr, yi, n);
}

inline FP_TYPE* fftshift(FP_TYPE* x, int n) {
  FP_TYPE* y = calloc(n, sizeof(FP_TYPE));
  int halfs = n / 2;
  int halfl = (n + 1) / 2;
  for(int i = 0; i < halfs; i ++)
    y[i] = x[i + halfl];
  for(int i = 0; i < halfl; i ++)
    y[i + halfs] = x[i];
  return y;
}

inline FP_TYPE* unwrap(FP_TYPE* x, int n) {
  FP_TYPE* y = calloc(n, sizeof(FP_TYPE));
  y[0] = x[0];
  for(int i = 1; i < n; i ++) {
    if(fabs(x[i] - x[i - 1]) > M_PI)
      y[i] = y[i - 1] + x[i] - (x[i - 1] + 2.0 * M_PI * (x[i] > x[i - 1] ? 1.0 : -1.0));
    else
      y[i] = y[i - 1] + x[i] - x[i - 1];
  }
  return y;
}

inline FP_TYPE* abscplx(FP_TYPE* xr, FP_TYPE* xi, int n) {
  FP_TYPE* y = calloc(n, sizeof(FP_TYPE));
  for(int i = 0; i < n; i ++)
    y[i] = sqrt(xr[i] * xr[i] + xi[i] * xi[i]);
  return y;
}

inline FP_TYPE* argcplx(FP_TYPE* xr, FP_TYPE* xi, int n) {
  FP_TYPE* y = calloc(n, sizeof(FP_TYPE));
  for(int i = 0; i < n; i ++)
    y[i] = atan2(xi[i], xr[i]);
  return y;
}

inline FP_TYPE linterp(FP_TYPE v1, FP_TYPE v2, FP_TYPE ratio) {
    return v1 + (v2 - v1) * ratio;
}

inline void complete_symm(FP_TYPE* x, int n) {
  if(n / 2 == (n + 1) / 2) // even
    x[n / 2] = x[n / 2 - 1];
  for(int i = n / 2 + 1; i < n; i ++)
    x[i] = x[n - i];
}

inline void complete_asymm(FP_TYPE* x, int n) {
  if(n / 2 == (n + 1) / 2) // even
    x[n / 2] = x[n / 2 - 1];
  for(int i = n / 2 + 1; i < n; i ++)
    x[i] = -x[n - i];
}

inline FP_TYPE* fir1(int order, FP_TYPE cutoff, char* type, char* window) {
  return llsm_winfir(order, cutoff / 2.0, type, window);
}

inline FP_TYPE* conv(FP_TYPE* x, FP_TYPE* h, int nx, int nh) {
  return llsm_convolution(x, h, nx, nh);
}

inline FP_TYPE* interp1(FP_TYPE* xi, FP_TYPE* yi, int ni, FP_TYPE* x, int nx) {
  return llsm_interp(xi, yi, ni, x, nx);
}

inline FP_TYPE* white_noise(FP_TYPE amplitude, int n) {
  FP_TYPE* y = calloc(n, sizeof(FP_TYPE));
  for(int i = 0; i < n; i ++)
    y[i] = ((FP_TYPE)rand() / RAND_MAX - 0.5) * amplitude * 2.0;
  return y;
}

inline FP_TYPE* moving_avg(FP_TYPE* x, int nx, int order) {
  FP_TYPE* h = boxcar(order);
  for(int i = 0; i < order; i ++) h[i] /= order;
  FP_TYPE* y = conv(x, h, nx, order);
  free(h);
  return y;
}

#endif

