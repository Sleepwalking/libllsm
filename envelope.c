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

#include "envelope.h"
#include "math-funcs.h"

void llsm_reduce_spectrum_depth(FP_TYPE* spectrum, int ns, int nhop, FP_TYPE minimum, FP_TYPE depth) {
  int nm = ns / nhop;
  FP_TYPE* smax = calloc(nm, sizeof(FP_TYPE));
  for(int i = 0; i < nm; i ++)
    smax[i] = maxfp(& spectrum[i * nhop], nhop);
  for(int i = 0; i < nm; i ++) {
    if(smax[i] < minimum) {
      int j = 0;
      while(j < nm) {
        if(i - j >= 0 && smax[i - j] > minimum) {
          smax[i] = smax[i - j];
          break;
        } else
        if(i + j < nm && smax[i + j] > minimum) {
          smax[i] = smax[i + j];
          break;
        } else
        if(i + j < 0 || i + j >= nm) {
          smax[i] = minimum;
          break;
        }
        j ++;
      }
    }
  }
  for(int i = 0; i < ns; i ++) {
    FP_TYPE base = smax[i / nhop] - depth;
    spectrum[i] = spectrum[i] < minimum ? base : spectrum[i];
  }
  free(smax);
}

FP_TYPE* llsm_true_envelope(FP_TYPE* spectrum, int ns, int order, int niter) {
  int nfft = ns * 2;
  FP_TYPE* y  = calloc(nfft, sizeof(FP_TYPE));
  FP_TYPE* ya = calloc(nfft, sizeof(FP_TYPE));
  FP_TYPE* c  = calloc(nfft, sizeof(FP_TYPE));
  FP_TYPE* c2 = calloc(nfft, sizeof(FP_TYPE));
  FP_TYPE fftbuff[65536];
  
  for(int i = 0; i < ns; i ++)
    y[i] = spectrum[i];
  complete_symm(y, nfft);
  for(int i = 0; i < nfft; i ++)
    ya[i] = y[i];
  
  for(int n = 0; n < niter; n ++) {
    if(n > 0)
      fft(c, NULL, y, NULL, nfft, fftbuff);
    for(int i = 0; i < nfft; i ++)
      y[i] = max(y[i], ya[i]);
    ifft(y, NULL, c2, NULL, nfft, fftbuff);
    if(n > 1)
      for(int i = 0; i < nfft; i ++)
        c[i] = c2[i] + (c2[i] - c[i]) * 3.0;
    else
      for(int i = 0; i < nfft; i ++)
        c[i] = c2[i];
    for(int i = order; i <= nfft - order; i ++)
      c[i] = 0;
  }
  
  free(y);
  free(ya);
  free(c2);
  return c;
}

FP_TYPE* llsm_wrap_freq(FP_TYPE fmin, FP_TYPE fmax, int n, FP_TYPE wrap_const) {
  FP_TYPE* freq = calloc(n + 1, sizeof(FP_TYPE));
  FP_TYPE wmin = 5000.0 * log(1 + fmin / wrap_const);
  FP_TYPE wmax = 5000.0 * log(1 + fmax / wrap_const);
  for(int i = 0; i <= n; i ++)
    freq[i] = wrap_const * (exp(((FP_TYPE)i / n * (wmax - wmin) + wmin) / 5000.0) - 1.0);
  return freq;
}

FP_TYPE* llsm_geometric_envelope(FP_TYPE* spectrum, int ns, int fs, FP_TYPE* freq, int nf) {
  FP_TYPE* env = calloc(nf, sizeof(FP_TYPE));
  int nfft = ns * 2;
  for(int i = 0; i < nf; i ++) {
    FP_TYPE fprev = i == 0 ? freq[0] - 50.0 : freq[i - 1];
    FP_TYPE fnext = freq[i + 1];
    FP_TYPE fcenter = freq[i];
    int idxl = max(0, floor((fcenter + fprev) / 2.0 / fs * nfft));
    int idxh = min(ns - 1, ceil((fcenter + fnext) / 2.0 / fs * nfft));
    // mean energy representation seems to introduce magnitude distortion
    /*
    env[i] = 0;
    for(int j = idxl; j <= idxh; j ++)
      env[i] += exp(spectrum[j] * 2.0);
    env[i] /= idxh - idxl + 1;
    env[i] = log(env[i]) / 2.0;
    */
    env[i] = maxfp(spectrum + idxl, idxh - idxl + 1);
  }
  return env;
}

FP_TYPE* llsm_spectrum_from_envelope(FP_TYPE* freq, FP_TYPE* magn, int nf, int ns, int fs) {
  FP_TYPE* x = calloc(ns, sizeof(FP_TYPE));
  int nfft = ns * 2;
  for(int i = 0; i < ns; i ++)
    x[i] = (double)i * fs / nfft;
  FP_TYPE* y = interp1(freq, magn, nf, x, ns);
  free(x);
  return y;
}

FP_TYPE* llsm_nonuniform_envelope(FP_TYPE* x, int nx, int* instant, int* winlen, int ni, int mode) {
  FP_TYPE* env = calloc(ni, sizeof(FP_TYPE));
  for(int i = 0; i < ni; i ++) {
    int t = instant[i];
    int l = winlen[i];
    env[i] = x[t];
    if(mode == 0)
      for(int j = 0; j < l; j ++) {
        int idx = t - l / 2 + j;
        if(idx >= 0 && idx < nx)
          env[i] = min(env[i], x[idx]);
      }
    else
      for(int j = 0; j < l; j ++) {
        int idx = t - l / 2 + j;
        if(idx >= 0 && idx < nx)
          env[i] = max(env[i], x[idx]);
      }
  }
  return env;
}

