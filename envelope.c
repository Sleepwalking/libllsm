/*
  libllsm - Low Level Speech Model
  ===
  Copyright (c) 2015-2017 Kanru Hua.

  libllsm is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  libllsm is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with libllsm. If not, see <http://www.gnu.org/licenses/>.
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

FP_TYPE* llsm_true_envelope(FP_TYPE* spectrum, int nfft, int order, int niter) {
  int ns = nfft / 2 + 1;
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

FP_TYPE* llsm_warp_freq(FP_TYPE fmin, FP_TYPE fmax, int n, FP_TYPE warp_const) {
  FP_TYPE* freq = calloc(n + 1, sizeof(FP_TYPE));
  FP_TYPE wmin = 5000.0 * log(1 + fmin / warp_const);
  FP_TYPE wmax = 5000.0 * log(1 + fmax / warp_const);
  for(int i = 0; i <= n; i ++)
    freq[i] = warp_const * (exp(((FP_TYPE)i / n * (wmax - wmin) + wmin) / 5000.0) - 1.0);
  return freq;
}

FP_TYPE* llsm_geometric_envelope(FP_TYPE* spectrum, int nfft, int fs, FP_TYPE* freq, int nf) {
  FP_TYPE* env = calloc(nf, sizeof(FP_TYPE));
  int ns = nfft / 2 + 1;
  for(int i = 0; i < nf; i ++) {
    FP_TYPE fprev = i == 0 ? freq[0] - 50.0 : freq[i - 1];
    FP_TYPE fnext = freq[i + 1];
    FP_TYPE fcenter = freq[i];
    int idxl = max(0, floor((fcenter + fprev) / 2.0 / fs * nfft));
    int idxh = ceil((fcenter + fnext) / 2.0 / fs * nfft);
    if(idxh >= ns)
      env[i] = env[i - 1];
    else
      env[i] = sumfp(spectrum + idxl, idxh - idxl + 1) / (idxh - idxl + 1);
  }
  return env;
}

FP_TYPE* llsm_spectrum_from_envelope(FP_TYPE* freq, FP_TYPE* magn, int nf, int nfft, int fs) {
  int ns = nfft / 2 + 1;
  FP_TYPE* x = calloc(ns, sizeof(FP_TYPE));
  for(int i = 0; i < ns; i ++)
    x[i] = (FP_TYPE)i * fs / nfft;
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

