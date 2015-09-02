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

#ifndef LLSM_COMMON
#define LLSM_COMMON

#include <malloc.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#define FP_TYPE double

#ifndef M_PI
  #define M_PI 3.1415926535897932385
#endif
#define EPS 0.0000000001
#define INF 9999999999

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

inline void** malloc2d(size_t m, size_t n, size_t size) {
  void** ret = calloc(m, sizeof(void*));
  for(size_t i = 0; i < m; i++)
    ret[i] = calloc(n, size);
  return ret;
}

#define copy2d(src, m, n, size) copy2d_((void**)src, m, n, size)
inline void** copy2d_(void** src, size_t m, size_t n, size_t size) {
  void** ret = malloc2d(m, n, sizeof(void*));
  for(size_t i = 0; i < m; i ++)
    memcpy(ret[i], src[i], n * size);
  return ret;
}

#define free2d(ptr, m) free2d_((void**)(ptr), m)
inline void free2d_(void** ptr, size_t m) {
  for(size_t i = 0; i < m; i ++)
    free(ptr[i]);
  free(ptr);
}

inline FP_TYPE* fetch_frame(FP_TYPE* x, int nx, int center, int nf) {
  FP_TYPE* y = calloc(nf, sizeof(FP_TYPE));
  for(int i = 0; i < nf; i ++) {
    int isrc = center + i - nf / 2;
    y[i] = (isrc >= 0 && isrc < nx) ? x[isrc] : 0;
  }
  return y;
}

#endif

