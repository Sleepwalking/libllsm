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

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "../external/matlabfunctions.h"
#include "../external/libpyin/pyin.h"
#include "../llsm.h"

double get_time() {
  struct timeval t;
  gettimeofday(&t, NULL);
  return (t.tv_sec + (t.tv_usec / 1000000.0)) * 1000.0;
}

int main(int argc, char** argv) {
  if(argc < 2) {
    fprintf(stderr, "Missing argument.\n");
    return 1;
  }
  
  int fs = 0;
  int nbit = 0;
  int nx = 0;
  double* x = wavread(argv[1], & fs, & nbit, & nx);

  pyin_paramters param = pyin_init(pow(2, ceil(log2(fs * 0.005))));
  param.fmin = 50.0;
  param.fmax = 800.0;
  param.trange = 24;
  param.nf = ceil(fs * 0.025);
  param.w = param.nf / 4;
  
  int nfrm = 0;
  double* f0 = pyin_analyze(param, x, nx, fs, & nfrm);

  llsm_parameters lparam = llsm_init();
  lparam.a_mvf = 11000;
  lparam.a_nhop = pow(2, ceil(log2(fs * 0.005)));
  llsm* model = llsm_analyze(lparam, x, nx, fs, f0, nfrm);

  lparam.s_fs = fs;
  
  double t0 = get_time();
  int ny;
  double* y = llsm_synthesize(lparam, model, & ny);
  double tspent = get_time() - t0;
  printf("%f ms, %fs/s\n", tspent, (double)nx / fs / tspent * 1000);
  
  wavwrite(y, ny, lparam.s_fs, nbit, "resynth.wav");

  llsm_delete(model);
  free(f0);
  free(x);
  free(y);

  return 0;
}

