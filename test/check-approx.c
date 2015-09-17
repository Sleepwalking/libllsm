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
#include <math.h>
#include "../external/matlabfunctions.h"
#include "../math-funcs.h"

#define check_single_var(func, stdfunc, lbound, ubound, wrap) \
  { \
    printf("Checking " #func " ...\n"); \
    FP_TYPE err_max = 0; \
    FP_TYPE err_sum = 0; \
    FP_TYPE err_sqr = 0; \
    FP_TYPE err_max_p = 0; \
    FP_TYPE err_sum_p = 0; \
    FP_TYPE step = (ubound - lbound) / 1000.0; \
    for(FP_TYPE x = lbound; x < ubound; x += step) { \
      FP_TYPE err = func(wrap(x)) - stdfunc(wrap(x)); \
      FP_TYPE err_p = fabs(err) / stdfunc(wrap(x)); \
      if(fabs(err) > err_max) err_max = fabs(err); \
      if(err_p > err_max_p) err_max_p = err_p; \
      err_sum += err; \
      err_sqr += err * err; \
      err_sum_p += err_p; \
    } \
    FP_TYPE mean = err_sum / 1000.0; \
    printf("Max error = %e, mean = %e, std = %e\n", err_max, mean, sqrt(err_sqr / 1000.0 - mean * mean)); \
    printf("Max relative error = %e, mean = %e\n", err_max_p, err_sum_p / 1000.0); \
  }

#define uniform(x) x

int main(int argc, char** argv) {
  check_single_var(fastcosfull, cos, -M_PI, M_PI, uniform);
  check_single_var(fastsinfull, sin, -M_PI, M_PI, uniform);
  check_single_var(fastlog, log, -20.0, 10.0, exp);
  check_single_var(fastexp, exp, -20.0, 10.0, uniform);

  { // we're not going to use fastatan2 due to its large approximation error as revealed below
    printf("Checking fastatan2 ...\n");
    const FP_TYPE lbound = -20.0;
    const FP_TYPE ubound = 20.0;
    FP_TYPE err_max = 0;
    FP_TYPE err_sum = 0;
    FP_TYPE err_sqr = 0;
    FP_TYPE step = (ubound - lbound) / 100.0;
    for(FP_TYPE x = lbound; x < ubound; x += step) {
      for(FP_TYPE y = lbound; y < ubound; y += step) {
        FP_TYPE err = fastatan2(y, x) - atan2(y, x);
        if(fabs(err) > err_max) err_max = fabs(err);
        err_sum += err;
        err_sqr += err * err;
      }
    }
    FP_TYPE mean = err_sum / 10000.0;
    printf("Max error = %e, mean = %e, std = %e\n", err_max, mean, sqrt(err_sqr / 10000.0 - mean * mean));
  }
}

