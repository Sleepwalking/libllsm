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

#include "math-funcs.h"
#include "filter-coef.h"

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

FP_TYPE* llsm_chebyfilt(FP_TYPE* x, int nx, FP_TYPE cutoff1, FP_TYPE cutoff2, char* type) {
  if(! strcmp(type, "bandpass")) {
    FP_TYPE* x1 = llsm_chebyfilt(x , nx, cutoff1, 0, "highpass");
    FP_TYPE* y  = llsm_chebyfilt(x1, nx, cutoff2, 0, "lowpass");
    free(x1);
    return y;
  }
  FP_TYPE* a, *b;
  int order = llsm_get_iir_filter(cutoff1, type, &a, &b);
  FP_TYPE* y = filtfilt(b, order, a, order, x, nx);

  free(a);
  free(b);
  return y;
}

