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

#ifndef LLSM_MFUNCS
#define LLSM_MFUNCS

#include "external/ciglet/ciglet.h"

#define LLSM_CHEBY_ORDER 5

FP_TYPE* llsm_chebyfilt(FP_TYPE* x, int nx, FP_TYPE cutoff1, FP_TYPE cutoff2, char* type);

static inline FP_TYPE* chebyfilt(FP_TYPE* x, int nx, FP_TYPE cutoff1, FP_TYPE cutoff2, char* type) {
  return llsm_chebyfilt(x, nx, cutoff1 / 2.0, cutoff2 / 2.0, type);
}

#endif

