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

#ifndef LLSM_ENV
#define LLSM_ENV

#include "common.h"

/*
  These two functions are not used in the current version but we might find them helpful some day.
  llsm_reduce_spectrum_depth: fill up low energy parts in a magnitude spectrum
  llsm_true_envelope: find the true envelope (Robel and Rodet, 2005) given a magnitude spectrum; returns cepstral coefficients
*/
void llsm_reduce_spectrum_depth(FP_TYPE* spectrum, int ns, int nhop, FP_TYPE minimum, FP_TYPE depth);
FP_TYPE* llsm_true_envelope(FP_TYPE* spectrum, int ns, int order, int niter);

/*
  llsm_wrap_freq: generate an array of wrapped frequencies from fmin to fmax; equivalent to mel-frequency
    spacing when wrap_const = 700.0
  llsm_geometric_envelope: obtain a geometric (piecewise linear) representation of spectrum by taking maximum
    values around freq
  llsm_spectrum_from_envelope: interpolate a piecewise linear function to generate a spectral envelope
  llsm_nonuniform_envelope: obtain a geometric representation of any signal by taking local maximum/minimum
    around nonuniformly sampled time instants; mode: 0 - minimum, 1 - maximum
*/
FP_TYPE* llsm_wrap_freq(FP_TYPE fmin, FP_TYPE fmax, int n, FP_TYPE wrap_const);
FP_TYPE* llsm_geometric_envelope(FP_TYPE* spectrum, int ns, int fs, FP_TYPE* freq, int nf);
FP_TYPE* llsm_spectrum_from_envelope(FP_TYPE* freq, FP_TYPE* magn, int nf, int ns, int fs);
FP_TYPE* llsm_nonuniform_envelope(FP_TYPE* x, int nx, int* instant, int* winlen, int ni, int mode);

#endif

