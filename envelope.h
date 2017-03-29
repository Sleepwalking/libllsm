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

#ifndef LLSM_ENV
#define LLSM_ENV

/*
  These two functions are not used in the current version but we might find them helpful some day.
  llsm_reduce_spectrum_depth: fill up low energy parts in a magnitude spectrum
  llsm_true_envelope: find the true envelope (Robel and Rodet, 2005) given a magnitude spectrum; returns cepstral coefficients
*/
void llsm_reduce_spectrum_depth(FP_TYPE* spectrum, int ns, int nhop, FP_TYPE minimum, FP_TYPE depth);
FP_TYPE* llsm_true_envelope(FP_TYPE* spectrum, int nfft, int order, int niter);

/*
  llsm_warp_freq: generate an array of warpped frequencies from fmin to fmax; equivalent to mel-frequency
    spacing when warp_const = 700.0
  llsm_geometric_envelope: obtain a geometric (piecewise linear) representation of spectrum by taking mean
    values around freq
  llsm_spectrum_from_envelope: interpolate a piecewise linear function to generate a spectral envelope
  llsm_nonuniform_envelope: obtain a geometric representation of any signal by taking local maximum/minimum
    around nonuniformly sampled time instants; mode: 0 - minimum, 1 - maximum
*/
FP_TYPE* llsm_warp_freq(FP_TYPE fmin, FP_TYPE fmax, int n, FP_TYPE warp_const);
FP_TYPE* llsm_geometric_envelope(FP_TYPE* spectrum, int nfft, int fs, FP_TYPE* freq, int nf);
FP_TYPE* llsm_spectrum_from_envelope(FP_TYPE* freq, FP_TYPE* magn, int nf, int nfft, int fs);
FP_TYPE* llsm_nonuniform_envelope(FP_TYPE* x, int nx, int* instant, int* winlen, int ni, int mode);

#endif

