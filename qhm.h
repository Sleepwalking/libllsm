/*
qhm.h - Quasi-Harmonic Model with Adaptive Iterative Refinement
===
Copyright (c) 2016 tuxzz. All rights reserved.
Developed by: tuxzz

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

#ifndef QHM
#define QHM

#include "llsm.h"
#include "complex.h"

#define QHM_LSMETHOD_QR 0
#define QHM_LSMETHOD_SVD 1

typedef struct {
  FP_TYPE* x; // original signal. shape = (nX, 1)
  FP_TYPE* window; // shape = (nX, 1)

  FP_TYPE f0;
  FP_TYPE maxcorr; // max correction per step (in Hz)
  int nx, fs, nhar;
  int ls_method;

  // run status
  int K, N; // (nHar * 2 + 1), (nX - 1) / 2
  FP_TYPE* fk_hat; // freq. of harmonics. (K, 1)
  FP_TYPE* n; // time vec. (1, nX)
  FP_TYPE* t; // arg of cplx exp, (K, 2N + 1)
  double complex* Ea; // matrix with cplx exp, (2K, 2N + 1)
  double complex* Ew; // E multiply by the window, (2N + 1, 2K)
  double complex* R; // E multiply by the window, (2K, 2K)
  double complex* windowed_x; // (nX, 1)
  double complex* lstsqb; // (2 * K, 1), contains ak and bk

} qhm_solve_status;

/*
* Switch of progress output, 1 by default.
*/
void qhm_progress(int v);

/*
* Least-Square F0 Adaptive Iteration Refinement(AIR)
*/
void qhm_air(llsm_parameters param, FP_TYPE* x, int nx, int fs,
  FP_TYPE* f0, int nhop, const char* wtype);

/*
* Least-Square Sinusoid parameter iteration anaysis.
*/
void qhm_analyze(llsm_parameters param, FP_TYPE* x, int nx, int fs, FP_TYPE* f0, int nhop,
  llsm_sinparam* sinu, const char* wtype);

/*
* QHM iterator functions
*/
qhm_solve_status* qhm_status_init(FP_TYPE* x, FP_TYPE* window, FP_TYPE f0, FP_TYPE maxcorr,
  int nx, int fs, int nhar, int lsmethod);
void qhm_status_reset(qhm_solve_status* status);
void delete_qhm_status(qhm_solve_status* status);
int qhm_iter(qhm_solve_status* status);

/*
* QHM result analysis functions.
*/
// work.size == K + 1 * nout
void qhm_synth(double complex* ak, FP_TYPE* fk_hat, FP_TYPE* synthrange, FP_TYPE* out,
  int K, int nout, int fs, double complex* work);
void qhm_synth_half(double complex* ak, FP_TYPE* fk_hat, FP_TYPE* synthrange, FP_TYPE* out,
  int K, int nout, int fs);

void qhm_harmonic_parameter_from_ak(double complex* ak, FP_TYPE* fk_hat, int K, int nhar,
  FP_TYPE* freq, FP_TYPE* ampl, FP_TYPE* phse);

#endif // QHM
