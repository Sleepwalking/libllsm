/*
qhm.c - Quasi-Harmonic Model with Adaptive Iterative Refinement
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

#include "qhm.h"
#include <stdio.h>
#include <math.h>
#include <cblas.h>
#include <lapacke.h>
#include "math-funcs.h"

typedef FP_TYPE*(*window_getter)(int);

static int isshowprogress = 1;

void qhm_progress(int v) {
  isshowprogress = v;
}

void qhm_air(llsm_parameters param, FP_TYPE* x, int nx, int fs, FP_TYPE* f0, int nhop, const char* wtype) {
  double window_periods = 2.5;
  window_getter window_func;
  if(! strcmp(wtype, "blackman_harris")) {
    window_func = blackman_harris;
    window_periods = 2.2;
  }
  else if(! strcmp(wtype, "hamming")) {
    window_func = hamming;
    window_periods = 1.5;
  }
  else if(! strcmp(wtype, "hanning")) {
    window_func = hanning;
    window_periods = 1.5;
  }
  else if(! strcmp(wtype, "blackman")) {
    window_func = blackman;
    window_periods = 1.65;
  }
  else {
    fprintf(stderr, "Unknown window type %s. Aborting...", wtype);
    abort();
  }

  FP_TYPE* s_window = hanning(param.a_nhop * 2);
  FP_TYPE* s_frame = calloc(sizeof(FP_TYPE), param.a_nhop * 2);
  FP_TYPE* s_range = calloc(sizeof(FP_TYPE), param.a_nhop * 2);
  FP_TYPE* srer_work = calloc(sizeof(FP_TYPE), param.a_nhop * 2);

  // F0 adaptive iterative refinement
  int lastoutputed = 0;
  double sumsrer = 0.0;
  int nvoiced = 0;
  for(int ihop = 0; ihop < nhop; ++ ihop) {
    if(f0[ihop] <= 0.0) continue;
    ++nvoiced;

    if(isshowprogress) {
      for(int i = 0; i < lastoutputed; ++ i)
        printf("\b");
      lastoutputed = printf("AIR: %d / %d", ihop, nhop);
      fflush(stdout);
    }

    // get comparison frame
    FP_TYPE* c_frame = fetch_frame(x, nx, ihop * param.a_nhop, param.a_nhop * 2);
    for(int i = 0; i < param.a_nhop * 2; ++ i) {
      s_range[i] = i - param.a_nhop;
      c_frame[i] *= s_window[i];
    }

    // do iterations
    FP_TYPE lastsrer = -INFINITY;
    for(int iiter = 0; iiter < param.a_maxairiter; ++ iiter) {
      // analysis
      int nhar = param.a_mvf / f0[ihop];
      int winlen = (int)(((FP_TYPE)fs) / f0[ihop] * window_periods) * 2 + 1;
      FP_TYPE* a_window = window_func(winlen);
      FP_TYPE* a_frame = fetch_frame(x, nx, ihop * param.a_nhop, winlen);
      qhm_solve_status *status = qhm_status_init(a_frame, a_window, f0[ihop], param.a_maxqhmcorr, winlen, fs, nhar, param.a_qhmlsmethod);
      int info = qhm_iter(status);
      if(info != 0) abort(); // error occurred

      // synthesis and calculate srer
      double complex* s_work = calloc(sizeof(double complex), (status->K + 1) * param.a_nhop * 2);
      qhm_synth(status->lstsqb, status->fk_hat, s_range, s_frame, status->K, param.a_nhop * 2, fs, s_work);
      for(int i = 0; i < param.a_nhop * 2; ++ i)
        s_frame[i] *= s_window[i];
      FP_TYPE currsrer = qhm_calc_srer(c_frame, s_frame, param.a_nhop * 2, srer_work);

      // update f0
      if(currsrer > lastsrer) {
        int nmean = nhar / 3 > 16 ? 16 : nhar / 3;
        FP_TYPE meanf0 = 0.0;
        for(int i = 0; i < nmean; ++ i)
          meanf0 += status->fk_hat[nhar + 1 + i] / (i + 1);
        meanf0 /= nmean;
        f0[ihop] = meanf0;
      }

      // check whether stop iteration
      if(currsrer - lastsrer < 0.1) {
        if(currsrer > lastsrer)
          lastsrer = currsrer;
        free(s_work);
        qhm_status_free(status);
        free(a_frame);
        free(a_window);
        break;
      }

      if(currsrer > lastsrer)
        lastsrer = currsrer;

      free(s_work);
      qhm_status_free(status);
      free(a_frame);
      free(a_window);
    }
    sumsrer += lastsrer;
    free(c_frame);
  }

  if(isshowprogress) {
    for(int i = 0; i < lastoutputed; ++ i)
      printf("\b");
    printf("AIR finished. Average SRER = %lf.\n", sumsrer / nvoiced);
  }

  free(srer_work);
  free(s_range);
  free(s_frame);
  free(s_window);
}

void qhm_analyze(llsm_parameters param, FP_TYPE* x, int nx, int fs, FP_TYPE* f0, int nhop, llsm_sinparam* sinu, const char* wtype) {
  double window_periods = 2.5;
  window_getter window_func;
  if(! strcmp(wtype, "blackman_harris")) {
    window_func = blackman_harris;
    window_periods = 2.2;
  }
  else if(! strcmp(wtype, "hamming")) {
    window_func = hamming;
    window_periods = 1.5;
  }
  else if(! strcmp(wtype, "hanning")) {
    window_func = hanning;
    window_periods = 1.5;
  }
  else if(! strcmp(wtype, "blackman")) {
    window_func = blackman;
    window_periods = 1.65;
  }
  else {
    fprintf(stderr, "Unknown window type %s. Aborting...", wtype);
    abort();
  }

  FP_TYPE* s_window = hanning(param.a_nhop * 2);
  FP_TYPE* s_frame = calloc(sizeof(FP_TYPE), param.a_nhop * 2);
  FP_TYPE* s_range = calloc(sizeof(FP_TYPE), param.a_nhop * 2);
  FP_TYPE* srer_work = calloc(sizeof(FP_TYPE), param.a_nhop * 2);

  int lastoutputed = 0;
  double sumsrer = 0.0;
  int nvoiced = 0;

  for(int ihop = 0; ihop < nhop; ++ ihop) {
    if(f0[ihop] <= 0.0) continue;
    ++nvoiced;

    if(isshowprogress) {
      for(int i = 0; i < lastoutputed; ++ i)
        printf("\b");
      lastoutputed = printf("QHM: %d / %d", ihop, nhop);
      fflush(stdout);
    }

    // get comparison frame and analysis frame
    FP_TYPE* c_frame = fetch_frame(x, nx, ihop * param.a_nhop, param.a_nhop * 2);
    for(int i = 0; i < param.a_nhop * 2; ++ i) {
      s_range[i] = i - param.a_nhop;
      c_frame[i] *= s_window[i];
    }

    int nhar = param.a_mvf / f0[ihop];
    int winlen = (int)(((FP_TYPE)fs) / f0[ihop] * window_periods) * 2 + 1;
    FP_TYPE* a_window = window_func(winlen);
    FP_TYPE* a_frame = fetch_frame(x, nx, ihop * param.a_nhop, winlen);
    qhm_solve_status *status = qhm_status_init(a_frame, a_window, f0[ihop], param.a_maxqhmcorr, winlen, fs, nhar, param.a_qhmlsmethod);
    FP_TYPE bestsrer = -INFINITY;
    double complex* bestak = calloc(sizeof(double complex), status->K);
    FP_TYPE* bestfk_hat = calloc(sizeof(FP_TYPE), status->K);

    // iteration
    for(int iiter = 0; iiter < param.a_maxqhmiter; ++ iiter) {
      // analysis
      int info = qhm_iter(status);
      if(info != 0) abort(); // error occurred

      // synthesis and calculate srer
      qhm_synth_half(status->lstsqb, status->fk_hat, s_range, s_frame, status->K, param.a_nhop * 2, fs);
      for(int i = 0; i < param.a_nhop * 2; ++ i)
        s_frame[i] *= s_window[i];
      FP_TYPE currsrer = qhm_calc_srer(c_frame, s_frame, param.a_nhop * 2, srer_work);
      // check whether update bestak
      if(currsrer > bestsrer) {
        bestsrer = currsrer;
        memcpy(bestfk_hat, status->fk_hat, status->K * sizeof(FP_TYPE));
        memcpy(bestak, status->lstsqb, status->K * sizeof(double complex));
      }
    }

    // export sinusoidal parameters
    if(isinf(bestsrer) || isnan(bestsrer)) {
      fprintf(stderr, "QHM: Failed to analyze at %d. Aborting...\n", ihop);
      abort();
    }
    qhm_harmonic_parameter_from_ak(bestak, bestfk_hat, status->K, param.a_nhar, sinu -> freq[ihop], sinu -> ampl[ihop], sinu -> phse[ihop]);
    sumsrer += bestsrer;

    free(bestfk_hat);
    free(bestak);
    qhm_status_free(status);
    free(a_frame);
    free(a_window);
    free(c_frame);
  }

  if(isshowprogress) {
    for(int i = 0; i < lastoutputed; ++ i)
      printf("\b");
    printf("QHM finished. Average SRER = %lf.\n", sumsrer / nvoiced);
  }

  free(srer_work);
  free(s_frame);
  free(s_range);
  free(s_window);
}

qhm_solve_status* qhm_status_init(FP_TYPE* x, FP_TYPE* window, FP_TYPE f0, FP_TYPE maxcorr, int nx, int fs, int nhar, char lsmethod) {
  qhm_solve_status *status = calloc(sizeof(qhm_solve_status), 1);

  status->x = x;
  status->window = window;
  status->f0 = f0;
  status->maxcorr = maxcorr;
  status->nx = nx;
  status->fs = fs;
  status->nhar = nhar;
  status->ls_method = lsmethod;

  status->K = status->nhar * 2 + 1;
  status->N = (status->nx - 1) / 2;
  status->fk_hat = calloc(sizeof(FP_TYPE), status->K);
  status->n = calloc(sizeof(FP_TYPE), nx);
  status->t = calloc(sizeof(FP_TYPE), status->K * nx);
  status->Ea = calloc(sizeof(double complex), 2 * status->K * nx);
  status->Ew = calloc(sizeof(double complex), nx * 2 * status->K);
  status->R = calloc(sizeof(double complex), 2 * status->K * 2 * status->K);
  status->windowed_x = calloc(sizeof(double complex), nx);
  status->lstsqb = calloc(sizeof(double complex), 2 * status->K);

  qhm_status_reset(status);

  return status;
}

void qhm_status_reset(qhm_solve_status* status) {
  for(int i = -status->nhar; i < status->nhar + 1; ++ i)
    status->fk_hat[i + status->nhar] = ((FP_TYPE)i) * status->f0;
  for(int i = -status->N; i < status->N + 1; ++ i)
    status->n[i + status->N] = (FP_TYPE)i;
}

void qhm_status_free(qhm_solve_status* status) {
  free(status->fk_hat);
  free(status->n);
  free(status->t);
  free(status->Ea);
  free(status->Ew);
  free(status->R);
  free(status->windowed_x);
  free(status->lstsqb);

  free(status);
}

int qhm_iter(qhm_solve_status* status) {
  int m, n, k;
  /* LS Analysis */
  // t = np.dot(2 * np.pi * fk_hat, n) / sr
  for(int row = 0; row < status->K; ++ row) {
    double a = 2.0 * M_PI * status->fk_hat[row] / (FP_TYPE)status->fs;
    for(int col = 0; col < status->nx; ++col)
      status->t[row * status->nx + col] = a * status->n[col];
  }

  // E = np.cos(t) + 1j * np.sin(t), mat with cplx exp
  // Ea = np.concatenate((E, np.tile(n, (K, 1)) * E), axis = 0)
  int base = status->K * status->nx;
  for(int i = 0; i < base; ++ i)
    status->Ea[i] = cos(status->t[i]) + sin(status->t[i]) * I;
  for(int i = 0; i < base; ++ i)
    status->Ea[i + base] = status->Ea[i] * status->n[i % status->nx];
  // Ew = np.tile(window, (1, 2 * K)) * Ea.T # multiply the window
  for(int row = 0; row < status->nx; ++ row)
    for(int col = 0; col < 2 * status->K; ++col)
      status->Ew[row * 2 * status->K + col] = status->window[row] * status->Ea[col * status->nx + row];

  // R = np.dot(Ew.T, Ew), compute the matrix to be inverted
  m = 2 * status->K;
  n = 2 * status->K;
  k = status->nx;
  double complex alpha = 1.0, beta = 0.0;
  cblas_zgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
    m, n, k, &alpha,
    status->Ew, m,
    status->Ew, n,
    &beta, status->R, n
  );
  // lstsqb = np.dot(Ew.T, x * window)
  for(int i = 0; i < status->nx; ++ i)
    status->windowed_x[i] = status->x[i] * status->window[i];
  m = 2 * status->K;
  n = 1;
  k = status->nx;
  cblas_zgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
    m, n, k, &alpha,
    status->Ew, m,
    status->windowed_x, n,
    &beta, status->lstsqb, n
  );

  // theta = sla.lstsq(R, lstsqb)[0]
  int rank;
  m = 2 * status->K;
  n = 2 * status->K;
  int nrhs = 1;
  double s[2 * status->K];
  int info;
  if(status->ls_method == 'Q')
    info = LAPACKE_zgels(LAPACK_ROW_MAJOR, 'N', m, n, nrhs, status->R, n, status->lstsqb, nrhs);
  else if(status->ls_method == 'S')
    info = LAPACKE_zgelsd(LAPACK_ROW_MAJOR, m, n, nrhs, status->R, n, status->lstsqb, nrhs, s, -1.0, &rank);
  else {
    fprintf(stderr, "Invalid ls_method.\n");
    abort();
  }

  if(info > 0) {
    fprintf(stderr, "Least-square calculation failed.\n");
    return info;
  }
  else if(info < 0) {
    fprintf(stderr, "Wrong parameter.\n");
    abort();
  }

  // corr
  // fk_hat += np.clip(sr / (2 * np.pi) * (ak.real * bk.imag - ak.imag * bk.real) / (np.abs(ak) ** 2), -maxCorr, maxCorr)
  base = status->K;
  for(int i = 0; i < status->K; ++i) {
    double absak = cabs(status->lstsqb[i]);
    double delta = status->fs / (2.0 * M_PI) * (creal(status->lstsqb[i]) * cimag(status->lstsqb[base + i]) - cimag(status->lstsqb[i]) * creal(status->lstsqb[base + i])) / (absak * absak);
    if(delta > status->maxcorr) delta = status->maxcorr;
    else if(delta < -status->maxcorr) delta = -status->maxcorr;
    status->fk_hat[i] += delta;
  }

  return 0;
}

void qhm_synth(double complex* ak, FP_TYPE* fk_hat, FP_TYPE* synthrange, FP_TYPE* out, int K, int nout, int fs, double complex* work) {
  int needfree = 0;
  if(work == NULL) {
    work = calloc(sizeof(double complex*), (K + 1) * nout);
    needfree = 1;
  }

  double complex* dotb = work; // (K, nout)
  double complex* outtemp = work + (K * nout); // (nout)

  // dotb = np.exp(np.dot(2j * np.pi * fk_hat, r) / sr)
  for(int row = 0; row < K; ++ row) {
    double complex a = (2 * M_PI * fk_hat[row]) * I;
    for(int col = 0; col < nout; ++ col)
      dotb[row * nout + col] = cexp(a * synthrange[col] / ((double complex)fs));
  }

  // out = np.dot(ak.T, work).real
  int m, n, k;
  m = 1;
  n = nout;
  k = K;
  double complex alpha = 1.0, beta = 0.0;
  cblas_zgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
    m, n, k, &alpha,
    ak, m,
    dotb, n,
    &beta, outtemp, n
  );

  for(int i = 0; i < nout; ++i)
    out[i] = creal(outtemp[i]);

  if(needfree) free(work);
}

FP_TYPE qhm_stdev(FP_TYPE* x, int n) {
  FP_TYPE mean = sumfp(x, n) / n;
  FP_TYPE sum_deviation = 0.0;
  for(int i = 0; i < n; ++ i)
    sum_deviation += (x[i] - mean) * (x[i] - mean);
  return sqrt(sum_deviation / n);
}

FP_TYPE qhm_calc_srer(FP_TYPE* x, FP_TYPE* y, int n, FP_TYPE* work) {
  int needfree = 0;
  if(work == NULL) {
    work = calloc(sizeof(FP_TYPE), n);
    needfree = 1;
  }

  for(int i = 0; i < n; ++ i)
    work[i] = x[i] - y[i];

  FP_TYPE orig_stdev = qhm_stdev(x, n);
  FP_TYPE delta_stdev = qhm_stdev(work, n);

  if(needfree) free(work);

  return log10(orig_stdev / delta_stdev) * 20.0;
}

void qhm_synth_half(double complex* ak, FP_TYPE* fk_hat, FP_TYPE* synthrange, FP_TYPE* out, int K, int nout, int fs) {
  int nhar = (K - 1) / 2;
  for(int i = 0; i < nout; ++ i)
    out[i] = 0.0;
  for(int ihar = 0; ihar < nhar; ++ ihar) {
    double freq = fk_hat[nhar + 1 + ihar];
    double ampl = cabs(ak[nhar + 1 + ihar]) * 2.0;
    double phse = carg(ak[nhar + 1 + ihar]);
    double fph = 2.0 * M_PI / ((double)fs) * freq;
    for(int i = 0; i < nout; ++ i)
      out[i] += cos(fph * ((double)synthrange[i]) + phse) * ampl;
  }
}

void qhm_harmonic_parameter_from_ak(double complex* ak, FP_TYPE* fk_hat, int K, int nhar, FP_TYPE* freq, FP_TYPE* ampl, FP_TYPE* phse)
{
  int akhar = (K - 1) / 2;
  int outhar = akhar < nhar ? akhar : nhar;

  for(int i = 0; i < outhar; ++i) {
    freq[i] = fk_hat[akhar + 1 + i];
    ampl[i] = cabs(ak[akhar + 1 + i]) * 2.0;
    phse[i] = carg(ak[akhar + 1 + i]);
  }

  FP_TYPE* phse_ = unwrap(phse, outhar);
  for(int i = 0; i < outhar; ++i)
    phse[i] = phse_[i];
  free(phse_);
}
