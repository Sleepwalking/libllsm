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

//#define ENABLE_DBGFUNCS

#include "llsm.h"
#include "qhm.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "math-funcs.h"
#include "envelope.h"
#include "external/matlabfunctions.h"

#ifdef ENABLE_DBGFUNCS
static void normalize_write(FP_TYPE* x, int nx, int fs, char* path, double gain) {
  FP_TYPE* normx = calloc(nx, sizeof(FP_TYPE));
  for(int i = 0; i < nx; i ++) normx[i] = x[i] * gain;
  wavwrite(normx, nx, fs, 16, path);
  free(normx);
}
#endif

llsm_parameters llsm_init(int nnosband) {
  llsm_parameters ret;
  ret.a_nhop = 256;
  ret.a_nhar = 80;
  ret.a_nhare = 4;
  ret.a_nnos = 54;
  ret.a_nosf = -1.0;
  ret.a_tfft = 0.04;
  ret.a_mvf = 8000.0;
  ret.a_noswrap = 5000.0;
  ret.a_nnosband = nnosband;
  ret.a_nosbandf = calloc(nnosband - 1, sizeof(FP_TYPE));
  ret.a_nosbandf[0] = 2000;
  ret.a_method = qhm;
  ret.a_maxairiter = 16;
  ret.a_maxqhmiter = 4;
  ret.a_targetsrer = 30.0;
  ret.a_maxqhmcorr = 20.0;
  for(int i = 1; i < nnosband - 1; i ++) ret.a_nosbandf[i] = ret.a_nosbandf[i - 1] * 2.0;
  ret.s_fs = 0;
  return ret;
}

void llsm_deinit(llsm_parameters dst) {
  if(dst.a_nosbandf != NULL)
    free(dst.a_nosbandf);
}

static void spectrogram_analyze(llsm_parameters param, FP_TYPE* x, int nx, int fs, FP_TYPE* f0, int nf0,
  int nfft, FP_TYPE* fftbuff, const char* wtype,
  FP_TYPE** spectrogram, FP_TYPE** phasegram, FP_TYPE** phasegram_d, FP_TYPE* minvec) {

  FP_TYPE* xbuff = calloc(nfft, sizeof(FP_TYPE));
  FP_TYPE* ybuffr = calloc(nfft, sizeof(FP_TYPE));
  FP_TYPE* ybuffi = calloc(nfft, sizeof(FP_TYPE));

  double window_periods = 2.5;
  if(! strcmp(wtype, "blackman_harris"))
  /*
    2.2 periods from both sides works well, though theoretically determined value is 3.0.
    As long as the sinusoids are resolved, shorter is better.
  */
    window_periods = 2.2;
  else if(! strcmp(wtype, "hamming"))
    window_periods = 1.5;
  else if(! strcmp(wtype, "hanning"))
    window_periods = 1.5;
  else if(! strcmp(wtype, "blackman"))
    window_periods = 1.65;

  for(int t = 0; t < nf0; t ++) {
    FP_TYPE resf = f0[t];
    int winlen = resf > 0 ? min(nfft, floor(fs / resf * window_periods) * 2) : param.a_nhop * 2;
    int tn = t * param.a_nhop;

    FP_TYPE* w = NULL;
    if(! strcmp(wtype, "blackman_harris"))
      w = blackman_harris(winlen);
    else if(! strcmp(wtype, "hamming"))
      w = hamming(winlen);
    else if(! strcmp(wtype, "hanning"))
      w = hanning(winlen);
    else if(! strcmp(wtype, "blackman"))
      w = blackman(winlen);
    FP_TYPE norm_factor = 2.0 / sumfp(w, winlen);

    FP_TYPE* xfrm, * xfrm_d;
    FP_TYPE* spec_magn, * spec_phse_, * spec_phse, * spec_phse_d;
    xfrm = xfrm_d = spec_magn = spec_phse_ = spec_phse = spec_phse_d = NULL;

    xfrm = fetch_frame(x, nx, tn, winlen);
    FP_TYPE mean_xfrm = sumfp(xfrm, winlen) / winlen;
    if(minvec != NULL) minvec[t] = minfp(xfrm, winlen);
    for(int i = 0; i < winlen; i ++)
      xfrm[i] -= mean_xfrm;
    if(phasegram_d) {
      xfrm_d = fetch_frame(x, nx, tn - 1, winlen);
      FP_TYPE mean_xfrm_d = sumfp(xfrm_d, winlen) / winlen;
      for(int i = 0; i < winlen; i ++)
        xfrm_d[i] -= mean_xfrm_d;
    }

    for(int i = 0; i < winlen; i ++) {
      xfrm[i] *= w[i];
      if(phasegram_d) xfrm_d[i] *= w[i];
    }

    // current frame
    memset(xbuff, 0, nfft * sizeof(FP_TYPE));
    for(int i = 0; i < winlen / 2; i ++) {
      xbuff[i] = xfrm[i + winlen / 2];
      xbuff[nfft - winlen / 2 + i] = xfrm[i];
    }
    fft(xbuff, NULL, ybuffr, ybuffi, nfft, fftbuff);

    spec_magn = abscplx(ybuffr, ybuffi, nfft / 2);
    if(phasegram) {
      spec_phse_ = argcplx(ybuffr, ybuffi, nfft / 2);
      spec_phse = unwrap(spec_phse_, nfft / 2); free(spec_phse_);
    }

    // delayed frame
    if(phasegram_d) {
      memset(xbuff, 0, nfft * sizeof(FP_TYPE));
      for(int i = 0; i < winlen / 2; i ++) {
        xbuff[i] = xfrm_d[i + winlen / 2];
        xbuff[nfft - winlen / 2 + i] = xfrm_d[i];
      }
      fft(xbuff, NULL, ybuffr, ybuffi, nfft, fftbuff);

      spec_phse_ = argcplx(ybuffr, ybuffi, nfft / 2);
      spec_phse_d = unwrap(spec_phse_, nfft / 2); free(spec_phse_);
    }

    for(int i = 0; i < nfft / 2; i ++) {
      spectrogram[t][i] = log_3(spec_magn[i] * norm_factor);
      if(isnan(spectrogram[t][i]) || isinf(spectrogram[t][i]))
        spectrogram[t][i] = -100.0;
      if(phasegram) phasegram[t][i] = spec_phse[i];
      if(phasegram_d) phasegram_d[t][i] = spec_phse_d[i];
    }

    free(spec_magn);
    if(spec_phse) free(spec_phse);
    if(spec_phse_d) free(spec_phse_d);
    free(xfrm);
    if(xfrm_d) free(xfrm_d);
    free(w);
  }

  free(xbuff);
  free(ybuffr);
  free(ybuffi);
}

static FP_TYPE* refine_f0(llsm_parameters param, int nfft, int fs, FP_TYPE* f0, int nf0,
  FP_TYPE** spectrogram, FP_TYPE** phasegram, FP_TYPE** phasegram_d) {

  FP_TYPE* rf0 = calloc(nf0, sizeof(FP_TYPE));
  for(int t = 0; t < nf0; t ++) {
    if(f0[t] <= 0) {
      rf0[t] = 0;
      continue;
    }
    int i_l = floor(f0[t] / fs * nfft * 0.7);
    int i_h = ceil(f0[t] / fs * nfft * 1.3);
    int i_max = round(f0[t] / fs * nfft);
    FP_TYPE peak = spectrogram[t][i_max];

    for(int i = i_l; i <= i_h; i ++)
      if(spectrogram[t][i] > peak) {
        peak = spectrogram[t][i];
        i_max = i;
      }

    FP_TYPE p   = phasegram  [t][i_max];
    FP_TYPE p_d = phasegram_d[t][i_max];
    p   -= floor(p   / 2.0 / M_PI) * 2.0 * M_PI;
    p_d -= floor(p_d / 2.0 / M_PI) * 2.0 * M_PI;
    if(p < p_d)
      p += 2.0 * M_PI;

    rf0[t] = (p - p_d) / 2.0 / M_PI * fs;
    if(fabs(rf0[t] / fs * nfft - i_max) > 1 || fabs(rf0[t] - f0[t]) > 10)
      rf0[t] = f0[t];
  }

  return rf0;
}

static FP_TYPE find_min_left(FP_TYPE* arr, int idx) {
  FP_TYPE ret = arr[idx];
  int left = idx - 100 > 0 ? idx - 100 : 0;
  while(idx > left && arr[idx - 1] <= arr[idx])
      ret = arr[-- idx];
  return ret;
}

static FP_TYPE find_min_right(FP_TYPE* arr, int idx, int length) {
  FP_TYPE ret = arr[idx];
  int right = idx + 100 > length ? length : idx + 100;
  while(idx < right - 1 && arr[idx + 1] <= arr[idx])
      ret = arr[++ idx];
  return ret;
}

int find_peak(FP_TYPE* arr, int l_idx, int u_idx, int direction) {
  FP_TYPE max = arr[l_idx] * direction;
  FP_TYPE max_idx = l_idx;
  for(int i = l_idx; i < u_idx; i ++)
    if(arr[i] * direction > arr[i - 1] * direction && arr[i] * direction > arr[i + 1] * direction)
      if(arr[i] * direction > max) {
        max = arr[i] * direction;
        max_idx = i;
      }
  return max_idx;
}

static int find_peak_cand(FP_TYPE* dst_cand, FP_TYPE* spectrum, int length,
  int l_idx, int u_idx, FP_TYPE threshold, int maxnum) {
  int n = 0;
  for(int i = l_idx; i <= u_idx; i ++)
    if(spectrum[i] >= spectrum[i - 1] && spectrum[i] > spectrum[i + 1]) {
      FP_TYPE l_height = spectrum[i] - find_min_left(spectrum, i);
      FP_TYPE r_height = spectrum[i] - find_min_right(spectrum, i, length);
      if(fmin(l_height, r_height) > threshold) {
        if(n == maxnum) { // in case that the number of peaks exceeds the buffer capacity
          n ++;
          break;
        } else
          dst_cand[n ++] = i;
      }
    }
  if(n > maxnum) // raise the threshold to select more prominent peaks
    return find_peak_cand(dst_cand, spectrum, length, l_idx, u_idx, threshold * 1.02, maxnum);
  return n;
}

int select_nearest_cand(FP_TYPE* cand_list, int cand_n, FP_TYPE center) {
  FP_TYPE dist = fabs(cand_list[0] - center);
  int ret = 0;
  for(int i = 0; i < cand_n; i ++)
    if(fabs(cand_list[i] - center) < dist) {
      dist = fabs(cand_list[i] - center);
      ret = i;
    }
  return ret;
}

void interp_peak(FP_TYPE* dst_freq, FP_TYPE* dst_ampl, FP_TYPE* spectrum, int n) {
  // quadratic interpolation
  FP_TYPE a, b, c, a1, a2, x;
  a = spectrum[n - 1];
  b = spectrum[n + 0];
  c = spectrum[n + 1];
  a1 = (a + c) / 2.0 - b;
  a2 = c - b - a1;
  x = - a2 / a1 * 0.5;

  x = (fabs(x) < 1.0) ? x : 0; // in case we get some x outside of [n-1, n+1]

  *dst_freq = (FP_TYPE)n + x;
  *dst_ampl = a1 * x * x + a2 * x + b;
  *dst_ampl = *dst_ampl > b + 0.1 ? b + 0.1 : *dst_ampl; // in prevention of numerical instability
}

static void find_harmonic_trajectory(llsm_parameters param, FP_TYPE** spectrogram, FP_TYPE** phasegram,
  int nfft, int fs, FP_TYPE* f0, int nf0, int idx, FP_TYPE* freq, FP_TYPE* ampl, FP_TYPE* phse) {
  const FP_TYPE hardev = 0.3;
  FP_TYPE cand_bin[20];
  int cand_bin_n = 0;
  for(int i = 0; i < nf0; i ++) {
    freq[i] = 0;
    ampl[i] = 0;
    phse[i] = 0;
    if(f0[i] < 50 || f0[i] * idx > param.a_mvf)
      continue;

    int l_idx = round(f0[i] * (idx - hardev) / fs * nfft);
    int u_idx = round(f0[i] * (idx + hardev) / fs * nfft);
    l_idx = l_idx < 1 ? 1 : l_idx;
    u_idx = u_idx > nfft / 2 - 1 ? nfft / 2 - 1 : u_idx;

    cand_bin_n = find_peak_cand(cand_bin, spectrogram[i], nfft / 2, l_idx, u_idx, 0.0, 20);
    int peak_bin;
    if(cand_bin_n > 0)
      peak_bin = cand_bin[select_nearest_cand(cand_bin, cand_bin_n, f0[i] * idx / fs * nfft)];
    else
      peak_bin = find_peak(spectrogram[i], l_idx, u_idx, 1);

    FP_TYPE peak_freq, peak_ampl;
    interp_peak(& peak_freq, & peak_ampl, spectrogram[i], peak_bin);
    freq[i] = peak_freq * fs / nfft;
    ampl[i] = exp_3(peak_ampl);
    phse[i] = linterp(phasegram[i][(int)peak_freq], phasegram[i][(int)peak_freq + 1], fmod(peak_freq, 1.0));

    if(isnan(ampl[i])) // peak amplitude goes nan if one of the bins = -INF
      ampl[i] = 0;
  }
}

static FP_TYPE* synth_sinusoid_frame(FP_TYPE* freq, FP_TYPE* ampl, FP_TYPE* phse, int nhar, int fs, int length) {
  FP_TYPE* x = calloc(length, sizeof(FP_TYPE));
  FP_TYPE* s = calloc(length, sizeof(FP_TYPE));
  for(int i = 0; i < nhar; i ++) {
    FP_TYPE h_f = freq[i];
    FP_TYPE h_a = ampl[i];
    FP_TYPE h_p = phse[i];
    if(i > 0 && h_f < 50.0) // reaching the end
      break;

    FP_TYPE tphffs = 2.0 * M_PI / fs * h_f;
    FP_TYPE c = 2.0 * cos_3(tphffs);
    s[0] = cos_3(tphffs * (- length / 2) + h_p) * h_a;
    s[1] = cos_3(tphffs * (- length / 2 + 1) + h_p) * h_a;
    x[0] += s[0];
    x[1] += s[1];
    for(int t = 2; t < length; t ++) {
      s[t] = c * s[t - 1] - s[t - 2];
      x[t] += s[t];
    }
  }
  free(s);
  return x;
}

static FP_TYPE* filter_noise(llsm_parameters param, llsm_conf conf, FP_TYPE* x, int nx,
  FP_TYPE** wrapped_spectrogram, int winsize) {
/*
  To suppress glitches caused by aliasing, we apply MLT sine window twice, before analysis and after synthesis respectively;
  Overlapping factor should be greater or equal to 4 for MLT sine window;
  To preserve time resolution, we actually lower the hopsize by half, meanwhile double sampling wrapped_spectrogram.
*/
  conf.nhop /= 2;
  int nfft = winsize * conf.nhop * 2;
  int ny = conf.nfrm * conf.nhop * 2 + nfft * 16;
  FP_TYPE* y = calloc(ny, sizeof(FP_TYPE));
  FP_TYPE* w = hanning(nfft);
  FP_TYPE* realbuff = calloc(nfft, sizeof(FP_TYPE));
  FP_TYPE* imagbuff = calloc(nfft, sizeof(FP_TYPE));
  FP_TYPE* yfrm = calloc(nfft, sizeof(FP_TYPE));
  FP_TYPE fftbuff[65536];
  FP_TYPE norm_factor = 0.5 * sumfp(w, nfft);
  FP_TYPE norm_factor_win = 0;
  {
    int i = 0;
    while(i < nfft) {
      norm_factor_win += w[i];
      i += conf.nhop;
    }
  }
  for(int j = 0; j < nfft; j ++)
    w[j] = sqrt(w[j]);

  FP_TYPE* freqwrap = llsm_wrap_freq(0, conf.nosf, conf.nnos, conf.noswrap);
  for(int i = 0; i < conf.nfrm * 2; i ++) {
    int t = i * conf.nhop;
    FP_TYPE* spec = llsm_spectrum_from_envelope(freqwrap, wrapped_spectrogram[i / 2], conf.nnos, nfft / 2, param.s_fs);
    FP_TYPE* xfrm = fetch_frame(x, nx, t, nfft);
    for(int j = 0; j < nfft; j ++)
      xfrm[j] *= w[j];
    fft(xfrm, NULL, realbuff, imagbuff, nfft, fftbuff);
    for(int j = 0; j < nfft / 2; j ++) {
      FP_TYPE a = exp_2(spec[j]) * norm_factor; // amplitude
      FP_TYPE p = atan2_2(imagbuff[j], realbuff[j]);
      realbuff[j] = a * cos_1(p);
      imagbuff[j] = a * sin_1(p);
    }
    complete_symm (realbuff, nfft);
    complete_asymm(imagbuff, nfft);
    ifft(realbuff, imagbuff, yfrm, NULL, nfft, fftbuff);
    for(int j = 0; j < nfft; j ++) {
      int idx = t + j - nfft / 2;
      if(idx >= 0)
        y[idx] += yfrm[j] * w[j] / norm_factor_win;
    }
    free(spec);
    free(xfrm);
  }

  free(w);
  free(yfrm);
  free(realbuff);
  free(imagbuff);
  free(freqwrap);
  return y;
}

static int* gen_ps_instants(FP_TYPE* f0, int nhop, int nf0, int nx, int fs, int* ni) {
  FP_TYPE t = 0;
  *ni = 0;
  int Ni = 100;
  int* instants = calloc(100, sizeof(int));
  while(t < (double)nx / fs) {
    int i = t * fs / nhop;
    instants[(*ni) ++] = t * fs;
    if((*ni) * 2 > Ni) {
      Ni *= 2;
      instants = realloc(instants, Ni * sizeof(int));
    }
    t += 1.0 / max(100, f0[i >= nf0 ? nf0 - 1 : i]);
  }
  return instants;
}

static void subtract_minimum_envelope(FP_TYPE* x, int nx, FP_TYPE* f0, int nhop, int nfrm, int fs) {
  int ninstant = 0;
  int* env_instants = gen_ps_instants(f0, nhop, nfrm, nx, fs, & ninstant);
  FP_TYPE* fenv_instants = calloc(ninstant, sizeof(FP_TYPE));
  int* env_winlen = calloc(ninstant, sizeof(int));
  for(int i = 0; i < ninstant; i ++) env_winlen[i] = fs / max(100, f0[i >= nfrm ? nfrm - 1 : i]) * 1.5;
  FP_TYPE* env_samples = llsm_nonuniform_envelope(x, nx, env_instants, env_winlen, ninstant, 0);
  FP_TYPE* iota = calloc(nx, sizeof(FP_TYPE));
  for(int i = 0; i < nx; i ++) iota[i] = i;
  for(int i = 0; i < ninstant; i ++) fenv_instants[i] = env_instants[i];
  FP_TYPE* env = interp1(fenv_instants, env_samples, ninstant, iota, nx);
  for(int i = 0; i < nx; i ++) x[i] -= env[i];

  free(env_instants);
  free(fenv_instants);
  free(env_winlen);
  free(env_samples);
  free(iota);
  free(env);
}

static FP_TYPE* stretch_static_noise(FP_TYPE* x, int nx, int ny, int novl) {
  FP_TYPE* y = calloc(ny, sizeof(FP_TYPE));
  for(int i = 0; i < min(nx, ny); i ++)
    y[i] = x[i];
  if(ny <= nx) return y;

  int top = nx;
  while(1) {
    for(int i = 0; i < novl; i ++) {
      y[top - novl + i] *= 1.0 - (FP_TYPE)i / novl;
      y[top - novl + i] += x[i] * (FP_TYPE)i / novl;
    }
    for(int i = 0; i < nx - novl; i ++) {
      if(top + i >= ny) return y;
      y[top + i] = x[i + novl];
    }
    top += nx - novl;
  }

  return y;
}

void llsm_delete(llsm* model) {
  if(model == NULL) return;
  int nfrm = model -> conf.nfrm;
  if(model -> conf.nosbandf != NULL) free(model -> conf.nosbandf);
  free(model -> f0);
  free2d(model -> noise, nfrm);
  free2d(model -> sinu -> freq, nfrm);
  free2d(model -> sinu -> ampl, nfrm);
  free2d(model -> sinu -> phse, nfrm);
  free(model -> sinu);

  for(int i = 0; i < model -> conf.nnosband; i ++) {
    free2d(model -> nosch[i] -> eenv -> freq, nfrm);
    free2d(model -> nosch[i] -> eenv -> ampl, nfrm);
    free2d(model -> nosch[i] -> eenv -> phse, nfrm);
    free(model -> nosch[i] -> eenv);
    if(model -> nosch[i] -> emin != NULL) free(model -> nosch[i] -> emin);
    free(model -> nosch[i]);
  }
  free(model -> nosch);
  free(model);
}

llsm* llsm_analyze(llsm_parameters param, FP_TYPE* x, int nx, int fs, FP_TYPE* f0, int nf0) {
  llsm* model = malloc(sizeof(llsm));
  model -> sinu = malloc(sizeof(llsm_sinparam));
  int nfft = pow(2, ceil(log2(fs * param.a_tfft)));
  double nosf = param.a_nosf > 0.0 ? param.a_nosf : fs / 2.0;
  FP_TYPE fftbuff[65536];

  // C2
  FP_TYPE** spectrogram = (FP_TYPE**)malloc2d(nf0, nfft / 2, sizeof(FP_TYPE));
  FP_TYPE** phasegram   = (FP_TYPE**)malloc2d(nf0, nfft / 2, sizeof(FP_TYPE));
  FP_TYPE** phasegram_d = (FP_TYPE**)malloc2d(nf0, nfft / 2, sizeof(FP_TYPE));

  spectrogram_analyze(param, x, nx, fs, f0, nf0, nfft, fftbuff, "blackman",
    spectrogram, phasegram, phasegram_d, NULL);

  // C3
  model -> f0 = refine_f0(param, nfft, fs, f0, nf0, spectrogram, phasegram, phasegram_d);
  FP_TYPE* rf0 = model -> f0;

  // C4
  model -> sinu -> freq = (FP_TYPE**)malloc2d(nf0, param.a_nhar, sizeof(FP_TYPE));
  model -> sinu -> ampl = (FP_TYPE**)malloc2d(nf0, param.a_nhar, sizeof(FP_TYPE));
  model -> sinu -> phse = (FP_TYPE**)malloc2d(nf0, param.a_nhar, sizeof(FP_TYPE));
  model -> sinu -> nfrm = nf0;
  model -> sinu -> nhar = param.a_nhar;
  FP_TYPE* tmpfreq = calloc(nf0, sizeof(FP_TYPE));
  FP_TYPE* tmpampl = calloc(nf0, sizeof(FP_TYPE));
  FP_TYPE* tmpphse = calloc(nf0, sizeof(FP_TYPE));

  if(param.a_method == qfft) {
    for(int i = 0; i < param.a_nhar; i ++) {
      find_harmonic_trajectory(param, spectrogram, phasegram, nfft, fs, rf0, nf0, i + 1, tmpfreq, tmpampl, tmpphse);
      for(int j = 0; j < nf0; j ++) {
        model -> sinu -> freq[j][i] = tmpfreq[j];
        model -> sinu -> ampl[j][i] = tmpampl[j];
        model -> sinu -> phse[j][i] = tmpphse[j];
      }
      if(i == 5)
        for(int j = 0; j < nf0; j ++) {
          FP_TYPE avg_f0 = 0;
          for(int k = 0; k < i; k ++)
            avg_f0 += model -> sinu -> freq[j][k] / (k + 1.0);
          avg_f0 /= i;
          if(fabs(rf0[j] - avg_f0) > 1)
            rf0[j] = avg_f0;
        }
    }
  }
  else if(param.a_method == qhm) {
    qhm_air(param, x, nx, fs, rf0, nf0, "blackman");
    qhm_analyze(param, x, nx, fs, rf0, nf0, model -> sinu, "blackman");
  }
  else {
    fprintf(stderr, "Invalid parameter->a_method. Aborting...\n");
    abort();
  }

  free2d(spectrogram, nf0);
  free2d(phasegram, nf0);
  free2d(phasegram_d, nf0);

  // C6
  FP_TYPE* resynth = calloc(nx + param.a_nhop, sizeof(FP_TYPE));
  int nresynth = 2 * param.a_nhop;
  FP_TYPE* resynth_window = hanning(nresynth);
  for(int i = 0; i < nf0; i ++) {
    int tn = i * param.a_nhop;
    FP_TYPE* resynth_frame = synth_sinusoid_frame(
      model -> sinu -> freq[i], model -> sinu -> ampl[i], model -> sinu -> phse[i],
      param.a_nhar, fs, nresynth);
    for(int j = 0; j < nresynth; j ++)
      if(tn + j - nresynth / 2 > 0)
        resynth[tn + j - nresynth / 2] += resynth_frame[j] * resynth_window[j];
    free(resynth_frame);
  }
  free(resynth_window);

  // C7
  for(int i = 0; i < nx; i ++)
    resynth[i] = x[i] - resynth[i];
  wavwrite(resynth, nx, fs, 16, "noise.wav");

  // C21
  FP_TYPE** noise_spectrogram = (FP_TYPE**)malloc2d(nf0, nfft / 2, sizeof(FP_TYPE));
  model -> noise = (FP_TYPE**)calloc(nf0, sizeof(FP_TYPE*));
  spectrogram_analyze(param, resynth, nx, fs, f0, nf0, nfft, fftbuff, "hanning",
    noise_spectrogram, NULL, NULL, NULL);

  // C22
  FP_TYPE* freqwrap = llsm_wrap_freq(0, nosf, param.a_nnos, param.a_noswrap);
  for(int t = 0; t < nf0; t ++) {
  /*
    // cut out the spectrum content around harmonics
    FP_TYPE resf = f0[t];
    int winlen = min(nfft, floor(fs / resf * 2.5) * 2);
    int wmlobe = ceil(1.3 * nfft / winlen);
    for(int i = 0; i < param.a_nhar; i ++) {
      FP_TYPE centerf = model -> sinu -> freq[t][i];
      if(centerf > fs / nfft) { // make sure the frequency is valid
        int lbin = max(0           , round(centerf / fs * nfft - wmlobe));
        int hbin = min(nfft / 2 - 1, round(centerf / fs * nfft + wmlobe));
        for(int j = lbin; j <= hbin; j ++)
          noise_spectrogram[t][j] = -30;
      }
    }*/
    model -> noise[t] = llsm_geometric_envelope(noise_spectrogram[t], nfft / 2, fs, freqwrap, param.a_nnos);
  }

  // remove the content above nyquist frequency
  for(int i = 0; i < param.a_nnos; i ++)
    if(freqwrap[i] >= fs / 2)
      for(int t = 0; t < nf0; t ++)
        model -> noise[t][i] = -100;

  free(freqwrap);
  free2d(noise_spectrogram, nf0);

  // for each noise channel
  model -> nosch = calloc(param.a_nnosband, sizeof(llsm_echannel*));
  for(int b = 0; b < param.a_nnosband; b ++) {
    model -> nosch[b] = malloc(sizeof(llsm_echannel));
    // CC2
    int filtord = LLSM_CHEBY_ORDER * 2;
    FP_TYPE* b_filtered = NULL;
    FP_TYPE favg = 0;
    if(b == 0) {
      b_filtered = chebyfilt(resynth, nx, param.a_nosbandf[b] / fs * 2.0, 0, "lowpass");
      favg = param.a_nosbandf[b] / 2.0;
      filtord -= LLSM_CHEBY_ORDER;
    } else
    if(b == param.a_nnosband - 1) {
      b_filtered = chebyfilt(resynth, nx, param.a_nosbandf[b - 1] / fs * 2.0, 0, "highpass");
      favg = (param.a_nosbandf[b - 1] + fs) / 2.0;
      filtord -= LLSM_CHEBY_ORDER;
    } else {
      b_filtered = chebyfilt(resynth, nx, param.a_nosbandf[b - 1] / fs * 2.0,
        param.a_nosbandf[b] / fs * 2.0, "bandpass");
      favg = (param.a_nosbandf[b - 1] + param.a_nosbandf[b]) / 2.0;
    }

    // CC3
    for(int i = 0; i < nx + LLSM_CHEBY_ORDER - 1; i ++)
      b_filtered[i] = b_filtered[i] * b_filtered[i];
    int mavgord = round(fs / favg * 2);
    FP_TYPE* b_env = moving_avg(b_filtered, nx, mavgord);
    for(int i = 0; i < nx; i ++)
      b_env[i] = sqrt(b_env[i]);
    free(b_filtered);

    // CC5
    FP_TYPE** b_spectrogram = (FP_TYPE**)malloc2d(nf0, nfft / 2, sizeof(FP_TYPE));
    FP_TYPE** b_phasegram   = (FP_TYPE**)malloc2d(nf0, nfft / 2, sizeof(FP_TYPE));
    model -> nosch[b] -> emin = calloc(nf0, sizeof(FP_TYPE));
    spectrogram_analyze(param, b_env + mavgord / 2, nx, fs, f0, nf0, nfft, fftbuff, "blackman",
      b_spectrogram, b_phasegram, NULL, model -> nosch[b] -> emin);
    free(b_env);

    // CC6
    model -> nosch[b] -> eenv = malloc(sizeof(llsm_sinparam));
    llsm_sinparam* b_eenv = model -> nosch[b] -> eenv;
    b_eenv -> nfrm = nf0;
    b_eenv -> nhar = param.a_nhare;
    b_eenv -> freq = (FP_TYPE**)malloc2d(nf0, param.a_nhar, sizeof(FP_TYPE));
    b_eenv -> ampl = (FP_TYPE**)malloc2d(nf0, param.a_nhar, sizeof(FP_TYPE));
    b_eenv -> phse = (FP_TYPE**)malloc2d(nf0, param.a_nhar, sizeof(FP_TYPE));
    for(int i = 0; i < param.a_nhare; i ++) {
      find_harmonic_trajectory(param, b_spectrogram, b_phasegram, nfft, fs, rf0, nf0, i + 1, tmpfreq, tmpampl, tmpphse);
      for(int j = 0; j < nf0; j ++) {
        b_eenv -> freq[j][i] = tmpfreq[j];
        b_eenv -> ampl[j][i] = tmpampl[j];
        b_eenv -> phse[j][i] = tmpphse[j];
      }
    }
    free2d(b_spectrogram, nf0);
    free2d(b_phasegram, nf0);
  }
  free(resynth);
  free(tmpfreq);
  free(tmpampl);
  free(tmpphse);

  // C8
  for(int i = 0; i < nf0; i ++) {
    FP_TYPE base = model -> sinu -> phse[i][0];
    if(rf0[i] <= 0.0) continue;
    for(int j = 0; j < param.a_nhar; j ++) {
      model -> sinu -> phse[i][j] -= base * model -> sinu -> freq[i][j] / rf0[i];
      model -> sinu -> phse[i][j] = fmod(model -> sinu -> phse[i][j] + 1001.0 * M_PI, M_PI * 2.0) - M_PI;
    }
    for(int b = 0; b < param.a_nnosband; b ++)
      for(int j = 0; j < param.a_nhare; j ++) {
        model -> nosch[b] -> eenv -> phse[i][j] -= base * model -> nosch[b] -> eenv -> freq[i][j] / rf0[i];
        model -> nosch[b] -> eenv -> phse[i][j] = fmod(model -> nosch[b] -> eenv -> phse[i][j] + 1001.0 * M_PI, M_PI * 2.0) - M_PI;
      }
  }

  model -> sinu -> nfrm = nf0;
  model -> conf.nfrm = nf0;
  model -> conf.nhop = param.a_nhop;
  model -> conf.thop = (FP_TYPE)param.a_nhop / fs;
  model -> conf.nhar = param.a_nhar;
  model -> conf.nhare = param.a_nhare;
  model -> conf.nnos = param.a_nnos;
  model -> conf.nosf = nosf;
  model -> conf.mvf = param.a_mvf;
  model -> conf.noswrap = param.a_noswrap;
  model -> conf.nosbandf = calloc(param.a_nnosband - 1, sizeof(FP_TYPE));
  model -> conf.nnosband = param.a_nnosband;
  for(int i = 0; i < param.a_nnosband - 1; i ++)
    model -> conf.nosbandf[i] = param.a_nosbandf[i];

  return model;
}

FP_TYPE* llsm_synthesize(llsm_parameters param, llsm* model, int* ny) {
  int nfrm = model -> conf.nfrm;
  int nhop = model -> conf.nhop;
  int fs = param.s_fs;
  int nfft = nhop * 4;
  int ola_factor = 4 / 2;

  *ny = nfrm * nhop + nfft;
  FP_TYPE* y_sin = calloc(*ny, sizeof(FP_TYPE));

  // D1
  FP_TYPE** sin_phse = (FP_TYPE**)malloc2d(nfrm * ola_factor, model -> conf.nhar, sizeof(FP_TYPE));
  FP_TYPE* sin_phse_sync = calloc(nfrm * ola_factor, sizeof(FP_TYPE));
  FP_TYPE phse0 = 0;
  for(int i = 1; i < nfrm * ola_factor; i ++) {
    FP_TYPE f0 = model -> f0[i / ola_factor];
    phse0 += f0 * nhop / ola_factor / fs * 2.0 * M_PI;
    sin_phse_sync[i] = fmod(phse0, 2.0 * M_PI) - M_PI;
    for(int j = 0; j < model -> conf.nhar; j ++)
      sin_phse[i][j] = sin_phse_sync[i] / f0 * model -> sinu -> freq[i / ola_factor][j]
        + model -> sinu -> phse[i / ola_factor][j];
  }

  // D2
  FP_TYPE* ola_window = hanning(nhop * 2);
# pragma omp parallel for
  for(int i = 0; i < nfrm * ola_factor; i ++) {
    int tn = i * nhop / ola_factor;
    if(model -> f0[i / ola_factor] <= 0.0) continue;

    FP_TYPE* sin_frame = synth_sinusoid_frame(
      model -> sinu -> freq[i / ola_factor], model -> sinu -> ampl[i / ola_factor], sin_phse[i],
      model -> conf.nhar, fs, nhop * 2);
#   pragma omp critical
    for(int j = 0; j < nhop * 2; j ++)
      if(tn + j - nhop > 0)
        y_sin[tn + j - nhop] += sin_frame[j] * ola_window[j] / ola_factor;
    free(sin_frame);
  }

  // DC1
  FP_TYPE* s = white_noise(1.0, fs); // one-second-long noise template
  FP_TYPE* noise_excitation = calloc(*ny, sizeof(FP_TYPE));

  // for each noise channel
# pragma omp parallel for
  for(int b = 0; b < model -> conf.nnosband; b ++) {
    llsm_echannel* b_channel = model -> nosch[b];
    FP_TYPE* b_env = calloc(*ny, sizeof(FP_TYPE));
    FP_TYPE* b_env_mix = calloc(*ny, sizeof(FP_TYPE));

    // D1
    FP_TYPE** b_phse = (FP_TYPE**)copy2d(b_channel -> eenv -> phse, nfrm, model -> conf.nhare, sizeof(FP_TYPE));
    for(int i = 0; i < nfrm; i ++) {
      FP_TYPE f0 = model -> f0[i];
      for(int j = 0; j < model -> conf.nhare; j ++)
        b_phse[i][j] = sin_phse_sync[i * ola_factor] / f0 * b_channel -> eenv -> freq[i][j] + b_channel -> eenv -> phse[i][j];
    }

    // D3
    for(int i = 0; i < nfrm; i ++) {
      int tn = i * nhop;
      if(model -> f0[i] <= 0.0) continue;
      FP_TYPE* b_env_frame = synth_sinusoid_frame(
        b_channel -> eenv -> freq[i], b_channel -> eenv -> ampl[i], b_phse[i],
        model -> conf.nhare, fs, nhop * 2);
      for(int j = 0; j < nhop * 2; j ++)
        if(tn + j - nhop > 0)
          b_env[tn + j - nhop] += b_env_frame[j] * ola_window[j];
      free(b_env_frame);
    }
    free2d(b_phse, nfrm);

    // DC2
    int filtord = LLSM_CHEBY_ORDER * 2;
    FP_TYPE* b_template = NULL;
    FP_TYPE* b_filtered = NULL;
    if(b == 0) {
      b_template = chebyfilt(s, fs, model -> conf.nosbandf[b] / fs * 2.0, 0, "lowpass");
      filtord -= LLSM_CHEBY_ORDER;
    } else
    if(b == model -> conf.nnosband - 1)
      b_template = chebyfilt(s, fs, model -> conf.nosbandf[b - 1] / fs * 2.0,
        model -> conf.nosf / fs * 2.0, "bandpass");
    else
      b_template = chebyfilt(s, fs, model -> conf.nosbandf[b - 1] / fs * 2.0,
        model -> conf.nosbandf[b] / fs * 2.0, "bandpass");
    b_filtered = stretch_static_noise(b_template, fs, *ny, 256);
    free(b_template);

    // DC3
    subtract_minimum_envelope(b_env, *ny, model -> f0, nhop, nfrm, fs);

    FP_TYPE* b_normalized = calloc(*ny, sizeof(FP_TYPE));
    for(int i = 0; i < nfrm; i ++) {
      FP_TYPE* hfrm = fetch_frame(b_filtered, *ny, i * nhop, nhop * 2);
      FP_TYPE* efrm = fetch_frame(b_env, *ny, i * nhop, nhop * 2);
      FP_TYPE havg = 0;
      for(int j = 0; j < nhop * 2; j ++)
        havg += hfrm[j] * hfrm[j];
      havg /= nhop * 2;
      for(int j = 0; j < nhop * 2; j ++)
        hfrm[j] *= sqrt(1.0 / (havg + EPS));

      for(int j = 0; j < nhop * 2; j ++)
        efrm[j] += b_channel -> emin[i];

      for(int j = 0; j < nhop * 2; j ++)
        if(i * nhop + j - nhop > 0) {
          b_normalized[i * nhop + j - nhop] += hfrm[j] * ola_window[j];
          b_env_mix   [i * nhop + j - nhop] += efrm[j] * ola_window[j];
        }

      free(hfrm);
      free(efrm);
    }
    /*
    char* writepath = strdup("/tmp/s .wav");
    writepath[6] = '0' + b;
    normalize_write(b_env_mix, *ny, fs, writepath, 1);
    free(writepath);
    */
    for(int i = 0; i < *ny; i ++)
      b_env_mix[i] = b_env_mix[i];

    // DC4
#   pragma omp critical
    for(int i = 0; i < *ny; i ++)
      noise_excitation[i] += b_normalized[i] * b_env_mix[i];

    free(b_filtered);
    free(b_normalized);
    free(b_env_mix);
    free(b_env);
  }
  free(sin_phse_sync);
  free2d(sin_phse, nfrm * ola_factor);
  free(s);

  // DC5
  FP_TYPE* y_nos = filter_noise(param, model -> conf, noise_excitation, *ny, model -> noise, 2);

  // DB2
  for(int i = 0; i < *ny; i ++)
    y_sin[i] += y_nos[i];

  free(noise_excitation);
  free(y_nos);

  free(ola_window);
  return y_sin;
}
