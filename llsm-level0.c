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

//#define DEBUG_EXPORT_WAV

#include "llsm.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "envelope.h"
#include "math-funcs.h"

#ifdef DEBUG_EXPORT_WAV
static void gain_write(FP_TYPE* x, int nx, int fs, char* path, FP_TYPE gain) {
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
  ret.a_f0refine = 1;
  ret.a_wsize = 2.2;
  ret.a_tfft = 0.04;
  ret.a_mvf = 8000.0;
  ret.a_noswarp = 15000.0;
  ret.a_nnosband = nnosband;
  ret.a_nosbandf = calloc(nnosband - 1, sizeof(FP_TYPE));
  ret.a_nosbandf[0] = 2000;
  for(int i = 1; i < nnosband - 1; i ++) ret.a_nosbandf[i] = ret.a_nosbandf[i - 1] * 2.0;
  ret.s_fs = 0;
  ret.s_noiseonly = 0;
  ret.s_n0 = 0;
  ret.s_n1 = 0;
  return ret;
}

void llsm_deinit(llsm_parameters dst) {
  if(dst.a_nosbandf != NULL)
    free(dst.a_nosbandf);
}

// compute magnitude, phase and delayed phase spectrogram on input signal
void spectrogram_analyze(FP_TYPE* x, int nx, int* t, int nfrm, int* winsize, int nfft, FP_TYPE fs,
  const char* wintype, FP_TYPE** spectrogram, FP_TYPE** phasegram, FP_TYPE** phasegram_d, FP_TYPE* dc) {
  
# pragma omp parallel for
  for(int n = 0; n < nfrm; n ++) {
    FP_TYPE fftbuff[65536];
    FP_TYPE inoutbuff[65536];
    FP_TYPE* xbuff = inoutbuff + 0;
    FP_TYPE* ybuffr = inoutbuff + nfft;
    FP_TYPE* ybuffi = inoutbuff + nfft * 2;
    int winlen = min(nfft, winsize[n]);
    int tn = t[n];
    
    FP_TYPE* w = NULL;
    if(! strcmp(wintype, "blackman_harris"))
      w = blackman_harris(winlen);
    else if(! strcmp(wintype, "hamming"))
      w = hamming(winlen);
    else if(! strcmp(wintype, "hanning"))
      w = hanning(winlen);
    else if(! strcmp(wintype, "mltsine"))
      w = mltsine(winlen);
    else if(! strcmp(wintype, "blackman"))
      w = blackman(winlen);
    FP_TYPE norm_factor = 2.0 / sumfp(w, winlen);

    FP_TYPE* xfrm, * xfrm_d;
    FP_TYPE* spec_magn, * spec_phse_, * spec_phse, * spec_phse_d;
    xfrm = xfrm_d = spec_magn = spec_phse_ = spec_phse = spec_phse_d = NULL;
    
    xfrm = fetch_frame(x, nx, tn, winlen);
    FP_TYPE mean_xfrm = sumfp(xfrm, winlen) / winlen;
    if(dc != NULL) dc[n] = minfp(xfrm, winlen);
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
      spectrogram[n][i] = log_3(spec_magn[i] * norm_factor);
      if(isnan(spectrogram[n][i]) || isinf(spectrogram[n][i]))
        spectrogram[n][i] = -100.0;
      if(phasegram) phasegram[n][i] = spec_phse[i];
      if(phasegram_d) phasegram_d[n][i] = spec_phse_d[i];
    }

    free(spec_magn);
    if(spec_phse) free(spec_phse);
    if(spec_phse_d) free(spec_phse_d);
    free(xfrm);
    if(xfrm_d) free(xfrm_d);
    free(w);
  }
}

// a simple routine for refining F0 trajectory
// It is recommended to use a dedicated external routine for F0 refinement instead.
static FP_TYPE* refine_f0(llsm_parameters param, int nfft, int fs, FP_TYPE* f0, int nf0,
  FP_TYPE** spectrogram, FP_TYPE** phasegram, FP_TYPE** phasegram_d) {
  
  FP_TYPE* rf0 = calloc(nf0, sizeof(FP_TYPE));
  if(param.a_f0refine == 0) {
    for(int t = 0; t < nf0; t ++)
      rf0[t] = f0[t];
    return rf0;
  }

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
    
    FP_TYPE peak_freq;
    qifft(spectrogram[t], i_max, & peak_freq);
    peak_freq *= (FP_TYPE)fs / nfft;
    
    FP_TYPE p   = phasegram  [t][i_max];
    FP_TYPE p_d = phasegram_d[t][i_max];
    p   -= floor(p   / 2.0 / M_PI) * 2.0 * M_PI;
    p_d -= floor(p_d / 2.0 / M_PI) * 2.0 * M_PI;
    if(p < p_d)
      p += 2.0 * M_PI;
    
    rf0[t] = (p - p_d) / 2.0 / M_PI * fs;
    if(fabs(rf0[t] / fs * nfft - i_max) > 1)
      rf0[t] = peak_freq;
  }
  
  // take out spurious fluctuations
  for(int t = 1; t < nf0 - 1; t ++) {
    if(rf0[t] <= 0) continue;
    if(rf0[t + 1] > 0 && fabs(rf0[t + 1] - rf0[t]) > 10) {
      if(rf0[t - 1] == 0) rf0[t] = rf0[t + 1];
      else rf0[t] = (rf0[t - 1] + rf0[t + 1]) / 2.0;
    }
  }
  
  return rf0;
}

// estimate harmonic parameters by peak-picking
static int harmonic_peak_picking(llsm_parameters param, FP_TYPE* spectrum, FP_TYPE* phasetrum,
  int nfft, int fs, int nhar, FP_TYPE f0, FP_TYPE* freq, FP_TYPE* ampl, FP_TYPE* phse) {
  const FP_TYPE hardev = 0.3;

  if(f0 < 10) return 0;
  for(int i = 1; i <= nhar; i ++) {
    if(f0 * i > param.a_mvf || f0 * i > fs / 2)
      return i - 1;
    
    int l_idx = round(f0 * (i - hardev) / fs * nfft);
    int u_idx = round(f0 * (i + hardev) / fs * nfft);
    l_idx = l_idx < 1 ? 1 : l_idx;
    u_idx = u_idx > nfft / 2 ? nfft / 2 : u_idx;

    int peak_bin = cig_find_peak(spectrum, l_idx, u_idx, 1);
    
    FP_TYPE peak_freq, peak_ampl;
    peak_ampl = qifft(spectrum, peak_bin, & peak_freq);
    freq[i - 1] = peak_freq * fs / nfft;
    ampl[i - 1] = exp_3(peak_ampl);
    phse[i - 1] = linterp(phasetrum[(int)peak_freq], phasetrum[(int)peak_freq + 1], fmod(peak_freq, 1.0));
    
    if(isnan(ampl[i - 1])) // peak amplitude goes nan if one of the bins = -INF
      ampl[i - 1] = 0;
  }
  return nhar;
}

// generate a bunch of pitch-synchronous time marks over the entire signal
static int* gen_ps_instants(FP_TYPE* f0, int nhop, int nf0, int nx, int fs, int* ni) {
  FP_TYPE t = 0;
  *ni = 0;
  int ci = 100;
  int* instants = calloc(100, sizeof(int));
  while(t < (FP_TYPE)nx / fs) {
    int i = t * fs / nhop;
    instants[(*ni) ++] = t * fs;
    if((*ni) * 2 > ci) {
      ci *= 2;
      instants = realloc(instants, ci * sizeof(int));
    }
    t += 1.0 / max(50, f0[i >= nf0 ? nf0 - 1 : i]);
  }
  return instants;
}

// reshape colored noise according to some time-varying PSD
static FP_TYPE* filter_noise(llsm_parameters param, llsm_conf conf, FP_TYPE* x, int nx,
  FP_TYPE** warpped_spectrogram) {
  const int nfade = 16;
  int nhop = conf.nhop;
  int nwin = nhop * 4;
  int nfft = pow(2, ceil(log2(nwin + 2 * nfade)));
  int ns = nfft / 2 + 1;
  int nfrm = conf.nfrm;
  int frm0 = max(0, (param.s_n0 - nfft) / nhop);
  int frm1 = param.s_n1 == 0 ? nfrm : min(nfrm, (param.s_n1 + nfft) / nhop);

  FP_TYPE* y = calloc(nx, sizeof(FP_TYPE));
  FP_TYPE* w = hanning(nwin);
  FP_TYPE* realbuff = calloc(nfft, sizeof(FP_TYPE));
  FP_TYPE* imagbuff = calloc(nfft, sizeof(FP_TYPE));
  FP_TYPE* xfrm = calloc(nfft, sizeof(FP_TYPE));
  FP_TYPE* yfrm = calloc(nfft, sizeof(FP_TYPE));
  FP_TYPE* fftbuff = calloc(nfft * 2, sizeof(FP_TYPE));
  FP_TYPE norm_factor = 0;
  int i = 0;
  while(i < nwin) {
    norm_factor += w[i];
    i += nhop;
  }
  FP_TYPE weight_factor = 0.5 * sumfp(w, nwin);

  FP_TYPE* freqwarp = llsm_warp_freq(0, conf.nosf, conf.nnos, conf.noswarp);
  for(int i = frm0; i < frm1; i ++) {
    int t = i * nhop;
    FP_TYPE* xwin = fetch_frame(x, nx, t, nwin);
    for(int j = 0; j < nwin; j ++) xwin[j] *= w[j];
    for(int j = 0; j < nfft; j ++) xfrm[j] = 0;
    for(int j = 0; j < nwin; j ++) xfrm[j - nwin / 2 + nfft / 2] = xwin[j];
    fft(xfrm, NULL, realbuff, imagbuff, nfft, fftbuff);
    
    // add E_EPS to avoid subnormal floats
    FP_TYPE* xA = abscplx(realbuff, imagbuff, ns);
    for(int j = 0; j < ns; j ++) xA[j] = xA[j] / weight_factor;
    FP_TYPE* xE = llsm_geometric_envelope(xA, nfft, param.s_fs, freqwarp, conf.nnos);
    for(int j = 0; j < conf.nnos; j ++)
      xE[j] = warpped_spectrogram[i][j] - log(xE[j] + M_EPS);
    FP_TYPE* spec = llsm_spectrum_from_envelope(freqwarp, xE, conf.nnos, nfft, param.s_fs);

    for(int j = 0; j < ns; j ++) {
      FP_TYPE jgain = exp_2(spec[j]) + M_EPS;
      realbuff[j] *= jgain;
      imagbuff[j] *= jgain;
    }
    free(xE);
    free(xA);
    
    complete_symm (realbuff, nfft);
    complete_asymm(imagbuff, nfft);
    ifft(realbuff, imagbuff, yfrm, NULL, nfft, fftbuff);
    for(int j = 0; j < nfade; j ++) {
      yfrm[j] *= (FP_TYPE)j / nfade;
      yfrm[nfft - j - 1] *= 1.0 - (FP_TYPE)j / nfade;
    }
    for(int j = 0; j < nfft; j ++) {
      int idx = t + j - nfft / 2;
      if(idx >= 0 && idx < nx)
        y[idx] += yfrm[j] / norm_factor;
    }
    free(spec);
    free(xwin);
  }
  
  free(w);
  free(xfrm); free(yfrm);
  free(realbuff); free(imagbuff);
  free(fftbuff);
  free(freqwarp);
  return y;
}

// subtract lower envelope from a periodic signal
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

// lengthen a noise signal by duplication
static FP_TYPE* stretch_static_noise(FP_TYPE* x, int nx, int ny, int novl) {
  FP_TYPE* y = calloc(ny, sizeof(FP_TYPE));
  for(int i = 0; i < min(nx, ny); i ++)
    y[i] = x[i];
  if(ny <= nx) return y;
  
  int top = nx;
  while(1) {
    for(int i = 0; i < novl; i ++) {
      FP_TYPE r = (FP_TYPE)i / novl;
      y[top - novl + i] *= 1.0 - r;
      y[top - novl + i] += x[i] * r;
      y[top - novl + i] /= sqrt(2 * r * (r - 1) + 1);
    }
    for(int i = 0; i < nx - novl; i ++) {
      if(top + i >= ny) return y;
      y[top + i] = x[i + novl];
    }
    top += nx - novl;
  }
  
  return y;
}

llsm_frame* llsm_create_frame(int nhar, int nhare, int nnos, int nchannel) {
  llsm_frame* ret = malloc(sizeof(llsm_frame));

  ret -> sinu = malloc(sizeof(llsm_sinframe));
  ret -> sinu -> ampl = calloc(nhar, sizeof(FP_TYPE));
  ret -> sinu -> phse = calloc(nhar, sizeof(FP_TYPE));
  ret -> sinu -> nhar = nhar;
  
  ret -> noise = malloc(sizeof(llsm_nosframe));
  ret -> noise -> eenv = calloc(nchannel, sizeof(llsm_sinframe*));
  ret -> noise -> emin = calloc(nchannel, sizeof(FP_TYPE));
  ret -> noise -> spec = calloc(nnos, sizeof(FP_TYPE));
  ret -> noise -> nchannel = nchannel;
  ret -> noise -> nnos = nnos;
  
  for(int i = 0; i < nchannel; i ++) {
    ret -> noise -> eenv[i] = malloc(sizeof(llsm_sinframe));
    llsm_sinframe* ieenv = ret -> noise -> eenv[i];
    ieenv -> ampl = calloc(nhare, sizeof(FP_TYPE));
    ieenv -> phse = calloc(nhare, sizeof(FP_TYPE));
    ieenv -> nhar = nhare;
  }
  
  return ret;
}

llsm_layer0* llsm_create_empty_layer0(llsm_conf conf) {
  llsm_layer0* ret = calloc(1, sizeof(llsm_layer0));
  ret -> conf = conf;
  ret -> conf.nosbandf = calloc(conf.nnosband - 1, sizeof(FP_TYPE));
  for(int i = 0; i < conf.nnosband - 1; i ++)
    ret -> conf.nosbandf[i] = conf.nosbandf[i];
  if(conf.nfrm > 0)
    ret -> frames = calloc(conf.nfrm, sizeof(llsm_frame*));
  else
    ret -> frames = NULL;
  return ret;
}

void llsm_copy_frame(llsm_frame* dst, llsm_frame* src) {
  llsm_copy_sinframe(dst -> sinu, src -> sinu);
  llsm_copy_nosframe(dst -> noise, src -> noise);
  dst -> f0 = src -> f0;
}

void llsm_copy_sinframe(llsm_sinframe* dst, llsm_sinframe* src) {
  if(dst -> nhar < src -> nhar) {
    dst -> ampl = realloc(dst -> ampl, src -> nhar * sizeof(FP_TYPE));
    dst -> phse = realloc(dst -> phse, src -> nhar * sizeof(FP_TYPE));
  }
  dst -> nhar = src -> nhar;
  for(int i = 0; i < dst -> nhar; i ++) {
    dst -> ampl[i] = src -> ampl[i];
    dst -> phse[i] = src -> phse[i];
  }
}

void llsm_copy_nosframe(llsm_nosframe* dst, llsm_nosframe* src) {
  if(dst -> nnos < src -> nnos)
    dst -> spec = realloc(dst -> spec, src -> nnos * sizeof(FP_TYPE));
  dst -> nnos = src -> nnos;
  for(int i = 0; i < dst -> nnos; i ++)
    dst -> spec[i] = src -> spec[i];
  if(dst -> nchannel < src -> nchannel) {
    dst -> emin = realloc(dst -> emin, src -> nchannel * sizeof(FP_TYPE));
    dst -> eenv = realloc(dst -> eenv, src -> nchannel * sizeof(llsm_sinframe*));
    for(int i = dst -> nchannel; i < src -> nchannel; i ++) {
      dst -> eenv[i] = malloc(sizeof(llsm_sinframe));
      dst -> eenv[i] -> nhar = 0;
      dst -> eenv[i] -> ampl = NULL;
      dst -> eenv[i] -> phse = NULL;
    }
  } else if(dst -> nchannel > src -> nchannel)
    for(int i = src -> nchannel; i < dst -> nchannel; i ++) {
      free(dst -> eenv[i] -> ampl);
      free(dst -> eenv[i] -> phse);
      free(dst -> eenv[i]);
    }
  dst -> nchannel = src -> nchannel;
  for(int i = 0; i < dst -> nchannel; i ++) {
    dst -> emin[i] = src -> emin[i];
    llsm_copy_sinframe(dst -> eenv[i], src -> eenv[i]);
  }
}

void llsm_delete_frame(llsm_frame* dst) {
  if(dst == NULL) return;
  if(dst -> sinu != NULL) {
    free(dst -> sinu -> ampl);
    free(dst -> sinu -> phse);
    free(dst -> sinu);
  }
  if(dst -> noise != NULL) {
    for(int i = 0; i < dst -> noise -> nchannel; i ++) {
      llsm_sinframe* ieenv = dst -> noise -> eenv[i];
      if(ieenv == NULL) continue;
      free(ieenv -> ampl);
      free(ieenv -> phse);
      free(ieenv);
    }
    free(dst -> noise -> eenv);
    free(dst -> noise -> emin);
    free(dst -> noise -> spec);
    free(dst -> noise);
  }
  free(dst);
}

void llsm_delete_layer0(llsm_layer0* dst) {
  if(dst == NULL) return;
  int nfrm = dst -> conf.nfrm;
  if(dst -> conf.nosbandf != NULL) free(dst -> conf.nosbandf);
  if(dst -> frames != NULL) {
    for(int i = 0; i < nfrm; i ++)
      llsm_delete_frame(dst -> frames[i]);
    free(dst -> frames);
  }
  free(dst);
}

void llsm_delete_output(llsm_output* dst) {
  if(dst == NULL) return;
  free(dst -> y);
  free(dst -> y_sin);
  free(dst -> y_nos);
  free(dst);
}

static void llsm_layer0_analyze_noise(llsm_parameters param, llsm_layer0* model, FP_TYPE* x,
    FP_TYPE* residual, int nx, int fs, int* analysis_instant, FP_TYPE* rf0, int nf0) {
  int nfft = pow(2, ceil(log2(fs * param.a_tfft)));
  int ns = nfft / 2 + 1;
  FP_TYPE nosf = param.a_nosf > 0.0 ? param.a_nosf : fs / 2.0;
  int* winsize = calloc(nf0, sizeof(int));
  
  FP_TYPE* tmpfreq = calloc(param.a_nhar, sizeof(FP_TYPE));
  FP_TYPE* tmpampl = calloc(param.a_nhar, sizeof(FP_TYPE));
  FP_TYPE* tmpphse = calloc(param.a_nhar, sizeof(FP_TYPE));

  for(int i = 0; i < nf0; i ++) winsize[i] = param.a_nhop * 4;
  
  // estimate time-varying PSD

  FP_TYPE** noise_spectrogram = (FP_TYPE**)malloc2d(nf0, ns, sizeof(FP_TYPE));
  spectrogram_analyze(residual, nx, analysis_instant, nf0, winsize, nfft, fs, "hanning",
    noise_spectrogram, NULL, NULL, NULL);
  
  FP_TYPE* freqwarp = llsm_warp_freq(0, nosf, param.a_nnos, param.a_noswarp);
  for(int t = 0; t < nf0; t ++) {
    for(int j = 0; j < ns; j ++)
      noise_spectrogram[t][j] = exp(noise_spectrogram[t][j]);
    FP_TYPE* wrpspec = llsm_geometric_envelope(noise_spectrogram[t], nfft, fs, freqwarp, param.a_nnos);
    for(int k = 0; k < param.a_nnos; k ++)
      model -> frames[t] -> noise -> spec[k] = log(wrpspec[k]);
    model -> frames[t] -> noise -> nnos = param.a_nnos;
    free(wrpspec);
  }
  free(freqwarp);
  free2d(noise_spectrogram, nf0);

  // extract subband envelopes

  for(int i = 0; i < nf0; i ++) {
    if(rf0[i] == 0) winsize[i] = model -> conf.nhop * 2;
    else winsize[i] = round(fs / rf0[i] * 1.8) * 2;
  }
  // for each noise channel
  for(int b = 0; b < param.a_nnosband; b ++) {
    int filtord = LLSM_CHEBY_ORDER * 2;
    FP_TYPE* b_filtered = NULL;
    FP_TYPE favg = 0;
    if(b == 0) {
      b_filtered = chebyfilt(residual, nx, param.a_nosbandf[b] / fs * 2.0, 0, "lowpass");
      favg = param.a_nosbandf[b] / 2.0;
      filtord -= LLSM_CHEBY_ORDER;
    } else
    if(b == param.a_nnosband - 1) {
      b_filtered = chebyfilt((param.a_nosbandf[b - 1] > 6000.0 ? x : residual), nx,
        param.a_nosbandf[b - 1] / fs * 2.0, 0, "highpass");
      favg = (param.a_nosbandf[b - 1] + fs) / 2.0;
      filtord -= LLSM_CHEBY_ORDER;
    } else {
      b_filtered = chebyfilt((param.a_nosbandf[b - 1] > 6000.0 ? x : residual), nx,
        param.a_nosbandf[b - 1] / fs * 2.0, param.a_nosbandf[b] / fs * 2.0, "bandpass");
      favg = (param.a_nosbandf[b - 1] + param.a_nosbandf[b]) / 2.0;
    }

    FP_TYPE* b_env = moving_rms(b_filtered, nx, fs / favg * 5);
    free(b_filtered);
    
#   ifdef DEBUG_EXPORT_WAV
    char* strdup(const char*);
    char* writepath = strdup("./env .wav");
    writepath[5] = '0' + b;
    gain_write(b_env, nx, fs, writepath, 10);
    free(writepath);
#   endif

    FP_TYPE** b_spectrogram = (FP_TYPE**)malloc2d(nf0, nfft / 2 + 1, sizeof(FP_TYPE));
    FP_TYPE** b_phasegram   = (FP_TYPE**)malloc2d(nf0, nfft / 2 + 1, sizeof(FP_TYPE));
    FP_TYPE* emin = calloc(nf0, sizeof(FP_TYPE));
    spectrogram_analyze(b_env, nx, analysis_instant, nf0, winsize, nfft, fs, "blackman",
      b_spectrogram, b_phasegram, NULL, emin);
    for(int i = 0; i < nf0; i ++)
      model -> frames[i] -> noise -> emin[b] = emin[i];
    free(b_env);
    free(emin);
    
    for(int i = 0; i < nf0; i ++) {
      if(rf0[i] <= 0.0) continue;
      llsm_sinframe* ieenv = model -> frames[i] -> noise -> eenv[b];
      int inhare = harmonic_peak_picking(param, b_spectrogram[i], b_phasegram[i],
        nfft, fs, param.a_nhare, rf0[i], tmpfreq, tmpampl, tmpphse);
      for(int k = 0; k < inhare; k ++) {
        ieenv -> ampl[k] = tmpampl[k];
        ieenv -> phse[k] = tmpphse[k];
      }
    }
    free2d(b_spectrogram, nf0);
    free2d(b_phasegram, nf0);
  }
  
  free(winsize);
  free(tmpfreq);
  free(tmpampl);
  free(tmpphse);
}

llsm_layer0* llsm_layer0_analyze(llsm_parameters param, FP_TYPE* x, int nx, int fs,
    FP_TYPE* f0, int nf0, FP_TYPE** xap) {
  llsm_layer0* model = malloc(sizeof(llsm_layer0));
  int nfft = pow(2, ceil(log2(fs * param.a_tfft)));
  int ns = nfft / 2 + 1;
  FP_TYPE nosf = param.a_nosf > 0.0 ? param.a_nosf : fs / 2.0;

  model -> conf.nfrm = nf0;
  model -> conf.nhop = param.a_nhop;
  model -> conf.thop = (FP_TYPE)param.a_nhop / fs;
  model -> conf.nhar = param.a_nhar;
  model -> conf.nhare = param.a_nhare;
  model -> conf.nnos = param.a_nnos;
  model -> conf.nosf = nosf;
  model -> conf.mvf = param.a_mvf;
  model -> conf.noswarp = param.a_noswarp;
  model -> conf.nosbandf = calloc(param.a_nnosband - 1, sizeof(FP_TYPE));
  model -> conf.nnosband = param.a_nnosband;
  for(int i = 0; i < param.a_nnosband - 1; i ++)
    model -> conf.nosbandf[i] = param.a_nosbandf[i];

  // compute spectrogram

  FP_TYPE** spectrogram = (FP_TYPE**)malloc2d(nf0, ns, sizeof(FP_TYPE));
  FP_TYPE** phasegram   = (FP_TYPE**)malloc2d(nf0, ns, sizeof(FP_TYPE));
  FP_TYPE** phasegram_d = (FP_TYPE**)malloc2d(nf0, ns, sizeof(FP_TYPE));

  int* winsize = calloc(nf0, sizeof(int));
  int* analysis_instant = calloc(nf0, sizeof(int));
  for(int i = 0; i < nf0; i ++) {
    if(f0[i] == 0) winsize[i] = model -> conf.nhop * 2;
    else winsize[i] = round(fs / f0[i] * param.a_wsize) * 2;
    analysis_instant[i] = model -> conf.nhop * i;
  }
  spectrogram_analyze(x, nx, analysis_instant, nf0, winsize, nfft, fs, "blackman", spectrogram, phasegram, phasegram_d, NULL);

  FP_TYPE* rf0 = refine_f0(param, nfft, fs, f0, nf0, spectrogram, phasegram, phasegram_d);

  // harmonic analysis

  model -> frames = calloc(nf0, sizeof(llsm_frame*));
  FP_TYPE* tmpfreq = calloc(param.a_nhar, sizeof(FP_TYPE));
  FP_TYPE* tmpampl = calloc(param.a_nhar, sizeof(FP_TYPE));
  FP_TYPE* tmpphse = calloc(param.a_nhar, sizeof(FP_TYPE));
  FP_TYPE** matfreq = calloc(nf0, sizeof(FP_TYPE*));
  
  for(int i = 0; i < nf0; i ++) {
    if(rf0[i] <= 0.0) {
      model -> frames[i] = llsm_create_frame(0, param.a_nhare, param.a_nnos, param.a_nnosband);
      model -> frames[i] -> f0 = 0.0;
      continue;
    }
    
    int inhar = harmonic_peak_picking(param, spectrogram[i], phasegram[i],
      nfft, fs, param.a_nhar, rf0[i], tmpfreq, tmpampl, tmpphse);
    if(inhar == 0) {
      model -> frames[i] = llsm_create_frame(0, param.a_nhare, param.a_nnos, param.a_nnosband);
      model -> frames[i] -> f0 = 0.0;
      continue;
    }

    matfreq[i] = malloc(inhar * sizeof(FP_TYPE));
    model -> frames[i] = llsm_create_frame(inhar, param.a_nhare, param.a_nnos, param.a_nnosband);
    model -> frames[i] -> f0 = rf0[i];

    for(int k = 0; k < inhar; k ++) {
      model -> frames[i] -> sinu -> ampl[k] = tmpampl[k];
      model -> frames[i] -> sinu -> phse[k] = tmpphse[k];
      matfreq[i][k] = tmpfreq[k];
    }
  }

  free2d(spectrogram, nf0);
  free2d(phasegram, nf0);
  free2d(phasegram_d, nf0);
  free(tmpfreq);
  free(tmpampl);
  free(tmpphse);

  // resynthesis of harmonic component

  FP_TYPE* resynth = calloc(nx + param.a_nhop, sizeof(FP_TYPE));
  int nrwin = 2 * param.a_nhop;
  FP_TYPE* rwin = hanning(nrwin);
  for(int i = 0; i < nf0; i ++) {
    int tn = i * param.a_nhop;
    FP_TYPE* rfrm = gensins(
      matfreq[i], model -> frames[i] -> sinu -> ampl, model -> frames[i] -> sinu -> phse,
      model -> frames[i] -> sinu -> nhar, fs, nrwin);
    for(int j = 0; j < nrwin; j ++)
      if(tn + j - nrwin / 2 > 0)
        resynth[tn + j - nrwin / 2] += rfrm[j] * rwin[j];
    free(rfrm);
  }
  free(rwin);
  free2d(matfreq, nf0);
# ifdef DEBUG_EXPORT_WAV
  wavwrite(resynth, nx, fs, 16, "sin.wav");
# endif
  
  // extract noise/aperiodic component by subtraction

  for(int i = 0; i < nx; i ++)
    resynth[i] = x[i] - resynth[i];
  if(xap != NULL) {
    *xap = calloc(nx, sizeof(FP_TYPE));
    for(int i = 0; i < nx; i ++)
      xap[0][i] = resynth[i];
  }

# ifdef DEBUG_EXPORT_WAV
  wavwrite(resynth, nx, fs, 16, "res.wav");
# endif

  // analyze noise component
  llsm_layer0_analyze_noise(param, model, x, resynth, nx, fs, analysis_instant, rf0, nf0);

  free(resynth);
  free(rf0);
  free(winsize);
  free(analysis_instant);

  return model;
}

void llsm_layer0_phaseshift(llsm_layer0* dst, FP_TYPE* phaseshift) {
  int nfrm = dst -> conf.nfrm;
  for(int i = 0; i < nfrm; i ++) {
    llsm_sinframe* isinu = dst -> frames[i] -> sinu;
    if(isinu -> nhar <= 0) continue;
    
    for(int j = 0; j < isinu -> nhar; j ++)
      isinu -> phse[j] = wrap(isinu -> phse[j] + phaseshift[i] * (1.0 + j));
    for(int b = 0; b < dst -> conf.nnosband; b ++) {
      llsm_sinframe* ieenv = dst -> frames[i] -> noise -> eenv[b];
      for(int j = 0; j < ieenv -> nhar; j ++)
        ieenv -> phse[j] = wrap(ieenv -> phse[j] + phaseshift[i] * (1.0 + j));
    }
  }
}

llsm_output* llsm_layer0_synthesize(llsm_parameters param, llsm_layer0* model) {
  int nfrm = model -> conf.nfrm;
  int nhop = model -> conf.nhop;
  int fs = param.s_fs;
  int nfft = nhop * 4;
  llsm_output* ret = malloc(sizeof(llsm_output));
  ret -> fs = fs;
  int* ny = & ret -> ny;
  int frm0 = max(0, (param.s_n0 - nfft) / nhop);
  int frm1 = param.s_n1 == 0 ? nfrm : min(nfrm, (param.s_n1 + nfft) / nhop);

  *ny = nfrm * nhop + nfft;
  FP_TYPE* y_sin = calloc(*ny, sizeof(FP_TYPE));
  FP_TYPE* f0 = calloc(nfrm, sizeof(FP_TYPE));
  for(int i = 0; i < nfrm; i ++)
    f0[i] = model -> frames[i] -> f0;
  
  // synthesize harmonic component
  FP_TYPE* ola_window = hanning(nhop * 2);
  if(! param.s_noiseonly) {
    for(int i = frm0; i < min(nfrm, frm1); i ++) {
      int tn = i * nhop;
      llsm_sinframe* sinu = model -> frames[i] -> sinu;
      FP_TYPE* freq = malloc(sinu -> nhar * sizeof(FP_TYPE));
      for(int k = 0; k < sinu -> nhar; k ++) freq[k] = (1.0 + k) * f0[i];
      FP_TYPE* sfrm = gensins(
        freq, model -> frames[i] -> sinu -> ampl, model -> frames[i] -> sinu -> phse,
        model -> frames[i] -> sinu -> nhar, fs, nhop * 2);
      for(int j = 0; j < nhop * 2; j ++)
        if(tn + j - nhop > 0)
          y_sin[tn + j - nhop] += sfrm[j] * ola_window[j];
      free(sfrm);
      free(freq);
    }
  }
  
  FP_TYPE* s = white_noise(1.0, fs); // one-second-long noise template
  FP_TYPE* noise_excitation = calloc(*ny, sizeof(FP_TYPE));

  // for each noise channel
  for(int b = 0; b < model -> conf.nnosband; b ++) {
    FP_TYPE* b_env = calloc(*ny, sizeof(FP_TYPE));
    FP_TYPE* b_env_mix = calloc(*ny, sizeof(FP_TYPE));
    
    // generate noise envelope
    for(int i = frm0; i < frm1; i ++) {
      if(f0[i] <= 0.0) continue;
      llsm_sinframe* ieenv = model -> frames[i] -> noise -> eenv[b];
      int tn = i * nhop;
      FP_TYPE ifreq[32];
      for(int k = 0; k < ieenv -> nhar; k ++) ifreq[k] = f0[i] * (k + 1.0);
      FP_TYPE* b_env_frame = gensins(ifreq, ieenv -> ampl, ieenv -> phse, ieenv -> nhar, fs, nhop * 2);
      for(int j = 0; j < nhop * 2; j ++)
        if(tn + j - nhop > 0)
          b_env[tn + j - nhop] += b_env_frame[j] * ola_window[j];
      free(b_env_frame);
    }
    
    subtract_minimum_envelope(b_env, *ny, f0, nhop, nfrm, fs);
    
    // generate band-limited noise
    int filtord = LLSM_CHEBY_ORDER * 2;
    FP_TYPE* b_template = NULL;
    FP_TYPE* b_filtered = NULL;
    if(b == 0) {
      b_template = chebyfilt(s, fs, model -> conf.nosbandf[b] / fs * 2.0, 0, "lowpass");
      filtord -= LLSM_CHEBY_ORDER;
    } else
    if(b == model -> conf.nnosband - 1) {
      b_template = chebyfilt(s, fs, model -> conf.nosbandf[b - 1] / fs * 2.0,
        model -> conf.nosf / fs * 2.0, "bandpass");
    } else {
      b_template = chebyfilt(s, fs, model -> conf.nosbandf[b - 1] / fs * 2.0,
        model -> conf.nosbandf[b] / fs * 2.0, "bandpass");
    }
    b_filtered = stretch_static_noise(b_template, fs, *ny, 256);
    free(b_template);
    
    // normalize filtered noise and recover envelope at DC
    FP_TYPE* b_normalized = calloc(*ny, sizeof(FP_TYPE));
    for(int i = frm0; i < frm1; i ++) {
      FP_TYPE* hfrm = fetch_frame(b_filtered, *ny, i * nhop, nhop * 2);
      FP_TYPE* efrm = fetch_frame(b_env, *ny, i * nhop, nhop * 2);
      FP_TYPE havg = 0;
      for(int j = 0; j < nhop * 2; j ++)
        havg += hfrm[j] * hfrm[j];
      havg /= nhop * 2;
      for(int j = 0; j < nhop * 2; j ++)
        hfrm[j] *= sqrt(1.0 / (havg + 1e-8)); // in case havg is very small
      
      if(model -> frames[i] -> f0 == 0)
        for(int j = 0; j < nhop * 2; j ++)
          efrm[j] += model -> frames[i] -> noise -> emin[b];
      for(int j = 0; j < nhop * 2; j ++) {
        efrm[j] = max(1e-7, efrm[j]);
      }
      
      for(int j = 0; j < nhop * 2; j ++)
        if(i * nhop + j - nhop > 0) {
          b_normalized[i * nhop + j - nhop] += hfrm[j] * ola_window[j];
          b_env_mix   [i * nhop + j - nhop] += efrm[j] * ola_window[j];
        }
      
      free(hfrm);
      free(efrm);
    }

    // modulate noise by envelope
    for(int i = 0; i < *ny; i ++)
      noise_excitation[i] += b_normalized[i] * b_env_mix[i];
    
    free(b_filtered);
    free(b_normalized);
    free(b_env_mix);
    free(b_env);
  }
  free(f0);
  free(s);

  // re-color excitation
  FP_TYPE** noisegram = calloc(nfrm, sizeof(FP_TYPE*));
  for(int i = 0; i < nfrm; i ++)
    noisegram[i] = model -> frames[i] -> noise -> spec;
  FP_TYPE* y_nos = filter_noise(param, model -> conf, noise_excitation, *ny, noisegram);
  free(noisegram);
# ifdef DEBUG_EXPORT_WAV
  wavwrite(noise_excitation, *ny, fs, 16, "synexc.wav");
  wavwrite(y_nos, *ny, fs, 16, "synres.wav");
  wavwrite(y_sin, *ny, fs, 16, "synsin.wav");
# endif
  free(noise_excitation);

  FP_TYPE* y_mix = calloc(*ny, sizeof(FP_TYPE));
  for(int i = 0; i < *ny; i ++)
    y_mix[i] = y_sin[i] + y_nos[i];

  ret -> y = y_mix;
  ret -> y_sin = y_sin;
  ret -> y_nos = y_nos;
  
  free(ola_window);
  return ret;
}
