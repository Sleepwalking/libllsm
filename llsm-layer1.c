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

#include "llsm.h"
#include <stdlib.h>
#include <string.h>
#include "envelope.h"
#include "math-funcs.h"

FP_TYPE* llsm_uniform_faxis(int nfft, int fs) {
  FP_TYPE* faxis = calloc(nfft / 2 + 1, sizeof(FP_TYPE));
  for(int i = 0; i < nfft / 2 + 1; i ++)
    faxis[i] = (FP_TYPE)i * fs / nfft;
  return faxis;
}

FP_TYPE* llsm_liprad(FP_TYPE* freq, int nf, FP_TYPE radius, FP_TYPE* dst_phaseresp) {
  FP_TYPE Rr = 128.0 / 9.0 / M_PI / M_PI;
  FP_TYPE Lr = 8.0 * radius / 100.0 / 3.0 / M_PI / 340.0;
  FP_TYPE* dst_resp_magn = calloc(nf, sizeof(FP_TYPE));
  
  for(int i = 0; i < nf; i ++) {
    cplx iresp = c_mul(c_cplx(0, 1),
      c_div(c_cplx(freq[i] * 2.0 * M_PI * Lr * Rr, 0),
            c_cplx(Rr, freq[i] * 2.0 * M_PI * Lr)
      ));
    dst_resp_magn[i] = c_abs(iresp);
    if(dst_phaseresp) dst_phaseresp[i] = c_arg(iresp);
  }
  
  return dst_resp_magn;
}

typedef struct {
  FP_TYPE** G;
  FP_TYPE* rd;
  int nhar;
  int nG;
} lfrdset;

// generate a group of LF model freq. responses for different Rd parameters
static lfrdset* generate_lfrd_responses(FP_TYPE* rd, int nrd, int nhar) {
  lfrdset* ret = malloc(sizeof(lfrdset));
  ret -> nG = nrd;
  ret -> nhar = nhar;
  ret -> G = calloc(nrd, sizeof(FP_TYPE*));
  ret -> rd = calloc(nrd, sizeof(FP_TYPE));
  FP_TYPE f0 = 200.0;
  FP_TYPE* freq = calloc(nhar, sizeof(FP_TYPE));
  for(int i = 0; i < nhar; i ++) freq[i] = f0 * (1.0 + i);
  for(int i = 0; i < nrd; i ++) {
    ret -> rd[i] = rd[i];
    lfmodel ilf = lfmodel_from_rd(rd[i], 1.0 / f0, 1.0);
    ret -> G[i] = lfmodel_spectrum(ilf, freq, nhar, NULL);
    for(int j = 0; j < nhar; j ++) {
      ret -> G[i][j] /= j + 1.0;
      ret -> G[i][j] *= ret -> G[i][j];
    }
  }
  free(freq);
  return ret;
}

static int find_extremum(FP_TYPE* arr, int n, int direction) {
  FP_TYPE extremum = arr[0] * direction;
  FP_TYPE extremum_idx = 0;
  for(int i = 0; i < n; i ++)
    if(arr[i] * direction > extremum) {
      extremum = arr[i] * direction;
      extremum_idx = i;
    }
  return extremum_idx;
}

// amplitude-only Rd parameter estimation (single-frame)
static FP_TYPE amprd(FP_TYPE* ampl, int nhar, FP_TYPE f0, lfrdset* lfresp,
  FP_TYPE* faxis, FP_TYPE* liprad, int nf) {
  int nhar_8k = min(lfresp -> nhar, min(8000.0 / f0, nhar));
  FP_TYPE* freq = calloc(nhar_8k, sizeof(FP_TYPE));
  FP_TYPE* ampl_sqr = calloc(nhar_8k, sizeof(FP_TYPE));
  for(int i = 0; i < nhar_8k; i ++)
    freq[i] = (i + 1.0) * f0;
  
  FP_TYPE faxis_upper = faxis[nf - 1] * 2 - faxis[nf - 2];
  FP_TYPE* ampl_lip = interp1u(faxis[0], faxis_upper, liprad, nf, freq, nhar_8k);
  for(int i = 0; i < nhar_8k; i ++)
    ampl_sqr[i] = max(1e-8, ampl[i] * ampl[i] / ampl_lip[i] / ampl_lip[i]);
  free(ampl_lip);
  
  FP_TYPE* dist_curve = calloc(lfresp -> nG, sizeof(FP_TYPE));
  FP_TYPE* G_normalized = calloc(nhar_8k, sizeof(FP_TYPE));
  for(int i = 0; i < lfresp -> nG; i ++) {
    FP_TYPE* G = lfresp -> G[i];
    FP_TYPE normalization = sumfp(ampl_sqr, nhar_8k) / sumfp(G, nhar_8k);
    for(int j = 0; j < nhar_8k; j ++) G_normalized[j] = G[j] * normalization;
    dist_curve[i] = itakura_saito(ampl_sqr, G_normalized, nhar_8k);
  }
  free(G_normalized);
  
  int valley = find_extremum(dist_curve, lfresp -> nG, -1);
  FP_TYPE rd_refined = lfresp -> rd[valley];
  if(valley > 0 && valley < lfresp -> nG - 1) {
    qifft(dist_curve, valley, & rd_refined);
    rd_refined = linterp(lfresp -> rd[(int)rd_refined],
      lfresp -> rd[(int)rd_refined + 1], fmod(rd_refined, 1.0));
  }
  
  free(dist_curve);
  free(freq);
  free(ampl_sqr);
  return rd_refined;
}

// amplitude-only Rd parameter estimation with smoothing
static FP_TYPE* llsm_analyze_rd_amponly(llsm_layer0* model, FP_TYPE* faxis, FP_TYPE* liprad, int nf) {
  int nfrm = model -> conf.nfrm;
  FP_TYPE* rd = calloc(nfrm, sizeof(FP_TYPE));
  const FP_TYPE rd_min = 0.02;
  const FP_TYPE rd_step = 0.04;
  const FP_TYPE rd_max = 3.0;
  const int nrd = (rd_max - rd_min) / rd_step + 1;
  FP_TYPE* range = calloc(nrd, sizeof(FP_TYPE));
  for(int i = 0; i < nrd; i ++) range[i] = rd_min + rd_step * i;
  lfrdset* lfresp = generate_lfrd_responses(range, nrd, 80);
  
# pragma omp parallel for
  for(int i = 0; i < nfrm; i ++) {
    FP_TYPE if0 = model -> frames[i] -> f0;
    if(if0 <= 0) continue;
    rd[i] = amprd(model -> frames[i] -> sinu -> ampl, model -> frames[i] -> sinu -> nhar,
      if0, lfresp, faxis, liprad, nf);
  }
  free(range);
  free(lfresp -> rd); free2d(lfresp -> G, lfresp -> nG);
  free(lfresp);
  
  FP_TYPE* rd_xi = calloc(nfrm + 2, sizeof(FP_TYPE));
  FP_TYPE* rd_yi = calloc(nfrm + 2, sizeof(FP_TYPE));
  FP_TYPE* rd_idx = calloc(nfrm + 1, sizeof(FP_TYPE));
  int rd_ni = 1;
  rd_xi[0] = 0;
  for(int i = 0; i < nfrm; i ++) {
    rd_idx[i] = i;
    if(rd[i] <= 0) continue;
    rd_xi[rd_ni] = i;
    rd_yi[rd_ni] = rd[i];
    rd_ni ++;
  }
  rd_yi[0] = rd_yi[1];
  rd_xi[rd_ni] = nfrm - 1;
  rd_yi[rd_ni] = rd_yi[rd_ni - 1];
  rd_ni ++;
  FP_TYPE* rd_filledup = interp1(rd_xi, rd_yi, rd_ni, rd_idx, nfrm);
  free(rd_idx);
  free(rd_xi);
  free(rd_yi);
  
  const int filtord = 5;
  FP_TYPE* avgfilt = boxcar(filtord);
  for(int i = 0; i < filtord; i ++) avgfilt[i] /= filtord;
  FP_TYPE* rd_filt = filtfilt(avgfilt, filtord, NULL, 1, rd_filledup, nfrm);
  for(int i = 0; i < filtord; i ++) {
    rd_filt[i] = rd_filt[filtord];
    rd_filt[nfrm - i - 1] = rd_filt[nfrm - filtord - 1];
  }
  for(int i = 0; i < nfrm; i ++)
    if(rd[i] > 0) rd[i] = rd_filt[i];
  free(avgfilt);
  free(rd_filledup);
  free(rd_filt);
  return rd;
}

// 3-period-long hanning windowed narrow-band spectrum
static FP_TYPE* calculate_harmonic_spectrum(FP_TYPE* ampl, int nhar, FP_TYPE f0, int nfft, int fs) {
  FP_TYPE* X = calloc(nfft / 2 + 1, sizeof(FP_TYPE));
  int nX = nfft / 2 + 1;
  int T = fs / f0 * 3.0;
  int width = ceil(f0 / fs * nfft * 1.5);
  for(int i = 0; i < nhar; i ++) {
    FP_TYPE ifreq = f0 * (1.0 + i);
    int center = round(ifreq / fs * nfft);
    for(int j = max(0, center - width); j < min(nX, center + width + 1); j ++) {
      FP_TYPE omega = ((FP_TYPE)j / nfft - ifreq / fs) * 2.0 * M_PI;
      FP_TYPE resp = 0.5 * safe_aliased_sinc(T, omega) +
        0.25 * safe_aliased_sinc(T, omega - 2.0 * M_PI / T) +
        0.25 * safe_aliased_sinc(T, omega + 2.0 * M_PI / T);
      X[j] = fmax(X[j], resp * ampl[i]);
      //printf("%f\n", X[j]);
    }
  }
  
  // normalization
  for(int i = 0; i < nfft / 2 + 1; i ++)
    X[i] *= f0 / fs;
  return X;
}

static FP_TYPE range_compress(FP_TYPE x) {
  if(x > -10) return x;
  return (x + 10.0) / 2 - 10.0;
}

static FP_TYPE range_decompress(FP_TYPE x) {
  if(x > -10) return x;
  return (x + 10.0) * 2 - 10.0;
}

FP_TYPE* llsm_harmonic_cheaptrick(FP_TYPE* ampl, int nhar, FP_TYPE f0, int nfft, FP_TYPE fs) {
  FP_TYPE* compampl = calloc(nhar, sizeof(FP_TYPE));
  FP_TYPE peak = log(maxfp(ampl, nhar)) + 1;
  for(int i = 0; i < nhar; i ++)
    compampl[i] = exp(range_compress(log(ampl[i]) - peak));

  FP_TYPE* X = calculate_harmonic_spectrum(compampl, nhar, f0, nfft, fs);
  FP_TYPE* full_spectrum = cig_spec2env(X, nfft, f0 / fs, nhar, NULL);
  free(X);
  for(int i = 0; i < nfft / 2 + 1; i ++)
    full_spectrum[i] = range_decompress(full_spectrum[i]) + peak;
  free(compampl);
  return full_spectrum;
}

FP_TYPE* llsm_harmonic_minphase(FP_TYPE* ampl, int nhar) {
  int nfft = max(64, pow(2, ceil(log2(nhar) + 2)));
  FP_TYPE* hari = calloc(nhar + 1, sizeof(FP_TYPE));
  FP_TYPE* hara = calloc(nhar + 1, sizeof(FP_TYPE));
  FP_TYPE* ffti = calloc(nfft / 2 + 1, sizeof(FP_TYPE));
  for(int i = 0; i < nhar; i ++) {
    hari[i + 1] = (i + 1.0) / (nhar + 1.0) * nfft / 2.0;
    hara[i + 1] = log(ampl[i] + 1e-10);
  }
  hara[0] = hara[1];
  for(int i = 0; i < nfft / 2 + 1; i ++)
    ffti[i] = i;
  
  FP_TYPE* hA = interp1u(0, hari[nhar] * 2 - hari[nhar - 1], hara, nhar + 1, ffti, nfft / 2 + 1);
  FP_TYPE* hP = minphase(hA, nfft);
  FP_TYPE* minP = interp1u(0, nfft / 2 + 1, hP, nfft / 2 + 1, hari, nhar + 1);
  for(int i = 1; i < nhar; i ++)
    minP[i - 1] = minP[i];
  free(hari);
  free(hara);
  free(ffti);
  free(hA);
  free(hP);
  return minP;
}

llsm_layer1* llsm_layer1_from_layer0(llsm_parameters param, llsm_layer0* mlv0, int nfft, int fs) {
  if(fs <= 0)
    fs = max(mlv0 -> conf.nosf, mlv0 -> conf.mvf);
  int nfrm = mlv0 -> conf.nfrm;
  int ns = nfft / 2 + 1;
  FP_TYPE* faxis = llsm_uniform_faxis(nfft, fs);
  FP_TYPE* freqwarp = llsm_warp_freq(0, mlv0 -> conf.nosf, mlv0 -> conf.nnos, mlv0 -> conf.noswarp);
  FP_TYPE faxis_upper = faxis[ns - 1] * 2 - faxis[ns - 2];
  
  llsm_layer1* mlv1 = malloc(sizeof(llsm_layer1));
  mlv1 -> nfrm = nfrm;
  mlv1 -> nfft = nfft;
  mlv1 -> fnyquist = fs / 2;
  mlv1 -> vt_resp_magn = calloc(nfrm, sizeof(FP_TYPE*));
  mlv1 -> vs_har_ampl = calloc(nfrm, sizeof(FP_TYPE*));
  mlv1 -> vs_har_phse = calloc(nfrm, sizeof(FP_TYPE*));
  mlv1 -> lip_resp_phse = calloc(ns, sizeof(FP_TYPE));
  mlv1 -> lip_resp_magn = llsm_liprad(faxis, ns, 1.5, mlv1 -> lip_resp_phse); // 1.5cm mouth radius
  mlv1 -> vs_rd = llsm_analyze_rd_amponly(mlv0, faxis, mlv1 -> lip_resp_magn, ns);

  FP_TYPE** spectrogram = (FP_TYPE**)malloc2d(ns, nfrm, sizeof(FP_TYPE)); // ns * nfrm transposed matrix

# pragma omp parallel for
  for(int i = 0; i < nfrm; i ++) {
    FP_TYPE har_freq[1024];
    FP_TYPE f0 = mlv0 -> frames[i] -> f0;
    FP_TYPE* resp_magn = NULL;
    if(f0 > 0) {
      int nhar = mlv0 -> frames[i] -> sinu -> nhar;
      for(int j = 0; j < nhar; j ++) har_freq[j + 1] = (1.0 + j) * f0;
      har_freq[nhar + 1] = fs / 2; // [0 f0 2f0 ... nhar*f0 fs/2]
      
      FP_TYPE* lip_har_magn = interp1u(faxis[0], faxis_upper, mlv1 -> lip_resp_magn, ns, har_freq + 1, nhar);
      FP_TYPE* har_ampl = calloc(nhar, sizeof(FP_TYPE)); // [f0 2f0 3f0 ... nhar*f0]
      for(int j = 0; j < nhar; j ++)
        har_ampl[j] = mlv0 -> frames[i] -> sinu -> ampl[j] / lip_har_magn[j];

      // magnitude response without lip contribution
      resp_magn = llsm_harmonic_cheaptrick(har_ampl, nhar, f0, nfft, fs);
      
      free(lip_har_magn);
      free(har_ampl);
    } else {
      resp_magn = llsm_spectrum_from_envelope(freqwarp, mlv0 -> frames[i] -> noise -> spec,
        mlv0 -> conf.nnos, nfft, fs);
      for(int j = 0; j < ns; j ++)
        resp_magn[j] -= log(mlv1 -> lip_resp_magn[(j < 10 ? 10 : j)]);
    }
    for(int j = 0; j < ns; j ++) spectrogram[j][i] = resp_magn[j];
    free(resp_magn);
  }
  free(freqwarp);
  
# pragma omp parallel for
  for(int i = 0; i < nfrm; i ++) {
    FP_TYPE har_freq[1024];
    FP_TYPE f0 = mlv0 -> frames[i] -> f0;
    if(f0 > 0) {

    int nhar = mlv0 -> frames[i] -> sinu -> nhar;
    for(int j = 0; j <= nhar; j ++) har_freq[j] = j * f0;
    har_freq[nhar + 1] = fs / 2; // [0 f0 2f0 ... nhar*f0 fs/2]
    
    FP_TYPE* lip_har_magn = interp1u(faxis[0], faxis_upper, mlv1 -> lip_resp_magn, ns, har_freq + 1, nhar);
    FP_TYPE* lip_har_phse = interp1u(faxis[0], faxis_upper, mlv1 -> lip_resp_phse, ns, har_freq + 1, nhar);
    FP_TYPE* har_ampl = calloc(nhar, sizeof(FP_TYPE)); // [f0 2f0 3f0 ... nhar*f0]
    FP_TYPE* har_phse = calloc(nhar, sizeof(FP_TYPE)); // [f0 2f0 3f0 ... nhar*f0]
    for(int j = 0; j < nhar; j ++) {
      har_ampl[j] = mlv0 -> frames[i] -> sinu -> ampl[j] / lip_har_magn[j];
      har_phse[j] = mlv0 -> frames[i] -> sinu -> phse[j] - lip_har_phse[j]; // phase response on harmonics, without lip contribution
    }
    free(lip_har_phse);
    
    FP_TYPE* resp_magn = calloc(ns, sizeof(FP_TYPE)); // magnitude response without lip contribution
    for(int j = 0; j < ns; j ++) resp_magn[j] = spectrogram[j][i];
    lfmodel lf = lfmodel_from_rd(mlv1 -> vs_rd[i], 1.0 / f0, 1.0);
    FP_TYPE* vs_har_magn_ = lfmodel_spectrum(lf, har_freq + 1, nhar, NULL); // glottal flow derivative
    for(int j = 1; j < nhar; j ++) vs_har_magn_[j] /= (1.0 + j) * vs_har_magn_[0]; // integration + normalization
    vs_har_magn_[0] = 1.0;
    FP_TYPE* vs_har_magn = calloc(nhar + 2, sizeof(FP_TYPE)); // [0 f0 2f0 ... nhar*f0 fs/2]
    for(int j = 0; j < nhar; j ++) vs_har_magn[j + 1] = log(vs_har_magn_[j]);
    vs_har_magn[0] = vs_har_magn[1]; // fake DC
    vs_har_magn[nhar + 1] = vs_har_magn[nhar];
    free(vs_har_magn_);
    FP_TYPE* vs_resp_magn = interp1(har_freq, vs_har_magn, nhar + 2, faxis, ns);
    
    for(int j = 0; j < ns; j ++)
      resp_magn[j] -= vs_resp_magn[j]; // recovered vocal-tract (log) magnitude response
    for(int j = 0; j < nhar; j ++)
      har_ampl[j] /= exp(vs_har_magn[j + 1]); // recovered vocal-tract (linear) magnitude response on harmonics
    /* // alternatively, compute phase from interpolated VT response
    FP_TYPE* phaseresp = minphase(resp_magn, nfft);
    FP_TYPE* har_vt_phse = interp1(faxis, phaseresp, ns, har_freq + 1, nhar);*/
    FP_TYPE* har_vt_phse = llsm_harmonic_minphase(har_ampl, nhar);
    for(int j = 0; j < nhar; j ++)
      har_phse[j] -= har_vt_phse[j]; // recovered glottal flow phase response on harmonics
    
    mlv1 -> vt_resp_magn[i] = calloc(ns, sizeof(FP_TYPE));
    mlv1 -> vs_har_ampl[i] = calloc(nhar, sizeof(FP_TYPE));
    mlv1 -> vs_har_phse[i] = calloc(nhar, sizeof(FP_TYPE));
    for(int j = 0; j < ns; j ++) mlv1 -> vt_resp_magn[i][j] = resp_magn[j];
    for(int j = 1; j < nhar; j ++) {
      mlv1 -> vs_har_ampl[i][j] = exp(vs_har_magn[j + 1]);
      mlv1 -> vs_har_phse[i][j] = har_phse[j];
    }
    mlv1 -> vs_har_ampl[i][0] = exp(vs_har_magn[1]);
    mlv1 -> vs_har_phse[i][0] = har_phse[0];
    
    free(har_ampl);
    free(har_phse);
    free(vs_har_magn);
    free(har_vt_phse);
    free(lip_har_magn);
    free(vs_resp_magn);
    free(resp_magn);
    } // f0 > 0
  }
  
  free2d(spectrogram, ns);
  free(faxis);

  return mlv1;
}

void llsm_delete_layer1(llsm_layer1* dst) {
  if(dst == NULL) return;
  free2d(dst -> vt_resp_magn, dst -> nfrm);
  free2d(dst -> vs_har_ampl, dst -> nfrm);
  free2d(dst -> vs_har_phse, dst -> nfrm);
  free(dst -> lip_resp_magn);
  free(dst -> lip_resp_phse);
  if(dst -> vs_rd) free(dst -> vs_rd);
  free(dst);
}
