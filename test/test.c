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

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "../external/libpyin/pyin.h"
#include "../math-funcs.h"
#include "../envelope.h"
#include "../llsm.h"

#define FP_TYPE float

char* strdup(const char* s);

double get_time() {
  struct timeval t;
  gettimeofday(&t, NULL);
  return (t.tv_sec + (t.tv_usec / 1000000.0)) * 1000.0;
}

int main(int argc, char** argv) {
  if(argc < 2) {
    fprintf(stderr, "Missing argument.\n");
    return 1;
  }
  
  int fs = 0;
  int nbit = 0;
  int nx = 0;
  int nfrm = 0;
  int nhop = 256;
  FP_TYPE* f0 = NULL;
  FP_TYPE* x = wavread(argv[1], & fs, & nbit, & nx);
  
  llsm_parameters lparam = llsm_init(4);
  lparam.a_nosbandf[0] = 2000;
  lparam.a_nosbandf[1] = 5000;
  lparam.a_nosbandf[2] = 9000;
  lparam.a_mvf = 12000;
  lparam.a_nnos = 192;
  lparam.a_nosf = fs / 2;
  lparam.a_nhop = nhop;
  lparam.a_nhar = 400;
  lparam.a_nhare = 5;
  
  if(argc < 4) {
    pyin_config param = pyin_init(nhop);
    param.fmin = 50.0;
    param.fmax = 900.0;
    param.trange = 24;
    param.nf = ceil(fs * 0.025);
    param.w = param.nf / 4;
    param.bias = 10;
    f0 = pyin_analyze(param, x, nx, fs, & nfrm);
  } else {
    // load f0 from binary file
    FILE* fp = fopen(argv[3], "rb");
    float tmp;
    fread(& tmp, 4, 1, fp); // nfrm
    nfrm = tmp;
    f0 = calloc(nfrm, sizeof(FP_TYPE));
    for(int i = 0; i < nfrm; i ++) {
      fread(& tmp, 4, 1, fp); // f0[i]
      f0[i] = tmp;
    }
    fread(& tmp, 4, 1, fp); // naxis[0]
    nhop = tmp;
    fread(& tmp, 4, 1, fp); // naxis[1]
    nhop = tmp - nhop;
    fclose(fp);
    lparam.a_f0refine = 0;
  }

  printf("Analyzing (lv0)...\n");
  llsm_layer0* model = llsm_layer0_analyze(lparam, x, nx, fs, f0, nfrm, NULL);
  
  FP_TYPE* phase0 = calloc(nfrm, sizeof(FP_TYPE));
  for(int i = 0; i < nfrm; i ++)
    phase0[i] = model -> frames[i] -> f0 > 0 ? -model -> frames[i] -> sinu -> phse[0] : 0;
  llsm_layer0_phaseshift(model, phase0);
  free(phase0);

  printf("Analyzing (lv1)...\n");
  llsm_layer1* model_lv1 = llsm_layer1_from_layer0(lparam, model, 2048, fs);

  printf("\nModifying...\n");
  FP_TYPE* faxis = llsm_uniform_faxis(2048, lparam.a_nosf * 2);
  FP_TYPE* freqwarp = llsm_warp_freq(0, lparam.a_nosf, lparam.a_nnos, lparam.a_noswarp);
  for(int i = 0; i < model -> conf.nfrm; i ++) {
    llsm_frame* iframe = model -> frames[i];
    if(iframe -> f0 <= 0) continue;
    const int nfft = 2048;
    
    FP_TYPE origf0 = iframe -> f0;
    iframe -> f0 *= argc > 2 ? atof(argv[2]) : 1.0;
    //iframe -> f0 += f0_vel;
    FP_TYPE freq[512];
    for(int j = 0; j < iframe -> sinu -> nhar; j ++) {
      freq[j] = iframe -> f0 * (j + 1.0);
      if(freq[j] > lparam.a_mvf) {
        iframe -> sinu -> nhar = j;
        break;
      }
    }
    FP_TYPE* newampl = interp1(faxis, model_lv1 -> vt_resp_magn[i], nfft / 2 + 1, freq, iframe -> sinu -> nhar);
    for(int j = 0; j < iframe -> sinu -> nhar; j ++) newampl[j] = exp(newampl[j]);
    FP_TYPE* phaseresp = minphase(model_lv1 -> vt_resp_magn[i], nfft);
    FP_TYPE* newphse = llsm_harmonic_minphase(newampl, iframe -> sinu -> nhar);
    FP_TYPE* lipampl = interp1(faxis, model_lv1 -> lip_resp_magn, nfft / 2 + 1, freq, iframe -> sinu -> nhar);
    FP_TYPE* lipphse = interp1(faxis, model_lv1 -> lip_resp_phse, nfft / 2 + 1, freq, iframe -> sinu -> nhar);

    for(int j = 0; j < iframe -> sinu -> nhar; j ++) {
      iframe -> sinu -> ampl[j] = model_lv1 -> vs_har_ampl[i][j] * newampl[j] * lipampl[j] * origf0 / iframe -> f0;
      iframe -> sinu -> phse[j] = model_lv1 -> vs_har_phse[i][j] + newphse[j] + lipphse[j];
    }

    free(newampl);
    free(newphse);
    free(phaseresp);
    free(lipampl);
    free(lipphse);
  }
  free(faxis); free(freqwarp);
  llsm_delete_layer1(model_lv1);
  
  phase0 = calloc(nfrm, sizeof(FP_TYPE));
  for(int i = 1; i < nfrm; i ++)
    if(model -> frames[i] -> f0 > 0)
      phase0[i] = phase0[i - 1] + model -> frames[i] -> f0 * nhop / fs * 2 * M_PI;
  llsm_layer0_phaseshift(model, phase0);
  free(phase0);

  lparam.s_fs = fs;

  printf("Synthesizing...\n");
  double t0 = get_time();
  llsm_output* out = llsm_layer0_synthesize(lparam, model);
  double tspent = get_time() - t0;
  printf("%f ms, %fs/s\n", tspent, (FP_TYPE)nx / fs / tspent * 1000);
  
  wavwrite(out -> y, out -> ny, lparam.s_fs, 16, "resynth.wav");

  llsm_deinit(lparam);
  llsm_delete_layer0(model);
  free(f0);
  free(x); // free(y);
  llsm_delete_output(out);

  return 0;
}

