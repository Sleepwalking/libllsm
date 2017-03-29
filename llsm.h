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

/*
  Naming Convention
    x         input time domain signal
    y         output time domain signal
    h         FIR filter coefficient
    nx/ny/... length of x/y/...
    f0        fundamental frequency
    fs        sampling frequency
*/

#ifndef LLSM
#define LLSM

/*
  llsm_sinframe: model parameters for a single frame of sinusoidal signal
*/
typedef struct {
  FP_TYPE* ampl; // linear
  FP_TYPE* phse; // rad
  int nhar;
} llsm_sinframe;

/*
  llsm_nosframe: model parameters for a single frame of noise signal
*/
typedef struct {
  llsm_sinframe** eenv; // noise energy envelope for each channel
  FP_TYPE* emin;        // minimum noise energy contour for each channel
  FP_TYPE* spec;        // noise spectral envelope
  int nchannel;
  int nnos;
} llsm_nosframe;

/*
  llsm_frame: model parameters for a single frame
*/
typedef struct {
  llsm_sinframe* sinu;  // harmonic component
  llsm_nosframe* noise; // stochastic component
  FP_TYPE f0;           // in Hz
} llsm_frame;

/*
  llsm_conf: model configurations
*/
typedef struct {
  int nfrm;             // number of frames
  int nhop;             // hop size in samples
  int nhar;             // number of harmonics
  int nhare;            // number of harmonics for noise energy
  int nnos;             // size of noise spectrum
  FP_TYPE noswarp;      // warpping factor for noise spectrum
  FP_TYPE mvf;          // maximum voiced frequency
  FP_TYPE nosf;         // upper bound of noise frequency
  FP_TYPE thop;         // hop size in seconds
  
  FP_TYPE* nosbandf;    // upper frequency of each noise band, excluding nosf;
                        //   example: [2000, 4000, 8000], size = nnosband - 1
  int nnosband;         // number of noise bands
} llsm_conf;

/*
  llsm_layer0: signal-domain model paramters
  llsm_layer1: pseudo-acoustic domain model paramters

  Note: llsm_layer1 is not complete on its own as it does not include the noise parameters.
    Only partial conversion from llsm_layer1 down to llsm_layer0 is possible unless noise parameters
    are retained.
*/
typedef struct {
  llsm_conf conf;       // configuration

  llsm_frame** frames;  // llsm -> frames[n] -> sinu  -> ampl[k] | phse[k] | f0 | nhar
                        //                      noise -> eenv[n] | nchannel -> ampl[k] | phse[k] | f0 | nhar
                        //                               emin[n]
                        //                               spec[k]
} llsm_layer0;

typedef struct {
  FP_TYPE** vt_resp_magn;     // vocal tract (log) magnitude response
  FP_TYPE** vs_har_ampl;      // harmonic amplitudes without vocal tract contribution
  FP_TYPE** vs_har_phse;      // harmonic phases without vocal tract contribution
  FP_TYPE*  vs_rd;            // Rd parameterization of LF model
  FP_TYPE*  lip_resp_magn;    // lip (linear) magnitude response
  FP_TYPE*  lip_resp_phse;    // lip phase response
  int nfrm;
  int nfft;
  int fnyquist;
} llsm_layer1;

/*
  llsm_parameters: configurations for analysis/synthesis processes

  note: there might be a naming confusion between llsm and llsm_parameters. The former is the model
    paramters (i.e. what we get after analyzing the speech); the later is a set of configurations for
    analysis/synthesis processes. We name so to make it consistent with libpyin (see pyin_paramters).
*/
typedef struct {
  // params for analysis
  int a_nhop;           // hop size in samples
  int a_nhar;           // number of harmonics
  int a_nhare;          // number of harmonics for noise energy
  int a_nnos;           // size of noise spectrum
  int a_f0refine;       // turn on/off f0 refinement
  FP_TYPE a_wsize;      // window size relative to period length
  FP_TYPE a_noswarp;    // warpping factor for noise spectrum
  FP_TYPE a_tfft;       // size of transformation in seconds
  FP_TYPE a_mvf;        // maximum voiced frequency
  FP_TYPE a_nosf;       // maximum noise frequency
  FP_TYPE* a_nosbandf;  // upper frequency of each noise band
  int a_nnosband;       // number of noise bands

  // params for synthesis
  int s_fs;             // sampling frequency
  int s_noiseonly;      // synthesize noise component only
  int s_n0;             // range (starts from)
  int s_n1;             // range (ends at)
} llsm_parameters;

typedef struct {
  // synthesis results
  FP_TYPE* y;
  FP_TYPE* y_sin;
  FP_TYPE* y_nos;
  int ny;
  FP_TYPE fs;
} llsm_output;

// === LLSM Layer 0 functions ===

/*
  llsm_create_frame: create a llsm parameter frame given number of harmonics,
    number of noise energy harmonics, size of noise spectrum, number of noise energy channels
  llsm_delete_frame: free a llsm parameter frame
*/
llsm_frame* llsm_create_frame(int nhar, int nhare, int nnos, int nchannel);
llsm_layer0* llsm_create_empty_layer0(llsm_conf conf);
void llsm_copy_frame(llsm_frame* dst, llsm_frame* src);
void llsm_copy_sinframe(llsm_sinframe* dst, llsm_sinframe* src);
void llsm_copy_nosframe(llsm_nosframe* dst, llsm_nosframe* src);
void llsm_delete_frame(llsm_frame* dst);

/*
  llsm_init: obtain a configuration set initialized with default values
*/
llsm_parameters llsm_init(int nnosband);
void llsm_deinit(llsm_parameters dst);

/*
  llsm_analyze
  llsm_synthesize           their names faithfully reflect what they do
  llsm_delete
*/
llsm_layer0* llsm_layer0_analyze(llsm_parameters param, FP_TYPE* x, int nx, int fs,
  FP_TYPE* f0, int nf0, FP_TYPE** xap);
llsm_output* llsm_layer0_synthesize(llsm_parameters param, llsm_layer0* model);
void llsm_delete_layer0(llsm_layer0* dst);
void llsm_delete_output(llsm_output* dst);

// === LLSM Layer 1 functions ===

/*
  llsm_uniform_faxis        generate a (nfft / 2 + 1)-sized array of uniformly spaced frequency values
  llsm_liprad               generate the approximated frequency response of lip radiation with radius
                              in cm (the default radius is 1.5 cm)
  llsm_harmonic_cheaptrick  estimate spectral envelope from harmonics using modified Cheaptrick algorithm
  llsm_harmonic_minphase    compute the phases of harmonics from amplitudes, assuming minimum phase
*/
FP_TYPE* llsm_uniform_faxis(int nfft, int fs);
FP_TYPE* llsm_liprad(FP_TYPE* freq, int nf, FP_TYPE radius, FP_TYPE* dst_phaseresp);
FP_TYPE* llsm_harmonic_cheaptrick(FP_TYPE* ampl, int nhar, FP_TYPE f0, int nfft, FP_TYPE fs);
FP_TYPE* llsm_harmonic_minphase(FP_TYPE* ampl, int nhar);

/*
  llsm_layer1_from_layer0   convert from layer 0 representation to layer 1 representation
  llsm_delete_layer1        deallocate memeory
*/
llsm_layer1* llsm_layer1_from_layer0(llsm_parameters param, llsm_layer0* model, int nfft, int fs);
void llsm_delete_layer1(llsm_layer1* dst);

void llsm_layer0_phaseshift(llsm_layer0* dst, FP_TYPE* phaseshift);

#endif
