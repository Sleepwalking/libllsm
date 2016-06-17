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

/*
  Naming Conventions
    x         input time domain signal
    y         output time domain signal
    h         FIR filter coefficient
    nx/ny/... length of x/y/...
    f0        fundamental frequency
    fs        sampling frequency

  Memory Structure for FP_TYPE**
    Unless otherwise specified, the first dimension is time and the second is frequency.
    For example,
      llsm_sinparam* speech_sin;
      speech_sin -> freq = (FP_TYPE**)malloc2d(500, 80, sizeof(FP_TYPE)); // allocate 500 frames, each containing 80 sinusoids
      speech_sin -> freq[25]; // pointer to sinusoid frequencies in the 25th frame
      speech_sin -> freq[25][40]; // frequency of the 40th sinusoid (in this case, also harmonic) in the 25th frame
*/

#ifndef LLSM
#define LLSM

#ifndef FP_TYPE
  #define FP_TYPE double
#endif

/*
  llsm_sinparam: model parameters for sinusoidal signals
*/
typedef struct {
  FP_TYPE** freq; // Hz
  FP_TYPE** ampl; // linear
  FP_TYPE** phse; // rad
  int nfrm;
  int nhar;
} llsm_sinparam;

/*
  llsm_echannel: model parameters for each noise energy channel
*/
typedef struct {
  llsm_sinparam* eenv;  // sinusoidal parameters for noise energy envelope
  FP_TYPE* emin;        // minimum noise energy
} llsm_echannel;

/*
  llsm_conf: model configurations
*/
typedef struct {
  int nfrm;             // number of frames
  int nhop;             // hop size in samples
  int nhar;             // number of harmonics
  int nhare;            // number of harmonics for noise energy
  int nnos;             // size of noise spectrum
  FP_TYPE noswrap;      // wrapping factor for noise spectrum
  FP_TYPE mvf;          // maximum voiced frequency
  FP_TYPE nosf;         // upper bound of noise frequency
  FP_TYPE thop;         // hop size in seconds

  FP_TYPE* nosbandf;    // upper frequency of each noise band, excluding nosf; example: [2000, 4000, 8000]
  int nnosband;         // number of noise bands
} llsm_conf;

/*
  llsm: Low Level Speech Model paramters
*/
typedef struct {
  llsm_conf conf;       // configuration
  FP_TYPE* f0;          // fundamental frequency (Hz)

  llsm_sinparam* sinu;  // sinusoidal parameters
  FP_TYPE** noise;      // wrapped noise spectrogram

  llsm_echannel** nosch;// noise energy channels
} llsm;

/*
  llsm_parameters: configurations for analysis/synthesis processes

  note: there might be a naming confusion between llsm and llsm_parameters. The former is the model
    paramters (i.e. what we get after analyzing the speech); the later is a set of configurations for
    analysis/synthesis processes. We name so to make it consistent with libpyin (see pyin_paramters).
*/

typedef enum {
  qfft = 0, // use fft spectrum peak picking on quasi-harmonic analysis
  qhm // use Least-Square
} llsm_harmonic_analysis_method;

typedef struct {
  // params for analysis
  int a_nhop;           // hop size in samples
  int a_nhar;           // number of harmonics
  int a_nhare;          // number of harmonics for noise energy
  int a_nnos;           // size of noise spectrum
  FP_TYPE a_noswrap;    // wrapping factor for noise spectrum
  FP_TYPE a_tfft;       // size of transformation in seconds
  FP_TYPE a_mvf;        // maximum voiced frequency
  FP_TYPE a_nosf;       // maximum noise frequency
  FP_TYPE* a_nosbandf;  // upper frequency of each noise band
  int a_nnosband;       // number of noise bands

  llsm_harmonic_analysis_method a_method; // method of quasi-harmonic analysis, qfft by default
  char a_qhmlsmethod; // method of LS. 'Q' remains QR/LR(faster). 'S' remains SVD(better quality, very slow).
  int a_maxairiter;  // maximum number of air iteration, 16 by default
  int a_maxqhmiter; // maximum number of qhm iteration, 4 by default
  FP_TYPE a_maxqhmcorr; // maximum frequency correction per iteration(in Hz)

  // params for synthesis
  int s_fs;             // sampling frequency
} llsm_parameters;

/*
  llsm_init: obtain a configuration set initialized with default values
*/
llsm_parameters llsm_init(int nnosband);
void llsm_deinit(llsm_parameters dst);

/*
  llsm_analyze
  llsm_synthesize       their names faithfully reflect what they do
  llsm_delete
*/
llsm* llsm_analyze(llsm_parameters param, FP_TYPE* x, int nx, int fs, FP_TYPE* f0, int nf0);
FP_TYPE* llsm_synthesize(llsm_parameters param, llsm* model, int* ny);
void llsm_delete(llsm* model);

#endif
