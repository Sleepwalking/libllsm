libllsm
===

Low Level Speech Model for high-quality speech analysis/synthesis

About
---

libllsm is a C library providing data structures and routines for analysis (parametrization), modification and synthesis of digital speech.

### The model

LLSM is a two-layer model of speech. The first layer (layer 0) is a signal-level parametrization that separately models harmonic and noise (aspiration and consonants) components; the second layer (layer 1) is an acoustic-level parametrization that further decomposes, in an approximated manner, the harmonic component by different parts of speech production system (glottis - vocal tract - lip).

* libllsm can be viewed as a lossy speech coder, but it differs from conventional speech coders in being suitable for modifications (rather than compression).

Compiling
---

1. `mkdir external build`
2. create symbolic link to [libpyin](https://github.com/Sleepwalking/libpyin) under `external`
3. run `make`, `make test`

Note: there's a macro named `FP_TYPE`, which is either float or double, that has to be specified as a complier flag (i.e., `-DFP_TYPE=float`).

The test speech `test/arctic_a0001.wav` is a sample taken from the CMU Arctic database.

How to use
---

`test/test.c` is a bare-bones example of doing pitch shifting with libllsm.

Once complied, run it with the following command.

`./build/llsm-test path-to-wav-file pitch-shift-ratio`

where `pitch-shift-ratio` is a positive number for scaling the fundamental frequency.

### Pitch shifting with libllsm: a code walk-through

The following code initializes analysis parameters for libllsm. The only argument of `llsm_init` is the number of bands for noise excitation (the way LLSM models noise excitation is basically a multi-band extension to [2]).
```c
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
```

Given some F0 estimation stored in a float-point array `f0`, we first call `llsm_layer0_analyze` on the input signal to obtain the layer 0 representation.
```c
  llsm_layer0* model = llsm_layer0_analyze(lparam, x, nx, fs, f0, nfrm, NULL);
```

The harmonic parameters in the layer 0 model are in absolute phase, which is somewhat inconvenient to manipulate (consider pitch shifting or interpolation). We apply a time shift to each frame so that the phases are made relative to the first harmonic. This is called Relative Phase Shift (RPS) [3].
```c
  FP_TYPE* phase0 = calloc(nfrm, sizeof(FP_TYPE));
  for(int i = 0; i < nfrm; i ++)
    phase0[i] = model -> frames[i] -> f0 > 0 ? -model -> frames[i] -> sinu -> phse[0] : 0;
  llsm_layer0_phaseshift(model, phase0);
  free(phase0);
```

Then we go from layer 0 up to layer 1 by calling `llsm_layer1_from_layer0`. The layer 1 model contains separate information about vocal tract and source. For pitch shifting we simply need to keep the layer 1 model intact and resample the vocal tract transfer function at scaled harmonic frequencies. Note that the layer 1 model is for harmonic component only, so the layer 0 model is still relevant and it should not be discarded at this point.
```c
  llsm_layer1* model_lv1 = llsm_layer1_from_layer0(lparam, model, 2048, fs);
```

Next in the for loop over frames, we first scale the F0 and make an array of harmonic frequencies.
```c
    llsm_frame* iframe = model -> frames[i];
    FP_TYPE origf0 = iframe -> f0;
    iframe -> f0 *= argc > 2 ? atof(argv[2]) : 1.0;
    FP_TYPE freq[512];
    for(int j = 0; j < iframe -> sinu -> nhar; j ++) {
      freq[j] = iframe -> f0 * (j + 1.0);
      if(freq[j] > lparam.a_mvf) {
        iframe -> sinu -> nhar = j;
        break;
      }
    }
```

Then as described above, vocal tract and lip frequency responses are subsampled at new frequencies; the vocal tract phase response is computed from amplitudes under minimum phase assumption.
```c
    FP_TYPE* newampl = interp1(faxis, model_lv1 -> vt_resp_magn[i], nfft / 2 + 1, freq, iframe -> sinu -> nhar);
    for(int j = 0; j < iframe -> sinu -> nhar; j ++) newampl[j] = exp(newampl[j]);
    FP_TYPE* newphse = llsm_harmonic_minphase(newampl, iframe -> sinu -> nhar);
    FP_TYPE* lipampl = interp1(faxis, model_lv1 -> lip_resp_magn, nfft / 2 + 1, freq, iframe -> sinu -> nhar);
    FP_TYPE* lipphse = interp1(faxis, model_lv1 -> lip_resp_phse, nfft / 2 + 1, freq, iframe -> sinu -> nhar);
```

These parts are re-combined by multiplication and the result is written back to layer 0. The `origf0 / iframe -> f0` term compensates for the amplitude gain due to change in number of harmonics within audible range.
```c
    for(int j = 0; j < iframe -> sinu -> nhar; j ++) {
      iframe -> sinu -> ampl[j] = model_lv1 -> vs_har_ampl[i][j] * newampl[j] * lipampl[j] * origf0 / iframe -> f0;
      iframe -> sinu -> phse[j] = model_lv1 -> vs_har_phse[i][j] + newphse[j] + lipphse[j];
    }
```

At this point pitch shifting is done on layer 0. But before synthesis we first need to recover the phase progression along time axis (which is the integral of F0).
```c
  phase0 = calloc(nfrm, sizeof(FP_TYPE));
  for(int i = 1; i < nfrm; i ++)
    if(model -> frames[i] -> f0 > 0)
      phase0[i] = phase0[i - 1] + model -> frames[i] -> f0 * nhop / fs * 2 * M_PI;
  llsm_layer0_phaseshift(model, phase0);
  free(phase0);
```

Finally we call `llsm_layer0_synthesize` to convert from layer 0 model back to signal. The output is a structure containing harmonic and noise components of the synthesized speech.
```c
  llsm_output* out = llsm_layer0_synthesize(lparam, model);
```

Licensing
---

libllsm is licensed under GPLv3.

I have a pending patent on LLSM-related technology. However the patent license is granted to libllsm users, free from royalty, under the terms of GPLv3.

Please contact the author for an alternatively licensed version primarily for commercial purposes.

Publications
---

Currently there's no publication directly associated with LLSM. However there is [a poster](https://khua5.web.engr.illinois.edu/writings/hua-spcc-poster.pdf) on the pseudo glottal inverse filtering method in layer 1 LLSM.

K. Hua, "Speech Analysis/Synthesis by Non-parametric Separation of Vocal Source and Tract Responses," presented at Speech Processing Courses in Crete, 2016.

The following are the major publications that LLSM draws inspiration from.

1. G. Degottex, P. Lanchantin, A. Roebel, and X. Rodet, "Mixed source model and its adapted vocal tract filter estimate for voice transformation and synthesis," Speech Communication, vol. 55, no. 2, pp. 278–294, 2013.

2. Y. Pantazis and Y. Stylianou, "Improving the modeling of the noise part in the harmonic plus noise model of speech," 2008 IEEE International Conference on Acoustics, Speech and Signal Processing, 2008.

3. I. Saratxaga, Hernáez I., D. Erro, E. Navas, and Sánchez J., "Simple representation of signal phase for harmonic speech models," Electronics Letters, vol. 45, no. 7, p. 381, 2009.
