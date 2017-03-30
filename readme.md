libllsm
===

Low Level Speech Model for high-quality speech analysis/synthesis

About
---

libllsm is a C library providing data structures and routines for analysis (parametrization), modification and synthesis of digital speech.

### The model

LLSM is a two-layer model of speech. The first layer (layer 0) is a signal-level parametrization that separately models harmonic and noise (aspiration and consonants) components; the second layer (layer 1) is an acoustic-level parametrization that further decomposes, in an approximated manner, the harmonic component by different parts of speech production system (glottis - vocal tract - lip).

* libllsm can be viewed as a lossy speech coder, but it differs from conventional speech coders in being suitable for modifications (rather than compression).

How to use
---

`test/test.c` is a bare-bones example of doing pitch shifting with libllsm.

Once complied, run it with the following command.

`./build/llsm-test path-to-wav-file pitch-shift-ratio`

where `pitch-shift-ratio` is a positive number for scaling the fundamental frequency.

Licensing
---

libllsm is licensed under GPLv3.

I have a pending patent on LLSM-related technology. However the patent license is granted to libllsm users, free from royalty, under the terms of GPLv3.

Please contact the author for an alternatively licensed version primarily for commercial purposes.
