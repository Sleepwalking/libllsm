libllsm
===

Low Level Speech Model for high quality speech analysis and synthesis.

LLSM is an extended *Harmonic + Noise Model* with turbulent noise modeling, shape-invariant phase preservation and a range of other techniques combined.

It is designed to be an **interpolatable**, **low-level** speech representation for **high quality** speech/singing synthesis, which means

* we're willing to sacrifice some storage consumption in exchange for better quality;
* we assume the noise/distortion in the input is negligible;
* data frames can be easily interpolated with little quality loss;
* the model deals with really low-level aspects in speech parameterization;
* where "low-level aspects" can be amplitudes of harmonics, f0, and spectral envelope; semantics of the context would be a counterexample;

License: [UIUC license](https://en.wikipedia.org/wiki/University_of_Illinois/NCSA_Open_Source_License)

Acknowledgement
---

The author would like to thank Professor Jont Allen and HSR Group at Univ. of Illinois for their kind and constructive advice.

Works cited
---

* Pantazis, Yannis, and Yannis Stylianou. "Improving the modeling of the noise part in the harmonic plus noise model of speech." Acoustics, Speech and Signal Processing, 2008. ICASSP 2008. IEEE International Conference on. IEEE, 2008.

* Röbel, Axel, and Xavier Rodet. "Efficient spectral envelope estimation and its application to pitch shifting and envelope preservation." Proc. DAFx. 2005.

* Serra, Xavier. "A system for sound analysis/transformation/synthesis based on a deterministic plus stochastic decomposition." Diss. Universitat Pompeu Fabra. 1989.

* Stylianou, Yannis. "Harmonic plus noise models for speech, combined with statistical methods, for speech and speaker modification." Diss. Ecole Nationale Supérieure des Télécommunications. 1996.
