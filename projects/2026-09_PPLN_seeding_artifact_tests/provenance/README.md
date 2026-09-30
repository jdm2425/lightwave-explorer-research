# PPLN coherent-seed artefact diagnostics

All runs use the GDD=0 pump/seed parameters from the matched PPLN study.

The tests are:

1. Exact-length 500 um single-entrance-seed controls for 20->24 um and
   24->20 um chirps.

2. First-half tests:
   - 20->22 um over exactly 250 um
   - 24->22 um over exactly 250 um

3. Distributed-seed tests:
   - 20->24 um over exactly 500 um
   - 24->20 um over exactly 500 um

For the distributed test the original 10 nJ coherent MIR seed is divided
equally between one injection per full QPM period. The first fractional seed
is present at the entrance through the normal interface pulse. A fresh
fractional seed is then added at the start of every subsequent QPM period with
LWE addPulse(). Thus the nominal injected pulse energies sum to 10 nJ.

This is a diagnostic model, not a model of vacuum-seeded OPG. The injected
fields remain deterministic and coherent. In particular, coherent field
interference means that "10 nJ total nominal injected energy" does not imply
that the field energy at the output must equal the arithmetic sum of injection
energies.

One injection per full QPM period was chosen rather than every half-period.
This avoids tying the injection itself to each sign flip of d_eff while still
sampling the seed throughout the grating.

The grating sequences in this project use an exact physical length. Domain
boundaries are generated from a linearly varying local period Lambda(z) using
the accumulated grating phase, and the final domain is truncated at exactly
250 or 500 um.
