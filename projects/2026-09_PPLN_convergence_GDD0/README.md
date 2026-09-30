# GDD=0 numerical convergence study

Reference physics:
- fixed QPM period: 21.75 um
- pump GDD: 0 fs^2
- all other optical/material parameters copied verbatim from the current
  PPLN_FIXED_21p75um production input.

Resolution tests change one numerical parameter at a time from:
- dt = 0.45 fs
- dx = 6 um
- dz = 0.25 um

Acceptance targets used by analysis:
- integrated 950-1050 nm energy error < 1%
- normalized spectral L1 error over 700-1600 nm < 2%

Do not use a coarser production grid until a combined-grid validation has
also been run after the one-at-a-time study.
