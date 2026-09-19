# Pump-chirp comparison: figure and report notes

## Intended main figure

x-axis:
- fixed PPLN poling period (um)

y-axis:
- integrated output energy in 950–1050 nm (pJ)

Curves:
- pump GDD = -1230 fs^2
- pump GDD = 0 fs^2

Do not include the former peak-wavelength secondary axis; it was removed because it did not provide a useful diagnostic for the scientific question.

## Numerical-resolution annotation

For the GDD = 0 production curve, the accepted numerical grid is:

- dt = 0.60 fs
- dx = 6 um
- dz = 0.50 um

Fine-grid validation:
- 21.75 um: 0.142348% target-band energy difference; 0.045162% spectral L1 difference
- 23.0 um: 0.942253% target-band energy difference; 0.041465% spectral L1 difference

If desired, show a conservative ±1% band around the GDD = 0 curve and describe it explicitly as a **numerical-resolution sensitivity envelope**, not a statistical uncertainty.

The -1230 fs^2 curve was generated on the finer numerical grid. This convergence study was performed explicitly for the GDD = 0 case. Do not claim a separately measured ±1% statistical uncertainty for the chirped curve from this study alone.

## Sampling note

The region 20–24 um is sampled at 0.125 um poling-period intervals. Outside this locally refined region the scan is coarser.

The horizontal spacing is sampling resolution, not an error bar on the physical poling period.

## Interpretation cautions

- The comparison isolates the effect of pump spectral phase/GDD within the deterministic LWE model.
- GDD changes temporal duration and peak field while leaving the initial spectral-energy magnitude unchanged.
- The model uses a coherent MIR seed as a surrogate for vacuum-initiated OPG.
- Absolute output should therefore be interpreted cautiously; relative trends between otherwise matched simulations are more robust.
