# PPLN fixed-period study: numerical convergence and uncertainty record

**Study date:** 2026-09-19  
**Purpose:** establish a defensible numerical grid for the fixed-period PPLN pump-chirp comparison before running the full `GDD = 0 fs^2` scan.

## Scientific quantity of interest

The primary reported observable is the **integrated output energy between 950 and 1050 nm**, evaluated by integrating the LWE spectral energy density on its native frequency grid.

A secondary convergence diagnostic is the normalized L1 difference between the test and reference spectra over **700–1600 nm**.

The convergence acceptance criteria were fixed before the final production grid was selected:

- integrated 950–1050 nm energy error: **< 1%**
- normalized spectral L1 error over 700–1600 nm: **< 2%**

These are numerical convergence criteria, not statistical confidence intervals.

## Reference numerical grid

The high-resolution reference grid was:

- `dt = 0.45 fs`
- `dx = 6 um`
- `dz = 0.25 um`
- time window = `2.9988 ps`
- transverse window = `768 um`
- crystal length = `500 um`

At `dt = 0.45 fs`, the time grid has approximately 6664 samples.  
At `dx = 6 um`, the transverse grid has 128 samples.  
At `dz = 0.25 um`, a 500 um crystal is propagated using approximately 2000 longitudinal steps.

The convergence study was performed for the **pump GDD = 0 fs^2** case, because the transform-limited pump has the larger peak field and was therefore treated as the more numerically demanding case.

## Phase 1: one-parameter-at-a-time convergence

All tests below were performed at a fixed poling period of **21.75 um**, near the maximum of the existing fixed-period response.

| Run | dt (fs) | dx (um) | dz (um) | 950–1050 nm energy (pJ) | Energy difference (%) | Spectral L1 difference 700–1600 nm (%) | Result |
|---|---:|---:|---:|---:|---:|---:|---|
| Reference | 0.45 | 6 | 0.25 | 1195.211173 | 0 | 0 | reference |
| dz = 0.50 um | 0.45 | 6 | 0.50 | 1195.261693 | 0.004227 | 0.004679 | pass |
| dz = 1.00 um | 0.45 | 6 | 1.00 | 1196.323814 | 0.093092 | 0.127842 | pass |
| dx = 8 um | 0.45 | 8 | 0.25 | 1196.052354 | 0.070379 | 0.034240 | pass |
| dx = 12 um | 0.45 | 12 | 0.25 | 1198.443765 | 0.270462 | 0.131972 | pass |
| dt = 0.50 fs | 0.50 | 6 | 0.25 | 1193.409186 | 0.150767 | 0.041372 | pass |
| dt = 0.60 fs | 0.60 | 6 | 0.25 | 1193.462012 | 0.146347 | 0.042113 | pass |

All individual coarsenings passed comfortably.

## Phase 2: aggressive combined grid

An aggressive combined grid was tested:

- `dt = 0.60 fs`
- `dx = 12 um`
- `dz = 1.00 um`

Results:

| Poling period (um) | Reference energy (pJ) | Candidate energy (pJ) | Energy difference (%) | Spectral L1 difference (%) | Result |
|---:|---:|---:|---:|---:|---|
| 21.75 | 1195.211173 | 1197.793761 | 0.216078 | 0.266336 | pass |
| 23.00 | 600.206145 | 607.894975 | 1.281032 | 0.247577 | **fail** |

The candidate was rejected because the 23 um integrated-energy difference exceeded the pre-declared 1% criterion.

## Phase 3: refinement of the combined grid

At 23 um:

- `dt = 0.60 fs, dx = 8 um, dz = 1.00 um`: energy difference = **1.07905%**, spectral L1 difference = **0.15362%** — fail.
- `dt = 0.60 fs, dx = 6 um, dz = 1.00 um`: energy difference = **1.008003%** — fail by the pre-declared criterion. The broadband spectral error was not recorded for this intermediate test.
- `dt = 0.60 fs, dx = 6 um, dz = 0.50 um`: energy difference = **0.942253%**, spectral L1 difference = **0.041465%** — pass.

The final candidate was then independently checked at 21.75 um:

- reference = **1195.211173 pJ**
- candidate = **1193.509819 pJ**
- energy difference = **0.142348%**
- spectral L1 difference = **0.045162%**
- result = **pass**

## Accepted production grid

The accepted production grid for the `GDD = 0 fs^2` fixed-period scan is therefore:

- `dt = 0.60 fs`
- `dx = 6 um`
- `dz = 0.50 um`

This grid passed the pre-declared convergence criteria at both:

- **21.75 um**, near the strongest target-band response; and
- **23.0 um**, on the upper flank of the phase-matching response.

## Numerical uncertainty / sensitivity statement

For the production grid, the largest observed discrepancy in integrated 950–1050 nm energy relative to the fine reference was **0.942%** across the two final validation points.

A conservative way to communicate this in figures and text is:

> The numerical discretization sensitivity of the integrated 950–1050 nm energy is below 1% at the two validated poling periods.

If an uncertainty band is shown on a plot, use **±1% as a numerical-resolution sensitivity envelope**, not as a statistical 1-sigma uncertainty or confidence interval.

The corresponding broadband spectral differences for the accepted grid were only **0.041–0.045%** over 700–1600 nm.

## What this convergence study does and does not establish

It establishes that the chosen numerical grid reproduces the fine-grid result for the selected observable at two representative poling periods within the adopted tolerance.

It does **not** establish a complete experimental uncertainty budget. In particular, the following are separate sources of uncertainty or model dependence and should not be folded into the ±1% numerical envelope unless they are independently quantified:

- pump pulse energy;
- pump spectral bandwidth and spectral shape;
- pump GDD and higher-order phase;
- beam waist and spatial overlap;
- seed pulse parameters;
- nonlinear coefficients and material dispersion model;
- crystal length and poling-period fabrication uncertainty;
- model-form assumptions.

## Important model limitation

The deterministic LWE calculations use a **coherent MIR seed** as a surrogate for vacuum-initiated OPG. Therefore the absolute simulated conversion/output should not automatically be interpreted as a quantitative prediction of spontaneous OPG yield. The fixed-period and chirp scans are most defensible as controlled comparisons within the specified deterministic model.

## Poling-period sampling

The scan uses a broad 17–30 um range and local refinement around the interesting region. Around 20–24 um the final sampling interval is **0.125 um**.

This is a parameter-sampling resolution, not a vertical error bar. If the optimum is quoted from the sampled maximum without fitting, report the sampled poling period and the 0.125 um sampling interval rather than inventing a statistical uncertainty.

## Report-ready methods wording

> Numerical convergence was assessed before performing the pump-chirp parameter scan. A fine reference grid with temporal, transverse and longitudinal step sizes of 0.45 fs, 6 um and 0.25 um, respectively, was compared with progressively coarser grids. Convergence was assessed using both the integrated 950–1050 nm output energy and the normalized spectral L1 difference over 700–1600 nm. Acceptance thresholds of <1% and <2%, respectively, were specified before selecting the production grid. The final grid, dt = 0.60 fs, dx = 6 um and dz = 0.50 um, differed from the fine-grid calculation by 0.14% at a 21.75 um poling period and 0.94% at 23.0 um, while the corresponding broadband spectral differences were below 0.05%. This grid was therefore used for the GDD = 0 fs^2 production scan.

## Figure-caption / uncertainty wording

Short form:

> Shaded uncertainty, where shown, indicates a conservative ±1% numerical-resolution sensitivity derived from comparison with a finer reference grid; it is not a statistical confidence interval.

Longer form:

> Numerical-resolution sensitivity was evaluated at 21.75 and 23.0 um by comparison with a finer grid (dt = 0.45 fs, dx = 6 um, dz = 0.25 um). The accepted production grid (dt = 0.60 fs, dx = 6 um, dz = 0.50 um) changed the integrated 950–1050 nm energy by 0.14% and 0.94%, respectively. A conservative ±1% envelope is therefore used where a numerical-resolution band is displayed.

## Reproducibility rule

Do not rely on chat history as the long-term record. The durable record should be the Git repository plus archived raw outputs.

Track in Git:

- all input generators;
- exact input files or deterministic generator configuration;
- `run_manifest.csv`;
- convergence analysis scripts;
- convergence summary CSVs;
- this methods note;
- final processed CSVs;
- final figures;
- figure/report notes;
- runtime/version provenance.

Raw LWE output ZIPs may remain outside Git if large, but should be retained in an immutable archive together with a checksum manifest.
