# Parameter comparison: reconstructed scan vs archived report-model input

Source: `reference/ChirpedPPLN_18_24um_n1230fs2_0d1uJ.txt` (LWE 2025.5).

| Parameter | Previous reconstructed scan | Archived report-model input | Revised scan |
|---|---:|---:|---:|
| Pump energy | 0.025 uJ | 0.100 uJ | 0.100 uJ |
| Pump centre | 374.74 THz (~800 nm) | 374.74 THz (~800 nm) | same |
| Pump bandwidth parameter | 46.84 THz | 25.00 THz | 25.00 THz |
| Pump SG order | 2 | 4 | 4 |
| Pump GDD | -500 fs^2 | -1230 fs^2 | -1230 fs^2 |
| Pump waist | 200 um | 200 um | 200 um |
| Seed energy | 10 nJ | 10 nJ | 10 nJ |
| Seed centre | 56.00 THz (~5.353 um) | 74.95 THz (~4.000 um) | 74.95 THz |
| Seed bandwidth parameter | 75.00 THz | 56.21 THz | 56.21 THz |
| Seed SG order | 2 | 10 | 10 |
| Seed GDD | 0 | 0 | 0 |
| Seed waist | 150 um | 200 um | 200 um |
| Grid width | 768 um | 768 um | same |
| dx | 6 um | 6 um | same |
| Time window | 2.9988 ps | 2.9988 ps | same |
| dt | 0.45 fs | 0.45 fs | same |
| dz | 0.25 um | 0.25 um | same |
| Material | MgO:LiNbO3 | MgO:LiNbO3 | same |

## Sequence observation

Despite its filename, the archived file does not chirp the domain size within
one propagation. It uses `i34` as a constant domain thickness inside the loop,
while Batch mode 34 scans parameter 34 (crystal thickness) from 9 to 12 um in
13 steps. This corresponds to uniform first-order QPM periods from 18.0 to
24.0 um in 0.5 um steps.

The revised scan therefore uses the same 13 full periods, with one input file
per period. It keeps the physical crystal length exactly 500 um for every run;
where necessary only the final ferroelectric domain is truncated at the exit
face.
