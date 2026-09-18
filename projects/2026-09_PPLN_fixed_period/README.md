# ESA PPLN fixed-period scan

This is the first project in the reusable Docker-backed LWE workspace.

## Scientific scan

- fixed first-order QPM periods: 17–30 µm in 1 µm increments
- crystal length: exactly 500 µm
- pump: ~800 nm, 0.025 µJ, -500 fs² GDD
- artificial coherent MIR surrogate seed retained from the report model
- target NIR metric: integrated 950–1050 nm output plus peak position in that band

The purpose is to isolate the effect of **fixed poling period**, instead of
inferring the trend only from different chirped-period ranges.

Important modelling limitation: the exact original LWE input used to generate
the report's existing chirped-PPLN figure was not available, so these inputs
reconstruct the stated report conditions using the archived MgO:LiNbO3
numerical conventions. Do not interpret absolute output energy as a direct
prediction of pump-only OPG.

## Run

First smoke test:

```zsh
lwe one inputs/PPLN_FIXED_17um.txt
```

Then:

```zsh
lwe run
```

or, after some completed runs:

```zsh
lwe run --remaining
```

Check:

```zsh
lwe status
```

Analyse:

```zsh
lwe analyze
```

Package for review:

```zsh
lwe pack --include-inputs
```
