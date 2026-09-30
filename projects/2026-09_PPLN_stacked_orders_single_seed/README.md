# Rerun both stacked PPLN orders with one MIR seed each

## What these inputs model

Each run launches a 100 nJ, zero-GDD pump and a **single 10 nJ coherent MIR
seed at the entrance**. The input contains precisely one `init()` and no
`addPulse()` call. The pump, seed and generated fields evolve continuously
through 500 µm of the first chirp and 500 µm of the second, on the same
numerical grid. Walk-off is included naturally; nothing is injected to undo
it at the midpoint.

| LWE input | Grating order | MIR seed injections |
| --- | --- | ---: |
| `STACK_NO_RESEED.txt` | 20→24→24→20 µm | 1 at entrance |
| `REVERSE_STACK_SINGLE_SEED.txt` | 24→20→20→24 µm | 1 at entrance |

The first input is the *exact archived input* for the previously plotted
stack. `generate_reverse_input.py` makes the second by concatenating the
validated exact-length 24→20 and 20→24 controls in the opposite order. It
verifies that concatenating those controls in their original order reproduces
the archived forward stack exactly. Each run has 92 nonlinear domains in
total and a physical length of 1 mm.

This is contiguous numerical propagation: there is no air gap, Fresnel loss,
alignment change or independently reset grating phase at the midpoint. These
are deterministic coherent-seed calculations, not a stochastic OPG model.

## Where the original output may be on your Mac

The original project installer defaulted to:

```text
~/LightwaveExplorer/lightwave-explorer-research/projects/2026-09_PPLN_two_stage_reseed/
```

Look in `inputs/STACK_NO_RESEED.txt` and `outputs/STACK_NO_RESEED.zip`.
The supplied report bundle contained only the extracted
`spectra_only/STACK_NO_RESEED_spectrum.dat`, not the full LWE ZIP. No
reverse-order output was present in the supplied bundles, although one
could exist elsewhere on your Mac. To search the likely projects root:

```zsh
find "$HOME/LightwaveExplorer/lightwave-explorer-research/projects" -type f \
  \( -iname '*STACK*NO*RESEED*' -o -iname '*REVERSE*STACK*' \
     -o -iname '*24*20*20*24*' \) -print
```

## Rerun both orders

Download this ZIP and unzip it as a *new* project. If your repository is in
the default location and the downloaded ZIP is in `~/Downloads`, run:

```zsh
project="$HOME/LightwaveExplorer/lightwave-explorer-research/projects/2026-09_PPLN_stacked_orders_single_seed"
mkdir -p "$project"
unzip "$HOME/Downloads/PPLN_stacked_order_rerun.zip" -d "$project"
cd "$project"
caffeinate -dimsu zsh run_both.sh
```

The script regenerates and verifies the reverse input, then runs the two LWE
files sequentially using the `lwe one` wrapper from your earlier projects.
If the wrapper is absent, it tries
`/Applications/LightwaveExplorer.app/Contents/MacOS/LightwaveExplorer`.
If the app is installed elsewhere, edit the `binary=` line in `run_both.sh`.
The script stops on a failed run and will not package incomplete results.

## Upload the results to this chat

After both runs succeed, this **one small file** appears in the project root:

```text
PPLN_single_seed_stacked_results.zip
```

Upload that ZIP here using the attachment button. It contains the *two input
files*, the *two raw `*_spectrum.dat` files* and a manifest confirming one
seed injection per run. The much larger LWE output ZIPs are not needed for
plotting; keep them on your Mac for provenance. If packaging fails, upload
both `outputs/*.zip` files and both `inputs/*.txt` files instead.

`package_results.py` also plots the available five spectra locally if NumPy
and Matplotlib are installed. All plotted curves are **solid**, with distinct
colours and thicker lines for the 1 mm stacks. Once the results are uploaded,
the completed figure can be checked for overlaps and finalised for the report.
