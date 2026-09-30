#!/bin/zsh
set -euo pipefail
cd "${0:A:h}"

# Regenerate and verify the reverse sequence from the archived exact-length
# single-crystal controls. Both runs have one init() and no addPulse().
python3 generate_reverse_input.py

runs=(STACK_NO_RESEED REVERSE_STACK_SINGLE_SEED)
if command -v lwe >/dev/null 2>&1; then
  for run in "${runs[@]}"; do
    echo "Running single-seed stack: $run"
    lwe one "inputs/${run}.txt"
  done
else
  binary="/Applications/LightwaveExplorer.app/Contents/MacOS/LightwaveExplorer"
  if [[ ! -x "$binary" ]]; then
    echo "LWE not found. Check the 'lwe' command or installed app path in README.md." >&2
    exit 1
  fi
  for run in "${runs[@]}"; do
    echo "Running single-seed stack: $run"
    "$binary" "inputs/${run}.txt"
  done
fi

python3 package_results.py
echo "Ready to upload: PPLN_single_seed_stacked_results.zip"
