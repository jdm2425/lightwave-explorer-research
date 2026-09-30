#!/bin/zsh
set -euo pipefail

runs=(
  SEEDTEST_FWD_20_24_500_SINGLE
  SEEDTEST_REV_24_20_500_SINGLE
  SEEDTEST_FWD_20_22_250_SINGLE
  SEEDTEST_REV_24_22_250_SINGLE
  SEEDTEST_FWD_20_24_500_DISTRIBUTED
  SEEDTEST_REV_24_20_500_DISTRIBUTED
)

for run in "${runs[@]}"; do
  echo
  echo "=================================================="
  echo "Running $run"
  echo "=================================================="
  lwe one "inputs/${run}.txt"
done

echo
echo "=================================================="
echo "Analysing"
echo "=================================================="
~/LightwaveExplorer/.venv/bin/python analysis/analyse_seeding_tests.py

echo
echo "All seeding-artifact diagnostics completed."
