#!/bin/zsh
set -euo pipefail

echo
echo "=== Fine-grid chirped validation ==="
lwe one inputs/VAL_CHIRP_GDD0_FINE.txt

echo
echo "=== Seed energy = 1 nJ ==="
lwe one inputs/VAL_FIXED_SEED1nJ.txt
lwe one inputs/VAL_CHIRP_SEED1nJ.txt

echo
echo "=== Seed delay = -100 fs ==="
lwe one inputs/VAL_FIXED_DELAYm100fs.txt
lwe one inputs/VAL_CHIRP_DELAYm100fs.txt

echo
echo "=== Seed delay = +100 fs ==="
lwe one inputs/VAL_FIXED_DELAYp100fs.txt
lwe one inputs/VAL_CHIRP_DELAYp100fs.txt

echo
echo "All seven validation simulations completed."
