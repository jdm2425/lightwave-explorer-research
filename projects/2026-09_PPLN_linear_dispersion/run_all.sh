#!/bin/zsh
set -euo pipefail

echo "=== Exact input field ==="
lwe one inputs/LINEAR_INPUT.txt

echo
echo "=== 500 um linear MgO:LiNbO3 propagation ==="
lwe one inputs/LINEAR_500um.txt

echo
echo "=== Analysis ==="
~/LightwaveExplorer/.venv/bin/python analysis/analyse_linear_dispersion.py
