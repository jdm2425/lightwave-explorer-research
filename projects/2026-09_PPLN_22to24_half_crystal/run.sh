#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"

python3 make_input.py
mkdir -p outputs
input="inputs/SEEDTEST_FWD_22_24_250_SINGLE.txt"

if command -v lwe >/dev/null 2>&1; then
    lwe one "$input"
else
    binary="/Applications/LightwaveExplorer.app/Contents/MacOS/LightwaveExplorer"
    if [[ ! -x "$binary" ]]; then
        echo "Lightwave Explorer was not found at $binary or as the 'lwe' command." >&2
        exit 1
    fi
    "$binary" "$input"
fi

python3 analyze.py
