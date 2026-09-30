#!/usr/bin/env python3
"""Extract both completed LWE spectra and make one small upload archive."""

from hashlib import sha256
import json
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED


ROOT = Path(__file__).resolve().parent
RUNS = ("STACK_NO_RESEED", "REVERSE_STACK_SINGLE_SEED")
manifest = {
    "description": "Two exact 500+500 um PPLN stacks, both with a single coherent MIR seed at the entrance",
    "runs": [],
}
payloads = []

for stem in RUNS:
    input_path = ROOT / "inputs" / f"{stem}.txt"
    inp = input_path.read_bytes()
    lines = inp.decode("utf-8").splitlines()
    seq = next(line.split(": ", 1)[1] for line in lines if line.startswith("Sequence: "))
    if seq.count("init()") != 1 or seq.count("addPulse(") != 0 or seq.count("nonlinear(") != 92:
        raise ValueError(f"{stem}: input is not the intended one-entrance-seed, 92-domain stack")
    zip_path = ROOT / "outputs" / f"{stem}.zip"
    if not zip_path.is_file():
        raise FileNotFoundError(f"Missing completed LWE output {zip_path}")
    with ZipFile(zip_path) as z:
        matches = [name for name in z.namelist()
                   if name.split("/")[-1] == f"{stem}_spectrum.dat"]
        if len(matches) != 1:
            raise ValueError(f"{zip_path}: expected one {stem}_spectrum.dat, found {len(matches)}")
        spec = z.read(matches[0])
    time_span = float(next(line.split(": ", 1)[1] for line in lines
                           if line.startswith("Time span (s): ")))
    dt = float(next(line.split(": ", 1)[1] for line in lines
                    if line.startswith("dt (s): ")))
    n_time = int(8 * round(time_span / (8 * dt)))
    expected_bytes = 3 * (n_time // 2 + 1) * 8
    if len(spec) != expected_bytes:
        raise ValueError(f"{stem}: {len(spec)} spectrum bytes; expected {expected_bytes}")

    (ROOT / "data" / f"{stem}.txt").write_bytes(inp)
    (ROOT / "data" / f"{stem}_spectrum.dat").write_bytes(spec)
    payloads.extend(((f"inputs/{stem}.txt", inp),
                     (f"spectra/{stem}_spectrum.dat", spec)))
    manifest["runs"].append({
        "name": stem,
        "init_count": 1,
        "midpoint_addPulse_count": 0,
        "nonlinear_domain_count": 92,
        "spectrum_bytes": len(spec),
        "spectrum_sha256": sha256(spec).hexdigest(),
    })

out = ROOT / "PPLN_single_seed_stacked_results.zip"
with ZipFile(out, "w", ZIP_DEFLATED, compresslevel=6) as z:
    z.writestr("manifest.json", json.dumps(manifest, indent=2) + "\n")
    for name, data in payloads:
        z.writestr(name, data)

print(f"Saved {out} ({out.stat().st_size:,} bytes)")
try:
    import numpy  # noqa: F401
    import matplotlib  # noqa: F401
except ImportError:
    print("Optional local plot skipped: numpy and matplotlib are required.")
else:
    import subprocess
    import sys

    subprocess.run([sys.executable, str(ROOT / "plot_figure12b_stacks.py")],
                   check=True, cwd=ROOT)
