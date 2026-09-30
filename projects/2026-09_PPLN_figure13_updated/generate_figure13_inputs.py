#!/usr/bin/env python3
"""Generate updated Figure-13 chirped-PPLN runs from the current fixed-period input.

The existing fixed 21.75 um production input is used as the source of truth for
all optical, material and numerical parameters. Only the poling profile and run
name/output path are changed.

This makes the new chirped-vs-fixed comparison genuinely apples-to-apples.
"""
from __future__ import annotations

from pathlib import Path
import csv
import re
import sys

PROJECT = Path(__file__).resolve().parent
INPUTS = PROJECT / "inputs"
INPUTS.mkdir(exist_ok=True)

if len(sys.argv) != 2:
    raise SystemExit(
        "Usage: python generate_figure13_inputs.py "
        "/path/to/PPLN_FIXED_21p75um.txt"
    )

TEMPLATE_PATH = Path(sys.argv[1]).expanduser().resolve()
if not TEMPLATE_PATH.exists():
    raise SystemExit(f"Template input not found: {TEMPLATE_PATH}")

template = TEMPLATE_PATH.read_text()

DESIGNS = [
    ("FIG13_CHIRP_18_24",    18.0, 24.0, "18–24 µm"),
    ("FIG13_CHIRP_19p4_22",  19.4, 22.0, "19.4–22 µm"),
    ("FIG13_CHIRP_20_24",    20.0, 24.0, "20–24 µm"),
    ("FIG13_CHIRP_17p5_20p4",17.5, 20.4, "17.5–20.4 µm"),
]

def get_value(text: str, label: str) -> str:
    prefix = label + ":"
    for line in text.splitlines():
        if line.startswith(prefix):
            return line.split(":", 1)[1].strip()
    raise RuntimeError(f"Could not find {label!r} in template")

def replace_line(text: str, label: str, value: str) -> str:
    pattern = rf"(?m)^{re.escape(label)}:.*$"
    out, n = re.subn(pattern, f"{label}: {value}", text, count=1)
    if n != 1:
        raise RuntimeError(f"Could not uniquely replace {label!r}")
    return out

def chirped_sequence(start_um: float, end_um: float) -> str:
    # This reproduces the linearly chirped APPLN construction used in the
    # archived Lightwave Explorer runs: the local full period changes linearly,
    # and d_eff is reversed after each half-period.
    return (
        "init()"
        f"set(00,{start_um:g})"
        f"set(01,{end_um:g})"
        "set(02,2*i34/(v00+v01)-1)"
        "set(03,(v01-v00)/v02)"
        "for(v02,04){"
        "nonlinear(d,d,d,(v00+v04*v03)/2,d)rotate(180)"
        "nonlinear(d,d,d,(v00+v04*v03)/2,d)rotate(180)"
        "}"
    )

# Record the source configuration before modifying anything.
snapshot_labels = [
    "Pulse energy 1 (J)", "Pulse energy 2 (J)",
    "Frequency 1 (Hz)", "Frequency 2 (Hz)",
    "Bandwidth 1 (Hz)", "Bandwidth 2 (Hz)",
    "SG order 1", "SG order 2",
    "GDD 1 (s^-2)", "GDD 2 (s^-2)",
    "Beamwaist 1 (m)", "Beamwaist 2 (m)",
    "Grid width (m)", "dx (m)",
    "Time span (s)", "dt (s)", "dz (m)",
]
(PROJECT / "provenance").mkdir(exist_ok=True)
with (PROJECT / "provenance" / "template_parameters.csv").open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["parameter", "value"])
    for label in snapshot_labels:
        w.writerow([label, get_value(template, label)])
    w.writerow(["source_template", str(TEMPLATE_PATH)])

rows = []
for run_name, start_um, end_um, label in DESIGNS:
    text = template
    text = replace_line(text, "Thickness (m)", "0.0005")
    text = replace_line(text, "Batch mode", "0")
    text = replace_line(text, "Batch destination", "1e-05")
    text = replace_line(text, "Batch steps", "1")
    text = replace_line(text, "Batch mode 2", "0")
    text = replace_line(text, "Batch destination 2", "1000")
    text = replace_line(text, "Batch steps 2", "1")
    text = replace_line(text, "Sequence", chirped_sequence(start_um, end_um))
    text = replace_line(text, "Output base path", f"outputs/{run_name}")

    path = INPUTS / f"{run_name}.txt"
    path.write_text(text)

    rows.append({
        "run_name": run_name,
        "design_label": label,
        "poling_start_um": start_um,
        "poling_end_um": end_um,
        "crystal_length_um": 500.0,
        "input_file": f"inputs/{run_name}.txt",
    })

with (PROJECT / "run_manifest.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys())
    w.writeheader()
    w.writerows(rows)

print("Generated updated Figure-13 chirped-PPLN inputs:")
for r in rows:
    print(f"  {r['run_name']:24s} {r['design_label']}")
print()
print("All pump/seed/grid parameters were copied from:")
print(f"  {TEMPLATE_PATH}")
print()
print("Only the poling profile and output path were changed.")
