#!/usr/bin/env python3
"""Generate the missing 24→20→20→24 µm, single-entrance-seed LWE run.

The geometry is taken directly from the exact-length single-crystal inputs
used for Figure 12(b). It is the order reversal of the archived forward stack.
"""

from pathlib import Path


ROOT = Path(__file__).resolve().parent
REF = ROOT / "reference_inputs"
OUT = ROOT / "inputs"
FWD = "SEEDTEST_FWD_20_24_500_SINGLE"
REV = "SEEDTEST_REV_24_20_500_SINGLE"
STACK = "STACK_NO_RESEED"
TARGET = "REVERSE_STACK_SINGLE_SEED"


def field(text, label):
    prefix = label + ": "
    hits = [line[len(prefix):] for line in text.splitlines()
            if line.startswith(prefix)]
    if len(hits) != 1:
        raise ValueError(f"Expected one {label!r}, found {len(hits)}")
    return hits[0]


def change(text, label, value):
    prefix = label + ": "
    lines = text.splitlines(keepends=True)
    hits = [i for i, line in enumerate(lines) if line.startswith(prefix)]
    if len(hits) != 1:
        raise ValueError(f"Expected one {label!r}, found {len(hits)}")
    i = hits[0]
    lines[i] = prefix + value + ("\n" if lines[i].endswith("\n") else "")
    return "".join(lines)


fwd = (REF / f"{FWD}.txt").read_text()
rev = (REF / f"{REV}.txt").read_text()
original = (REF / f"{STACK}.txt").read_text()

seq_fwd = field(fwd, "Sequence")
seq_rev = field(rev, "Sequence")
seq_stack = field(original, "Sequence")
assert seq_fwd.startswith("init()") and seq_rev.startswith("init()")
assert seq_stack == seq_fwd + seq_rev[6:], (
    "Archived forward stack no longer matches its two source crystals"
)
assert all("addPulse(" not in s for s in (seq_fwd, seq_rev, seq_stack))
assert seq_fwd.count("nonlinear(") == seq_rev.count("nonlinear(") == 46
assert seq_stack.count("nonlinear(") == 92
assert float(field(original, "Thickness (m)")) == 0.001
assert float(field(original, "Pulse energy 1 (J)")) == 1e-7
assert float(field(original, "Pulse energy 2 (J)")) == 1e-8

seq_new = seq_rev + seq_fwd[6:]
assert seq_new.count("init()") == 1 and "addPulse(" not in seq_new
assert seq_new.count("nonlinear(") == 92

updated = change(original, "Sequence", seq_new)
updated = change(updated, "Output base path", f"outputs/{TARGET}")
OUT.mkdir(exist_ok=True)
target = OUT / f"{TARGET}.txt"
target.write_text(updated)
print(f"Wrote {target}")
print("1 mm total; 24→20→20→24 µm; one 10 nJ MIR entrance seed")
