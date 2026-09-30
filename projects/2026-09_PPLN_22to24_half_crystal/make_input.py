#!/usr/bin/env python3
"""Create the exact-length 22→24 µm PPLN half-crystal LWE input."""

from math import expm1, isclose
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent
STEM = "SEEDTEST_FWD_22_24_250_SINGLE"
SOURCE = ROOT / "reference_inputs" / "SEEDTEST_FWD_20_24_500_SINGLE.txt"


def get(text, key):
    matches = [line.split(": ", 1)[1] for line in text.splitlines()
               if line.startswith(key + ": ")]
    if len(matches) != 1:
        raise ValueError("Expected one field: " + key)
    return matches[0]


def set_field(text, key, value):
    pattern = r"^" + re.escape(key) + r": .*?$"
    result, count = re.subn(pattern, key + ": " + str(value), text, flags=re.M)
    if count != 1:
        raise ValueError("Expected one field to replace: " + key)
    return result


def domains(start_um, end_um, length_um):
    """Domain cuts at integer values of integral[0,z] 2/Lambda(z) dz.

    The same construction reproduces the archived 20→24 µm, 500 µm
    sequence, including the shortened final domain.
    """
    slope = (end_um - start_um) / length_um
    cuts = [0.0]
    n = 1
    while True:
        z = start_um / slope * expm1(n * slope / 2)
        if z >= length_um:
            break
        cuts.append(z)
        n += 1
    cuts.append(float(length_um))
    return [right - left for left, right in zip(cuts, cuts[1:])]


source = SOURCE.read_text()
assert float(get(source, "Pulse energy 1 (J)")) == 1e-7
assert float(get(source, "Pulse energy 2 (J)")) == 1e-8
assert float(get(source, "GDD 1 (s^-2)")) == 0.0
assert float(get(source, "GDD 2 (s^-2)")) == 0.0
assert float(get(source, "Thickness (m)")) == 0.0005

old = list(map(float, re.findall(r"nonlinear\(d,d,d,([\d.]+),d\)",
                                get(source, "Sequence"))))
reference = domains(20, 24, 500)
assert len(old) == len(reference) == 46
assert all(isclose(x, y, abs_tol=1e-8) for x, y in zip(old, reference))

segments = domains(22, 24, 250)
assert len(segments) == 22 and isclose(sum(segments), 250, abs_tol=1e-10)
sequence = "init()" + "".join(
    "nonlinear(d,d,d,{:.12g},d)rotate(180)".format(x) for x in segments
)
assert sequence.count("init()") == 1 and "addPulse(" not in sequence

result = set_field(source, "Thickness (m)", "0.00025")
result = set_field(result, "Sequence", sequence)
result = set_field(result, "Output base path", "outputs/" + STEM)
out = ROOT / "inputs" / (STEM + ".txt")
out.parent.mkdir(exist_ok=True)
out.write_text(result)
print("Generated {} ({} domains, {:.6f} µm total)".format(
    out, len(segments), sum(segments)))
