#!/usr/bin/env python3
"""Integrate LWE field-3 spectral energy over the two target half-bands."""

from array import array
import csv
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import sys

ROOT = Path(__file__).resolve().parent
STEM = "SEEDTEST_FWD_22_24_250_SINGLE"
C = 299792458.0


def field(text, key):
    matches = [line.split(": ", 1)[1] for line in text.splitlines()
               if line.startswith(key + ": ")]
    if len(matches) != 1:
        raise ValueError("Expected one " + key)
    return matches[0]


input_path = ROOT / "inputs" / (STEM + ".txt")
input_bytes = input_path.read_bytes()
input_text = input_bytes.decode("utf-8")
dt = float(field(input_text, "dt (s)"))
span = float(field(input_text, "Time span (s)"))
n_time = int(8 * round(span / (8 * dt)))
n_freq = n_time // 2 + 1
zip_path = ROOT / "outputs" / (STEM + ".zip")
if not zip_path.exists():
    sys.exit("Missing LWE output: " + str(zip_path))
with ZipFile(zip_path) as z:
    matches = [n for n in z.namelist()
               if Path(n).name == STEM + "_spectrum.dat"]
    if len(matches) != 1:
        sys.exit("Expected one field-spectrum file in " + str(zip_path))
    spectrum_bytes = z.read(matches[0])

raw = array("d")
raw.frombytes(spectrum_bytes)
if len(raw) != 3 * n_freq:
    sys.exit("Unexpected spectrum length: {} vs expected {}".format(
        len(raw), 3 * n_freq))
field3 = raw[2 * n_freq:]
frequencies = [i / (n_time * dt) for i in range(n_freq)]


def band_energy_pj(low_nm, high_nm):
    indices = [i for i, f in enumerate(frequencies)
               if C / (high_nm * 1e-9) <= f <= C / (low_nm * 1e-9)]
    return sum((field3[i] + field3[j]) * (frequencies[j] - frequencies[i]) / 2
               for i, j in zip(indices, indices[1:])) * 1e12


rows = [
    (950, 1000, band_energy_pj(950, 1000)),
    (1000, 1050, band_energy_pj(1000, 1050)),
    (950, 1050, band_energy_pj(950, 1050)),
]
for low, high, value in rows:
    print("{}–{} nm: {:.3f} pJ".format(low, high, value))
print("The two half-band trapezoids need not sum exactly to the direct full-band integral;")
print("the discrete frequency samples leave a small interval around 1000 nm between them.")

result_dir = ROOT / "results"
result_dir.mkdir(exist_ok=True)
csv_path = result_dir / (STEM + "_band_energies.csv")
with csv_path.open("w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["case", "pump_GDD_fs2", "lower_nm", "upper_nm", "energy_pJ"])
    for low, high, value in rows:
        writer.writerow([STEM, 0, low, high, "{:.12g}".format(value)])

package = ROOT / (STEM + "_results.zip")
with ZipFile(package, "w", ZIP_DEFLATED) as z:
    z.writestr("inputs/" + input_path.name, input_bytes)
    z.writestr("spectra/" + STEM + "_spectrum.dat", spectrum_bytes)
    z.write(csv_path, "results/" + csv_path.name)
print("Ready to upload: " + str(package))
