#!/usr/bin/env python3
from pathlib import Path
import csv
import zipfile
import numpy as np
import matplotlib.pyplot as plt

from publication_style import apply_publication_style, TWO_PANEL_FIGSIZE

apply_publication_style()

ROOT = Path(__file__).resolve().parents[1]
INPUT_RUN = "LINEAR_INPUT"
OUTPUT_RUN = "LINEAR_500um"


def read_value(path, label):
    prefix = label + ":"
    for line in Path(path).read_text().splitlines():
        if line.startswith(prefix):
            return float(line.split(":", 1)[1].strip())
    raise RuntimeError(f"Missing {label!r} in {path}")


def load_ext(run):
    inp = ROOT / "inputs" / f"{run}.txt"
    out = ROOT / "outputs" / f"{run}.zip"

    dt = read_value(inp, "dt (s)")
    T = read_value(inp, "Time span (s)")
    dx = read_value(inp, "dx (m)")
    width = read_value(inp, "Grid width (m)")

    nt = int(8 * round(T / (8 * dt)))
    nx = int(8 * round(width / (8 * dx)))

    with zipfile.ZipFile(out) as z:
        member = f"{run}_Ext.dat"
        if member not in z.namelist():
            matches = [x for x in z.namelist() if x.endswith("_Ext.dat")]
            if len(matches) != 1:
                raise RuntimeError(f"{run}: could not uniquely find Ext.dat")
            member = matches[0]
        raw = np.frombuffer(z.read(member), dtype=np.float64)

    expected = nt * nx * 2
    if raw.size < expected:
        raise RuntimeError(
            f"{run}: Ext.dat contains {raw.size} doubles; expected >= {expected}"
        )

    arr = raw[:expected].reshape((nt, nx, 2), order="F")
    ex = arr[:, :, 0]

    # Cylindrical grid: average the two points closest to r=0.
    j0 = nx // 2
    trace = 0.5 * (ex[:, j0 - 1] + ex[:, j0])

    return dt, trace


def analytic_band_envelope(trace, dt, f_lo, f_hi):
    n = trace.size
    f = np.fft.fftfreq(n, dt)
    E = np.fft.fft(trace)

    mask = (f > f_lo) & (f < f_hi)

    # Positive-frequency analytic representation.
    analytic = np.fft.ifft(2.0 * E * mask)
    intensity = np.abs(analytic) ** 2
    return intensity


def fwhm(t, y):
    y = np.asarray(y, float)
    y = y / y.max()
    i0 = int(np.argmax(y))

    # Search locally outwards from the main peak.
    left = i0
    while left > 0 and y[left] >= 0.5:
        left -= 1
    right = i0
    while right < len(y) - 1 and y[right] >= 0.5:
        right += 1

    def crossing(i, j):
        x1, x2 = t[i], t[j]
        y1, y2 = y[i], y[j]
        if y2 == y1:
            return 0.5 * (x1 + x2)
        return x1 + (0.5 - y1) * (x2 - x1) / (y2 - y1)

    tl = crossing(left, left + 1)
    tr = crossing(right - 1, right)
    return tr - tl


dt_in, Ein = load_ext(INPUT_RUN)
dt_out, Eout = load_ext(OUTPUT_RUN)
if abs(dt_in - dt_out) > 1e-30:
    raise RuntimeError("Input/output dt mismatch")

dt = dt_in
n = Ein.size

f_seed = read_value(ROOT / "inputs" / f"{INPUT_RUN}.txt", "Frequency 2 (Hz)")
f_pump = read_value(ROOT / "inputs" / f"{INPUT_RUN}.txt", "Frequency 1 (Hz)")
split = 0.5 * (f_seed + f_pump)

# Very broad, deliberately non-overlapping passbands.
seed_lo, seed_hi = 5e12, split
pump_lo, pump_hi = split, 0.5 / dt * 0.98

profiles = {}
for state, E in [("input", Ein), ("output", Eout)]:
    profiles[(state, "pump")] = analytic_band_envelope(E, dt, pump_lo, pump_hi)
    profiles[(state, "seed")] = analytic_band_envelope(E, dt, seed_lo, seed_hi)

# Set t=0 at the input pump peak.
pump_in_peak_i = int(np.argmax(profiles[("input", "pump")]))
t = (np.arange(n) - pump_in_peak_i) * dt

rows = []
for state in ("input", "output"):
    for pulse in ("pump", "seed"):
        y = profiles[(state, pulse)]
        peak_t = t[int(np.argmax(y))]
        width = fwhm(t, y)
        rows.append({
            "state": state,
            "pulse": pulse,
            "peak_time_fs": peak_t * 1e15,
            "intensity_FWHM_fs": width * 1e15,
        })

def get(state, pulse, key):
    return next(r[key] for r in rows if r["state"] == state and r["pulse"] == pulse)

rel_in = get("input", "seed", "peak_time_fs") - get("input", "pump", "peak_time_fs")
rel_out = get("output", "seed", "peak_time_fs") - get("output", "pump", "peak_time_fs")
walkoff = rel_out - rel_in

summary = [
    ["input_pump_FWHM_fs", get("input", "pump", "intensity_FWHM_fs")],
    ["output_pump_FWHM_fs", get("output", "pump", "intensity_FWHM_fs")],
    ["input_seed_FWHM_fs", get("input", "seed", "intensity_FWHM_fs")],
    ["output_seed_FWHM_fs", get("output", "seed", "intensity_FWHM_fs")],
    ["input_seed_minus_pump_peak_delay_fs", rel_in],
    ["output_seed_minus_pump_peak_delay_fs", rel_out],
    ["measured_linear_walkoff_fs", walkoff],
]

processed = ROOT / "processed"
figures = ROOT / "figures"
processed.mkdir(exist_ok=True)
figures.mkdir(exist_ok=True)

with (processed / "linear_dispersion_summary.csv").open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["quantity", "value"])
    w.writerows(summary)

with (processed / "linear_dispersion_pulse_metrics.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys())
    w.writeheader()
    w.writerows(rows)

fig, axes = plt.subplots(1, 2, figsize=TWO_PANEL_FIGSIZE)

for ax, pulse, title in [
    (axes[0], "pump", "800 nm pump"),
    (axes[1], "seed", "4 µm MIR seed"),
]:
    yin = profiles[("input", pulse)]
    yout = profiles[("output", pulse)]
    yin = yin / yin.max()
    yout = yout / yout.max()

    # Wide enough for the strongly dispersed seed.
    m = (t * 1e15 >= -700) & (t * 1e15 <= 700)
    ax.plot(t[m] * 1e15, yin[m], label="Input")
    ax.plot(t[m] * 1e15, yout[m], label="After 500 µm")
    ax.set_title(title)
    ax.set_xlabel("Time (fs)")
    ax.grid(True)

axes[0].set_ylabel("Normalised intensity")
axes[1].legend()
fig.tight_layout()
fig.savefig(figures / "linear_dispersion_500um.pdf")
fig.savefig(figures / "linear_dispersion_500um.png")
plt.close(fig)

print("Linear-propagation result")
print("=========================")
for k, v in summary:
    print(f"{k:42s} {v:10.3f}")

print()
print("Saved:")
print("  processed/linear_dispersion_summary.csv")
print("  processed/linear_dispersion_pulse_metrics.csv")
print("  figures/linear_dispersion_500um.pdf")
print("  figures/linear_dispersion_500um.png")
