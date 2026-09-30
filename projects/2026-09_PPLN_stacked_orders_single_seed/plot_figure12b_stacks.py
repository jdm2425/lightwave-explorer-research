"""Compare Figure 12(b) spectra with single-seeded stacked PPLN runs.

All curves use the LWE field-3 spectral energy density, converted from J/Hz
to pJ/nm. The optional reverse stack appears when its matching input and
spectrum files are placed in data/ as REVERSE_STACK_SINGLE_SEED.*.
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
C = 299_792_458.0

mpl.rcParams.update({
    "figure.figsize": (3.4, 2.6),
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "font.family": "sans-serif",
    "font.size": 7,
    "axes.labelsize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 5.6,
    "lines.linewidth": 1.5,
    "axes.linewidth": 0.8,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "grid.linestyle": ":",
    "grid.linewidth": 0.5,
    "grid.alpha": 0.6,
    "legend.frameon": False,
    "pdf.fonttype": 42,
})


def input_value(path, label):
    prefix = label + ":"
    for line in path.read_text().splitlines():
        if line.startswith(prefix):
            return float(line.partition(":")[2].strip())
    raise ValueError(f"Missing {label} in {path}")


def load_spectrum(stem):
    inp = DATA / (stem + ".txt")
    spec = DATA / (stem + "_spectrum.dat")
    time_span = input_value(inp, "Time span (s)")
    dt = input_value(inp, "dt (s)")
    n_time = int(8 * round(time_span / (8 * dt)))
    f = np.fft.rfftfreq(n_time, dt)
    raw = np.fromfile(spec, dtype=np.float64)
    if raw.size != 3 * f.size:
        raise ValueError(f"{spec}: {raw.size} values, expected {3 * f.size}")
    sf = raw.reshape((f.size, 3), order="F")[:, 2]

    mask = f > 0
    wavelength_m = C / f[mask]
    wavelength_nm = wavelength_m * 1e9
    # Convert J/Hz into pJ/nm.
    density = sf[mask] * C / wavelength_m**2 * 1e3
    order = np.argsort(wavelength_nm)

    band = (f >= C / (1050e-9)) & (f <= C / (950e-9))
    energy_pj = np.trapz(sf[band], f[band]) * 1e12
    return wavelength_nm[order], density[order], energy_pj


cases = [
    ("Fixed 21.75 µm · 0.5 mm", "PPLN_FIXED_21p75um", "#4D4D4D", 1.3),
    ("20→24 µm · 0.5 mm", "SEEDTEST_FWD_20_24_500_SINGLE", "#009E73", 1.3),
    ("24→20 µm · 0.5 mm", "SEEDTEST_REV_24_20_500_SINGLE", "#0072B2", 1.3),
    ("20→24→24→20 µm · 1 mm", "STACK_NO_RESEED", "#D55E00", 1.9),
]

reverse_stem = "REVERSE_STACK_SINGLE_SEED"
if all((DATA / (reverse_stem + suffix)).is_file()
       for suffix in (".txt", "_spectrum.dat")):
    cases.append(("24→20→20→24 µm · 1 mm", reverse_stem, "#CC79A7", 1.9))

fig, ax = plt.subplots()
ax.axvspan(950, 1050, color="0.75", alpha=0.18, zorder=0)
for label, stem, color, linewidth in cases:
    lam, density, energy = load_spectrum(stem)
    shown = (lam >= 900) & (lam <= 1125)
    ax.plot(lam[shown], density[shown], label=label,
            color=color, linestyle="-", linewidth=linewidth)
    print(f"{label}: {energy:.1f} pJ in 950–1050 nm")

ax.set_xlim(900, 1125)
ax.set_ylim(bottom=0)
ax.set_xlabel("Wavelength (nm)")
ax.set_ylabel("Spectral energy density (pJ/nm)")
ax.grid(True)
ax.legend(loc="upper right", ncol=1, handlelength=2.5)
fig.tight_layout()
for suffix in ("png", "pdf"):
    fig.savefig(ROOT / f"Figure12b_with_single_seed_stack.{suffix}",
                bbox_inches="tight", pad_inches=0.03)
plt.close(fig)
