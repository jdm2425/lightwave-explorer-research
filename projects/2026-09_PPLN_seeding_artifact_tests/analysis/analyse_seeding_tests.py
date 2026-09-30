#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import csv
import zipfile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from publication_style import apply_publication_style, SINGLE_COLUMN, TWO_PANEL_FIGSIZE

apply_publication_style()

C = 299_792_458.0
ROOT = Path(__file__).resolve().parents[1]
FIXED = ROOT.parent / "2026-09_PPLN_fixed_period_GDD0"
LEGACY = ROOT.parent / "2026-09_PPLN_figure13_updated"


def read_value(path: Path, label: str) -> float:
    prefix = label + ":"
    for line in path.read_text().splitlines():
        if line.startswith(prefix):
            return float(line.split(":", 1)[1].strip())
    raise RuntimeError(f"Could not find {label!r} in {path}")


def load(run: str, inp: Path, out: Path):
    T = read_value(inp, "Time span (s)")
    dt = read_value(inp, "dt (s)")
    nt = int(8 * round(T / (8 * dt)))
    nf = nt // 2 + 1
    f = np.fft.rfftfreq(nt, dt)

    with zipfile.ZipFile(out) as z:
        exact = f"{run}_spectrum.dat"
        if exact in z.namelist():
            member = exact
        else:
            matches = [x for x in z.namelist() if x.endswith("_spectrum.dat")]
            if len(matches) != 1:
                raise RuntimeError(f"{run}: spectrum not unique: {matches}")
            member = matches[0]
        raw = np.frombuffer(z.read(member), dtype=np.float64)

    if raw.size != 3 * nf:
        raise RuntimeError(f"{run}: {raw.size} doubles; expected {3*nf}")

    sf = raw.reshape((nf, 3), order="F")[:, 2]
    return f, sf


def load_project_run(run):
    return load(
        run,
        ROOT / "inputs" / f"{run}.txt",
        ROOT / "outputs" / f"{run}.zip",
    )


def wavelength_density_pj_nm(f, sf):
    m = f > 0
    f = f[m]
    sf = sf[m]
    lm = C / f
    lam_nm = lm * 1e9
    # J/Hz -> J/m -> pJ/nm gives an overall factor 1e3.
    sl = sf * C / lm**2 * 1e3
    order = np.argsort(lam_nm)
    return lam_nm[order], sl[order]


def integrate_band_pj(f, sf, lo_nm, hi_nm):
    flo = C / (hi_nm * 1e-9)
    fhi = C / (lo_nm * 1e-9)
    m = (f >= flo) & (f <= fhi)
    return np.trapezoid(sf[m], f[m]) * 1e12


def metrics(run, f, sf):
    lam, sl = wavelength_density_pj_nm(f, sf)
    target = (lam >= 950) & (lam <= 1050)
    x = lam[target]
    y = sl[target]

    total = integrate_band_pj(f, sf, 950, 1050)
    blue = integrate_band_pj(f, sf, 950, 1000)
    red = integrate_band_pj(f, sf, 1000, 1050)
    peak_i = int(np.argmax(y))
    centroid = np.trapezoid(x * y, x) / np.trapezoid(y, x)

    return {
        "run_name": run,
        "E_950_1050_pJ": total,
        "E_950_1000_pJ": blue,
        "E_1000_1050_pJ": red,
        "red_blue_ratio": red / blue if blue else np.nan,
        "peak_wavelength_950_1050_nm": x[peak_i],
        "spectral_centroid_950_1050_nm": centroid,
        "peak_density_pJ_per_nm": y[peak_i],
        "equivalent_width_nm": total / y[peak_i] if y[peak_i] else np.nan,
    }


runs = [
    "SEEDTEST_FWD_20_24_500_SINGLE",
    "SEEDTEST_REV_24_20_500_SINGLE",
    "SEEDTEST_FWD_20_22_250_SINGLE",
    "SEEDTEST_REV_24_22_250_SINGLE",
    "SEEDTEST_FWD_20_24_500_DISTRIBUTED",
    "SEEDTEST_REV_24_20_500_DISTRIBUTED",
]

data = {}
rows = []
for run in runs:
    f, sf = load_project_run(run)
    data[run] = (f, sf)
    rows.append(metrics(run, f, sf))

processed = ROOT / "processed"
figures = ROOT / "figures"
processed.mkdir(exist_ok=True)
figures.mkdir(exist_ok=True)

df = pd.DataFrame(rows)
df.to_csv(processed / "seeding_artifact_metrics.csv", index=False)


def plot_curve(ax, run, label, **kwargs):
    f, sf = data[run]
    lam, sl = wavelength_density_pj_nm(f, sf)
    m = (lam >= 900) & (lam <= 1125)
    ax.plot(lam[m], sl[m], label=label, **kwargs)


# 1. Exact-length direction comparison including fixed-period reference.
fig, ax = plt.subplots(figsize=SINGLE_COLUMN)

fixed_run = "PPLN_FIXED_21p75um"
ff, fs = load(
    fixed_run,
    FIXED / "inputs" / f"{fixed_run}.txt",
    FIXED / "outputs" / f"{fixed_run}.zip",
)
lam, sl = wavelength_density_pj_nm(ff, fs)
m = (lam >= 900) & (lam <= 1125)
ax.plot(lam[m], sl[m], label="Fixed 21.75 µm")

plot_curve(ax, "SEEDTEST_FWD_20_24_500_SINGLE", "Chirped 20→24 µm")
plot_curve(ax, "SEEDTEST_REV_24_20_500_SINGLE", "Chirped 24→20 µm")

ax.axvspan(950, 1050, alpha=0.12)
ax.axvline(1000, linestyle=":", linewidth=0.8)
ax.set_xlim(900, 1125)
ax.set_xlabel("Wavelength (nm)")
ax.set_ylabel("Spectral energy density (pJ/nm)")
ax.grid(True)
ax.legend()
fig.tight_layout()
fig.savefig(figures / "exact_chirp_direction_comparison.pdf")
fig.savefig(figures / "exact_chirp_direction_comparison.png")
plt.close(fig)


# 2. First-half crystal diagnostic.
fig, axes = plt.subplots(1, 2, figsize=TWO_PANEL_FIGSIZE, sharex=True)

plot_curve(
    axes[0], "SEEDTEST_FWD_20_24_500_SINGLE",
    "20→24 µm, 500 µm"
)
plot_curve(
    axes[0], "SEEDTEST_FWD_20_22_250_SINGLE",
    "20→22 µm, first 250 µm"
)
axes[0].set_title("Increasing chirp")

plot_curve(
    axes[1], "SEEDTEST_REV_24_20_500_SINGLE",
    "24→20 µm, 500 µm"
)
plot_curve(
    axes[1], "SEEDTEST_REV_24_22_250_SINGLE",
    "24→22 µm, first 250 µm"
)
axes[1].set_title("Decreasing chirp")

for ax in axes:
    ax.axvspan(950, 1050, alpha=0.12)
    ax.axvline(1000, linestyle=":", linewidth=0.8)
    ax.set_xlim(900, 1125)
    ax.set_xlabel("Wavelength (nm)")
    ax.grid(True)
    ax.legend()

axes[0].set_ylabel("Spectral energy density (pJ/nm)")
fig.tight_layout()
fig.savefig(figures / "half_crystal_test.pdf")
fig.savefig(figures / "half_crystal_test.png")
plt.close(fig)


# 3. Distributed-seed diagnostic.
fig, axes = plt.subplots(1, 2, figsize=TWO_PANEL_FIGSIZE, sharex=True)

plot_curve(
    axes[0], "SEEDTEST_FWD_20_24_500_SINGLE",
    "Single entrance seed"
)
plot_curve(
    axes[0], "SEEDTEST_FWD_20_24_500_DISTRIBUTED",
    "Seed distributed by period"
)
axes[0].set_title("20→24 µm")

plot_curve(
    axes[1], "SEEDTEST_REV_24_20_500_SINGLE",
    "Single entrance seed"
)
plot_curve(
    axes[1], "SEEDTEST_REV_24_20_500_DISTRIBUTED",
    "Seed distributed by period"
)
axes[1].set_title("24→20 µm")

for ax in axes:
    ax.axvspan(950, 1050, alpha=0.12)
    ax.axvline(1000, linestyle=":", linewidth=0.8)
    ax.set_xlim(900, 1125)
    ax.set_xlabel("Wavelength (nm)")
    ax.grid(True)
    ax.legend()

axes[0].set_ylabel("Spectral energy density (pJ/nm)")
fig.tight_layout()
fig.savefig(figures / "distributed_seed_test.pdf")
fig.savefig(figures / "distributed_seed_test.png")
plt.close(fig)


def row(run):
    return df.loc[df.run_name == run].iloc[0]


comparisons = []

for direction, full, half in [
    ("forward", "SEEDTEST_FWD_20_24_500_SINGLE", "SEEDTEST_FWD_20_22_250_SINGLE"),
    ("reverse", "SEEDTEST_REV_24_20_500_SINGLE", "SEEDTEST_REV_24_22_250_SINGLE"),
]:
    a, b = row(full), row(half)
    comparisons.append({
        "test": "half_vs_full",
        "direction": direction,
        "comparison": "half/full",
        "target_energy_ratio": b.E_950_1050_pJ / a.E_950_1050_pJ,
        "blue_energy_ratio": b.E_950_1000_pJ / a.E_950_1000_pJ,
        "red_energy_ratio": b.E_1000_1050_pJ / a.E_1000_1050_pJ,
        "full_red_blue_ratio": a.red_blue_ratio,
        "test_red_blue_ratio": b.red_blue_ratio,
    })

for direction, single, dist in [
    ("forward", "SEEDTEST_FWD_20_24_500_SINGLE", "SEEDTEST_FWD_20_24_500_DISTRIBUTED"),
    ("reverse", "SEEDTEST_REV_24_20_500_SINGLE", "SEEDTEST_REV_24_20_500_DISTRIBUTED"),
]:
    a, b = row(single), row(dist)
    comparisons.append({
        "test": "distributed_vs_single",
        "direction": direction,
        "comparison": "distributed/single",
        "target_energy_ratio": b.E_950_1050_pJ / a.E_950_1050_pJ,
        "blue_energy_ratio": b.E_950_1000_pJ / a.E_950_1000_pJ,
        "red_energy_ratio": b.E_1000_1050_pJ / a.E_1000_1050_pJ,
        "full_red_blue_ratio": a.red_blue_ratio,
        "test_red_blue_ratio": b.red_blue_ratio,
    })

cdf = pd.DataFrame(comparisons)
cdf.to_csv(processed / "seeding_artifact_comparisons.csv", index=False)


# Optional comparison to the legacy loop-based chirps if those outputs are present.
legacy_rows = []
for direction, legacy_run, exact_run in [
    ("forward", "FIG13_CHIRP_20_24_GDD0", "SEEDTEST_FWD_20_24_500_SINGLE"),
    ("reverse", "FIG13_CHIRP_24_20_GDD0", "SEEDTEST_REV_24_20_500_SINGLE"),
]:
    linp = LEGACY / "inputs" / f"{legacy_run}.txt"
    lout = LEGACY / "outputs" / f"{legacy_run}.zip"
    if linp.exists() and lout.exists():
        lf, ls = load(legacy_run, linp, lout)
        lm = metrics(legacy_run, lf, ls)
        em = row(exact_run)
        legacy_rows.append({
            "direction": direction,
            "legacy_run": legacy_run,
            "exact_run": exact_run,
            "legacy_E_950_1050_pJ": lm["E_950_1050_pJ"],
            "exact_E_950_1050_pJ": em.E_950_1050_pJ,
            "exact_over_legacy_energy": em.E_950_1050_pJ / lm["E_950_1050_pJ"],
            "legacy_red_blue_ratio": lm["red_blue_ratio"],
            "exact_red_blue_ratio": em.red_blue_ratio,
        })

if legacy_rows:
    pd.DataFrame(legacy_rows).to_csv(
        processed / "legacy_vs_exact_chirp_geometry.csv", index=False
    )


print()
print("Absolute metrics")
print("================")
print(df.to_string(index=False))

print()
print("Diagnostic comparisons")
print("======================")
print(cdf.to_string(index=False))

if legacy_rows:
    print()
    print("Legacy loop geometry vs exact 500 µm geometry")
    print("=============================================")
    print(pd.DataFrame(legacy_rows).to_string(index=False))

print()
print("Saved:")
print("  processed/seeding_artifact_metrics.csv")
print("  processed/seeding_artifact_comparisons.csv")
if legacy_rows:
    print("  processed/legacy_vs_exact_chirp_geometry.csv")
print("  figures/exact_chirp_direction_comparison.pdf")
print("  figures/half_crystal_test.pdf")
print("  figures/distributed_seed_test.pdf")
