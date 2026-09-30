#!/usr/bin/env python3
"""Replot Figure 13 and compare the best fixed period against 20–24 um chirped PPLN."""
from __future__ import annotations

from pathlib import Path
import argparse
import zipfile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from publication_style import apply_publication_style, SINGLE_COLUMN, WIDE_FIGSIZE

apply_publication_style()

C = 299_792_458.0
TARGET_MIN_NM = 950.0
TARGET_MAX_NM = 1050.0

ROOT = Path(__file__).resolve().parents[1]


def read_value(path: Path, label: str) -> float:
    prefix = label + ":"
    for line in path.read_text().splitlines():
        if line.startswith(prefix):
            return float(line.split(":", 1)[1].strip())
    raise ValueError(f"Could not find {label!r} in {path}")


def load_spectrum(input_path: Path, zip_path: Path, run_name: str):
    T = read_value(input_path, "Time span (s)")
    dt = read_value(input_path, "dt (s)")
    nt = int(8 * round(T / (8 * dt)))
    nf = nt // 2 + 1
    f = np.fft.rfftfreq(nt, dt)

    if not zip_path.exists():
        raise FileNotFoundError(zip_path)

    with zipfile.ZipFile(zip_path) as z:
        exact = f"{run_name}_spectrum.dat"
        if exact in z.namelist():
            member = exact
        else:
            matches = [n for n in z.namelist() if n.endswith("_spectrum.dat")]
            if len(matches) != 1:
                raise RuntimeError(f"{run_name}: spectrum not unique: {matches}")
            member = matches[0]
        raw = np.frombuffer(z.read(member), dtype=np.float64)

    expected = 3 * nf
    if raw.size != expected:
        raise ValueError(
            f"{run_name}: spectrum has {raw.size} doubles; expected {expected}"
        )

    sf = raw.reshape((nf, 3), order="F")[:, 2]
    return f, sf


def wavelength_density(freq_hz, sf_j_per_hz):
    positive = freq_hz > 0
    f = freq_hz[positive]
    sf = sf_j_per_hz[positive]
    lam_m = C / f

    # J/Hz -> J/m -> nJ/nm. 1 m = 1e9 nm and 1 J = 1e9 nJ,
    # so the numerical conversion factors cancel.
    s_lambda_nj_per_nm = sf * C / lam_m**2

    lam_nm = lam_m * 1e9
    order = np.argsort(lam_nm)
    return lam_nm[order], s_lambda_nj_per_nm[order]


def band_energy_nj(freq_hz, sf):
    lo = C / (TARGET_MAX_NM * 1e-9)
    hi = C / (TARGET_MIN_NM * 1e-9)
    m = (freq_hz >= lo) & (freq_hz <= hi)
    return np.trapezoid(sf[m], freq_hz[m]) * 1e9


def shape_metrics(freq_hz, sf):
    lam, sl = wavelength_density(freq_hz, sf)
    m = (lam >= TARGET_MIN_NM) & (lam <= TARGET_MAX_NM)
    x = lam[m]
    y = sl[m]

    E = band_energy_nj(freq_hz, sf)
    peak = float(np.nanmax(y))
    eq_width = E / peak if peak > 0 else np.nan
    occupancy = eq_width / (TARGET_MAX_NM - TARGET_MIN_NM)

    mean = float(np.nanmean(y))
    cv = float(np.nanstd(y) / mean) if mean > 0 else np.nan

    half = 0.5 * peak
    above = x[y >= half]
    width50 = float(above.max() - above.min()) if len(above) >= 2 else 0.0

    return {
        "integrated_950_1050_nJ": E,
        "peak_density_nJ_per_nm": peak,
        "equivalent_width_nm": eq_width,
        "spectral_occupancy": occupancy,
        "coefficient_of_variation": cv,
        "width_above_50pct_peak_nm": width50,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fixed-project",
        type=Path,
        default=ROOT.parent / "2026-09_PPLN_fixed_period",
    )
    parser.add_argument(
        "--fixed-run",
        default="PPLN_FIXED_21p75um",
    )
    args = parser.parse_args()

    manifest = pd.read_csv(ROOT / "run_manifest.csv")
    spectra = {}
    metrics = []

    for _, row in manifest.iterrows():
        run = row["run_name"]
        inp = ROOT / row["input_file"]
        z = ROOT / "outputs" / f"{run}.zip"
        try:
            f, sf = load_spectrum(inp, z, run)
        except FileNotFoundError:
            print(f"Missing: {z}")
            continue

        spectra[run] = (row["design_label"], f, sf)
        rec = {"case": row["design_label"], "run_name": run, "geometry": "chirped"}
        rec.update(shape_metrics(f, sf))
        metrics.append(rec)

    fixed_project = args.fixed_project.resolve()
    fixed_run = args.fixed_run
    fixed_input = fixed_project / "inputs" / f"{fixed_run}.txt"
    fixed_zip = fixed_project / "outputs" / f"{fixed_run}.zip"

    if not fixed_input.exists() or not fixed_zip.exists():
        raise SystemExit(
            "Fixed reference not found. Expected:\\n"
            f"  {fixed_input}\\n"
            f"  {fixed_zip}"
        )

    ff, fs = load_spectrum(fixed_input, fixed_zip, fixed_run)
    fixed_label = "Fixed 21.75 µm"
    fixed_metrics = {"case": fixed_label, "run_name": fixed_run, "geometry": "fixed"}
    fixed_metrics.update(shape_metrics(ff, fs))
    metrics.append(fixed_metrics)

    processed = ROOT / "processed"
    figures = ROOT / "figures"
    processed.mkdir(exist_ok=True)
    figures.mkdir(exist_ok=True)

    mdf = pd.DataFrame(metrics)
    mdf.to_csv(processed / "figure13_metrics.csv", index=False)

    # ---- Full updated Figure 13: four chirped designs + fixed reference ----
    fig, ax = plt.subplots(figsize=WIDE_FIGSIZE)

    for run in manifest["run_name"]:
        if run not in spectra:
            continue
        label, f, sf = spectra[run]
        lam, sl = wavelength_density(f, sf)
        m = (lam >= 650) & (lam <= 1250)
        ax.plot(lam[m], sl[m], label=label)

    lam, sl = wavelength_density(ff, fs)
    m = (lam >= 650) & (lam <= 1250)
    ax.plot(lam[m], sl[m], linestyle="--", label=fixed_label)

    ax.axvspan(TARGET_MIN_NM, TARGET_MAX_NM, alpha=0.12)
    ax.set_xlim(650, 1250)
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Spectral energy density (nJ/nm)")
    ax.grid(True)
    ax.legend(ncol=2)
    fig.tight_layout()
    fig.savefig(figures / "Figure13_chirped_designs_plus_fixed.png")
    fig.savefig(figures / "Figure13_chirped_designs_plus_fixed.pdf")
    plt.close(fig)

    # ---- Direct scientific comparison: best fixed vs 20-24 chirped ----
    target_run = "FIG13_CHIRP_20_24"
    if target_run in spectra:
        _, cf, cs = spectra[target_run]

        fig, ax = plt.subplots(figsize=SINGLE_COLUMN)
        lam_f, sl_f = wavelength_density(ff, fs)
        lam_c, sl_c = wavelength_density(cf, cs)

        mf = (lam_f >= 900) & (lam_f <= 1125)
        mc = (lam_c >= 900) & (lam_c <= 1125)

        ax.plot(lam_f[mf], sl_f[mf], label="Fixed 21.75 µm")
        ax.plot(lam_c[mc], sl_c[mc], label="Chirped 20–24 µm")
        ax.axvspan(TARGET_MIN_NM, TARGET_MAX_NM, alpha=0.12)

        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Spectral energy density (nJ/nm)")
        ax.grid(True)
        ax.legend()
        fig.tight_layout()
        fig.savefig(figures / "Figure13_fixed_vs_chirped_20_24.png")
        fig.savefig(figures / "Figure13_fixed_vs_chirped_20_24.pdf")
        plt.close(fig)

    print()
    print(mdf.to_string(index=False))
    print()
    print("Saved:")
    print("  processed/figure13_metrics.csv")
    print("  figures/Figure13_chirped_designs_plus_fixed.png/.pdf")
    if target_run in spectra:
        print("  figures/Figure13_fixed_vs_chirped_20_24.png/.pdf")


if __name__ == "__main__":
    main()
