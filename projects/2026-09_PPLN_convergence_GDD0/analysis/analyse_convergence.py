#!/usr/bin/env python3
"""Analyse one-at-a-time resolution convergence for the GDD=0 PPLN case."""
from __future__ import annotations

from pathlib import Path
import csv
import zipfile

import matplotlib.pyplot as plt
from publication_style import apply_publication_style, SINGLE_COLUMN, WIDE_FIGSIZE

apply_publication_style()
import numpy as np
import pandas as pd

C = 299_792_458.0
TARGET_MIN_NM = 950.0
TARGET_MAX_NM = 1050.0
SHAPE_MIN_NM = 700.0
SHAPE_MAX_NM = 1600.0

ENERGY_TOL_PCT = 1.0
SHAPE_TOL_PCT = 2.0

ROOT = Path(__file__).resolve().parents[1]


def read_value(path: Path, label: str) -> float:
    prefix = label + ":"
    for line in path.read_text().splitlines():
        if line.startswith(prefix):
            return float(line.split(":", 1)[1].strip())
    raise ValueError(f"Could not find {label!r} in {path}")


def grid_from_input(input_path: Path):
    time_span = read_value(input_path, "Time span (s)")
    dt = read_value(input_path, "dt (s)")
    grid_width = read_value(input_path, "Grid width (m)")
    dx = read_value(input_path, "dx (m)")
    dz = read_value(input_path, "dz (m)")

    nt = int(8 * round(time_span / (8 * dt)))
    nf = nt // 2 + 1
    freq = np.fft.rfftfreq(nt, dt)
    nspace = int(round(grid_width / dx))

    return {
        "time_span": time_span,
        "dt": dt,
        "dx": dx,
        "dz": dz,
        "nt": nt,
        "nf": nf,
        "freq": freq,
        "nspace": nspace,
    }


def read_total_spectrum(run_name: str, nf: int):
    zpath = ROOT / "outputs" / f"{run_name}.zip"
    if not zpath.exists():
        raise FileNotFoundError(zpath)

    with zipfile.ZipFile(zpath) as zf:
        member = f"{run_name}_spectrum.dat"
        if member not in zf.namelist():
            matches = [n for n in zf.namelist() if n.endswith("_spectrum.dat")]
            if len(matches) != 1:
                raise RuntimeError(
                    f"{run_name}: could not uniquely locate spectrum: {matches}"
                )
            member = matches[0]
        raw = np.frombuffer(zf.read(member), dtype=np.float64)

    expected = 3 * nf
    if raw.size != expected:
        raise ValueError(
            f"{run_name}: read {raw.size} doubles, expected {expected}"
        )

    return raw.reshape((nf, 3), order="F")[:, 2]


def integrate_band(freq, spec, min_nm, max_nm):
    fmin = C / (max_nm * 1e-9)
    fmax = C / (min_nm * 1e-9)
    mask = (freq >= fmin) & (freq <= fmax)
    return np.trapezoid(spec[mask], freq[mask])


def spectral_l1_error_pct(ref_f, ref_s, test_f, test_s):
    # Compare on the test grid over the shared 700-1600 nm interval.
    fmin = C / (SHAPE_MAX_NM * 1e-9)
    fmax = C / (SHAPE_MIN_NM * 1e-9)
    lo = max(fmin, ref_f.min(), test_f.min())
    hi = min(fmax, ref_f.max(), test_f.max())

    mask = (test_f >= lo) & (test_f <= hi)
    f = test_f[mask]
    s_test = test_s[mask]

    s_ref = np.interp(f, ref_f, ref_s)
    denom = np.trapezoid(np.abs(s_ref), f)
    if denom == 0:
        return np.nan

    return 100.0 * np.trapezoid(np.abs(s_test - s_ref), f) / denom


def main():
    manifest = pd.read_csv(ROOT / "run_manifest.csv")
    data = {}

    missing = []
    for _, row in manifest.iterrows():
        run = row["run_name"]
        input_path = ROOT / row["input_file"]
        grid = grid_from_input(input_path)
        try:
            spec = read_total_spectrum(run, grid["nf"])
        except FileNotFoundError:
            missing.append(run)
            continue
        data[run] = (row, grid, spec)

    if missing:
        print("Missing runs:")
        for x in missing:
            print("  ", x)
        print()

    if "CONV_ref" not in data:
        raise SystemExit("Reference output CONV_ref is required before analysis.")

    ref_row, ref_grid, ref_spec = data["CONV_ref"]
    ref_f = ref_grid["freq"]
    ref_energy = integrate_band(
        ref_f, ref_spec, TARGET_MIN_NM, TARGET_MAX_NM
    )

    records = []
    for _, mrow in manifest.iterrows():
        run = mrow["run_name"]
        if run not in data:
            continue

        row, grid, spec = data[run]
        energy = integrate_band(
            grid["freq"], spec, TARGET_MIN_NM, TARGET_MAX_NM
        )
        energy_err = 100.0 * abs(energy - ref_energy) / abs(ref_energy)
        shape_err = spectral_l1_error_pct(
            ref_f, ref_spec, grid["freq"], spec
        )

        is_ref = run == "CONV_ref"
        passes = is_ref or (
            energy_err < ENERGY_TOL_PCT
            and shape_err < SHAPE_TOL_PCT
        )

        records.append({
            "run_name": run,
            "purpose": row["purpose"],
            "dt_fs": float(row["dt_fs"]),
            "dx_um": float(row["dx_um"]),
            "dz_um": float(row["dz_um"]),
            "Ntime": grid["nt"],
            "Nspace": grid["nspace"],
            "integrated_950_1050_pJ": energy * 1e12,
            "energy_error_pct": energy_err,
            "spectral_L1_error_700_1600_pct": shape_err,
            "passes_1pct_energy_and_2pct_shape": passes,
        })

    df = pd.DataFrame(records)

    order = [
        "CONV_ref",
        "CONV_dz_0p50",
        "CONV_dz_1p00",
        "CONV_dx_8",
        "CONV_dx_12",
        "CONV_dt_0p50",
        "CONV_dt_0p60",
    ]
    rank = {name: i for i, name in enumerate(order)}
    df["_order"] = df["run_name"].map(rank)
    df = df.sort_values("_order").drop(columns="_order")

    outdir = ROOT / "processed"
    plotdir = ROOT / "plots"
    outdir.mkdir(exist_ok=True)
    plotdir.mkdir(exist_ok=True)

    csv_path = outdir / "convergence_summary.csv"
    df.to_csv(csv_path, index=False)

    # One compact convergence diagnostic figure.
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    x = np.arange(len(df))
    ax.plot(
        x,
        df["energy_error_pct"],
        marker="o",
        label="950–1050 nm energy error",
    )
    ax.plot(
        x,
        df["spectral_L1_error_700_1600_pct"],
        marker="s",
        label="700–1600 nm spectral L1 error",
    )
    ax.axhline(ENERGY_TOL_PCT, linestyle="--", alpha=0.6, label="1% energy tolerance")
    ax.axhline(SHAPE_TOL_PCT, linestyle=":", alpha=0.6, label="2% shape tolerance")
    ax.set_xticks(x)
    ax.set_xticklabels(df["run_name"], rotation=35, ha="right")
    ax.set_ylabel("Difference from reference (%)")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()

    png = plotdir / "convergence_summary.png"
    pdf = plotdir / "convergence_summary.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    print()
    print(f"Reference 950–1050 nm energy: {ref_energy * 1e12:.6g} pJ")
    print()
    print(df.to_string(index=False))
    print()
    print("Acceptance rule:")
    print(f"  energy error < {ENERGY_TOL_PCT:g}%")
    print(f"  spectral L1 error over 700–1600 nm < {SHAPE_TOL_PCT:g}%")
    print()
    print("Saved:")
    print(f"  {csv_path.relative_to(ROOT)}")
    print(f"  {png.relative_to(ROOT)}")
    print(f"  {pdf.relative_to(ROOT)}")
    print()
    print(
        "Do not change the production grid yet. "
        "Use these one-at-a-time results to choose candidate dt/dx/dz, "
        "then run one combined-grid validation."
    )


if __name__ == "__main__":
    main()
