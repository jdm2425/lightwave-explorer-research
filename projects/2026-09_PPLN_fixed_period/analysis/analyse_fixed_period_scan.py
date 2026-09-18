#!/usr/bin/env python3
"""Analyse the LWE fixed-period PPLN scan and make the proposed ESA figure.

Outputs
-------
fixed_period_scan_summary.csv
PPLN_fixed_period_summary.png
PPLN_fixed_period_summary.pdf

The left y-axis is the integrated output energy in 950--1050 nm.
The right y-axis is the peak wavelength within that band.

Important: the energy integral is performed on LWE's native frequency axis.
For the peak-wavelength metric, the spectrum is converted to a wavelength
spectral density using the Jacobian |df/dlambda| = c/lambda^2 before taking
its maximum.
"""
from __future__ import annotations

from pathlib import Path
import zipfile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

C = 299_792_458.0
TIME_SPAN_S = 2.9988e-12
DT_S = 4.5e-16
NT = int(8 * round(TIME_SPAN_S / (8 * DT_S)))
NF = NT // 2 + 1
FREQ_HZ = np.fft.rfftfreq(NT, DT_S)

LAM_MIN_M = 950e-9
LAM_MAX_M = 1050e-9
FMIN_HZ = C / LAM_MAX_M
FMAX_HZ = C / LAM_MIN_M
BAND = (FREQ_HZ >= FMIN_HZ) & (FREQ_HZ <= FMAX_HZ)

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "outputs"
PROCESSED = BASE / "processed"
PLOTS = BASE / "plots"
PROCESSED.mkdir(exist_ok=True)
PLOTS.mkdir(exist_ok=True)
MANIFEST = pd.read_csv(BASE / "run_manifest.csv")


def locate_spectrum(run_name: str):
    """Find an LWE *_spectrum.dat either raw or inside a run ZIP."""
    candidates = [
        OUT / f"{run_name}_spectrum.dat",
        BASE / f"{run_name}_spectrum.dat",
    ]
    for path in candidates:
        if path.exists():
            return ("raw", path, None)

    recursive = list(OUT.rglob(f"{run_name}_spectrum.dat"))
    if recursive:
        return ("raw", recursive[0], None)

    zip_candidates = [OUT / f"{run_name}.zip", BASE / f"{run_name}.zip"]
    zip_candidates += list(OUT.rglob(f"{run_name}.zip"))
    for zpath in zip_candidates:
        if not zpath.exists():
            continue
        with zipfile.ZipFile(zpath) as zf:
            members = [n for n in zf.namelist() if n.endswith("_spectrum.dat")]
            if members:
                # Prefer exact basename when present.
                exact = [n for n in members if Path(n).name == f"{run_name}_spectrum.dat"]
                return ("zip", zpath, (exact or members)[0])

    raise FileNotFoundError(
        f"Could not find {run_name}_spectrum.dat in {OUT}. "
        "Place the LWE outputs in the outputs/ folder."
    )


def read_total_spectrum(run_name: str) -> np.ndarray:
    mode, path, member = locate_spectrum(run_name)
    if mode == "raw":
        a = np.fromfile(path, dtype=np.float64)
    else:
        with zipfile.ZipFile(path) as zf:
            a = np.frombuffer(zf.read(member), dtype=np.float64)

    expected = 3 * NF
    if a.size != expected:
        raise ValueError(
            f"{run_name}: spectrum contains {a.size} doubles; expected {expected}. "
            "This usually means the LWE time grid or batch dimensions differ from the input files."
        )

    # LWE layout used by the existing project analysis:
    # frequency index x [polarisation 1, polarisation 2, total], Fortran order.
    return a.reshape((NF, 3), order="F")[:, 2]


def band_metrics(s_f: np.ndarray):
    f = FREQ_HZ[BAND]
    sf = s_f[BAND]

    # Physical band energy: native-frequency integral.
    energy_j = np.trapezoid(sf, f)

    # Convert J/Hz -> J/m before asking for the peak on a wavelength axis.
    lam_m = C / f
    s_lambda = sf * C / lam_m**2
    peak_i = int(np.nanargmax(s_lambda))
    peak_nm = lam_m[peak_i] * 1e9

    # Also retain the raw-frequency-density peak for diagnostics.
    peak_raw_nm = (C / f[int(np.nanargmax(sf))]) * 1e9
    return energy_j, peak_nm, peak_raw_nm


def main():
    records = []
    missing = []

    for _, row in MANIFEST.sort_values("period_um").iterrows():
        run = row["run_name"]
        try:
            s = read_total_spectrum(run)
        except FileNotFoundError:
            missing.append(run)
            continue

        energy_j, peak_nm, peak_raw_nm = band_metrics(s)
        records.append({
            "period_um": float(row["period_um"]),
            "integrated_950_1050_J": energy_j,
            "integrated_950_1050_pJ": energy_j * 1e12,
            "peak_wavelength_nm": peak_nm,
            "peak_wavelength_raw_frequency_density_nm": peak_raw_nm,
        })

    if missing:
        print("Missing runs:")
        for run in missing:
            print("  ", run)
        print()

    if not records:
        raise SystemExit("No LWE spectra were found. Run the simulations first.")

    df = pd.DataFrame(records).sort_values("period_um")
    df.to_csv(PROCESSED / "fixed_period_scan_summary.csv", index=False)

    fig, ax1 = plt.subplots(figsize=(7.2, 4.6))
    line1, = ax1.plot(
        df["period_um"],
        df["integrated_950_1050_pJ"],
        marker="o",
        label="Integrated 950-1050 nm output",
    )
    ax1.set_xlabel("Fixed poling period (um)")
    ax1.set_ylabel("Integrated 950-1050 nm output (pJ)")
    ax1.grid(True, alpha=0.25)

    ax2 = ax1.twinx()
    line2, = ax2.plot(
        df["period_um"],
        df["peak_wavelength_nm"],
        marker="s",
        linestyle="--",
        label="Peak wavelength in target band",
    )
    ax2.set_ylabel("Peak wavelength in 950-1050 nm band (nm)")
    ax2.set_ylim(945, 1055)

    ax1.legend([line1, line2], [line1.get_label(), line2.get_label()], loc="best")
    fig.tight_layout()
    fig.savefig(PLOTS / "PPLN_fixed_period_summary.png", dpi=300, bbox_inches="tight")
    fig.savefig(PLOTS / "PPLN_fixed_period_summary.pdf", bbox_inches="tight")
    plt.close(fig)

    print(df.to_string(index=False))
    print()
    print("Saved:")
    print("  processed/fixed_period_scan_summary.csv")
    print("  plots/PPLN_fixed_period_summary.png")
    print("  plots/PPLN_fixed_period_summary.pdf")


if __name__ == "__main__":
    main()
