#!/usr/bin/env python3
"""Overlay spectra from the report-matched fixed-period PPLN scan.

Run from the project root, e.g.
    ~/LightwaveExplorer/.venv/bin/python analysis/overlay_fixed_period_spectra.py

Prerequisite:
    Run `lwe spectrum` for at least one completed production run so that an
    exact zero-propagation input spectrum exists under `.lwe_diagnostics/`.

Outputs:
    processed/diagnostics/PPLN_fixed_period_output_overlay.png/.pdf
    processed/diagnostics/PPLN_fixed_period_excess_overlay.png/.pdf
    processed/diagnostics/PPLN_fixed_period_overlay_metrics.csv
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import zipfile

import matplotlib.pyplot as plt
from publication_style import apply_publication_style, WIDE_FIGSIZE

apply_publication_style()
import numpy as np

C = 299_792_458.0


def read_input_value(path: Path, label: str) -> float:
    prefix = label + ":"
    for line in path.read_text().splitlines():
        if line.startswith(prefix):
            return float(line.split(":", 1)[1].strip())
    raise ValueError(f"Could not find {label!r} in {path}")


def load_total_spectrum(zip_path: Path, stem: str, nf: int) -> np.ndarray:
    member = f"{stem}_spectrum.dat"
    with zipfile.ZipFile(zip_path) as zf:
        if member not in zf.namelist():
            matches = [n for n in zf.namelist() if n.endswith("_spectrum.dat")]
            if len(matches) == 1:
                member = matches[0]
            else:
                raise FileNotFoundError(
                    f"Could not uniquely locate a spectrum in {zip_path}; "
                    f"found {matches}"
                )
        raw = np.frombuffer(zf.read(member), dtype=np.float64)

    expected = 3 * nf
    if raw.size != expected:
        raise ValueError(
            f"{zip_path.name}: {raw.size} doubles, expected {expected}"
        )
    return raw.reshape((nf, 3), order="F")[:, 2]


def find_input_diagnostic(root: Path, preferred_run: str | None) -> tuple[Path, str]:
    diag_dir = root / ".lwe_diagnostics" / "outputs"

    if preferred_run:
        stem = f"{preferred_run}_INPUT"
        z = diag_dir / f"{stem}.zip"
        if z.exists():
            return z, stem

    candidates = sorted(diag_dir.glob("*_INPUT.zip"))
    if not candidates:
        raise FileNotFoundError(
            "No exact input-spectrum diagnostic exists.\n"
            "Run, for example:\n"
            "  lwe spectrum PPLN_FIXED_18um --no-show"
        )

    z = candidates[0]
    return z, z.stem


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-nm", type=float, default=700.0)
    parser.add_argument("--max-nm", type=float, default=1600.0)
    parser.add_argument("--target-min-nm", type=float, default=950.0)
    parser.add_argument("--target-max-nm", type=float, default=1050.0)
    parser.add_argument(
        "--input-run",
        default="PPLN_FIXED_18um",
        help="Run whose cached init()-only spectrum should be used as the common input.",
    )
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    root = Path.cwd().resolve()
    manifest_path = root / "run_manifest.csv"
    if not manifest_path.exists():
        raise SystemExit("Run this script from the PPLN project root.")

    with manifest_path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    rows.sort(key=lambda r: float(r["period_um"]))

    first_input = root / rows[0]["input_file"]
    time_span = read_input_value(first_input, "Time span (s)")
    dt = read_input_value(first_input, "dt (s)")
    nt = int(8 * round(time_span / (8 * dt)))
    nf = nt // 2 + 1

    freq = np.fft.rfftfreq(nt, dt)
    keep = freq > 0
    f = freq[keep]

    input_zip, input_stem = find_input_diagnostic(root, args.input_run)
    s_in_f = load_total_spectrum(input_zip, input_stem, nf)[keep]

    lam_m = C / f
    lam_nm = lam_m * 1e9
    conversion = C / lam_m**2 * 1e3  # J/Hz -> pJ/nm

    order = np.argsort(lam_nm)
    lam_nm = lam_nm[order]
    s_in_lam = (s_in_f * conversion)[order]

    plot_mask = (lam_nm >= args.min_nm) & (lam_nm <= args.max_nm)
    x = lam_nm[plot_mask]
    input_plot = s_in_lam[plot_mask]

    target_f_min = C / (args.target_max_nm * 1e-9)
    target_f_max = C / (args.target_min_nm * 1e-9)
    target_f = (f >= target_f_min) & (f <= target_f_max)

    completed = []
    metrics = []

    for row in rows:
        run = row["run_name"]
        period = float(row["period_um"])
        zpath = root / "outputs" / f"{run}.zip"
        if not zpath.exists():
            print(f"Skipping missing output: {run}")
            continue

        s_out_f = load_total_spectrum(zpath, run, nf)[keep]
        s_out_lam = (s_out_f * conversion)[order]
        excess_lam = np.maximum(s_out_lam - s_in_lam, 0.0)

        positive_excess_pj = (
            np.trapezoid(
                np.maximum(s_out_f[target_f] - s_in_f[target_f], 0.0),
                f[target_f],
            ) * 1e12
        )
        net_target_pj = (
            np.trapezoid(
                s_out_f[target_f] - s_in_f[target_f],
                f[target_f],
            ) * 1e12
        )

        target_lam = (
            (lam_nm >= args.target_min_nm)
            & (lam_nm <= args.target_max_nm)
        )
        if np.any(target_lam):
            target_excess = excess_lam[target_lam]
            if np.nanmax(target_excess) > 0:
                peak_excess_nm = float(
                    lam_nm[target_lam][np.nanargmax(target_excess)]
                )
            else:
                peak_excess_nm = np.nan
        else:
            peak_excess_nm = np.nan

        completed.append((period, run, s_out_lam, excess_lam))
        metrics.append(
            {
                "period_um": period,
                "run_name": run,
                "positive_excess_950_1050_pJ": positive_excess_pj,
                "net_output_minus_input_950_1050_pJ": net_target_pj,
                "peak_positive_excess_in_target_nm": peak_excess_nm,
            }
        )

    if not completed:
        raise SystemExit("No production output ZIPs were found.")

    outdir = root / "processed" / "diagnostics"
    outdir.mkdir(parents=True, exist_ok=True)

    # ---------- Raw input/output overlay ----------
    fig, ax = plt.subplots(figsize=WIDE_FIGSIZE)
    positive = input_plot[input_plot > 0]
    floor_seed = positive.max() if positive.size else 1.0

    all_positive = [positive]
    for _, _, s_out_lam, _ in completed:
        y = s_out_lam[plot_mask]
        all_positive.append(y[y > 0])
    all_positive = [a for a in all_positive if a.size]
    global_max = max(float(np.nanmax(a)) for a in all_positive) if all_positive else floor_seed
    floor = global_max * 1e-12

    ax.semilogy(
        x,
        np.maximum(input_plot, floor),
        linestyle="--",
        
        label="Input",
    )
    for period, _, s_out_lam, _ in completed:
        ax.semilogy(
            x,
            np.maximum(s_out_lam[plot_mask], floor),
            
            label=f"{period:g} µm",
        )

    ax.axvspan(args.target_min_nm, args.target_max_nm, alpha=0.10)
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Spectral energy density (pJ/nm)")
    ax.set_xlim(args.min_nm, args.max_nm)
    ax.grid(True, which="both")
    ax.legend(ncol=2)
    fig.tight_layout()

    output_png = outdir / "PPLN_fixed_period_output_overlay.png"
    output_pdf = outdir / "PPLN_fixed_period_output_overlay.pdf"
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    fig.savefig(output_pdf, bbox_inches="tight")
    if args.show:
        plt.show()
    else:
        plt.close(fig)

    # ---------- Positive nonlinear excess overlay ----------
    fig, ax = plt.subplots(figsize=WIDE_FIGSIZE)
    excess_values = []
    for _, _, _, excess_lam in completed:
        y = excess_lam[plot_mask]
        excess_values.append(y[y > 0])

    excess_values = [a for a in excess_values if a.size]
    excess_max = max(float(np.nanmax(a)) for a in excess_values) if excess_values else 1.0
    excess_floor = excess_max * 1e-10

    for period, _, _, excess_lam in completed:
        y = excess_lam[plot_mask]
        ax.semilogy(
            x,
            np.where(y > excess_floor, y, np.nan),
            
            label=f"{period:g} µm",
        )

    ax.axvspan(args.target_min_nm, args.target_max_nm, alpha=0.10)
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Positive output - input (pJ/nm)")
    ax.set_xlim(args.min_nm, args.max_nm)
    ax.grid(True, which="both")
    ax.legend(ncol=2)
    fig.tight_layout()

    excess_png = outdir / "PPLN_fixed_period_excess_overlay.png"
    excess_pdf = outdir / "PPLN_fixed_period_excess_overlay.pdf"
    fig.savefig(excess_png, dpi=300, bbox_inches="tight")
    fig.savefig(excess_pdf, bbox_inches="tight")
    if args.show:
        plt.show()
    else:
        plt.close(fig)

    metrics_path = outdir / "PPLN_fixed_period_overlay_metrics.csv"
    with metrics_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=metrics[0].keys())
        writer.writeheader()
        writer.writerows(metrics)

    print()
    print(f"Plotted {len(completed)} completed runs.")
    print("Saved:")
    print(f"  {output_png.relative_to(root)}")
    print(f"  {output_pdf.relative_to(root)}")
    print(f"  {excess_png.relative_to(root)}")
    print(f"  {excess_pdf.relative_to(root)}")
    print(f"  {metrics_path.relative_to(root)}")


if __name__ == "__main__":
    main()
