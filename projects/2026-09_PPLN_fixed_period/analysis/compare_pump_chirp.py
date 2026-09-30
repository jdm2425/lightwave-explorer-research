#!/usr/bin/env python3
"""Compare baseline (-1230 fs^2) and zero-GDD fixed-period PPLN scans."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from publication_style import apply_publication_style, SINGLE_COLUMN

apply_publication_style()


def main():
    here = Path.cwd().resolve()

    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", type=Path, default=here)
    parser.add_argument(
        "--zero-chirp",
        type=Path,
        default=here.parent / "2026-09_PPLN_fixed_period_GDD0",
    )
    args = parser.parse_args()

    baseline = args.baseline.resolve()
    zero = args.zero_chirp.resolve()

    base_csv = baseline / "processed" / "fixed_period_scan_summary.csv"
    zero_csv = zero / "processed" / "fixed_period_scan_summary.csv"

    if not base_csv.exists():
        raise SystemExit(f"Missing {base_csv}\nRun `lwe analyze` in the baseline project first.")
    if not zero_csv.exists():
        raise SystemExit(f"Missing {zero_csv}\nRun `lwe analyze` in the zero-GDD project first.")

    a = pd.read_csv(base_csv)[["period_um", "integrated_950_1050_nJ"]].rename(
        columns={"integrated_950_1050_nJ": "baseline_minus1230fs2_nJ"}
    )
    b = pd.read_csv(zero_csv)[["period_um", "integrated_950_1050_nJ"]].rename(
        columns={"integrated_950_1050_nJ": "zero_chirp_nJ"}
    )

    df = pd.merge(a, b, on="period_um", how="inner").sort_values("period_um")
    if df.empty:
        raise SystemExit("No common poling periods found between scans.")

    df["difference_zero_minus_baseline_nJ"] = (
        df["zero_chirp_nJ"] - df["baseline_minus1230fs2_nJ"]
    )
    df["ratio_zero_over_baseline"] = np.where(
        df["baseline_minus1230fs2_nJ"] != 0,
        df["zero_chirp_nJ"] / df["baseline_minus1230fs2_nJ"],
        np.nan,
    )

    outdir = baseline / "processed" / "chirp_comparison"
    plotdir = baseline / "plots"
    outdir.mkdir(parents=True, exist_ok=True)
    plotdir.mkdir(parents=True, exist_ok=True)

    csv_path = outdir / "pump_GDD_comparison.csv"
    df.to_csv(csv_path, index=False)

    # IMPORTANT: fixed publication-size canvas. Do not override with a large
    # presentation-style figsize; that makes the 7 pt typography look tiny.
    fig, ax = plt.subplots(figsize=SINGLE_COLUMN)

    ax.plot(
        df["period_um"],
        df["baseline_minus1230fs2_nJ"],
        marker="o",
        label="Pump GDD = -1230 fs$^2$",
    )
    ax.plot(
        df["period_um"],
        df["zero_chirp_nJ"],
        marker="s",
        label="Pump GDD = 0 fs$^2$",
    )

    ax.set_xlabel("Fixed poling period ($\\mu$m)")
    ax.set_ylabel("Integrated 950-1050 nm output (nJ)")
    ax.grid(True)
    ax.legend(loc="best")

    fig.tight_layout()

    png = plotdir / "PPLN_pump_GDD_comparison.png"
    pdf = plotdir / "PPLN_pump_GDD_comparison.pdf"
    fig.savefig(png)
    fig.savefig(pdf)
    plt.close(fig)

    print(df.to_string(index=False))
    print()
    print("Saved:")
    print(f"  {csv_path}")
    print(f"  {png}")
    print(f"  {pdf}")


if __name__ == "__main__":
    main()
