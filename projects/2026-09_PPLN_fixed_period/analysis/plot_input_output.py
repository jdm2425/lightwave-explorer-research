#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import zipfile

import matplotlib.pyplot as plt
import numpy as np

C = 299_792_458.0


def read_value(path: Path, label: str) -> float:
    prefix = label + ":"
    for line in path.read_text().splitlines():
        if line.startswith(prefix):
            return float(line.split(":", 1)[1].strip())
    raise ValueError(f"Could not find {label!r} in {path}")


def read_total_spectrum(zip_path: Path, member_stem: str, nf: int) -> np.ndarray:
    member = f"{member_stem}_spectrum.dat"
    with zipfile.ZipFile(zip_path) as zf:
        if member not in zf.namelist():
            matches = [n for n in zf.namelist() if n.endswith("_spectrum.dat")]
            raise FileNotFoundError(
                f"{member} not found in {zip_path}. Available: {matches}"
            )
        raw = np.frombuffer(zf.read(member), dtype=np.float64)

    expected = 3 * nf
    print(f"{member_stem}: read {raw.size} doubles (expected {expected})")
    if raw.size != expected:
        raise ValueError(
            f"Unexpected spectrum length: {raw.size}, expected {expected}"
        )

    return raw.reshape((nf, 3), order="F")[:, 2]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("run_name")
    ap.add_argument("--input-file", required=True)
    ap.add_argument("--output-zip", required=True)
    ap.add_argument("--input-zip", required=True)
    ap.add_argument("--input-member-stem", required=True)
    ap.add_argument("--min-nm", type=float, default=500.0)
    ap.add_argument("--max-nm", type=float, default=6000.0)
    ap.add_argument("--target-min-nm", type=float, default=950.0)
    ap.add_argument("--target-max-nm", type=float, default=1050.0)
    ap.add_argument("--no-show", action="store_true")
    args = ap.parse_args()

    root = Path.cwd().resolve()
    input_file = Path(args.input_file)
    if not input_file.is_absolute():
        input_file = (root / input_file).resolve()

    time_span = read_value(input_file, "Time span (s)")
    dt = read_value(input_file, "dt (s)")

    nt = int(8 * round(time_span / (8 * dt)))
    nf = nt // 2 + 1
    freq = np.fft.rfftfreq(nt, dt)

    s_in_f = read_total_spectrum(Path(args.input_zip), args.input_member_stem, nf)
    s_out_f = read_total_spectrum(Path(args.output_zip), args.run_name, nf)

    positive = freq > 0
    f = freq[positive]
    s_in_f = s_in_f[positive]
    s_out_f = s_out_f[positive]

    target_f_min = C / (args.target_max_nm * 1e-9)
    target_f_max = C / (args.target_min_nm * 1e-9)
    band = (f >= target_f_min) & (f <= target_f_max)

    input_target_pj = np.trapezoid(s_in_f[band], f[band]) * 1e12
    output_target_pj = np.trapezoid(s_out_f[band], f[band]) * 1e12
    positive_excess_pj = (
        np.trapezoid(np.maximum(s_out_f[band] - s_in_f[band], 0.0), f[band])
        * 1e12
    )
    net_target_pj = output_target_pj - input_target_pj

    lam_m = C / f
    lam_nm = lam_m * 1e9
    conversion = C / lam_m**2 * 1e3  # J/Hz -> pJ/nm

    s_in_lam = s_in_f * conversion
    s_out_lam = s_out_f * conversion
    s_excess_lam = np.maximum(s_out_lam - s_in_lam, 0.0)

    order = np.argsort(lam_nm)
    lam_nm = lam_nm[order]
    s_in_lam = s_in_lam[order]
    s_out_lam = s_out_lam[order]
    s_excess_lam = s_excess_lam[order]

    pmask = (lam_nm >= args.min_nm) & (lam_nm <= args.max_nm)
    x = lam_nm[pmask]
    yin = s_in_lam[pmask]
    yout = s_out_lam[pmask]
    yex = s_excess_lam[pmask]

    positives = np.concatenate([yin[yin > 0], yout[yout > 0]])
    floor = positives.max() * 1e-12 if positives.size else 1e-30

    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.semilogy(x, np.maximum(yin, floor), "--", lw=1.5, label="Input")
    ax.semilogy(x, np.maximum(yout, floor), lw=1.7, label="After PPLN")
    ax.semilogy(
        x,
        np.where(yex > floor, yex, np.nan),
        lw=1.2,
        alpha=0.8,
        label="Positive output - input",
    )
    ax.axvspan(
        args.target_min_nm,
        args.target_max_nm,
        alpha=0.12,
        label=f"{args.target_min_nm:g}-{args.target_max_nm:g} nm target",
    )
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Spectral energy density (pJ/nm)")
    ax.set_xlim(args.min_nm, args.max_nm)
    ax.grid(True, which="both", alpha=0.2)
    ax.legend()
    fig.tight_layout()

    outdir = root / "processed" / "diagnostics"
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / f"{args.run_name}_input_vs_output_log.png"
    fig.savefig(outpath, dpi=300, bbox_inches="tight")

    print()
    print(f"{args.target_min_nm:g}-{args.target_max_nm:g} nm:")
    print(f"  input energy:              {input_target_pj:.6g} pJ")
    print(f"  output energy:             {output_target_pj:.6g} pJ")
    print(f"  positive spectral excess:  {positive_excess_pj:.6g} pJ")
    print(f"  net output - input:        {net_target_pj:.6g} pJ")
    print()
    print(f"Saved: {outpath.relative_to(root)}")

    if args.no_show:
        plt.close(fig)
    else:
        plt.show()


if __name__ == "__main__":
    main()
