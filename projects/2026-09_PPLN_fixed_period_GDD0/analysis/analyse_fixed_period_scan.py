#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import zipfile
import matplotlib.pyplot as plt
from publication_style import apply_publication_style, SINGLE_COLUMN

apply_publication_style()
import numpy as np
import pandas as pd

C = 299_792_458.0
LAM_MIN_M = 950e-9
LAM_MAX_M = 1050e-9

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "outputs"
PROCESSED = BASE / "processed"
PLOTS = BASE / "plots"
PROCESSED.mkdir(exist_ok=True)
PLOTS.mkdir(exist_ok=True)
MANIFEST = pd.read_csv(BASE / "run_manifest.csv")


def read_input_value(path: Path, label: str) -> float:
    prefix = label + ":"
    for line in path.read_text().splitlines():
        if line.startswith(prefix):
            return float(line.split(":", 1)[1].strip())
    raise ValueError(f"Could not find {label!r} in {path}")


def frequency_grid_for_run(input_path: Path) -> np.ndarray:
    time_span_s = read_input_value(input_path, "Time span (s)")
    dt_s = read_input_value(input_path, "dt (s)")
    nt = int(8 * round(time_span_s / (8 * dt_s)))
    return np.fft.rfftfreq(nt, dt_s)


def locate_spectrum(run_name: str):
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
                exact = [n for n in members if Path(n).name == f"{run_name}_spectrum.dat"]
                return ("zip", zpath, (exact or members)[0])

    raise FileNotFoundError(f"Could not find {run_name}_spectrum.dat in {OUT}.")


def read_total_spectrum(run_name: str, nf: int) -> np.ndarray:
    mode, path, member = locate_spectrum(run_name)
    if mode == "raw":
        a = np.fromfile(path, dtype=np.float64)
    else:
        with zipfile.ZipFile(path) as zf:
            a = np.frombuffer(zf.read(member), dtype=np.float64)

    expected = 3 * nf
    if a.size != expected:
        raise ValueError(
            f"{run_name}: spectrum contains {a.size} doubles; expected {expected}. "
            "The input file used for analysis must match the grid used to produce the output."
        )
    return a.reshape((nf, 3), order="F")[:, 2]


def main():
    records = []
    missing = []

    for _, row in MANIFEST.sort_values("period_um").iterrows():
        run = row["run_name"]
        input_path = BASE / row["input_file"] if "input_file" in row and isinstance(row["input_file"], str) else BASE / "inputs" / f"{run}.txt"

        freq_hz = frequency_grid_for_run(input_path)
        try:
            spectrum = read_total_spectrum(run, len(freq_hz))
        except FileNotFoundError:
            missing.append(run)
            continue

        band = (freq_hz >= C / LAM_MAX_M) & (freq_hz <= C / LAM_MIN_M)
        energy_j = np.trapezoid(spectrum[band], freq_hz[band])

        records.append({
            "period_um": float(row["period_um"]),
            "integrated_950_1050_J": energy_j,
            "integrated_950_1050_nJ": energy_j * 1e9,
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

    fig, ax = plt.subplots(figsize=SINGLE_COLUMN)
    ax.plot(df["period_um"], df["integrated_950_1050_nJ"], marker="o", markersize=3.5, linewidth=1.4)
    ax.set_xlabel("Fixed poling period (µm)")
    ax.set_ylabel("Integrated 950–1050 nm output (nJ)")
    ax.grid(True)
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
