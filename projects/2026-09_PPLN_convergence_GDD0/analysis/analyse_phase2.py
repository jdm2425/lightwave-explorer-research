#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import zipfile
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


def load(run_name: str):
    inp = ROOT / "inputs" / f"{run_name}.txt"
    if not inp.exists():
        raise FileNotFoundError(inp)

    time_span = read_value(inp, "Time span (s)")
    dt = read_value(inp, "dt (s)")
    nt = int(8 * round(time_span / (8 * dt)))
    nf = nt // 2 + 1
    freq = np.fft.rfftfreq(nt, dt)

    zpath = ROOT / "outputs" / f"{run_name}.zip"
    if not zpath.exists():
        raise FileNotFoundError(zpath)

    with zipfile.ZipFile(zpath) as zf:
        member = f"{run_name}_spectrum.dat"
        if member not in zf.namelist():
            matches = [n for n in zf.namelist() if n.endswith("_spectrum.dat")]
            if len(matches) != 1:
                raise RuntimeError(f"{run_name}: spectrum not unique: {matches}")
            member = matches[0]
        raw = np.frombuffer(zf.read(member), dtype=np.float64)

    expected = 3 * nf
    if raw.size != expected:
        raise ValueError(f"{run_name}: {raw.size} doubles, expected {expected}")

    spec = raw.reshape((nf, 3), order="F")[:, 2]
    return freq, spec


def band_energy(freq, spec):
    fmin = C / (TARGET_MAX_NM * 1e-9)
    fmax = C / (TARGET_MIN_NM * 1e-9)
    m = (freq >= fmin) & (freq <= fmax)
    return np.trapezoid(spec[m], freq[m])


def l1_shape(ref_f, ref_s, test_f, test_s):
    fmin = C / (SHAPE_MAX_NM * 1e-9)
    fmax = C / (SHAPE_MIN_NM * 1e-9)
    lo = max(fmin, ref_f.min(), test_f.min())
    hi = min(fmax, ref_f.max(), test_f.max())

    m = (test_f >= lo) & (test_f <= hi)
    f = test_f[m]
    t = test_s[m]
    r = np.interp(f, ref_f, ref_s)

    denom = np.trapezoid(np.abs(r), f)
    return 100.0 * np.trapezoid(np.abs(t - r), f) / denom


def compare(label, ref_name, test_name):
    try:
        rf, rs = load(ref_name)
        tf, ts = load(test_name)
    except FileNotFoundError:
        return None

    re = band_energy(rf, rs)
    te = band_energy(tf, ts)
    eerr = 100.0 * abs(te - re) / abs(re)
    serr = l1_shape(rf, rs, tf, ts)
    passes = (eerr < ENERGY_TOL_PCT) and (serr < SHAPE_TOL_PCT)

    return {
        "comparison": label,
        "reference_run": ref_name,
        "test_run": test_name,
        "reference_950_1050_pJ": re * 1e12,
        "test_950_1050_pJ": te * 1e12,
        "energy_error_pct": eerr,
        "spectral_L1_error_700_1600_pct": serr,
        "passes": passes,
    }


def main():
    results = []

    r1 = compare(
        "combined grid at 21.75 um",
        "CONV_ref",
        "CONV_combined_21p75",
    )
    if r1 is not None:
        results.append(r1)

    r2 = compare(
        "combined grid at 23.0 um",
        "CONV_spot23_ref",
        "CONV_spot23_combined",
    )
    if r2 is not None:
        results.append(r2)

    if not results:
        raise SystemExit(
            "No complete phase-2 comparison pair found. "
            "Run CONV_combined_21p75 first."
        )

    df = pd.DataFrame(results)
    out = ROOT / "processed" / "convergence_phase2_summary.csv"
    df.to_csv(out, index=False)

    print(df.to_string(index=False))
    print()
    print("Pass criteria:")
    print("  target-band energy error < 1%")
    print("  spectral L1 error over 700–1600 nm < 2%")
    print()
    print(f"Saved: {out.relative_to(ROOT)}")

    if results[0]["passes"]:
        print()
        print("21.75 um combined grid PASSES.")
        if len(results) == 1:
            print("Next run the 23 um reference and combined spot-check.")
        elif results[1]["passes"]:
            print("23 um independent spot-check also PASSES.")
            print(
                "Validated production grid: "
                "dt=0.60 fs, dx=12 um, dz=1.00 um."
            )
    else:
        print()
        print("21.75 um combined grid FAILS. Do not use it for production.")


if __name__ == "__main__":
    main()
