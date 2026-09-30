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

    return freq, raw.reshape((nf, 3), order="F")[:, 2]


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
    test = test_s[m]
    ref = np.interp(f, ref_f, ref_s)

    denom = np.trapezoid(np.abs(ref), f)
    return 100.0 * np.trapezoid(np.abs(test - ref), f) / denom


def compare(period, reference, candidate):
    try:
        rf, rs = load(reference)
        cf, cs = load(candidate)
    except FileNotFoundError:
        return None

    er = band_energy(rf, rs)
    ec = band_energy(cf, cs)
    eerr = 100.0 * abs(ec - er) / abs(er)
    serr = l1_shape(rf, rs, cf, cs)

    return {
        "period_um": period,
        "reference_run": reference,
        "candidate_run": candidate,
        "reference_950_1050_pJ": er * 1e12,
        "candidate_950_1050_pJ": ec * 1e12,
        "energy_error_pct": eerr,
        "spectral_L1_error_700_1600_pct": serr,
        "passes": eerr < ENERGY_TOL_PCT and serr < SHAPE_TOL_PCT,
    }


def main():
    results = []

    r23 = compare(23.0, "CONV_spot23_ref", "CONV_candidate_23")
    if r23 is not None:
        results.append(r23)

    rpeak = compare(21.75, "CONV_ref", "CONV_candidate_21p75")
    if rpeak is not None:
        results.append(rpeak)

    if not results:
        raise SystemExit("Run CONV_candidate_23 first.")

    df = pd.DataFrame(results)
    path = ROOT / "processed" / "convergence_phase3_summary.csv"
    df.to_csv(path, index=False)

    print(df.to_string(index=False))
    print()
    print("Pass criteria:")
    print("  target-band energy error < 1%")
    print("  spectral L1 error over 700–1600 nm < 2%")
    print()
    print(f"Saved: {path.relative_to(ROOT)}")

    if r23 is not None:
        print()
        if r23["passes"]:
            print("23.0 um candidate grid PASSES.")
            if rpeak is None:
                print("Next run CONV_candidate_21p75.")
        else:
            print("23.0 um candidate grid FAILS.")
            print("Do not use dt=0.60 fs, dx=8 um, dz=1.00 um for production.")

    if r23 is not None and rpeak is not None:
        print()
        if r23["passes"] and rpeak["passes"]:
            print(
                "Both validation points PASS. "
                "Use dt=0.60 fs, dx=8 um, dz=1.00 um for the GDD=0 production scan."
            )
        else:
            print("At least one validation point failed; keep refining the candidate grid.")


if __name__ == "__main__":
    main()
