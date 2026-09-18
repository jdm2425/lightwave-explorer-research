#!/usr/bin/env python3
"""Generate LWE 2025.2 inputs for a fixed-poling-period PPLN scan.

Purpose
-------
Reproduce the ESA PPLN comparison conditions while changing only the
uniform QPM period.  One 500 um crystal is generated for each period from
17 to 30 um (1 um increments by default).

The nonlinear sequence is explicitly unrolled so that the propagated
crystal length is exactly 500 um even when 500 um is not an integer number
of poling periods.  The final ferroelectric domain is truncated at the end
face of the crystal, as it would be for a real cut crystal.
"""
from __future__ import annotations

from pathlib import Path
import csv
import math

# -----------------------------------------------------------------------------
# USER-REVIEWABLE SETTINGS
# These match the values stated in the ESA report for the Figure-13 comparison.
# -----------------------------------------------------------------------------
C = 299_792_458.0

PERIODS_UM = [float(x) for x in range(17, 31)]  # 17, 18, ..., 30 um
CRYSTAL_LENGTH_UM = 500.0
DZ_M = 2.5e-7

# 800 nm pump, ~100 nm bandwidth, 0.025 uJ, -500 fs^2 GDD.
PUMP_CENTRE_HZ = 374.74e12
PUMP_BANDWIDTH_HZ = 46.84e12
PUMP_ENERGY_J = 0.025e-6
PUMP_GDD_S2 = -500.0e-30
PUMP_WAIST_M = 200e-6

# Artificial coherent MIR seed used in the report as a surrogate for
# vacuum-initiated OPG: 56 THz centre, 75 THz bandwidth, 10 nJ.
SEED_CENTRE_HZ = 56.0e12
SEED_BANDWIDTH_HZ = 75.0e12
SEED_ENERGY_J = 10e-9
SEED_GDD_S2 = 0.0
SEED_WAIST_M = 150e-6

# Same numerical grid / material block as the archived Chimera PPLN LWE 2025.2
# input used elsewhere in the project.
GRID_WIDTH_M = 0.000768
DX_M = 6e-6
TIME_SPAN_S = 2.9988e-12
DT_S = 4.5e-16

MATERIAL_BLOCK = """Material name: MgO:LiNbO3
Sellmeier reference: O. Gayer, et al., Appl. Phys. B 91, 343-348 (2008)
Chi2 reference: Nikogosian, Nonlinear Optical Crystals - A Complete Survey, values from stoichimetric, except for d33, reduced to 26, based on https://arxiv.org/abs/2111.13796
Chi3 reference: https://arxiv.org/abs/2111.13796
5.65089,0.118417,0,-0.043728,89.6158,0,-117.722,0,0,0,0,0,0,-0.0197,0,0,0,0,0,0,0,0
5.7484,0.098175,0,-0.0407384,188.917,0,-156.75,0,0,0,0,0,0,-0.0132,0,0,0,0,0,0,0,0
5.7484,0.098175,0,-0.0407384,188.917,0,-156.75,0,0,0,0,0,0,-0.0132,0,0,0,0,0,0,0,0
Code version: 2025.2"""

BASE = Path(__file__).resolve().parent
INPUT_DIR = BASE / "inputs"
OUTPUT_DIR = BASE / "outputs"
INPUT_DIR.mkdir(exist_ok=True)
OUTPUT_DIR.mkdir(exist_ok=True)


def fixed_qpm_sequence(period_um: float, length_um: float) -> tuple[str, int, float]:
    """Return an exactly-length-matched alternating-domain LWE sequence.

    For first-order QPM the sign of d_eff is reversed every Lambda/2.  We use
    the same rotate(180) convention as the existing Chimera PPLN model.
    """
    half = period_um / 2.0
    remaining = length_um
    segments: list[float] = []

    # Explicitly fill the crystal, truncating only the last domain if needed.
    while remaining > 1e-10:
        seg = min(half, remaining)
        segments.append(seg)
        remaining -= seg

    pieces = ["init()"]
    for i, seg in enumerate(segments):
        pieces.append(f"nonlinear(d,d,d,{seg:.12g},d)")
        if i != len(segments) - 1:
            pieces.append("rotate(180)")

    total = sum(segments)
    if not math.isclose(total, length_um, rel_tol=0.0, abs_tol=1e-8):
        raise RuntimeError(f"Sequence length mismatch: {total} um vs {length_um} um")

    return "".join(pieces), len(segments), segments[-1]


def lwe_text(run_name: str, period_um: float) -> tuple[str, int, float]:
    sequence, n_domains, last_domain_um = fixed_qpm_sequence(period_um, CRYSTAL_LENGTH_UM)

    # Relative output paths work when run_all.zsh is launched from this folder.
    output_base = f"outputs/{run_name}"

    text = f"""Pulse energy 1 (J): {PUMP_ENERGY_J:.16g}
Pulse energy 2 (J): {SEED_ENERGY_J:.16g}
Frequency 1 (Hz): {PUMP_CENTRE_HZ:.16g}
Frequency 2 (Hz): {SEED_CENTRE_HZ:.16g}
Bandwidth 1 (Hz): {PUMP_BANDWIDTH_HZ:.16g}
Bandwidth 2 (Hz): {SEED_BANDWIDTH_HZ:.16g}
SG order 1: 2
SG order 2: 2
CEP 1 (rad): 0
CEP 2 (rad): 0
Delay 1 (s): 0
Delay 2 (s): 0
GDD 1 (s^-2): {PUMP_GDD_S2:.16g}
GDD 2 (s^-2): {SEED_GDD_S2:.16g}
TOD 1 (s^-3): 0
TOD 2 (s^-3): 0
Phase material 1 index: 2
Phase material 2 index: 2
Phase material thickness 1 (mcr.): 0
Phase material thickness 2 (mcr.): 0
Beam mode placeholder: 0
Beamwaist 1 (m): {PUMP_WAIST_M:.16g}
Beamwaist 2 (m): {SEED_WAIST_M:.16g}
x offset 1 (m): 0
x offset 2 (m): 0
y offset 1 (m): 0
y offset 2 (m): 0
z offset 1 (m): 0
z offset 2 (m): 0
NC angle 1 (rad): 0
NC angle 2 (rad): 0
NC angle phi 1 (rad): 0
NC angle phi 2 (rad): 0
Polarization 1 (rad): 0
Polarization 2 (rad): 0
Circularity 1: 0
Circularity 2: 0
Material index: 12
Alternate material index: 0
Crystal theta (rad): 1.5707963267949
Crystal phi (rad): 0
Grid width (m): {GRID_WIDTH_M:.16g}
Grid height (m): 0
dx (m): {DX_M:.16g}
Time span (s): {TIME_SPAN_S:.16g}
dt (s): {DT_S:.16g}
Thickness (m): {CRYSTAL_LENGTH_UM * 1e-6:.16g}
dz (m): {DZ_M:.16g}
Nonlinear absorption parameter: 0
Initial carrier density (m^-3): 0
Band gap (eV): 6
Effective mass (relative): 0
Drude gamma (Hz): 0
Propagation mode: 1
Batch mode: 0
Batch destination: 1e-05
Batch steps: 1
Batch mode 2: 0
Batch destination 2: 1000
Batch steps 2: 1
Sequence: {sequence}
Fitting: 
Fitting mode: 0
Output base path: {output_base}
Field 1 from file type: 0
Field 2 from file type: 0
Field 1 file path: 
Field 2 file path: 
Fitting reference file path: 
{MATERIAL_BLOCK}
"""
    return text, n_domains, last_domain_um


def main() -> None:
    rows = []
    for period_um in PERIODS_UM:
        tag = f"{period_um:g}".replace(".", "p")
        run_name = f"PPLN_FIXED_{tag}um"
        text, n_domains, last_domain_um = lwe_text(run_name, period_um)
        path = INPUT_DIR / f"{run_name}.txt"
        path.write_text(text)
        rows.append({
            "run_name": run_name,
            "period_um": period_um,
            "crystal_length_um": CRYSTAL_LENGTH_UM,
            "half_period_um": period_um / 2.0,
            "domain_segments": n_domains,
            "final_domain_length_um": last_domain_um,
            "pump_energy_uJ": PUMP_ENERGY_J * 1e6,
            "pump_gdd_fs2": PUMP_GDD_S2 / 1e-30,
            "pump_centre_nm": C / PUMP_CENTRE_HZ * 1e9,
            "seed_centre_um": C / SEED_CENTRE_HZ * 1e6,
            "seed_energy_nJ": SEED_ENERGY_J * 1e9,
            "input_file": f"inputs/{run_name}.txt",
        })

    with (BASE / "run_manifest.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    print(f"Generated {len(rows)} LWE inputs in {INPUT_DIR}")
    print(f"Fixed crystal length: {CRYSTAL_LENGTH_UM:.1f} um")
    print(f"Periods: {PERIODS_UM[0]:.1f} to {PERIODS_UM[-1]:.1f} um")
    print("Review run_manifest.csv, then execute: lwe run")


if __name__ == "__main__":
    main()
