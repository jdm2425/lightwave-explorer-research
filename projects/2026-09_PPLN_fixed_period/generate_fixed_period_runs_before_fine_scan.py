#!/usr/bin/env python3
"""Generate report-matched fixed-period PPLN inputs for Lightwave Explorer.

Optical input-field parameters are taken directly from the archived LWE 2025.5
input `ChirpedPPLN_18_24um_n1230fs2_0d1uJ.txt`.

Each generated run uses one uniform first-order QPM period and an exactly
500 um physical crystal length. The final ferroelectric domain is truncated
at the output face where necessary so crystal length does not vary with period.
"""
from pathlib import Path
import csv
import math

C = 299_792_458.0

PERIODS_UM = [17.0 + 0.5*i for i in range(27)]
CRYSTAL_LENGTH_UM = 500.0
DZ_M = 2.5e-7

PUMP_CENTRE_HZ = 374.74e12
PUMP_BANDWIDTH_HZ = 25.0e12
PUMP_SG_ORDER = 4
PUMP_ENERGY_J = 1.0e-7
PUMP_GDD_S2 = -1.23e-27
PUMP_WAIST_M = 200e-6

SEED_CENTRE_HZ = 74.95e12
SEED_BANDWIDTH_HZ = 56.21e12
SEED_SG_ORDER = 10
SEED_ENERGY_J = 1.0e-8
SEED_GDD_S2 = 0.0
SEED_WAIST_M = 200e-6

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
Code version: 2025.5"""

BASE = Path(__file__).resolve().parent
INPUT_DIR = BASE / "inputs"
OUTPUT_DIR = BASE / "outputs"
INPUT_DIR.mkdir(exist_ok=True)
OUTPUT_DIR.mkdir(exist_ok=True)


def fixed_qpm_sequence(period_um: float, length_um: float):
    half = period_um / 2.0
    remaining = length_um
    segments = []

    while remaining > 1e-10:
        seg = min(half, remaining)
        segments.append(seg)
        remaining -= seg

    pieces = ["init()"]
    for i, seg in enumerate(segments):
        pieces.append(f"nonlinear(d,d,d,{seg:.12g},d)")
        if i != len(segments)-1:
            pieces.append("rotate(180)")

    if not math.isclose(sum(segments), length_um, rel_tol=0, abs_tol=1e-8):
        raise RuntimeError("Generated sequence does not equal requested crystal length.")

    return "".join(pieces), len(segments), segments[-1]


def make_input(run_name: str, period_um: float):
    sequence, n_domains, last_domain_um = fixed_qpm_sequence(
        period_um, CRYSTAL_LENGTH_UM
    )

    text = f"""Pulse energy 1 (J): {PUMP_ENERGY_J:.16g}
Pulse energy 2 (J): {SEED_ENERGY_J:.16g}
Frequency 1 (Hz): {PUMP_CENTRE_HZ:.16g}
Frequency 2 (Hz): {SEED_CENTRE_HZ:.16g}
Bandwidth 1 (Hz): {PUMP_BANDWIDTH_HZ:.16g}
Bandwidth 2 (Hz): {SEED_BANDWIDTH_HZ:.16g}
SG order 1: {PUMP_SG_ORDER}
SG order 2: {SEED_SG_ORDER}
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
Thickness (m): {CRYSTAL_LENGTH_UM*1e-6:.16g}
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
Output base path: outputs/{run_name}
Field 1 from file type: 0
Field 2 from file type: 0
Field 1 file path: 
Field 2 file path: 
Fitting reference file path: 
{MATERIAL_BLOCK}
"""
    return text, n_domains, last_domain_um


def main():
    rows = []
    for period_um in PERIODS_UM:
        tag = f"{period_um:g}".replace(".", "p")
        run_name = f"PPLN_FIXED_{tag}um"
        text, n_domains, last_domain_um = make_input(run_name, period_um)
        (INPUT_DIR / f"{run_name}.txt").write_text(text)

        rows.append({
            "run_name": run_name,
            "period_um": period_um,
            "crystal_length_um": CRYSTAL_LENGTH_UM,
            "half_period_um": period_um/2,
            "domain_segments": n_domains,
            "final_domain_length_um": last_domain_um,
            "pump_energy_uJ": PUMP_ENERGY_J*1e6,
            "pump_gdd_fs2": PUMP_GDD_S2/1e-30,
            "pump_centre_nm": C/PUMP_CENTRE_HZ*1e9,
            "pump_bandwidth_THz": PUMP_BANDWIDTH_HZ/1e12,
            "pump_sg_order": PUMP_SG_ORDER,
            "seed_centre_um": C/SEED_CENTRE_HZ*1e6,
            "seed_bandwidth_THz": SEED_BANDWIDTH_HZ/1e12,
            "seed_sg_order": SEED_SG_ORDER,
            "seed_energy_nJ": SEED_ENERGY_J*1e9,
            "input_file": f"inputs/{run_name}.txt",
        })

    with (BASE / "run_manifest.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    print(f"Generated {len(rows)} report-matched fixed-period inputs.")
    print(f"Periods: {PERIODS_UM[0]:g} to {PERIODS_UM[-1]:g} um in 0.5 um steps")
    print(f"Pump: {PUMP_ENERGY_J*1e6:g} uJ, {C/PUMP_CENTRE_HZ*1e9:.3f} nm, "
          f"{PUMP_BANDWIDTH_HZ/1e12:g} THz bandwidth parameter, SG{PUMP_SG_ORDER}, "
          f"{PUMP_GDD_S2/1e-30:g} fs^2")
    print(f"Seed: {SEED_ENERGY_J*1e9:g} nJ, {C/SEED_CENTRE_HZ*1e6:.4f} um, "
          f"{SEED_BANDWIDTH_HZ/1e12:g} THz bandwidth parameter, SG{SEED_SG_ORDER}")


if __name__ == "__main__":
    main()
