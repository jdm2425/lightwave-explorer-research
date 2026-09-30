#!/usr/bin/env python3
"""Generate exact-length chirped-PPLN seeding diagnostic simulations.

The grating is defined by a linearly varying local full QPM period Lambda(z).
Domain boundaries satisfy

    integral_0^z 2 / Lambda(z') dz' = integer,

so each boundary corresponds to a pi shift of the grating phase. The final
domain is truncated exactly at the requested physical crystal length.

This avoids the integer-loop truncation in the earlier exploratory chirped
sequence and makes forward/reverse crystals exactly equal in length.
"""
from __future__ import annotations

from pathlib import Path
import csv
import math
import re
import sys

if len(sys.argv) != 3:
    raise SystemExit("Usage: generate_inputs.py TEMPLATE.txt PROJECT_DIR")

template_path = Path(sys.argv[1]).expanduser().resolve()
root = Path(sys.argv[2]).expanduser().resolve()
text0 = template_path.read_text()

(root / "inputs").mkdir(parents=True, exist_ok=True)
(root / "provenance").mkdir(parents=True, exist_ok=True)


def get(label: str) -> float:
    prefix = label + ":"
    for line in text0.splitlines():
        if line.startswith(prefix):
            return float(line.split(":", 1)[1].strip())
    raise RuntimeError(f"Could not find {label!r}")


def replace(text: str, label: str, value) -> str:
    pat = rf"(?m)^{re.escape(label)}:.*$"
    out, n = re.subn(pat, f"{label}: {value}", text, count=1)
    if n != 1:
        raise RuntimeError(f"Could not uniquely replace {label!r}")
    return out


def fmt(x: float) -> str:
    return f"{x:.12g}"


def chirped_domains(length_um: float, start_um: float, end_um: float):
    """Return exact domain lengths for Lambda(z) linear in z."""
    a = (end_um - start_um) / length_um

    boundaries = [0.0]
    m = 1
    while True:
        if abs(a) < 1e-15:
            z = 0.5 * m * start_um
        else:
            # integral 2 dz/Lambda(z) = m
            z = (start_um / a) * (math.exp(0.5 * a * m) - 1.0)

        if z >= length_um - 1e-10:
            break
        boundaries.append(z)
        m += 1

    boundaries.append(length_um)
    domains = [
        boundaries[i + 1] - boundaries[i]
        for i in range(len(boundaries) - 1)
    ]
    return boundaries, domains


# addPulse() sequence arguments use GUI units rather than SI units.
seed_args = dict(
    frequency_thz=get("Frequency 2 (Hz)") / 1e12,
    bandwidth_thz=get("Bandwidth 2 (Hz)") / 1e12,
    sg_order=int(round(get("SG order 2"))),
    cep_pi=get("CEP 2 (rad)") / math.pi,
    delay_fs=get("Delay 2 (s)") / 1e-15,
    gdd_fs2=get("GDD 2 (s^-2)") / 1e-30,
    tod_fs3=get("TOD 2 (s^-3)") / 1e-45,
    phase_material=int(round(get("Phase material 2 index"))),
    phase_thickness_um=get("Phase material thickness 2 (mcr.)"),
    beamwaist_um=get("Beamwaist 2 (m)") / 1e-6,
    x0_um=get("x offset 2 (m)") / 1e-6,
    y0_um=get("y offset 2 (m)") / 1e-6,
    z0_um=get("z offset 2 (m)") / 1e-6,
    angle_x_deg=math.degrees(get("NC angle 2 (rad)")),
    angle_y_deg=math.degrees(get("NC angle phi 2 (rad)")),
    polarization_deg=math.degrees(get("Polarization 2 (rad)")),
    circularity=get("Circularity 2"),
    material_index=int(round(get("Material index"))),
    theta_deg=math.degrees(get("Crystal theta (rad)")),
    phi_deg=math.degrees(get("Crystal phi (rad)")),
)
seed_total_j = get("Pulse energy 2 (J)")


def add_pulse_command(energy_j: float) -> str:
    a = seed_args
    vals = [
        energy_j,
        a["frequency_thz"],
        a["bandwidth_thz"],
        a["sg_order"],
        a["cep_pi"],
        a["delay_fs"],
        a["gdd_fs2"],
        a["tod_fs3"],
        a["phase_material"],
        a["phase_thickness_um"],
        a["beamwaist_um"],
        a["x0_um"],
        a["y0_um"],
        a["z0_um"],
        a["angle_x_deg"],
        a["angle_y_deg"],
        a["polarization_deg"],
        a["circularity"],
        a["material_index"],
        a["theta_deg"],
        a["phi_deg"],
    ]
    return "addPulse(" + ",".join(fmt(float(v)) for v in vals) + ")"


def build_sequence(domains, distributed: bool):
    # Two domains make one full QPM period. A final truncated pair still counts
    # as one spatial seeding interval.
    n_periods = math.ceil(len(domains) / 2)
    e_each = seed_total_j / n_periods

    commands = ["init()"]

    for i, dz_um in enumerate(domains):
        # The first fractional seed is already created by init() through Pulse 2.
        # Re-inject at the start of every subsequent full grating period.
        if distributed and i > 0 and i % 2 == 0:
            commands.append(add_pulse_command(e_each))

        commands.append(f"nonlinear(d,d,d,{fmt(dz_um)},d)")
        commands.append("rotate(180)")

    return "".join(commands), n_periods, e_each


cases = [
    {
        "run_name": "SEEDTEST_FWD_20_24_500_SINGLE",
        "start_um": 20.0, "end_um": 24.0, "length_um": 500.0,
        "seeding": "single_entrance",
    },
    {
        "run_name": "SEEDTEST_REV_24_20_500_SINGLE",
        "start_um": 24.0, "end_um": 20.0, "length_um": 500.0,
        "seeding": "single_entrance",
    },
    {
        "run_name": "SEEDTEST_FWD_20_22_250_SINGLE",
        "start_um": 20.0, "end_um": 22.0, "length_um": 250.0,
        "seeding": "single_entrance",
    },
    {
        "run_name": "SEEDTEST_REV_24_22_250_SINGLE",
        "start_um": 24.0, "end_um": 22.0, "length_um": 250.0,
        "seeding": "single_entrance",
    },
    {
        "run_name": "SEEDTEST_FWD_20_24_500_DISTRIBUTED",
        "start_um": 20.0, "end_um": 24.0, "length_um": 500.0,
        "seeding": "distributed_per_period",
    },
    {
        "run_name": "SEEDTEST_REV_24_20_500_DISTRIBUTED",
        "start_um": 24.0, "end_um": 20.0, "length_um": 500.0,
        "seeding": "distributed_per_period",
    },
]

manifest_rows = []
geometry_rows = []
injection_rows = []

for case in cases:
    boundaries, domains = chirped_domains(
        case["length_um"], case["start_um"], case["end_um"]
    )
    distributed = case["seeding"] == "distributed_per_period"
    sequence, n_periods, e_each = build_sequence(domains, distributed)

    text = text0
    text = replace(text, "Thickness (m)", fmt(case["length_um"] * 1e-6))
    text = replace(text, "Batch mode", "0")
    text = replace(text, "Batch steps", "1")
    text = replace(text, "Batch mode 2", "0")
    text = replace(text, "Batch steps 2", "1")
    text = replace(text, "Sequence", sequence)
    text = replace(text, "Output base path", f"outputs/{case['run_name']}")

    if distributed:
        # init() supplies injection #1 at the crystal entrance. Further
        # injections are added by addPulse() in the sequence.
        text = replace(text, "Pulse energy 2 (J)", fmt(e_each))

    (root / "inputs" / f"{case['run_name']}.txt").write_text(text)

    manifest_rows.append({
        **case,
        "n_domains": len(domains),
        "n_seed_injections": n_periods if distributed else 1,
        "seed_energy_per_injection_nJ":
            e_each * 1e9 if distributed else seed_total_j * 1e9,
        "nominal_total_seed_energy_nJ": seed_total_j * 1e9,
        "input_file": f"inputs/{case['run_name']}.txt",
    })

    z = 0.0
    for i, dz in enumerate(domains):
        z0 = z
        z += dz
        geometry_rows.append({
            "run_name": case["run_name"],
            "domain_index": i,
            "z_start_um": z0,
            "z_end_um": z,
            "domain_length_um": dz,
            "local_period_at_start_um":
                case["start_um"]
                + (case["end_um"] - case["start_um"]) * z0 / case["length_um"],
        })

    if distributed:
        # One injection at the entrance, then one every two domains.
        for period_idx in range(n_periods):
            domain_idx = 2 * period_idx
            z_inj = boundaries[domain_idx] if domain_idx < len(boundaries) else case["length_um"]
            injection_rows.append({
                "run_name": case["run_name"],
                "injection_index": period_idx,
                "z_um": z_inj,
                "energy_nJ": e_each * 1e9,
                "source": "interface_init" if period_idx == 0 else "addPulse",
            })

with (root / "run_manifest.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=manifest_rows[0].keys())
    w.writeheader()
    w.writerows(manifest_rows)

with (root / "provenance" / "grating_geometry.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=geometry_rows[0].keys())
    w.writeheader()
    w.writerows(geometry_rows)

with (root / "provenance" / "seed_injection_schedule.csv").open("w", newline="") as f:
    if injection_rows:
        w = csv.DictWriter(f, fieldnames=injection_rows[0].keys())
        w.writeheader()
        w.writerows(injection_rows)

print("Generated:")
for r in manifest_rows:
    print(
        f"  {r['run_name']}: {r['length_um']:.0f} um, "
        f"{r['n_domains']} domains, {r['n_seed_injections']} seed injection(s)"
    )
print()
print("The distributed cases use one injection per FULL QPM period.")
print("The nominal injected pulse energies sum to the original 10 nJ seed.")
print("Because addPulse adds coherent fields, output field energy is not constrained")
print("to equal the arithmetic sum of those nominal injection energies.")
