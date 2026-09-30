from pathlib import Path
import zipfile

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from publication_style import apply_publication_style

apply_publication_style()

c = 299792458.0

ROOT = Path(__file__).resolve().parents[1]
FIXED_NEG = ROOT.parent / "2026-09_PPLN_fixed_period"
FIXED_ZERO = ROOT.parent / "2026-09_PPLN_fixed_period_GDD0"

TARGET_MIN = 950.0
TARGET_MAX = 1050.0


def value(path, label):
    for line in Path(path).read_text().splitlines():
        if line.startswith(label + ":"):
            return float(line.split(":", 1)[1])
    raise RuntimeError(label)


def load(run, input_path, zip_path):
    dt = value(input_path, "dt (s)")
    T = value(input_path, "Time span (s)")

    NT = int(8 * round(T / (8 * dt)))
    NF = NT // 2 + 1
    f = np.fft.rfftfreq(NT, dt)

    with zipfile.ZipFile(zip_path) as z:
        exact = f"{run}_spectrum.dat"

        if exact in z.namelist():
            member = exact
        else:
            matches = [
                n for n in z.namelist()
                if n.endswith("_spectrum.dat")
            ]
            if len(matches) != 1:
                raise RuntimeError(
                    f"{run}: cannot uniquely find spectrum"
                )
            member = matches[0]

        raw = np.frombuffer(
            z.read(member),
            dtype=np.float64
        )

    if raw.size != 3 * NF:
        raise RuntimeError(
            f"{run}: got {raw.size} doubles, expected {3*NF}"
        )

    S_f = raw.reshape((NF, 3), order="F")[:, 2]

    return f, S_f


def to_wavelength(f, S_f):
    m = f > 0

    f = f[m]
    S_f = S_f[m]

    lam_m = c / f
    lam_nm = lam_m * 1e9

    # J/Hz -> nJ/nm
    S_lam = S_f * c / lam_m**2

    order = np.argsort(lam_nm)

    return lam_nm[order], S_lam[order]


def integrated_energy(f, S_f):
    fmin = c / (TARGET_MAX * 1e-9)
    fmax = c / (TARGET_MIN * 1e-9)

    m = (f >= fmin) & (f <= fmax)

    return np.trapezoid(
        S_f[m], f[m]
    ) * 1e9


def metrics(f, S_f):
    lam, S = to_wavelength(f, S_f)

    m = (
        (lam >= TARGET_MIN)
        & (lam <= TARGET_MAX)
    )

    lam = lam[m]
    S = S[m]

    energy = integrated_energy(f, S_f)

    peak = np.max(S)

    equivalent_width = energy / peak

    occupancy = equivalent_width / (
        TARGET_MAX - TARGET_MIN
    )

    mean = np.mean(S)

    flatness_cv = np.std(S) / mean

    above_half = lam[S >= 0.5 * peak]

    if len(above_half) > 1:
        width_50 = (
            above_half.max()
            - above_half.min()
        )
    else:
        width_50 = 0.0

    return {
        "integrated_950_1050_nJ": energy,
        "peak_density_nJ_per_nm": peak,
        "equivalent_width_nm": equivalent_width,
        "spectral_occupancy": occupancy,
        "flatness_CV": flatness_cv,
        "width_above_50pct_peak_nm": width_50,
    }


cases = [
    {
        "label": "−1230 fs²",

        "fixed_run": "PPLN_FIXED_21p75um",
        "fixed_input":
            FIXED_NEG
            / "inputs/PPLN_FIXED_21p75um.txt",
        "fixed_zip":
            FIXED_NEG
            / "outputs/PPLN_FIXED_21p75um.zip",

        "chirped_run": "FIG13_CHIRP_20_24",
        "chirped_input":
            ROOT
            / "inputs/FIG13_CHIRP_20_24.txt",
        "chirped_zip":
            ROOT
            / "outputs/FIG13_CHIRP_20_24.zip",
    },

    {
        "label": "0 fs²",

        "fixed_run": "PPLN_FIXED_21p75um",
        "fixed_input":
            FIXED_ZERO
            / "inputs/PPLN_FIXED_21p75um.txt",
        "fixed_zip":
            FIXED_ZERO
            / "outputs/PPLN_FIXED_21p75um.zip",

        "chirped_run":
            "FIG13_CHIRP_20_24_GDD0",
        "chirped_input":
            ROOT
            / "inputs/FIG13_CHIRP_20_24_GDD0.txt",
        "chirped_zip":
            ROOT
            / "outputs/FIG13_CHIRP_20_24_GDD0.zip",
    },
]


fig, axes = plt.subplots(
    1,
    2,
    figsize=(6.8, 2.8),
    sharex=True,
)

rows = []

for ax, case in zip(axes, cases):

    ff, Sf = load(
        case["fixed_run"],
        case["fixed_input"],
        case["fixed_zip"],
    )

    fc, Sc = load(
        case["chirped_run"],
        case["chirped_input"],
        case["chirped_zip"],
    )

    lm_f, Sl_f = to_wavelength(ff, Sf)
    lm_c, Sl_c = to_wavelength(fc, Sc)

    mf = (
        (lm_f >= 900)
        & (lm_f <= 1125)
    )

    mc = (
        (lm_c >= 900)
        & (lm_c <= 1125)
    )

    ax.plot(
        lm_f[mf],
        Sl_f[mf],
        label="Fixed 21.75 µm",
    )

    ax.plot(
        lm_c[mc],
        Sl_c[mc],
        label="Chirped 20–24 µm",
    )

    ax.axvspan(
        950,
        1050,
        alpha=0.12,
    )

    ax.set_title(
        f"Pump GDD = {case['label']}"
    )

    ax.set_xlabel("Wavelength (nm)")
    ax.grid(True)

    fixed_metrics = metrics(ff, Sf)
    chirped_metrics = metrics(fc, Sc)

    rows.append({
        "pump_GDD": case["label"],
        "geometry": "Fixed 21.75 µm",
        **fixed_metrics,
    })

    rows.append({
        "pump_GDD": case["label"],
        "geometry": "Chirped 20–24 µm",
        **chirped_metrics,
    })

axes[0].set_ylabel(
    "Spectral energy density (nJ/nm)"
)

axes[1].legend()

fig.tight_layout()

figures = ROOT / "figures"
processed = ROOT / "processed"

figures.mkdir(exist_ok=True)
processed.mkdir(exist_ok=True)

fig.savefig(
    figures
    / "Figure13_fixed_vs_chirped_two_GDD.pdf"
)

fig.savefig(
    figures
    / "Figure13_fixed_vs_chirped_two_GDD.png"
)

plt.close(fig)

df = pd.DataFrame(rows)

df.to_csv(
    processed
    / "fixed_vs_chirped_two_GDD_metrics.csv",
    index=False,
)

print(df.to_string(index=False))

print()
print("Saved:")
print(
    "  figures/"
    "Figure13_fixed_vs_chirped_two_GDD.pdf"
)
print(
    "  processed/"
    "fixed_vs_chirped_two_GDD_metrics.csv"
)
