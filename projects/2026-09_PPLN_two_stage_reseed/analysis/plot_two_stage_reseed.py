from pathlib import Path
import zipfile

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

# -----------------------------------------------------
# Publication style
# -----------------------------------------------------

fontsize = 7

mpl.rcParams.update({
    "figure.figsize": (3.4, 2.6),
    "figure.dpi": 300,
    "savefig.dpi": 300,

    "font.family": "sans-serif",
    "font.size": fontsize,
    "axes.labelsize": fontsize,
    "axes.titlesize": fontsize,
    "xtick.labelsize": fontsize - 1,
    "ytick.labelsize": fontsize - 1,
    "legend.fontsize": fontsize - 1,

    "lines.linewidth": 1.5,
    "lines.markersize": 4,

    "axes.linewidth": 0.8,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,

    "grid.linestyle": ":",
    "grid.linewidth": 0.5,
    "grid.alpha": 0.6,

    "legend.frameon": False,

    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.03,
})

c = 299792458.0

ROOT = Path(__file__).resolve().parents[1]

FIRST = (
    ROOT.parent
    / "2026-09_PPLN_seeding_artifact_tests"
)


# -----------------------------------------------------
# Helpers
# -----------------------------------------------------

def value(path, label):

    for line in Path(path).read_text().splitlines():

        if line.startswith(label + ":"):

            return float(
                line.split(":", 1)[1].strip()
            )

    raise RuntimeError(
        f"Could not find {label} in {path}"
    )


def load(run, input_path, zip_path):

    dt = value(input_path, "dt (s)")
    T = value(input_path, "Time span (s)")

    NT = int(
        8 * round(T / (8 * dt))
    )

    NF = NT // 2 + 1

    f = np.fft.rfftfreq(
        NT,
        dt,
    )

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
                    f"{run}: spectrum not unique"
                )

            member = matches[0]

        raw = np.frombuffer(
            z.read(member),
            dtype=np.float64,
        )

    expected = 3 * NF

    if raw.size != expected:

        raise RuntimeError(
            f"{run}: got {raw.size} doubles; "
            f"expected {expected}"
        )

    S_f = raw.reshape(
        (NF, 3),
        order="F",
    )[:, 2]

    return f, S_f


def wavelength_density(f, S_f):

    m = f > 0

    f = f[m]
    S_f = S_f[m]

    lam_m = c / f
    lam_nm = lam_m * 1e9

    # J/Hz -> pJ/nm
    S_lambda = (
        S_f
        * c
        / lam_m**2
        * 1e3
    )

    order = np.argsort(lam_nm)

    return (
        lam_nm[order],
        S_lambda[order],
    )


def integrate(f, S_f, min_nm, max_nm):

    fmin = c / (max_nm * 1e-9)
    fmax = c / (min_nm * 1e-9)

    m = (
        (f >= fmin)
        & (f <= fmax)
    )

    return (
        np.trapezoid(
            S_f[m],
            f[m],
        )
        * 1e12
    )


# -----------------------------------------------------
# Runs
# -----------------------------------------------------

cases = [

    (
        "After first 20→24 µm stage",

        "SEEDTEST_FWD_20_24_500_SINGLE",

        FIRST
        / "inputs"
        / "SEEDTEST_FWD_20_24_500_SINGLE.txt",

        FIRST
        / "outputs"
        / "SEEDTEST_FWD_20_24_500_SINGLE.zip",
    ),

    (
        "Stacked, no fresh seed",

        "STACK_NO_RESEED",

        ROOT
        / "inputs"
        / "STACK_NO_RESEED.txt",

        ROOT
        / "outputs"
        / "STACK_NO_RESEED.zip",
    ),

    (
        "Stacked, fresh midpoint seed",

        "STACK_FRESH_SEED",

        ROOT
        / "inputs"
        / "STACK_FRESH_SEED.txt",

        ROOT
        / "outputs"
        / "STACK_FRESH_SEED.zip",
    ),
]


# -----------------------------------------------------
# Plot + metrics
# -----------------------------------------------------

fig, ax = plt.subplots()

rows = []

for label, run, inp, out in cases:

    f, S_f = load(
        run,
        inp,
        out,
    )

    lam, S_lambda = wavelength_density(
        f,
        S_f,
    )

    m = (
        (lam >= 900)
        & (lam <= 1125)
    )

    ax.plot(
        lam[m],
        S_lambda[m],
        label=label,
    )

    blue = integrate(
        f,
        S_f,
        950,
        1000,
    )

    red = integrate(
        f,
        S_f,
        1000,
        1050,
    )

    rows.append({
        "case": label,
        "E_950_1000_pJ": blue,
        "E_1000_1050_pJ": red,
        "E_950_1050_pJ": blue + red,
        "red_blue_ratio": red / blue,
    })


# Target region
ax.axvspan(
    950,
    1050,
    alpha=0.12,
)

ax.axvline(
    1000,
    linestyle=":",
    linewidth=0.8,
)

ax.set_xlim(
    900,
    1125,
)

ax.set_xlabel(
    "Wavelength (nm)"
)

ax.set_ylabel(
    "Spectral energy density (pJ/nm)"
)

ax.grid(True)

ax.legend()

fig.tight_layout()


# -----------------------------------------------------
# Save
# -----------------------------------------------------

figures = ROOT / "figures"
processed = ROOT / "processed"

figures.mkdir(exist_ok=True)
processed.mkdir(exist_ok=True)

fig.savefig(
    figures
    / "two_stage_reseed_comparison.pdf"
)

fig.savefig(
    figures
    / "two_stage_reseed_comparison.png"
)

plt.close(fig)


df = pd.DataFrame(rows)

df.to_csv(
    processed
    / "two_stage_reseed_metrics.csv",
    index=False,
)


print()
print(df.to_string(index=False))

print()

no_seed = df[
    df["case"]
    == "Stacked, no fresh seed"
].iloc[0]

fresh = df[
    df["case"]
    == "Stacked, fresh midpoint seed"
].iloc[0]


print("Effect of fresh midpoint seed")
print("=" * 40)

print(
    "950–1050 nm energy ratio: "
    f"{fresh.E_950_1050_pJ / no_seed.E_950_1050_pJ:.4f}"
)

print(
    "950–1000 nm energy ratio: "
    f"{fresh.E_950_1000_pJ / no_seed.E_950_1000_pJ:.4f}"
)

print(
    "1000–1050 nm energy ratio: "
    f"{fresh.E_1000_1050_pJ / no_seed.E_1000_1050_pJ:.4f}"
)

print(
    "Red/blue without reseed: "
    f"{no_seed.red_blue_ratio:.4f}"
)

print(
    "Red/blue with fresh seed: "
    f"{fresh.red_blue_ratio:.4f}"
)

print()

print("Saved:")
print(
    "  figures/"
    "two_stage_reseed_comparison.pdf"
)
print(
    "  processed/"
    "two_stage_reseed_metrics.csv"
)
