from pathlib import Path
import zipfile

import numpy as np
import matplotlib.pyplot as plt

from publication_style import apply_publication_style, SINGLE_COLUMN

apply_publication_style()

c = 299792458.0

ROOT = Path(__file__).resolve().parents[1]
FIXED = ROOT.parent / "2026-09_PPLN_fixed_period_GDD0"


def value(path, label):
    for line in Path(path).read_text().splitlines():
        if line.startswith(label + ":"):
            return float(line.split(":", 1)[1])
    raise RuntimeError(f"Could not find {label}")


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
                    f"{run}: could not uniquely identify spectrum"
                )

            member = matches[0]

        raw = np.frombuffer(
            z.read(member),
            dtype=np.float64
        )

    if raw.size != 3 * NF:
        raise RuntimeError(
            f"{run}: got {raw.size} doubles; expected {3*NF}"
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
    S_lambda = S_f * c / lam_m**2

    order = np.argsort(lam_nm)

    return lam_nm[order], S_lambda[order]


def integrate_band(f, S_f, min_nm, max_nm):
    fmin = c / (max_nm * 1e-9)
    fmax = c / (min_nm * 1e-9)

    m = (f >= fmin) & (f <= fmax)

    return np.trapezoid(
        S_f[m],
        f[m]
    ) * 1e9


cases = [
    (
        "Fixed 21.75 µm",
        "PPLN_FIXED_21p75um",
        FIXED / "inputs/PPLN_FIXED_21p75um.txt",
        FIXED / "outputs/PPLN_FIXED_21p75um.zip",
    ),
    (
        "Chirped 20→24 µm",
        "FIG13_CHIRP_20_24_GDD0",
        ROOT / "inputs/FIG13_CHIRP_20_24_GDD0.txt",
        ROOT / "outputs/FIG13_CHIRP_20_24_GDD0.zip",
    ),
    (
        "Chirped 24→20 µm",
        "FIG13_CHIRP_24_20_GDD0",
        ROOT / "inputs/FIG13_CHIRP_24_20_GDD0.txt",
        ROOT / "outputs/FIG13_CHIRP_24_20_GDD0.zip",
    ),
]


fig, ax = plt.subplots(figsize=SINGLE_COLUMN)

print()
print("Spectral comparison")
print("=" * 72)

for label, run, inp, out in cases:

    f, S_f = load(run, inp, out)

    lam, S_lam = to_wavelength(f, S_f)

    m = (
        (lam >= 900)
        & (lam <= 1125)
    )

    ax.plot(
        lam[m],
        S_lam[m],
        label=label,
    )

    total = integrate_band(
        f, S_f, 950, 1050
    )

    blue = integrate_band(
        f, S_f, 950, 1000
    )

    red = integrate_band(
        f, S_f, 1000, 1050
    )

    ratio = red / blue

    print(label)
    print(
        f"  E(950–1050)  = {total:.4f} nJ"
    )
    print(
        f"  E(950–1000)  = {blue:.4f} nJ"
    )
    print(
        f"  E(1000–1050) = {red:.4f} nJ"
    )
    print(
        f"  red/blue      = {ratio:.4f}"
    )
    print()


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

ax.set_xlim(900, 1125)

ax.set_xlabel("Wavelength (nm)")
ax.set_ylabel(
    "Spectral energy density (nJ/nm)"
)

ax.grid(True)
ax.legend()

fig.tight_layout()

figures = ROOT / "figures"
figures.mkdir(exist_ok=True)

fig.savefig(
    figures / "chirp_direction_comparison.pdf"
)

fig.savefig(
    figures / "chirp_direction_comparison.png"
)

plt.close(fig)

print("Saved:")
print("  figures/chirp_direction_comparison.pdf")
print("  figures/chirp_direction_comparison.png")
