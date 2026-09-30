"""Publication plotting style for Jack's Laser Physics PhD figures."""
import matplotlib as mpl

FONTSIZE = 7
SINGLE_COLUMN = (3.4, 2.6)
WIDE_FIGSIZE = (6.8, 3.6)

def apply_publication_style():
    mpl.rcParams.update({
        "figure.figsize": SINGLE_COLUMN,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "font.family": "sans-serif",
        "font.size": FONTSIZE,
        "axes.labelsize": FONTSIZE,
        "axes.titlesize": FONTSIZE,
        "xtick.labelsize": FONTSIZE - 1,
        "ytick.labelsize": FONTSIZE - 1,
        "legend.fontsize": FONTSIZE - 1,
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
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })
