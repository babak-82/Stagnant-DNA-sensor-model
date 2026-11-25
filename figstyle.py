"""Shared figure style for the paper figures (one theme, one palette)."""
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

INK, MUTED, GRID, PDE = "#1F2937", "#6B7280", "#E5E7EB", "#111827"
STAG = {5: "#8FC1E3", 15: "#4A9CCF", 45: "#0072B2", 120: "#004A75"}
WM = {5: "#F0B58A", 15: "#E88A4C", 45: "#D55E00", 120: "#8F3F00"}
GROW, FIXED, BESTFIXED = "#009E73", "#6B7280", "#E69F00"
BETA = {0.01: "#A08BD6", 0.1: "#7C5CC4", 1.0: "#522A9E", 10.0: "#240C5E"}
RATIO = {0.25: "#6A51A3", 1.0: "#1B9E77", 4.0: "#C98A00"}
CM = 1 / 2.54
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "figures_v2")


def setup():
    available = {font.name for font in mpl.font_manager.fontManager.ttflist}
    preferred = "Inter" if "Inter" in available else "DejaVu Sans"
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": [preferred, "Arial", "Helvetica", "DejaVu Sans"],
        "mathtext.fontset": "dejavusans" if preferred == "DejaVu Sans" else "custom",
        **({"mathtext.rm": "Inter", "mathtext.it": "Inter:italic", "mathtext.bf": "Inter:bold"} if preferred == "Inter" else {}),
        "mathtext.default": "it",
        "font.size": 8, "axes.labelsize": 8.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7,
        "axes.linewidth": 0.7, "axes.edgecolor": INK, "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": INK, "ytick.color": INK, "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3, "xtick.minor.size": 1.8, "ytick.minor.size": 1.8,
        "xtick.major.width": 0.7, "ytick.major.width": 0.7, "xtick.minor.width": 0.5, "ytick.minor.width": 0.5,
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
        "grid.color": GRID, "grid.linewidth": 0.5, "axes.axisbelow": True,
        "lines.linewidth": 1.8, "legend.frameon": False,
        "figure.dpi": 150, "savefig.dpi": 600, "savefig.facecolor": "white",
        "pdf.fonttype": 42, "svg.fonttype": "none",
    })


def panel_letter(ax, s, x=-0.12, y=1.04):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=9.5, fontweight="bold", color=INK, va="bottom", ha="left")


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=600)
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    fig.savefig(os.path.join(OUT, name + ".svg"))
    fig.savefig(os.path.join(OUT, name + "_preview.png"), dpi=150)
    plt.close(fig)


DEEPSEA = LinearSegmentedColormap.from_list(
    "deepsea", ["#DEF5E5", "#60CEAC", "#3497A9", "#395D9C", "#382A54", "#0B0405"], N=256)
