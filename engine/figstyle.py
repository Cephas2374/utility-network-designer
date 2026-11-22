"""
House style for all figures of the paper series: Source Sans 3, recessive chrome,
utility identity colours (APWA uniform colour code, validated all-pairs for colour
vision deficiency) with a line pattern as the secondary encoding.
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager

FONTS = Path(__file__).resolve().parent / "fonts"
TEXT_WIDTH_IN = 5.47  # text width of the papers (13.9 cm)

INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SURFACE = "#ffffff"
STREET = "#c9c7bf"
STREET_MINOR = "#e1e0d9"
BUILDING = "#d9d7cf"
WATER_BODY = "#dbe8f5"

UTIL = {"water": "#2a78d6", "electricity": "#e34948", "gas": "#eda100"}
UTIL_DARK = {"water": "#184f95", "electricity": "#a82a2a", "gas": "#8a5d00"}  # labels in the utility hue family
DASH = {"water": "solid", "electricity": (0, (6, 2)), "gas": (0, (2, 1.5))}
LABEL = {"water": "Water", "electricity": "Electricity", "gas": "Gas"}

SEQ_BLUE = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
STATUS = {"good": "#0ca30c", "warning": "#fab219", "serious": "#ec835a", "critical": "#d03b3b"}
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]


def setup():
    for f in FONTS.glob("*.ttf"):
        font_manager.fontManager.addfont(str(f))
    mpl.rcParams.update({
        "font.family": "Source Sans 3",
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "axes.edgecolor": AXIS,
        "axes.labelcolor": INK2,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "axes.grid": False,
        "grid.color": GRID,
        "grid.linewidth": 0.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03,
        "pdf.fonttype": 42,
        "legend.frameon": False,
        "lines.solid_capstyle": "round",
    })


def map_axes(ax, bounds=None, pad=60):
    ax.set_aspect("equal")
    ax.set_axis_off()
    if bounds is not None:
        x0, y0, x1, y1 = bounds
        ax.set_xlim(x0 - pad, x1 + pad)
        ax.set_ylim(y0 - pad, y1 + pad)


def halo(color="white", w=2.0):
    import matplotlib.patheffects as pe
    return [pe.withStroke(linewidth=w, foreground=color)]


def legend_below(ax, handles, ncol=None, y=-0.02, **kw):
    """Legend in one row (or a few) under a map panel, outside the drawing."""
    return ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, y), ncol=ncol or len(handles),
                     fontsize=6.5, handlelength=1.4, columnspacing=1.2, handletextpad=0.4, borderaxespad=0, **kw)


def not_a_map(ax):
    """Exclude an equal-aspect panel (e.g. in degrees) from the scale bar / north arrow."""
    ax.set_aspect(1.0000001)


def save(fig, path, furnish=True, furnish_axes=None):
    """Save PDF + PNG.  furnish_axes: list of map axes that get their own scale bar and north arrow."""
    from .map_furniture import furnish as add_furniture, is_map
    if furnish_axes is not None:
        maps = [ax for ax in fig.axes if is_map(ax)]
        for target in furnish_axes:
            hidden = [ax for ax in maps if ax is not target]
            for ax in hidden:
                ax.set_aspect(1.0000001)
            add_furniture(fig)
            for ax in hidden:
                ax.set_aspect(1.0)
    elif furnish:
        add_furniture(fig)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".pdf"))
    fig.savefig(path.with_suffix(".png"), dpi=200)
    plt.close(fig)
