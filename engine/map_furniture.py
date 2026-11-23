"""Scale bar and north arrow for the map figures (projected axes in metres, EPSG:32733).

furnish(fig) finds the map axes of a figure (equal aspect, axis switched off), renders the figure and
places a scale bar and a north arrow where they cover least of each map (corners preferred), on a white background so they
stay legible where a corner is not completely empty. Charts and colour bars are left untouched.
Sizes are scaled with the figure width so that they print at about 7 pt once the figure is shrunk to the
text width of the paper (13.9 cm).
"""
import numpy as np
from matplotlib.patches import FancyBboxPatch, Rectangle

TEXT_WIDTH_IN = 5.47
NICE = [10, 20, 25, 50, 100, 200, 250, 500, 1000, 2000, 2500, 5000]
GRID = [(fx / 8, fy / 6) for fx in range(9) for fy in range(7)]     # candidate positions, (0|1, 0|1) = corners


def is_map(ax):
    return ax.get_aspect() == 1.0 and not ax.axison and ax.get_label() != "<colorbar>"


def _box(ax_px, pos, w, h, m):
    """Pixel box (x0, y0, x1, y1) of size w*h at fractional position pos=(fx, fy) inside the axes, margin m."""
    x0, y0, x1, y1 = ax_px
    bx = x0 + m + pos[0] * (x1 - x0 - 2 * m - w)
    by = y0 + m + pos[1] * (y1 - y0 - 2 * m - h)
    return bx, by, bx + w, by + h


def _grow(box, m):
    return box[0] - m, box[1] - m, box[2] + m, box[3] + m


def _corner(pos):
    return pos[0] in (0, 1) and pos[1] in (0, 1)


def _ink(img, box, H):
    """Share of non-white pixels in a box given in display coordinates (origin bottom left)."""
    x0, y0, x1, y1 = (int(round(v)) for v in box)
    sub = img[max(H - y1, 0):max(H - y0, 0), max(x0, 0):max(x1, 0), :3]
    # darkness rather than a plain pixel count, so pale background streets weigh less than map content
    return 1.0 if sub.size == 0 else float((255 - sub.min(axis=2)).mean() / 255)


def furnish(fig, corners=None):
    """Add a scale bar and a north arrow to every map axes of fig.

    corners: optional ((fx, fy), (fx, fy)) positions of the scale bar and north arrow, to override the
    automatic choice ((0, 0) = lower left, (1, 1) = upper right).
    """
    maps = [ax for ax in fig.axes if is_map(ax)]
    if not maps:
        return
    k = fig.get_figwidth() / TEXT_WIDTH_IN          # how much the figure is shrunk in the paper
    fs, lw = 6.8 * k, 0.6 * k
    if len(maps) > 1:                                # panels share one extent: furnish the first only, smaller
        maps, fs, lw = maps[:1], fs * 0.8, lw * 0.8
    fig.canvas.draw()
    img = np.asarray(fig.canvas.buffer_rgba())
    H = img.shape[0]
    for ax in maps:
        bb = ax.get_window_extent()
        ax_px = (bb.x0, bb.y0, bb.x1, bb.y1)
        xa, xb = ax.get_xlim()
        m_per_px = (xb - xa) / bb.width
        pt = fig.dpi / 72
        # bar lengths long enough for the "0" and end labels, short enough to fit in a corner
        lengths = [n for n in NICE if n / m_per_px >= 3.2 * fs * pt and n <= 0.3 * (xb - xa)][:3] \
            or [max([n for n in NICE if n <= 0.3 * (xb - xa)] or [NICE[0]])]
        sb_size = lambda n: (n / m_per_px + 3.6 * fs * pt, 3.2 * fs * pt)
        na_w, na_h = 2.4 * fs * pt, 4.2 * fs * pt
        margin = 0.6 * fs * pt
        if corners:
            sb_c, na_c = corners
            length = lengths[0]
        else:
            # least map covered (box plus a clear margin); ties go to the shorter bar, then to a corner
            score = lambda n, c: (_ink(img, _grow(_box(ax_px, c, *sb_size(n), margin), margin), H)
                                  + 0.004 * lengths.index(n) + 0.01 * (not _corner(c)))
            length, sb_c = min(((n, c) for n in lengths for c in GRID), key=lambda nc: score(*nc))
            sb_box = _grow(_box(ax_px, sb_c, *sb_size(length), margin), margin)
            rest = [c for c in GRID if not _overlap(_box(ax_px, c, na_w, na_h, margin), sb_box)]
            na_c = min(rest, key=lambda c: _ink(img, _grow(_box(ax_px, c, na_w, na_h, margin), margin), H)
                       + 0.01 * (not _corner(c)))
        _scalebar(ax, _box(ax_px, sb_c, *sb_size(length), margin), length, fs, lw)
        _north(ax, _box(ax_px, na_c, na_w, na_h, margin), fs, lw)


def _overlap(a, b):
    return a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]


def _to_data(ax, box):
    inv = ax.transData.inverted()
    (x0, y0), (x1, y1) = inv.transform([(box[0], box[1]), (box[2], box[3])])
    return x0, y0, x1, y1


def _background(ax, x0, y0, x1, y1, lw):
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle="round,pad=0,rounding_size=0",
                                fc="white", ec="#BBBBBB", lw=lw * .6, alpha=1, zorder=50, clip_on=False))


def _scalebar(ax, box, length, fs, lw):
    x0, y0, x1, y1 = _to_data(ax, box)
    _background(ax, x0, y0, x1, y1, lw)
    h = (y1 - y0) * 0.16
    bx = x0 + ((x1 - x0) - length) / 2
    by = y0 + (y1 - y0) * 0.22
    for i in range(4):
        ax.add_patch(Rectangle((bx + i * length / 4, by), length / 4, h, fc="black" if i % 2 == 0 else "white",
                               ec="black", lw=lw, zorder=51, clip_on=False))
    unit, div = ("km", 1000) if length >= 1000 else ("m", 1)
    fmt = lambda v: f"{v / div:g}"
    for d, t in [(0, "0"), (length, f"{fmt(length)} {unit}")]:
        ax.text(bx + d, by + h + (y1 - y0) * 0.08, t, ha="center", va="bottom", fontsize=fs, zorder=52, clip_on=False)


def _north(ax, box, fs, lw):
    x0, y0, x1, y1 = _to_data(ax, box)
    _background(ax, x0, y0, x1, y1, lw)
    cx, hh = (x0 + x1) / 2, y1 - y0
    ax.annotate("", xy=(cx, y0 + hh * 0.62), xytext=(cx, y0 + hh * 0.12), zorder=52, annotation_clip=False,
                arrowprops=dict(arrowstyle="-|>,head_length=0.45,head_width=0.22", lw=lw * 1.4, color="black",
                                mutation_scale=fs * 1.6))
    ax.text(cx, y0 + hh * 0.66, "N", ha="center", va="bottom", fontsize=fs * 1.1, fontweight="bold", zorder=52,
            clip_on=False)
