"""
Result figures of paper P1 (from results/gombe/model.pkl).

    python engine/figures_p1_results.py [networks|rac|methods]
"""
import sys
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from shapely.geometry import LineString, Point, box

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from engine import figstyle as fs  # noqa: E402
from engine import store  # noqa: E402
from engine import validate  # noqa: E402
from engine.config import DATA, RESULTS, ROOT  # noqa: E402
from engine.crossings import Corridor  # noqa: E402

OUT = ROOT / "papers" / "P1_graph_formulation" / "figures"
F = DATA / "gombe_inputs.gpkg"


def load():
    return store.load(RESULTS / "model.pkl")


def read(layer):
    return gpd.read_file(F, layer=layer)


def net_lines(N, C):
    segs, data = [], []
    for a, b, d in N.edges(data=True):
        ba, bb = validate.design_base(a), validate.design_base(b)
        if d.get("eid") is not None and C.G.has_edge(ba, bb):
            g = C.geometry(ba, bb)
        else:
            g = LineString([N.nodes[a]["xy"], N.nodes[b]["xy"]])
        segs.append(np.asarray(g.coords))
        data.append(d)
    return segs, data


def base_map(ax, bnd, water, edges):
    from engine.figures_p1 import river_patch
    river_patch(ax, water, bnd)
    ax.add_collection(LineCollection([np.asarray(g.coords) for g in edges.geometry], colors="#e4e2db",
                                     linewidths=0.35, zorder=1))
    fs.map_axes(ax, bnd, pad=60)


def _aspect(b, pad):
    return (b[3] - b[1] + 2 * pad) / (b[2] - b[0] + 2 * pad)


# -----------------------------------------------------------------------------
def fig_networks():
    m = load()
    G, C = m["G"], m["G"].graph["corridor"]
    edges, water, gombe = read("street_edges"), read("water"), read("commune_boundary")
    bnd = gombe.total_bounds
    W = fs.TEXT_WIDTH_IN
    hm = W * _aspect(bnd, 60)
    gap = 0.62
    H = 3 * hm + 3 * gap
    fig = plt.figure(figsize=(W, H))
    axs = [fig.add_axes([0, 1 - (i + 1) * (hm + gap) / H + gap * 0.62 / H, 1, hm / H]) for i in range(3)]

    # (a) water: mains by diameter, inflows, loops
    ax = axs[0]
    base_map(ax, bnd, water, edges)
    segs, data = net_lines(m["Nw"], C)
    D = np.array([d.get("diameter_mm", 50) for d in data])
    lw = np.interp(D, [50, 100, 200, 400, 600], [0.35, 0.6, 1.0, 1.6, 2.2])
    ax.add_collection(LineCollection(segs, colors=fs.UTIL["water"], linewidths=lw, zorder=3))
    loop = [s for s, d in zip(segs, data) if d.get("loop_closure")]
    ax.add_collection(LineCollection(loop, colors=fs.UTIL_DARK["water"], linewidths=0.9, zorder=4,
                                     linestyles=(0, (2, 1))))
    for i, s in enumerate(m["wi"]["sources"]):
        x, y = G.nodes[s]["xy"]
        ax.plot(x, y, marker="v", ms=6.5, color=fs.UTIL_DARK["water"], mec="white", mew=0.7, zorder=6)
    wi = m["wi"]
    fs.legend_below(ax, [Line2D([], [], color=fs.UTIL["water"], lw=0.5, label="main 50 mm"),
                         Line2D([], [], color=fs.UTIL["water"], lw=1.6, label="400 mm"),
                         Line2D([], [], color=fs.UTIL["water"], lw=2.2, label="600 mm"),
                         Line2D([], [], color=fs.UTIL_DARK["water"], lw=0.9, ls=(0, (2, 1)), label=f"loop closure ({len(wi['augmentation']['faces'])} faces)"),
                         Line2D([], [], marker="v", ls="", color=fs.UTIL_DARK["water"], ms=6, label=f"inflow ({wi['K']})")],
                    y=-0.01)

    # (b) electricity: LV forest, transformers, MV rings
    ax = axs[1]
    base_map(ax, bnd, water, edges)
    segs, data = net_lines(m["Fe"], C)
    lv_tint = "#f2a5a4"  # lighter step of the electricity hue, so the 20 kV rings stand out
    ax.add_collection(LineCollection(segs, colors=lv_tint, linewidths=0.5, zorder=3))
    er = m["er"]
    mv = []
    for r in er["rings"]:
        for leg in r["legs"]:
            pts = [G.nodes[x]["xy"] for x in leg]
            if len(pts) > 1:
                mv.append(np.asarray(pts))
    ax.add_collection(LineCollection(mv, colors=fs.UTIL_DARK["electricity"], linewidths=1.0, zorder=4,
                                     linestyles=fs.DASH["electricity"]))
    tr = np.array([G.nodes[s]["xy"] for s, t in er["transformers"].items() if not t.get("mv_customer")])
    mvc = np.array([G.nodes[s]["xy"] for s, t in er["transformers"].items() if t.get("mv_customer")])
    ax.scatter(tr[:, 0], tr[:, 1], s=6, marker="s", color=fs.UTIL_DARK["electricity"], linewidths=0, zorder=6)
    if len(mvc):
        ax.scatter(mvc[:, 0], mvc[:, 1], s=14, marker="D", facecolor="white", edgecolor=fs.UTIL_DARK["electricity"],
                   linewidths=0.8, zorder=6)
    x, y = G.nodes[m["infeed"]]["xy"]
    ax.plot(x, y, marker="*", ms=9, color=fs.INK, mec="white", mew=0.6, zorder=7)
    i = er["info"]
    fs.legend_below(ax, [Line2D([], [], color=lv_tint, lw=1.0, label="LV feeder"),
                         Line2D([], [], color=fs.UTIL_DARK["electricity"], lw=0.8, ls=fs.DASH["electricity"],
                                label=f"20 kV ring ({i['mv_rings']})"),
                         Line2D([], [], marker="s", ls="", color=fs.UTIL_DARK["electricity"], ms=3.5,
                                label=f"transformer ({i['n_transformers']})"),
                         Line2D([], [], marker="D", ls="", mfc="white", mec=fs.UTIL_DARK["electricity"], ms=4,
                                label=f"MV customer ({i['n_mv_customers']})"),
                         Line2D([], [], marker="*", ls="", color=fs.INK, ms=7, label="MV infeed")], y=-0.01)

    # (c) gas: LP zones, MP backbone, regulators
    ax = axs[2]
    base_map(ax, bnd, water, edges)
    segs, data = net_lines(m["LP"], C)
    ax.add_collection(LineCollection(segs, colors=fs.UTIL["gas"], linewidths=0.45, zorder=3))
    segs, data = net_lines(m["MP"], C)
    ax.add_collection(LineCollection(segs, colors=fs.UTIL_DARK["gas"], linewidths=1.3, zorder=4,
                                     linestyles=fs.DASH["gas"]))
    reg = np.array([G.nodes[s]["xy"] for s in m["gr"]["sites"]])
    ax.scatter(reg[:, 0], reg[:, 1], s=20, marker="^", color=fs.UTIL_DARK["gas"], edgecolor="white", linewidths=0.6,
               zorder=6)
    x, y = G.nodes[m["gate"]]["xy"]
    ax.plot(x, y, marker="*", ms=9, color=fs.INK, mec="white", mew=0.6, zorder=7)
    gi = m["gr"]["info"]
    fs.legend_below(ax, [Line2D([], [], color=fs.UTIL["gas"], lw=0.8, label="LP 100 mbar"),
                         Line2D([], [], color=fs.UTIL_DARK["gas"], lw=1.3, ls=fs.DASH["gas"], label="MP 4 bar"),
                         Line2D([], [], marker="^", ls="", color=fs.UTIL_DARK["gas"], ms=5,
                                label=f"district regulator ({gi['n_regulators']})"),
                         Line2D([], [], marker="*", ls="", color=fs.INK, ms=7, label="city gate (hypothetical)")],
                    y=-0.01)
    from engine.figures_p1 import _title
    for ax, t in zip(axs, ("a  Water: mains by diameter, loops closed on street faces",
                           "b  Electricity: LV shortest-path forests of the transformer cells and 20 kV open rings",
                           "c  Gas (prospective): LP trees of the regulator cells and MP backbone")):
        _title(fig, ax, t)
    fs.save(fig, OUT / "fig_networks", furnish_axes=[axs[0]])


def fig_rac():
    """Right-angle-crossing drawing of the three networks around one junction (lanes to scale)."""
    m = load()
    G, S, rows = m["G"], m["S"], m["rows"]
    nets = {"gas": [m["LP"], m["MP"]], "electricity": [m["Fe"]], "water": [m["Nw"]]}
    cor = Corridor(S, G, nets, rows, m["packs"], m["cr"])
    cor.sigma, cor.dual = m["sigma"], m["dual"]
    # junction with degree 4 where all three utilities are present and many services
    cand = []
    for v, hs in cor.rot.items():
        if len(hs) != 4:
            continue
        es = [e for e, _ in hs]
        if all(any(e in cor.E[u] for e in es) for u in cor.U):
            n = sum(len(cor.taps.get(e, [])) for e in es)
            if all(S.edges[e][3]["length"] > 60 for e in es):
                cand.append((n, v))
    v = max(cand)[1]
    cx, cy = S.pos[v]
    r = 42
    win = (cx - r, cy - r, cx + r, cy + r)
    wb = box(*win)
    fb = read("buildings")
    W = fs.TEXT_WIDTH_IN
    fig = plt.figure(figsize=(W, W * 0.78))
    ax = fig.add_axes([0, 0.1, 0.72, 0.88])
    axl = fig.add_axes([0.74, 0.1, 0.26, 0.88])
    axl.set_axis_off()
    fb[fb.intersects(wb)].plot(ax=ax, color="#e7e5de", edgecolor="#c9c7bf", linewidth=0.4, zorder=1)
    # carriageways
    for e in S.edges:
        g = S.edge_geom(e)
        if g.intersects(wb.buffer(30)):
            w = cor.width(e)
            poly = g.buffer(w / 2, cap_style="flat")
            gpd.GeoSeries([poly]).plot(ax=ax, color="#f1f0ec", edgecolor="#b9b7ae", linewidth=0.5, zorder=2)
            ax.plot(*np.asarray(g.coords).T, color="#c9c7bf", lw=0.5, ls=(0, (6, 4)), zorder=3)
    # mains offset into their lanes, trimmed at the junction disk
    disk = Point(cx, cy).buffer(max(cor.width(e) for e, _ in cor.rot[v]) / 2 + 8.0)
    for u in cor.U:
        for e in cor.E[u]:
            g = S.edge_geom(e)
            if not g.intersects(wb.buffer(60)):
                continue
            for side in cor.sides_of(u, e):
                off = cor.offset(u, e, side)
                oc = g.offset_curve(off * side)
                oc = oc.difference(disk) if oc.intersects(disk) else oc
                for part in getattr(oc, "geoms", [oc]):
                    if not part.is_empty:
                        ax.plot(*np.asarray(part.coords).T, color=fs.UTIL[u], lw=1.4, ls=fs.DASH[u],
                                solid_capstyle="butt", zorder=5)
    # junction connections: ports on the disk, joined along the arc that avoids the longest
    # unoccupied run; a carriageway is crossed by a leg from its right-side port to its left-side
    # port at the same distance from v, hence perpendicular to the street
    hs = cor.rot[v]
    kk = len(hs)
    Rd = max(cor.width(e) for e, _ in hs) / 2 + 8.0  # ports outside every carriageway of the junction

    def port(e, from_start, out_side, u):
        g = S.edge_geom(e)
        s_ = min(Rd, g.length / 2)
        s_ = s_ if from_start else g.length - s_
        p = np.array(g.interpolate(s_).coords[0])
        t = np.subtract(g.interpolate(min(s_ + 0.5, g.length)).coords[0], g.interpolate(max(s_ - 0.5, 0)).coords[0])
        t = t / np.linalg.norm(t) * (1 if from_start else -1)  # outward direction
        n = np.array([-t[1], t[0]])
        stored = out_side if from_start else -out_side
        return dict(e=e, side=out_side, xy=p + n * out_side * cor.offset(u, e, stored), t=t)

    def corner(r1, r2):
        """Meeting point of the two lane lines (through the ports, along the streets)."""
        A = np.array([r1["t"], -r2["t"]]).T
        if abs(np.linalg.det(A)) < 1e-6:
            return None
        lam = np.linalg.solve(A, r2["xy"] - r1["xy"])
        c = r1["xy"] + lam[0] * r1["t"]
        return c if np.linalg.norm(c - np.array([cx, cy])) < 1.6 * Rd else None

    for u in cor.U:
        occ = cor.occupied(u, v)
        if len(occ) < 2:
            continue
        best, start = -1, None
        for i in sorted(occ):
            g = 0
            while (i + 1 + g) % kk not in occ:
                g += 1
            if g > best:
                best, start = g, i
        first = (start + best + 1) % kk
        recs = []
        for step in range(kk - best):
            sec = (first + step) % kk
            e_a, fs_a = hs[sec]
            e_b, fs_b = hs[(sec + 1) % kk]
            last = step == kk - best - 1
            has_a = e_a in cor.port[u].get(v, ()) and any((x if fs_a else -x) > 0 for x in cor.sides_of(u, e_a))
            has_b = e_b in cor.port[u].get(v, ()) and any((x if fs_b else -x) < 0 for x in cor.sides_of(u, e_b))
            if has_a or step > 0:
                recs.append(port(e_a, fs_a, +1, u))  # left of h_sec (a port, or the end of a leg)
            if has_b or not last:
                recs.append(port(e_b, fs_b, -1, u))  # right of h_{sec+1} (a port, or the start of a leg)
        pts = []
        for r1, r2 in zip(recs[:-1], recs[1:]):
            pts.append(r1["xy"])
            if r1["e"] != r2["e"]:  # corner inside a sector
                c = corner(r1, r2)
                if c is not None:
                    pts.append(c)
        pts.append(recs[-1]["xy"])
        ax.plot(*np.array(pts).T, color=fs.UTIL[u], lw=1.4, ls=fs.DASH[u], zorder=5)

    # services in the window: perpendicular from the main lane to the facade
    for rr in rows:
        if rr.get("host") or rr["dist"] <= 0:
            continue
        T = Point(rr["path"][0])
        if not wb.contains(T):
            continue
        e = rr["eid"]
        line = S.edge_geom(e)
        tvec = np.subtract(line.interpolate(min(rr["s"] + 0.5, line.length)).coords[0],
                           line.interpolate(max(rr["s"] - 0.5, 0)).coords[0])
        tvec = tvec / np.linalg.norm(tvec)
        n = np.array([-tvec[1], tvec[0]])
        for k, u in enumerate(cor.U):
            if e not in cor.E[u]:
                continue
            su = rr["side"] if rr["side"] in cor.sides_of(u, e) else -rr["side"]
            off = cor.offset(u, e, su)
            start = np.array(T.coords[0]) + n * off * su
            end = np.array(rr["facade"]) + tvec * (k - 1) * 0.7
            col = fs.UTIL[u]
            ax.plot([start[0], end[0]], [start[1], end[1]], color=col, lw=0.6, zorder=6, alpha=0.9)
            if su != rr["side"]:
                mid = np.array(T.coords[0])
                ax.plot(*mid, marker="o", ms=3.0, mfc="none", mec=fs.INK, mew=0.6, zorder=8)
    fs.map_axes(ax, win, pad=0)
    # legend panel: lanes and rules
    y = 0.95
    for u in ("gas", "electricity", "water"):
        axl.plot([0.02, 0.22], [y, y], color=fs.UTIL[u], lw=1.6, ls=fs.DASH[u], transform=axl.transAxes)
        axl.text(0.27, y, fs.LABEL[u] + " main", fontsize=7, va="center", transform=axl.transAxes)
        y -= 0.07
    axl.plot([0.02, 0.22], [y, y], color=fs.INK2, lw=0.6, transform=axl.transAxes)
    axl.text(0.27, y, "service (perpendicular)", fontsize=7, va="center", transform=axl.transAxes)
    y -= 0.07
    axl.plot([0.12], [y], marker="o", ms=4, mfc="none", mec=fs.INK, transform=axl.transAxes)
    axl.text(0.27, y, "service crossing the road", fontsize=7, va="center", transform=axl.transAxes)
    y -= 0.11
    txt = ("Lanes from the kerb outwards:\ngas, electricity, water\n(offsets: kerb 0.5 m, then the\n"
           "separations 0.5 m, 0.5 m)\n\nEvery crossing is perpendicular:\nthe drawing is RAC (rule C1).")
    axl.text(0.02, y, txt, fontsize=6.6, va="top", color=fs.INK2, transform=axl.transAxes, linespacing=1.35)
    from engine.figures_p1 import _title
    _title(fig, ax, "Junction of the three networks drawn in their lanes")
    fs.save(fig, OUT / "fig_rac", furnish_axes=[ax])


def fig_methods():
    import json
    met = json.loads((RESULTS / "metrics.json").read_text())
    m = load()
    fig, axs = plt.subplots(1, 3, figsize=(fs.TEXT_WIDTH_IN, 2.0), gridspec_kw=dict(wspace=0.55))
    # (a) total length of the three tree methods on the water terminal set
    ax = axs[0]
    keys = [("steiner", "Steiner"), ("mst", "MST (pruned)"), ("spt", "shortest-path")]
    L = [met["methods"][k]["length_m"] / 1000 for k, _ in keys]
    ax.barh(range(3), L, color=fs.CAT[0], height=0.55)
    ax.set_yticks(range(3))
    ax.set_yticklabels([n for _, n in keys], fontsize=7)
    ax.invert_yaxis()
    for i, x in enumerate(L):
        ax.text(x + 1, i, f"{x:.1f}", va="center", fontsize=6.5, color=fs.INK)
    ax.set_xlabel("main length (km)")
    ax.set_xlim(0, max(L) * 1.22)
    ax.set_title("a  Tree methods (water taps)", fontsize=8)
    # (b) degree distribution of the mains before and after splitting (water)
    ax = axs[1]
    deg = np.bincount([d for _, d in m["Nw"].degree()], minlength=6)[1:5]
    ax.bar(np.arange(1, 5), np.where(deg > 0, deg, np.nan), color=fs.CAT[0], width=0.6)
    ax.set_yscale("log")
    ax.set_ylim(100, deg.max() * 3)
    for k, n in zip(range(1, 5), deg):
        y = n * 1.15 if n > 0 else 130  # a zero bar is labelled at the baseline
        ax.text(k, y, f"{n:,}", ha="center", va="bottom", fontsize=6.3)
    ax.set_xticks(range(1, 5))
    ax.set_xlabel("degree in the water mains")
    ax.set_title("b  T-junctions after splitting", fontsize=8)
    # (c) cyclomatic number per network
    ax = axs[2]
    mus = [met["water"]["mu"], met["electricity"]["mu_closed"], met["electricity"]["mu_with_ties"],
           met["gas"]["mu_lp"] + met["gas"]["mu_mp"]]
    labs = ["water", "LV (ties open)", "LV (ties closed)", "gas"]
    cols = [fs.UTIL["water"], fs.UTIL["electricity"], fs.UTIL_DARK["electricity"], fs.UTIL["gas"]]
    ax.barh(range(4), mus, color=cols, height=0.55)
    ax.set_yticks(range(4))
    ax.set_yticklabels(labs, fontsize=7)
    ax.invert_yaxis()
    for i, x in enumerate(mus):
        ax.text(x + 3, i, str(x), va="center", fontsize=6.5)
    ax.set_xlim(0, max(mus) * 1.3 + 5)
    ax.set_xlabel(r"cyclomatic number $\mu$")
    ax.set_title("c  Loops = bounded faces", fontsize=8)
    fs.save(fig, OUT / "fig_methods", furnish=False)


FIGS = dict(networks=fig_networks, rac=fig_rac, methods=fig_methods)

if __name__ == "__main__":
    fs.setup()
    for n in sys.argv[1:] or list(FIGS):
        print("figure", n)
        FIGS[n]()
