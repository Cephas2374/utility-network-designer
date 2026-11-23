"""
Figures of paper P2 (from results/gombe/model.pkl and results/gombe/water_p2.pkl).

    python engine/figures_p2.py [network|n1|dma|fire|charts]
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from engine import figstyle as fs  # noqa: E402
from engine import store  # noqa: E402
from engine.config import RESULTS, ROOT  # noqa: E402
from engine.figures_p1_results import _aspect, base_map, net_lines, read  # noqa: E402

OUT = ROOT / "papers" / "P2_water" / "figures"
COL = {"looped": fs.UTIL["water"], "tree": fs.CAT[1]}  # validated pair (dataviz validator, light)
NAME = {"looped": "looped", "tree": "shortest-path forest"}


def load():
    m = store.load(RESULTS / "model.pkl")
    p = store.load(RESULTS / "water_p2.pkl")
    return m, p


def frame():
    edges, water, gombe = read("street_edges"), read("water"), read("commune_boundary")
    return gombe.total_bounds, water, edges


def map_stack(n, gap=0.62):
    """n full-width map panels stacked vertically, as in P1."""
    bnd, water, edges = frame()
    W = fs.TEXT_WIDTH_IN
    hm = W * _aspect(bnd, 60)
    H = n * hm + n * gap
    fig = plt.figure(figsize=(W, H))
    axs = [fig.add_axes([0, 1 - (i + 1) * (hm + gap) / H + gap * 0.62 / H, 1, hm / H]) for i in range(n)]
    for ax in axs:
        base_map(ax, bnd, water, edges)
    return fig, axs


def panel_label(ax, s):
    ax.text(0.0, 1.0, s, transform=ax.transAxes, fontsize=9, fontweight="bold", va="top", color=fs.INK)


def inflows(ax, N, sources):
    for s in sources:
        x, y = N.nodes[s]["xy"]
        ax.plot(x, y, marker="v", ms=6.5, color=fs.UTIL_DARK["water"], mec="white", mew=0.7, zorder=8)


def binned(values, bounds, colors):
    cmap = ListedColormap(colors)
    return cmap, BoundaryNorm(bounds, cmap.N)


def bin_legend(ax, bounds, colors, fmt, title, y=-0.01, extra=()):
    hs = [Line2D([], [], color=c, lw=2.4, label=fmt(a, b)) for a, b, c in zip(bounds[:-1], bounds[1:], colors)]
    leg = fs.legend_below(ax, list(extra) + hs, y=y, title=title, title_fontsize=6.5)
    leg._legend_box.align = "left"


# -----------------------------------------------------------------------------
def fig_network():
    m = store.load(RESULTS / "model.pkl")
    N, C, wi = m["Nw"], m["G"].graph["corridor"], m["wi"]
    fig, (a, b) = map_stack(2)
    segs, data = net_lines(N, C)
    D = np.array([d.get("diameter_mm", 50) for d in data])
    lw = np.interp(D, [50, 100, 200, 400, 600], [0.35, 0.6, 1.0, 1.6, 2.2])
    a.add_collection(LineCollection(segs, colors=fs.UTIL["water"], linewidths=lw, zorder=3))
    loop = [s for s, d in zip(segs, data) if d.get("loop_closure")]
    a.add_collection(LineCollection(loop, colors=fs.UTIL_DARK["water"], linewidths=0.9, zorder=4,
                                    linestyles=(0, (2, 1))))
    inflows(a, N, wi["sources"])
    panel_label(a, "a")
    fs.legend_below(a, [Line2D([], [], color=fs.UTIL["water"], lw=0.5, label="main 50 mm"),
                        Line2D([], [], color=fs.UTIL["water"], lw=1.6, label="400 mm"),
                        Line2D([], [], color=fs.UTIL["water"], lw=2.2, label="600 mm"),
                        Line2D([], [], color=fs.UTIL_DARK["water"], lw=0.9, ls=(0, (2, 1)), label="loop closure"),
                        Line2D([], [], marker="v", ls="", color=fs.UTIL_DARK["water"], ms=6, label=f"inflow ({wi['K']})")])
    # (b) peak-hour pressure at the vertices that carry demand
    b.add_collection(LineCollection(segs, colors="#b9cfe9", linewidths=0.4, zorder=3))
    v = [x for x in N if N.nodes[x].get("q_req", 0) > 0 or x in m["G"].graph["tap_vertex"].values()]
    xy = np.array([N.nodes[x]["xy"] for x in v])
    pr = np.array([N.nodes[x]["pressure_mwc"] for x in v])
    bounds = [20, 25, 30, 35, 40, 50, 80]
    colors = fs.SEQ_BLUE[:6][::-1]  # dark = low pressure (the critical end)
    cmap, norm = binned(pr, bounds, colors)
    o = np.argsort(-pr)
    b.scatter(xy[o, 0], xy[o, 1], c=pr[o], cmap=cmap, norm=norm, s=1.6, linewidths=0, zorder=5)
    inflows(b, N, wi["sources"])
    panel_label(b, "b")
    bin_legend(b, bounds, colors, lambda x, y: f"{x}–{y}", "pressure at the taps, peak hour (m)")
    fs.save(fig, OUT / "fig_network", furnish_axes=[a, b])


# -----------------------------------------------------------------------------
def fig_tree_path():
    """Why the Steiner tree fails: the supply path of its worst tap (K = 7) against the path of the
    same tap in the shortest-path forest on the looped design's inflows."""
    import copy

    import geopandas as gpd

    from engine import water, water_analysis as wa
    from engine.config import DATA
    from engine.dem import SRTM
    from engine.run_gombe import MAJOR
    m = store.load(RESULTS / "model.pkl")
    G, S, b, rows, packs = m["G"], m["S"], m["b"], m["rows"], m["packs"]
    C = G.graph["corridor"]
    bi = b.set_index("bldg_id")
    served = [r["bldg_id"] for r in rows]
    dem = SRTM(DATA / "S05E015.hgt")
    elev = lambda xy: float(dem.sample([xy])[0])
    boundary = read("commune_boundary").geometry.iloc[0]
    wp = gpd.read_file(DATA / "gombe_inputs.gpkg", layer="water")
    river = wp.loc[wp.area.idxmax()].geometry
    major = lambda v: any(G.edges[v, u]["road_class"] in MAJOR for u in G.neighbors(v))
    cands = water.inflow_candidates(G, boundary, river, major)
    terminals = sorted({G.graph["tap_vertex"][x] for x in served}, key=str)
    st = copy.deepcopy(packs)
    st["water"]["graph"]["bridgeless_trunk"] = False
    src7 = water.choose_inflows(G, cands, 7, elev, terminals, packs["water"]["capacity"]["source_pressure_mwc"])[0]
    Ns, _ = water.design_water(G, S, served, bi, st, src7, elev, print)
    sp = copy.deepcopy(st)
    sp["water"]["graph"]["tree"] = "spt"
    src3 = m["wi"]["sources"]
    Nf, _ = water.design_water(G, S, served, bi, sp, src3, elev, print)
    tapset = set(terminals)
    wst = wa.worst_supply_path(Ns, src7, tapset)
    w = min(tapset & set(Ns), key=lambda v: Ns.nodes[v]["pressure_mwc"])

    def best_path(N, src):
        paths = [nx.shortest_path(N, s, w, weight="length") for s in src if nx.has_path(N, s, w)]
        return max(paths, key=lambda p: N.nodes[p[0]]["head_m"] - (N.nodes[p[0]]["head_m"] - N.nodes[w]["head_m"]))

    fig, axs = map_stack(2, gap=0.2)  # one legend, under the lower panel
    for ax, N, src, lab, title in ((axs[0], Ns, src7, "a", "Steiner tree, 7 inflows"),
                                   (axs[1], Nf, src3, "b", f"shortest-path forest, {len(src3)} inflows")):
        segs, data = net_lines(N, C)
        ax.add_collection(LineCollection(segs, colors="#b9cfe9", linewidths=0.45, zorder=3))
        path = best_path(N, src)
        pe = {frozenset(e) for e in zip(path[:-1], path[1:])}
        hl = [s for s, e in zip(segs, N.edges()) if frozenset(e) in pe]
        ax.add_collection(LineCollection(hl, colors=fs.UTIL_DARK["water"], linewidths=2.0, zorder=5))
        inflows(ax, N, src)
        x, y = N.nodes[w]["xy"]
        ax.plot(x, y, marker="o", ms=6, mfc=fs.STATUS["critical"], mec="white", mew=0.8, zorder=9)
        L = sum(N.edges[e]["length"] for e in zip(path[:-1], path[1:]))
        p = N.nodes[w]["pressure_mwc"]
        panel_label(ax, lab)
        ax.text(1.0, 1.0, f"{title}: supply path {L / 1000:.1f} km, tap pressure {p:.1f} m", transform=ax.transAxes,
                ha="right", va="top", fontsize=7.5, color=fs.INK2)
    fs.legend_below(axs[1], [Line2D([], [], color="#b9cfe9", lw=1, label="mains"),
                             Line2D([], [], color=fs.UTIL_DARK["water"], lw=2, label="supply path of the tap"),
                             Line2D([], [], marker="o", ls="", mfc=fs.STATUS["critical"], mec="white", ms=6,
                                    label=f"worst tap of the Steiner tree (z = {wst['z']:.0f} m)"),
                             Line2D([], [], marker="v", ls="", color=fs.UTIL_DARK["water"], ms=6, label="inflow")])
    fs.save(fig, OUT / "fig_tree_path", furnish_axes=list(axs))


# -----------------------------------------------------------------------------
def fig_n1():
    """Demand not served when each analysed segment is closed: looped design and tree."""
    m, p = load()
    C = m["G"].graph["corridor"]
    out = p["out"]
    designs = {"looped": (m["Nw"], m["wi"]["sources"]), "tree": (p["Nt"], p["ti"]["sources"])}
    fig, axs = map_stack(2)
    bounds = [0, 0.5, 1, 2, 5, 10, 100]
    colors = fs.SEQ_BLUE[1:7]
    cmap, norm = binned(None, bounds, colors)
    for ax, (name, (N, src)), lab in zip(axs, designs.items(), "ab"):
        r = out[name]
        key = {str(v): v for v in N}
        seg_of = {key[k]: s for k, s in r["seg_of"].items() if k in key}
        uns = {x["segment"]: 100 * (1 - x["served_fraction"]) for x in r["n_minus_one"]}
        segs, data = net_lines(N, C)
        ends = [(a, b) for a, b in N.edges()]
        col, lw = [], []
        for (a, b), d in zip(ends, data):
            s = seg_of.get(a) if seg_of.get(a) == seg_of.get(b) else None
            if s in uns:
                col.append(cmap(norm(uns[s])))
                lw.append(1.4)
            else:
                col.append("#cfd8e3")
                lw.append(0.4)
        o = np.argsort(lw)
        ax.add_collection(LineCollection([segs[i] for i in o], colors=[col[i] for i in o],
                                         linewidths=[lw[i] for i in o], zorder=3))
        inflows(ax, N, src)
        panel_label(ax, lab)
        ax.text(1.0, 1.0, f"{NAME[name]} ({r['K']} inflows)", transform=ax.transAxes, ha="right", va="top",
                fontsize=7.5, color=fs.INK2)
        bin_legend(ax, bounds, colors, lambda x, y: f"{x:g}–{y:g}" if y < 100 else f">{x:g}",
                   "demand not served when the segment is closed (% of total)",
                   extra=[Line2D([], [], color="#cfd8e3", lw=1, label="not analysed (< 150 mm)")])
    fs.save(fig, OUT / "fig_n1", furnish_axes=list(axs))


# -----------------------------------------------------------------------------
def fig_dma():
    m, p = load()
    N, C = m["Nw"], m["G"].graph["corridor"]
    dm = p["dma"]
    dma = dm["dma"]
    # colour the districts by a greedy colouring of the quotient graph: neighbours never share a hue
    col_of = nx.coloring.greedy_color(dm["quotient"], strategy="largest_first")
    hues = fs.CAT[:max(col_of.values()) + 1]
    fig, (ax,) = map_stack(1)
    segs, data = net_lines(N, C)
    ends = list(N.edges())
    closed = {frozenset(e) for e in dm["closed"]}
    inlets = {frozenset(e) for e in dm["inlets"]}
    inner = [(s, hues[col_of[dma[a]]]) for s, (a, b) in zip(segs, ends) if dma[a] == dma[b]]
    ax.add_collection(LineCollection([s for s, _ in inner], colors=[c for _, c in inner], linewidths=0.9, zorder=3))
    for e, s in zip(ends, segs):
        f = frozenset(e)
        if f in closed or f in inlets:
            x, y = s.mean(axis=0)
            if f in inlets:
                ax.plot(x, y, marker="o", ms=4.5, mfc="white", mec=fs.INK, mew=1.0, zorder=7)
            else:
                ax.plot(x, y, marker="x", ms=3.2, color=fs.INK2, mew=0.8, zorder=6)
    # district numbers at the vertex nearest the centroid of each district
    for d in set(dma.values()):
        pts = np.array([N.nodes[v]["xy"] for v in N if dma[v] == d])
        c = pts[np.argmin(((pts - pts.mean(0)) ** 2).sum(1))]
        ax.text(c[0], c[1], str(d + 1), fontsize=6.5, fontweight="bold", ha="center", va="center", color=fs.INK,
                path_effects=fs.halo(), zorder=9)
    inflows(ax, N, p["sources"])
    fs.legend_below(ax, [Line2D([], [], color=fs.INK2, lw=0.9, label="mains, coloured by DMA"),
                         Line2D([], [], marker="o", ls="", mfc="white", mec=fs.INK, ms=4.5, label="open inlet (bulk meter)"),
                         Line2D([], [], marker="x", ls="", color=fs.INK2, ms=3.5, label="closed boundary valve"),
                         Line2D([], [], marker="v", ls="", color=fs.UTIL_DARK["water"], ms=6, label="inflow")])
    fs.save(fig, OUT / "fig_dma", furnish_axes=[ax])


# -----------------------------------------------------------------------------
def fig_fire():
    m, p = load()
    N, C = m["Nw"], m["G"].graph["corridor"]
    fire = p["out"]["fire"]["all"]
    fig, (ax,) = map_stack(1)
    segs, _ = net_lines(N, C)
    ax.add_collection(LineCollection(segs, colors="#b9cfe9", linewidths=0.4, zorder=3))
    xy = np.array([N.nodes[f["node"]]["xy"] for f in fire])
    q = np.array([f["capacity_ls"] for f in fire])
    bounds = [0, 15, 30, 50, 75, 100.1]
    colors = fs.SEQ_BLUE[1:6][::-1]  # dark = small capacity (the critical end)
    cmap, norm = binned(q, bounds, colors)
    o = np.argsort(-q)
    ax.scatter(xy[o, 0], xy[o, 1], c=q[o], cmap=cmap, norm=norm, s=5, linewidths=0.3, edgecolors="white", zorder=5)
    inflows(ax, N, p["sources"])
    bin_legend(ax, bounds, colors, lambda x, y: f"{x:g}–{y:g}" if y < 100 else f"$\\geq${x:g}",
               "fire capacity of the hydrant at 15 m residual (L/s)")
    fs.save(fig, OUT / "fig_fire", furnish_axes=[ax])


# -----------------------------------------------------------------------------
def fig_charts():
    """(a) ECDF of demand not served per closure, (b) sensitivity, (c) EPANET cross-check."""
    m, p = load()
    out = p["out"]
    fig, axs = plt.subplots(1, 3, figsize=(fs.TEXT_WIDTH_IN, 2.05), gridspec_kw=dict(wspace=0.55))
    ax = axs[0]
    for name in ("looped", "tree"):
        u = np.sort([100 * (1 - x["served_fraction"]) for x in out[name]["n_minus_one"]])
        y = 1 - np.arange(len(u)) / len(u)
        ax.step(np.maximum(u, 1e-3), y, where="post", color=COL[name], lw=1.6, label=NAME[name])
    ax.set_xscale("log")
    ax.set_xlim(1e-3, 100)
    ax.set_xlabel("demand not served (%)")
    ax.set_ylabel("share of closures exceeding")
    ax.legend(loc="lower left", fontsize=6.3, handlelength=1.2)
    ax.set_title("a", loc="left")
    # (b) sensitivity of the minimum tap pressure to demand and roughness (indexed to the peak, C = 120)
    ax = axs[1]
    sens = out["sensitivity"]
    dem = [s for s in sens if s["kind"] == "demand"]
    rough = [s for s in sens if s["kind"] == "C"]
    ax.plot([s["value"] for s in dem], [s["min_pressure"] for s in dem], color=fs.UTIL["water"], lw=1.6,
            marker="o", ms=3.5, label="demand factor")
    ax2x = [s["value"] / 120 for s in rough]
    ax.plot(ax2x, [s["min_pressure"] for s in rough], color=fs.INK2, lw=1.6, marker="s", ms=3.5,
            label="roughness $C/120$")
    ax.axhline(20, color=fs.STATUS["critical"], lw=0.8, ls=(0, (3, 2)))
    ax.text(1.62, 20.6, "20 m", color=fs.INK2, fontsize=6.3, ha="right")
    ax.set_xlabel("factor")
    ax.set_ylabel("minimum tap pressure (m)")
    ax.legend(loc="lower left", fontsize=6.3, handlelength=1.2)
    ax.set_title("b", loc="left")
    # (c) heads of the engine against EPANET 2.2
    ax = axs[2]
    ep = out["epanet"]
    h1, h2 = np.array(ep["head_engine"]), np.array(ep["head_epanet"])
    ax.plot(h2, h1, ls="", marker="o", ms=1.5, color=fs.UTIL["water"], alpha=0.5, mew=0)
    lo, hi = min(h1.min(), h2.min()), max(h1.max(), h2.max())
    ax.plot([lo, hi], [lo, hi], color=fs.INK2, lw=0.6)
    ax.set_xlabel("head, EPANET 2.2 (m)")
    ax.set_ylabel("head, engine (m)")
    ax.text(0.04, 0.96, f"max $|\\Delta h|$ = {1000 * ep['max_head_diff_m']:.1f} mm", transform=ax.transAxes, va="top",
            fontsize=6.3, color=fs.INK2)
    ax.set_title("c", loc="left")
    for a in axs:
        a.grid(True, axis="y")
    fs.save(fig, OUT / "fig_charts", furnish=False)


if __name__ == "__main__":
    fs.setup()
    which = sys.argv[1:] or ["network", "n1", "dma", "fire", "charts"]
    for w in which:
        print(w)
        globals()[f"fig_{w}"]()
