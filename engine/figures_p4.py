"""
Figures of paper P4: engineering drawings generated from the design (longitudinal profile, cross-section,
junction plan, single-line diagram), the street-cut map and the charts.

    python engine/figures_p4.py [profile|section|junction|sld|cuts|charts]
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Polygon, Rectangle
from shapely.geometry import Point

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from engine import figstyle as fs, vertical  # noqa: E402
from engine import store  # noqa: E402
from engine.config import RESULTS, ROOT, load_rules  # noqa: E402
from engine.figures_p2 import binned, bin_legend, map_stack  # noqa: E402

OUT = ROOT / "papers" / "P4_corridors" / "figures"
RES = RESULTS


def configure(results=None, out=None):
    """Read the design from another results folder and write the drawings elsewhere (QGIS plugin, P5)."""
    global RES, OUT
    RES = Path(results) if results else RES
    OUT = Path(out) if out else OUT
U3 = ("gas", "electricity", "water")
NAME = {"p1": "kerb order of P1 (gas, electricity, water)", "compatible": "compatible order (gas, water, electricity)"}
LEGACY = ["gas", "electricity", "water"]


def title(p, lay):
    """Panel title: the wording of the paper for the design of P1, else the kerb order of the layout."""
    kerb = p["keep"][lay]["kerb"]
    if lay == "compatible" or list(kerb) == LEGACY:
        return NAME[lay]
    return "design: kerb order " + ", ".join(kerb)


def layouts(p):
    return [lay for lay in ("p1", "compatible") if lay in p["keep"]]


def load():
    m = store.load(RES / "model.pkl")
    p = store.load(RES / "corridor_p4.pkl")
    m["cr"] = dict(m["cr"], vertical=load_rules()[1]["vertical"])
    return m, p


def sample_edge(m, p):
    """The street edge whose water main is lowered most in the layout of P1, among edges of 80-250 m
    that carry all three utilities on one side in both layouts."""
    S = m["S"]
    best, key = None, -1
    for line, (s, z) in p["keep"]["p1"]["prof"].items():
        _, u, e, side = line
        if u != "water":
            continue
        L = S.edges[e][3]["length"]
        if not 80 <= L <= 250:
            continue
        if not all(("main", w, e, side) in p["keep"]["p1"]["prof"] or
                   p["keep"]["p1"]["sigma"].get(w, {}).get(e) == side for w in U3):
            continue
        low = float(np.sum(z > 0.8 + 1e-6))
        if low > key:
            best, key = (e, side), low
    if best is None:  # another district: the most lowered water main of any length, else any water main
        mains = [(float(np.sum(z > 0.8 + 1e-6)), line) for line, (s, z) in p["keep"]["p1"]["prof"].items()
                 if line[1] == "water"]
        if mains:
            _, (_, _, e, side) = max(mains)
            return e, side
        e = next(iter(vertical.corridor(m).E["water"]))
        return e, p["keep"]["p1"]["sigma"]["water"][e]
    return best


# -----------------------------------------------------------------------------
def fig_profile():
    """Longitudinal profile of the water and gas mains along the sample edge, both layouts."""
    m, p = load()
    e, side = sample_edge(m, p)
    L = m["S"].edges[e][3]["length"]
    lays = layouts(p)
    fig, axs = plt.subplots(len(lays), 1, figsize=(fs.TEXT_WIDTH_IN, 1.95 * len(lays)), sharex=True,
                            gridspec_kw=dict(hspace=0.35), squeeze=False)
    axs = axs[:, 0]
    d0 = {u: m["packs"][u]["constraints"]["min_burial_depth_m"] for u in U3}
    for ax, lay in zip(axs, lays):
        k = p["keep"][lay]
        ax.axhline(0, color=fs.INK, lw=0.8)
        ax.fill_between([0, L], 0, -0.06, color=fs.BUILDING, lw=0)
        for u in ("electricity", "water", "gas"):
            line = ("main", u, e, side)
            col = fs.UTIL[u]
            if line in k["prof"]:
                s, z = k["prof"][line]
                ax.plot(s, -z, color=col, lw=1.6, label=f"{u} main")
                ax.plot([0, L], [-d0[u], -d0[u]], color=col, lw=0.6, ls=(0, (3, 2)))
                pts = k["req"].get(line, [])
                if pts:
                    xs, ds = zip(*pts)
                    ax.plot(xs, -np.array(ds), ls="", marker="v", ms=2.6, color=fs.UTIL_DARK[u], alpha=0.8)
                    mx = max(ds)
                    ax.plot([0, L], [-mx, -mx], color=col, lw=0.6, ls=(0, (1, 1.5)))
            else:
                ax.plot([0, L], [-d0[u], -d0[u]], color=col, lw=1.6, label=f"{u} main")
        ax.set_ylim(-1.75, 0.15)
        ax.set_ylabel("depth of crown (m)")
        ax.set_yticks([0, -0.5, -1.0, -1.5])
        ax.set_yticklabels(["0", "0.5", "1.0", "1.5"])
        ax.set_title(title(p, lay), loc="left", fontsize=8)
        ax.grid(True, axis="y")
    axs[-1].set_xlabel(f"position along the street (m); street edge {e}, side {'left' if side > 0 else 'right'}")
    axs[0].legend(handles=[Line2D([], [], color=fs.UTIL[u], lw=1.6, label=f"{u} main") for u in ("electricity", "water", "gas")]
                  + [Line2D([], [], color=fs.INK2, lw=0.6, ls=(0, (3, 2)), label="standard depth"),
                     Line2D([], [], color=fs.INK2, lw=0.6, ls=(0, (1, 1.5)), label="uniform lowering"),
                     Line2D([], [], color=fs.INK2, ls="", marker="v", ms=3, label="required at a crossing")],
                  loc="lower center", bbox_to_anchor=(0.5, 1.12), ncol=3, fontsize=6.3)
    fs.save(fig, OUT / "fig_profile", furnish=False)


# -----------------------------------------------------------------------------
def fig_section():
    """Cross-section through one tap of the sample edge, both layouts, with service dips."""
    m, p = load()
    e, side = sample_edge(m, p)
    cr = m["cr"]
    S = m["S"]
    cls = S.edges[e][3]["road_class"]
    rc = cr["road_classes"].get(cls, cr["road_classes"]["residential"])
    w2 = rc["width_m"] / 2
    taps = [r for r in m["rows"] if r["eid"] == e and not r.get("host") and r["side"] == side]
    # a typical facade: the tap whose building is nearest to 8 m beyond the verge
    r = min(taps, key=lambda r: abs(r["dist"] - (w2 + rc["verge_m"] + 8)))
    lays = layouts(p)
    fig, axs = plt.subplots(len(lays), 1, figsize=(fs.TEXT_WIDTH_IN, 2.3 * len(lays)), sharex=True,
                            gridspec_kw=dict(hspace=0.35), squeeze=False)
    axs = axs[:, 0]
    vr = cr["vertical"]
    Dm = vertical.main_dims(m, vr)
    for ax, lay in zip(axs, lays):
        k = p["keep"][lay]
        cor = vertical.corridor(m)  # the layout's sides and kerb order
        cor.kerb = k["kerb"]
        cor.order = {u: i for i, u in enumerate(cor.kerb)}
        cor.sigma, cor.dual = k["sigma"], k["dual"]
        prof = vertical.analyse(cor, m, m["packs"], cr)["prof"]  # services included
        # road surface, carriageway, building (the section is drawn towards the building side)
        x0 = -w2 - 1.5
        ax.add_patch(Rectangle((-w2, 0), 2 * w2, 0.15, color="#bdbdbd", lw=0))
        ax.plot([x0, r["dist"]], [0, 0], color=fs.INK, lw=0.8)
        ax.add_patch(Rectangle((r["dist"], -0.3), 1.5, 2.0, color=fs.BUILDING, lw=0))
        ax.text(r["dist"] + 0.75, 0.7, "building", rotation=90, ha="center", va="center", fontsize=6, color=fs.INK2)
        ax.text(0, 0.3, f"carriageway {2 * w2:.0f} m", ha="center", fontsize=6, color=fs.INK2)
        ax.plot([w2, w2], [0, 0.15], color=fs.INK, lw=0.8)
        ax.text(w2 + 0.1, 0.3, "kerb", fontsize=6, color=fs.INK2)
        for u in cor.U:
            if e not in cor.E[u]:
                continue
            for sd in cor.sides_of(u, e):
                x = cor.offset(u, e, sd) * (1 if sd == side else -1)
                line = ("main", u, e, sd)
                z = m["packs"][u]["constraints"]["min_burial_depth_m"]
                if line in k["prof"]:
                    s, zz = k["prof"][line]
                    z = float(np.interp(r["s"], s, zz))
                D = Dm[u][e] if u != "electricity" else vr["outer_dimension_m"]["electricity_main"]
                D = max(D, 0.05)
                ax.add_patch(Circle((x, -z - D / 2), D / 2, color=fs.UTIL[u], zorder=4))
                if sd == side:  # the service to the building, along its least profile
                    zs = m["packs"][u]["constraints"]["min_burial_depth_m"]
                    sv = ("svc", u, r["bldg_id"], e)
                    if sv in prof:
                        s_, z_ = prof[sv]
                        ax.plot(x + s_, -z_, color=fs.UTIL[u], lw=1.0, zorder=3)
                    else:
                        ax.plot([x, r["dist"]], [-zs, -zs], color=fs.UTIL[u], lw=1.0, zorder=3)
        ax.set_xlim(x0, r["dist"] + 1.8)
        ax.set_ylim(-1.7, 0.5)
        ax.set_aspect(4.0)  # depth exaggerated four times
        ax.set_title(title(p, lay), loc="left", fontsize=7.5)
        ax.set_ylabel("depth (m)")
        ax.set_yticks([0, -0.5, -1.0, -1.5])
        ax.set_yticklabels(["0", "0.5", "1.0", "1.5"])
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.grid(True, axis="y")
    axs[-1].set_xlabel("distance from the street axis (m); depth exaggerated four times")
    fs.legend_below(axs[-1], [Line2D([], [], marker="o", ls="-", lw=1, color=fs.UTIL[u], ms=5, label=f"{u} main and service")
                             for u in U3], y=-0.42)
    fs.save(fig, OUT / "fig_section", furnish=False)


# -----------------------------------------------------------------------------
def fig_junction():
    """Plan of a junction of four streets: carriageways, the three mains on their sides and the legs
    that cross the carriageways (P1, Proposition 8)."""
    m, p = load()
    cor = vertical.corridor(m)
    S = m["S"]
    cr = m["cr"]
    rec = dict(services={}, road={}, road_service={}, road_junction={}, inter={}, inter_service=0,
               inter_junction=0, c2_shifted={}, c2_unavoidable={}, legs=[])
    from collections import defaultdict
    rec = {k: (defaultdict(int) if isinstance(v, dict) else v) for k, v in rec.items()}
    rec["road"] = defaultdict(lambda: defaultdict(int))
    rec["legs"] = []
    best, nb = None, -1
    for v, hs in cor.rot.items():
        if len(hs) != 4:
            continue
        rec["legs"] = []
        cor.junction(v, rec)
        n = len({u for u, _, _ in rec["legs"]})
        lens = [S.edges[e][3]["length"] for e, _ in hs]
        if n > nb and min(lens) > 40 and all(any(e in cor.E[u] for u in cor.U) for e, _ in hs):
            best, nb, legs = v, n, list(rec["legs"])
    if best is None:  # another district: the junction of degree >= 3 with the most legs
        for v, hs in cor.rot.items():
            if len(hs) < 3:
                continue
            rec["legs"] = []
            cor.junction(v, rec)
            n = len(rec["legs"])
            if n > nb:
                best, nb, legs = v, n, list(rec["legs"])
    v = best
    e0 = cor.rot[v][0][0]
    g0 = S.edge_geom(e0)
    center = Point(g0.coords[0] if S.edges[e0][3]["start"] == v else g0.coords[-1])
    fig, ax = plt.subplots(figsize=(fs.TEXT_WIDTH_IN, 3.6))
    R = 38.0
    for e, _ in cor.rot[v]:
        g = S.edge_geom(e)
        if g.distance(center) > 1:
            continue
        g = g.intersection(center.buffer(R))
        w2 = cor.width(e) / 2
        ax.add_patch(Polygon(np.asarray(g.buffer(w2, cap_style=2).exterior.coords), color="#e6e6e6", lw=0))
        disk = center.buffer(cor.width(e) / 2 + cor.verge(e) + cr["corridor_layout"]["junction_disk_margin_m"])
        for u in cor.U:
            if e not in cor.E[u]:
                continue
            for sd in cor.sides_of(u, e):
                # the main ends where the legs of the junction disk take over
                gg = (S.edge_geom(e).offset_curve(sd * cor.offset(u, e, sd)).intersection(center.buffer(R))
                      .difference(disk))
                if gg.is_empty:
                    continue
                parts = getattr(gg, "geoms", [gg])
                for part in parts:
                    xy = np.asarray(part.coords)
                    ax.plot(xy[:, 0], xy[:, 1], color=fs.UTIL[u], lw=1.4, ls=fs.DASH[u], solid_capstyle="butt")
    margin = cr["corridor_layout"]["junction_disk_margin_m"]
    for u, vv, e in legs:
        g = S.edge_geom(e)
        Lg = g.length
        Rl = min(Lg / 2, cor.width(e) / 2 + cor.verge(e) + margin)
        at_start = S.edges[e][3]["start"] == v
        s = Rl if at_start else Lg - Rl
        a, b = np.array(g.interpolate(max(0, s - 0.5)).coords[0]), np.array(g.interpolate(min(Lg, s + 0.5)).coords[0])
        t = (b - a) / np.linalg.norm(b - a)
        nrm = np.array([-t[1], t[0]])
        c = np.array(g.interpolate(s).coords[0])
        o = cor.lane[u] + cor.width(e) / 2
        ax.plot(*np.array([c - o * nrm, c + o * nrm]).T, color=fs.UTIL[u], lw=1.6, ls=fs.DASH[u])
    ax.plot(center.x, center.y, marker="o", ms=3, color=fs.INK)
    ax.set_aspect("equal")
    ax.set_xlim(center.x - R, center.x + R)
    ax.set_ylim(center.y - R * 0.75, center.y + R * 0.75)
    ax.set_axis_off()
    ax.plot([center.x - R + 3, center.x - R + 13], [center.y - R * 0.7] * 2, color=fs.INK, lw=1.2)
    ax.text(center.x - R + 8, center.y - R * 0.7 + 1, "10 m", ha="center", fontsize=6.5)
    fs.legend_below(ax, [Line2D([], [], color=fs.UTIL[u], lw=1.4, ls=fs.DASH[u], label=f"{u} main / leg") for u in U3]
                    + [Polygon([[0, 0]], color="#e6e6e6", label="carriageway")])
    fs.save(fig, OUT / "fig_junction", furnish=False)


# -----------------------------------------------------------------------------
def fig_sld():
    """Single-line diagram of one 20 kV open ring: infeed, ring main units, transformers, open point."""
    m, p = load()
    er = m["er"]
    ring = min(er["rings"], key=lambda r: abs(len(set(r["tour"][1:-1])) - 12))
    subs = list(dict.fromkeys(ring["tour"][1:-1]))
    tr = er["transformers"]
    n = len(subs)
    k = ring["nop_leg"]
    fig, ax = plt.subplots(figsize=(fs.TEXT_WIDTH_IN, 3.0))
    xs = np.linspace(0.6, n - 0.4, n)
    y0, y1 = 2.6, 0.0
    ax.plot([-0.4, -0.4], [y0, y1 - 0.0], color=fs.INK, lw=1.2)
    ax.plot([-0.9, 0.1], [y0 + 0.35] * 2, color=fs.INK, lw=3)
    ax.text(-0.4, y0 + 0.55, "20 kV infeed", ha="center", fontsize=6.5)
    ax.plot([-0.4, -0.4], [y0 + 0.35, y0], color=fs.INK, lw=1.2)
    # feeder A along the top to the open point, feeder B along the bottom back
    half = (n + 1) // 2
    pos = []
    for i in range(n):
        if i < half:
            pos.append((xs[i] * n / half * 0.5 if False else 0.4 + i * (n - 0.8) / max(1, half - 1) * 0.98, y0))
        else:
            j = n - 1 - i
            pos.append((0.4 + j * (n - 0.8) / max(1, n - half - 1 if n - half > 1 else 1) * 0.98, y1))
    pts = [(-0.4, y0)] + pos + [(-0.4, y1)]
    for i in range(len(pts) - 1):
        a, b = pts[i], pts[i + 1]
        open_leg = i == k
        if a[1] != b[1]:
            mid = (max(a[0], b[0]) + 0.5, None)
            path = [a, (mid[0], a[1]), (mid[0], b[1]), b]
        else:
            path = [a, b]
        xy = np.array(path)
        ax.plot(xy[:, 0], xy[:, 1], color=fs.UTIL["electricity"], lw=1.2)
        if open_leg:
            c = xy[len(xy) // 2 - 1:len(xy) // 2 + 1].mean(axis=0)
            ax.plot([c[0] - 0.12, c[0] + 0.12], [c[1] + 0.18, c[1] - 0.18], color=fs.INK, lw=1.6)
            ax.text(c[0] + 0.15, c[1] + 0.2, "open point", fontsize=6, color=fs.INK)
    for s, (x, y) in zip(subs, pos):
        ax.add_patch(Rectangle((x - 0.17, y - 0.12), 0.34, 0.24, facecolor="white", edgecolor=fs.INK, lw=0.8, zorder=5))
        dy = -0.5 if y == y0 else 0.5
        ax.plot([x, x], [y + np.sign(dy) * 0.12, y + dy * 0.6], color=fs.INK, lw=0.8)
        for off in (0.0, 0.11):
            ax.add_patch(Circle((x, y + dy * (0.7 + off * 2.2)), 0.1, facecolor="none", edgecolor=fs.INK, lw=0.7, zorder=5))
        t = tr[s]
        lab = f"{t['rating_kva']:.0f}" if "mv_customer" not in t else f"MV {t['kva']:.0f}"
        ax.text(x, y + dy * 1.25, lab, ha="center", va="center", fontsize=5.5, color=fs.INK2)
    ax.set_xlim(-1.2, n + 0.8)
    ax.set_ylim(-0.4, 3.6)
    ax.set_axis_off()
    ax.text(n + 0.6, y0 + 0.55, f"ring of {n} substations, {ring['length_m'] / 1000:.1f} km; labels: kVA rating",
            ha="right", fontsize=6.5, color=fs.INK2)
    fs.save(fig, OUT / "fig_sld", furnish=False)


# -----------------------------------------------------------------------------
def fig_cuts():
    """Street edges coloured by the connections that lose two or three utilities when the street is
    dug up."""
    m, p = load()
    S = m["S"]
    fig, (ax,) = map_stack(1)
    rows = {r["eid"]: r for r in p["cuts"]}
    segs, vals = [], []
    for e, r in rows.items():
        segs.append(np.asarray(S.edge_geom(e).coords))
        vals.append(r["two"])
    vals = np.array(vals)
    bounds = [0, 1, 10, 50, 200, 2000]
    colors = ["#d9d7cf", "#cfc3ef", "#a08fe0", "#6250d6", "#3b2e8f"]
    cmap, norm = binned(vals, bounds, colors)
    o = np.argsort(vals)
    ax.add_collection(LineCollection([segs[i] for i in o], colors=[cmap(norm(vals[i])) for i in o],
                                     linewidths=[0.5 + 0.35 * np.log10(1 + vals[i]) for i in o], zorder=3))
    bin_legend(ax, bounds, colors, lambda x, y: "0" if x == 0 else f"{x}–{y - 1}",
               "connections losing at least two utilities when the street is dug up")
    fs.save(fig, OUT / "fig_cuts", furnish_axes=[ax])


# -----------------------------------------------------------------------------
def fig_charts():
    """(a) extra excavation for the lowered mains by policy and layout; (b) lowered main length
    against the largest gradient; (c) expected connections losing each utility per street cut."""
    m, p = load()
    o = p["out"]
    fig, axs = plt.subplots(1, 3, figsize=(fs.TEXT_WIDTH_IN, 2.2), gridspec_kw=dict(wspace=0.6))
    ax = axs[0]
    labs, cone, unif, stack = [], [], [], []
    for lay in ("p1", "compatible"):
        t = o[lay]["grad_0.1"]["totals"]
        labs.append("P1" if lay == "p1" else "compat.")
        cone.append(sum(t[u]["mains_extra_m3"] for u in ("water", "gas")) / 1000)
        unif.append(sum(t[u]["mains_uniform_m3"] for u in ("water", "gas")) / 1000)
        stack.append(sum(t[u]["stack_extra_m3"] for u in ("water", "gas")) / 1000)
    x = np.arange(2)
    for i, (vals, col, lab) in enumerate(((cone, fs.UTIL["water"], "least profile"), (unif, fs.CAT[1], "edge lowered"),
                                          (stack, fs.INK2, "depth stack"))):
        ax.bar(x + (i - 1) * 0.26, vals, width=0.24, color=col, label=lab)
    ax.set_xticks(x)
    ax.set_xticklabels(labs)
    ax.set_ylabel("extra excavation (1000 m³)")
    ax.legend(fontsize=5.8, loc="upper center", bbox_to_anchor=(0.5, -0.16), handlelength=1, ncol=1)
    ax.set_title("a", loc="left")
    ax = axs[1]
    for lay, col in (("p1", fs.CAT[1]), ("compatible", fs.UTIL["water"])):
        g = [0.05, 0.1, 0.2]
        ax.plot(g, [o[lay][f"grad_{x}"]["totals"]["water"]["mains_lowered_m"] / 1000 for x in g], color=col, lw=1.6,
                marker="o", ms=3.5, label=("P1" if lay == "p1" else "compatible") + ", water")
        ax.plot(g, [o[lay][f"grad_{x}"]["totals"]["gas"]["mains_lowered_m"] / 1000 for x in g], color=col, lw=1.2,
                ls=(0, (3, 2)), marker="s", ms=3, label=("P1" if lay == "p1" else "compatible") + ", gas")
    ax.set_xscale("log")
    from matplotlib.ticker import NullFormatter
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xticks([0.05, 0.1, 0.2])
    ax.set_xticklabels(["1:20", "1:10", "1:5"])
    ax.set_xlabel("largest gradient")
    ax.set_ylabel("main lowered (km)")
    ax.legend(fontsize=5.5, loc="upper right", handlelength=1.4)
    ax.set_title("b", loc="left")
    ax = axs[2]
    ex = o["street_cuts"]["expected"]
    keys = [("water", "water"), ("electricity", "electr."), ("gas", "gas"), ("two", "≥ 2"), ("three", "all 3")]
    ax.barh(range(len(keys)), [ex[k] for k, _ in keys],
            color=[fs.UTIL["water"], fs.UTIL["electricity"], fs.UTIL["gas"], fs.CAT[6], fs.CAT[6]])
    ax.set_yticks(range(len(keys)))
    ax.set_yticklabels([l for _, l in keys])
    ax.invert_yaxis()
    ax.set_xlabel("connections per cut")
    ax.set_title("c", loc="left")
    for a in axs:
        a.grid(True, axis="y" if a is not axs[2] else "x")
    fs.save(fig, OUT / "fig_charts", furnish=False)


if __name__ == "__main__":
    fs.setup()
    for w in sys.argv[1:] or ["profile", "section", "junction", "sld", "cuts", "charts"]:
        print(w)
        globals()[f"fig_{w}"]()
