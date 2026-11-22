"""
Validation - one independent check per graph rule (rules/graph_rules.py).

The checks recompute each property from the designed networks; they do not rely on how
the networks were built.  Every result row: rule id, utility, value, limit, pass.
"""
import math
from collections import defaultdict

import networkx as nx
import numpy as np
from shapely.geometry import LineString

from . import design
from .planar import crossing_pairs


def _row(rule, utility, quantity, value, limit, ok, note=""):
    return dict(rule=rule, utility=utility, quantity=quantity, value=value, limit=limit, ok=bool(ok), note=note)


def is_tap(x):
    return isinstance(x, tuple) and x[0] == "t"


def drawing(N, Gc):
    """Line drawing of a network on the street axes (split links as short segments)."""
    lines = []
    for a, b, d in N.edges(data=True):
        if d.get("eid") is not None and Gc.has_edge(design_base(a), design_base(b)):
            g = Gc_geometry(Gc, design_base(a), design_base(b))
        else:
            g = LineString([N.nodes[a]["xy"], N.nodes[b]["xy"]])
        if g.length > 0:
            lines.append(g)
    return lines


def design_base(x):
    if isinstance(x, tuple) and x[0] == "x":
        return x[1]
    return x


def Gc_geometry(Gc, a, b):
    return Gc.graph["corridor"].geometry(a, b)


def check_network(u, nets, Gc, packs, taps, sources, corridor, extra=None, exclude_terminals=()):
    """Checks G1-G12 for utility u.  nets: list of graphs of u (checked one by one: the two gas
    pressure tiers are separate plane graphs); sources: supply vertices."""
    p = packs[u]
    topo, gr = p["topology"], p["graph"]
    rows = []
    tv = Gc.graph["tap_vertex"]
    tapset = {tv[r["bldg_id"]] for r in taps if r["bldg_id"] not in exclude_terminals}
    own = [r for r in taps if not r.get("host") and r["dist"] > 0]
    hosted = [r for r in taps if r.get("host") and r["dist"] > 0]

    # G1 plane drawing of each network with the services it feeds; a shared service meets its
    # host at a tee, which is a vertex of the service tree
    for k, N in enumerate(nets):
        lines, kinds, ids = [], [], []
        for g in drawing(N, Gc):
            lines.append(g)
            kinds.append("main")
            ids.append(None)
        for r in own + hosted:
            if tv[r["bldg_id"]] in N:
                lines.append(LineString(r["path"]))
                kinds.append("service")
                ids.append(r)
        bad = []
        for i, j in crossing_pairs(lines):
            a_, b_ = ids[i], ids[j]
            if a_ is not None and b_ is not None and (
                    a_.get("host") == b_["bldg_id"] or b_.get("host") == a_["bldg_id"]
                    or (a_.get("host") and a_.get("host") == b_.get("host"))):
                continue
            bad.append((i, j))
        kinds_bad = defaultdict(int)
        for i, j in bad:
            kinds_bad["-".join(sorted((kinds[i], kinds[j])))] += 1
        n_svc = sum(1 for x in kinds if x == "service")
        rows.append(_row("G1", u, f"crossing pairs in the drawing (net {k})", len(bad), 0, not bad,
                         f"{dict(kinds_bad)}; {n_svc} services"))

    # G2 planarity (combinatorial) and the bound |E| <= 3|V| - 6
    for k, N in enumerate(nets):
        planar = nx.check_planarity(N)[0]
        V, E = N.number_of_nodes(), N.number_of_edges()
        rows.append(_row("G2", u, f"planar (net {k})", int(planar), 1, planar and (V < 3 or E <= 3 * V - 6),
                         f"|V|={V}, |E|={E}"))

    # G3 degree with services (a tap carries one service edge)
    worst, worst_src = 0, 0
    lim_src = gr.get("source_max_degree", topo["max_branching"])
    for N in nets:
        for x in N:
            deg = N.degree(x) + (1 if x in tapset else 0)
            if x in sources:
                worst_src = max(worst_src, deg)
            else:
                worst = max(worst, deg)
    rows.append(_row("G3", u, "max degree (junctions, with services)", worst, topo["max_branching"],
                     worst <= topo["max_branching"]))
    rows.append(_row("G3", u, "max degree (supply points)", worst_src, lim_src, worst_src <= lim_src))

    # G4 every terminal connected to a supply point
    H = nx.compose_all(nets) if len(nets) > 1 else nets[0]
    comp_has_src = {}
    for c in nx.connected_components(H):
        has = any(s in c for s in sources)
        comp_has_src.update({x: has for x in c})
    lost = sum(1 for t in tapset if not comp_has_src.get(t, False))
    rows.append(_row("G4", u, "terminals without supply", lost, 0, lost == 0))

    # G5 loop policy
    policy = gr["loop_policy"]
    for k, N in enumerate(nets):
        mu = N.number_of_edges() - N.number_of_nodes() + nx.number_connected_components(N)
        if policy in ("forest", "radial_with_ties"):
            rows.append(_row("G5", u, f"cyclomatic number mu (net {k}, closed)", mu, 0, mu == 0))
        else:
            rows.append(_row("G5", u, f"cyclomatic number mu (net {k})", mu, ">=0", True))

    # G8 segment length
    Lmax = max(d["length"] for N in nets for _, _, d in N.edges(data=True))
    rows.append(_row("G8", u, "longest edge (m)", round(Lmax, 1), topo["max_segment_length_m"],
                     Lmax <= topo["max_segment_length_m"] + 1e-6))

    # G9 services: perpendicular, length (avoidable violations only), distinct taps
    perp = np.mean([r["perpendicular"] for r in taps if not r.get("host")])
    rows.append(_row("G9", u, "share of perpendicular service taps", round(float(perp), 4), "-", True,
                     "the others leave the tap obliquely inside the face"))
    lim = p["routing"]["max_service_distance"]
    over = [r for r in taps if corridor.service_length(u, r) > lim + 1e-9]
    unavoid = [r for r in over if corridor.min_feasible_service_length(u, r) > lim + 1e-9]
    lens = np.array([corridor.service_length(u, r) for r in taps])
    av = corridor.avoidability_cache
    cur, low = av["current"]["service"].get(u, 0), av["minimum"]["service"].get(u, 0)
    rows.append(_row("G9", u, "services over the limit above the achievable minimum", cur - low, 0, cur == low,
                     f"limit {lim} m; {len(over)} over the limit, {low} at the minimum over all side / dual choices, "
                     f"{len(unavoid)} farther than the limit from every lane (private plot length); "
                     f"P95 {np.percentile(lens, 95):.1f} m"))
    per_edge = defaultdict(list)
    for r in taps:
        if not r.get("host"):
            per_edge[r["eid"]].append(r["s"])
    gaps = [b - a for ss in per_edge.values() for a, b in zip(sorted(ss)[:-1], sorted(ss)[1:])]
    mg = min(gaps) if gaps else math.inf
    rows.append(_row("G9", u, "min tap spacing (m)", round(mg, 3), gr["min_tap_spacing_m"],
                     mg >= gr["min_tap_spacing_m"] - 1e-6,
                     f"{len(hosted)} shared services use the tap of a neighbour"))

    # G12 mains on street corridors
    off = sum(1 for N in nets for _, _, d in N.edges(data=True) if d.get("eid") is None and not d.get("split_link"))
    rows.append(_row("G12", u, "main edges off the street graph", off, 0, off == 0))
    if extra:
        rows.extend(extra)
    return rows


def check_spacing(u, N, spacing, points, rule="G11", what="components"):
    """Largest along-run distance from any point of the mains to the nearest component."""
    runs = design.runs(N)
    by_run = defaultdict(list)
    for ri, s, e, off in points:
        by_run[ri].append(s)
    worst = 0.0
    for ri, path in enumerate(runs):
        L = sum(N.edges[a, b]["length"] for a, b in zip(path[:-1], path[1:]))
        ss = sorted(by_run.get(ri, []))
        if not ss:
            continue
        g = max([ss[0]] + [(b - a) / 2 for a, b in zip(ss[:-1], ss[1:])] + [L - ss[-1]])
        worst = max(worst, g)
    return _row(rule, u, f"max distance to nearest {what} along a run (m)", round(worst, 1), spacing / 2,
                worst <= spacing / 2 + 1e-6)
