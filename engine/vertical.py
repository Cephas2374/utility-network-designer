"""
Vertical design of the street corridor (paper P4): every crossing of a service or junction leg with the
main of another utility, the depth it requires, and the least profiles of mains and services.

Depths are to the crown.  At a crossing the utility listed earlier in the crossing depth order passes
below: its crown must lie at least (depth of the other + its outer dimension + clearance).  With a total
depth order the requirements form a DAG (an element is only constrained by shallower ones), so the
utilities are processed from the shallowest: each main or service gets the least profile that meets its
point requirements with a bounded gradient, the upper envelope of cones (Proposition in P4).
"""
from collections import defaultdict

import numpy as np

from .config import vertical_clearance
from .crossings import Corridor


def corridor(m, kerb_order=None, log=print):
    """The corridor of the design (sides of P1), or a new layout with another kerb order, optimised
    like P1 (side assignment and joint repair)."""
    nets = {"gas": [m["LP"], m["MP"]], "electricity": [m["Fe"]], "water": [m["Nw"]]}
    cr = dict(m["cr"])
    if kerb_order is None:
        cor = Corridor(m["S"], m["G"], nets, m["rows"], m["packs"], cr)
        cor.sigma, cor.dual = m["sigma"], m["dual"]
        return cor
    cr["corridor_layout"] = dict(cr["corridor_layout"], kerb_order=list(kerb_order))
    cor = Corridor(m["S"], m["G"], nets, m["rows"], m["packs"], cr)
    cor.optimise(log=log)
    cor.joint_repair(log=log)
    return cor


def main_dims(m, vr):
    """Outer dimension of the main of each utility on each street edge (largest pipe on the edge)."""
    D = {"electricity": defaultdict(lambda: vr["outer_dimension_m"]["electricity_main"])}
    for u, graphs in (("water", [m["Nw"]]), ("gas", [m["LP"], m["MP"]])):
        d = defaultdict(float)
        for N in graphs:
            for a, b, x in N.edges(data=True):
                if x.get("eid") is not None:
                    d[x["eid"]] = max(d[x["eid"]], x.get("diameter_mm", 0) / 1000)
        D[u] = d
    return D


def envelope(L, d0, pts, grad, step=0.25):
    """Least profile on [0, L] with z >= d0, z(s_i) >= d_i and |z'| <= grad: max(d0, d_i - grad |s - s_i|).
    Returns (grid, z)."""
    n = max(2, int(np.ceil(L / step)) + 1)
    s = np.linspace(0.0, L, n)
    z = np.full(n, d0)
    for si, di in pts:
        if di > d0:
            z = np.maximum(z, di - grad * np.abs(s - si))
    return s, z


def stats(s, z, d0, width):
    if len(s) < 2:
        return 0.0, 0.0, 0.0
    ds = np.diff(s)
    ex = np.maximum(z - d0, 0.0)
    area = float(np.sum((ex[:-1] + ex[1:]) / 2 * ds))
    low = float(np.sum(ds[(ex[:-1] > 1e-6) | (ex[1:] > 1e-6)]))
    return area * width, low, float(ex.max())


def analyse(cor, m, packs, cr):
    """All crossings of the corridor with their requirements, and the least profiles.  Returns totals
    per utility and policy, and the profiles of the mains (for drawings)."""
    vr = cr["vertical"]
    grad = vr["max_gradient"]
    order = cr["crossing_depth_order"]  # first = deepest at a crossing
    rank = {u: i for i, u in enumerate(order)}  # smaller rank passes below
    d_std = {u: packs[u]["constraints"]["min_burial_depth_m"] for u in cor.U}
    road_cover = cr["road_crossing_cover_m"]
    Dm = main_dims(m, vr)
    Ds = {u: vr["outer_dimension_m"][f"{u}_service"] for u in cor.U}
    margin = cr["corridor_layout"]["junction_disk_margin_m"]
    S = cor.S

    # -- crossing events -------------------------------------------------------------------
    # a "line" is either a main (u, e, side) or a service (u, bldg) or a leg (u, v, e)
    events = []  # (upper line, lower line, position on upper, position on lower)
    for e in set(cor.taps) | set().union(*cor.E.values()):
        for r in cor.taps.get(e, []):
            sb = r["side"]
            for u in cor.U:
                if e not in cor.E[u]:
                    continue
                su = sb if sb in cor.sides_of(u, e) else -sb
                road = su != sb
                ou = cor.offset(u, e, su)
                for w in cor.U:
                    if w == u:
                        continue
                    for sw in cor.sides_of(w, e):
                        if road:
                            hit = (sw == su and cor.order[w] < cor.order[u]) or sw == sb
                        else:
                            hit = sw == su and cor.order[w] > cor.order[u] and cor.offset(w, e, sw) < r["dist"]
                        if not hit:
                            continue
                        ow = cor.offset(w, e, sw)
                        x = abs(ou - ow) if sw == su else ou + ow  # position along the service
                        events.append((("svc", u, r["bldg_id"], e), ("main", w, e, sw), x, r["s"],
                                       "road" if road else "same"))
    rec = dict(services=defaultdict(int), road=defaultdict(lambda: defaultdict(int)), road_service=defaultdict(int),
               road_junction=defaultdict(int), inter=defaultdict(int), inter_service=0, inter_junction=0,
               c2_shifted=defaultdict(int), c2_unavoidable=defaultdict(int), legs=[])
    for v in cor.rot:
        cor.junction(v, rec)
    for u, v, e in rec["legs"]:
        L = S.edges[e][3]["length"]
        R = min(L / 2, cor.width(e) / 2 + cor.verge(e) + margin)
        s = R if S.edges[e][3]["start"] == v else L - R
        for w in cor.U:
            if w == u or cor.order[w] >= cor.order[u] or e not in cor.port[w].get(v, ()):
                continue
            for sw in cor.sides_of(w, e):
                events.append((("leg", u, v, e), ("main", w, e, sw), 0.0, s, "leg"))

    # -- requirements, from the shallowest utility down ------------------------------------------
    def depth0(line):
        kind, u = line[0], line[1]
        return road_cover[u] if kind == "leg" else d_std[u]

    def dim(line):
        kind, u = line[0], line[1]
        return Dm[u][line[2]] if kind == "main" else Ds[u]

    def length(line):
        kind, u = line[0], line[1]
        if kind == "main":
            return S.edges[line[2]][3]["length"]
        if kind == "svc":
            return max(0.5, cor.service_length(u, cor.byid[line[2]]))
        return cor.width(line[3])  # a leg spans the carriageway

    req = defaultdict(list)  # line -> [(position, required crown depth)]
    prof = {}  # line -> (grid, z)
    by_u = defaultdict(list)
    for ev in events:
        by_u[ev[0][1]].append(ev)
        by_u[ev[1][1]].append(ev)

    def final_depth(line, pos):
        if line not in prof:
            return depth0(line)
        s, z = prof[line]
        return float(np.interp(pos, s, z))

    counts = defaultdict(lambda: defaultdict(int))
    for u in sorted(cor.U, key=lambda x: -rank[x]):  # shallowest first
        lines_u = set()
        for a, b, xa, xb, kind in by_u[u]:
            # a is a service or leg, b a main; positions xa on a, xb on b
            ua, ub = a[1], b[1]
            if ua == u and rank[ua] < rank[ub]:  # u's service or leg passes below the main of ub
                need = final_depth(b, xb) + dim(b) + vertical_clearance(cr, ua, ub)
                req[a].append((xa, need))
                lines_u.add(a)
                counts[u]["service_below"] += 1
                counts[u][f"service_below_{kind}"] += 1
                counts[u][f"service_below_{kind}_needs"] += need > depth0(a) + 1e-9
            elif ub == u and rank[ub] < rank[ua]:  # u's main passes below the service or leg of ua
                need = final_depth(a, xa) + dim(a) + vertical_clearance(cr, ua, ub)
                req[b].append((xb, need))
                lines_u.add(b)
                counts[u]["main_below"] += 1
                counts[u][f"main_below_{kind}"] += 1
                counts[u][f"main_below_{kind}_needs"] += need > depth0(b) + 1e-9
        for line in lines_u:
            prof[line] = envelope(length(line), depth0(line), req[line], grad)

    # -- totals ------------------------------------------------------------------------------------
    extra = vr["trench_extra_width_m"]
    out = {}
    for u in cor.U:
        t = dict(main_crossings=counts[u]["main_below"], service_crossings=counts[u]["service_below"],
                 by_kind={k: v for k, v in counts[u].items() if k.count("_") >= 2},
                 mains_lowered=0, mains_lowered_m=0.0, mains_extra_m3=0.0, mains_uniform_m3=0.0,
                 mains_max_m=0.0, services_lowered=0, services_extra_m3=0.0, legs_lowered=0)
        for line, (s, z) in prof.items():
            if line[1] != u:
                continue
            d0 = depth0(line)
            if line[0] == "main":
                w = dim(line) + extra
                v, low, mx = stats(s, z, d0, w)
                if mx > 1e-6:
                    t["mains_lowered"] += 1
                    t["mains_lowered_m"] += low
                    t["mains_extra_m3"] += v
                    t["mains_max_m"] = max(t["mains_max_m"], mx)
                    t["mains_uniform_m3"] += mx * length(line) * w  # whole edge lowered to its deepest need
            else:
                v, low, mx = stats(s, z, d0, vr["service_trench_width_m"])
                if mx > 1e-6:
                    t["services_lowered" if line[0] == "svc" else "legs_lowered"] += 1
                    t["services_extra_m3"] += v
        out[u] = t
    # a consistent stack of standard depths, with no requirement at any crossing
    stack, z = {}, None
    for u in sorted(cor.U, key=lambda x: -rank[x]):
        dmax = max(max(Dm[u].values(), default=0.0), Ds[u]) if u != "electricity" else vr["outer_dimension_m"]["electricity_main"]
        stack[u] = d_std[u] if z is None else max(d_std[u], z)
        z = stack[u] + dmax + max(vertical_clearance(cr, u, w) for w in cor.U if w != u)
    for u in cor.U:
        vol = std = 0.0
        for e in cor.E[u]:
            L = S.edges[e][3]["length"] * len(cor.sides_of(u, e))
            w = dim(("main", u, e)) + extra
            vol += L * w * (stack[u] - d_std[u])
            std += L * w * (d_std[u] + dim(("main", u, e)))
        out[u]["stack_depth_m"] = stack[u]
        out[u]["stack_extra_m3"] = vol
        out[u]["base_m3"] = std
    return dict(totals=out, events=len(events), legs=len(rec["legs"]), prof=prof, req=dict(req),
                inter=dict(rec["inter"]))
