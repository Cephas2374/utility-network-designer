"""
Gas distribution design (prospective scenario for Kinshasa): district regulators by
graph-Voronoi refinement with hazard clearance, LP Steiner trees inside the cells,
an MP Steiner tree from the city gate, Renouard sizing on trees (flows are unique),
valves and emergency shut-offs.
"""
import math
from collections import defaultdict

import networkx as nx
import numpy as np

from . import design

P_ATM = 1.01325  # bar


def renouard_lp_drop(Q, D_mm, L_m, s, c=23200):
    """P1 - P2 [mbar] = c s L[m] Q^1.82 D^-4.82   (Q m3/h at standard conditions, D internal mm)."""
    return c * s * L_m * Q ** 1.82 * D_mm ** -4.82


def renouard_mp_p2(P1_abs, Q, D_mm, L_m, s, c=48600):
    """P2 from P1^2 - P2^2 [bar abs^2] = c s L[km] Q^1.82 D^-4.82."""
    v = P1_abs ** 2 - c * s * (L_m / 1000) * Q ** 1.82 * D_mm ** -4.82
    return math.sqrt(v) if v > 0 else 0.0


def velocity(Q, D_mm, P_abs):
    """Mean velocity (m/s) of a standard flow Q (m3/h) at absolute pressure P_abs (bar)."""
    A = math.pi * (D_mm / 1000) ** 2 / 4
    return Q / 3600 * (P_ATM / P_abs) / A


def design_gas(Gc, taps_by_bldg, bldg, packs, city_gate, hazard_geoms, log):
    g = packs["gas"]
    cap, topo, gr = g["capacity"], g["topology"], g["graph"]
    reg = g["components"]["regulators"]
    s_rel = cap["relative_density"]
    sdr = cap["pipe_sdr"]
    ID = {od: od * (1 - 2 / sdr) for od in g["standard_sizes"]}
    tapv = {b: Gc.graph["tap_vertex"][b] for b in taps_by_bldg}
    demand = {}
    for b, t in tapv.items():
        q = float(bldg.loc[b, "gas_peak_m3h"])
        d = demand.setdefault(t, {"Q": 0.0, "w": 0.0})
        d["Q"] += q
        d["w"] += q
    load_fn = lambda recs: sum(r["Q"] for r in recs)

    # district regulator stations: not within the hazard clearance of schools / hospitals
    from shapely import STRtree
    from shapely.geometry import Point
    tree = STRtree(hazard_geoms) if len(hazard_geoms) else None
    clear = g["safety"]["min_clearance_from_hazards_m"]

    def candidate(v):
        if isinstance(v, tuple):
            return False
        if tree is None:
            return True
        p = Point(Gc.nodes[v]["xy"])
        return len(tree.query(p, predicate="dwithin", distance=clear)) == 0

    sites, dist, pred, owner = design.voronoi_refine(Gc, demand, load_fn, reg["max_capacity_m3_h"],
                                                     gr["lp_zone_max_reach_m"], candidate)
    sites = [s for s in sites if any(owner[t] == s for t in demand)]
    sites, n_removed, n_conflicts = design.enforce_spacing(Gc, sites, demand, load_fn, reg["max_capacity_m3_h"],
                                                           gr["lp_zone_max_reach_m"], topo["min_regulator_spacing_m"])
    dist, pred, owner = design.multi_source_dijkstra(Gc, sites)

    # LP: one Steiner tree per cell, on the subgraph induced by the cell (connected)
    cells = defaultdict(list)
    for t in demand:
        cells[owner[t]].append(t)
    members = defaultdict(set)
    for v, o in owner.items():
        members[o].add(v)
    LP = nx.Graph()
    for s in sites:
        H = Gc.subgraph(members[s])
        T = design.steiner_tree(H, set(cells[s]) | {s}) if len(cells[s]) > 0 else nx.Graph()
        T.add_node(s)
        T = design.copy_attrs(T, Gc)
        for a, b in T.edges():
            T.edges[a, b]["zone"] = s
            T.edges[a, b]["tier"] = "LP"
        LP = nx.union(LP, T) if not set(LP) & set(T) else nx.compose(LP, T)
    for x in LP:
        LP.nodes[x]["zone"] = owner.get(x)

    # MP: Steiner tree from the city gate to the regulators
    MP = design.copy_attrs(design.steiner_tree(Gc, set(sites) | {city_gate}), Gc)
    for a, b in MP.edges():
        MP.edges[a, b]["tier"] = "MP"

    # vertex splitting in each tier (G3); a regulator joins the two tiers
    n_split = design.split_vertices(LP, topo["max_branching"])[0] + design.split_vertices(MP, topo["max_branching"])[0]

    # LP sizing: each arc carries its downstream demand (unique in a tree)
    lp_res = _size_lp(LP, sites, cells, demand, g, ID, s_rel)
    mp_res = _size_mp(MP, city_gate, sites, cells, demand, g, ID, s_rel)

    info = dict(n_regulators=len(sites), total_peak_m3h=sum(r["Q"] for r in demand.values()),
                lp_length_m=sum(d["length"] for _, _, d in LP.edges(data=True)),
                mp_length_m=sum(d["length"] for _, _, d in MP.edges(data=True)),
                n_split=n_split, lp=lp_res, mp=mp_res,
                mu_lp=LP.number_of_edges() - LP.number_of_nodes() + nx.number_connected_components(LP),
                mu_mp=MP.number_of_edges() - MP.number_of_nodes() + nx.number_connected_components(MP),
                max_zone_reach_m=max(max(dist[t] for t in cells[s]) for s in sites),
                max_regulator_load_m3h=max(load_fn([demand[t] for t in cells[s]]) for s in sites),
                min_regulator_spacing_m=_min_spacing(Gc, sites), spacing_removed=n_removed,
                spacing_conflicts=n_conflicts)
    log(f"gas: {len(sites)} regulators, LP {info['lp_length_m']/1000:.1f} km, MP {info['mp_length_m']/1000:.1f} km")
    return LP, MP, dict(info=info, sites=sites, owner=owner)


def _min_spacing(Gc, sites):
    if len(sites) < 2:
        return math.inf
    xy = np.array([Gc.nodes[s]["xy"] for s in sites])
    d = np.sqrt(((xy[:, None, :] - xy[None, :, :]) ** 2).sum(-1))
    d[np.diag_indices(len(sites))] = np.inf
    return float(d.min())


def _size_lp(LP, sites, cells, demand, g, ID, s_rel):
    cap = g["capacity"]
    p_out = g["components"]["regulators"]["outlet_pressure_bar"] * 1000  # mbar
    p_min = cap["min_pressure_bar"] * 1000
    grad = cap["max_pressure_drop_per_100m_mbar"] / 100  # mbar per m
    vmax = cap["max_velocity_lp_ms"]
    sizes = [od for od in g["standard_sizes"] if od >= 32]
    parent, order, down = design.orient(LP, sites, {t: {"Q": demand[t]["Q"]} for t in demand})
    p = {}
    worst = math.inf
    qd_max = 0.0
    for v in order:
        if parent[v] is None:
            p[v] = p_out
            continue
        u = parent[v]
        e = LP.edges[u, v]
        Q = down[v]["Q"]
        od = next((d for d in sizes if renouard_lp_drop(Q, ID[d], 1.0, s_rel) <= grad
                   and velocity(Q, ID[d], P_ATM + p_out / 1000) <= vmax), sizes[-1])
        dp = renouard_lp_drop(Q, ID[od], e["length"], s_rel)
        e.update(diameter_mm=od, flow_m3h=Q, drop_mbar=dp, gradient_mbar_per_100m=dp / e["length"] * 100)
        p[v] = p[u] - dp
        qd_max = max(qd_max, Q / ID[od])
    for v in LP:
        LP.nodes[v]["pressure_mbar"] = p.get(v, p_out)
    tp = [p[t] for t in demand if t in p]
    return dict(min_pressure_mbar=min(tp), taps_below_min=sum(x < p_min for x in tp),
                max_gradient_mbar_per_100m=max(d.get("gradient_mbar_per_100m", 0) for _, _, d in LP.edges(data=True)),
                max_q_over_d=qd_max)


def _size_mp(MP, gate, sites, cells, demand, g, ID, s_rel):
    cap = g["capacity"]
    P0 = cap["inlet_pressure_bar"] + P_ATM
    pmin = cap["min_mp_pressure_bar"] + P_ATM
    vmax = cap["max_velocity_mp_ms"]
    sizes = [od for od in g["standard_sizes"] if od >= 63]
    load = {s: {"Q": sum(demand[t]["Q"] for t in cells[s])} for s in sites}
    parent, order, down = design.orient(MP, [gate], load)
    size = {}
    for v in order:
        if parent[v] is None:
            continue
        Q = down[v]["Q"]
        size[(parent[v], v)] = next((d for d in sizes if velocity(Q, ID[d], pmin) <= vmax), sizes[-1])

    def pressures():
        P = {gate: P0}
        for v in order:
            u = parent[v]
            if u is None:
                continue
            P[v] = renouard_mp_p2(P[u], down[v]["Q"], ID[size[(u, v)]], MP.edges[u, v]["length"], s_rel)
        return P

    for _ in range(200):
        P = pressures()
        low = [s for s in sites if P[s] < pmin]
        if not low:
            break
        x, path = min(low, key=P.get), []
        while parent[x] is not None:
            path.append((parent[x], x))
            x = parent[x]
        up = [a for a in path if size[a] < sizes[-1]]
        if not up:
            break
        a = max(up, key=lambda a: down[a[1]]["Q"] ** 1.82 * MP.edges[a]["length"] / ID[size[a]] ** 4.82)
        size[a] = sizes[sizes.index(size[a]) + 1]
    P = pressures()
    for (u, v), od in size.items():
        MP.edges[u, v].update(diameter_mm=od, flow_m3h=down[v]["Q"],
                              velocity_ms=velocity(down[v]["Q"], ID[od], P[v]))
    for v in MP:
        MP.nodes[v]["pressure_bar_g"] = P.get(v, P0) - P_ATM
    return dict(min_regulator_inlet_bar_g=min(P[s] for s in sites) - P_ATM,
                max_velocity_ms=max(d.get("velocity_ms", 0) for _, _, d in MP.edges(data=True)))
