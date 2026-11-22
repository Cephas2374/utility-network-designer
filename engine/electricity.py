"""
Electricity distribution design: transformer siting by graph-Voronoi refinement,
LV shortest-path forest, vertex splitting, voltage-drop cable sizing, normally-open
LV ties between adjacent cells, and 20 kV open rings by sweep + Christofides.
"""
import math
from collections import defaultdict

import networkx as nx
import numpy as np

from . import design
from .demand import diversity_factor, diversity_factor_monotone


def _g(n, ed):
    """Diversity factor of n units: the rule-pack table, or its monotone majorant when the demand
    pack sets monotone_diversity (paper P3; P1 used the table)."""
    f = diversity_factor_monotone if ed.get("monotone_diversity") else diversity_factor
    return f(n, ed["diversity_factors"])


def _kva(recs, ed, pf):
    P = sum(r["P"] for r in recs)
    n = sum(r["units"] for r in recs)
    return _g(n, ed) * P / pf if recs else 0.0


def design_electricity(Gc, taps_by_bldg, bldg, packs, mv_infeed, service_len, log):
    e = packs["electricity"]
    cap, topo, gr, ed = e["capacity"], e["topology"], e["graph"], e["demand"]
    comp = e["components"]["transformers"]
    pf = cap["power_factor"]
    rating_max = max(comp["standard_ratings_kva"])
    load_cap = rating_max * e["validation"]["max_transformer_loading_percent"] / 100

    tapv = {b: Gc.graph["tap_vertex"][b] for b in taps_by_bldg}
    # loads per tap (a shared service feeds all its buildings); taps above the MV threshold get
    # their own substation on the MV ring
    agg = {}
    for b, t in tapv.items():
        d = agg.setdefault(t, {"P": 0.0, "units": 0, "w": 0.0, "svc": 0.0, "bldgs": []})
        d["P"] += float(bldg.loc[b, "elec_kw"])
        d["units"] += int(bldg.loc[b, "units"])
        d["w"] += float(bldg.loc[b, "elec_kw"])
        d["svc"] += float(bldg.loc[b, "elec_service_kw"]) / pf
        d["bldgs"].append(b)
    svc_kva = {t: d["svc"] for t, d in agg.items()}
    mv_taps = [t for t, d in agg.items() if d["svc"] >= cap["mv_customer_min_kva"]]
    mv_customers = [b for t in mv_taps for b in agg[t]["bldgs"]]
    demand = {t: {k: d[k] for k in ("P", "units", "w")} for t, d in agg.items() if t not in set(mv_taps)}
    load_fn = lambda recs: _kva(recs, ed, pf)
    candidate = lambda v: not isinstance(v, tuple) or v[0] == "t"  # junctions or in front of a plot

    # 1-3. transformer sites (graph-Voronoi cells within capacity and reach, G10), LV
    #      shortest-path forest, cable sizing; a cell whose voltage drop cannot be met with the
    #      largest cable is split towards its worst customer and the cells are recomputed.
    init, rounds = None, 0
    while True:
        rounds += 1
        sites, dist, pred, owner = design.voronoi_refine(Gc, demand, load_fn, load_cap,
                                                         topo["transformer_max_reach_m"], candidate, init=init)
        sites = [s for s in sites if any(owner[t] == s for t in demand)]
        F = design.forest_from_pred(Gc, pred, list(demand))
        for x in F:
            F.nodes[x]["cell"] = owner[x]
        for u, v in F.edges():
            F.edges[u, v]["cell"] = owner[u]
        exempt = {s: gr["source_max_degree"] for s in sites}
        n_split, n_links = design.split_vertices(F, topo["max_branching"], exempt=exempt)
        for x in F:
            F.nodes[x].setdefault("cell", owner.get(F.nodes[x].get("split_of", x)))
        for u, v, d in F.edges(data=True):
            d.setdefault("cell", F.nodes[u]["cell"])
        n_sub = 0  # LV runs between taps are far below the span limit; checked in validation

        cells = defaultdict(list)
        for t in demand:
            cells[owner[t]].append(t)
        trafos, worst_drop, infeasible, split_at = {}, 0.0, 0, []
        for s in sites:
            res = _size_cell(F, s, cells[s], demand, e, service_len, bldg, tapv)
            kva = load_fn([demand[t] for t in cells[s]])
            rating = next((r for r in comp["standard_ratings_kva"]
                           if kva <= r * e["validation"]["max_transformer_loading_percent"] / 100), rating_max)
            trafos[s] = dict(kva=kva, rating_kva=rating, loading=kva / rating, n_customers=len(cells[s]),
                             max_drop_pct=res["max_drop_pct"], max_reach_m=max(dist[t] for t in cells[s]),
                             feeders=F.degree(s) if s in F else 0)
            worst_drop = max(worst_drop, res["max_drop_pct"])
            infeasible += res["infeasible"]
            if res["infeasible"] and len(cells[s]) > 1:
                split_at.append(design._walk_back(Gc, pred, res["worst"], dist[res["worst"]] / 2,
                                                  candidate, sites))
        if not split_at or rounds >= 25:
            break
        init = sites + split_at
        log(f"  electricity round {rounds}: {len(split_at)} cells split for voltage drop")

    # 4. normally-open LV ties between adjacent cells (shortest Gc edge between two LV trees)
    ties = {}
    for a, b, d in Gc.edges(data=True):
        oa, ob = owner.get(a), owner.get(b)
        if oa is None or ob is None or oa == ob or a not in F or b not in F:
            continue
        key = tuple(sorted((oa, ob), key=str))
        if key not in ties or d["length"] < ties[key][2]["length"]:
            ties[key] = (a, b, d)

    # 5. MV open rings: sweep by angle around the infeed, Christofides tour per group;
    #    MV customer substations are ring nodes like the distribution transformers
    for v in mv_taps:
        k = svc_kva[v]
        rating = next((r for r in comp["standard_ratings_kva"] if k <= r * 0.8), math.ceil(k / 0.8 / 500) * 500)
        trafos[v] = dict(kva=k, rating_kva=rating, loading=k / rating, n_customers=len(agg[v]["bldgs"]),
                         max_drop_pct=0.0, max_reach_m=0.0, feeders=0, mv_customer=agg[v]["bldgs"][0],
                         mv_customer_bldgs=agg[v]["bldgs"])
    ring_nodes = sites + mv_taps
    rings = _mv_rings(Gc, mv_infeed, ring_nodes, trafos, gr["mv_ring_max_kva"])

    lv_len = sum(d["length"] for _, _, d in F.edges(data=True))
    lvt = [t for t in trafos.values() if "mv_customer" not in t]
    info = dict(n_transformers=len(sites), n_mv_customers=len(mv_taps),
                mv_customer_kva=sum(svc_kva[t] for t in mv_taps),
                installed_kva=sum(t["rating_kva"] for t in lvt),
                diversified_kva=sum(t["kva"] for t in lvt),
                max_loading=max(t["loading"] for t in lvt),
                lv_length_m=lv_len, n_split=n_split, n_split_links=n_links, n_subdivided=n_sub,
                max_drop_pct=worst_drop, taps_drop_infeasible=infeasible,
                max_reach_m=max(t["max_reach_m"] for t in lvt),
                max_feeders=max(t["feeders"] for t in lvt),
                n_lv_ties=len(ties), mv_rings=len(rings), siting_rounds=rounds,
                mv_length_m=sum(r["length_m"] for r in rings),
                mu_closed=F.number_of_edges() - F.number_of_nodes() + nx.number_connected_components(F),
                mu_with_ties=F.number_of_edges() + len(ties) - F.number_of_nodes() + 1)
    log(f"electricity: {len(sites)} transformers, LV {lv_len/1000:.1f} km, MV {info['mv_length_m']/1000:.1f} km, "
        f"max drop {worst_drop:.2f} %")
    return F, dict(info=info, transformers=trafos, ties=ties, rings=rings, sites=sites, owner=owner)


def _size_cell(F, s, taps, demand, e, service_len, bldg, tapv_inv_unused):
    """Size the LV cables of one transformer tree: ampacity first, then upsize along the
    worst path until the voltage drop at every customer is within the limit."""
    cap, ed = e["capacity"], e["demand"]
    V, pf = cap["nominal_voltage_v"], cap["power_factor"]
    sinphi = math.sqrt(1 - pf ** 2)
    rho, X = cap["conductor_resistivity_ohmm"], cap["conductor_reactance_ohm_per_km"] / 1000
    amp = {float(k): v * cap["ampacity_derating_factor"] for k, v in e["cable_capacity"].items()}
    sizes = [a for a in e["standard_sizes"] if a >= cap["min_main_size_mm2"]]
    limit = cap["max_voltage_drop_percent"]

    T = F.subgraph(nx.node_connected_component(F, s)).copy() if s in F else nx.Graph()
    loads = {t: demand[t] for t in taps}
    parent, order, down = design.orient(T, [s], {t: {"P": loads[t]["P"], "n": loads[t]["units"]} for t in taps})
    arcs = [(parent[v], v) for v in order if parent[v] is not None]
    I = {}
    for p, v in arcs:
        P = down[v]["P"] * _g(int(down[v]["n"]), ed) * 1000
        I[(p, v)] = P / (math.sqrt(3) * V * pf)
    size = {}
    for arc in arcs:
        size[arc] = next((a for a in sizes if amp[a] >= I[arc]), sizes[-1])

    def drops():
        cum = {s: 0.0}
        for p, v in arcs:
            L = F.edges[p, v]["length"]
            dv = math.sqrt(3) * I[(p, v)] * L * (rho / (size[(p, v)] * 1e-6) * pf + X * sinphi)
            cum[v] = cum[p] + dv / V * 100
        return cum

    infeasible = 0
    svc_size = {}
    sd = lambda t: _service_drop(t, loads, service_len, e, svc_size.get(t))
    for _ in range(600):
        cum = drops()
        worst = max(taps, key=lambda t: cum[t] + sd(t))
        if cum[worst] + sd(worst) <= limit:
            break
        path, x = [], worst
        while parent[x] is not None:
            path.append((parent[x], x))
            x = parent[x]
        up = [a for a in path if size[a] < sizes[-1]]
        # candidates: the arcs of the worst path and the service cable of the worst customer
        cand = [(I[a] * F.edges[a]["length"] / size[a], "arc", a) for a in up]
        A_svc = svc_size.get(worst) or _service_size(worst, loads, e)
        if A_svc < e["standard_sizes"][-1]:
            cand.append((_svc_current(worst, loads, e) * service_len.get(worst, 10.0) / A_svc, "svc", worst))
        if not cand:
            infeasible += 1
            break
        _, kind, a = max(cand, key=lambda c: c[0])
        if kind == "arc":
            size[a] = sizes[sizes.index(size[a]) + 1]
        else:
            std = e["standard_sizes"]
            svc_size[a] = std[std.index(A_svc) + 1]
    cum = drops()
    for t, A in svc_size.items():
        F.nodes[t]["service_mm2"] = A
    for (p, v), A in size.items():
        F.edges[p, v].update(size_mm2=A, current_a=I[(p, v)], hierarchy=None)
    for t in taps:
        F.nodes[t]["drop_pct"] = cum[t] + sd(t)
    total = down[s]["P"] if s in down else 0.0
    for p, v in arcs:
        F.edges[p, v]["hierarchy"] = "primary" if down[v]["P"] >= 0.3 * total else "secondary"
    worst = max(taps, key=lambda t: F.nodes[t]["drop_pct"])
    return {"max_drop_pct": F.nodes[worst]["drop_pct"], "infeasible": infeasible, "worst": worst}


def _svc_current(t, loads, e):
    cap = e["capacity"]
    P = loads[t]["P"] * _g(loads[t]["units"], e["demand"]) * 1000
    return P / (math.sqrt(3) * cap["nominal_voltage_v"] * cap["power_factor"])


def _service_size(t, loads, e):
    """Smallest service cable meeting the ampacity of the customer load."""
    cap = e["capacity"]
    I = _svc_current(t, loads, e)
    amp = {float(k): v * cap["ampacity_derating_factor"] for k, v in e["cable_capacity"].items()}
    return next((a for a in e["standard_sizes"] if a >= cap["min_service_size_mm2"] and amp[a] >= I),
                e["standard_sizes"][-1])


def _service_drop(t, loads, service_len, e, A=None):
    """Voltage drop (%) in the service cable of size A (default: ampacity size), length from the tap."""
    cap = e["capacity"]
    V, pf = cap["nominal_voltage_v"], cap["power_factor"]
    I = _svc_current(t, loads, e)
    A = A or _service_size(t, loads, e)
    L = service_len.get(t, 10.0)
    rho, X = cap["conductor_resistivity_ohmm"], cap["conductor_reactance_ohm_per_km"] / 1000
    return math.sqrt(3) * I * L * (rho / (A * 1e-6) * pf + X * math.sqrt(1 - pf ** 2)) / V * 100


def _mv_rings(Gc, infeed, sites, trafos, ring_kva):
    """Sweep heuristic: transformers sorted by angle around the infeed, grouped up to the ring
    capacity; each group is closed as a cycle through the infeed (Christofides on the metric
    closure) and opened at the point that balances the two half-feeders."""
    x0, y0 = Gc.nodes[infeed]["xy"]
    ang = {s: math.atan2(Gc.nodes[s]["xy"][1] - y0, Gc.nodes[s]["xy"][0] - x0) for s in sites}
    order = sorted(sites, key=ang.get)
    # start the sweep at the largest angular gap so that groups are compact sectors
    gaps = [((ang[order[(i + 1) % len(order)]] - ang[order[i]]) % (2 * math.pi), i) for i in range(len(order))]
    start = (max(gaps)[1] + 1) % len(order)
    order = order[start:] + order[:start]
    groups, cur, kva = [], [], 0.0
    for s in order:
        r = trafos[s]["rating_kva"]
        if cur and kva + r > ring_kva:
            groups.append(cur)
            cur, kva = [], 0.0
        cur.append(s)
        kva += r
    if cur:
        groups.append(cur)

    rings = []
    for g in groups:
        nodes = [infeed] + g
        K = nx.Graph()
        paths = {}
        for i, a in enumerate(nodes):
            d, p = nx.single_source_dijkstra(Gc, a, weight="length")
            for b in nodes[i + 1:]:
                K.add_edge(a, b, weight=d[b])
                paths[(a, b)] = p[b]
                paths[(b, a)] = p[b][::-1]
        if len(nodes) == 2:
            tour = [infeed, g[0], infeed]
        else:
            tour = nx.approximation.traveling_salesman_problem(K, cycle=True, method=nx.approximation.christofides)
        # rotate so that the tour starts and ends at the infeed
        i = tour.index(infeed)
        tour = tour[i:-1] + tour[:i] + [infeed] if tour[0] != infeed else tour
        # the TSP wrapper expands ties of the metric closure through intermediate nodes, so a
        # substation on the shortest path between two others can appear twice; keep the first visit
        # (the leg between its neighbours runs through it anyway, so no length changes; paper P3)
        tour = list(dict.fromkeys(tour[:-1])) + [infeed]
        legs = [paths[(a, b)] for a, b in zip(tour[:-1], tour[1:])]
        length = sum(K.edges[a, b]["weight"] for a, b in zip(tour[:-1], tour[1:]))
        # normally-open point: between the transformers where the cumulative kVA passes half
        tot = sum(trafos[s]["rating_kva"] for s in g)
        acc, nop = 0.0, 0
        for j, s in enumerate(tour[1:-1]):
            acc += trafos[s]["rating_kva"]
            if acc >= tot / 2:
                nop = j + 1
                break
        rings.append(dict(transformers=g, tour=tour, legs=legs, length_m=length, kva=tot, nop_leg=nop))
    return rings
