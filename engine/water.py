"""
Water distribution design: Steiner trunk, face-based loops, vertex splitting,
Global Gradient Algorithm (Todini-Pilati) hydraulics with Hazen-Williams,
pipe sizing, fire-flow check, valves and hydrants.
"""
import math

import networkx as nx
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

from . import design

HW_K = 10.67
HW_N = 1.852


def hw_r(L, D, C):
    """Hazen-Williams resistance: h = r |Q|^0.852 Q  (SI: m, m3/s)."""
    return HW_K * L / (C ** HW_N * D ** 4.87)


def capacity_at_gradient(D, C, grad):
    """Flow (m3/s) of a pipe of diameter D at hydraulic gradient grad (m/m)."""
    return (grad * C ** HW_N * D ** 4.87 / HW_K) ** (1 / HW_N)


def inflow_count(total_q, packs):
    """Lower bound on the number of inflows: K = ceil(Q / q_max(D_max, gradient limit))."""
    cap = packs["water"]["capacity"]
    D = max(packs["water"]["standard_sizes"]) / 1000
    qmax = capacity_at_gradient(D, cap["hazen_williams_c"], cap["max_pressure_drop_per_km_m"] / 1000)
    return max(1, math.ceil(total_q / qmax)), qmax


def inflow_candidates(Gc, boundary, river, cands_filter, river_buffer=300.0, ring_tol=40.0):
    """Street vertices on major roads within ring_tol of the land boundary (river bank excluded)."""
    from shapely.geometry import Point
    ring = boundary.exterior
    out = []
    for v, d in Gc.nodes(data=True):
        if not d.get("street_vertex") or not cands_filter(v):
            continue
        p = Point(d["xy"])
        if ring.distance(p) <= ring_tol and (river is None or river.distance(p) > river_buffer):
            out.append(v)
    if not out:
        raise ValueError("no inflow candidates on the boundary")
    return out


def choose_inflows(Gc, cands, K, elev, terminals, p_src, g=0.003):
    """K inflows maximising the worst head margin (greedy k-center on the metric graph):

        margin(t) = max_s [ z_s + p_src - z_t - g d_N(s, t) ]

    with g an assumed trunk gradient (m/m).  Since p(t) <= max_s H_s - z(t), a node with
    margin < p_min for every candidate set needs a separate pressure zone.
    Returns (inflows, worst margin).
    """
    terminals = list(terminals)
    zt = np.array([float(elev(Gc.nodes[t]["xy"])) for t in terminals])
    rows = {}
    for c in cands:
        d = nx.single_source_dijkstra_path_length(Gc, c, weight="length")
        zc = float(elev(Gc.nodes[c]["xy"]))
        rows[c] = zc + p_src - zt - g * np.array([d.get(t, math.inf) for t in terminals])
    best = np.full(len(terminals), -np.inf)
    chosen = []
    while len(chosen) < min(K, len(cands)):
        c = max((c for c in cands if c not in chosen),
                key=lambda c: (np.maximum(best, rows[c]).min(), np.maximum(best, rows[c]).mean()))
        chosen.append(c)
        best = np.maximum(best, rows[c])
    return chosen, float(best.min())


class GGA:
    """Global Gradient Algorithm for a looped network with fixed-head sources."""

    def __init__(self, N, sources, head0, demand, C):
        self.nodes = list(N.nodes())
        self.edges = list(N.edges())
        idx = {v: i for i, v in enumerate(self.nodes)}
        self.fixed = [idx[s] for s in sources]
        free = [i for i in range(len(self.nodes)) if i not in set(self.fixed)]
        self.free = free
        fpos = {i: j for j, i in enumerate(free)}
        hpos = {i: j for j, i in enumerate(self.fixed)}
        m = len(self.edges)
        r12, c12, v12, r10, c10, v10 = [], [], [], [], [], []
        for k, (a, b) in enumerate(self.edges):
            for v, sgn in ((a, -1.0), (b, 1.0)):
                i = idx[v]
                if i in fpos:
                    r12.append(k); c12.append(fpos[i]); v12.append(sgn)
                else:
                    r10.append(k); c10.append(hpos[i]); v10.append(sgn)
        self.A12 = sp.csr_matrix((v12, (r12, c12)), shape=(m, len(free)))
        self.A10 = sp.csr_matrix((v10, (r10, c10)), shape=(m, len(self.fixed)))
        self.H0 = np.array([head0[self.nodes[i]] for i in self.fixed])
        self.q = np.array([demand.get(self.nodes[i], 0.0) for i in free])
        self.L = np.array([N.edges[e]["length"] for e in self.edges])
        self.C = C
        self.idx = idx

    def solve(self, D, Q0=None, tol=1e-6, max_iter=60):
        r = hw_r(self.L, D, self.C)
        Q = np.full(len(self.edges), 1e-3) if Q0 is None else Q0.copy()
        A21 = self.A12.T.tocsr()
        a10h = self.A10 @ self.H0
        for it in range(max_iter):
            aq = np.maximum(np.abs(Q), 1e-7)
            F1 = r * aq ** (HW_N - 1) * Q
            Dd = HW_N * r * aq ** (HW_N - 1)
            Dinv = 1.0 / Dd
            # Newton step on  r|Q|^0.852 Q + A12 H + A10 H0 = 0,  A21 Q = q
            if it == 0:
                H = np.zeros(len(self.free))
            F1 = F1 + self.A12 @ H + a10h
            F2 = A21 @ Q - self.q
            M = (A21 @ sp.diags(Dinv) @ self.A12).tocsc()
            dH = spsolve(M, F2 - A21 @ (Dinv * F1))
            dQ = -Dinv * (F1 + self.A12 @ dH)
            H = H + dH
            Q = Q + dQ
            if np.max(np.abs(dQ)) < tol:
                break
        self.iterations = it + 1
        heads = {self.nodes[i]: H[j] for j, i in enumerate(self.free)}
        heads.update({self.nodes[i]: self.H0[j] for j, i in enumerate(self.fixed)})
        hf = r * np.abs(Q) ** (HW_N - 1) * Q
        return Q, heads, hf


def design_water(Gc, street, taps_by_bldg, bldg, packs, sources, elev, log):
    """Design the water network.  Returns (N, info) with sizes and hydraulic results on N."""
    w = packs["water"]
    cap, topo, gr = w["capacity"], w["topology"], w["graph"]
    tapv = {b: Gc.graph["tap_vertex"][b] for b in taps_by_bldg}
    q_tap = {}
    for b, t in tapv.items():  # a shared service carries the demand of all its buildings
        q_tap[t] = q_tap.get(t, 0.0) + bldg.loc[b, "water_peak_ls"] / 1000.0
    terminals = set(q_tap) | set(sources)

    # 1. Steiner tree on the weighted corridor graph, or (tree variant of P2) the shortest-path
    #    forest from the inflows: every tap on its shortest street path from the nearest inflow
    if gr.get("tree") == "spt":
        H = Gc.copy()
        H.add_edges_from((("super", s) for s in sources), length=0.0, weight=0.0)
        T, _ = design.shortest_path_tree(H, "super", terminals, weight="length")
        T.remove_node("super")
        N = design.copy_attrs(T, Gc)
    else:
        N = design.copy_attrs(design.steiner_tree(Gc, terminals), Gc)
    tree_len = sum(d["length"] for _, _, d in N.edges(data=True))

    # 2. trunk = arcs carrying at least the trunk share of the total demand
    parent, order, down = design.orient(N, sources, {v: {"q": q} for v, q in q_tap.items()})
    total = sum(q_tap.values())
    trunk = [(parent[v], v) for v in order if parent[v] is not None
             and down[v]["q"] >= gr["trunk_demand_share"] * total / len(sources)]

    # 3. face-based loops around trunk bridges (G6)
    aug = design.augment_faces(N, Gc, street, trunk, gr["max_loop_chord_m"]) if gr["bridgeless_trunk"] else {}

    # 4. vertex splitting (G3) and subdivision (G8)
    n_split, n_links = design.split_vertices(N, topo["max_branching"])
    n_sub = design.subdivide(N, Gc, street, topo["max_segment_length_m"])

    # 5. hydraulics and sizing
    z = {v: float(elev(N.nodes[v]["xy"])) for v in N}
    head0 = {s: z[s] + cap["source_pressure_mwc"] for s in sources}
    for v in N:
        N.nodes[v]["z"] = z[v]
    sizes = [d / 1000 for d in w["standard_sizes"]]
    solver = GGA(N, sources, head0, q_tap, cap["hazen_williams_c"])
    # first pass: tree flows at the design velocity
    parent, order, down = design.orient(N, sources, {v: {"q": q} for v, q in q_tap.items()})
    Qt = np.array([down[b]["q"] if parent.get(b) == a else (down[a]["q"] if parent.get(a) == b else 0.0)
                   for a, b in solver.edges])
    D = np.array([_size_for_velocity(abs(q), sizes, cap["design_velocity_ms"]) for q in Qt])
    grad_max = cap["max_pressure_drop_per_km_m"] / 1000
    hmax = max(head0.values())
    Q = None
    for it in range(600):
        Q, H, hf = solver.solve(D, Q)
        v = np.abs(Q) / (math.pi * D ** 2 / 4)
        grad = np.abs(hf) / solver.L
        up = (v > cap["max_velocity_ms"]) | (grad > grad_max * 1.0001)
        # nodes that no pipe size can bring to the minimum pressure (static head bound)
        feasible = [t for t in q_tap if hmax - z[t] >= cap["min_pressure_mwc"]]
        low = [t for t in feasible if H[t] - z[t] < cap["min_pressure_mwc"]]
        if low:
            # upsize the pipes on the least-loss path to the worst feasible node
            worst = min(low, key=lambda t: H[t] - z[t])
            up |= _path_mask(N, solver, sources, worst, hf)
        can = np.array([d < sizes[-1] for d in D])
        if not (up & can).any():
            break
        D = np.where(up & can, [_next(d, sizes) for d in D], D)
    # feasible iff no velocity / gradient violation and no feasible node below the minimum
    feasible_design = not up.any()
    viol = dict(velocity=int((v > cap["max_velocity_ms"]).sum()),
                gradient=int((grad > grad_max * 1.0001).sum()), low_pressure=len(low),
                at_max_size=int(((v > cap["max_velocity_ms"]) & ~can).sum()))

    # 6. fire flow at the hydrant node with the lowest pressure
    v = np.abs(Q) / (math.pi * D ** 2 / 4)
    for k, (a, b) in enumerate(solver.edges):
        e = N.edges[a, b]
        e.update(diameter_mm=round(D[k] * 1000), flow_ls=float(Q[k] * 1000), velocity_ms=float(v[k]),
                 headloss_m=float(hf[k]), gradient_m_per_km=float(abs(hf[k]) / solver.L[k] * 1000))
    for x in N:
        N.nodes[x]["head_m"] = float(H[x])
        N.nodes[x]["pressure_mwc"] = float(H[x] - z[x])

    _hierarchy(N, sources, q_tap, gr["trunk_demand_share"] * total / len(sources))
    comps = _components(N, w)
    fire = _fire_check(N, solver, D, Q, q_tap, comps["hydrants"], z, w, sources)

    info = dict(sources=sources, total_peak_ls=total * 1000, steiner_length_m=tree_len,
                length_m=sum(d["length"] for _, _, d in N.edges(data=True)),
                augmentation=aug, n_split=n_split, n_split_links=n_links, n_subdivided=n_sub,
                gga_iterations=solver.iterations, sizing_passes=it + 1,
                min_pressure_mwc=min(N.nodes[t]["pressure_mwc"] for t in q_tap),
                taps_below_min=sum(N.nodes[t]["pressure_mwc"] < cap["min_pressure_mwc"] - 1e-6 for t in q_tap),
                taps_static_infeasible=sum(hmax - z[t] < cap["min_pressure_mwc"] for t in q_tap),
                source_elevations_m=[z[s] for s in sources],
                max_pressure_mwc=max(N.nodes[t]["pressure_mwc"] for t in N),
                max_velocity_ms=float(v.max()), max_gradient_m_per_km=float((np.abs(hf) / solver.L).max() * 1000),
                feasible=bool(feasible_design), violations=viol, components=comps, fire=fire,
                mu=N.number_of_edges() - N.number_of_nodes() + nx.number_connected_components(N))
    log(f"water: {info['length_m']/1000:.1f} km, mu={info['mu']}, p_min={info['min_pressure_mwc']:.1f} m")
    return N, info


def design_water_auto(Gc, street, taps_by_bldg, bldg, packs, cands, elev, log, k_max=10):
    """Add inflows (K <- K + 1, rule G4 repair) until the sized network meets all hydraulic rules."""
    K, qmax = inflow_count(bldg.loc[list(taps_by_bldg), "water_peak_ls"].sum() / 1000, packs)
    K0 = K
    terminals = sorted({Gc.graph["tap_vertex"][b] for b in taps_by_bldg}, key=str)
    p_src = packs["water"]["capacity"]["source_pressure_mwc"]
    while True:
        sources, margin = choose_inflows(Gc, cands, K, elev, terminals, p_src)
        N, info = design_water(Gc, street, taps_by_bldg, bldg, packs, sources, elev, log)
        info.update(K_lower_bound=K0, K=K, q_max_per_inflow_ls=qmax * 1000, head_margin_m=margin)
        log(f"  K={K}: feasible={info['feasible']} {info['violations']}")
        if info["feasible"] or K >= min(k_max, len(cands)):
            return N, info
        K += 1


def _size_for_velocity(q, sizes, vmax):
    for d in sizes:
        if q <= vmax * math.pi * d * d / 4:
            return d
    return sizes[-1]


def _next(d, sizes):
    for s in sizes:
        if s > d + 1e-9:
            return s
    return d


def _path_mask(N, solver, sources, target, hf):
    """Edges on the path from the nearest source to target that minimises the head loss."""
    idx = {e: k for k, e in enumerate(solver.edges)}
    H = nx.Graph()
    for k, (a, b) in enumerate(solver.edges):
        H.add_edge(a, b, w=abs(hf[k]) + 1e-9)
    path = None
    best = math.inf
    for s in sources:
        try:
            L, p = nx.single_source_dijkstra(H, s, target, weight="w")
        except nx.NetworkXNoPath:
            continue
        if L < best:
            best, path = L, p
    mask = np.zeros(len(solver.edges), dtype=bool)
    if path:
        for a, b in zip(path[:-1], path[1:]):
            k = idx.get((a, b), idx.get((b, a)))
            mask[k] = True
    return mask


def _hierarchy(N, sources, q_tap, trunk_q):
    for a, b, d in N.edges(data=True):
        d["hierarchy"] = "trunk" if abs(d.get("flow_ls", 0)) / 1000 >= trunk_q else "secondary"


def _components(N, w):
    """Valves at every branch of a tee and every 200 m; hydrants every 100 m on mains >= 100 mm."""
    comp = w["components"]
    valves = design.place_along_runs(N, comp["valves"]["required_spacing_m"])
    # every branch of a tee gets a valve (required_at_branches); taps are tees with the service
    branch_valves = [(v, u) for v in N if N.degree(v) == 3 for u in N.neighbors(v)]
    dmin = w["capacity"]["fire_flow_min_diameter_mm"]
    hyd = design.place_along_runs(N, comp["hydrants"]["spacing_m"],
                                  select=lambda d: d.get("diameter_mm", 0) >= dmin)
    return {"valves": valves, "branch_valves": branch_valves, "hydrants": hyd}


def _fire_check(N, solver, D, Q, q_tap, hydrants, z, w, sources):
    """Peak hour plus the fire flow at the hydrant of lowest pressure."""
    if not hydrants:
        return {}
    cap = w["capacity"]
    cand = []
    for ri, s, (a, b), off in hydrants:
        node = a if off < N.edges[a, b]["length"] / 2 else b
        cand.append(node)
    worst = min(cand, key=lambda x: N.nodes[x]["pressure_mwc"])
    q = dict(q_tap)
    q[worst] = q.get(worst, 0.0) + cap["fire_flow_ls"] / 1000
    s2 = GGA(N, sources, {s: z[s] + cap["source_pressure_mwc"] for s in sources}, q, cap["hazen_williams_c"])
    Q2, H2, _ = s2.solve(D, Q)
    return {"node": str(worst), "pressure_normal_mwc": N.nodes[worst]["pressure_mwc"],
            "residual_mwc": float(H2[worst] - z[worst]),
            "min_pressure_other_mwc": float(min(H2[t] - z[t] for t in q_tap if t != worst))}
