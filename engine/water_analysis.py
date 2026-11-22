"""
Analysis of the designed water network (paper P2): pressure-driven hydraulics, valve
segments, N-1 segment failures, District Metered Areas, fire flow at every hydrant,
sensitivity, bill of quantities, EPANET export and cross-check.
"""
import math
import subprocess
from collections import defaultdict
from pathlib import Path

import networkx as nx
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

from . import design
from .water import GGA, HW_N, hw_r

ROOT = Path(__file__).resolve().parents[1]
EPANET_PY = ROOT / "dev" / ".venv_epanet" / "Scripts" / "python.exe"


# -- hydraulic model on a (possibly reduced) network ---------------------------------
class Hydraulics:
    """Demand-driven (GGA) and pressure-driven solutions of a sized network."""

    def __init__(self, N, sources, packs):
        self.N = N
        self.sources = list(sources)
        cap = packs["water"]["capacity"]
        self.C = cap["hazen_williams_c"]
        self.p_src = cap["source_pressure_mwc"]
        self.p_req = cap["min_pressure_mwc"]
        self.z = {v: N.nodes[v]["z"] for v in N}
        self.q_req = {v: N.nodes[v].get("q_req", 0.0) for v in N}
        self._base = self._setup(N)  # the full network's solver, reused by every call on it

    def _setup(self, H):
        """Solver of H restricted to the components that contain a source."""
        keep = set()
        for c in nx.connected_components(H):
            if any(s in c for s in self.sources):
                keep |= c
        Hk = H.subgraph(keep)
        solver = GGA(Hk, [s for s in self.sources if s in keep], {s: 0.0 for s in self.sources}, {}, self.C)
        D = np.array([Hk.edges[e]["diameter_mm"] / 1000 for e in solver.edges])
        return dict(solver=solver, D=D, keep=keep, Q=None)

    def closed(self, nodes=(), edges=()):
        """Solver context of the full network with some vertices and pipes out of service, cut from
        the full network's incidence matrices (no graph copy).  Only the part still joined to a
        source is kept; the flows of the full network are the warm start."""
        base = self._base
        b = base["solver"]
        if "_ends" not in base:
            ea = np.array([b.idx[a] for a, _ in b.edges])
            eb = np.array([b.idx[v] for _, v in b.edges])
            base["_ends"] = (ea, eb, {e: k for k, e in enumerate(b.edges)})
        ea, eb, eidx = base["_ends"]
        n = len(b.nodes)
        alive = np.ones(n, dtype=bool)
        alive[[b.idx[v] for v in nodes if v in b.idx]] = False
        open_ = alive[ea] & alive[eb]
        for a, v in edges:
            k = eidx.get((a, v), eidx.get((v, a)))
            if k is not None:
                open_[k] = False
        adj = sp.coo_matrix((np.ones(open_.sum()), (ea[open_], eb[open_])), shape=(n, n))
        _, lab = connected_components(adj, directed=False)
        src_lab = {lab[i] for i in b.fixed if alive[i]}
        keep_mask = alive & np.isin(lab, list(src_lab))
        rows = np.flatnonzero(open_ & keep_mask[ea])
        fcols = [j for j, i in enumerate(b.free) if keep_mask[i]]
        hcols = [j for j, i in enumerate(b.fixed) if keep_mask[i]]
        s = GGA.__new__(GGA)
        s.nodes, s.idx, s.C = b.nodes, b.idx, b.C
        s.edges = [b.edges[k] for k in rows]
        s.free = [b.free[j] for j in fcols]
        s.fixed = [b.fixed[j] for j in hcols]
        s.A12 = b.A12[rows][:, fcols]
        s.A10 = b.A10[rows][:, hcols]
        s.L = b.L[rows]
        keep = {b.nodes[i] for i in np.flatnonzero(keep_mask)}
        Q0 = None if base["Q"] is None else base["Q"][rows]
        return dict(solver=s, D=base["D"][rows], keep=keep, Q=Q0)

    def _solve(self, ctx, demand, head_add=0.0, C=None):
        """DDA in a solver context; the last flows of the context are the warm start."""
        solver = ctx["solver"]
        solver.H0 = np.array([self.z[solver.nodes[i]] + self.p_src + head_add for i in solver.fixed])
        solver.q = np.array([demand.get(solver.nodes[i], 0.0) for i in solver.free])
        solver.C = C or self.C
        Q, heads, hf = solver.solve(ctx["D"], ctx["Q"], max_iter=80)
        ctx["Q"] = Q
        return heads, dict(zip(solver.edges, Q)), ctx["keep"]

    def _ctx(self, H):
        """H: None (full network), a solver context from closed(), or a subgraph of the network."""
        if isinstance(H, dict):
            return H
        return self._base if H is None or H is self.N else self._setup(H)

    def dda(self, H=None, scale=1.0, extra=None, head_add=0.0, C=None):
        dem = {v: q * scale for v, q in self.q_req.items()}
        for v, q in (extra or {}).items():
            dem[v] = dem.get(v, 0.0) + q
        heads, flows, keep = self._solve(self._ctx(H), dem, head_add, C)
        p = {v: heads[v] - self.z[v] for v in keep}
        return p, flows, keep, dem

    def pda(self, H=None, scale=1.0, extra=None, head_add=0.0, C=None, iters=40, relax=0.5, tol=1e-6):
        """Pressure-driven demand (Wagner): d = q (p / p_req)^0.5 for 0 < p < p_req, by a damped
        fixed point around the demand-driven solver.  Returns served demand per node."""
        ctx = self._ctx(H)
        req = {v: q * scale for v, q in self.q_req.items()}
        for v, q in (extra or {}).items():
            req[v] = req.get(v, 0.0) + q
        dem = dict(req)
        keep = set()
        for _ in range(iters):
            heads, flows, keep = self._solve(ctx, dem, head_add, C)
            p = {v: heads[v] - self.z[v] for v in keep}
            new = {}
            change = 0.0
            for v, q in req.items():
                if v not in keep or q == 0:
                    new[v] = 0.0
                    continue
                s = min(1.0, max(0.0, p[v] / self.p_req)) ** 0.5
                d = (1 - relax) * dem.get(v, q) + relax * q * s
                change = max(change, abs(d - dem.get(v, 0.0)))
                new[v] = d
            dem = new
            if change < tol:
                break
        served = sum(dem[v] for v in keep if v in dem)
        return dict(pressure=p, served=dem, keep=keep, served_total=served, required_total=sum(req.values()))


def worst_supply_path(N, sources, taps):
    """The tap of lowest pressure in a sized design and its best supply path: the path from the
    inflow that loses least head on the way (the only path in a tree), with its length, head drop,
    diameters and the straight-line distance it replaces."""
    w = min(taps, key=lambda v: N.nodes[v]["pressure_mwc"])
    best = None
    for s in sources:
        try:
            path = nx.shortest_path(N, s, w, weight="length")
        except nx.NetworkXNoPath:
            continue
        drop = N.nodes[s]["head_m"] - N.nodes[w]["head_m"]
        if best is None or N.nodes[s]["head_m"] - drop > best[0]:
            best = (N.nodes[s]["head_m"] - drop, s, path, drop)
    _, s, path, drop = best
    es = list(zip(path[:-1], path[1:]))
    L = sum(N.edges[e]["length"] for e in es)
    D = np.array([N.edges[e]["diameter_mm"] for e in es])
    Lmax = sum(N.edges[e]["length"] for e in es if N.edges[e]["diameter_mm"] == D.max())
    return dict(pressure=N.nodes[w]["pressure_mwc"], z=N.nodes[w]["z"], source_head=N.nodes[s]["head_m"],
                static_bound=max(N.nodes[x]["head_m"] for x in sources) - N.nodes[w]["z"],
                path_length_m=L, head_drop_m=drop, gradient_m_per_km=1000 * drop / L,
                d_max_mm=int(D.max()), share_at_d_max=Lmax / L,
                straight_m=float(math.dist(N.nodes[s]["xy"], N.nodes[w]["xy"])))


# -- valves and segments ----------------------------------------------------------------
def valve_edges(N, comps):
    """Edges that carry a valve: every branch of a tee (at the tee) and the interval valves."""
    ve = set()
    for v, u in comps["branch_valves"]:
        ve.add(frozenset((v, u)))
    for ri, s, (a, b), off in comps["valves"]:
        ve.add(frozenset((a, b)))
    return ve


def segments(N, valves, taps):
    """Valve segments: components of the network after cutting every valve edge.

    A valve edge belongs half to each side.  Returns segment id per vertex, and per segment
    its length, customers and the valve edges that close it.
    """
    H = N.copy()
    H.remove_edges_from([tuple(e) for e in valves])
    seg_of = {}
    for i, c in enumerate(nx.connected_components(H)):
        for v in c:
            seg_of[v] = i
    info = defaultdict(lambda: dict(length=0.0, customers=0, valves=set(), nodes=set()))
    for a, b, d in N.edges(data=True):
        if frozenset((a, b)) in valves:
            for x in (a, b):
                info[seg_of[x]]["length"] += d["length"] / 2
                info[seg_of[x]]["valves"].add(frozenset((a, b)))
        else:
            info[seg_of[a]]["length"] += d["length"]
    for v, sid in seg_of.items():
        info[sid]["nodes"].add(v)
        info[sid]["customers"] += taps.get(v, 0)
    return seg_of, dict(info)


def segment_graph(N, seg_of, valves, sources):
    """Segment graph: one vertex per valve segment, an edge per valve between two segments, and a
    super source "w" joined to the segments that hold an inflow."""
    S = nx.Graph()
    S.add_nodes_from(set(seg_of.values()))
    for e in valves:
        a, b = tuple(e)
        if seg_of[a] != seg_of[b]:
            S.add_edge(seg_of[a], seg_of[b])
    for s in sources:
        S.add_edge("w", seg_of[s])
    return S


def predicted_cut_off(S, segs):
    """Customers cut off (outside the closed segment) when each segment is closed, read from the
    segment graph: those of the segments that lose every path to "w" (Proposition: segments)."""
    out = {}
    cut_vertices = set(nx.articulation_points(S))
    for sid in segs:
        if sid not in cut_vertices:
            out[sid] = 0
            continue
        H = S.copy()
        H.remove_node(sid)
        alive = nx.node_connected_component(H, "w")
        out[sid] = sum(segs[t]["customers"] for t in H if t not in alive)
    return out


def n_minus_one(hyd, segs, taps, select, log=print):
    """Pressure-driven N-1 over the selected valve segments: demand not served and customers
    cut off (outside the closed segment) when each segment is closed."""
    base = hyd.pda()
    out = []
    items = [(sid, s) for sid, s in segs.items() if select(sid, s)]
    for k, (sid, s) in enumerate(items):
        r = hyd.pda(hyd.closed(nodes=s["nodes"]))
        in_seg = s["customers"]
        lost = sum(taps.get(v, 0) for v in hyd.N if v not in s["nodes"] and v not in r["keep"])
        out.append(dict(segment=sid, length=s["length"], customers_in_segment=in_seg, customers_cut_off=lost,
                        served_fraction=r["served_total"] / base["required_total"]))
        if (k + 1) % 50 == 0:
            log(f"    N-1: {k + 1}/{len(items)}")
    return out


# -- District Metered Areas --------------------------------------------------------------
def supply_forest(N, sources):
    """Peak-hour supply forest: every vertex takes as parent the neighbour that feeds it the largest
    flow.  Vertices are attached in order of decreasing head and a parent is always attached
    earlier (heads fall strictly along a positive flow), so the parent pointers form a forest rooted
    at the inflows.  Flows are read from the sized network (flow_ls, signed along N.edges())."""
    inflow = defaultdict(dict)  # inflow[v][u] = flow from u into v (m3/s)
    for a, b, d in N.edges(data=True):
        q = d.get("flow_ls", 0.0) / 1000
        inflow[b][a] = q
        inflow[a][b] = -q
    parent = {s: None for s in sources}
    pending = sorted((v for v in N if v not in parent), key=lambda v: -N.nodes[v]["head_m"])
    while pending:
        rest = []
        for v in pending:
            cand = [u for u in N.neighbors(v) if u in parent]
            if cand:
                parent[v] = max(cand, key=lambda u: inflow[v].get(u, 0.0))
            else:
                rest.append(v)
        if len(rest) == len(pending):
            raise ValueError("vertices not joined to an inflow")
        pending = rest
    return parent


def tree_partition(parent, weight, cap):
    """Kundu-Misra bottom-up partition of a vertex-weighted forest into the fewest connected parts of
    weight at most cap (every single weight <= cap): at each vertex, while the weight of the vertex
    and its uncut child subtrees exceeds cap, the heaviest child subtree is cut off as a part.
    Returns part id per vertex; the roots' parts hold the inflows."""
    children = defaultdict(list)
    for v, p in parent.items():
        if p is not None:
            children[p].append(v)
    order, stack = [], [v for v, p in parent.items() if p is None]
    while stack:  # pre-order; reversed it is a post-order
        v = stack.pop()
        order.append(v)
        stack.extend(children[v])
    load, cut = {}, set()
    for v in reversed(order):
        kids = sorted(children[v], key=lambda c: -load[c])
        total = weight.get(v, 0) + sum(load[c] for c in kids)
        for c in kids:
            if total <= cap:
                break
            cut.add(c)
            total -= load[c]
        load[v] = total
    part, nxt = {}, 0
    for v in order:  # a vertex starts a part if it is a root or was cut, else joins its parent's
        if parent[v] is None or v in cut:
            part[v] = nxt
            nxt += 1
        else:
            part[v] = part[parent[v]]
    return part


def dmas(N, sources, taps, target, hyd, log=print):
    """DMAs as the fewest connected parts of the peak-hour supply forest with at most `target`
    connections (Kundu-Misra); open inlets = a maximum spanning tree of the quotient graph rooted at
    the DMAs that hold an inflow (largest pipe on each tree edge), every other boundary pipe closed;
    inlets are added while a DMA lacks pressure."""
    parent = supply_forest(N, sources)
    dma = tree_partition(parent, taps, target)
    # quotient graph
    Qg = nx.Graph()
    Qg.add_nodes_from(set(dma.values()))
    boundary = defaultdict(list)
    for a, b, d in N.edges(data=True):
        if dma[a] != dma[b]:
            key = tuple(sorted((dma[a], dma[b])))
            boundary[key].append((a, b))
    for key, es in boundary.items():
        best = max(es, key=lambda e: N.edges[e]["diameter_mm"])
        Qg.add_edge(*key, w=-N.edges[best]["diameter_mm"], pipe=best)
    roots = {dma[s] for s in sources}
    # maximum-capacity spanning forest rooted at the source DMAs: join the roots by a super node
    T = Qg.copy()
    T.add_node("root")
    for r in roots:
        T.add_edge("root", r, w=-1e9)
    F = nx.minimum_spanning_tree(T, weight="w")
    inlets = {tuple(sorted((a, b), key=str)): F.edges[a, b]["pipe"] for a, b in F.edges() if "root" not in (a, b)}
    closed = [e for key, es in boundary.items() for e in es if e not in inlets.values()]

    def evaluate(closed_edges):
        return hyd.pda(hyd.closed(edges=closed_edges))

    res = evaluate(closed)
    added = 0
    for _ in range(30):
        low = [v for v in res["keep"] if hyd.q_req.get(v, 0) > 0 and res["pressure"][v] < hyd.p_req - 1e-6]
        lost = [v for v in N if hyd.q_req.get(v, 0) > 0 and v not in res["keep"]]
        bad = low + lost
        if not bad:
            break
        # open the closed boundary pipe that is largest among those touching the worst DMA
        worst = max(set(dma[v] for v in bad), key=lambda d_: sum(1 for v in bad if dma[v] == d_))
        cand = [e for e in closed if worst in (dma[e[0]], dma[e[1]])]
        if not cand:
            break
        e = max(cand, key=lambda e: N.edges[e]["diameter_mm"])
        closed.remove(e)
        inlets[("extra", added)] = e
        added += 1
        res = evaluate(closed)
    sizes = defaultdict(int)
    for v, d_ in dma.items():
        sizes[d_] += taps.get(v, 0)
    log(f"  DMAs: {len(set(dma.values()))}, inlets {len(inlets)}, closed boundary pipes {len(closed)}")
    return dict(dma=dma, inlets=list(inlets.values()), closed=closed, quotient=Qg, forest=F, roots=roots,
                sizes=dict(sizes), extra_inlets=added, served_fraction=res["served_total"] / res["required_total"],
                min_pressure=min(res["pressure"][v] for v in res["keep"] if hyd.q_req.get(v, 0) > 0))


# -- fire flow at every hydrant -------------------------------------------------------------
def fire_all(hyd, hydrant_nodes, q_fire, log=print):
    out = []
    for k, h in enumerate(hydrant_nodes):
        p, flows, keep, dem = hyd.dda(extra={h: q_fire})
        others = [p[v] for v in keep if hyd.q_req.get(v, 0) > 0 and v != h]
        out.append(dict(node=h, residual=p[h], min_other=min(others)))
        if (k + 1) % 200 == 0:
            log(f"    fire: {k + 1}/{len(hydrant_nodes)}")
    return out


def fire_capacity(hyd, h, p_res, q_hi=0.1, tol=5e-4):
    """Largest flow (m3/s) that hydrant h delivers on top of the peak demand with a residual pressure of
    at least p_res there.  The residual decreases with the drawn flow (Proposition: monotone heads),
    so bisection applies.  Returns q_hi if the hydrant still holds p_res at q_hi."""
    res = lambda q: hyd.dda(extra={h: q})[0][h]
    if res(0.0) < p_res:
        return 0.0
    if res(q_hi) >= p_res:
        return q_hi
    lo, hi = 0.0, q_hi
    while hi - lo > tol:
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if res(mid) >= p_res else (lo, mid)
    return lo


def hydrant_nodes(N, comps):
    nodes = []
    for ri, s, (a, b), off in comps["hydrants"]:
        nodes.append(a if off < N.edges[a, b]["length"] / 2 else b)
    return list(dict.fromkeys(nodes))


# -- bill of quantities -----------------------------------------------------------------------
def boq(N, packs, comps):
    w = packs["water"]
    costs = w["materials"]["costs_per_m"][w["materials"]["default"]]
    known = sorted((int(k), v) for k, v in costs.items())
    unit = lambda d: np.interp(d, [k for k, _ in known], [v for _, v in known])
    by_d = defaultdict(float)
    for a, b, d in N.edges(data=True):
        by_d[d["diameter_mm"]] += d["length"]
    rows = [dict(diameter_mm=D, length_m=L, unit_cost=float(unit(D)), cost=float(unit(D) * L)) for D, L in sorted(by_d.items())]
    pipe = sum(r["cost"] for r in rows)
    trench = sum(by_d.values()) * w["materials"]["trench_cost_per_m"]
    return dict(rows=rows, pipe_cost=pipe, trench_cost=trench, total=pipe + trench,
                n_valves=len(comps["valves"]) + len(comps["branch_valves"]), n_hydrants=len(comps["hydrants"]))


# -- EPANET export and cross-check --------------------------------------------------------
def write_inp(N, sources, packs, path, demand=None):
    cap = packs["water"]["capacity"]
    ids = {v: f"J{i}" for i, v in enumerate(N)}
    for i, s in enumerate(sources):
        ids[s] = f"R{i}"
    lines = ["[TITLE]", "Gombe water network (engine export)", "", "[JUNCTIONS]", ";ID Elev Demand"]
    for v in N:
        if v in sources:
            continue
        q = (demand or {}).get(v, N.nodes[v].get("q_req", 0.0)) * 1000
        lines.append(f"{ids[v]} {N.nodes[v]['z']:.4f} {q:.6f}")
    lines += ["", "[RESERVOIRS]", ";ID Head"]
    for s in sources:
        lines.append(f"{ids[s]} {N.nodes[s]['z'] + cap['source_pressure_mwc']:.4f}")
    lines += ["", "[PIPES]", ";ID N1 N2 Length Diameter Roughness MinorLoss Status"]
    edges = list(N.edges())
    for k, (a, b) in enumerate(edges):
        d = N.edges[a, b]
        lines.append(f"P{k} {ids[a]} {ids[b]} {d['length']:.4f} {d['diameter_mm']:.1f} {cap['hazen_williams_c']} 0 Open")
    lines += ["", "[OPTIONS]", "Units LPS", "Headloss H-W", "Accuracy 0.000001", "Trials 200", "",
              "[TIMES]", "Duration 0", "", "[REPORT]", "Status No", "", "[END]"]
    Path(path).write_text("\n".join(lines))
    return ids, edges


def epanet_check(N, sources, packs, workdir):
    """Run EPANET 2.2 (WNTR) on the exported network and compare heads and flows with the GGA."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    inp = workdir / "gombe_water.inp"
    ids, edges = write_inp(N, sources, packs, inp)
    script = ROOT / "dev" / "epanet_run.py"
    script.write_text(
        "import sys, wntr\n"
        "wn = wntr.network.WaterNetworkModel(sys.argv[1])\n"
        "res = wntr.sim.EpanetSimulator(wn).run_sim(file_prefix=sys.argv[2])\n"
        "res.node['head'].iloc[0].to_csv(sys.argv[2] + '_head.csv')\n"
        "res.link['flowrate'].iloc[0].to_csv(sys.argv[2] + '_flow.csv')\n")
    prefix = str((workdir / "epanet").resolve())
    pr = subprocess.run([str(EPANET_PY), str(script), str(inp.resolve()), prefix], capture_output=True, text=True,
                        cwd=str(workdir.resolve()))
    if pr.returncode != 0:
        raise RuntimeError(pr.stderr[-2000:])
    import pandas as pd
    h = pd.read_csv(prefix + "_head.csv", index_col=0).iloc[:, 0]
    f = pd.read_csv(prefix + "_flow.csv", index_col=0).iloc[:, 0]
    # pipes are exported in the order and orientation of N.edges(), the order used by the GGA,
    # so EPANET's flow from N1 to N2 has the sign of the engine's flow
    dh = [abs(h[ids[v]] - N.nodes[v]["head_m"]) for v in N]
    df = [abs(f[f"P{k}"] - N.edges[a, b]["flow_ls"] / 1000) for k, (a, b) in enumerate(edges)]
    return dict(max_head_diff_m=float(max(dh)), mean_head_diff_m=float(np.mean(dh)),
                max_flow_diff_ls=float(max(df) * 1000), n_nodes=len(dh), n_pipes=len(df),
                head_engine=[N.nodes[v]["head_m"] for v in N], head_epanet=[float(h[ids[v]]) for v in N])
