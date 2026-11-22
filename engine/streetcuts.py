"""
One excavation, three utilities (paper P4): when the street edge e is dug up, every main laid along it
is cut.  For each utility, the connections that lose supply:

  water        the valve segments that hold a pipe of e are closed (P2); lost: their connections and
               those cut off from every inflow in the segment graph
  electricity  LV: the connections below a cut arc of their transformer tree; MV: the substations a
               trench cut separates from the infeed (P3), with all their connections
  gas          LP: the connections below a cut arc of their zone tree; MP: the zones of the regulators
               below a cut arc of the medium-pressure tree

Connections are tap vertices, weighted by the buildings they serve.
"""
from collections import defaultdict, deque

from . import electricity_analysis as ea, water_analysis as wa


def _intervals(T, roots):
    """Pre-order index and subtree end of every vertex of a forest, children of the given roots."""
    pre, end = {}, {}
    order = []
    for r in roots:
        if r not in T or r in pre:
            continue
        stack = [(r, None, False)]
        while stack:
            v, p, done = stack.pop()
            if done:
                end[v] = len(order)
                continue
            pre[v] = len(order)
            order.append(v)
            stack.append((v, p, True))
            for w in T.neighbors(v):
                if w != p and w not in pre:
                    stack.append((w, v, False))
    return pre, end, order


def tree_losses(T, roots, weight):
    """For every street edge (eid) of the forest T: the tap vertices below a cut arc on it."""
    pre, end, order = _intervals(T, roots)
    by_eid = defaultdict(list)
    for a, b, d in T.edges(data=True):
        e = d.get("eid")
        if e is None or a not in pre or b not in pre:
            continue
        child = b if pre[b] > pre[a] and end[b] <= end[a] else a  # the end farther from the root
        by_eid[e].append((pre[child], end[child]))
    out = {}
    for e, iv in by_eid.items():
        iv.sort()
        lost, hi = set(), -1
        for lo, h in iv:
            if h <= hi:
                continue
            for i in range(max(lo, hi), h):
                if order[i] in weight:
                    lost.add(order[i])
            hi = max(hi, h)
        out[e] = lost
    return out


def water_losses(Nw, wi, taps):
    """Street edge -> tap vertices without water when the segments holding its pipes are closed."""
    ve = wa.valve_edges(Nw, wi["components"])
    seg_of, segs = wa.segments(Nw, ve, taps)
    Sg = wa.segment_graph(Nw, seg_of, ve, wi["sources"])
    nodes_of = defaultdict(set)
    for v, s in seg_of.items():
        nodes_of[s].add(v)
    closed_by = defaultdict(set)
    for a, b, d in Nw.edges(data=True):
        e = d.get("eid")
        if e is not None:
            closed_by[e] |= {seg_of[a], seg_of[b]} if frozenset((a, b)) in ve else {seg_of[a]}
    out = {}
    for e, closed in closed_by.items():
        seen, q = {"w"}, deque(["w"])
        while q:
            x = q.popleft()
            for y in Sg.neighbors(x):
                if y not in seen and y not in closed:
                    seen.add(y)
                    q.append(y)
        out[e] = {v for s in Sg if s != "w" and s not in seen for v in nodes_of[s] if v in taps}
    return out


def electricity_losses(G, F, er, infeed, cell_taps, mv_tap_of):
    """Street edge -> tap vertices without electricity: LV subtrees and MV trench cuts."""
    lv_sites = [s for s, t in er["transformers"].items() if "mv_customer" not in t]
    weight = {t for taps in cell_taps.values() for t in taps}
    out = defaultdict(set, tree_losses(F, lv_sites, weight))
    for r in er["rings"]:
        tc = ea.trench_cuts(r, er["transformers"], G, infeed)
        subs = list(dict.fromkeys(r["tour"][1:-1]))
        for f, (n, _) in tc["per_edge"].items():
            if not n:
                continue
            a, b = tuple(f)
            e = G.edges[a, b].get("eid")
            # recompute which substations (cheap: n > 0 only), from the multigraph of the ring
            lost = _ring_lost(r, f, infeed, subs)
            for s in lost:
                out[e] |= set(cell_taps.get(s, ())) | set(mv_tap_of.get(s, ()))
    return dict(out)


def _ring_lost(ring, f, infeed, subs):
    adj = defaultdict(list)
    for leg in ring["legs"]:
        for a, b in zip(leg[:-1], leg[1:]):
            g = frozenset((a, b))
            if g != f:
                adj[a].append(b)
                adj[b].append(a)
    seen, q = {infeed}, deque([infeed])
    while q:
        x = q.popleft()
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                q.append(y)
    return [s for s in subs if s not in seen]


def gas_losses(LP, MP, gr, gate, zone_taps):
    """Street edge -> tap vertices without gas: LP subtrees below a cut arc, and the whole zones of
    the regulators below a cut arc of the MP tree."""
    regs = list(gr["sites"])
    weight = {t for taps in zone_taps.values() for t in taps}
    out = defaultdict(set, tree_losses(LP, regs, weight))
    mp = tree_losses(MP, [gate], set(regs))
    for e, lost_regs in mp.items():
        for s in lost_regs:
            out[e] |= set(zone_taps.get(s, ()))
    return dict(out)


def street_cuts(m, demand_lv, mv_taps, gate):
    """Per street edge: the tap vertices that lose water, electricity and gas when it is dug up."""
    G, rows = m["G"], m["rows"]
    tv = G.graph["tap_vertex"]
    taps = defaultdict(int)
    for r in rows:
        taps[tv[r["bldg_id"]]] += 1
    W = water_losses(m["Nw"], m["wi"], taps)
    er = m["er"]
    cell_taps = defaultdict(list)
    for t in demand_lv:
        cell_taps[er["owner"][t]].append(t)
    mv_tap_of = {t: [t] for t in mv_taps}
    E = electricity_losses(G, m["Fe"], er, m["infeed"], cell_taps, mv_tap_of)
    zone_taps = defaultdict(list)
    for t in taps:
        z = m["gr"]["owner"].get(t)
        if z is not None:
            zone_taps[z].append(t)
    Gs = gas_losses(m["LP"], m["MP"], m["gr"], gate, zone_taps)
    return dict(water=W, electricity=E, gas=Gs, taps=dict(taps))
