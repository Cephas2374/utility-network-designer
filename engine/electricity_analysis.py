"""
Analysis of the designed electricity network (paper P3): monotone diversified loads, LV losses and
cable quantities, transformer outages restored through the normally-open LV ties, 20 kV open rings
(open point, losses, N-1), trench cuts and ring re-routing, and load growth.
"""
import copy
import math
from collections import defaultdict, deque

import networkx as nx

from . import design
from .electricity import _g, _service_drop, _size_cell


# -- loads, as in design_electricity ------------------------------------------------------
def loads(G, rows, bldg, packs):
    """LV loads per tap vertex (P kW, units), service lengths, and the taps of MV customers."""
    e = packs["electricity"]
    pf = e["capacity"]["power_factor"]
    agg = {}
    for r in rows:
        b = r["bldg_id"]
        t = G.graph["tap_vertex"][b]
        d = agg.setdefault(t, {"P": 0.0, "units": 0, "svc": 0.0})
        d["P"] += float(bldg.loc[b, "elec_kw"])
        d["units"] += int(bldg.loc[b, "units"])
        d["svc"] += float(bldg.loc[b, "elec_service_kw"]) / pf
    mv = {t for t, d in agg.items() if d["svc"] >= e["capacity"]["mv_customer_min_kva"]}
    demand = {t: {"P": d["P"], "units": d["units"]} for t, d in agg.items() if t not in mv}
    svc_len = {G.graph["tap_vertex"][r["bldg_id"]]: r["dist"] for r in rows if not r.get("host")}
    return demand, svc_len, mv


class LV:
    """Electrical model of the LV forest: currents and voltage drops on any tree of it."""

    def __init__(self, F, er, demand, svc_len, packs):
        self.F, self.er, self.demand, self.svc_len = F, er, demand, svc_len
        self.e = packs["electricity"]
        cap = self.e["capacity"]
        self.V, self.pf = cap["nominal_voltage_v"], cap["power_factor"]
        self.sin = math.sqrt(1 - self.pf ** 2)
        self.rho, self.X = cap["conductor_resistivity_ohmm"], cap["conductor_reactance_ohm_per_km"] / 1000
        self.amp = {float(k): v * cap["ampacity_derating_factor"] for k, v in self.e["cable_capacity"].items()}
        self.sizes = [a for a in self.e["standard_sizes"] if a >= cap["min_main_size_mm2"]]
        self.owner = er["owner"]
        self.sites = [s for s, t in er["transformers"].items() if "mv_customer" not in t]
        self.cell_taps = defaultdict(list)
        for t in demand:
            self.cell_taps[self.owner[t]].append(t)
        # service drops do not depend on the feeding transformer
        self.svc = {t: _service_drop(t, demand, svc_len, self.e, F.nodes[t].get("service_mm2") if t in F else None)
                    for t in demand}

    def current(self, P_kw, n, scale=1.0):
        return scale * P_kw * _g(int(n), self.e["demand"]) * 1000 / (math.sqrt(3) * self.V * self.pf)

    def r_per_m(self, A):
        return self.rho / (A * 1e-6)

    def drop_pct(self, I, L, A):
        return math.sqrt(3) * I * L * (self.r_per_m(A) * self.pf + self.X * self.sin) / self.V * 100

    def solve(self, T, root, taps, size_of, scale=1.0):
        """Currents and drops on tree T fed at root, for the loads of taps; size_of(arc) -> mm2."""
        parent, order, down = design.orient(T, [root], {t: {"P": self.demand[t]["P"], "n": self.demand[t]["units"]}
                                                       for t in taps})
        I, drop = {}, {root: 0.0}
        for v in order:
            p = parent[v]
            if p is None:
                continue
            I[(p, v)] = self.current(down[v]["P"], down[v]["n"], scale)
            drop[v] = drop[p] + self.drop_pct(I[(p, v)], T.edges[p, v]["length"], size_of((p, v)))
        cust = {t: drop[t] + scale * self.svc[t] for t in taps if t in drop}
        return dict(parent=parent, order=order, down=down, I=I, drop=drop, cust=cust)

    def cell_tree(self, s):
        return self.F.subgraph(nx.node_connected_component(self.F, s))

    def size(self, a):
        return self.F.edges[a]["size_mm2"]


# -- diversity: monotone design loads ----------------------------------------------------------
def diversity_inversions(lv):
    """Arcs whose design current is smaller than that of a child arc (stepwise diversity table)."""
    inv, n_arcs, worst = 0, 0, 0.0
    for s in lv.sites:
        if not lv.cell_taps[s] or s not in lv.F:
            continue
        r = lv.solve(lv.cell_tree(s), s, lv.cell_taps[s], lv.size)
        kids = defaultdict(list)
        for (p, v), i in r["I"].items():
            kids[p].append(i)
        for (p, v), i in r["I"].items():
            n_arcs += 1
            c = max(kids.get(v, [0.0]))
            if c > i + 1e-9:
                inv += 1
                worst = max(worst, c / i - 1)
    return dict(arcs=n_arcs, inversions=inv, worst_excess=worst)


def resize_monotone(lv, packs, bldg):
    """Re-size the LV cables of every cell with the monotone majorant of the diversity factor
    (same transformers and trees).  Returns the re-sized forest and the changes."""
    p2 = copy.deepcopy(packs)
    p2["electricity"]["demand"]["monotone_diversity"] = True
    F2 = lv.F.copy()
    for s in lv.sites:
        if lv.cell_taps[s]:
            _size_cell(F2, s, lv.cell_taps[s], lv.demand, p2["electricity"], lv.svc_len, bldg, None)
    up = sum(1 for a, b, d in F2.edges(data=True) if "size_mm2" in d and d["size_mm2"] > lv.F.edges[a, b]["size_mm2"])
    down = sum(1 for a, b, d in F2.edges(data=True) if "size_mm2" in d and d["size_mm2"] < lv.F.edges[a, b]["size_mm2"])
    return F2, dict(upsized=up, downsized=down, cost_before=cable_cost(lv.F, packs), cost_after=cable_cost(F2, packs))


def cable_cost(F, packs):
    e = packs["electricity"]
    costs = {float(k): v for k, v in e["materials"]["costs_per_m"]["XLPE_Cu"].items()}
    by = defaultdict(float)
    for a, b, d in F.edges(data=True):
        if "size_mm2" in d:
            by[d["size_mm2"]] += d["length"]
    cable = sum(costs[A] * L for A, L in by.items())
    return dict(cable=cable, trench=sum(by.values()) * e["materials"]["trench_cost_per_m"], by_size=dict(by))


def lv_losses(lv, scale=1.0):
    """Peak I^2 R losses (kW) of the LV mains at the diversified design currents."""
    loss = 0.0
    for s in lv.sites:
        if not lv.cell_taps[s] or s not in lv.F:
            continue
        r = lv.solve(lv.cell_tree(s), s, lv.cell_taps[s], lv.size, scale)
        for (p, v), i in r["I"].items():
            loss += 3 * i ** 2 * lv.r_per_m(lv.size((p, v))) * lv.F.edges[p, v]["length"]
    return loss / 1000


# -- transformer outage restored through an LV tie ----------------------------------------------
def tie_restoration(lv, packs):
    """For every distribution transformer: can its whole cell be fed from a neighbour through one
    normally-open LV tie?  Cell tree + tie + neighbour tree is a tree (Proposition), so the check of
    loading, ampacity and voltage drop is exact."""
    cap = packs["electricity"]["capacity"]
    lim_drop = cap["emergency_voltage_drop_percent"]
    lim_load = cap["emergency_transformer_loading_percent"] / 100
    trafos = lv.er["transformers"]
    nbrs = defaultdict(list)
    for (oa, ob), (a, b, d) in lv.er["ties"].items():
        x, y = (a, b) if lv.owner.get(a) == oa else (b, a)
        nbrs[oa].append((ob, x, y, d))
        nbrs[ob].append((oa, y, x, d))
    out = []
    for s in lv.sites:
        taps = lv.cell_taps[s]
        if not taps or s not in lv.F:
            continue
        best = None
        for s2, x_mine, x_theirs, d in nbrs.get(s, []):
            if s2 not in trafos or "mv_customer" in trafos[s2] or s2 not in lv.F:
                continue
            T = nx.Graph(lv.cell_tree(s))
            T.add_edges_from(lv.cell_tree(s2).edges(data=True))
            T.add_edge(x_theirs, x_mine, length=d["length"], tie=True)
            both = taps + lv.cell_taps[s2]
            kva = (sum(lv.demand[t]["P"] for t in both) * _g(sum(lv.demand[t]["units"] for t in both),
                                                             lv.e["demand"]) / lv.pf)
            # the tie cable is sized for the current it carries
            n_tie = sum(lv.demand[t]["units"] for t in taps)
            I_tie = lv.current(sum(lv.demand[t]["P"] for t in taps), n_tie)
            A_tie = next((A for A in lv.sizes if lv.amp[A] >= I_tie), None)
            size_of = lambda a: (A_tie or lv.sizes[-1]) if T.edges[a].get("tie") else T.edges[a]["size_mm2"]
            r = lv.solve(T, s2, both, size_of)
            over_amp = A_tie is None or any(i > lv.amp[size_of(a)] for a, i in r["I"].items())
            drop = max(r["cust"].values())
            load = kva / trafos[s2]["rating_kva"]
            ok = load <= lim_load and not over_amp and drop <= lim_drop
            cand = dict(neighbour=s2, ok=ok, loading=load, ampacity_ok=not over_amp, max_drop=drop)
            if best is None or (ok, -drop) > (best["ok"], -best["max_drop"]):
                best = cand
        out.append(dict(site=s, customers=len(taps), kva=trafos[s]["kva"], n_ties=len(nbrs.get(s, [])),
                        **(best or dict(neighbour=None, ok=False, loading=None, ampacity_ok=None, max_drop=None))))
    return out


def feeder_restoration(lv, packs, upsize=False, req=None):
    """Partial restoration: when transformer s fails, T_s - s splits into its feeders (one subtree per
    LV board way).  Each feeder is back-fed through a tie on it by a neighbouring transformer, if the
    neighbour's tree with all the feeders it has taken (still a tree) keeps the emergency loading,
    the ampacities and the emergency voltage drop.  Feeders are placed largest first, each on the
    first tie that stays feasible.  With upsize, cables that the back-feed would overload are given
    the smallest standard size that carries it (recorded in req: arc -> required mm2, the largest over
    all outages), and the drop is checked with those sizes."""
    req = {} if req is None else req
    cap = packs["electricity"]["capacity"]
    lim_drop = cap["emergency_voltage_drop_percent"]
    lim_load = cap["emergency_transformer_loading_percent"] / 100
    trafos = lv.er["transformers"]
    ties_at = defaultdict(list)  # vertex of the failed cell -> [(neighbour site, its tie vertex, tie data)]
    for (oa, ob), (a, b, d) in lv.er["ties"].items():
        ties_at[a].append((lv.owner.get(b), b, d))
        ties_at[b].append((lv.owner.get(a), a, d))
    out = []
    for s in lv.sites:
        taps = lv.cell_taps[s]
        if not taps or s not in lv.F:
            continue
        Ts = lv.cell_tree(s)
        H = nx.Graph(Ts)
        H.remove_node(s)
        feeders = []
        for comp in nx.connected_components(H):
            ft = [t for t in taps if t in comp]
            if ft:
                feeders.append((sum(lv.demand[t]["P"] for t in ft), comp, ft))
        feeders.sort(key=lambda x: -x[0])
        taken = defaultdict(list)  # neighbour -> [(feeder nodes, feeder taps, own vertex, their vertex, tie)]
        restored = 0
        why = defaultdict(int)  # binding constraint of every rejected (feeder, tie) trial
        for P, comp, ft in feeders:
            cands = [(s2, x, y, d) for x in comp for (s2, y, d) in ties_at.get(x, [])
                     if s2 in trafos and "mv_customer" not in trafos[s2] and s2 in lv.F and s2 != s]
            for s2, x, y, d in sorted(cands, key=lambda c: c[3]["length"]):
                trial = taken[s2] + [(comp, ft, x, y, d)]
                T = nx.Graph(lv.cell_tree(s2))
                ties = set()
                load_taps = list(lv.cell_taps[s2])
                for c_, ft_, x_, y_, d_ in trial:
                    T.add_edges_from(Ts.subgraph(c_).edges(data=True))
                    T.add_edge(y_, x_, length=d_["length"], tie=True)
                    ties.add(frozenset((x_, y_)))
                    load_taps += ft_
                if not nx.is_tree(T):
                    why["not a tree"] += 1
                    continue
                kva = (sum(lv.demand[t]["P"] for t in load_taps)
                       * _g(sum(lv.demand[t]["units"] for t in load_taps), lv.e["demand"]) / lv.pf)
                if kva > lim_load * trafos[s2]["rating_kva"]:
                    why["loading"] += 1
                    continue
                base = lambda a: lv.sizes[-1] if T.edges[a].get("tie") else T.edges[a]["size_mm2"]
                r = lv.solve(T, s2, load_taps, base)
                over = {a: i for a, i in r["I"].items() if i > lv.amp[base(a)]}
                new = {}
                if over:
                    if not upsize:
                        why["ampacity"] += 1
                        continue
                    for a, i in over.items():
                        A = next((A for A in lv.sizes if lv.amp[A] >= i), None)
                        if A is None:
                            break
                        new[a] = A
                    if len(new) < len(over):
                        why["ampacity"] += 1
                        continue
                    r = lv.solve(T, s2, load_taps, lambda a: new.get(a, new.get(a[::-1], base(a))))
                if max(r["cust"].values()) > lim_drop:
                    why["drop"] += 1
                    continue
                for a, A in new.items():
                    if not T.edges[a].get("tie"):
                        k = frozenset(a)
                        req[k] = max(req.get(k, 0), A)
                taken[s2] = trial
                restored += len(ft)
                break
        out.append(dict(site=s, customers=len(taps), feeders=len(feeders), restored=restored,
                        neighbours=sum(1 for v in taken.values() if v), why=dict(why),
                        feeder_kva=[P * _g(len(ft), lv.e["demand"]) / lv.pf for P, _, ft in feeders]))
    return out


# -- 20 kV open rings -----------------------------------------------------------------------------
def ring_feeders(ring, trafos, packs, nop=None, scale=1.0):
    """Currents, losses and drops of an open ring with the open point on leg `nop` (between tour[nop]
    and tour[nop+1]); legs before it feed tour[1..nop] from the infeed, legs after it the rest."""
    e = packs["electricity"]
    V, pf = e["capacity"]["mv_voltage_v"], e["capacity"]["power_factor"]
    c = e["mv_cable"]
    R, X = c["r_ohm_per_km"] / 1000, c["x_ohm_per_km"] / 1000
    tour = ring["tour"]
    n = len(tour) - 2
    S = [scale * trafos[s]["kva"] * first for s, first in visits(tour)]  # diversified load, kVA, once per substation
    L = ring["leg_m"]
    k = ring["nop_leg"] if nop is None else nop
    I = [0.0] * (n + 1)  # current on leg i (A)
    for i in range(k):  # feeder A: legs 0..k-1 serve tour[1..k]
        I[i] = sum(S[i:k]) * 1000 / (math.sqrt(3) * V)
    for i in range(k + 1, n + 1):  # feeder B: legs n..k+1 serve tour[k+1..n]
        I[i] = sum(S[k:i]) * 1000 / (math.sqrt(3) * V)
    loss = sum(3 * I[i] ** 2 * R * L[i] for i in range(n + 1)) / 1000
    du = lambda legs: sum(math.sqrt(3) * I[i] * L[i] * (R * pf + X * math.sqrt(1 - pf ** 2)) for i in legs) / V * 100
    return dict(I=I, loss_kw=loss, drop_a=du(range(k)), drop_b=du(range(k + 1, n + 1)), k=k, S=S, L=L)


def visits(tour):
    """(substation, 1 at its first visit else 0) along a ring tour.  The tours of the P1 design can
    pass a substation twice where it lies on the shortest path between two others (a tie of the
    Christofides expansion); its load counts once."""
    seen = set()
    out = []
    for s in tour[1:-1]:
        out.append((s, 0 if s in seen else 1))
        seen.add(s)
    return out


def ring_analysis(rings, trafos, packs, G, scale=1.0):
    """Per ring: open point by load balance (design) and by losses, and the N-1 case of the whole ring
    fed from one end (Proposition: open ring)."""
    e = packs["electricity"]
    amp = e["mv_cable"]["ampacity_a"]
    out = []
    for r in rings:
        rr = dict(r, leg_m=[leg_length(G, leg) for leg in r["legs"]])
        n = len(r["tour"]) - 2
        bal = ring_feeders(rr, trafos, packs, scale=scale)
        alls = [ring_feeders(rr, trafos, packs, nop=k, scale=scale) for k in range(n + 1)]
        opt = min(alls, key=lambda x: x["loss_kw"])
        full_a, full_b = alls[n], alls[0]  # open at the last leg: all fed from A; open at leg 0: all from B
        out.append(dict(n=len(set(r["tour"][1:-1])), visits=n, kva=sum(bal["S"]), length_m=sum(rr["leg_m"]),
                        nop_balance=bal["k"], nop_loss=opt["k"],
                        loss_balance_kw=bal["loss_kw"], loss_opt_kw=opt["loss_kw"],
                        I_max_normal=max(bal["I"]), I_full=max(full_a["I"]),
                        drop_normal=max(bal["drop_a"], bal["drop_b"]),
                        drop_full=max(full_a["drop_a"], full_b["drop_b"]),
                        n1_ok=max(full_a["I"]) <= amp and max(full_a["drop_a"], full_b["drop_b"])
                        <= e["capacity"]["max_mv_voltage_drop_percent"]))
    return out


def leg_length(G, leg):
    return sum(G.edges[a, b]["length"] for a, b in zip(leg[:-1], leg[1:]))


# -- trench cuts ------------------------------------------------------------------------------------
def trench_cuts(ring, trafos, G, infeed):
    """Cut every street edge of the ring in turn (all ring cables in it fail, the open point closed):
    the substations no longer joined to the infeed, and the length-weighted expectation."""
    occ = defaultdict(int)
    adj = defaultdict(list)  # multigraph of cable runs: vertex -> [(neighbour, street edge)]
    for leg in ring["legs"]:
        for a, b in zip(leg[:-1], leg[1:]):
            f = frozenset((a, b))
            occ[f] += 1
            adj[a].append((b, f))
            adj[b].append((a, f))
    subs = list(dict.fromkeys(ring["tour"][1:-1]))
    res = {}
    for f in occ:
        seen, q = {infeed}, deque([infeed])
        while q:
            x = q.popleft()
            for y, g in adj[x]:
                if g != f and y not in seen:
                    seen.add(y)
                    q.append(y)
        lost = [s for s in subs if s not in seen]
        res[f] = (len(lost), sum(trafos[s]["kva"] for s in lost))
    Ls = {f: G.edges[tuple(f)]["length"] for f in occ}
    tot = sum(Ls.values())
    return dict(street_m=tot, cable_m=sum(Ls[f] * occ[f] for f in occ),
                shared_m=sum(Ls[f] for f in occ if occ[f] > 1),
                exposed_m=sum(Ls[f] for f in occ if res[f][0] > 0),
                exp_lost_kva=sum(Ls[f] * res[f][1] for f in occ) / tot,
                exp_lost_subs=sum(Ls[f] * res[f][0] for f in occ) / tot,
                max_lost_subs=max(v[0] for v in res.values()), per_edge=res)


def spokes(G, infeed, first, last):
    """Cheapest pair of edge-disjoint street paths infeed -> first and infeed -> last: a min-cost flow
    of value 2 from the infeed to a sink joined to first and last (Suurballe).  None if the street
    graph has no such pair (a bridge separates them from the infeed, Menger)."""
    D = nx.DiGraph()
    for a, b, d in G.edges(data=True):
        w = int(round(d["length"] * 100))
        D.add_edge(a, b, capacity=1, weight=w)
        D.add_edge(b, a, capacity=1, weight=w)
    D.add_edge(first, "_sink", capacity=1, weight=0)
    D.add_edge(last, "_sink", capacity=1, weight=0)
    D.nodes[infeed]["demand"] = -2
    D.nodes["_sink"]["demand"] = 2
    try:
        flow = nx.network_simplex(D)[1]
    except nx.NetworkXUnfeasible:
        return None
    # path decomposition: follow unused arcs with flow from the infeed to the sink, twice
    succ = {a: [b for b, f in nbr.items() if f > 0] for a, nbr in flow.items()}
    paths = {}
    for _ in range(2):
        p, x = [infeed], infeed
        while x != "_sink":
            x = succ[x].pop()
            p.append(x)
        paths[p[-2]] = p[:-1]
    return [paths[first], paths[last]]


def two_edge_component(G, infeed):
    """Vertices joined to the infeed by two edge-disjoint paths (Menger): the infeed's component of G
    without its bridges."""
    H = G.copy()
    H.remove_edges_from(nx.bridges(G))
    return nx.node_connected_component(H, infeed)


def second_infeed(G, infeed, cands, subs, trafos):
    """Bridge bound with a second infeed: a substation can be on a trench-safe ring if it has two
    edge-disjoint paths to either infeed (Menger on a super-source joined to both).  Every candidate is
    tried; returns the one that leaves the fewest substations exposed, and that count and load."""
    best = None
    for c in cands:
        if c == infeed:
            continue
        H = G.copy()
        H.add_edge("_src", infeed, length=0.0)
        H.add_edge("_src", c, length=0.0)
        safe = two_edge_component(H, "_src")
        lost = [s for s in subs if s not in safe]
        key = (len(lost), sum(trafos[s]["kva"] for s in lost))
        if best is None or key < best[0]:
            best = (key, c)
    (n, kva), c = best
    return dict(site=c, behind_bridge=n, kva_behind=kva, candidates=len(cands) - 1)


def reroute(ring, G, trafos, infeed=None, penalty=10.0, inner=True):
    """The same tour with the two spokes from the infeed laid as the cheapest edge-disjoint pair
    (spokes), then the legs between substations laid one by one on shortest street paths that avoid
    the streets already used by this ring; a street is reused only if no path avoids it (then at
    `penalty` times its length).  The open point is put back at the load balance."""
    used = set()
    tour = ring["tour"]
    sp = spokes(G, tour[0], tour[1], tour[-2]) if infeed is not None and len(tour) > 3 else None
    if sp:
        used |= {frozenset(e) for p in sp for e in zip(p[:-1], p[1:])}

    def w_avoid(a, b, d):
        return None if frozenset((a, b)) in used else d["length"]

    def w_pen(a, b, d):
        return d["length"] * (1 + penalty) if frozenset((a, b)) in used else d["length"]

    legs = []
    n = len(tour) - 1
    for i, (a, b) in enumerate(zip(tour[:-1], tour[1:])):
        if sp and i == 0:
            p = sp[0]
        elif sp and i == n - 1:
            p = sp[1][::-1]
        elif not inner:  # spokes only: keep the design's leg
            p = ring["legs"][i]
        else:
            try:
                p = nx.shortest_path(G, a, b, weight=w_avoid)
            except nx.NetworkXNoPath:
                p = nx.shortest_path(G, a, b, weight=w_pen)
        legs.append(p)
        used |= {frozenset(e) for e in zip(p[:-1], p[1:])}
    tot = sum(trafos[s]["kva"] * f for s, f in visits(tour))
    acc, nop = 0.0, 0
    for j, (s, f) in enumerate(visits(tour)):
        acc += trafos[s]["kva"] * f
        if acc >= tot / 2:
            nop = j + 1
            break
    length = sum(leg_length(G, leg) for leg in legs)
    return dict(ring, legs=legs, length_m=length, nop_leg=nop)


# -- load growth ------------------------------------------------------------------------------------
def growth_headroom(lv, packs, rings_out):
    """Load factor at which each element first breaks a rule.  Currents scale linearly with the load
    (diversity depends on the unit count only), so every drop and loading does as well."""
    e = packs["electricity"]
    lim = e["capacity"]["max_voltage_drop_percent"]
    load_cap = e["validation"]["max_transformer_loading_percent"] / 100
    trafos = lv.er["transformers"]
    h_tr = [load_cap * t["rating_kva"] / t["kva"] for s, t in trafos.items() if "mv_customer" not in t and t["kva"] > 0]
    h_tr_rating = [t["rating_kva"] / t["kva"] for s, t in trafos.items() if "mv_customer" not in t and t["kva"] > 0]
    h_drop = [lim / max(1e-9, max(F_drop for F_drop in [lv.F.nodes[t].get("drop_pct", 0.0) for t in lv.cell_taps[s]]))
              for s in lv.sites if lv.cell_taps[s]]
    amp = e["mv_cable"]["ampacity_a"]
    h_mv = [amp / r["I_full"] for r in rings_out]
    h_mv_drop = [e["capacity"]["max_mv_voltage_drop_percent"] / r["drop_full"] for r in rings_out]
    return dict(transformer_80=h_tr, transformer_100=h_tr_rating, lv_drop=h_drop, mv_n1_current=h_mv,
                mv_n1_drop=h_mv_drop)
