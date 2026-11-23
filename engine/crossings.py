"""
Corridor layout and crossings (rules C1-C6).

Drawing convention
  * Every main of utility u runs on one side of each street edge e (sigma(u, e) = +1 left of
    the edge direction, -1 right), or on both sides (dual main), in the lane of u:
        offset = w(e)/2 + kerb offset + cumulative separations in the kerb order
    of the rule pack (gas, water, electricity from the kerb outwards since P4; P1 used gas,
    electricity, water).  Mains are parallel to the axis, so
    two mains never cross along an edge.
  * A service leaves its main perpendicularly at the tap; if the building is on the other
    side it crosses the carriageway perpendicularly (C1 by construction).
  * At a junction v with incident half-edges h_0..h_{k-1} in counter-clockwise order, the
    corner sectors S_0..S_{k-1} lie between consecutive half-edges.  A main of u on edge h_i
    lies in S_i (left of the outward direction) or in S_{i-1} (right).  Joining the occupied
    sectors O_u needs a crossing leg over every radius h_j except those of one gap: if the
    longest run of unoccupied sectors (possibly empty) has g sectors, u crosses k - (g + 1)
    carriageways at v, the minimum possible (Proposition in P1).  Legs are drawn at the mouth
    of the crossed street, where only the mains parallel to that street are met.
  * A leg or service of u crosses the main of w iff w lies between u and the point it heads
    for on the cross-section (interval stabbing).
"""
from collections import defaultdict


from .config import separation, vertical_clearance

UTILS = ("gas", "electricity", "water")  # the utilities (the lane order is the rule pack's kerb_order)


def base_vertex(x):
    if isinstance(x, tuple) and x[0] == "x":
        return x[1]
    return x


class Corridor:
    def __init__(self, street, Gc, nets, taps, packs, cr, served_by=None):
        """nets: utility -> list of network graphs (their edges carry 'eid')."""
        self.S, self.Gc, self.packs, self.cr = street, Gc, packs, cr
        self.rc = cr["road_classes"]
        lay = cr["corridor_layout"]
        self.kerb = lay["kerb_order"]
        # lanes: cumulative pairwise separations (consecutive pairs) from the kerb offset
        self.lane = {}
        off = lay["min_kerb_offset_m"]
        for i, u in enumerate(self.kerb):
            if i:
                off += separation(packs, self.kerb[i - 1], u)
            self.lane[u] = off
        self.order = {u: i for i, u in enumerate(self.kerb)}
        self.kerb_off = lay["min_kerb_offset_m"]
        self.ccost = cr["crossing_cost_equivalent_m"]
        self.w_inter = cr["optimization"]["inter_utility_crossing_weight"]
        self.pen = cr["optimization"]["penalties"]

        # street edges used by each utility and the junction ports
        self.E = {u: set() for u in nets}
        self.port = {u: defaultdict(set) for u in nets}  # u -> v -> {eid}
        for u, graphs in nets.items():
            for N in graphs:
                for a, b, d in N.edges(data=True):
                    e = d.get("eid")
                    if e is None:
                        continue
                    self.E[u].add(e)
                    st, en = street.edges[e][3]["start"], street.edges[e][3]["end"]
                    for x in (a, b):
                        bx = base_vertex(x)
                        if bx == st or bx == en:
                            self.port[u][bx].add(e)
        self.U = list(nets)
        # services: building taps per edge (a shared service follows its host)
        self.byid = {r["bldg_id"]: r for r in taps}
        self.taps = defaultdict(list)
        self.hosted = [r for r in taps if r.get("host")]
        for r in taps:
            if not r.get("host"):
                self.taps[r["eid"]].append(r)
        self.hosted_on = defaultdict(list)  # eid of the root host -> (root row, extra branch length)
        for r in self.hosted:
            extra, h = r["dist"], self.byid[r["host"]]
            while h.get("host"):
                extra += h["dist"]
                h = self.byid[h["host"]]
            self.hosted_on[h["eid"]].append((h, extra))
        self.served_by = served_by or {u: None for u in self.U}
        # rotation: half-edges at each vertex
        self.rot = {v: [(eid, frm) for _, eid, frm in hs] for v, hs in street.rotation().items()}
        self.inc = defaultdict(set)
        for v, hs in self.rot.items():
            for eid, _ in hs:
                self.inc[eid].add(v)
        # initial sides: majority side of the buildings on the edge
        self.sigma = {u: {} for u in self.U}
        self.dual = {u: {} for u in self.U}
        for u in self.U:
            for e in self.E[u]:
                sides = [r["side"] for r in self.taps.get(e, [])]
                self.sigma[u][e] = 1 if sum(sides) >= 0 else -1
                self.dual[u][e] = False

    # -- geometry helpers -------------------------------------------------------
    def width(self, e):
        return self.rc.get(self.S.edges[e][3]["road_class"], self.rc["residential"])["width_m"]

    def cls(self, e):
        c = self.S.edges[e][3]["road_class"]
        return c if c in self.ccost else "residential"

    def sides_of(self, u, e):
        if e not in self.E[u]:
            return set()
        return {1, -1} if self.dual[u][e] else {self.sigma[u][e]}

    def present(self, e, side):
        """Utilities with a main on this side of e, in kerb order."""
        return [w for w in self.kerb if side in self.sides_of(w, e)]

    def offset(self, u, e, side):
        """Distance of u's main from the street axis: compact lanes of the utilities present on
        that side, kerb offset plus the pairwise separations of consecutive present utilities."""
        pres = self.present(e, side)
        off = self.width(e) / 2 + self.kerb_off
        for a, b in zip(pres[:-1], pres[1:]):
            if a == u:
                break
            off += separation(self.packs, a, b)
        return off

    def lane_need(self, e, side):
        """Verge width used by the utilities on one side (outermost offset minus half-width)."""
        pres = self.present(e, side)
        if not pres:
            return 0.0
        return self.offset(pres[-1], e, side) - self.width(e) / 2

    def verge(self, e):
        return self.rc.get(self.S.edges[e][3]["road_class"], self.rc["residential"])["verge_m"]

    # -- cross-section cost of an edge (services) ---------------------------------
    def section(self, e, record=None, penalty_only=False, counts=None):
        cost = pen = 0.0
        for r in self.taps.get(e, []):
            sb, dist = r["side"], r["dist"]
            for u in self.U:
                if e not in self.E[u]:
                    continue
                su = sb if sb in self.sides_of(u, e) else -sb
                road = su != sb
                if road:
                    cost += self.ccost[self.cls(e)]
                # utilities crossed by the service of u (interval on the signed cross-section)
                for w in self.U:
                    if w == u:
                        continue
                    for sw in self.sides_of(w, e):
                        if road:
                            hit = (sw == su and self.order[w] < self.order[u]) or sw == sb
                        else:
                            hit = sw == su and self.order[w] > self.order[u] and self.offset(w, e, sw) < dist
                        if hit:
                            cost += self.w_inter
                            if record is not None:
                                record["inter"][tuple(sorted((u, w)))] += 1
                                record["inter_service"] += 1
                # G9 service length (soft rule)
                L = self._svc_len(u, r, su)
                if L > self.packs[u]["routing"]["max_service_distance"]:
                    pen += self.pen["service_over_limit"]
                    if counts is not None:
                        counts["service"][u] += 1
                # C5 gas clearance to the building it runs along
                if u == "gas" and not road and dist - self.offset("gas", e, su) < \
                        self.packs["gas"]["constraints"]["min_separation_from_buildings_m"]:
                    pen += self.pen["gas_building_clearance"]
                    if counts is not None:
                        counts["gas"] += 1
                if record is not None:
                    record["services"][u] += 1
                    if road:
                        record["road"][u][self.cls(e)] += 1
                        record["road_service"][u] += 1
                        self._junction_distance(e, r, record, u)
        for h, extra in self.hosted_on.get(e, []):
            for u in self.U:
                if e in self.E[u]:
                    su = h["side"] if h["side"] in self.sides_of(u, e) else -h["side"]
                    if self._svc_len(u, h, su) + extra > self.packs[u]["routing"]["max_service_distance"]:
                        pen += self.pen["service_over_limit"]
                        if counts is not None:
                            counts["service"][u] += 1
        for side in (1, -1):
            if self.lane_need(e, side) > self.verge(e) + 1e-9:
                pen += self.pen["lane_outside_verge"]
                if counts is not None:
                    counts["verge"] += 1
        if penalty_only:
            return pen
        cost += pen
        for u in self.U:
            if e in self.E[u] and self.dual[u][e]:
                cost += self.dual_length(u, e)  # the second main, in metres of main
        return cost

    def dual_length(self, u, e):
        """Length of the second main of u on e: fed from the nearer end of the edge, it reaches the
        farthest tap on the side of the first main's opposite."""
        L = self.S.edges[e][3]["length"]
        far = [r["s"] for r in self.taps.get(e, []) if r["side"] != self.sigma[u][e]]
        if not far:
            return L
        return min(max(far), L - min(far))

    def _svc_len(self, u, r, su):
        off = self.offset(u, r["eid"], su)
        return abs(r["dist"] - off) if su == r["side"] else r["dist"] + off

    def _junction_distance(self, e, r, record, u):
        lim = self.cr["geometry"]["min_distance_from_junction_m"]
        L = self.S.edges[e][3]["length"]
        st, en = self.S.edges[e][3]["start"], self.S.edges[e][3]["end"]
        near = []
        if self.S.G.degree(st) >= 3 and r["s"] < lim:
            near.append("start")
        if self.S.G.degree(en) >= 3 and L - r["s"] < lim:
            near.append("end")
        if near:
            if L >= 2 * lim:
                record["c2_shifted"][u] += 1  # crossing leg moved to lim from the junction
            else:
                record["c2_unavoidable"][u] += 1

    # -- junction cost -----------------------------------------------------------
    def occupied(self, u, v):
        hs = self.rot.get(v, [])
        k = len(hs)
        occ = set()
        for i, (e, from_start) in enumerate(hs):
            if e not in self.port[u].get(v, ()):
                continue
            for s in self.sides_of(u, e):
                out = s if from_start else -s  # side w.r.t. the outward direction from v
                occ.add(i if out > 0 else (i - 1) % k)
        return occ

    def junction(self, v, record=None):
        hs = self.rot.get(v, [])
        k = len(hs)
        if k < 2:
            return 0.0
        cost = 0.0
        legs = {}
        for u in self.U:
            occ = self.occupied(u, v)
            if len(occ) <= 1:
                legs[u] = set()
                continue
            radii_cost = [self.ccost[self.cls(hs[j][0])] for j in range(k)]
            # best gap: run of unoccupied sectors j..j+g-1 (g >= 0) between occupied ones; the
            # radii h_j .. h_{j+g} are not crossed.  Radius h_j separates S_{j-1} and S_j.
            best, best_keep = -1.0, None
            for i in sorted(occ):
                g = 0
                while (i + 1 + g) % k not in occ:
                    g += 1
                keep = [(i + 1 + t) % k for t in range(g + 1)]
                val = sum(radii_cost[j] for j in keep)
                if val > best:
                    best, best_keep = val, keep
            crossed = set(range(k)) - set(best_keep)
            legs[u] = crossed
            cost += sum(radii_cost[j] for j in crossed)
            if record is not None:
                record["road_junction"][u] += len(crossed)
                for j in crossed:
                    record["road"][u][self.cls(hs[j][0])] += 1
                    if "legs" in record:  # leg positions, for the vertical design of P4
                        record["legs"].append((u, v, hs[j][0]))
        # inter-utility crossings of the legs: a leg of u over radius h_j meets every main of w
        # parallel to h_j that lies between u and the kerb (w before u in the kerb order), on each side
        for u, crossed in legs.items():
            for j in crossed:
                e = hs[j][0]
                for w in self.U:
                    if w == u or self.order[w] >= self.order[u] or e not in self.port[w].get(v, ()):
                        continue
                    n = len(self.sides_of(w, e))
                    cost += self.w_inter * n
                    if record is not None:
                        record["inter"][tuple(sorted((u, w)))] += n
                        record["inter_junction"] += n
        return cost

    # -- optimisation --------------------------------------------------------------
    def local(self, e):
        return self.section(e) + sum(self.junction(v) for v in self.inc[e])

    def optimise(self, sweeps=None, log=print):
        sweeps = sweeps or self.cr["optimization"]["side_assignment_iterations"]
        hist = [self.total()]
        for it in range(sweeps):
            changed = 0
            for u in self.U:
                allow_dual = lambda e: self.rc.get(self.S.edges[e][3]["road_class"], {}).get("dual_mains", False)
                for e in sorted(self.E[u]):
                    cur = (self.sigma[u][e], self.dual[u][e])
                    opts = [(1, False), (-1, False)] + ([(cur[0], True)] if allow_dual(e) else [])
                    best, best_c = cur, None
                    for o in opts:
                        self.sigma[u][e], self.dual[u][e] = o
                        c = self.local(e)
                        if best_c is None or c < best_c - 1e-9:
                            best, best_c = o, c
                    self.sigma[u][e], self.dual[u][e] = best
                    changed += best != cur
            hist.append(self.total())
            log(f"  side assignment sweep {it + 1}: cost {hist[-1]:.0f}, changes {changed}")
            if not changed:
                break
        return hist

    def violations(self, e):
        """Rule violations on edge e, in priority order: (gas clearance, verge, services)."""
        c = {"gas": 0, "verge": 0, "service": defaultdict(int)}
        self.section(e, penalty_only=True, counts=c)
        return (c["gas"], c["verge"], sum(c["service"].values())), dict(c["service"])

    OPTIONS = ((1, False), (-1, False), (1, True), (-1, True))

    def best_choice(self, e):
        """Lexicographically best joint side / dual choice of the utilities on e:
        fewest gas clearance violations, then verge overflows, then long services, then cost.
        Leaves the state unchanged; returns (choice, violations, service counts)."""
        import itertools
        us = [u for u in self.U if e in self.E[u]]
        cur = {u: (self.sigma[u][e], self.dual[u][e]) for u in us}
        best = None
        for combo in itertools.product(self.OPTIONS, repeat=len(us)):
            for u, o in zip(us, combo):
                self.sigma[u][e], self.dual[u][e] = o
            v, svc = self.violations(e)
            key = (v, self.local(e))
            if best is None or key < best[0]:
                best = (key, dict(zip(us, combo)), v, svc)
        for u, o in cur.items():
            self.sigma[u][e], self.dual[u][e] = o
        return best[1], best[2], best[3]

    def joint_repair(self, log=print):
        """Edges that break a rule get their lexicographically best joint choice."""
        edges = [e for e in set(self.taps).union(*self.E.values()) if any(self.violations(e)[0])]
        improved = 0
        for e in edges:
            before = self.violations(e)[0]
            choice, v, _ = self.best_choice(e)
            for u, o in choice.items():
                self.sigma[u][e], self.dual[u][e] = o
            improved += v < before
        log(f"  joint repair: {len(edges)} edges with violations, {improved} improved")
        return len(edges), improved

    def avoidability(self):
        """Per-edge exact test: current violations against the lexicographic minimum over all
        joint choices on the edge.  Returns the counts and the number of edges above the minimum."""
        cur = {"gas": 0, "verge": 0, "service": defaultdict(int)}
        low = {"gas": 0, "verge": 0, "service": defaultdict(int)}
        above = 0
        for e in set(self.taps).union(*self.E.values()):
            v, svc = self.violations(e)
            if not any(v):
                continue
            _, vmin, svc_min = self.best_choice(e)
            above += v > vmin
            for d, vv, ss in ((cur, v, svc), (low, vmin, svc_min)):
                d["gas"] += vv[0]
                d["verge"] += vv[1]
                for u, n in ss.items():
                    d["service"][u] += n
        f = lambda d: {"gas": d["gas"], "verge": d["verge"], "service": dict(d["service"])}
        return {"current": f(cur), "minimum": f(low), "edges_above_minimum": above}

    def total(self):
        edges = set(self.taps).union(*self.E.values())
        return sum(self.section(e) for e in edges) + sum(self.junction(v) for v in self.rot)

    def service_length(self, u, r):
        """Length of the service of u to a building: main lane to facade (perpendicular legs);
        a shared service adds its branch to the service of its host."""
        if r.get("host"):
            h = self.byid[r["host"]]
            while h.get("host"):
                h = self.byid[h["host"]]
            return self.service_length(u, h) + r["dist"]
        e = r["eid"]
        su = r["side"] if r["side"] in self.sides_of(u, e) else -r["side"]
        return self._svc_len(u, r, su)

    def gas_side_allowed(self, e, side):
        """Gas may run on a side only if every building there keeps the clearance from the gas
        lane at its innermost position (kerb offset)."""
        gmin = self.packs["gas"]["constraints"]["min_separation_from_buildings_m"]
        lane = self.width(e) / 2 + self.kerb_off
        return all(r["dist"] - lane >= gmin for r in self.taps.get(e, []) if r["side"] == side)

    def min_feasible_service_length(self, u, r):
        """Shortest service over the main positions the other rules allow: same side at the
        outer verge edge (for gas only on an allowed side), or the other side at the kerb."""
        if r.get("host"):
            h = self.byid[r["host"]]
            while h.get("host"):
                h = self.byid[h["host"]]
            return self.min_feasible_service_length(u, h) + r["dist"]
        e = r["eid"]
        w2 = self.width(e) / 2
        off_max = w2 + self.kerb_off
        for a, b in zip(self.kerb[:-1], self.kerb[1:]):
            if a == u:
                break
            off_max += separation(self.packs, a, b)
        opts = []
        if u != "gas" or self.gas_side_allowed(e, r["side"]):
            opts.append(max(0.0, r["dist"] - off_max))
        if u != "gas" or self.gas_side_allowed(e, -r["side"]):
            opts.append(r["dist"] + w2 + self.kerb_off)
        return min(opts) if opts else float("inf")

    def min_service_length(self, r):
        """Shortest service any side choice allows: main at the outer edge of the verge."""
        if r.get("host"):
            h = self.byid[r["host"]]
            while h.get("host"):
                h = self.byid[h["host"]]
            return self.min_service_length(h) + r["dist"]
        e = r["eid"]
        return max(0.0, r["dist"] - (self.width(e) / 2 + self.verge(e)))

    # -- report ----------------------------------------------------------------------
    def report(self, rail=None, water=None):
        rec = dict(services=defaultdict(int), road=defaultdict(lambda: defaultdict(int)),
                   road_service=defaultdict(int), road_junction=defaultdict(int),
                   inter=defaultdict(int), inter_service=0, inter_junction=0,
                   c2_shifted=defaultdict(int), c2_unavoidable=defaultdict(int))
        for e in set(self.taps) | set().union(*self.E.values()):
            self.section(e, rec)
        for v in self.rot:
            self.junction(v, rec)
        out = {k: (dict(v) if isinstance(v, defaultdict) else v) for k, v in rec.items()}
        out["road"] = {u: dict(d) for u, d in rec["road"].items()}
        out["inter"] = {f"{a}-{b}": n for (a, b), n in rec["inter"].items()}
        # C3: casings and trenchless crossings
        cas = defaultdict(int)
        trench = defaultdict(int)
        for u, d in rec["road"].items():
            for c, n in d.items():
                if self.rc[c]["casing_required"]:
                    cas[u] += n
                if self.rc[c]["trenchless"]:
                    trench[u] += n
        out["casings"], out["trenchless"] = dict(cas), dict(trench)
        # C4: local deepening at inter-utility crossings (verge depths of the packs)
        deep = {}
        for pair, n in rec["inter"].items():
            a, b = pair
            da = self.packs[a]["constraints"]["min_burial_depth_m"]
            db = self.packs[b]["constraints"]["min_burial_depth_m"]
            need = vertical_clearance(self.cr, a, b)
            deep[f"{a}-{b}"] = dict(crossings=n, depth_gap_m=round(abs(da - db), 3), clearance_m=need,
                                    deepening_points=n if abs(da - db) < need - 1e-9 else 0,
                                    lowering_m=round(max(0.0, need - abs(da - db)), 3))
        out["deepening"] = deep
        # C5: lanes inside the verge; gas clearance to buildings
        sides_over = 0
        for e in set().union(*self.E.values()):
            for side in (1, -1):
                sides_over += self.lane_need(e, side) > self.verge(e) + 1e-9
        out["lane_outside_verge_sides"] = sides_over
        gmin = self.packs["gas"]["constraints"]["min_separation_from_buildings_m"]
        gas_bad = gas_tight = 0
        if "gas" in self.U:
            for e in self.E["gas"]:
                for r in self.taps.get(e, []):
                    tight = r["dist"] - (self.width(e) / 2 + self.kerb_off) < gmin
                    gas_tight += tight
                    if r["side"] in self.sides_of("gas", e) and r["dist"] - self.offset("gas", e, r["side"]) < gmin:
                        gas_bad += 1
        out["gas_building_clearance_violations"] = gas_bad
        avoidable = 0
        if "gas" in self.U:
            for e in self.E["gas"]:
                for side in self.sides_of("gas", e):
                    if not self.gas_side_allowed(e, side) and self.gas_side_allowed(e, -side):
                        avoidable += 1
        out["gas_on_blocked_side_avoidable"] = avoidable
        out["gas_edges_blocked_both_sides"] = sum(
            1 for e in self.E.get("gas", ()) if not self.gas_side_allowed(e, 1) and not self.gas_side_allowed(e, -1))
        out["gas_buildings_within_clearance_of_kerb"] = gas_tight
        out["kerb_order"] = list(self.kerb)
        out["avoidability"] = self.avoidability_cache = self.avoidability()
        out["kerb_offset_m"] = self.kerb_off
        out["dual_mains"] = {u: sum(self.dual[u].values()) for u in self.U}
        out["dual_mains_length_m"] = {u: sum(self.dual_length(u, e) for e, d in self.dual[u].items() if d)
                                      for u in self.U}
        # rail and watercourse crossings of the mains
        for name, geoms in (("rail_crossings", rail), ("watercourse_crossings", water)):
            if geoms is None:
                continue
            from shapely import STRtree
            t = STRtree(geoms)
            out[name] = {u: sum(len(t.query(self.S.edge_geom(e), predicate="intersects")) > 0 for e in self.E[u])
                         for u in self.U}
        return out
