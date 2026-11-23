"""
Corridor graph - service taps on the street edges and the subdivided graph
Gc on which the mains are designed.

Each served building gets one tap, preferably on its V8 access street (the street
stroke of its address): the foot of the shortest segment from the footprint to the
street edge.  When the foot lies inside the edge the segment is perpendicular to
the street, so the service is a right-angle connection (rule G9).  The building's
street side is the sign of that segment relative to the edge direction.

Gc is G0 subdivided at the taps.  Subdivision inserts degree-2 vertices only,
so Gc is plane, has the same faces as G0 and the same maximum degree.
"""
import math
from collections import defaultdict

import networkx as nx
import numpy as np
from shapely import STRtree, shortest_line
from shapely.geometry import LineString, Point
from shapely.ops import substring

END_CLEARANCE_M = 1.0  # Taps keep this distance from the edge end vertices


def _beyond(inter, point, tol=0.05):
    """True if some part of the intersection lies farther than tol from the contact point."""
    return not inter.is_empty and not inter.difference(point.buffer(tol)).is_empty


def _tangent(line, s, h=0.5):
    a = line.interpolate(max(s - h, 0.0)).coords[0]
    b = line.interpolate(min(s + h, line.length)).coords[0]
    t = np.subtract(b, a)
    n = np.linalg.norm(t)
    return t / n if n else np.array([1.0, 0.0])


class Taps:
    """One tap per served building: edge, arc position, side, distance to the facade.

    1. Admissible edges: the service segment (foot -> facade) must not cross another street
       edge, i.e. it stays in the face of the building (rule G1).  The V8 access street is
       preferred; otherwise the nearest admissible edge is used.
    2. Edges too short for their taps at the minimum spacing pass the buildings farthest
       away to their next admissible edge (rule G9).
    3. Taps are spread to the minimum spacing.
    4. Uncrossing: two crossing services on the same edge side exchange their taps.  By the
       triangle inequality each exchange shortens the total service length, so the process
       stops, and then no two services of the edge side cross.
    """

    def __init__(self, buildings, street, min_spacing=1.0, search_m=60.0):
        self.rows = []
        self.min_spacing = min_spacing
        edges = street.edges
        self.edges = edges
        by_stroke = defaultdict(list)
        for eid, (u, v, k, d) in edges.items():
            by_stroke[d["stroke_id"]].append(eid)
        eids = list(edges)
        geoms = [edges[e][3]["geometry"] for e in eids]
        self.tree, self.eids, self.geoms = STRtree(geoms), eids, geoms
        self.stats = dict(access_changed=0, non_planar=0, overflow_moved=0, uncross_swaps=0)

        for b in buildings.itertuples():
            if not b.served:
                continue
            fp = b.geometry
            stroke = sorted(by_stroke.get(b.stroke_id, []), key=lambda e: edges[e][3]["geometry"].distance(fp))
            dmin = geoms[int(self.tree.nearest(fp))].distance(fp)
            near = [eids[k] for k in self.tree.query(fp.buffer(dmin + search_m))]
            near = sorted(set(near) - set(stroke), key=lambda e: edges[e][3]["geometry"].distance(fp))
            options = [self._option(e, fp) for e in stroke + near]
            adm = [o for o in options if o["admissible"]]
            if adm:
                pick = adm[0]
                self.stats["access_changed"] += pick["eid"] not in stroke
            else:
                pick = options[0]
                self.stats["non_planar"] += 1
            self.rows.append(dict(bldg_id=b.bldg_id, **pick, alternatives=adm[1:4], fp=fp, options=adm[1:6]))
        self._overflow()
        self._space()
        self._uncross()
        self._repair()
        self._retap()
        if self.stats.get("service_crossings_left"):
            self._exchange_fixed()  # only when crossings remain, so designs without any are unchanged
        for r in self.rows:
            r.pop("alternatives", None)
            r.pop("options", None)
            r.pop("fp", None)

    def _option(self, eid, fp):
        line = self.edges[eid][3]["geometry"]
        sl = shortest_line(line, fp)
        P, Q = sl.coords[0], sl.coords[-1]
        s = line.project(Point(P))
        t = _tangent(line, s)
        w = np.subtract(Q, P)
        dist = float(np.linalg.norm(w))
        side = 1 if t[0] * w[1] - t[1] * w[0] >= 0 else -1
        perp = dist < 1e-6 or abs(np.dot(t, w)) / dist < math.sin(math.radians(10))
        adm = True
        if dist > 1e-6:
            for k in self.tree.query(sl, predicate="intersects"):
                if self.eids[k] == eid:
                    continue
                inter = sl.intersection(self.geoms[k])
                if _beyond(inter, Point(P)):
                    adm = False
                    break
        return dict(eid=eid, s=s, s0=s, side=side, dist=dist, perpendicular=bool(perp), foot=P, facade=Q,
                    admissible=adm)

    def _capacity(self, eid):
        L = self.edges[eid][3]["geometry"].length
        c = min(END_CLEARANCE_M, L / 3)
        return int((L - 2 * c) // self.min_spacing) + 1

    def _overflow(self, rounds=10):
        for _ in range(rounds):
            per = defaultdict(list)
            for r in self.rows:
                per[r["eid"]].append(r)
            moved = 0
            for eid in list(per):
                rs = per[eid]
                extra = len(rs) - self._capacity(eid)
                if extra <= 0:
                    continue
                for r in sorted(rs, key=lambda r: -r["dist"]):
                    if extra <= 0:
                        break
                    for alt in r["alternatives"]:
                        if len(per[alt["eid"]]) < self._capacity(alt["eid"]):
                            rs.remove(r)
                            per[alt["eid"]].append(r)
                            r["alternatives"] = [a for a in r["alternatives"] if a is not alt]
                            r.update({k: v for k, v in alt.items()})
                            extra -= 1
                            moved += 1
                            break
            self.stats["overflow_moved"] += moved
            if not moved:
                break

    def _space(self):
        """Distinct taps: at least min_spacing apart and END_CLEARANCE_M from the ends (G9)."""
        for eid, rs in self.by_edge(key="s0").items():
            L = self.edges[eid][3]["geometry"].length
            lo, hi = min(END_CLEARANCE_M, L / 3), L - min(END_CLEARANCE_M, L / 3)
            n = len(rs)
            step = min(self.min_spacing, (hi - lo) / max(n - 1, 1)) if n > 1 else 0
            prev = -np.inf
            for r in rs:
                r["s"] = min(max(r["s0"], lo, prev + step), hi)
                prev = r["s"]
            nxt = np.inf
            for r in reversed(rs):
                r["s"] = min(r["s"], nxt - step)
                nxt = r["s"]

    def _uncross(self, max_pass=60):
        groups = defaultdict(list)
        for r in self.rows:
            groups[(r["eid"], r["side"])].append(r)
        for (eid, side), rs in groups.items():
            if len(rs) < 2:
                continue
            line = self.edges[eid][3]["geometry"]
            seg = lambda r: LineString([line.interpolate(r["s"]).coords[0], r["facade"]])
            for _ in range(max_pass):
                segs = [seg(r) for r in rs]
                swapped = False
                for i in range(len(rs)):
                    for j in range(i + 1, len(rs)):
                        if segs[i].crosses(segs[j]):
                            rs[i]["s"], rs[j]["s"] = rs[j]["s"], rs[i]["s"]
                            segs[i], segs[j] = seg(rs[i]), seg(rs[j])
                            swapped = True
                            self.stats["uncross_swaps"] += 1
                if not swapped:
                    break
        for r in self.rows:
            r["s_shift"] = abs(r["s"] - r["s0"])
            line = self.edges[r["eid"]][3]["geometry"]
            T = line.interpolate(r["s"]).coords[0]
            w = np.subtract(r["facade"], T)
            d = float(np.linalg.norm(w))
            t = _tangent(line, r["s"])
            r["perpendicular"] = bool(d < 1e-6 or abs(np.dot(t, w)) / d < math.sin(math.radians(10)))

    # -- final repair: facade from the tap, global uncrossing, street clearance ------
    def _tap_point(self, r):
        return self.edges[r["eid"]][3]["geometry"].interpolate(r["s"]).coords[0]

    def _set_facade(self, r):
        T = Point(self._tap_point(r))
        r["facade"] = shortest_line(T, r["fp"]).coords[-1]

    def _service(self, r):
        return LineString(r["path"]) if r.get("path") else LineString([self._tap_point(r), r["facade"]])

    def _dogleg(self, r, offset=1.0):
        """Right-angle service: off the main at the tap, along the frontage (offset towards the
        building) to the foot of the building's perpendicular, then into the plot."""
        line = self.edges[r["eid"]][3]["geometry"]
        sl = shortest_line(line, r["fp"])
        P, Q = sl.coords[0], sl.coords[-1]
        sP = line.project(Point(P))
        T = self._tap_point(r)
        t = _tangent(line, sP)
        n = np.array([-t[1], t[0]]) * r["side"]
        a = np.add(T, n * offset)
        b = np.add(P, n * offset)
        if np.linalg.norm(np.subtract(Q, P)) <= offset:
            return [T, tuple(a), tuple(Q)] if np.linalg.norm(np.subtract(a, Q)) > 0 else [T, tuple(Q)]
        return [T, tuple(a), tuple(b), tuple(Q)]

    def _hits_street(self, r):
        g = self._service(r)
        if g.length == 0:
            return False
        T = Point(g.coords[0])
        for k in self.tree.query(g, predicate="intersects"):
            if self.eids[k] == r["eid"]:
                inter = g.intersection(self.geoms[k])
                if not _beyond(inter, T):
                    continue
                return True  # crosses its own edge again (bend)
            inter = g.intersection(self.geoms[k])
            if _beyond(inter, T):
                return True
        return False

    def _repair(self, rounds=8):
        """Facade = nearest footprint point to the final tap; then exchange taps of crossing
        services (any two, each exchange shortens the total length) and of services that
        clip another street, as long as both services stay admissible."""
        from .planar import crossing_pairs
        for r in self.rows:
            self._set_facade(r)
        by_slot = self.by_edge()
        for _ in range(rounds):
            changed = 0
            # services clipping another street: try the nearest slots of the same edge
            for r in [r for r in self.rows if self._hits_street(r)]:
                slots = sorted((o for o in by_slot[r["eid"]] if o is not r), key=lambda o: abs(o["s"] - r["s0"]))
                for o in slots[:12]:
                    r["s"], o["s"] = o["s"], r["s"]
                    self._set_facade(r); self._set_facade(o)
                    if not self._hits_street(r) and not self._hits_street(o):
                        changed += 1
                        break
                    r["s"], o["s"] = o["s"], r["s"]
                    self._set_facade(r); self._set_facade(o)
            # crossing services: exchange the two taps (also across edges of one face)
            rows = [r for r in self.rows if r["dist"] > 0]
            segs = [self._service(r) for r in rows]
            for i, j in crossing_pairs([g for g in segs]):
                a, b = rows[i], rows[j]
                old = (a["eid"], a["s"], a["side"], b["eid"], b["s"], b["side"])
                before = self._service(a).length + self._service(b).length
                a["eid"], a["s"], b["eid"], b["s"] = b["eid"], b["s"], a["eid"], a["s"]
                self._set_facade(a); self._set_facade(b)
                ok = (not self._hits_street(a) and not self._hits_street(b)
                      and not self._service(a).crosses(self._service(b))
                      and self._service(a).length + self._service(b).length < before - 1e-9)
                if ok:
                    for x in (a, b):
                        x["side"] = self._side_of(x)
                    changed += 1
                    self.stats["uncross_swaps"] += 1
                else:
                    a["eid"], a["s"], a["side"], b["eid"], b["s"], b["side"] = old
                    self._set_facade(a); self._set_facade(b)
            by_slot = self.by_edge()
            if not changed:
                break
        for r in self.rows:
            r["path"] = [self._tap_point(r), r["facade"]]
        self._share_services(radius=40.0)
        if self.stats["service_crossings_left"] or self.stats["street_clipping_left"]:
            self._share_services(radius=80.0)
        for r in self.rows:
            r["dist"] = float(self._service(r).length)
            if r.get("host"):
                continue
            r["s_shift"] = abs(r["s"] - r["s0"])
            line = self.edges[r["eid"]][3]["geometry"]
            w = np.subtract(r["facade"], self._tap_point(r))
            d = float(np.linalg.norm(w))
            t = _tangent(line, r["s"])
            r["perpendicular"] = bool(d < 1e-6 or abs(np.dot(t, w)) / d < math.sin(math.radians(10)))
        clip, pairs = self._conflicts()
        self.stats["service_crossings_left"] = len(pairs)
        self.stats["street_clipping_left"] = len(clip)

    def _retap(self):
        """Last resort for a service that still conflicts: tap it on its next admissible street
        edge, at the foot of its perpendicular, if that slot keeps the tap spacing and the new
        service crosses nothing."""
        clip, pairs = self._conflicts()
        todo = {id(r): r for r in clip}
        for a, b in pairs:
            for r in sorted((a, b), key=lambda r: -self._service(r).length):
                if not r.get("host") and not any(o.get("host") == r["bldg_id"] for o in self.rows):
                    todo[id(r)] = r
                    break
        moved = 0
        for r in todo.values():
            others = [o for o in self.rows if o is not r and (o["dist"] > 0 or o.get("host"))]
            segs = [self._service(o) for o in others]
            stree = STRtree(segs)
            for opt in r.get("options", []):
                taken = [o["s"] for o in self.rows if o["eid"] == opt["eid"] and not o.get("host") and o is not r]
                if any(abs(t - opt["s0"]) < self.min_spacing for t in taken):
                    continue
                old = {k: r[k] for k in ("eid", "s", "s0", "side", "facade", "path", "perpendicular", "dist")}
                r.update(eid=opt["eid"], s=opt["s0"], s0=opt["s0"], side=opt["side"], perpendicular=opt["perpendicular"])
                self._set_facade(r)
                r["path"] = [self._tap_point(r), r["facade"]]
                g = self._service(r)
                ok = not self._hits_street(r) and not any(
                    not g.intersection(segs[k]).is_empty for k in stree.query(g, predicate="intersects"))
                if ok:
                    r["dist"] = float(g.length)
                    moved += 1
                    break
                r.update(old)
        clip, pairs = self._conflicts()
        self.stats.update(retapped=moved, service_crossings_left=len(pairs), street_clipping_left=len(clip))

    def _exchange_fixed(self, rounds=500):
        """Last step, for crossing pairs that survive every repair: the exchange of the uncrossing lemma
        (P1) itself, with the facade points kept.  _repair recomputes the facade of each building after
        an exchange (nearest footprint point to its new tap), which can make the exchanged services cross
        again; with the facades kept, the lemma guarantees a shorter pair.  As in P1, an exchange is
        accepted if it shortens the pair and neither service then clips a street; the total length of
        the services strictly decreases, so the loop stops (and rounds bounds it in any case).

        Added 2026-11-23 (P5, Bukavu): removes crossings that _repair rejected in dense blocks.  It only
        acts on pairs left after _retap, so designs with none (Gombe, P1-P4) are unchanged."""
        moved = 0
        for _ in range(rounds):
            clip, pairs = self._conflicts()
            if not pairs:
                break
            progress = False
            for a, b in pairs:
                if a.get("host") or b.get("host") or any(o.get("host") in (a["bldg_id"], b["bldg_id"]) for o in self.rows):
                    continue  # services others hang on keep their taps
                keep = {id(r): {k: r[k] for k in ("eid", "s", "side", "path", "dist", "perpendicular")} for r in (a, b)}
                before = self._service(a).length + self._service(b).length
                a["eid"], a["s"], b["eid"], b["s"] = b["eid"], b["s"], a["eid"], a["s"]
                for r in (a, b):
                    r["path"] = [self._tap_point(r), r["facade"]]
                    r["dist"] = float(self._service(r).length)
                ok = (a["dist"] + b["dist"] < before - 1e-9 and not self._hits_street(a) and not self._hits_street(b))
                if ok:
                    for r in (a, b):
                        r["side"] = self._side_of(r)
                    moved += 1
                    progress = True
                    break  # pairs changed: recompute them
                for r in (a, b):
                    r.update(keep[id(r)])
            if not progress:
                break
        clip, pairs = self._conflicts()
        self.stats.update(exchanged_fixed=moved, service_crossings_left=len(pairs), street_clipping_left=len(clip))

    @staticmethod
    def _hangs_on(row, bldg_id, byid):
        """True if row's service hangs, directly or through other branches, on bldg_id (or is it)."""
        seen = set()
        x = row
        while x is not None and x["bldg_id"] not in seen:
            if x["bldg_id"] == bldg_id:
                return True
            seen.add(x["bldg_id"])
            x = byid.get(x.get("host")) if x.get("host") else None
        return False

    @staticmethod
    def _tee_related(a, b):
        """Services that meet at a tee by construction (host and branch, or two branches of one host)."""
        return (a.get("host") == b["bldg_id"] or b.get("host") == a["bldg_id"]
                or (a.get("host") is not None and a.get("host") == b.get("host")))

    def _conflicts(self):
        """Rows whose service clips a street, and crossing pairs of services (tees excluded)."""
        from .planar import crossing_pairs
        live = [r for r in self.rows if r["dist"] > 0 or r.get("host")]
        segs = [self._service(r) for r in live]
        clip = [r for r in live if not r.get("host") and self._hits_street(r)]
        pairs = [(live[i], live[j]) for i, j in crossing_pairs(segs) if not self._tee_related(live[i], live[j])]
        return clip, pairs

    def _share_services(self, radius=40.0, rounds=12):
        """Shared services, iterated to a fixpoint: a building whose service clips a street or
        crosses another service branches off the nearest service of a neighbour (a tee; it shares
        the neighbour's tap).  A branch is accepted only if it crosses no street and no service
        drawn at that moment, so every re-attachment removes a crossing and adds none: the number
        of crossings strictly decreases and the loop stops."""
        self.stats["shared"] = self.stats.get("shared", 0)
        for _ in range(rounds):
            clip, pairs = self._conflicts()
            if not clip and not pairs:
                break
            dependents = defaultdict(int)
            for r in self.rows:
                if r.get("host"):
                    dependents[r["host"]] += 1
            bad = {id(r): r for r in clip if not dependents[r["bldg_id"]]}
            for a, b in pairs:
                # re-attach the longer service, unless other buildings hang on it
                order = sorted((a, b), key=lambda r: -self._service(r).length)
                pick = next((r for r in order if not dependents[r["bldg_id"]]), None)
                if pick is not None:
                    bad[id(pick)] = pick
            moved = 0
            for r in bad.values():
                live = [o for o in self.rows if (o["dist"] > 0 or o.get("host")) and o is not r]
                segs = [self._service(o) for o in live]
                stree = STRtree(segs)
                best = None
                byid = {o["bldg_id"]: o for o in self.rows}
                for k in stree.query(r["fp"].buffer(radius)):
                    host = live[k]
                    if self._hangs_on(host, r["bldg_id"], byid):
                        continue  # keeps the host relation acyclic (a forest of service trees)
                    br = shortest_line(segs[k], r["fp"])
                    if br.length == 0 or br.length > radius:
                        continue
                    if self._branch_clear(br, host, stree, segs, live) and (best is None or br.length < best[1].length):
                        best = (host, br)
                if best is None:
                    continue
                host, br = best
                r.update(eid=host["eid"], s=host["s"], side=host["side"], host=host["bldg_id"],
                         path=list(br.coords), facade=br.coords[-1], perpendicular=host["perpendicular"],
                         dist=float(br.length))
                moved += 1
            self.stats["shared"] += moved
            if not moved:
                break
        clip, pairs = self._conflicts()
        self.stats["service_crossings_left"] = len(pairs)
        self.stats["street_clipping_left"] = len(clip)

    def _branch_clear(self, br, host, stree, segs, live):
        """A branch must not cross a street or any service other than its host (at the tee)."""
        tee = Point(br.coords[0])
        for k in self.tree.query(br, predicate="intersects"):
            inter = br.intersection(self.geoms[k])
            if _beyond(inter, tee):
                return False
        for k in stree.query(br, predicate="intersects"):
            other = live[k]
            if other is host:
                continue
            inter = br.intersection(segs[k])
            if _beyond(inter, tee):
                return False
        return True

    def _side_of(self, r):
        line = self.edges[r["eid"]][3]["geometry"]
        t = _tangent(line, r["s"])
        w = np.subtract(r["facade"], self._tap_point(r))
        return 1 if t[0] * w[1] - t[1] * w[0] >= 0 else -1

    def by_edge(self, key="s"):
        out = defaultdict(list)
        for r in self.rows:
            out[r["eid"]].append(r)
        for rs in out.values():
            rs.sort(key=lambda r: r[key])
        return out


class CorridorGraph:
    """G0 subdivided at the taps: a simple plane graph with vertex positions and edge arcs."""

    def __init__(self, street, taps, weight_factor=None):
        self.street = street
        self.G = nx.Graph()
        self.tap_vertex = {}
        self.G.graph["tap_vertex"] = self.tap_vertex
        per_edge = taps.by_edge()
        seen_pairs = set()
        for eid, (u0, v0, k, d) in street.edges.items():
            line = d["geometry"]
            u, v = d["start"], d["end"]
            stops = [(0.0, u)]
            for r in per_edge.get(eid, []):
                if r.get("host"):
                    continue  # shared service: uses the tap of its host
                tv = ("t", r["bldg_id"])
                stops.append((r["s"], tv))
                self.tap_vertex[r["bldg_id"]] = tv
            pair = (min(u, v), max(u, v))
            if len(stops) == 1 and (u == v or pair in seen_pairs):
                stops.append((line.length / 2, ("m", eid)))        # keep Gc simple (no loops / multi-edges)
            if u == v and len(stops) == 2:
                stops.append((line.length * 0.75, ("m2", eid)))
            stops.append((line.length, v))
            seen_pairs.add(pair)
            f = weight_factor(d) if weight_factor else 1.0
            for (sa, a), (sb, b) in zip(stops[:-1], stops[1:]):
                L = max(sb - sa, 1e-3)
                self.G.add_edge(a, b, eid=eid, s0=sa, s1=sb, length=L, weight=L * f,
                                road_class=d["road_class"])
            for s, x in stops:
                self.G.nodes[x]["xy"] = line.interpolate(s).coords[0]
                self.G.nodes[x]["street_vertex"] = not isinstance(x, tuple)
        hosted = {r["bldg_id"]: r["host"] for r in taps.rows if r.get("host")}
        for b, h in hosted.items():
            seen = {b}
            while h in hosted:  # a host may itself be hosted (the relation is acyclic)
                if h in seen:
                    raise ValueError(f"cyclic shared services at {b}")
                seen.add(h)
                h = hosted[h]
            self.tap_vertex[b] = self.tap_vertex[h]
        self.hosted = hosted

    def geometry(self, a, b):
        d = self.G.edges[a, b]
        line = self.street.edge_geom(d["eid"])
        g = substring(line, d["s0"], d["s1"])
        # orient from a to b
        if Point(g.coords[0]).distance(Point(self.G.nodes[a]["xy"])) > Point(g.coords[-1]).distance(Point(self.G.nodes[a]["xy"])):
            g = LineString(list(g.coords)[::-1])
        return g
