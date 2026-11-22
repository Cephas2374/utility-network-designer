"""
Graph design tools shared by the three utilities.

  steiner_tree       Mehlhorn 2-approximation of the Steiner tree on Gc
  shortest_path_tree union of shortest paths from a root (radial feeders)
  prune_leaves       remove non-terminal leaves
  split_vertices     vertex splitting to max degree (rule G3), rotation preserving
  orient             arborescence from the sources with downstream sums
  augment_faces      face-based loop augmentation of the trunk (rule G6)
  site_supply_points capacitated graph-Voronoi siting with reach (rule G10)
  place_along_runs   components at arc-length spacing (rule G11)
"""
import math
from collections import defaultdict, deque

import networkx as nx
import numpy as np
from networkx.algorithms.approximation import steiner_tree as nx_steiner

SPLIT_LINK_M = 0.5  # Length of the pipe/cable piece between two tees after vertex splitting


# -- trees ---------------------------------------------------------------------
def steiner_tree(Gc, terminals, weight="weight"):
    """Mehlhorn's 2(1-1/t) approximation of the minimum Steiner tree (Gc must be connected on terminals)."""
    T = nx_steiner(Gc, list(terminals), weight=weight, method="mehlhorn")
    return prune_leaves(nx.Graph(T), set(terminals))


def shortest_path_tree(Gc, root, terminals, weight="weight", allowed=None):
    """Union of shortest root-terminal paths (a tree: shortest paths from one root form an arborescence)."""
    H = Gc.subgraph(allowed) if allowed is not None else Gc
    pred, dist = nx.dijkstra_predecessor_and_distance(H, root, weight=weight)
    T = nx.Graph()
    T.add_node(root)
    for t in terminals:
        x = t
        while x != root and not T.has_edge(x, pred[x][0]):
            p = pred[x][0]
            T.add_edge(p, x, **Gc.edges[p, x])
            x = p
    return T, dist


def prune_leaves(T, keep):
    """Iteratively delete leaves that are not terminals or sources."""
    leaves = deque(v for v in T if T.degree(v) <= 1 and v not in keep)
    while leaves:
        v = leaves.popleft()
        if v not in T or v in keep or T.degree(v) > 1:
            continue
        nb = list(T.neighbors(v))
        T.remove_node(v)
        for u in nb:
            if T.degree(u) <= 1 and u not in keep:
                leaves.append(u)
    return T


def copy_attrs(N, Gc):
    for a, b in N.edges():
        if "length" not in N.edges[a, b] and Gc.has_edge(a, b):
            N.edges[a, b].update(Gc.edges[a, b])
    for v in N:
        N.nodes[v].update(Gc.nodes[v])
    return N


# -- vertex splitting (G3) -------------------------------------------------------
def split_vertices(N, max_degree, exempt=None):
    """Replace every vertex of degree k > D by a path of k-D+1 vertices of degree <= D.

    The incident edges are taken in their cyclic (rotation) order and given to the
    new vertices in consecutive blocks, so the split is drawn inside a small disk
    around the old vertex without crossings: planarity is preserved.
    Returns the number of vertices split and of link edges added.
    """
    exempt = exempt or {}
    n_split = n_links = 0
    for v in [x for x in N if N.degree(x) > exempt.get(x, max_degree)]:
        D = exempt.get(v, max_degree)
        k = N.degree(v)
        x0, y0 = N.nodes[v]["xy"]
        nbrs = sorted(N.neighbors(v), key=lambda w: _angle(N, v, w))
        # first block keeps D-1 edges, middle blocks D-2, last block D-1 (path of tees)
        blocks, i = [], 0
        sizes = []
        rest = k
        while rest > 0:
            cap = D - 1 if not sizes or rest <= D - 1 else D - 2
            take = min(cap, rest)
            sizes.append(take)
            rest -= take
        if len(sizes) == 1:
            continue
        attrs = dict(N.nodes[v])
        edata = {w: dict(N.edges[v, w]) for w in nbrs}
        N.remove_node(v)
        names = [v] + [("x", v, j) for j in range(1, len(sizes))]
        for j, (name, size) in enumerate(zip(names, sizes)):
            block = nbrs[i:i + size]
            i += size
            ang = np.mean([_angle_xy((x0, y0), N.nodes[w]["xy"]) for w in block])
            r = SPLIT_LINK_M * j  # tees placed on a short polyline inside the junction disk
            N.add_node(name, **{**attrs, "xy": (x0 + r * math.cos(ang), y0 + r * math.sin(ang)),
                                "split_of": v})
            for w in block:
                N.add_edge(name, w, **edata[w])
        for a, b in zip(names[:-1], names[1:]):
            N.add_edge(a, b, length=SPLIT_LINK_M, weight=SPLIT_LINK_M, split_link=True, eid=None)
            n_links += 1
        n_split += 1
    return n_split, n_links


def _angle_xy(p, q):
    return math.atan2(q[1] - p[1], q[0] - p[0])


def _angle(N, v, w):
    return _angle_xy(N.nodes[v]["xy"], N.nodes[w]["xy"])


# -- edge subdivision (G8) -------------------------------------------------------
def subdivide(N, Gc, street, max_len):
    """Subdivide every edge longer than max_len into equal pieces (degree-2 vertices).

    Subdivision keeps the drawing, so G1 (plane), G2 and G3 (degree) are unchanged.
    """
    count = 0
    for a, b in [e for e in N.edges() if N.edges[e]["length"] > max_len]:
        d = dict(N.edges[a, b])
        n = math.ceil(d["length"] / max_len)
        s0, s1 = d.get("s0"), d.get("s1")
        line = street.edge_geom(d["eid"]) if d.get("eid") is not None else None
        # orient the arc so that s runs from a to b
        if line is not None:
            pa = np.array(N.nodes[a]["xy"])
            if np.linalg.norm(np.array(line.interpolate(s0).coords[0]) - pa) > \
                    np.linalg.norm(np.array(line.interpolate(s1).coords[0]) - pa):
                s0, s1 = s1, s0
        N.remove_edge(a, b)
        prev = a
        for i in range(1, n + 1):
            x = b if i == n else ("s", a, b, i)
            if i < n:
                if line is not None:
                    xy = line.interpolate(s0 + (s1 - s0) * i / n).coords[0]
                else:
                    pa, pb = np.array(N.nodes[a]["xy"]), np.array(N.nodes[b]["xy"])
                    xy = tuple(pa + (pb - pa) * i / n)
                N.add_node(x, xy=xy, street_vertex=False)
            sa = None if s0 is None else s0 + (s1 - s0) * (i - 1) / n
            sb = None if s0 is None else s0 + (s1 - s0) * i / n
            piece = {**d, "length": d["length"] / n, "weight": d.get("weight", d["length"]) / n,
                     "s0": min(sa, sb) if sa is not None else None, "s1": max(sa, sb) if sa is not None else None}
            N.add_edge(prev, x, **piece)
            prev = x
        count += 1
    return count


# -- orientation and downstream sums -------------------------------------------
def orient(T, roots, load=None):
    """Arborescence of a forest from its roots (BFS); downstream sums of node loads.

    load: dict vertex -> dict of quantities.  Returns (parent, order, down) where
    down[v] is the sum of loads in the subtree of v.  In a tree the flow on the arc
    (parent(v), v) equals down[v] by conservation, and it is the unique solution.
    """
    parent, order = {}, []
    for r in roots:
        if r not in T or r in parent:
            continue
        parent[r] = None
        q = deque([r])
        while q:
            x = q.popleft()
            order.append(x)
            for y in T.neighbors(x):
                if y not in parent:
                    parent[y] = x
                    q.append(y)
    down = {}
    if load is not None:
        keys = {k for d in load.values() for k in d}
        for v in reversed(order):
            acc = {k: load.get(v, {}).get(k, 0.0) for k in keys}
            down[v] = acc
        for v in reversed(order):
            p = parent[v]
            if p is not None:
                for k, val in down[v].items():
                    down[p][k] += val
    return parent, order, down


# -- face-based loop augmentation (G6) ------------------------------------------
def augment_faces(N, Gc, street, trunk_edges, max_added_m):
    """Close loops around trunk bridges by adding street-face boundaries.

    For a bounded face f of G0 let W(f) be its boundary walk without dangling
    edges.  Every edge used once by a closed walk lies on a cycle, so after adding
    W(f) no edge of W(f) is a bridge.  Greedy: repeatedly add the face that covers
    the most remaining trunk bridges per metre of new main, with at most
    max_added_m of new main per face.  Trunk bridges that are bridges of G0 can
    never be covered (they lie on no cycle of G0) and are reported.
    """
    ef = street.edge_faces()
    face_edges = defaultdict(set)
    for eid, (fl, fr) in ef.items():
        if fl == fr:
            continue  # dangling edge: same face on both sides
        for f in (fl, fr):
            if f >= 0:
                face_edges[f].add(eid)
    sub = defaultdict(list)  # street edge -> its Gc edges
    for a, b, d in Gc.edges(data=True):
        if d.get("eid") is not None:
            sub[d["eid"]].append((a, b))

    def added(f):
        return [(a, b) for e in face_edges[f] for a, b in sub[e] if not N.has_edge(a, b)]

    chosen, added_len = [], 0.0
    g0_bridge_eids = _street_bridge_eids(street)
    while True:
        br = {tuple(sorted(e, key=str)) for e in nx.bridges(N)}
        tb = [e for e in trunk_edges if tuple(sorted(e, key=str)) in br and N.has_edge(*e)]
        if not tb:
            break
        tb_eids = defaultdict(int)
        for a, b in tb:
            tb_eids[N.edges[a, b].get("eid")] += 1
        best, best_score = None, 0.0
        cand = {f for e in tb_eids for f in ef.get(e, ()) if f >= 0 and f not in chosen}
        for f in cand:
            new = added(f)
            L = sum(Gc.edges[a, b]["length"] for a, b in new)
            if L > max_added_m:
                continue
            cover = sum(c for e, c in tb_eids.items() if e in face_edges[f])
            score = cover / (L + 1.0)
            if cover and score > best_score:
                best, best_score = f, score
        if best is None:
            break
        for a, b in added(best):
            N.add_edge(a, b, **Gc.edges[a, b], loop_closure=True)
            N.nodes[a].update(Gc.nodes[a])
            N.nodes[b].update(Gc.nodes[b])
            added_len += Gc.edges[a, b]["length"]
        chosen.append(best)
    br = {tuple(sorted(e, key=str)) for e in nx.bridges(N)}
    left = [e for e in trunk_edges if N.has_edge(*e) and tuple(sorted(e, key=str)) in br]
    forced = [e for e in left if N.edges[e].get("eid") in g0_bridge_eids]
    # bridges whose every incident face would need more than max_added_m of new main
    too_long = 0
    for a, b in left:
        e = N.edges[a, b].get("eid")
        if e in g0_bridge_eids:
            continue
        fs_ = [f for f in ef.get(e, ()) if f >= 0]
        if all(sum(Gc.edges[x, y]["length"] for x, y in added(f)) > max_added_m for f in fs_):
            too_long += 1
    return {"faces": chosen, "added_m": added_len, "trunk_bridges_left": len(left),
            "trunk_bridges_in_G0_bridges": len(forced), "trunk_bridges_face_too_long": too_long}


def _street_bridge_eids(street):
    H = nx.Graph()
    count = defaultdict(int)
    for eid, (u, v, k, d) in street.edges.items():
        count[(min(u, v), max(u, v))] += 1
        H.add_edge(u, v)
    out = set()
    br = {(min(a, b), max(a, b)) for a, b in nx.bridges(H)}
    for eid, (u, v, k, d) in street.edges.items():
        p = (min(u, v), max(u, v))
        if p in br and count[p] == 1 and u != v:
            out.add(eid)
    return out


# -- graph-Voronoi siting (G10) ---------------------------------------------------
def multi_source_dijkstra(Gc, sites, weight="length"):
    """Distances, predecessors and owner (nearest site) of every vertex.

    Ties are broken by the order of the sites, so every vertex has one owner and the
    predecessor of a vertex has the same owner: each graph-Voronoi cell induces a
    connected subgraph, and the predecessor arcs form a shortest-path forest.
    """
    import heapq
    dist, pred, owner = {}, {}, {}
    import itertools
    tie = itertools.count()
    heap = [(0.0, i, next(tie), s, None) for i, s in enumerate(sites)]
    heapq.heapify(heap)
    while heap:
        d, i, _, v, p = heapq.heappop(heap)
        if v in dist:
            continue
        dist[v], pred[v], owner[v] = d, p, sites[i]
        for w, e in Gc[v].items():
            if w not in dist:
                heapq.heappush(heap, (d + e[weight], i, next(tie), w, v))
    return dist, pred, owner


def voronoi_refine(Gc, demand, load_fn, capacity, reach, candidate, init=None, relocate=2, max_sites=2000):
    """Supply points whose graph-Voronoi cells meet the capacity and reach rules.

    demand:    terminal -> load record ('w' = weight for the 1-median)
    load_fn:   list of load records -> design load of a cell
    candidate: vertex -> bool, admissible supply-point locations
    A cell that is overloaded or reaches beyond `reach` is split by opening a new site on
    the path to its farthest terminal, at reach/2 from that terminal; then every site is
    moved to the 1-median of its cell when the cell stays valid.  Returns (sites, dist, pred, owner).
    """
    terms = list(demand)
    if not init:
        init = greedy_cover(Gc, demand, load_fn, capacity, reach, candidate)
    sites = list(dict.fromkeys(init))

    def cells_of(owner):
        cells = defaultdict(list)
        for t in terms:
            cells[owner[t]].append(t)
        return cells

    def bad_cells(dist, owner):
        """Cells that break a rule and can still be split (more than one terminal)."""
        out = []
        for s, ts in cells_of(owner).items():
            far = max(ts, key=lambda t: dist[t])
            if (dist[far] > reach or load_fn([demand[t] for t in ts]) > capacity) and len(ts) > 1:
                out.append((s, far))
        return out

    for _round in range(relocate + 1):
        while True:
            dist, pred, owner = multi_source_dijkstra(Gc, sites)
            bad = bad_cells(dist, owner)
            if not bad or len(sites) >= max_sites:
                break
            for s, far in bad:
                sites.append(_walk_back(Gc, pred, far, reach / 2, candidate, sites))
            sites = list(dict.fromkeys(sites))
        if _round == relocate:
            break
        moved = False
        cells = cells_of(owner)
        new = []
        for s in sites:
            ts = cells.get(s)
            if not ts:
                continue
            m = _one_median(Gc, s, ts, demand, reach, candidate, owner)
            moved |= m != s
            new.append(m)
        if not moved:
            break
        trial = list(dict.fromkeys(new))
        d2, p2, o2 = multi_source_dijkstra(Gc, trial)
        if len(bad_cells(d2, o2)) <= len(bad_cells(dist, owner)):
            sites = trial
    dist, pred, owner = multi_source_dijkstra(Gc, sites)
    return sites, dist, pred, owner


def enforce_spacing(Gc, sites, demand, load_fn, capacity, reach, spacing):
    """Remove supply points closer than `spacing` (straight line) to another one, as long as the
    remaining graph-Voronoi cells still meet capacity and reach.  Returns (sites, removed, conflicts)."""
    sites = list(sites)
    removed, conflicts = 0, 0
    tried = set()
    while True:
        xy = np.array([Gc.nodes[s]["xy"] for s in sites])
        d = np.sqrt(((xy[:, None, :] - xy[None, :, :]) ** 2).sum(-1))
        d[np.diag_indices(len(sites))] = np.inf
        pairs = sorted(((d[i, j], i, j) for i in range(len(sites)) for j in range(i + 1, len(sites))
                        if d[i, j] < spacing and (sites[i], sites[j]) not in tried))
        if not pairs:
            break
        _, i, j = pairs[0]
        tried.add((sites[i], sites[j]))
        ok = False
        for k in (j, i):
            trial = sites[:k] + sites[k + 1:]
            dist, pred, owner = multi_source_dijkstra(Gc, trial)
            cells = defaultdict(list)
            for t in demand:
                cells[owner[t]].append(t)
            if all(max(dist[t] for t in ts) <= reach and load_fn([demand[t] for t in ts]) <= capacity
                   for ts in cells.values()):
                sites, ok = trial, True
                removed += 1
                break
        if not ok:
            ok = _relocate_apart(Gc, sites, (i, j), demand, load_fn, capacity, reach, spacing)
        conflicts += not ok
    return sites, removed, conflicts


def _relocate_apart(Gc, sites, pair, demand, load_fn, capacity, reach, spacing, tries=40):
    """Move one site of a too-close pair inside its own cell to a vertex at least `spacing` from
    every other site, if the cells then still meet capacity and reach (sites edited in place)."""
    dist, pred, owner = multi_source_dijkstra(Gc, sites)
    for k in pair[::-1]:
        s = sites[k]
        others = np.array([Gc.nodes[x]["xy"] for x in sites if x != s])
        cell = [v for v, o in owner.items() if o == s and not isinstance(v, tuple)]
        cell.sort(key=lambda v: dist[v])
        for v in cell[:tries]:
            if np.sqrt(((others - np.array(Gc.nodes[v]["xy"])) ** 2).sum(1)).min() < spacing:
                continue
            trial = sites[:k] + [v] + sites[k + 1:]
            d2, p2, o2 = multi_source_dijkstra(Gc, trial)
            cells = defaultdict(list)
            for t in demand:
                cells[o2[t]].append(t)
            if all(max(d2[t] for t in ts) <= reach and load_fn([demand[t] for t in ts]) <= capacity
                   for ts in cells.values()):
                sites[k] = v
                return True
    return False


def greedy_cover(Gc, demand, load_fn, capacity, reach, candidate):
    """Lazy greedy set cover of the terminals by network balls B(c, reach), each ball taking
    its nearest uncovered terminals up to the capacity (ln n approximation of the minimum
    number of supply points).  Every terminal ends within `reach` of a chosen site, hence
    within `reach` of its nearest site: the graph-Voronoi assignment keeps the reach."""
    import heapq
    balls = {}
    for c in (v for v in Gc if candidate(v)):
        d = nx.single_source_dijkstra_path_length(Gc, c, cutoff=reach, weight="length")
        ts = sorted((t for t in d if t in demand), key=d.get)
        if ts:
            balls[c] = ts
    uncovered = set(demand)

    def gain(c):
        take, recs = [], []
        for t in balls[c]:
            if t not in uncovered:
                continue
            if recs and load_fn(recs + [demand[t]]) > capacity:
                break
            take.append(t)
            recs.append(demand[t])
        return take

    heap = [(-len(ts), i, c) for i, (c, ts) in enumerate(balls.items())]
    heapq.heapify(heap)
    sites = []
    while uncovered and heap:
        neg, i, c = heapq.heappop(heap)
        take = gain(c)
        if not take:
            continue
        if heap and len(take) < -heap[0][0]:
            heapq.heappush(heap, (-len(take), i, c))  # stale score: re-insert
            continue
        sites.append(c)
        uncovered -= set(take)
    return sites


def _farthest_points(Gc, terms, k, candidate):
    heavy = terms[0]
    pts = [_nearest_candidate(Gc, heavy, candidate)]
    d = nx.single_source_dijkstra_path_length(Gc, pts[0], weight="length")
    dmin = {t: d.get(t, math.inf) for t in terms}
    while len(pts) < k:
        far = max(terms, key=lambda t: dmin[t])
        c = _nearest_candidate(Gc, far, candidate)
        pts.append(c)
        d = nx.single_source_dijkstra_path_length(Gc, c, weight="length")
        dmin = {t: min(dmin[t], d.get(t, math.inf)) for t in terms}
    return pts


def _nearest_candidate(Gc, v, candidate):
    for x, _ in nx.single_source_dijkstra_path_length(Gc, v, weight="length").items():
        if candidate(x):
            return x
    return v


def _walk_back(Gc, pred, far, back, candidate, sites):
    """Vertex on the shortest path from the cell's site to `far`, about `back` metres from far."""
    x, run = far, 0.0
    best = None
    while pred.get(x) is not None:
        p = pred[x]
        run += Gc.edges[p, x]["length"]
        x = p
        if candidate(x) and x not in sites:
            best = x
            if run >= back:
                break
    if best is None:
        best = _nearest_candidate(Gc, far, lambda v: candidate(v) and v not in sites)
    return best


def _one_median(Gc, s, ts, demand, reach, candidate, owner, pool_size=15):
    """Admissible vertex of the cell minimising the load-weighted distance sum, within reach of all."""
    w = np.array([demand[t].get("w", 1.0) for t in ts])
    xy = np.array([Gc.nodes[t]["xy"] for t in ts])
    c = (xy * w[:, None]).sum(0) / w.sum()
    pool = [v for v, o in owner.items() if o == s and candidate(v)]
    pool.sort(key=lambda v: (Gc.nodes[v]["xy"][0] - c[0]) ** 2 + (Gc.nodes[v]["xy"][1] - c[1]) ** 2)
    pool = pool[:pool_size] + [s]
    best, best_cost = s, math.inf
    for v in pool:
        d = nx.single_source_dijkstra_path_length(Gc, v, cutoff=reach, weight="length")
        if any(t not in d for t in ts):
            continue
        cost = float(sum(wi * d[t] for wi, t in zip(w, ts)))
        if cost < best_cost:
            best, best_cost = v, cost
    return best


def forest_from_pred(Gc, pred, keep):
    """Shortest-path forest restricted to the paths that reach the vertices in keep."""
    F = nx.Graph()
    for t in keep:
        x = t
        F.add_node(x, **Gc.nodes[x])
        while pred.get(x) is not None and not F.has_edge(pred[x], x):
            p = pred[x]
            F.add_node(p, **Gc.nodes[p])
            F.add_edge(p, x, **Gc.edges[p, x])
            x = p
    return F


# -- arc-length placement (G11) -------------------------------------------------
def runs(N):
    """Maximal paths whose interior vertices have degree 2 (the 'runs' of the mains)."""
    ends = {v for v in N if N.degree(v) != 2}
    seen = set()
    out = []
    for a in ends:
        for b in N.neighbors(a):
            if (a, b) in seen:
                continue
            path = [a, b]
            seen.add((a, b)); seen.add((b, a))
            while path[-1] not in ends:
                x = path[-1]
                nxt = [y for y in N.neighbors(x) if y != path[-2]]
                if not nxt:
                    break
                y = nxt[0]
                seen.add((x, y)); seen.add((y, x))
                path.append(y)
                if y == a:
                    break
            out.append(path)
    # pure cycles without degree != 2 vertices
    for comp in nx.connected_components(N):
        if all(N.degree(v) == 2 for v in comp):
            c = list(nx.cycle_basis(N.subgraph(comp))[0])
            out.append(c + [c[0]])
    return out


def place_along_runs(N, spacing, length="length", select=None):
    """Points every `spacing` along each run, at (i + 1/2) L / n with n = ceil(L / spacing).

    Every point of a run is then within spacing/2 of a component along the run.
    Returns a list of (run index, arc position, edge (a, b), offset on edge).
    """
    pts = []
    for ri, path in enumerate(runs(N)):
        edges = list(zip(path[:-1], path[1:]))
        if select is not None and not all(select(N.edges[a, b]) for a, b in edges):
            continue
        lens = [N.edges[a, b][length] for a, b in edges]
        L = sum(lens)
        if L <= 0:
            continue
        n = max(1, math.ceil(L / spacing))
        cum = np.concatenate([[0], np.cumsum(lens)])
        for i in range(n):
            s = (i + 0.5) * L / n
            j = min(int(np.searchsorted(cum, s, side="right") - 1), len(edges) - 1)
            pts.append((ri, s, edges[j], s - cum[j]))
    return pts
