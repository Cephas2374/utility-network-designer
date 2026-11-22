"""
Plane street graph - noding, rotation system, faces, dual graph, Euler check,
planar indices and the geometric plane-embedding test.

G0 = (V0, E0) is built from the noded street edges: every crossing of two
streets is a vertex, so the drawing of G0 is a plane graph (bridges and
tunnels excepted).  For a plane graph with c components the Euler relation is
    |V| - |E| + F_b = c              (F_b bounded faces)
and the cyclomatic number mu = |E| - |V| + c equals F_b: every independent
cycle of a plane graph is a bounded face.
"""
import math
from collections import defaultdict

import networkx as nx
import numpy as np
from shapely import STRtree, polygonize
from shapely.geometry import LineString, Point
from shapely.ops import substring

ROUND = 2  # Vertex snapping precision (cm): noded street data share exact end points


def node_key(xy):
    return (round(xy[0], ROUND), round(xy[1], ROUND))


def bearing(p, q):
    """Bearing of p->q in degrees, clockwise from grid north, in [0, 360)."""
    return math.degrees(math.atan2(q[0] - p[0], q[1] - p[1])) % 360.0


def leaving_angle(line, from_start, probe=5.0):
    """Mathematical angle (rad, CCW from east) of a line leaving its start or end point."""
    L = line.length
    d = min(probe, L / 2) if L > 0 else 0
    if from_start:
        p, q = line.coords[0], line.interpolate(d).coords[0]
    else:
        p, q = line.coords[-1], line.interpolate(L - d).coords[0]
    return math.atan2(q[1] - p[1], q[0] - p[0])


class StreetGraph:
    """Plane street graph with geometry, rotation system, faces and dual."""

    def __init__(self, edges_gdf):
        self.G = nx.MultiGraph()
        self.pos = {}
        ids = {}
        for i, row in enumerate(edges_gdf.itertuples()):
            line = row.geometry
            if line is None or line.is_empty or line.length == 0:
                continue
            a, b = node_key(line.coords[0]), node_key(line.coords[-1])
            for k in (a, b):
                if k not in ids:
                    ids[k] = len(ids)
                    self.pos[ids[k]] = k
            u, v = ids[a], ids[b]
            self.G.add_edge(u, v, key=i, eid=i, start=u, end=v, geometry=line, length=line.length,
                            road_class=getattr(row, "road_class", "residential") or "residential",
                            stroke_id=getattr(row, "stroke_id", -1), name=getattr(row, "name", None))
        self.edges = {d["eid"]: (u, v, k, d) for u, v, k, d in self.G.edges(keys=True, data=True)}
        self._rotation = None
        self._faces = None

    # -- basic structure ------------------------------------------------------
    @property
    def n(self):
        return self.G.number_of_nodes()

    @property
    def m(self):
        return self.G.number_of_edges()

    def components(self):
        return nx.number_connected_components(self.G)

    def edge_geom(self, eid):
        return self.edges[eid][3]["geometry"]

    # -- rotation system ------------------------------------------------------
    def rotation(self):
        """Cyclic CCW order of half-edges (eid, leaves_from_start) around each vertex."""
        if self._rotation is None:
            rot = defaultdict(list)
            for eid, (u, v, k, d) in self.edges.items():
                g = d["geometry"]
                rot[d["start"]].append((leaving_angle(g, True), eid, True))
                rot[d["end"]].append((leaving_angle(g, False), eid, False))
            self._rotation = {x: sorted(h) for x, h in rot.items()}
        return self._rotation

    # -- faces and dual --------------------------------------------------------
    def faces(self):
        """Bounded faces of the drawing (shapely polygonize of the noded edges)."""
        if self._faces is None:
            geoms = [d["geometry"] for _, _, _, d in self.edges.values()]
            self._faces = list(polygonize(geoms).geoms)
        return self._faces

    def edge_faces(self, eps=0.05):
        """For every edge: (face left, face right) w.r.t. its stored direction; -1 = outer face."""
        faces = self.faces()
        tree = STRtree(faces)
        out = {}
        for eid, (u, v, k, d) in self.edges.items():
            g = d["geometry"]
            s = g.length / 2
            p = np.array(g.interpolate(s).coords[0])
            q = np.array(g.interpolate(min(s + 0.5, g.length)).coords[0])
            if np.allclose(p, q):
                q = np.array(g.interpolate(max(s - 0.5, 0)).coords[0])
                t = p - q
            else:
                t = q - p
            t = t / (np.linalg.norm(t) or 1)
            nl = np.array([-t[1], t[0]])
            sides = []
            for sgn in (+1, -1):
                pt = Point(*(p + sgn * eps * nl))
                hit = [int(j) for j in tree.query(pt, predicate="within")]
                sides.append(hit[0] if hit else -1)
            out[eid] = tuple(sides)
        return out

    def dual(self):
        """Planar dual: faces as vertices (outer face = -1), one dual edge per primal edge."""
        D = nx.MultiGraph()
        D.add_nodes_from(range(-1, len(self.faces())))
        for eid, (fl, fr) in self.edge_faces().items():
            D.add_edge(fl, fr, eid=eid)
        return D

    # -- indices ---------------------------------------------------------------
    def metrics(self):
        V, E, C = self.n, self.m, self.components()
        Fb = len(self.faces())
        mu = E - V + C
        deg = dict(self.G.degree())
        bridges = list(nx.bridges(nx.Graph(self.G)))
        # parallel edges are never bridges; nx.Graph collapses them, so keep only true bridges
        bridge_set = {(min(a, b), max(a, b)) for a, b in bridges
                      if self.G.number_of_edges(a, b) == 1}
        bridge_len = sum(d["length"] for u, v, k, d in self.G.edges(keys=True, data=True)
                         if (min(u, v), max(u, v)) in bridge_set)
        total = sum(d["length"] for _, _, _, d in self.edges.values())
        H, phi, hist = self.orientation_entropy()
        ra = self.right_angle_share()
        return {
            "V": V, "E": E, "C": C, "F_bounded": Fb, "euler_residual": V - E + Fb - C,
            "mu": mu, "alpha_meshedness": (E - V + C) / (2 * V - 5) if V > 2 else 0.0,
            "beta": E / V, "gamma": E / (3 * (V - 2)) if V > 2 else 0.0,
            "mean_degree": 2 * E / V,
            "degree_hist": {k: sum(1 for x in deg.values() if x == k) for k in sorted(set(deg.values()))},
            "n_bridges": len(bridge_set), "bridge_length_m": bridge_len, "total_length_m": total,
            "orientation_entropy": H, "orientation_order_phi": phi, "orientation_hist": hist,
            "right_angle_share": ra,
            "planarity_test": nx.check_planarity(nx.Graph(self.G))[0],
        }

    def orientation_entropy(self, bins=36):
        """Boeing (2019) street orientation entropy on length-weighted segment bearings."""
        w = np.zeros(bins)
        width = 360.0 / bins
        for _, _, _, d in self.edges.values():
            c = list(d["geometry"].coords)
            for p, q in zip(c[:-1], c[1:]):
                L = math.dist(p, q)
                if L == 0:
                    continue
                for b in (bearing(p, q), (bearing(p, q) + 180) % 360):
                    w[int(((b + width / 2) % 360) // width)] += L
        p = w / w.sum()
        H = float(-(p[p > 0] * np.log(p[p > 0])).sum())
        Hmax, Hg = math.log(bins), math.log(4)
        phi = 1 - ((H - Hg) / (Hmax - Hg)) ** 2
        return H, phi, p.tolist()

    def right_angle_share(self, tol=10.0):
        """Share of consecutive junction angles (deg >= 3) within 90 +- tol degrees."""
        good = tot = 0
        for x, hs in self.rotation().items():
            if len(hs) < 3:
                continue
            ang = [h[0] for h in hs]
            for i in range(len(ang)):
                a = math.degrees((ang[(i + 1) % len(ang)] - ang[i]) % (2 * math.pi))
                tot += 1
                good += abs(a - 90) <= tol
        return good / tot if tot else 0.0


# -- geometric plane-embedding test (rule G1) -----------------------------------
def crossing_pairs(lines, tol=1e-6):
    """Pairs (i, j) of lines whose drawings meet outside common end points.

    A drawing is plane iff this list is empty.
    """
    tree = STRtree(lines)
    left, right = tree.query(lines, predicate="intersects")
    bad = []
    for i, j in zip(left, right):
        if i >= j:
            continue
        a, b = lines[i], lines[j]
        inter = a.intersection(b)
        if inter.is_empty:
            continue
        ends = [Point(a.coords[0]), Point(a.coords[-1])]
        shared = [e for e in ends if e.distance(Point(b.coords[0])) < tol or e.distance(Point(b.coords[-1])) < tol]
        rest = inter
        for e in shared:
            rest = rest.difference(e.buffer(tol * 10))
        if not rest.is_empty:
            bad.append((int(i), int(j)))
    return bad


def split_line(line, distances):
    """Split a LineString at the given arc-length distances (sorted, strictly inside)."""
    cuts = [0.0] + sorted(distances) + [line.length]
    return [substring(line, a, b) for a, b in zip(cuts[:-1], cuts[1:])]
