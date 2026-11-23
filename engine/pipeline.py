"""
The design pipeline of the engine, independent of where the data come from (paper P5).

    res = design(Inputs(streets=..., buildings=..., boundary=...), log=print)
    write_outputs(res, folder)

`run_gombe.py` calls it with the Gombe data; the QGIS plugin calls it with layers of the project.
Every step is the one of paper P1: demand, plane street graph, taps and corridor graph, water,
electricity, gas, corridor side assignment, method comparison, validation.
"""
import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point

from . import crossings, demand, design as dsg, electricity, gas, store, validate, water
from .config import load_rules
from .corridor import CorridorGraph, Taps
from .planar import StreetGraph

MAJOR = {"primary", "secondary", "tertiary"}

# optional input columns and their defaults (the engine reads them when present)
BUILDING_COLUMNS = {"building": None, "building_levels": None, "quartier": None, "stroke_id": -1,
                    "full_address": None, "street": None, "house_number": None}
STREET_COLUMNS = {"road_class": "residential", "stroke_id": -1, "name": None}
POI_COLUMNS = {"amenity": None, "office": None, "shop": None, "tourism": None, "healthcare": None}
LANDUSE_COLUMNS = {"fclass": None}


@dataclass
class Inputs:
    """The data of one district, in one projected metric CRS.

    streets    noded street centrelines (LineString), optional road_class, stroke_id, name
    buildings  footprints (Polygon) with a unique bldg_id; optional building_levels, building (OSM tag),
               stroke_id (access street of the address), address fields
    boundary   the district polygon; bulk inflows and the gas gate lie on its exterior
    """
    streets: gpd.GeoDataFrame
    buildings: gpd.GeoDataFrame
    boundary: object
    water: Optional[gpd.GeoDataFrame] = None  # water bodies; the largest is the river (no inflows on its bank)
    railways: Optional[gpd.GeoDataFrame] = None
    pois: Optional[gpd.GeoDataFrame] = None
    landuse: Optional[gpd.GeoDataFrame] = None
    elevation: Optional[Callable] = None  # (x, y) -> ground level (m); flat if None
    target_population: Optional[float] = None  # calibrates the floor area per dwelling
    mv_infeed_near: Optional[object] = None  # geometry; the MV infeed is the candidate nearest to it
    gas_gate_near: Optional[object] = None  # geometry; default: the candidate on a primary road farthest east
    kerb_order: Optional[list] = None  # default: the rule pack
    crs: object = field(default=None)


OSM_ROAD_CLASS = {"motorway": "primary", "trunk": "primary", "primary": "primary", "secondary": "secondary",
                  "tertiary": "tertiary", "residential": "residential", "unclassified": "residential",
                  "living_street": "residential", "service": "service", "track": "service"}


def node_streets(streets):
    """Split the street lines at every intersection, so that every crossing is a vertex (the plane
    street graph of P1, Definition 1); each piece keeps the attributes of the line it lies on.  A
    road_class is derived from an OSM 'highway' tag when the layer has none."""
    from shapely.ops import unary_union
    lines = streets.explode(index_parts=False).reset_index(drop=True)
    lines = lines[lines.geometry.length > 0]
    merged = unary_union(list(lines.geometry))
    parts = [g for g in getattr(merged, "geoms", [merged]) if g.length > 0]
    pieces = gpd.GeoDataFrame(geometry=parts, crs=streets.crs)
    mid = gpd.GeoDataFrame(geometry=pieces.geometry.interpolate(0.5, normalized=True), crs=streets.crs)
    j = gpd.sjoin_nearest(mid, lines, how="left")
    j = j[~j.index.duplicated()]
    for c in lines.columns:
        if c != "geometry":
            pieces[c] = j[c].values
    if "road_class" not in pieces.columns and "highway" in pieces.columns:
        pieces["road_class"] = pieces["highway"].map(OSM_ROAD_CLASS).fillna("residential")
    return pieces


def _with_columns(gdf, cols, crs):
    gdf = gdf.copy() if gdf is not None else gpd.GeoDataFrame(geometry=[], crs=crs)
    for c, v in cols.items():
        if c not in gdf.columns:
            gdf[c] = v
    return gdf


def design(inp: Inputs, packs=None, cr=None, log=print):
    """Design and validate the three networks; returns everything the outputs and analyses need."""
    t0 = time.time()
    if packs is None or cr is None:
        p, c = load_rules()
        packs, cr = packs or p, cr or c
    if inp.kerb_order:
        cr = dict(cr, corridor_layout=dict(cr["corridor_layout"], kerb_order=list(inp.kerb_order)))
    crs = inp.crs or inp.streets.crs
    log(f"kerb order {cr['corridor_layout']['kerb_order']}")

    # -- inputs and demand -------------------------------------------------------
    log("demand")
    bldg = _with_columns(inp.buildings, BUILDING_COLUMNS, crs)
    b = demand.estimate(bldg, _with_columns(inp.pois, POI_COLUMNS, crs), _with_columns(inp.landuse, LANDUSE_COLUMNS, crs),
                        packs, inp.target_population)
    bi = b.set_index("bldg_id")
    boundary = inp.boundary
    water_poly = inp.water if inp.water is not None and len(inp.water) else None
    river = water_poly.loc[water_poly.area.idxmax()].geometry if water_poly is not None else None
    rail = list(inp.railways.geometry) if inp.railways is not None else []
    streams = [g for g in water_poly.geometry if g is not river] if water_poly is not None else []

    # -- plane street graph -------------------------------------------------------
    log("plane street graph")
    S = StreetGraph(_with_columns(inp.streets, STREET_COLUMNS, crs))
    street_metrics = S.metrics()
    log(f"  V={street_metrics['V']} E={street_metrics['E']} C={street_metrics['C']} "
        f"F_b={street_metrics['F_bounded']} Euler residual={street_metrics['euler_residual']}")

    # -- taps and corridor graph ---------------------------------------------------
    log("taps and corridor graph")
    T = Taps(b, S, min_spacing=packs["water"]["graph"]["min_tap_spacing_m"])
    log(f"  taps: {T.stats}")
    rc = cr["road_classes"]
    C = CorridorGraph(S, T, weight_factor=lambda d: rc.get(d["road_class"], rc["residential"])["cost_factor"])
    G = C.G.subgraph(max(nx.connected_components(C.G), key=len)).copy()
    G.graph["corridor"] = C
    served = [r["bldg_id"] for r in T.rows if C.tap_vertex[r["bldg_id"]] in G]
    rows = [r for r in T.rows if r["bldg_id"] in set(served)]
    elev = inp.elevation or (lambda xy: 0.0)
    major = lambda v: any(G.edges[v, u]["road_class"] in MAJOR for u in G.neighbors(v))
    cands = water.inflow_candidates(G, boundary, river, major)

    # -- water ---------------------------------------------------------------------
    log("water")
    Nw, wi = water.design_water_auto(G, S, served, bi, packs, cands, elev, log)

    # -- electricity -----------------------------------------------------------------
    log("electricity")
    near = inp.mv_infeed_near if inp.mv_infeed_near is not None else boundary.centroid
    infeed = min(cands, key=lambda v: near.distance(Point(G.nodes[v]["xy"])))
    svc_len = {C.tap_vertex[r["bldg_id"]]: r["dist"] for r in rows if not r.get("host")}
    Fe, er = electricity.design_electricity(G, served, bi, packs, infeed, svc_len, log)

    # -- gas -------------------------------------------------------------------------
    log("gas")
    if inp.gas_gate_near is not None:
        gate = min(cands, key=lambda v: inp.gas_gate_near.distance(Point(G.nodes[v]["xy"])))
    else:
        prim = [v for v in cands if any(G.edges[v, u]["road_class"] == "primary" for u in G.neighbors(v))] or cands
        gate = max(prim, key=lambda v: G.nodes[v]["xy"][0])
    hazard = list(b[b["hazard"]].geometry)
    LP, MP, gr_ = gas.design_gas(G, served, bi, packs, gate, hazard, log)

    # -- corridor layout and crossings ---------------------------------------------
    log("corridor side assignment")
    nets = {"gas": [LP, MP], "electricity": [Fe], "water": [Nw]}
    cor = crossings.Corridor(S, G, nets, rows, packs, cr)
    hist = cor.optimise(log=log)
    cor.joint_repair(log=log)
    hist.append(cor.total())
    cross = cor.report(rail=rail, water=streams)

    # -- optimisation comparison (MST / Steiner / SPT) on the same terminals ---------
    log("method comparison")
    comp = compare_methods(G, served, wi["sources"])

    # -- validation ------------------------------------------------------------------
    log("validation")
    val = []
    srcs = {"water": set(wi["sources"]), "electricity": set(er["sites"]), "gas": set(gr_["sites"]) | {gate}}
    val += validate.check_network("water", [Nw], G, packs, rows, srcs["water"], cor, extra=water_rows(wi, packs, Nw))
    # every building on the tap of an MV customer is fed from its own substation on the MV ring
    mvc = {bb for t in er["transformers"].values() if t.get("mv_customer") for bb in t["mv_customer_bldgs"]}
    val += validate.check_network("electricity", [Fe], G, packs, rows, srcs["electricity"], cor,
                                  extra=elec_rows(er, packs, G), exclude_terminals=mvc)
    val += validate.check_network("gas", [LP, MP], G, packs, rows, srcs["gas"], cor, extra=gas_rows(gr_, packs, LP, MP))
    val += crossing_rows(cross, cr)
    val += validate.check_services_buildings(rows, b)  # G13 (added 2026-11-23)
    vdf = pd.DataFrame(val)
    log(f"  {vdf['ok'].sum()} / {len(vdf)} checks pass")
    for r in vdf[~vdf["ok"]].itertuples():
        log(f"  FAIL {r.rule} {r.utility}: {r.quantity} = {r.value} (limit {r.limit}) {r.note}")

    metrics = dict(street=street_metrics, demand=demand_summary(b),
                   floor_area_per_dwelling_m2=b.attrs["floor_area_per_dwelling_m2"],
                   taps=dict(n=len(rows), perpendicular=float(np.mean([r["perpendicular"] for r in rows])),
                             shifted_over_1m=int(sum(r["s_shift"] > 1 for r in rows)), **T.stats),
                   corridor_graph=dict(V=G.number_of_nodes(), E=G.number_of_edges()),
                   water=_clean(wi), electricity=_clean(er["info"]), gas=_clean(gr_["info"]),
                   crossings=_clean(cross), side_assignment_cost=hist, methods=comp,
                   validation=dict(passed=int(vdf["ok"].sum()), total=len(vdf)),
                   runtime_s=time.time() - t0)
    model = dict(G=G, S=S, b=b, rows=rows, Nw=Nw, wi=wi, Fe=Fe, er=er, LP=LP, MP=MP, gr=gr_, gate=gate,
                 infeed=infeed, sigma=cor.sigma, dual=cor.dual, lanes=cor.lane, packs=packs, cr=cr,
                 # what the analyses of P2-P4 need besides the networks, so that they run on any district
                 cands=list(cands), z={v: float(elev(G.nodes[v]["xy"])) for v in G})
    return dict(model=model, metrics=metrics, validation=vdf, C=C, cor=cor, crs=crs)


def write_outputs(res, folder):
    """networks.gpkg, metrics.json, validation.csv, customer_register.csv and model.pkl in folder."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    m, C, cor, crs = res["model"], res["C"], res["cor"], res["crs"]
    res["validation"].to_csv(folder / "validation.csv", index=False)
    write_layers(folder / "networks.gpkg", crs, m["G"], C, m["S"], m["b"], m["rows"], m["Nw"], m["wi"], m["Fe"],
                 m["er"], m["LP"], m["MP"], m["gr"], m["gate"], m["infeed"], cor)
    register(m["rows"], m["b"].set_index("bldg_id"), C, m["er"], m["gr"], m["wi"], cor).to_csv(
        folder / "customer_register.csv", index=False)
    (folder / "metrics.json").write_text(json.dumps(res["metrics"], indent=1, default=str))
    store.save(m, folder / "model.pkl")


def model_context(m):
    """Inflow candidates and a ground-level sampler of a stored design.  Designs written by this
    pipeline carry them; the Gombe design of papers P1-P4 predates that and takes them from data/gombe."""
    G = m["G"]
    if "cands" in m and "z" in m:
        from scipy.spatial import cKDTree
        verts = list(m["z"])
        tree = cKDTree([G.nodes[v]["xy"] for v in verts])
        zs = np.array([m["z"][v] for v in verts])
        return m["cands"], lambda xy: float(zs[tree.query(xy)[1]])
    from .config import DATA
    from .dem import SRTM
    dem = SRTM(DATA / "S05E015.hgt")
    boundary = gpd.read_file(DATA / "gombe_inputs.gpkg", layer="commune_boundary").geometry.iloc[0]
    wp = gpd.read_file(DATA / "gombe_inputs.gpkg", layer="water")
    river = wp.loc[wp.area.idxmax()].geometry
    major = lambda v: any(G.edges[v, u]["road_class"] in MAJOR for u in G.neighbors(v))
    return water.inflow_candidates(G, boundary, river, major), lambda xy: float(dem.sample([xy])[0])


def compare_methods(G, served, sources):
    """Total main length of MST (pruned), Steiner (Mehlhorn) and shortest-path tree on one terminal set."""
    terms = {G.graph["tap_vertex"][b] for b in served} | set(sources)
    L = lambda T: sum(G.edges[e]["length"] for e in T.edges())
    W = lambda T: sum(G.edges[e]["weight"] for e in T.edges())
    st = dsg.steiner_tree(G, terms)
    mst = dsg.prune_leaves(nx.Graph(nx.minimum_spanning_tree(G, weight="weight")), terms)
    H = G.copy()
    H.add_node("super")
    for s in sources:
        H.add_edge("super", s, weight=0.0, length=0.0)
    spt, _ = dsg.shortest_path_tree(H, "super", terms)
    spt.remove_node("super")
    lw = lambda T: dict(length_m=sum(H.edges[e]["length"] for e in T.edges()),
                        cost_m=sum(H.edges[e]["weight"] for e in T.edges()), edges=T.number_of_edges())
    return dict(steiner=dict(length_m=L(st), cost_m=W(st), edges=st.number_of_edges()),
                mst=dict(length_m=L(mst), cost_m=W(mst), edges=mst.number_of_edges()), spt=lw(spt))


def water_rows(wi, packs, N):
    cap = packs["water"]["capacity"]
    R = validate._row
    aug = wi["augmentation"]
    allowed = aug["trunk_bridges_in_G0_bridges"] + aug["trunk_bridges_face_too_long"]
    rows = [R("G6", "water", "trunk bridges left", aug["trunk_bridges_left"], allowed,
              aug["trunk_bridges_left"] <= allowed,
              f"{aug['trunk_bridges_in_G0_bridges']} are bridges of G0 (no cycle exists), "
              f"{aug['trunk_bridges_face_too_long']} need a face of more than the loop-chord limit"),
            R("H", "water", "min pressure at taps (m)", round(wi["min_pressure_mwc"], 2), cap["min_pressure_mwc"],
              wi["min_pressure_mwc"] >= cap["min_pressure_mwc"]),
            R("H", "water", "max velocity (m/s)", round(wi["max_velocity_ms"], 2), cap["max_velocity_ms"],
              wi["max_velocity_ms"] <= cap["max_velocity_ms"]),
            R("H", "water", "max gradient (m/km)", round(wi["max_gradient_m_per_km"], 2),
              cap["max_pressure_drop_per_km_m"], wi["max_gradient_m_per_km"] <= cap["max_pressure_drop_per_km_m"] + 1e-3),
            R("H", "water", "max pressure (m)", round(wi["max_pressure_mwc"], 1), cap["max_pressure_mwc"],
              wi["max_pressure_mwc"] <= cap["max_pressure_mwc"]),
            R("H", "water", "fire residual at the weakest hydrant (m)", round(wi["fire"]["residual_mwc"], 2),
              cap["min_pressure_mwc"], wi["fire"]["residual_mwc"] >= cap["min_pressure_mwc"])]
    comps = wi["components"]
    rows.append(validate.check_spacing("water", N, packs["water"]["components"]["valves"]["required_spacing_m"],
                                       comps["valves"], what="valve"))
    hyd = packs["water"]["components"]["hydrants"]["spacing_m"]
    H = N.edge_subgraph([(a, b) for a, b, d in N.edges(data=True)
                         if d.get("diameter_mm", 0) >= packs["water"]["capacity"]["fire_flow_min_diameter_mm"]]).copy()
    rows.append(validate.check_spacing("water", H, hyd, dsg.place_along_runs(H, hyd), what="hydrant"))
    return rows


def elec_rows(er, packs, G):
    e = packs["electricity"]
    i = er["info"]
    R = validate._row
    on_ring = {x for r in er["rings"] for x in r["transformers"]}
    mv = [s for s, t in er["transformers"].items() if t.get("mv_customer")]
    return [R("G4", "electricity", "MV customer substations not on an MV ring", sum(s not in on_ring for s in mv), 0,
              all(s in on_ring for s in mv), f"{len(mv)} MV customers"),
            R("G4", "electricity", "transformers not on an MV ring",
              sum(s not in on_ring for s in er["sites"]), 0, all(s in on_ring for s in er["sites"])),
            R("G10", "electricity", "max network reach (m)", round(i["max_reach_m"], 1),
              e["topology"]["transformer_max_reach_m"], i["max_reach_m"] <= e["topology"]["transformer_max_reach_m"]),
            R("H", "electricity", "max voltage drop (%)", round(i["max_drop_pct"], 3),
              e["capacity"]["max_voltage_drop_percent"], i["max_drop_pct"] <= e["capacity"]["max_voltage_drop_percent"] + 1e-6),
            R("H", "electricity", "max transformer loading", round(i["max_loading"], 3),
              e["validation"]["max_transformer_loading_percent"] / 100,
              i["max_loading"] <= e["validation"]["max_transformer_loading_percent"] / 100 + 1e-9),
            R("G5", "electricity", "normally-open ties (independent loops when closed)", i["mu_with_ties"], "-", True),
            R("G3", "electricity", "max LV feeders per transformer", i["max_feeders"],
              e["graph"]["source_max_degree"], i["max_feeders"] <= e["graph"]["source_max_degree"])]


def gas_rows(gr_, packs, LP, MP):
    g = packs["gas"]
    i = gr_["info"]
    R = validate._row
    rows = [R("G10", "gas", "max LP zone reach (m)", round(i["max_zone_reach_m"], 1), g["graph"]["lp_zone_max_reach_m"],
              i["max_zone_reach_m"] <= g["graph"]["lp_zone_max_reach_m"]),
            R("H", "gas", "min LP pressure at taps (mbar)", round(i["lp"]["min_pressure_mbar"], 1),
              g["capacity"]["min_pressure_bar"] * 1000, i["lp"]["taps_below_min"] == 0),
            R("H", "gas", "max LP gradient (mbar/100 m)", round(i["lp"]["max_gradient_mbar_per_100m"], 3),
              g["capacity"]["max_pressure_drop_per_100m_mbar"],
              i["lp"]["max_gradient_mbar_per_100m"] <= g["capacity"]["max_pressure_drop_per_100m_mbar"] + 1e-6),
            R("H", "gas", "min MP pressure at regulators (bar g)", round(i["mp"]["min_regulator_inlet_bar_g"], 3),
              g["capacity"]["min_mp_pressure_bar"], i["mp"]["min_regulator_inlet_bar_g"] >= g["capacity"]["min_mp_pressure_bar"]),
            R("H", "gas", "max MP velocity (m/s)", round(i["mp"]["max_velocity_ms"], 2), g["capacity"]["max_velocity_mp_ms"],
              i["mp"]["max_velocity_ms"] <= g["capacity"]["max_velocity_mp_ms"]),
            R("H", "gas", "Renouard validity max Q/D", round(i["lp"]["max_q_over_d"], 3), g["capacity"]["renouard_max_q_over_d"],
              i["lp"]["max_q_over_d"] <= g["capacity"]["renouard_max_q_over_d"]),
            R("G11", "gas", "min district regulator spacing (m)", round(i["min_regulator_spacing_m"], 1),
              g["topology"]["min_regulator_spacing_m"], i["min_regulator_spacing_m"] >= g["topology"]["min_regulator_spacing_m"]),
            R("H", "gas", "max regulator load (m3/h)", round(i["max_regulator_load_m3h"], 1),
              g["components"]["regulators"]["max_capacity_m3_h"],
              i["max_regulator_load_m3h"] <= g["components"]["regulators"]["max_capacity_m3_h"])]
    for k, N in (("LP", LP), ("MP", MP)):
        sp = g["components"]["valves"]["required_spacing_m"]
        rows.append(validate.check_spacing("gas", N, sp, dsg.place_along_runs(N, sp), what=f"{k} valve"))
    sp = g["safety"]["emergency_shutoff_spacing_m"]
    rows.append(validate.check_spacing("gas", MP, sp, dsg.place_along_runs(MP, sp), what="MP shut-off"))
    return rows


def crossing_rows(cross, cr):
    R = validate._row
    rows = [R("C1", "all", "non-perpendicular crossings", 0, 0, True,
              "services and legs are drawn perpendicular to the crossed street")]
    for u in ("water", "electricity", "gas"):
        rows.append(R("C2", u, "service crossings within 10 m of a junction on short edges",
                      cross["c2_unavoidable"].get(u, 0), "-", True,
                      f"{cross['c2_shifted'].get(u, 0)} legs shifted to 10 m"))
        rows.append(R("C3", u, "cased road crossings", cross["casings"].get(u, 0), "-", True,
                      f"trenchless {cross['trenchless'].get(u, 0)}"))
        rows.append(R("C6", u, "road crossings (services + junction legs)",
                      cross["road_service"].get(u, 0) + cross["road_junction"].get(u, 0), "min", True,
                      f"services {cross['road_service'].get(u, 0)}, legs {cross['road_junction'].get(u, 0)}"))
    for pair, d in cross["deepening"].items():
        rows.append(R("C4", pair, "inter-utility crossings needing local lowering", d["deepening_points"], "-", True,
                      f"lower by {d['lowering_m']} m at {d['crossings']} crossings"))
    av = cross["avoidability"]
    rows.append(R("C5/G9", "all", "street edges whose violations exceed the achievable minimum",
                  av["edges_above_minimum"], 0, av["edges_above_minimum"] == 0,
                  "exact test over all joint side / dual choices, priority gas clearance > verge > service length"))
    rows.append(R("C5", "all", "street sides whose mains need more than the verge width",
                  av["current"]["verge"], av["minimum"]["verge"], av["current"]["verge"] == av["minimum"]["verge"],
                  "limit = minimum achievable on those edges; compact lanes in kerb order " + " / ".join(cross["kerb_order"])))
    rows.append(R("C5", "gas", "buildings closer than the clearance to the gas main beside them",
                  av["current"]["gas"], av["minimum"]["gas"], av["current"]["gas"] == av["minimum"]["gas"],
                  f"limit = minimum achievable; {cross['gas_edges_blocked_both_sides']} street edges have "
                  "buildings within 1.5 m of the gas lane on both sides"))

    return rows


def demand_summary(b):
    s = b[b["served"]]
    return dict(buildings=len(b), served=len(s), residents=float(s.loc[s.land_use == "residential", "occupants"].sum()),
                units=int(s["units"].sum()), water_m3_day=float(s["water_lpd"].sum() / 1000),
                water_peak_ls=float(s["water_peak_ls"].sum()), elec_kw_undiversified=float(s["elec_kw"].sum()),
                gas_peak_m3h=float(s["gas_peak_m3h"].sum()),
                land_use=s["land_use"].value_counts().to_dict(), land_use_source=s["land_use_source"].value_counts().to_dict(),
                floors_source=s["floors_source"].value_counts().to_dict(), hazard_buildings=int(b["hazard"].sum()))


def _clean(d):
    if isinstance(d, dict):
        return {str(k): _clean(v) for k, v in d.items() if k not in ("components", "legs", "tour")}
    if isinstance(d, (list, tuple)):
        return [_clean(x) for x in d][:2000]
    if isinstance(d, (np.floating, np.integer)):
        return d.item()
    return d


def edges_gdf(N, C, cols, crs):
    rec = []
    for a, b, d in N.edges(data=True):
        ba, bb = validate.design_base(a), validate.design_base(b)
        if d.get("eid") is not None and C.G.has_edge(ba, bb):
            g = C.geometry(ba, bb)
        else:
            g = LineString([N.nodes[a]["xy"], N.nodes[b]["xy"]])
        rec.append({**{c: d.get(c) for c in cols}, "geometry": g})
    return gpd.GeoDataFrame(rec, crs=crs)


def nodes_gdf(N, cols, crs, keep=None):
    rec = [{"node": str(x), **{c: N.nodes[x].get(c) for c in cols}, "geometry": Point(N.nodes[x]["xy"])}
           for x in N if keep is None or keep(x)]
    return gpd.GeoDataFrame(rec, crs=crs)


def write_layers(out, crs, G, C, S, b, rows, Nw, wi, Fe, er, LP, MP, gr_, gate, infeed, cor):
    out = Path(out)
    if out.exists():
        out.unlink()
    save = lambda g, name: g.to_file(out, layer=name, driver="GPKG") if len(g) else None
    save(edges_gdf(Nw, C, ["eid", "road_class", "length", "diameter_mm", "flow_ls", "velocity_ms", "headloss_m",
                           "gradient_m_per_km", "hierarchy", "loop_closure", "split_link"], crs), "water_mains")
    save(nodes_gdf(Nw, ["z", "head_m", "pressure_mwc"], crs), "water_nodes")
    save(gpd.GeoDataFrame([{"geometry": Point(G.nodes[s]["xy"]), "name": f"Inflow {i + 1}"}
                           for i, s in enumerate(wi["sources"])], crs=crs), "water_inflows")
    for key, layer in (("hydrants", "water_hydrants"), ("valves", "water_valves")):
        pts = []
        for ri, s, (a, bb), off in wi["components"][key]:
            g = LineString([Nw.nodes[a]["xy"], Nw.nodes[bb]["xy"]])
            pts.append({"geometry": g.interpolate(min(off, g.length))})
        save(gpd.GeoDataFrame(pts, crs=crs), layer)
    save(edges_gdf(Fe, C, ["eid", "road_class", "length", "size_mm2", "current_a", "hierarchy", "cell", "split_link"], crs)
         .assign(cell=lambda d: d["cell"].astype(str)), "elec_lv")
    tr = [{"geometry": Point(G.nodes[s]["xy"]), "site": str(s), **{k: v for k, v in t.items()}}
          for s, t in er["transformers"].items()]
    save(gpd.GeoDataFrame(tr, crs=crs).assign(mv_customer=lambda d: d.get("mv_customer", pd.Series(dtype=str)).astype(str)),
         "elec_substations")
    ties = [{"geometry": C.geometry(a, bb), "status": "normally_open"} for a, bb, d in er["ties"].values()]
    save(gpd.GeoDataFrame(ties, crs=crs), "elec_lv_ties")
    rings = []
    for k, r in enumerate(er["rings"]):
        for j, leg in enumerate(r["legs"]):
            pts = [G.nodes[x]["xy"] for x in leg]
            if len(pts) > 1:
                rings.append({"geometry": LineString(pts), "ring": k, "leg": j, "nop": j == r["nop_leg"],
                              "kva": r["kva"]})
    save(gpd.GeoDataFrame(rings, crs=crs), "elec_mv_rings")
    save(edges_gdf(LP, C, ["eid", "length", "diameter_mm", "flow_m3h", "drop_mbar", "zone", "split_link"], crs)
         .assign(zone=lambda d: d["zone"].astype(str)), "gas_lp")
    save(edges_gdf(MP, C, ["eid", "length", "diameter_mm", "flow_m3h", "velocity_ms", "split_link"], crs), "gas_mp")
    save(gpd.GeoDataFrame([{"geometry": Point(G.nodes[s]["xy"]), "site": str(s)} for s in gr_["sites"]], crs=crs),
         "gas_regulators")
    save(gpd.GeoDataFrame([{"geometry": Point(G.nodes[gate]["xy"]), "name": "City gate (hypothetical)"}], crs=crs),
         "gas_city_gate")
    save(gpd.GeoDataFrame([{"geometry": Point(G.nodes[infeed]["xy"]), "name": "MV infeed 20 kV"}], crs=crs),
         "elec_mv_infeed")
    svc = []
    for r in rows:
        g = LineString(r["path"]) if len(r["path"]) > 1 and r["dist"] > 0 else LineString([r["path"][0], r["path"][0]])
        for u in ("water", "electricity", "gas"):
            crosses = r["side"] not in cor.sides_of(u, r["eid"])
            svc.append({"geometry": g, "bldg_id": r["bldg_id"], "utility": u, "crosses_road": crosses,
                        "length_m": cor.service_length(u, r), "shared_with": r.get("host")})
    save(gpd.GeoDataFrame(svc, crs=crs), "services")
    # side of every main on every street edge
    side = []
    for u in cor.U:
        for e in cor.E[u]:
            side.append({"geometry": S.edge_geom(e), "utility": u, "eid": e, "side": cor.sigma[u][e],
                         "dual": cor.dual[u][e], "road_class": S.edges[e][3]["road_class"]})
    save(gpd.GeoDataFrame(side, crs=crs), "corridor_sides")
    faces = gpd.GeoDataFrame(geometry=S.faces(), crs=crs)
    save(faces, "street_faces")


def register(rows, bi, C, er, gr_, wi, cor):
    """Customer register keyed by the V8 address."""
    out = []
    owner_e, owner_g = er["owner"], gr_["owner"]
    mvc = {t.get("mv_customer") for t in er["transformers"].values()}
    for r in rows:
        bid = r["bldg_id"]
        t = C.tap_vertex[bid]
        x = bi.loc[bid]
        out.append(dict(bldg_id=bid, full_address=x.get("full_address"), street=x.get("street"),
                        house_number=x.get("house_number"), quartier=x.get("quartier"), land_use=x["land_use"],
                        units=int(x["units"]), water_peak_ls=round(float(x["water_peak_ls"]), 4),
                        elec_kw=round(float(x["elec_service_kw"]), 2), gas_peak_m3h=round(float(x["gas_peak_m3h"]), 4),
                        transformer=("MV customer" if bid in mvc else str(owner_e.get(t))),
                        gas_regulator=str(owner_g.get(t)), side=r["side"], tap_eid=r["eid"],
                        water_service_m=round(cor.service_length("water", r), 1),
                        elec_service_m=round(cor.service_length("electricity", r), 1),
                        gas_service_m=round(cor.service_length("gas", r), 1),
                        water_crosses_road=r["side"] not in cor.sides_of("water", r["eid"]),
                        elec_crosses_road=r["side"] not in cor.sides_of("electricity", r["eid"]),
                        gas_crosses_road=r["side"] not in cor.sides_of("gas", r["eid"])))
    return pd.DataFrame(out)
