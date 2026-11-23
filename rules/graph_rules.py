"""
Graph Rules - The topology rules of the utility rule packs written as
properties of graphs embedded in the plane.

Notation
  G0 = (V0, E0)      street graph, noded at every crossing; a plane graph
  N_u = (V_u, E_u)   network of utility u (water, electricity, gas), drawn
                     along the corridors of G0; services are pendant edges
  S_u, T_u           source vertices and demand terminals (building taps)
  c(N)               number of connected components
  mu(N)              cyclomatic number |E| - |V| + c(N)  (independent loops)
  F(N)               faces of the plane drawing; Euler: |V| - |E| + F = 1 + c
  Delta(N)           maximum vertex degree
  d_N(x, y)          shortest-path (network) distance, edge weight = length

Each rule names the rule-pack parameter it reads, so the values stay in
water_rules / electricity_rules / gas_rules / crossing_rules.  The engine
(engine/validate.py) implements one check per rule id.
"""


GRAPH_RULES = [
    {
        "id": "G1",
        "name": "plane_embedding",
        "statement": "N_u is a plane graph: the drawings of two edges meet only at a common end vertex.",
        "formal": "for all e != f in E_u: draw(e) ∩ draw(f) ⊆ ends(e) ∩ ends(f)",
        "parameter": None,
        "applies_to": ["water", "electricity", "gas"],
        "check": "check_plane",
        "repair": "node the crossing (insert a vertex) or re-route along the street graph",
    },
    {
        "id": "G2",
        "name": "planarity_bound",
        "statement": "N_u is planar (Kuratowski); hence |E_u| <= 3|V_u| - 6, and the street graph G0 is planar as well.",
        "formal": "N_u has no K5 or K3,3 subdivision; checked combinatorially (left-right / Boyer-Myrvold test)",
        "parameter": None,
        "applies_to": ["water", "electricity", "gas"],
        "check": "check_planarity",
        "repair": None,
    },
    {
        "id": "G3",
        "name": "max_degree",
        "statement": "Every vertex has degree at most the T-junction limit.",
        "formal": "Delta(N_u) <= topology.max_branching (= 3)",
        "parameter": "topology.max_branching",
        "applies_to": ["water", "electricity", "gas"],
        "check": "check_max_degree",
        "repair": "vertex splitting: replace a degree-k vertex by a path of k-2 degree-3 vertices, keeping the cyclic order of the rotation system (planarity preserved)",
    },
    {
        "id": "G4",
        "name": "source_connectivity",
        "statement": "Every terminal is connected to a source; every component holds at least one source.",
        "formal": "for all t in T_u there is s in S_u with d_N(s, t) < inf",
        "parameter": None,
        "applies_to": ["water", "electricity", "gas"],
        "check": "check_connectivity",
        "repair": "connect the component to the nearest source component by a shortest corridor path",
    },
    {
        "id": "G5",
        "name": "loop_policy",
        "statement": "The cyclomatic number obeys the loop policy. In a plane network the loops are the bounded faces: mu = F - 1.",
        "formal": "gas: mu(N) = 0 (forest);  electricity: mu(closed) = 0, ties are normally open;  water: mu(N) >= 1 allowed",
        "parameter": "topology.allow_loops",
        "applies_to": ["water", "electricity", "gas"],
        "check": "check_loops",
        "repair": "gas: delete the heaviest edge of each fundamental cycle; electricity: open one tie per cycle",
    },
    {
        "id": "G6",
        "name": "bridgeless_trunk",
        "statement": "For a looped utility the trunk subgraph has no bridge, so one pipe failure cannot cut the trunk.",
        "formal": "water with optimization.include_redundancy: bridges(N_primary) = ∅  (2-edge-connected)",
        "parameter": "optimization.include_redundancy",
        "applies_to": ["water"],
        "check": "check_bridgeless_trunk",
        "repair": "add the street chord that closes the shortest face around each bridge (edge augmentation)",
    },
    {
        "id": "G7",
        "name": "hierarchy_monotone",
        "statement": "Oriented from the sources, the hierarchy level never rises along the flow; services are pendant edges.",
        "formal": "for arcs (a->b), (b->c): level(b,c) >= level(a,b);  deg(building) = 1",
        "parameter": "topology.hierarchy_levels",
        "applies_to": ["water", "electricity", "gas"],
        "check": "check_hierarchy",
        "repair": "re-label by downstream demand share",
    },
    {
        "id": "G8",
        "name": "segment_length",
        "statement": "No edge is longer than the maximum segment; long edges are subdivided (degree-2 vertices keep G1-G3).",
        "formal": "for all e in E_u: len(e) <= topology.max_segment_length_m",
        "parameter": "topology.max_segment_length_m",
        "applies_to": ["water", "electricity", "gas"],
        "check": "check_segment_length",
        "repair": "edge subdivision",
    },
    {
        "id": "G9",
        "name": "service_geometry",
        "statement": "A service edge is perpendicular to the main at its tap, no longer than the service limit, and taps are distinct.",
        "formal": "angle(service, main) = 90° ± tol;  len(service) <= routing.max_service_distance;  |tap_i - tap_j| >= min spacing",
        "parameter": "routing.max_service_distance",
        "applies_to": ["water", "electricity", "gas"],
        "check": "check_services",
        "repair": "move the tap along the main (arc-length shift)",
    },
    {
        "id": "G10",
        "name": "supply_reach",
        "statement": "Each terminal lies within the reach of its supply point, measured along the network (graph-Voronoi cell radius).",
        "formal": "electricity: d_N(transformer(t), t) <= topology.transformer_max_reach_m;  gas LP zone: pressure >= min at all terminals",
        "parameter": "topology.transformer_max_reach_m",
        "applies_to": ["electricity", "gas"],
        "check": "check_reach",
        "repair": "add a supply point and recompute the graph-Voronoi partition",
    },
    {
        "id": "G11",
        "name": "component_spacing",
        "statement": "Components are placed at the prescribed arc-length spacing along the mains (metric graph).",
        "formal": "max gap between consecutive valves / hydrants / shut-offs along any main <= spacing",
        "parameter": "components.*.spacing_m",
        "applies_to": ["water", "electricity", "gas"],
        "check": "check_component_spacing",
        "repair": "insert a component at the arc-length position",
    },
    {
        "id": "G12",
        "name": "corridor_only",
        "statement": "Mains run along street edges; only service edges leave the street graph.",
        "formal": "every non-service edge of N_u is a sub-arc of an edge of G0 (offset to its side)",
        "parameter": "routing.prefer_roads",
        "applies_to": ["water", "electricity", "gas"],
        "check": "check_corridor",
        "repair": None,
    },
    {
        # Added 2026-11-23 (P5, Bukavu test): services are drawn as straight segments from the tap to
        # the nearest facade point; G1 checks them against each other and C-rules against streets, but
        # nothing checked them against the other buildings they pass.  On the design of P1-P4 for Gombe
        # 4,309 of 12,378 services cross another footprint.  Reported, not yet repaired.
        "id": "G13",
        "name": "services_clear_of_buildings",
        "statement": "A service line does not cross the footprint of a building other than the one it serves.",
        "formal": "for every service s(b): draw(s(b)) ∩ int(footprint(b')) = ∅ for all b' != b",
        "parameter": None,
        "applies_to": ["water", "electricity", "gas"],
        "check": "check_services_buildings",
        "repair": "route the service around the footprint (bent service along a passage) or share a neighbour's service",
    },
    {
        "id": "C1",
        "name": "right_angle_crossings",
        "statement": "Every crossing of a network edge with a street edge or another utility is perpendicular (RAC drawing).",
        "formal": "for every crossing point: |angle - 90°| <= crossing.geometry.angle_tolerance_deg",
        "parameter": "crossing.geometry.angle_tolerance_deg",
        "applies_to": ["water", "electricity", "gas"],
        "check": "check_crossing_angles",
        "repair": "insert a dog-leg so that the crossing leg is perpendicular",
    },
    {
        "id": "C2",
        "name": "crossing_away_from_junctions",
        "statement": "Service crossings keep a minimum distance from street vertices of degree >= 3.",
        "formal": "for service crossings x: min_{v in V0, deg v >= 3} |x - v| >= min_distance_from_junction_m",
        "parameter": "crossing.geometry.min_distance_from_junction_m",
        "applies_to": ["water", "electricity", "gas"],
        "check": "check_crossing_junction_distance",
        "repair": "move the tap along the main",
    },
    {
        "id": "C3",
        "name": "crossing_treatment",
        "statement": "Road crossings meet the cover depth and casing rules of the road class; rail and watercourse crossings are trenchless and cased.",
        "formal": "cover(x) >= road_crossing_cover_m[u];  casing(x) if road_classes[class].casing_required",
        "parameter": "crossing.road_crossing_cover_m",
        "applies_to": ["water", "electricity", "gas"],
        "check": "check_crossing_treatment",
        "repair": None,
    },
    {
        "id": "C4",
        "name": "vertical_clearance",
        "statement": "Where two utilities cross, their depth difference is at least the clearance; otherwise one is locally lowered.",
        "formal": "|depth_u(x) - depth_v(x)| >= vertical_clearance_m[u, v]",
        "parameter": "crossing.vertical_clearance_m",
        "applies_to": ["water", "electricity", "gas"],
        "check": "check_vertical_clearance",
        "repair": "local deepening of the lower-ranked utility",
    },
    {
        "id": "C5",
        "name": "horizontal_separation",
        "statement": "Parallel mains keep the pairwise separations of the rule packs; gas keeps its distance to buildings.",
        "formal": "|offset_u(e) - offset_v(e)| >= min_separation(u, v) on every shared corridor edge e",
        "parameter": "constraints.min_separation_*",
        "applies_to": ["water", "electricity", "gas"],
        "check": "check_horizontal_separation",
        "repair": "re-assign sides / offsets in the corridor",
    },
    {
        "id": "C6",
        "name": "minimum_road_crossings",
        "statement": "The number of carriageway crossings is minimised by choosing the side of each main (dual mains on primary roads).",
        "formal": "min over side assignments sigma: sum_services [side(building) != sigma(edge)] + penalty * side switches",
        "parameter": "crossing.optimization.minimise_road_crossings",
        "applies_to": ["water", "electricity", "gas"],
        "check": "report_road_crossings",
        "repair": None,
    },
]


def get_graph_rules(utility_type=None):
    """
    Get graph-theoretic rules, optionally filtered by utility

    Args:
        utility_type: 'water', 'electricity', 'gas' or None for all

    Returns:
        List of rule dictionaries
    """
    if utility_type is None:
        return list(GRAPH_RULES)
    return [r for r in GRAPH_RULES if utility_type in r["applies_to"]]
