"""
Crossing Rules - Road, rail, watercourse and inter-utility crossings,
and the lateral layout of utilities inside a street corridor.

The rules complement the per-utility rule packs (water_rules, electricity_rules,
gas_rules). Horizontal separations and burial depths are NOT repeated here: they
are read from the utility packs ('constraints' namespace) so that each value has
a single source.

Graph-theoretic reading (see graph_rules.py):
  * A crossing is a point where the drawing of one network edge meets the
    drawing of a street edge (road crossing) or of another utility's edge
    (inter-utility crossing).  Every single-utility network must be a plane
    graph (no self-crossings); the union of the utility networks with the
    street graph is a right-angle-crossing (RAC) drawing: every crossing is
    perpendicular within the tolerance below.
  * The number of road crossings is the crossing number of the network
    drawing against the street drawing and is minimised by the side
    assignment of mains (corridor.py / crossings.py).
"""


def get_crossing_rules():
    """
    Get rule pack for crossings and corridor layout

    Returns:
        Dictionary with crossing rules and parameters
    """
    return {
        "rule_set": "crossings",
        "display_name": "Crossings and Corridor Layout",
        "description": "Rules for networks crossing roads, railways, watercourses and other utilities",

        # Geometry of every crossing (roads, rails, other utilities)
        "geometry": {
            "crossing_angle_deg": 90,  # Crossings are perpendicular (RAC drawing)
            "angle_tolerance_deg": 10,  # Accepted deviation: 80-100 degrees
            "min_distance_from_junction_m": 10,  # Service crossings kept out of intersections
            "min_distance_from_fitting_m": 1.0,  # No crossing over valves, joints, chambers
            "max_crossings_per_service": 1,  # A service line crosses at most one carriageway
            "allow_diagonal_junction_crossing": False  # Never cross an intersection diagonally
        },

        # Road classes: carriageway width and crossing treatment.
        # Widths are planning defaults for Gombe's OSM classes.
        "road_classes": {
            "primary":     {"width_m": 20.0, "verge_m": 3.0, "casing_required": True,  "trenchless": True,  "dual_mains": True,  "cost_factor": 3.0},
            "secondary":   {"width_m": 14.0, "verge_m": 2.5, "casing_required": True,  "trenchless": False, "dual_mains": True, "cost_factor": 2.0},
            "tertiary":    {"width_m": 10.0, "verge_m": 2.0, "casing_required": False, "trenchless": False, "dual_mains": True, "cost_factor": 1.5},
            "residential": {"width_m": 8.0,  "verge_m": 1.5, "casing_required": False, "trenchless": False, "dual_mains": True, "cost_factor": 1.0},
            "service":     {"width_m": 6.0,  "verge_m": 1.0, "casing_required": False, "trenchless": False, "dual_mains": True, "cost_factor": 1.0}
        },

        # Minimum cover below the road surface at a carriageway crossing (m).
        # Deeper than the verge burial depths of the utility packs because of traffic loads.
        "road_crossing_cover_m": {
            "water": 1.0,
            "electricity": 0.8,
            "gas": 1.2
        },

        # Casing (sleeve / duct) at road crossings
        "casing": {
            "extend_beyond_carriageway_m": 0.5,  # Sleeve protrudes past each kerb
            "diameter_ratio": 1.5  # Casing diameter / carrier diameter
        },

        # Railway crossings
        "railways": {
            "avoid": True,
            "cost_factor": 10.0,  # Routing penalty on corridor edges crossing a railway
            "cover_m": 1.5,  # Cover below rail level
            "casing_required": True,
            "trenchless": True
        },

        # Watercourse crossings (rivers, canals, drains)
        "watercourses": {
            "avoid": True,
            "cost_factor": 20.0,
            "allowed_methods": ["bridge_attachment", "trenchless"],
            "cover_below_bed_m": 1.5
        },

        # Inter-utility crossings: minimum vertical clearance (m) between the
        # outer surfaces of two utilities that cross.  Symmetric.
        "vertical_clearance_m": {
            ("water", "gas"): 0.3,
            ("water", "electricity"): 0.3,
            ("gas", "electricity"): 0.3,
            ("water", "sewer"): 0.5
        },
        # When two utilities cross and their burial depths differ by less than the
        # clearance, the deeper-ranked one is locally lowered.  Ranking (deepest first):
        "crossing_depth_order": ["gas", "water", "electricity"],

        # vertical design at crossings (paper P4); depths are to the crown of the pipe or cable
        "vertical": {
            "outer_dimension_m": {  # pipes: their nominal diameter (from the design) is used for mains
                "electricity_main": 0.08,  # LV cable with its bedding layer
                "electricity_service": 0.03,
                "water_service": 0.04,  # 32 mm PE service with tolerance
                "gas_service": 0.032
            },
            "max_gradient": 0.1,  # largest change of depth per metre of main or service (1:10)
            "trench_extra_width_m": 0.4,  # trench width = outer dimension + this (pipes), cables 0.4 m
            "service_trench_width_m": 0.3
        },

        # Lateral layout of a street corridor.  On each side of the street axis the
        # utilities present keep this order from the kerb outwards; the offsets are
        # the cumulative pairwise separations of the utility packs.  Gas nearest the
        # kerb (farthest from buildings), water nearest the buildings:
        #   gas -0.5 m- electricity -0.5 m- water, and gas-water = 1.0 m.
        "corridor_layout": {
            # Lanes from the kerb outwards, by decreasing depth (paper P4, Proposition "compatible orders"):
            # a service then never forces a main on its own side down.  P1 used gas, electricity,
            # water (reproduce with run_gombe.py --legacy-kerb-order).
            "kerb_order": ["gas", "water", "electricity"],
            "min_kerb_offset_m": 0.5,  # First utility at least 0.5 m behind the kerb
            "junction_disk_margin_m": 2.0  # Junction disk radius = half road width + verge + margin
        },

        # Cost of one carriageway crossing, as an equivalent length of main (m).
        # Used to decide the side of each main and whether a second (dual) main
        # on the other side of the street is cheaper than the crossings it saves.
        "crossing_cost_equivalent_m": {
            "primary": 40.0,  # Trenchless, cased, traffic management
            "secondary": 20.0,
            "tertiary": 10.0,
            "residential": 6.0,
            "service": 4.0
        },

        # Optimisation of crossings
        "optimization": {
            "minimise_road_crossings": True,
            "inter_utility_crossing_weight": 0.5,  # Weight of a utility-utility crossing vs a road crossing
            "side_assignment_iterations": 10,  # Local search sweeps of the side assignment
            # Rule penalties of the side assignment (local search), in metres of main.  The final
            # joint repair and the validation use the strict priority gas clearance > verge >
            # service length > cost (lexicographic), so these weights only guide the search.
            "penalties": {
                "service_over_limit": 1000.0,  # G9: service longer than routing.max_service_distance
                "lane_outside_verge": 1000.0,  # C5: utilities on one side need more than the verge width
                "gas_building_clearance": 5000.0  # C5 (safety, outranks G9): gas main too close to a building
            }
        }
    }
