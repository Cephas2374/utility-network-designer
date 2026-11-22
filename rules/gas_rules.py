"""
Gas Network Rules - Configuration for gas distribution systems
"""


def get_gas_rules():
    """
    Get rule pack for gas distribution networks

    Returns:
        Dictionary with gas network rules and parameters
    """
    return {
        "utility_type": "gas",
        "display_name": "Gas Distribution",
        "description": "Rules and parameters for gas distribution network design",

        # Topology constraints
        "topology": {
            "max_branching": 3,  # Maximum derivations per node
            "allow_loops": False,  # Generally tree-structured for gas
            "hierarchy_levels": ["transmission", "distribution", "service"],
            "max_segment_length_m": 400,
            "min_regulator_spacing_m": 500,
            "safety_zone_buffer_m": 10  # Safety buffer around gas lines
        },

        # Graph-theoretic constraints (rules/graph_rules.py)
        "graph": {
            "loop_policy": "forest",  # mu = 0 in both the MP and the LP networks
            "two_stage": True,  # MP backbone (inlet pressure) -> district regulators -> LP zones
            "lp_zone_max_reach_m": 1000,  # Network radius of an LP zone around its district regulator
            "source_max_degree": 3,
            "min_tap_spacing_m": 1.0
        },

        # Capacity and pressure parameters
        "capacity": {
            "inlet_pressure_bar": 4.0,  # Medium pressure inlet (bar)
            "min_pressure_bar": 0.02,  # Minimum delivery pressure (bar) - low pressure
            "min_mp_pressure_bar": 1.0,  # Minimum MP pressure at a district regulator inlet
            "max_pressure_drop_per_100m_mbar": 1.0,  # Maximum pressure drop per 100m
            "gas_density_kg_m3": 0.8,  # Natural gas density
            "relative_density": 0.6,  # Density relative to air (Renouard formulas)
            "gas_viscosity_pas": 1.1e-5,  # Gas dynamic viscosity
            "compressibility_factor": 0.998,  # Gas compressibility
            "renouard_mp_coefficient": 48600,  # P1^2 - P2^2 [bar abs^2] = c s L[km] Q^1.82 D[mm]^-4.82
            "renouard_lp_coefficient": 23200,  # P1 - P2 [mbar] = c s L[m] Q^1.82 D[mm]^-4.82
            "renouard_max_q_over_d": 150,  # Validity limit Q[m3/h] / D[mm]
            "pipe_sdr": 11,  # PE100 SDR 11: internal diameter = OD (1 - 2/SDR)
            "max_velocity_mp_ms": 20,  # Velocity limits (noise, erosion)
            "max_velocity_lp_ms": 6
        },

        # Demand estimation parameters
        "demand": {
            "heating_load_m3_per_apt": 5.0,  # m³/day per apartment for heating
            "cooking_load_m3_per_apt": 0.5,  # m³/day per apartment for cooking
            "floor_area_per_apartment_m2": 100,  # Gross floor area of one dwelling unit
            "seasonal_factor": 1.5,  # Winter vs summer demand ratio
            "peak_hour_factor": 1.8,  # Peak hour multiplier
            "building_multipliers": {
                "residential": 1.0,
                "commercial": 0.8,  # May have central heating
                "industrial": 3.0,
                "public": 1.5,
                "mixed": 1.2
            },
            "climate_adjustment": 1.0  # Adjustment for local climate
        },

        # Physical constraints
        "constraints": {
            "max_slope_percent": 10,
            "min_burial_depth_m": 1.0,  # Deeper burial for safety
            "max_burial_depth_m": 2.5,
            "min_separation_from_water_m": 1.0,
            "min_separation_from_electricity_m": 0.5,
            "min_separation_from_buildings_m": 1.5,  # Safety clearance
            "min_separation_from_ignition_sources_m": 5.0  # Safety requirement
        },

        # Standard pipe sizes (mm)
        "standard_sizes": [25, 32, 40, 50, 63, 75, 90, 110, 125, 160, 200, 250, 315],
        "service_sizes": [20, 25, 32],

        # Materials and costs
        "materials": {
            "default": "PE100",  # Polyethylene PE100
            "options": ["PE80", "PE100", "Steel"],
            "costs_per_m": {
                "PE100": {"25": 10, "32": 15, "40": 18, "50": 22, "63": 30, "75": 38, "90": 50, "110": 62,
                          "125": 80, "160": 105, "200": 140, "250": 185, "315": 240},
                "PE80": {"32": 12, "63": 25, "90": 42, "125": 70},
                "Steel": {"50": 60, "100": 120, "150": 180}
            },
            "trench_cost_per_m": 30
        },

        # Components
        "components": {
            "regulators": {
                "types": ["district_regulator", "service_regulator"],
                "required_at_pressure_zones": True,
                "inlet_pressure_bar": 4.0,
                "outlet_pressure_bar": 0.1,
                "max_capacity_m3_h": 1500  # Capacity of one district regulator station
            },
            "valves": {
                "types": ["ball_valve", "gate_valve", "emergency_shutoff"],
                "required_spacing_m": 300  # Max interval along a main
            },
            "meters": {
                "types": ["diaphragm", "rotary", "turbine"],
                "required_at_service": True
            },
            "vents": {
                "required_at_high_points": True,
                "spacing_m": 500
            }
        },

        # Routing preferences
        "routing": {
            "prefer_roads": True,
            "road_corridor_weight": 1.0,
            "off_road_weight": 5.0,  # Higher penalty due to safety concerns
            "max_direct_distance": 300,
            "max_service_distance": 30,  # Shorter service connections
            "avoid_dense_areas": True  # Prefer less populated routes
        },

        # Safety requirements
        "safety": {
            "emergency_shutoff_spacing_m": 500,
            "pressure_monitoring_required": True,
            "leak_detection_required": True,
            "hazard_zones": ["schools", "hospitals", "public_buildings"],
            "min_clearance_from_hazards_m": 50,
            "hazard_clearance_applies_to": ["district_regulator"],  # Above-ground stations, not buried pipe
            "cathodic_protection_required": True  # For steel pipes
        },

        # Optimization parameters
        "optimization": {
            "method": "mst",  # Simpler tree structure for gas
            "objective": "safety_cost",  # Balance safety and cost
            "include_redundancy": False,  # Generally no redundancy in gas
            "pressure_zone_segmentation": True
        },

        # Validation rules
        "validation": {
            "check_pressure": True,
            "check_pressure_drop": True,
            "check_safety_clearances": True,
            "check_topology": True,
            "enforce_no_loops": True,
            "check_regulator_placement": True
        }
    }
