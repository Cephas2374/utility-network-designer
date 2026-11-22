"""
Water Network Rules - Configuration for water distribution systems
"""


def get_water_rules():
    """
    Get rule pack for water distribution networks

    Returns:
        Dictionary with water network rules and parameters
    """
    return {
        "utility_type": "water",
        "display_name": "Water Distribution",
        "description": "Rules and parameters for water distribution network design",

        # Topology constraints
        "topology": {
            "max_branching": 3,  # Maximum derivations per node
            "allow_loops": True,  # Allow loop structures for reliability
            "hierarchy_levels": ["trunk", "secondary", "service"],
            "max_segment_length_m": 500,  # Maximum segment length before intermediate node
            "min_valve_spacing_m": 200,  # Isolation valve interval along mains (max gap between valves)
            "hydrant_spacing_m": 100  # Hydrant interval along mains for fire protection
        },

        # Graph-theoretic constraints (rules/graph_rules.py)
        "graph": {
            "loop_policy": "looped",  # mu >= 0; loops are the bounded faces of the plane network
            "bridgeless_trunk": True,  # Trunk must be 2-edge-connected (G6)
            "trunk_demand_share": 0.10,  # Edges carrying >= 10% of total demand form the trunk
            "max_loop_chord_m": 400,  # Longest street chord added to close a trunk loop
            "source_max_degree": 3,  # Sources obey the T-junction limit as well
            "inflow_single_connection": True,  # Each bulk inflow arrives through one transmission main (<= largest size)
            "min_tap_spacing_m": 1.0  # Distinct service taps along a main (G9)
        },

        # Capacity and hydraulic parameters
        "capacity": {
            "source_pressure_mwc": 50,  # Source pressure in meters water column
            "min_pressure_mwc": 20,  # Minimum required pressure at demand nodes
            "max_pressure_mwc": 80,  # Maximum static pressure (pipe class and leakage)
            "max_velocity_ms": 2.5,  # Maximum flow velocity (m/s)
            "design_velocity_ms": 1.0,  # Economic velocity used for the first sizing pass
            "hazen_williams_c": 120,  # Hazen-Williams roughness coefficient
            "fire_flow_ls": 15,  # Fire flow requirement (L/s)
            "fire_flow_min_diameter_mm": 100,  # Hydrants only on mains of at least this size
            "fire_residual_mwc": 15,  # Residual pressure at the hydrant while it delivers the fire flow (DVGW W 405: 1.5 bar)
            "peak_hour_factor": 2.5,  # Peak hour demand multiplier
            "max_pressure_drop_per_km_m": 10,  # Maximum pressure drop per km
            "solver": "gga",  # Global Gradient Algorithm (Todini-Pilati) for looped networks
            "solver_tolerance_ls": 1e-4,  # Convergence tolerance on flows
            "solver_max_iterations": 50
        },

        # Demand estimation parameters
        "demand": {
            "daily_consumption_lpd": 150,  # Liters per person per day
            "occupants_per_room": 1.5,  # Average occupants per room
            "rooms_per_apartment": 3,  # Default when the building has no room count
            "floor_area_per_apartment_m2": 100,  # Gross floor area of one dwelling unit
            "peak_factor": 2.5,  # Peak demand factor
            "land_use_multipliers": {
                "residential": 1.0,
                "commercial": 1.5,
                "industrial": 2.0,
                "public": 1.8,
                "mixed": 1.3,
                "agricultural": 0.5
            },
            "seasonal_factor": 1.2  # Summer vs winter demand variation
        },

        # Physical constraints
        "constraints": {
            "max_slope_percent": 15,  # Maximum acceptable slope
            "min_burial_depth_m": 0.8,  # Minimum pipe burial depth
            "max_burial_depth_m": 3.0,  # Maximum practical burial depth
            "min_separation_from_sewer_m": 3.0,  # Separation from sewer lines
            "min_separation_from_gas_m": 1.0,  # Separation from gas lines
            "min_separation_from_electricity_m": 0.5,  # Separation from power cables
            "crossing_clearance_m": 0.3  # Vertical clearance for crossings
        },

        # Standard pipe sizes (mm)
        "standard_sizes": [50, 75, 100, 150, 200, 250, 300, 400, 500, 600],
        "service_sizes": [20, 25, 32, 40, 50],  # Service connection pipe sizes (mm)

        # Materials and costs
        "materials": {
            "default": "HDPE",
            "options": ["HDPE", "PVC", "DI", "Steel"],
            "costs_per_m": {  # Approximate costs per meter (local currency)
                "HDPE": {"25": 8, "50": 20, "75": 30, "100": 40, "150": 70, "200": 100, "250": 140,
                         "300": 185, "400": 290, "500": 420, "600": 560},
                "PVC": {"50": 15, "100": 30, "150": 50, "200": 75},
                "DI": {"100": 80, "150": 120, "200": 170, "300": 250}
            },
            "trench_cost_per_m": 25  # Excavation and reinstatement in a verge (local currency)
        },

        # Components
        "components": {
            "valves": {
                "types": ["gate", "butterfly", "check"],
                "required_spacing_m": 200,  # Max interval along a main
                "required_at_branches": True  # A valve on each branch of a T-junction
            },
            "hydrants": {
                "types": ["pillar", "underground"],
                "spacing_m": 100,
                "required_pressure_mwc": 25
            },
            "chambers": {
                "types": ["valve_chamber", "junction_chamber", "meter_chamber"],
                "required_at_branches": True
            },
            "meters": {
                "types": ["bulk", "customer"],
                "required_at_service": True,
                "bulk_meter_at_dma_inlets": True
            }
        },

        # Routing preferences
        "routing": {
            "prefer_roads": True,
            "road_corridor_weight": 1.0,
            "off_road_weight": 3.0,  # Penalty for routing off roads
            "max_direct_distance": 500,  # Maximum direct connection distance
            "max_service_distance": 100  # Maximum service connection length
        },

        # Optimization parameters
        "optimization": {
            "method": "steiner",  # 'mst', 'steiner', 'shortest_path'
            "objective": "cost",  # 'cost', 'length', 'reliability'
            "include_redundancy": True,  # Add redundant paths for reliability
            "dma_segmentation": True,  # Create District Metered Areas
            "dma_target_connections": 1500  # Connections per DMA (graph partition of the dual)
        },

        # Validation rules
        "validation": {
            "check_pressure": True,
            "check_velocity": True,
            "check_fire_flow": True,
            "check_topology": True,
            "allow_zero_pressure": False
        }
    }
