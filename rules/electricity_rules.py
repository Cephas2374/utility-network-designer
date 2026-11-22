"""
Electricity Network Rules - Configuration for electricity distribution systems
"""


def get_electricity_rules():
    """
    Get rule pack for electricity distribution networks

    Returns:
        Dictionary with electricity network rules and parameters
    """
    return {
        "utility_type": "electricity",
        "display_name": "Electricity Distribution",
        "description": "Rules and parameters for electricity distribution network design",

        # Topology constraints
        "topology": {
            "max_branching": 3,  # Maximum branches per node
            "allow_loops": True,  # Allow ring/loop configurations for reliability
            "hierarchy_levels": ["primary", "secondary", "service"],
            "max_segment_length_m": 300,  # Maximum span between poles
            "transformer_max_reach_m": 300,  # Maximum reach from transformer
            "protection_zone_coordination": True
        },

        # Graph-theoretic constraints (rules/graph_rules.py)
        "graph": {
            "loop_policy": "radial_with_ties",  # Closed network is a forest; ties are normally open
            "lv_ties_between_cells": True,  # One normally-open LV tie per adjacent pair of transformer cells
            "mv_ring": True,  # MV backbone closed as an open ring through the transformers
            "mv_ring_max_kva": 10000,  # Installed transformer capacity per MV open ring (feeder pair)
            "source_max_degree": 6,  # LV feeders per transformer board (source exempt from max_branching)
            "min_tap_spacing_m": 1.0
        },

        # Capacity and electrical parameters
        "capacity": {
            "nominal_voltage_v": 400,  # Low voltage 3-phase (400V)
            "mv_voltage_v": 11000,  # Medium voltage (11kV)
            "power_factor": 0.9,  # Typical power factor
            "max_voltage_drop_percent": 5,  # Maximum voltage drop (%)
            "conductor_resistivity_ohmm": 1.724e-8,  # Copper resistivity
            "conductor_reactance_ohm_per_km": 0.08,  # LV cable series reactance
            "ambient_temperature_c": 30,  # Design ambient temperature
            "max_conductor_temperature_c": 70,  # Maximum operating temperature
            "ampacity_derating_factor": 0.9,  # Grouping / soil-temperature derating of the table below
            "min_main_size_mm2": 25,  # Smallest LV main cable
            "mv_customer_min_kva": 400,  # Customers above this load get their own substation on the MV ring
            "min_service_size_mm2": 10,  # Smallest service cable
            # restoration through a normally-open LV tie (paper P3)
            "emergency_voltage_drop_percent": 10,  # EN 50160: supply voltage within +-10 % of nominal
            "emergency_transformer_loading_percent": 100,  # pick-up within the nameplate rating
            "max_mv_voltage_drop_percent": 5  # MV ring fed from one end after a fault
        },

        # 20 kV ring cable: 3 x 1 x 240 mm2 Al XLPE 12/20 kV in trefoil, buried (IEC 60502-2).
        # No SNEL cable specification is available; the values are derived from the standards:
        "mv_cable": {
            "type": "3x1x240 Al XLPE 12/20 kV",
            # IEC 60228: R_DC(20 C) = 0.125 ohm/km; at 90 C: 0.125 (1 + 0.00403 * 70) = 0.160 ohm/km
            "r_ohm_per_km": 0.160,
            "x_ohm_per_km": 0.11,  # single-core cables in trefoil
            "ampacity_base_a": 365,  # catalogue value, ground at 20 C, one circuit
            "ground_temp_factor": 0.93,  # IEC 60502-2 correction for ground at 30 C (Kinshasa)
            "ampacity_a": 340,  # 365 * 0.93
            "cost_per_m": 85  # local currency units, cable only (trench as LV)
        },

        # Demand estimation parameters
        "demand": {
            "base_load_kw_per_apt": 3.5,  # Base load per apartment (kW)
            "load_per_m2": 0.02,  # Load per square meter (kW/m²)
            "floor_area_per_apartment_m2": 100,  # Gross floor area of one dwelling unit
            "diversity_factors": {  # Diversity based on number of units
                "1-5": 0.9,
                "6-10": 0.8,
                "11-20": 0.7,
                "21-50": 0.6,
                "50+": 0.5
            },
            "building_multipliers": {
                "residential": 1.0,
                "commercial": 2.0,
                "industrial": 4.0,
                "public": 1.5,
                "mixed": 1.5,
                "data_center": 10.0
            },
            "time_of_use_factor": 1.3,  # Peak vs average demand
            "future_growth_factor": 1.25  # Design for 25% future growth
        },

        # Physical constraints
        "constraints": {
            "max_slope_percent": 20,  # More flexible for overhead lines
            "min_burial_depth_m": 0.7,  # For underground cables
            "max_burial_depth_m": 1.5,
            "min_separation_from_water_m": 0.5,
            "min_separation_from_gas_m": 0.5,
            "min_clearance_overhead_m": 5.5,  # Overhead line clearance
            "pole_spacing_m": 40  # Typical pole spacing for overhead
        },

        # Standard cable sizes (mm²)
        "standard_sizes": [1.5, 2.5, 4, 6, 10, 16, 25, 35, 50, 70, 95, 120, 150, 185, 240, 300],

        # Cable current capacity (A) - simplified table
        "cable_capacity": {
            1.5: 19, 2.5: 26, 4: 35, 6: 46, 10: 63,
            16: 85, 25: 112, 35: 138, 50: 168, 70: 213,
            95: 258, 120: 299, 150: 344, 185: 392, 240: 461, 300: 530
        },

        # Materials and costs
        "materials": {
            "default": "XLPE",  # Cross-linked polyethylene cable
            "options": ["XLPE", "PVC", "Overhead_ACSR"],
            "conductor_types": ["Copper", "Aluminum"],
            "costs_per_m": {
                "XLPE_Cu": {"10": 15, "16": 22, "25": 35, "35": 48, "50": 70, "70": 95, "95": 130,
                            "120": 165, "150": 205, "185": 250, "240": 320, "300": 400},
                "XLPE_Al": {"16": 12, "35": 28, "70": 55, "120": 100},
                "Overhead": {"35": 8, "50": 12, "95": 20}
            },
            "trench_cost_per_m": 20
        },

        # Components
        "components": {
            "transformers": {
                "types": ["pad_mounted", "pole_mounted", "ground_mounted"],
                "standard_ratings_kva": [25, 50, 100, 160, 250, 315, 400, 500, 630, 800, 1000],
                "typical_service_radius_m": 300,
                "efficiency_percent": 98
            },
            "switchgear": {
                "types": ["circuit_breaker", "fuse", "recloser"],
                "required_at_branches": True,
                "coordination_required": True
            },
            "meters": {
                "types": ["single_phase", "three_phase", "smart_meter"],
                "required_at_service": True
            },
            "poles": {
                "types": ["wood", "concrete", "steel"],
                "typical_height_m": 10,
                "spacing_m": 40
            },
            "joints": {
                "types": ["straight", "tee"],
                "max_cable_drum_length_m": 300  # A straight joint at least every drum length
            }
        },

        # Routing preferences
        "routing": {
            "prefer_roads": True,
            "road_corridor_weight": 1.0,
            "off_road_weight": 2.5,
            "max_direct_distance": 300,
            "max_service_distance": 50,
            "underground_vs_overhead": "mixed",  # 'underground', 'overhead', 'mixed'
            "avoid_crossings": True
        },

        # Protection and reliability
        "protection": {
            "fault_level_calculation": True,
            "protection_coordination": True,
            "grounding_required": True,
            "lightning_protection": True,
            "backup_feeders": False,  # Redundant feeders
            "automatic_reclosing": True
        },

        # Optimization parameters
        "optimization": {
            "method": "steiner",
            "objective": "voltage_stability",  # 'cost', 'voltage_stability', 'reliability'
            "include_redundancy": True,  # Ring mains for reliability
            "load_balancing": True,  # Balance loads across phases
            "feeder_zone_segmentation": True
        },

        # Validation rules
        "validation": {
            "check_voltage_drop": True,
            "check_cable_capacity": True,
            "check_transformer_loading": True,
            "check_protection_coordination": True,
            "check_topology": True,
            "check_phase_balance": False,  # Simplified - would need phase modeling
            "max_transformer_loading_percent": 80
        }
    }
