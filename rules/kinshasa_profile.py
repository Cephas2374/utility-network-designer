"""
Kinshasa Profile - Local overrides of the rule packs for Gombe (Kinshasa, DR Congo)

Only values that the local context changes are overridden; everything else
keeps the defaults of water_rules / electricity_rules / gas_rules.  Each
override carries its reason so the paper can list them (Table: local profile).
"""


def get_kinshasa_profile():
    """
    Get local overrides for Kinshasa

    Returns:
        Dictionary {utility_type: {namespace: {key: value}}} and a list of reasons
    """
    overrides = {
        "water": {
            # REGIDESO supplies central Kinshasa, including Gombe, from the
            # Ngaliema treatment plant west of the commune.
            "sources": {"name": "Ngaliema (REGIDESO) inflow", "side": "west"},
            "demand": {"seasonal_factor": 1.0},  # Tropical climate, no summer peak in the design day
        },
        "common": {
            # Resident population used to calibrate the floor area per dwelling
            "target_population": 89080,
            "target_population_source": "Gombe commune, 2015 estimate",
        },
        "electricity": {
            # SNEL is merging its 6.6 kV and 30 kV systems into 20 kV.
            "capacity": {"mv_voltage_v": 20000, "ambient_temperature_c": 30,
                         "conductor_reactance_ohm_per_km": 0.08},
            "demand": {"load_per_m2": 0.05},  # Air-conditioned offices and shops in the business district
        },
        "gas": {
            # No piped gas network exists in Kinshasa: gas is a prospective
            # scenario (cooking only, no space heating in a tropical climate).
            "demand": {"heating_load_m3_per_apt": 0.0, "seasonal_factor": 1.0},
            "capacity": {"relative_density": 0.6},
            "sources": {"name": "Hypothetical city gate", "side": "east"},
        },
    }

    reasons = [
        ("water", "sources", "Gombe is supplied from the Ngaliema plant (REGIDESO), west of the commune"),
        ("common", "target_population", "89,080 residents (2015 est.): calibrates the floor area per dwelling"),
        ("water", "demand.seasonal_factor", "1.0: equatorial climate, no seasonal design peak"),
        ("electricity", "capacity.mv_voltage_v", "20 kV: SNEL target MV level in Kinshasa (6.6/30 kV being merged)"),
        ("electricity", "capacity.conductor_reactance_ohm_per_km", "0.08 ohm/km: LV cable reactance for the voltage drop"),
        ("electricity", "demand.load_per_m2", "50 W/m2: air-conditioned non-residential floor area"),
        ("gas", "demand.heating_load_m3_per_apt", "0: no space heating in Kinshasa; cooking load only"),
        ("gas", "demand.seasonal_factor", "1.0: no heating season"),
        ("gas", "sources", "No existing network: city gate placed at the eastern (port/industrial) edge, hypothetical"),
    ]
    return {"overrides": overrides, "reasons": reasons}


def apply_profile(rules, utility_type, profile=None):
    """
    Return a copy of a rule pack with the Kinshasa overrides applied

    Args:
        rules: Rule pack dictionary
        utility_type: 'water', 'electricity' or 'gas'
        profile: Profile from get_kinshasa_profile() (default: loaded)
    """
    import copy

    profile = profile or get_kinshasa_profile()
    merged = copy.deepcopy(rules)
    for namespace, values in profile["overrides"].get(utility_type, {}).items():
        merged.setdefault(namespace, {}).update(values)
    return merged
