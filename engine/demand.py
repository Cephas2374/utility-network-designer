"""
Demand estimation per building, from footprints, OSM tags, POIs and the
V8 addresses.

Floors come from `building:levels` where mapped, otherwise from the footprint
area class (Table "demand assumptions" of the paper).  One dwelling unit is
`floor_area_per_apartment_m2` of gross floor area; non-residential floor area
is converted to unit equivalents the same way, so the electrical diversity
factor can be applied to the number of units downstream of every network edge.
"""
import re
from functools import lru_cache

import geopandas as gpd
import numpy as np
import pandas as pd

from .config import MIN_SERVED_AREA_M2

PUBLIC_TAGS = {"school", "hospital", "university", "college", "government", "public", "civic", "church",
               "mosque", "temple", "kindergarten", "clinic"}
COMMERCIAL_TAGS = {"commercial", "retail", "office", "hotel", "supermarket", "kiosk"}
RESIDENTIAL_TAGS = {"house", "apartments", "residential", "terrace", "detached", "semidetached_house", "bungalow"}
HAZARD_TAGS = {"school", "hospital", "university", "college", "kindergarten", "clinic"}

PUBLIC_AMENITY = {"school", "hospital", "university", "college", "kindergarten", "doctors", "clinic",
                  "place_of_worship", "townhall", "police", "courthouse", "embassy", "library"}
PUBLIC_OFFICE = {"government", "diplomatic", "administrative", "ngo", "foreign_national_agency"}
HAZARD_AMENITY = {"school", "hospital", "university", "college", "kindergarten", "clinic", "doctors"}
BUSINESS_QUARTIERS = {"Commerce", "Gare"}

# Gross floor area of one unit equivalent for non-residential uses (m2)
NONRES_UNIT_M2 = 100.0

# Floors from the footprint when building:levels is missing (planning assumption);
# industrial sheds and warehouses are taken as single-storey.
FLOOR_CLASSES = [(150, 1), (500, 2), (1500, 3), (np.inf, 4)]


def parse_levels(v):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return np.nan
    m = re.match(r"\s*(\d+(\.\d+)?)", str(v))
    return float(m.group(1)) if m else np.nan


def default_floors(area):
    for limit, floors in FLOOR_CLASSES:
        if area < limit:
            return floors
    return FLOOR_CLASSES[-1][1]


def classify_land_use(b, pois, landuse):
    """Land use per building: tag, then POI inside the footprint, then land-use polygon, then quartier."""
    lu = pd.Series("residential", index=b.index, dtype=object)
    src = pd.Series("default", index=b.index, dtype=object)

    by_quartier = b["quartier"].isin(BUSINESS_QUARTIERS)
    lu[by_quartier], src[by_quartier] = "mixed", "quartier"

    if len(landuse):
        cent = gpd.GeoDataFrame(geometry=b.geometry.representative_point(), index=b.index, crs=b.crs)
        j = gpd.sjoin(cent, landuse[["fclass", "geometry"]], predicate="within", how="inner")
        j = j[~j.index.duplicated()]
        m = {"commercial": "commercial", "retail": "commercial", "industrial": "industrial",
             "residential": "residential"}
        for i, f in j["fclass"].items():
            if f in m:
                lu[i], src[i] = m[f], "landuse"

    if len(pois):
        p = pois.copy()
        p["cls"] = None
        p.loc[p["amenity"].isin(PUBLIC_AMENITY) | p["office"].isin(PUBLIC_OFFICE), "cls"] = "public"
        p.loc[p["cls"].isna() & (p["shop"].notna() | p["office"].notna() | p["tourism"].isin(["hotel", "motel", "guest_house", "hostel"])
                                 | p["amenity"].isin(["restaurant", "bank", "cafe", "bar", "fast_food", "pharmacy", "fuel"])),
              "cls"] = "commercial"
        p["hazard"] = p["amenity"].isin(HAZARD_AMENITY) | p["healthcare"].notna()
        p = p[p["cls"].notna() | p["hazard"]]
        j = gpd.sjoin(p[["cls", "hazard", "geometry"]], b[["geometry"]].assign(geometry=b.geometry.buffer(2.0)),
                      predicate="within", how="inner")
        for i, grp in j.groupby("index_right"):
            cls = grp["cls"].dropna()
            if len(cls):
                lu[i] = "public" if (cls == "public").any() else "commercial"
                src[i] = "poi"
        hazard_poi = set(j.loc[j["hazard"], "index_right"])
    else:
        hazard_poi = set()

    tag = b["building"].fillna("yes")
    lu[tag.isin(PUBLIC_TAGS)], src[tag.isin(PUBLIC_TAGS)] = "public", "tag"
    lu[tag.isin(COMMERCIAL_TAGS)], src[tag.isin(COMMERCIAL_TAGS)] = "commercial", "tag"
    lu[tag == "industrial"], src[tag == "industrial"] = "industrial", "tag"
    lu[tag.isin(RESIDENTIAL_TAGS)], src[tag.isin(RESIDENTIAL_TAGS)] = "residential", "tag"
    hazard = tag.isin(HAZARD_TAGS) | b.index.isin(hazard_poi)
    return lu, src, hazard


def diversity_factor(n, table):
    """IEC-style diversity factor for n units from the rule-pack table {'1-5': 0.9, ..., '50+': 0.5}."""
    for key, f in table.items():
        if key.endswith("+"):
            if n >= int(key[:-1]):
                return f
        else:
            a, b = (int(x) for x in key.split("-"))
            if a <= n <= b:
                return f
    return min(table.values())


def diversity_factor_monotone(n, table):
    """Least factor g*(n) >= g(n) for which the diversified load n g*(n) never decreases with n:
    g*(n) = max_{m <= n} m g(m) / n (the least non-decreasing majorant of n g(n), divided by n).
    The stepwise table gives 10 units 8.0 unit loads but 11 units only 7.7 (paper P3)."""
    n = max(1, int(n))
    return _monotone_load(n, tuple(table.items())) / n


@lru_cache(maxsize=None)
def _monotone_load(n, items):
    table = dict(items)
    return max(m * diversity_factor(m, table) for m in range(1, n + 1))


def residential_units(floor_area, area_per_unit):
    return np.maximum(1, np.round(floor_area / area_per_unit)).astype(int)


def calibrate_unit_area(floor_area, persons_per_unit, target, lo=20.0, hi=2000.0):
    """Floor area per dwelling (m2) such that the modelled residents equal the target (bisection).

    The modelled population is non-increasing in the unit area, so bisection converges.
    """
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        pop = residential_units(floor_area, mid).sum() * persons_per_unit
        lo, hi = (mid, hi) if pop > target else (lo, mid)
    return 0.5 * (lo + hi)


def estimate(buildings, pois, landuse, packs, target_population=None):
    """Return the buildings GeoDataFrame with demand columns for the three utilities.

    If target_population is given, the floor area per dwelling is calibrated so that the
    residents of the served residential buildings equal it; the value is stored in b.attrs.
    """
    b = buildings.copy()
    b["area_m2"] = b.geometry.area
    lev = b["building_levels"].map(parse_levels)
    b["floors_source"] = np.where(lev.notna(), "osm", "area_class")
    b["land_use"], b["land_use_source"], b["hazard"] = classify_land_use(b, pois, landuse)
    guess = b["area_m2"].map(default_floors).where(b["land_use"] != "industrial", 1)
    b["floors"] = lev.fillna(guess).clip(lower=1).astype(int)
    b["floor_area_m2"] = b["area_m2"] * b["floors"]
    b["served"] = b["area_m2"] >= MIN_SERVED_AREA_M2

    wd = packs["water"]["demand"]
    ed = packs["electricity"]["demand"]
    gd = packs["gas"]["demand"]

    # Dwelling units (residential) or unit equivalents (other uses)
    persons = wd["rooms_per_apartment"] * wd["occupants_per_room"]
    res = (b["land_use"] == "residential") & b["served"]
    area_per_unit = wd["floor_area_per_apartment_m2"]
    if target_population:
        area_per_unit = calibrate_unit_area(b.loc[res, "floor_area_m2"].values, persons, target_population)
    b.attrs["floor_area_per_dwelling_m2"] = area_per_unit
    b["units"] = np.where(b["land_use"] == "residential", residential_units(b["floor_area_m2"], area_per_unit),
                          residential_units(b["floor_area_m2"], NONRES_UNIT_M2))
    b["occupants"] = b["units"] * persons

    # Water: average day (L/day) and peak hour (L/s)
    m_w = b["land_use"].map(wd["land_use_multipliers"]).fillna(1.0)
    b["water_lpd"] = b["occupants"] * wd["daily_consumption_lpd"] * m_w * wd["seasonal_factor"]
    b["water_peak_ls"] = b["water_lpd"] / 86400.0 * wd["peak_factor"]

    # Electricity: undiversified design load (kW); diversity is applied on the network
    m_e = b["land_use"].map(ed["building_multipliers"]).fillna(1.0)
    res = b["land_use"] == "residential"
    b["elec_kw"] = np.where(res, b["units"] * ed["base_load_kw_per_apt"] * m_e,
                            b["floor_area_m2"] * ed["load_per_m2"]) * ed["future_growth_factor"]
    b["elec_service_kw"] = b["elec_kw"] * b["units"].map(lambda n: diversity_factor(n, ed["diversity_factors"]))

    # Gas: day (m3/day) and peak hour (m3/h)
    m_g = b["land_use"].map(gd["building_multipliers"]).fillna(1.0)
    b["gas_m3d"] = (b["units"] * (gd["heating_load_m3_per_apt"] + gd["cooking_load_m3_per_apt"])
                    * gd["seasonal_factor"] * gd["climate_adjustment"] * m_g)
    b["gas_peak_m3h"] = b["gas_m3d"] / 24.0 * gd["peak_hour_factor"]

    for c in ["water_lpd", "water_peak_ls", "elec_kw", "elec_service_kw", "gas_m3d", "gas_peak_m3h"]:
        b.loc[~b["served"], c] = 0.0
    return b
