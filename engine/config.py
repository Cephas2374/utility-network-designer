"""
Configuration - rule packs, crossing rules, graph rules and the local profile
merged into one dictionary per utility.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from rules.water_rules import get_water_rules  # noqa: E402
from rules.electricity_rules import get_electricity_rules  # noqa: E402
from rules.gas_rules import get_gas_rules  # noqa: E402
from rules.crossing_rules import get_crossing_rules  # noqa: E402
from rules.graph_rules import get_graph_rules  # noqa: E402
from rules.kinshasa_profile import get_kinshasa_profile, apply_profile  # noqa: E402

UTILITIES = ("water", "electricity", "gas")
DATA = ROOT / "data" / "gombe"
RESULTS = ROOT / "results" / "gombe"
FIGURES = ROOT / "paper_final" / "figures"
CRS = 32733

# Building footprints below this area (kiosks, sheds, annexes) get no service.
MIN_SERVED_AREA_M2 = 15.0


def load_rules(profile=True):
    """Return {utility: merged rule pack} plus the crossing rules."""
    packs = {"water": get_water_rules(), "electricity": get_electricity_rules(), "gas": get_gas_rules()}
    if profile:
        prof = get_kinshasa_profile()
        packs = {u: apply_profile(p, u, prof) for u, p in packs.items()}
    return packs, get_crossing_rules()


def separation(packs, u, v):
    """Minimum horizontal separation between utilities u and v (max of both packs)."""
    if u == v:
        return 0.0
    keys = {"water": "water", "electricity": "electricity", "gas": "gas"}
    a = packs[u]["constraints"].get(f"min_separation_from_{keys[v]}_m", 0.0)
    b = packs[v]["constraints"].get(f"min_separation_from_{keys[u]}_m", 0.0)
    return max(a, b)


def vertical_clearance(crossing, u, v):
    table = crossing["vertical_clearance_m"]
    return table.get((u, v), table.get((v, u), 0.3))


__all__ = ["UTILITIES", "DATA", "RESULTS", "FIGURES", "CRS", "MIN_SERVED_AREA_M2", "load_rules",
           "separation", "vertical_clearance", "get_graph_rules"]
