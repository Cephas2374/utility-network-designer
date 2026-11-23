"""Regression of the QGIS plugin against the papers: the 'design' algorithm on the Gombe layers with the
settings of papers P1-P4 (lane order of P1, the papers' population, the native-grid SRTM clip), then
every leaf of metrics.json and every rule check compared with results/gombe.

    "C:\\Program Files\\QGIS 3.44.12\\bin\\python-qgis-ltr.bat" dev\\run_plugin_regression.py [--dem-mode MODE]

MODE: native (default; data/test/gombe_dem.tif, smoothed + bilinear, as fixed 2026-11-22)
      raw    (same raster, no smoothing)
Writes results/plugin_regression[_MODE]/ and prints the comparison.
"""
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
QGIS = Path(os.environ.get("QGIS_PREFIX_PATH", r"C:\Program Files\QGIS 3.44.12\apps\qgis-ltr"))
sys.path.append(str(QGIS / "python" / "plugins"))
sys.path.insert(0, str(ROOT / "qgis_plugin"))

from qgis.core import QgsApplication, QgsProcessingFeedback, QgsRasterLayer  # noqa: E402

app = QgsApplication([], False)
app.initQgis()
import processing  # noqa: E402
from processing.core.Processing import Processing  # noqa: E402

Processing.initialize()
import utility_network_graph  # noqa: E402,F401
from utility_network_graph.provider import UtilityNetworkProvider  # noqa: E402

provider = UtilityNetworkProvider()  # keep a reference, or Python collects it
QgsApplication.processingRegistry().addProvider(provider)
# the rule pack bundled in the plugin (lib/ is on the path); the project root must not shadow it
from rules.kinshasa_profile import get_kinshasa_profile  # noqa: E402


class Feedback(QgsProcessingFeedback):
    def pushInfo(self, info):  # noqa: N802
        print(f"[{time.strftime('%H:%M:%S')}]", info, flush=True)

    def reportError(self, error, fatalError=False):  # noqa: N802,N803
        print("ERROR:", error, flush=True)


def diff(x, y, path="", out=None):
    out = [] if out is None else out
    if isinstance(x, dict) and isinstance(y, dict):
        for k in set(x) | set(y):
            diff(x.get(k), y.get(k), f"{path}/{k}", out)
    elif isinstance(x, list) and isinstance(y, list) and len(x) == len(y):
        for i, (a, b) in enumerate(zip(x, y)):
            diff(a, b, f"{path}[{i}]", out)
    elif isinstance(x, (int, float)) and isinstance(y, (int, float)) and not isinstance(x, bool):
        if abs(x - y) > 1e-9 * max(1.0, abs(x), abs(y)):
            out.append((path, x, y))
    elif x != y:
        out.append((path, x, y))
    return out


mode = sys.argv[sys.argv.index("--dem-mode") + 1] if "--dem-mode" in sys.argv else "native"
gpkg = ROOT / "data" / "test" / "gombe_plugin_test.gpkg"
lyr = lambda name: f"{gpkg}|layername={name}"
out = ROOT / "results" / ("plugin_regression" if mode == "native" else f"plugin_regression_{mode}")
dem = QgsRasterLayer(str(ROOT / "data" / "test" / "gombe_dem.tif"), "gombe_dem")
pop = float(get_kinshasa_profile()["overrides"]["common"]["target_population"])
t = time.time()
processing.run("utilitynetworkgraph:design", {
    "STREETS": lyr("streets"), "BUILDINGS": lyr("buildings"), "BOUNDARY": lyr("boundary"), "WATER": lyr("water"),
    "RAILWAYS": lyr("railways"), "POIS": lyr("pois"), "LANDUSE": lyr("landuse"), "DEM": dem,
    "DEM_SMOOTH": mode != "raw", "INFEED": lyr("mv_infeed"), "GATE": lyr("gas_gate"), "POPULATION": pop,
    "KERB_ORDER": 1, "NODE": True, "OUTPUT": str(out)}, feedback=Feedback())
print(f"design: {time.time() - t:.0f} s")

ref = json.loads((ROOT / "results" / "gombe" / "metrics.json").read_text())
new = json.loads((out / "metrics.json").read_text())
TIMING = {"runtime_s"}  # wall-clock time, not a result of the design
import pandas as pd  # noqa: E402
_v = pd.read_csv(out / "validation.csv")
if "G13" in set(_v["rule"]):  # G13 (2026-11-23) postdates the stored design: count the 78 original checks
    _o = _v[_v["rule"] != "G13"]
    new["validation"] = dict(passed=int(_o["ok"].sum()), total=len(_o))
d = diff({k: v for k, v in ref.items() if k not in TIMING}, {k: v for k, v in new.items() if k not in TIMING})
leaves = lambda x: sum(leaves(v) for v in x.values()) if isinstance(x, dict) else (
    sum(leaves(v) for v in x) if isinstance(x, list) else 1)
import pandas as pd  # noqa: E402
va, vb = pd.read_csv(ROOT / "results" / "gombe" / "validation.csv"), pd.read_csv(out / "validation.csv")
# G13 (added 2026-11-23) postdates the stored design of the papers: compare the 78 original checks
g13 = vb[vb["rule"] == "G13"].to_dict("records")
vb = vb[vb["rule"] != "G13"].reset_index(drop=True)
res = dict(mode=mode, K=new["water"]["K"], passed=new["validation"]["passed"], total=new["validation"]["total"],
           leaves=leaves({k: v for k, v in ref.items() if k not in TIMING}), differing=len(d), validation_identical=bool(va.equals(vb)), g13=g13, first=d[:15])
(out / "regression.json").write_text(json.dumps(res, indent=1, default=str))
print(json.dumps(res, indent=1, default=str))
app.exitQgis()
