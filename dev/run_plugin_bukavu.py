"""Transferability test of P5: the three algorithms of the QGIS plugin on northern Ibanda (Bukavu), with
default settings only (no population, infeed or gate given; rule-pack lane order).

    "C:\\Program Files\\QGIS 3.44.12\\bin\\python-qgis-ltr.bat" dev\\run_plugin_bukavu.py [--design-only | --analyses-only]

Inputs: data/test/bukavu_ibanda_north.gpkg, data/test/bukavu_dem.tif (dev/make_bukavu_district.py).
Writes results/bukavu_ibanda_north/ (design, analyses, drawings/) and timings.json.
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


class Feedback(QgsProcessingFeedback):
    def pushInfo(self, info):  # noqa: N802
        print(f"[{time.strftime('%H:%M:%S')}]", info, flush=True)

    def reportError(self, error, fatalError=False):  # noqa: N802,N803
        print("ERROR:", error, flush=True)


gpkg = ROOT / "data" / "test" / "bukavu_ibanda_north.gpkg"
lyr = lambda name: f"{gpkg}|layername={name}"
out = ROOT / "results" / "bukavu_ibanda_north"
dem = QgsRasterLayer(str(ROOT / "data" / "test" / "bukavu_dem.tif"), "bukavu_dem")
tf = out / "timings.json"
T = json.loads(tf.read_text()) if ("--analyses-only" in sys.argv and tf.exists()) else {}
t = time.time()
if "--analyses-only" not in sys.argv:
  processing.run("utilitynetworkgraph:design", {
    "STREETS": lyr("streets"), "BUILDINGS": lyr("buildings"), "BOUNDARY": lyr("boundary"), "WATER": lyr("water"),
    "RAILWAYS": lyr("railways"), "POIS": lyr("pois"), "LANDUSE": lyr("landuse"), "DEM": dem, "DEM_SMOOTH": True,
    "POPULATION": 0.0, "KERB_ORDER": 0, "NODE": True, "OUTPUT": str(out)}, feedback=Feedback())
  T["design_s"] = round(time.time() - t)
  print(f"design: {T['design_s']} s", flush=True)
(out / "timings.json").write_text(json.dumps(T, indent=1))
if "--design-only" not in sys.argv:
    for key, params in [("analyse_electricity_corridor_s", dict(WATER=False, ELECTRICITY=True, CORRIDOR=True)),
                        ("analyse_water_s", dict(WATER=True, ELECTRICITY=False, CORRIDOR=False))]:
        t = time.time()
        processing.run("utilitynetworkgraph:analyse", {"FOLDER": str(out), **params}, feedback=Feedback())
        T[key] = round(time.time() - t)
        print(f"{key}: {T[key]} s", flush=True)
        (out / "timings.json").write_text(json.dumps(T, indent=1))
    t = time.time()
    processing.run("utilitynetworkgraph:drawings", {"FOLDER": str(out), "OUTPUT": str(out / "drawings")},
                   feedback=Feedback())
    T["drawings_s"] = round(time.time() - t)
    (out / "timings.json").write_text(json.dumps(T, indent=1))
print(json.dumps(T), flush=True)
app.exitQgis()
