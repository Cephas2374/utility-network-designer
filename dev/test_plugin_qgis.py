"""Headless test of the QGIS plugin (paper P5): register the provider and run design, analyses and
drawings on the test district.  Run with QGIS's Python:

    "C:\\Program Files\\QGIS 3.44.12\\bin\\python-qgis-ltr.bat" dev\\test_plugin_qgis.py [--analyses]
"""
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
QGIS = Path(os.environ.get("QGIS_PREFIX_PATH", r"C:\Program Files\QGIS 3.44.12\apps\qgis-ltr"))
sys.path.append(str(QGIS / "python" / "plugins"))  # the Processing framework
sys.path.insert(0, str(ROOT / "qgis_plugin"))

from qgis.core import QgsApplication, QgsProcessingFeedback  # noqa: E402

app = QgsApplication([], False)
app.initQgis()
import processing  # noqa: E402
from processing.core.Processing import Processing  # noqa: E402

Processing.initialize()
import utility_network_graph  # noqa: E402,F401  (puts lib/ on the path)
from utility_network_graph.provider import UtilityNetworkProvider  # noqa: E402

provider = UtilityNetworkProvider()
QgsApplication.processingRegistry().addProvider(provider)
print("algorithms:", [a.id() for a in provider.algorithms()])


class Feedback(QgsProcessingFeedback):
    def pushInfo(self, info):  # noqa: N802 (QGIS API)
        print("   ", info, flush=True)


gpkg = ROOT / "data" / "test" / "test_district.gpkg"
out = ROOT / "results" / "test_district"
t = time.time()
res = processing.run("utilitynetworkgraph:design", {
    "STREETS": f"{gpkg}|layername=streets", "BUILDINGS": f"{gpkg}|layername=buildings",
    "BOUNDARY": f"{gpkg}|layername=boundary", "POPULATION": 0, "KERB_ORDER": 0, "NODE": True,
    "OUTPUT": str(out)}, feedback=Feedback())
print(f"design: {time.time() - t:.0f} s ->", res)
if "--analyses" in sys.argv:
    t = time.time()
    processing.run("utilitynetworkgraph:analyse", {"FOLDER": str(out), "WATER": True, "ELECTRICITY": True,
                                                   "CORRIDOR": True}, feedback=Feedback())
    print(f"analyses: {time.time() - t:.0f} s")
    t = time.time()
    processing.run("utilitynetworkgraph:drawings", {"FOLDER": str(out), "OUTPUT": str(out / "drawings")},
                   feedback=Feedback())
    print(f"drawings: {time.time() - t:.0f} s ->", sorted(p.name for p in (out / "drawings").glob("*.pdf")))
app.exitQgis()
