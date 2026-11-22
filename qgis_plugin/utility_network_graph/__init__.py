"""QGIS plugin: the graph-theoretic utility network engine as Processing algorithms (paper P5)."""
import os
import sys

# the engine and the rule packs are bundled in lib/ by qgis_plugin/build_plugin.py
LIB = os.path.join(os.path.dirname(__file__), "lib")
if LIB not in sys.path:
    sys.path.insert(0, LIB)


def classFactory(iface):  # noqa: N802 (QGIS API)
    from .plugin import UtilityNetworkGraphPlugin
    return UtilityNetworkGraphPlugin(iface)
