"""Processing provider of the plugin."""
import os

from qgis.core import QgsProcessingProvider
from qgis.PyQt.QtGui import QIcon

from .algorithms import AnalyseAlgorithm, DesignAlgorithm, DrawingsAlgorithm


class UtilityNetworkProvider(QgsProcessingProvider):
    def id(self):
        return "utilitynetworkgraph"

    def name(self):
        return "Utility network designer (graph engine)"

    def icon(self):
        return QIcon(os.path.join(os.path.dirname(__file__), "icon.png"))

    def loadAlgorithms(self):  # noqa: N802 (QGIS API)
        for alg in (DesignAlgorithm(), AnalyseAlgorithm(), DrawingsAlgorithm()):
            self.addAlgorithm(alg)
