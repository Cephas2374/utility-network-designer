"""Plugin class: registers the Processing provider, and adds a menu and a toolbar that open the three
algorithms (they are also in the Processing toolbox, group 'Utility networks')."""
import os

from qgis.core import QgsApplication
from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtWidgets import QAction

from .provider import UtilityNetworkProvider

MENU = "&Utility Network Designer"
TOOLS = (("design", "Design utility networks"),
         ("analyse", "Analyse reliability and crossings"),
         ("drawings", "Engineering drawings"))


class UtilityNetworkGraphPlugin:
    def __init__(self, iface):
        self.iface = iface
        self.provider = None
        self.actions = []
        self.toolbar = None

    def initProcessing(self):  # noqa: N802 (QGIS API)
        self.provider = UtilityNetworkProvider()
        QgsApplication.processingRegistry().addProvider(self.provider)

    def initGui(self):  # noqa: N802 (QGIS API)
        self.initProcessing()
        icon = QIcon(os.path.join(os.path.dirname(__file__), "icon.png"))
        self.toolbar = self.iface.addToolBar("Utility Network Designer")
        self.toolbar.setObjectName("UtilityNetworkDesignerToolbar")
        for name, title in TOOLS:
            action = QAction(icon, title, self.iface.mainWindow())
            action.setToolTip(f"Utility Network Designer: {title}")
            action.triggered.connect(lambda checked=False, n=name: self.open(n))
            self.iface.addPluginToMenu(MENU, action)
            self.toolbar.addAction(action)
            self.actions.append(action)

    def open(self, name):
        """The dialog of one algorithm, as the Processing toolbox opens it."""
        import processing
        processing.execAlgorithmDialog(f"utilitynetworkgraph:{name}")

    def unload(self):
        for action in self.actions:
            self.iface.removePluginMenu(MENU, action)
        if self.toolbar is not None:
            self.toolbar.deleteLater()
        if self.provider is not None:
            QgsApplication.processingRegistry().removeProvider(self.provider)
