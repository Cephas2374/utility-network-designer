"""
Processing algorithms of the plugin: design, analyses and drawings.  Each converts the QGIS inputs to
the data structures of the engine (engine/pipeline.py) and calls it; the engine is the code of the
papers P1-P4, bundled unchanged in lib/.
"""
import os

from qgis.core import (QgsFeatureRequest, QgsProcessing, QgsProcessingAlgorithm, QgsProcessingContext,
                       QgsProcessingException, QgsProcessingOutputFolder, QgsProcessingParameterBoolean,
                       QgsProcessingParameterEnum, QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterFile, QgsProcessingParameterFolderDestination,
                       QgsProcessingParameterNumber, QgsProcessingParameterRasterLayer, QgsCoordinateTransform,
                       QgsPointXY, QgsVectorLayer)

try:  # QGIS >= 3.36
    from qgis.core import Qgis
    LINE, POLYGON, POINT = (Qgis.ProcessingSourceType.VectorLine, Qgis.ProcessingSourceType.VectorPolygon,
                            Qgis.ProcessingSourceType.VectorPoint)
    DOUBLE = Qgis.ProcessingNumberParameterType.Double
    FOLDER = Qgis.ProcessingFileParameterBehavior.Folder
except AttributeError:  # QGIS 3.28-3.34
    LINE, POLYGON, POINT = QgsProcessing.TypeVectorLine, QgsProcessing.TypeVectorPolygon, QgsProcessing.TypeVectorPoint
    DOUBLE = QgsProcessingParameterNumber.Double
    FOLDER = QgsProcessingParameterFile.Folder

KERB_ORDERS = [("rule pack: gas, water, electricity (compatible with the depth order, P4)", None),
               ("paper P1: gas, electricity, water", ["gas", "electricity", "water"])]


def _value(v):
    """QGIS attribute -> Python value (NULL -> None)."""
    if v is None:
        return None
    if hasattr(v, "isNull") and v.isNull():
        return None
    return v


def to_gdf(source, crs, context):
    """A Processing feature source as a GeoDataFrame in the CRS of the design."""
    import geopandas as gpd
    from shapely import wkb
    if source is None:
        return None
    req = QgsFeatureRequest().setDestinationCrs(crs, context.transformContext())
    names = [f.name() for f in source.fields()]
    rows, geoms = [], []
    for f in source.getFeatures(req):
        g = f.geometry()
        if g is None or g.isNull() or g.isEmpty():
            continue
        geoms.append(wkb.loads(bytes(g.asWkb())))
        rows.append({n: _value(v) for n, v in zip(names, f.attributes())})
    return gpd.GeoDataFrame(rows, geometry=geoms, crs=crs.toWkt())


def dem_sampler(layer, crs, context, smooth=True):
    """Ground level (x, y) -> m from a raster layer, sampled as the papers sample SRTM: the whole band is
    read once, smoothed by a 3x3 mean (optional) and interpolated bilinearly between pixel centres
    (engine.dem.GridDEM).  Use the DEM on its native grid (e.g. an SRTM clip in EPSG:4326) to reproduce
    the papers; a DEM warped to another CRS has already been resampled once."""
    import numpy as np
    from engine.dem import GridDEM
    prov = layer.dataProvider()
    ext, w, h = layer.extent(), layer.width(), layer.height()
    block = prov.block(1, ext, w, h)
    z = np.array([[block.value(r, c) for c in range(w)] for r in range(h)], dtype=float)
    nodata = prov.sourceNoDataValue(1) if prov.sourceHasNoDataValue(1) else None
    dx, dy = ext.width() / w, -ext.height() / h
    tr = QgsCoordinateTransform(crs, layer.crs(), context.transformContext())

    def to_grid(xs, ys):
        pts = [tr.transform(QgsPointXY(float(x), float(y))) for x, y in zip(xs, ys)]
        return np.array([p.x() for p in pts]), np.array([p.y() for p in pts])

    g = GridDEM(z, ext.xMinimum() + dx / 2, ext.yMaximum() + dy / 2, dx, dy, to_grid=to_grid, nodata=nodata,
                smooth=3 if smooth else 0)
    return lambda xy: float(g.sample([xy])[0])


def _logger(feedback):
    def log(msg):
        if feedback.isCanceled():
            raise QgsProcessingException("cancelled")
        feedback.pushInfo(str(msg))
    return log


BANNER = os.path.join(os.path.dirname(__file__), "banner.png").replace("\\", "/")


def _help(text):
    """Help panel of an algorithm: the banner of the plugin above the text."""
    return f'<p><img src="file:///{BANNER}" width="330"/></p><p>{text}</p>'  # Qt rich text: pixel widths only


class _Base(QgsProcessingAlgorithm):
    def icon(self):
        from qgis.PyQt.QtGui import QIcon
        return QIcon(os.path.join(os.path.dirname(__file__), "icon.png"))

    def group(self):
        return "Utility networks"

    def groupId(self):  # noqa: N802 (QGIS API)
        return "utilitynetworks"

    def createInstance(self):  # noqa: N802 (QGIS API)
        return type(self)()


# -----------------------------------------------------------------------------------------------
class DesignAlgorithm(_Base):
    """Design the water, electricity and gas networks of a district and validate every rule (P1)."""

    def name(self):
        return "design"

    def displayName(self):  # noqa: N802 (QGIS API)
        return "Design utility networks"

    def shortHelpString(self):  # noqa: N802 (QGIS API)
        return _help("Designs water, electricity and gas distribution networks on the plane street graph of the "
                "district and checks every rule of the rule packs (paper P1). Streets are noded at every "
                "crossing. All layers are taken to the CRS of the streets, which must be projected in metres. "
                "Writes networks.gpkg, validation.csv, customer_register.csv, metrics.json and model.pkl to the "
                "output folder and loads the network layers.")

    def initAlgorithm(self, config=None):  # noqa: N802 (QGIS API)
        P = QgsProcessingParameterFeatureSource
        self.addParameter(P("STREETS", "Street centrelines (road_class or OSM highway tag optional)", [LINE]))
        self.addParameter(P("BUILDINGS", "Building footprints (bldg_id, building_levels optional)", [POLYGON]))
        self.addParameter(P("BOUNDARY", "District boundary", [POLYGON]))
        self.addParameter(P("WATER", "Water bodies (no inflow on the bank of the largest)", [POLYGON], optional=True))
        self.addParameter(P("RAILWAYS", "Railways", [LINE], optional=True))
        self.addParameter(P("POIS", "Points of interest (OSM amenity, shop, office ...)", [POINT], optional=True))
        self.addParameter(P("LANDUSE", "Land use (OSM fclass)", [POLYGON], optional=True))
        self.addParameter(QgsProcessingParameterRasterLayer("DEM", "Ground elevation", optional=True))
        self.addParameter(QgsProcessingParameterBoolean(
            "DEM_SMOOTH", "Smooth the DEM (3x3 mean, as in the papers; damps roof and tree noise)", True))
        self.addParameter(P("INFEED", "MV infeed near (first feature)", [POINT, POLYGON], optional=True))
        self.addParameter(P("GATE", "Gas city gate near (first feature)", [POINT, POLYGON], optional=True))
        self.addParameter(QgsProcessingParameterNumber("POPULATION", "Population to calibrate the dwellings (0: none)",
                                                       DOUBLE, defaultValue=0, minValue=0))
        self.addParameter(QgsProcessingParameterEnum("KERB_ORDER", "Lane order from the kerb",
                                                     [k for k, _ in KERB_ORDERS], defaultValue=0))
        self.addParameter(QgsProcessingParameterBoolean("NODE", "Node the street layer at every crossing", True))
        self.addParameter(QgsProcessingParameterFolderDestination("OUTPUT", "Output folder"))

    def processAlgorithm(self, parameters, context, feedback):  # noqa: N802 (QGIS API)
        try:
            from engine import pipeline
        except ImportError as exc:
            raise QgsProcessingException(f"the engine needs numpy, scipy, pandas, shapely, geopandas and networkx: {exc}")
        streets_src = self.parameterAsSource(parameters, "STREETS", context)
        crs = streets_src.sourceCrs()
        if crs.isGeographic():
            raise QgsProcessingException("use a projected CRS in metres for the streets (e.g. the UTM zone)")
        src = lambda key: self.parameterAsSource(parameters, key, context)
        gdf = lambda key: to_gdf(src(key), crs, context)
        streets = gdf("STREETS")
        if self.parameterAsBool(parameters, "NODE", context):
            feedback.pushInfo("noding the streets")
            streets = pipeline.node_streets(streets)
        elif "road_class" not in streets.columns and "highway" in streets.columns:
            streets["road_class"] = streets["highway"].map(pipeline.OSM_ROAD_CLASS).fillna("residential")
        buildings = gdf("BUILDINGS")
        if "bldg_id" not in buildings.columns or buildings["bldg_id"].isna().any() or not buildings["bldg_id"].is_unique:
            buildings["bldg_id"] = [f"B{i:06d}" for i in range(len(buildings))]
        buildings = buildings.explode(index_parts=False).drop_duplicates("bldg_id").reset_index(drop=True)
        boundary = gdf("BOUNDARY").geometry.union_all()
        first = lambda key: (gdf(key).geometry.iloc[0] if src(key) is not None and src(key).featureCount() else None)
        elevation = None
        dem = self.parameterAsRasterLayer(parameters, "DEM", context)
        if dem is not None:
            elevation = dem_sampler(dem, crs, context, self.parameterAsBool(parameters, "DEM_SMOOTH", context))
        pop = self.parameterAsDouble(parameters, "POPULATION", context)
        inp = pipeline.Inputs(streets=streets, buildings=buildings, boundary=boundary, water=gdf("WATER"),
                              railways=gdf("RAILWAYS"), pois=gdf("POIS"), landuse=gdf("LANDUSE"), elevation=elevation,
                              target_population=pop or None, mv_infeed_near=first("INFEED"), gas_gate_near=first("GATE"),
                              kerb_order=KERB_ORDERS[self.parameterAsEnum(parameters, "KERB_ORDER", context)][1],
                              crs=streets.crs)
        res = pipeline.design(inp, log=_logger(feedback))
        out = self.parameterAsString(parameters, "OUTPUT", context)
        feedback.pushInfo(f"writing {out}")
        pipeline.write_outputs(res, out)
        gpkg = os.path.join(out, "networks.gpkg")
        import pyogrio
        for name, *_ in pyogrio.list_layers(gpkg):
            layer = QgsVectorLayer(f"{gpkg}|layername={name}", name, "ogr")
            if layer.isValid():
                context.temporaryLayerStore().addMapLayer(layer)
                context.addLayerToLoadOnCompletion(layer.id(), QgsProcessingContext.LayerDetails(name, context.project(),
                                                                                                    "OUTPUT"))
        v = res["metrics"]["validation"]
        feedback.pushInfo(f"{v['passed']} of {v['total']} rule checks pass")
        return {"OUTPUT": out}



# -----------------------------------------------------------------------------------------------
class AnalyseAlgorithm(_Base):
    """The analyses of papers P2-P4 on a design folder."""

    def name(self):
        return "analyse"

    def displayName(self):  # noqa: N802 (QGIS API)
        return "Analyse reliability and crossings"

    def shortHelpString(self):  # noqa: N802 (QGIS API)
        return _help("Runs on the output folder of 'Design utility networks': water reliability (valve segments, "
                "pressure-driven N-1, district metered areas, fire flow; paper P2), electricity (LV ties, MV open "
                "rings, trench cuts, load growth; paper P3), and the vertical design of the corridor with the "
                "losses of a street cut (paper P4). Results are written to the same folder as JSON.")

    def initAlgorithm(self, config=None):  # noqa: N802 (QGIS API)
        self.addParameter(QgsProcessingParameterFile("FOLDER", "Design folder (with model.pkl)", behavior=FOLDER))
        self.addParameter(QgsProcessingParameterBoolean("WATER", "Water reliability (P2, slow)", True))
        self.addParameter(QgsProcessingParameterBoolean("ELECTRICITY", "Electricity reliability (P3)", True))
        self.addParameter(QgsProcessingParameterBoolean("CORRIDOR", "Vertical design and street cuts (P4)", True))
        self.addOutput(QgsProcessingOutputFolder("OUTPUT", "Design folder"))

    def processAlgorithm(self, parameters, context, feedback):  # noqa: N802 (QGIS API)
        folder = self.parameterAsFile(parameters, "FOLDER", context)
        if not os.path.exists(os.path.join(folder, "model.pkl")):
            raise QgsProcessingException("the folder holds no model.pkl; run 'Design utility networks' first")
        log = _logger(feedback)
        if self.parameterAsBool(parameters, "WATER", context):
            from engine import run_water_p2
            run_water_p2.main(results=folder, log=log)
        if self.parameterAsBool(parameters, "ELECTRICITY", context):
            from engine import run_electricity_p3
            run_electricity_p3.main(results=folder, log=log)
        if self.parameterAsBool(parameters, "CORRIDOR", context):
            from engine import run_corridor_p4
            run_corridor_p4.main(results=folder, log=log)
        return {"OUTPUT": folder}



# -----------------------------------------------------------------------------------------------
class DrawingsAlgorithm(_Base):
    """Engineering drawings of P4 from a design folder."""

    def name(self):
        return "drawings"

    def displayName(self):  # noqa: N802 (QGIS API)
        return "Engineering drawings"

    def shortHelpString(self):  # noqa: N802 (QGIS API)
        return _help("Draws, from a design folder: the longitudinal profile of the mains of the street with the most "
                "lowering, a cross-section through one of its taps, the plan of a junction of four streets and the "
                "single-line diagram of one MV open ring (paper P4), as PDF and PNG. Runs the corridor analysis "
                "first if it is missing.")

    def initAlgorithm(self, config=None):  # noqa: N802 (QGIS API)
        self.addParameter(QgsProcessingParameterFile("FOLDER", "Design folder (with model.pkl)", behavior=FOLDER))
        self.addParameter(QgsProcessingParameterFolderDestination("OUTPUT", "Drawings folder"))

    def processAlgorithm(self, parameters, context, feedback):  # noqa: N802 (QGIS API)
        import matplotlib
        matplotlib.use("Agg")
        folder = self.parameterAsFile(parameters, "FOLDER", context)
        out = self.parameterAsString(parameters, "OUTPUT", context)
        if not os.path.exists(os.path.join(folder, "corridor_p4.pkl")):
            from engine import run_corridor_p4
            run_corridor_p4.main(results=folder, log=_logger(feedback))
        from engine import figstyle, figures_p4
        figures_p4.configure(results=folder, out=out)
        figstyle.setup()
        for name in ("profile", "section", "junction", "sld"):
            feedback.pushInfo(f"drawing {name}")
            getattr(figures_p4, f"fig_{name}")()
        return {"OUTPUT": out}

