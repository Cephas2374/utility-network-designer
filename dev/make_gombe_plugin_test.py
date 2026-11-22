"""The Gombe data of the papers as plugin inputs: one GeoPackage with a layer per parameter of the
'Design utility networks' dialog, and the elevation as a GeoTIFF beside it.

    python dev/make_gombe_plugin_test.py                      # the GeoPackage
    python-qgis-ltr dev/make_gombe_plugin_test.py --dem       # the GeoTIFF (needs GDAL of QGIS)

Writes data/test/gombe_plugin_test.gpkg (WGS 84 / UTM 33S, EPSG:32733) and data/test/gombe_dem.tif
(the SRTM tile clipped on its native grid, EPSG:4326).
"""
import pickle
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))  # the stored design unpickles engine objects
DATA = ROOT / "data" / "gombe"
OUT = ROOT / "data" / "test"


def gpkg():
    import geopandas as gpd
    from shapely.geometry import Point
    f = DATA / "gombe_inputs.gpkg"
    out = OUT / "gombe_plugin_test.gpkg"
    OUT.mkdir(parents=True, exist_ok=True)
    if out.exists():
        out.unlink()
    save = lambda g, name: (g.to_file(out, layer=name, driver="GPKG"), print(f"  {name:<10} {len(g):>6}"))
    save(gpd.read_file(f, layer="street_edges"), "streets")
    save(gpd.read_file(f, layer="buildings"), "buildings")
    save(gpd.read_file(f, layer="commune_boundary"), "boundary")
    save(gpd.read_file(f, layer="water"), "water")
    save(gpd.read_file(f, layer="railways"), "railways")
    pois = gpd.read_file(f, layer="pois")
    save(pois.drop(columns=[c for c in pois.columns if c == "tags"]), "pois")
    save(gpd.read_file(f, layer="landuse"), "landuse")
    # the MV infeed and the gas gate of the design of the papers (results/gombe/model.pkl)
    m = pickle.load(open(ROOT / "results" / "gombe" / "model.pkl", "rb"))
    G, crs = m["G"], gpd.read_file(f, layer="streets", rows=1).crs
    save(gpd.GeoDataFrame([{"name": "MV infeed 20 kV (papers)"}], geometry=[Point(G.nodes[m["infeed"]]["xy"])], crs=crs),
         "mv_infeed")
    save(gpd.GeoDataFrame([{"name": "Gas city gate (papers, hypothetical)"}],
                          geometry=[Point(G.nodes[m["gate"]]["xy"])], crs=crs), "gas_gate")
    print(out)


def dem():
    """The SRTM tile clipped on its native grid (EPSG:4326, 1 arc-second, no resampling).

    Fix 2026-11-22: this used to warp SRTM to UTM 33S at 30 m (bilinear), so the plugin sampled a surface
    resampled once more than the papers' one; with the raw pixel read, the plugin designed K=4 and passed
    77/78 on Gombe.  On the native grid, with the plugin's 3x3-mean + bilinear sampling (engine.dem.GridDEM),
    the plugin reproduces the papers (dev/run_plugin_regression.py)."""
    from osgeo import gdal
    import geopandas as gpd
    b = gpd.read_file(OUT / "gombe_plugin_test.gpkg", layer="boundary").to_crs(4326).total_bounds
    pad = 0.01  # degrees (~1.1 km)
    out = OUT / "gombe_dem.tif"
    for f in (out, out.with_suffix(".tif.aux.xml")):
        f.unlink(missing_ok=True)
    gdal.Translate(str(out), str(DATA / "S05E015.hgt"), projWin=(b[0] - pad, b[3] + pad, b[2] + pad, b[1] - pad),
                   format="GTiff", creationOptions=["COMPRESS=DEFLATE"])
    print(out)


if __name__ == "__main__":
    dem() if "--dem" in sys.argv else gpkg()
