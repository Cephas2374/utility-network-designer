"""Ground-level error of the old plugin DEM route (raw pixel of SRTM warped to UTM 30 m) against the papers' surface
(3x3-smoothed SRTM, bilinear) at every vertex of the stored Gombe design.  Needs GDAL (QGIS Python):

    "C:\Program Files\QGIS 3.44.12\bin\python-qgis-ltr.bat" dev\measure_dem_old_route.py

Writes results/plugin_regression/dem_old_route.json.
"""
import sys, pickle, json
from pathlib import Path
ROOT = Path(r"D:\Super Computer Prof Mommert\utility_network")
sys.path.insert(0, str(ROOT))
import numpy as np
from osgeo import gdal
from engine.dem import SRTM
from engine.config import DATA
s = SRTM(DATA / "S05E015.hgt")
m = pickle.load(open(ROOT / "results/gombe/model.pkl", "rb")); G = m["G"]
xy = np.array([G.nodes[v]["xy"] for v in G])
ref = s.sample(xy)
ds = gdal.Open(str(ROOT / "data/test/gombe_dem_utm30_old.tif")); a = ds.ReadAsArray().astype(float); gt = ds.GetGeoTransform()
c = np.floor((xy[:,0]-gt[0])/gt[1]).astype(int); r = np.floor((xy[:,1]-gt[3])/gt[5]).astype(int)
old = a[r, c]; d = np.abs(old - ref)
res = dict(n=len(xy), mean_abs=round(float(d.mean()),2), max_abs=round(float(d.max()),1))
out = ROOT / "results/plugin_regression"; out.mkdir(parents=True, exist_ok=True)
(out / "dem_old_route.json").write_text(json.dumps(res)); print(res)
