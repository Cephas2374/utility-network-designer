"""A small test district in the shape of generic input data (paper P5, plugin test): central Gombe in a
1.4 km square, streets as whole un-noded street lines, buildings with only an id and their levels, and
the boundary.  Writes data/test/test_district.gpkg."""
import sys
from pathlib import Path

import geopandas as gpd
from shapely.geometry import box

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from engine.config import DATA, ROOT  # noqa: E402

f = DATA / "gombe_inputs.gpkg"
b = gpd.read_file(f, layer="buildings")
commune = gpd.read_file(f, layer="commune_boundary").geometry.iloc[0]
c = b.geometry.union_all().centroid
sq = box(c.x - 700, c.y - 700, c.x + 700, c.y + 700).intersection(commune)
out = ROOT / "data" / "test" / "test_district.gpkg"
out.parent.mkdir(parents=True, exist_ok=True)
if out.exists():
    out.unlink()
streets = gpd.read_file(f, layer="streets")  # whole streets: not noded at their crossings
streets = gpd.clip(streets[["road_class", "name", "geometry"]], sq).explode(index_parts=False)
streets = streets[streets.geometry.length > 1]
streets.to_file(out, layer="streets", driver="GPKG")
bb = b[b.geometry.representative_point().within(sq)][["bldg_id", "building_levels", "geometry"]]
bb.to_file(out, layer="buildings", driver="GPKG")
gpd.GeoDataFrame(geometry=[sq], crs=b.crs).to_file(out, layer="boundary", driver="GPKG")
print(f"{out}: {len(streets)} street lines, {len(bb)} buildings, {sq.area / 1e6:.2f} km2")
