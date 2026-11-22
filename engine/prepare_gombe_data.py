"""
Extract the Gombe (Kinshasa, DR Congo) inputs used by the utility-network model
from the GeoAdresse RDC project into  data/gombe/.

Only the layers the model needs are copied, all reprojected to WGS 84 / UTM 33S
(EPSG:32733), the metric CRS used in the addressing papers on Gombe.

    python engine/prepare_gombe_data.py
"""
import gzip
import json
import shutil
from pathlib import Path

import geopandas as gpd
import pandas as pd

SRC = Path(r"D:\GeoAdresse_RDC_Project")
ADDR = SRC / "04_Kinshasa_Gombe_Data" / "data" / "gombe_addresses_v8.gpkg"
BLDG = SRC / "02_Model_Code_Data" / "data" / "african_test" / "gombe" / "gombe_buildings_osm_google.gpkg"
FACES = SRC / "06_Paper_Gombe_Oriented_Segments" / "results" / "gombe_oriented_segments.gpkg"
WATER = SRC / "06_Paper_Gombe_Oriented_Segments" / "data" / "gombe_water_osm.gpkg"
OSM = SRC / "04_Kinshasa_Gombe_Data" / "data" / "kinshasa_osm.gpkg"
GEOFABRIK = SRC / "04_Kinshasa_Gombe_Data" / "data" / "kinshasa_osm_geofabrik_layers.gpkg"
BOUNDS = SRC / "04_Kinshasa_Gombe_Data" / "data" / "kinshasa_boundaries.gpkg"
SRTM = SRC / "02_Model_Code_Data" / "data" / "global_dem" / "srtm" / "S05E015.hgt.gz"

OUT = Path(__file__).resolve().parents[1] / "data" / "gombe"
CRS = 32733


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    gpkg = OUT / "gombe_inputs.gpkg"
    if gpkg.exists():
        gpkg.unlink()
    log = {}

    def save(gdf, layer):
        gdf.to_file(gpkg, layer=layer, driver="GPKG")
        log[layer] = len(gdf)
        print(f"  {layer:<18} {len(gdf):>6}")

    commune = gpd.read_file(ADDR, layer="commune_boundary").to_crs(CRS)
    save(commune, "commune_boundary")
    aoi = commune.geometry.iloc[0]
    aoi_wgs = commune.to_crs(4326).geometry.iloc[0].buffer(0.01)
    bbox = tuple(aoi_wgs.bounds)

    # Buildings: fused OSM + Google Open Buildings footprints, with the
    # addresses of the V8 addressing model: access street, street side,
    # arc-length position along the street stroke and full address.
    b = gpd.read_file(BLDG, layer="buildings_utm33s")
    a = gpd.read_file(ADDR, layer="buildings", ignore_geometry=True)
    address_fields = ["bldg_id", "quartier", "block_code", "street", "street_code", "house_number",
                      "full_address", "side", "s_m", "perp_dist", "stroke_id", "access_street",
                      "assign_method", "quality"]
    b = b.merge(a[address_fields], on="bldg_id", how="left")
    b = b.drop(columns=["qa_overlap"], errors="ignore")
    save(b, "buildings")
    pts = gpd.read_file(ADDR, layer="address_points").to_crs(CRS)
    save(pts[address_fields + ["geometry"]], "address_points")

    edges = gpd.read_file(ADDR, layer="street_edges").to_crs(CRS)
    save(edges, "street_edges")
    save(gpd.read_file(ADDR, layer="streets").to_crs(CRS), "streets")
    save(gpd.read_file(ADDR, layer="blocks").to_crs(CRS), "blocks")
    save(gpd.read_file(ADDR, layer="quartiers").to_crs(CRS), "quartiers")
    save(gpd.read_file(FACES, layer="faces").to_crs(CRS), "faces")

    water = gpd.read_file(WATER).to_crs(CRS)
    save(water[water.intersects(aoi.buffer(1500))], "water")

    # Points of interest: schools, hospitals, public buildings are the
    # hazard zones of the gas rule pack and set land-use multipliers.
    pois = gpd.read_file(OSM, layer="pois", bbox=bbox).to_crs(CRS)
    save(pois[pois.intersects(aoi)], "pois")
    landuse = gpd.read_file(GEOFABRIK, layer="landuse_a", bbox=bbox).to_crs(CRS)
    save(landuse[landuse.intersects(aoi)], "landuse")
    rail = gpd.read_file(GEOFABRIK, layer="railways", bbox=bbox).to_crs(CRS)
    save(rail[rail.intersects(aoi.buffer(200))], "railways")

    # Context for the location map (kept in WGS 84).
    ctx = OUT / "kinshasa_context.gpkg"
    if ctx.exists():
        ctx.unlink()
    gpd.read_file(BOUNDS, layer="province").to_file(ctx, layer="province", driver="GPKG")
    gpd.read_file(BOUNDS, layer="admin_level_7").to_file(ctx, layer="communes", driver="GPKG")

    # SRTM 1-arc-second tile covering Gombe (lat -5..-4, lon 15..16).
    dem = OUT / "S05E015.hgt"
    with gzip.open(SRTM, "rb") as fi, open(dem, "wb") as fo:
        shutil.copyfileobj(fi, fo)
    log["dem"] = dem.name

    (OUT / "manifest.json").write_text(json.dumps(log, indent=2))
    print("written to", OUT)


if __name__ == "__main__":
    main()
