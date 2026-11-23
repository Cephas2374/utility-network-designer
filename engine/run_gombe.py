"""
Run the graph-theoretic utility network model on Gombe (Kinshasa) and write the results.

    python engine/run_gombe.py                       # kerb order of the rule pack (gas, water, electricity)
    python engine/run_gombe.py --legacy-kerb-order   # the design of papers P1-P4 (gas, electricity, water)
    python engine/run_gombe.py --out DIR             # write to DIR instead of results/gombe

Outputs: networks.gpkg, metrics.json, validation.csv, customer_register.csv, model.pkl (for the
figures), run.log.  The papers P1-P4 read results/gombe; a run whose kerb order differs from the stored
design's refuses to overwrite it unless --force is given.  The design itself is engine/pipeline.py.
"""
import sys
import time
from pathlib import Path

import geopandas as gpd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from engine.config import CRS, DATA, RESULTS, load_rules  # noqa: E402
from engine import store  # noqa: E402
from engine.dem import SRTM  # noqa: E402
from engine.pipeline import MAJOR, Inputs, design, write_outputs  # noqa: E402,F401  (MAJOR: used by other scripts)
from rules.kinshasa_profile import get_kinshasa_profile  # noqa: E402

LEGACY_KERB_ORDER = ["gas", "electricity", "water"]  # the lane order of paper P1


def gombe_inputs(kerb_order=None):
    """The Gombe data of data/gombe as pipeline inputs."""
    f = DATA / "gombe_inputs.gpkg"
    read = lambda layer: gpd.read_file(f, layer=layer)
    dem = SRTM(DATA / "S05E015.hgt")
    com = gpd.read_file(DATA / "kinshasa_context.gpkg", layer="communes").to_crs(CRS)
    lingwala = com[com["name"].str.contains("Lingwala")].geometry.iloc[0]  # the MV infeed faces Lingwala
    return Inputs(streets=read("street_edges"), buildings=read("buildings"),
                  boundary=read("commune_boundary").geometry.iloc[0], water=read("water"),
                  railways=read("railways"), pois=read("pois"), landuse=read("landuse"),
                  elevation=lambda xy: float(dem.sample([xy])[0]),
                  target_population=get_kinshasa_profile()["overrides"]["common"]["target_population"],
                  mv_infeed_near=lingwala, kerb_order=kerb_order, crs=CRS)


def main():
    out = Path(sys.argv[sys.argv.index("--out") + 1]) if "--out" in sys.argv else RESULTS
    out.mkdir(parents=True, exist_ok=True)
    packs, cr = load_rules()
    kerb = LEGACY_KERB_ORDER if "--legacy-kerb-order" in sys.argv else cr["corridor_layout"]["kerb_order"]
    stored = out / "model.pkl"
    if stored.exists() and "--force" not in sys.argv:
        old = store.load(stored)["cr"]["corridor_layout"]["kerb_order"]
        if list(old) != list(kerb):
            sys.exit(f"{out} holds a design with kerb order {old} (papers P1-P4); this run uses {kerb}.  "
                     f"Use --legacy-kerb-order to reproduce it, or --force to overwrite it.")
    logf = open(out / "run.log", "w", encoding="utf-8")

    def log(msg):
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        print(line)
        logf.write(line + "\n")
        logf.flush()

    t0 = time.time()
    res = design(gombe_inputs(kerb), packs, cr, log=log)
    log("writing outputs")
    write_outputs(res, out)
    log(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
