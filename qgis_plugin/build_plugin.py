"""
Build the QGIS plugin: bundle the engine and the rule packs into the plugin folder and zip it.

    python qgis_plugin/build_plugin.py

Writes qgis_plugin/utility_network_graph/lib/{engine,rules} and qgis_plugin/utility_network_graph.zip,
which QGIS installs with Plugins > Manage and Install Plugins > Install from ZIP.  The bundled engine is
the code of the papers, copied unchanged; only the scripts that produce the papers' figures and numbers,
which need the paper folders, are left out.
"""
import shutil
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PLUGIN = Path(__file__).resolve().parent / "utility_network_graph"
SKIP = {"make_numbers_p1.py", "make_numbers_p2.py", "make_numbers_p3.py", "make_numbers_p4.py",
        "figures_p1.py", "prepare_gombe_data.py"}


def main():
    lib = PLUGIN / "lib"
    if lib.exists():
        shutil.rmtree(lib)
    (lib / "engine").mkdir(parents=True)
    for f in (ROOT / "engine").glob("*.py"):
        if f.name not in SKIP:
            shutil.copy2(f, lib / "engine" / f.name)
    shutil.copytree(ROOT / "engine" / "fonts", lib / "engine" / "fonts")
    shutil.copytree(ROOT / "rules", lib / "rules", ignore=shutil.ignore_patterns("__pycache__"))
    icon = PLUGIN / "icon.png"
    if not icon.exists():
        from PIL import Image
        Image.open(ROOT / "icon.png").convert("RGBA").resize((64, 64)).save(icon)
    out = PLUGIN.parent / "utility_network_graph.zip"
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for f in PLUGIN.rglob("*"):
            if f.is_file() and "__pycache__" not in f.parts:
                z.write(f, f.relative_to(PLUGIN.parent))
    n = sum(1 for _ in zipfile.ZipFile(out).namelist())
    print(f"{out} ({out.stat().st_size / 1e6:.1f} MB, {n} files)")


if __name__ == "__main__":
    main()
