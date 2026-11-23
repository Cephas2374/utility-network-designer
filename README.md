# Utility Network Designer (graph engine) — v2.1.0

A QGIS Processing plugin that designs, validates and analyses the **water, electricity and gas
distribution networks** of a district on its plane street graph, from three layers: streets, buildings
and a boundary. It runs the unchanged engine of the paper series on Gombe (Kinshasa) by
C.B. Mwimangire and B. Mukengere Espoir.

## Algorithms (Processing toolbox → *Utility network designer (graph engine)*)
1. **Design utility networks**: water (looped, sized, pressures), LV/MV electricity (transformers, open
   rings), gas (district regulators), corridor lanes, and every design rule checked (validation table).
2. **Analyse reliability and crossings**: valve segments, N-1, district metered areas, fire flow; LV ties,
   MV rings, trench cuts; vertical design of the corridor and losses of a street cut.
3. **Engineering drawings**: longitudinal profile, cross-section, junction plan, MV single-line diagram.

## Install
```
python qgis_plugin/build_plugin.py      # bundles engine/ and rules/ into the plugin, writes the ZIP
```
QGIS → Plugins → Manage and Install Plugins → Install from ZIP → `qgis_plugin/utility_network_graph.zip`.
Requires QGIS ≥ 3.28 (tested: 3.44 LTR, Windows 11) with numpy, scipy, pandas, shapely, geopandas,
networkx, matplotlib (shipped with the QGIS Windows installers). EPANET 2.2 is optional (hydraulic cross-check).

## Reproduce the published results
```
python engine/run_gombe.py --legacy-kerb-order --out results/gombe_regression       # runner of the papers
"C:\Program Files\QGIS 3.44.12\bin\python-qgis-ltr.bat" dev/run_plugin_regression.py   # the plugin inside QGIS
```
Both compare every metric and rule check with `results/gombe/` (the design of the papers).
`dev/test_plugin_qgis.py` runs the three algorithms on the unnoded test district `data/test/test_district.gpkg`;
`dev/run_plugin_bukavu.py` runs them on northern Ibanda, Bukavu (`data/test/bukavu_ibanda_north.gpkg`).

## Data
`data/gombe` and `data/test`: © OpenStreetMap contributors (ODbL 1.0), Google Open Buildings v3 (CC BY 4.0),
SRTM 1 arc-second (NASA/USGS, public domain). See `DATA_LICENSE.md`.

## Cite
See `CITATION.cff`.

## Licence
Code: MIT (`LICENSE`).
