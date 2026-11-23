# Graph-theoretic utility network engine

QGIS-independent core of the Utility Network Designer. It designs water, electricity and gas
networks on the plane street graph of a district and checks every rule of the rule packs
(`../rules/`) as a property of plane graphs (`../rules/graph_rules.py`).

## Run on Gombe (Kinshasa)

```
python engine/prepare_gombe_data.py   # once: extract the inputs from D:\GeoAdresse_RDC_Project
python engine/run_gombe.py            # design + validation  -> results/gombe/  (~10 min)
python engine/figures_p1.py           # P1 figures from the inputs
python engine/figures_p1_results.py   # P1 figures from results/gombe/model.pkl
papers/P1_graph_formulation/build.ps1 # numbers.tex + PDF
python engine/run_water_p2.py         # P2 water analyses on model.pkl -> water_p2.json/.pkl (~25 min)
python engine/figures_p2.py           # P2 figures
papers/P2_water/build.ps1             # numbers.tex (make_numbers_p2.py) + PDF
```

## Modules

| Module | Role | Rules |
|---|---|---|
| `config.py` | rule packs + crossing rules + Kinshasa profile | – |
| `planar.py` | plane street graph: rotation, faces, dual, Euler, indices, crossing test | G1, G2 |
| `demand.py` | demand per building, population calibration | – |
| `corridor.py` | taps (face-respecting, spaced, uncrossed, shared services) and corridor graph G_c | G1, G9 |
| `design.py` | Steiner / SPT, vertex splitting, subdivision, tree flows, face augmentation, graph-Voronoi siting, arc-length placement | G3–G8, G10, G11 |
| `water.py` | inflows (head margin), GGA hydraulics, sizing, valves, hydrants, fire flow | H, G6, G11 |
| `electricity.py` | transformer cells, LV forests, cable sizing, ties, 20 kV open rings | H, G10 |
| `gas.py` | district regulators, LP / MP trees, Renouard sizing | H, G10, G11 |
| `crossings.py` | lanes, side assignment, junction and cross-section crossings, joint repair, exact avoidability | C1–C6 |
| `water_analysis.py` | P2: pressure-driven hydraulics, valve segments, segment graph, N-1, supply-forest DMAs, fire capacity, BoQ, EPANET check | H, G11 |
| `validate.py` | one independent check per rule | all |
| `figstyle.py`, `map_furniture.py` | house style of all figures (Source Sans 3, APWA utility colours) | – |

Outputs (`results/gombe/`): `networks.gpkg` (all layers, open in QGIS), `metrics.json`,
`validation.csv`, `customer_register.csv` (every customer under its V8 address), `model.pkl`, `run.log`;
P2: `water_p2.json`, `water_p2.pkl`, `run_water_p2.log`, `epanet/` (EPANET input and results).

## Changes for P2 (2026-11-01)

The P1 design and numbers are unchanged (`make_numbers_p1.py` output verified identical).

- **N-1 speed** (`Hydraulics.closed`): reduced networks are cut from the full network's incidence
  matrices with a warm start, instead of a graph copy and a new GGA per case (~5 s -> ~1 s per
  closure; identical to 3e-8 m, `dev/test_closed_ctx.py`). Described in P2 Section 4.1.
- **Tree comparison**: the Steiner tree without loop closures fails the 20 m rule for every
  K = 3..10 inflows (worst tap at the end of a 7.2 km path already at 600 mm, `dev/diag_tree_k.py`),
  so sizing cannot repair it. `design_water` got the option `graph.tree = "spt"` (shortest-path forest
  from the inflows), which is feasible on the 3 inflows of the looped design and is the tree of the
  P2 comparison (`dev/test_tree_spt.py`). Described in P2 Sections 4.5 and 5.
- **DMAs**: graph-Voronoi cells with an unbounded reach took 62 min (a full Dijkstra from every
  vertex in `greedy_cover`) and fragmented into 61-180 cells for a target of 9. Replaced by the
  Kundu-Misra partition of the peak-hour supply forest (`supply_forest`, `tree_partition`): fewest
  parts, < 1 s (`dev/test_dma_p2.py`). Proposition "Fewest districts" in P2 Section 3.5.
- **Fire capacity**: new rule key `capacity.fire_residual_mwc = 15` (DVGW W 405); the run adds
  rule keys missing from the pickled design packs (`setdefault`, design values kept).
- **Segment graph** check of the hydraulic N-1 cut-off counts (P2 Proposition "Cut-off customers").

## Changes for P3 (2026-11-01)

The P1 design and numbers are unchanged; the design code changes below are options or do not
change any length.

- **Diversity**: the stepwise table makes n g(n) drop at n = 11, 21, 51 (446 of 12,936 LV arcs carry
  less design current than a child arc). `demand.diversity_factor_monotone` (least monotone majorant)
  and the opt-in flag `demand.monotone_diversity` in `electricity._g`; P1 used the table.
- **Ring tours**: NetworkX's Christofides wrapper expands ties of the metric closure through a node
  already visited, so 13 of 24 P1 tours visit a substation twice (16 visits; no length change).
  `electricity._mv_rings` now keeps the first visit; the analyses count each load once (`visits`).
- **New rule keys** (P3): `capacity.emergency_voltage_drop_percent` (10, EN 50160),
  `capacity.emergency_transformer_loading_percent`, `capacity.max_mv_voltage_drop_percent`, and the
  section `mv_cable` (240 mm2 Al XLPE 12/20 kV planning values); added to the pickled design packs
  with `setdefault` by `run_electricity_p3.py`.
- **Restoration**: whole cells through one tie (26/361 pass) replaced by feeder-by-feeder back-feed
  (`feeder_restoration`), optionally with tie-path cables up-sized (`upsize=True`).
- **Ring re-routing**: sequential legs (287 km of cable) replaced by min-cost-flow spokes (`spokes`,
  Suurballe) plus legs that avoid the ring's streets (`reroute`); `two_edge_component` gives the
  bridge bound (Menger).

- **MV cable data** (2026-11-01, user decision: no SNEL specification available): R' derived from
  IEC 60228 (0.125 ohm/km at 20 C -> 0.160 at 90 C), X' 0.11 ohm/km (trefoil), ampacity 365 A x 0.93
  (IEC 60502-2, ground at 30 C) = 340 A.  **Infeed**: P1's single modelled infeed is kept; the bridge
  bound is also reported with the best second infeed of the 143 boundary candidates (`second_infeed`).

```
python engine/run_electricity_p3.py   # P3 analyses -> elec_p3.json/.pkl (~4 min)
python engine/figures_p3.py           # P3 figures
papers/P3_electricity/build.ps1       # numbers.tex (make_numbers_p3.py) + PDF
```

## Changes for P5, the plugin and the book (2026-11-01)

- **Lane order**: the rule pack's kerb order is now gas, water, electricity (P4).  `run_gombe.py
  --legacy-kerb-order` reproduces the design of P1-P4; the runner refuses to overwrite a stored design
  made with another order unless `--force`.
- **Pipeline** (`pipeline.py`): `design(Inputs)` and `write_outputs(result, folder)`; `run_gombe.py` is a
  thin Gombe adapter.  Verified: a full rerun reproduces all 347 numbers of metrics.json (0 differences)
  and the 78 validation rows (`results/gombe_regression`).  The stored design now carries the inflow
  candidates and the ground levels (`model_context`), so the analyses run on any district; they take
  `--results DIR`.  `node_streets` nodes ordinary street layers; OSM `highway` maps to road classes.
- **QGIS plugin** (`qgis_plugin/utility_network_graph`, built by `qgis_plugin/build_plugin.py` into a ZIP):
  Processing provider with design, analyses and drawings; tested headless in QGIS 3.44
  (`dev/test_plugin_qgis.py`) on `data/test/test_district.gpkg` (`dev/make_test_district.py`);
  screenshots by `dev/screenshots_qgis.py`.
- **Brand** (`branding/make_brand.py`, `make_portrait.py`): original artwork in the colours of the DRC flag
  (= water, gas, electricity), drawn from the Gombe design; no emblems or institutional imagery.
- **Book** (`book/`): `build_book.py` assembles P1-P4 as chapters (labels and number macros prefixed per
  chapter, floats released), `make_numbers_book.py` the software numbers, `build.ps1` the PDF,
  `cover/make_cover.py` the full-wrap cover sized from the page count.

## Changes for P4 (2026-11-01)

The P1 design and numbers are unchanged.

- **Vertical design** (`vertical.py`): every crossing of a service or junction leg with another
  utility's main is located (`Corridor.junction` now records leg positions when the record has a
  `legs` list) and the least depth profiles are computed from the shallowest utility down (cone
  envelopes).  P1 only counted "deepening points" from standard depths and clearance, without outer
  dimensions, and did not design the lowering.
- **New rule section** `crossing_rules.vertical`: outer dimensions of cables and services, largest
  gradient (1:10), trench widths.
- **Layout B**: kerb order gas, water, electricity (compatible with the depth order), sides
  re-optimised as in P1 (`vertical.corridor(m, kerb_order=...)`).
- **Street cuts** (`streetcuts.py`): water through P2's segment graph, electricity through LV trees and
  P3's ring trench cuts, gas through LP trees and MP zones, from the city gate stored in model.pkl.

```
python engine/run_corridor_p4.py      # P4 analyses -> corridor_p4.json/.pkl (~3 min)
python engine/figures_p4.py           # P4 figures and drawings
papers/P4_corridors/build.ps1         # numbers.tex (make_numbers_p4.py) + PDF
```

## Fix 2026-11-22 — plugin elevation route (plugin v2.0.1, P5 SoftwareX)

- **Symptom.** The QGIS plugin on the Gombe layers designed K=4 water inflows and passed 77/78 checks
  (fire residual 17.1 m < 20 m), against K=3 and 78/78 in P1–P4, even with the P1 lane order and the
  papers' population.
- **Cause.** Elevation. The papers sample SRTM on its native grid, smoothed by a 3x3 mean, bilinearly
  (`dem.SRTM`). The plugin read the raw pixel under each vertex (`provider.sample`) of a DEM warped to
  UTM 33S at 30 m. Over the 14,375 vertices of the design, ground levels differed by 0.65 m on average
  and 7.2 m at most (`dev/measure_dem_old_route.py`); with three inflows five nodes fell below the
  minimum pressure.
- **Fix.** `dem.GridDEM` samples any raster as `dem.SRTM` does (identical to 6e-12 m on the native grid);
  the plugin reads the band once and samples through it (parameter "Smooth the DEM", default on);
  `data/test/gombe_dem.tif` is now the SRTM tile clipped on its native grid (EPSG:4326), the old warped
  file is kept as `gombe_dem_utm30_old.tif`. `dem.SRTM` itself is unchanged, so P1–P4 are unaffected.
- **Check.** `dev/run_plugin_regression.py` (QGIS Python): K=3, 78/78, all 347 metrics identical
  (run time excluded), validation table identical.
- **Also found.** A regression script must keep a reference to the Processing provider it registers,
  or Python collects it and QGIS reports "Error creating algorithm from createInstance()".

## Fix 2026-11-23 — Bukavu transferability test (P5)

- **Tap exchange (corridor.Taps._exchange_fixed).** The first plugin run on northern Ibanda (Bukavu) left
  8 crossing service pairs (G1 failed for water, electricity and gas). In every pair, the exchange of P1's
  uncrossing lemma with the facade points kept shortens the pair and removes the crossing; `_repair`
  recomputed each facade (nearest footprint point to the new tap) before testing and then rejected the
  exchange. New last step: the lemma's own exchange (accept if shorter and no street clipped; total length
  strictly decreases, so it stops), run only when crossings survive `_retap`. Bukavu: 9 exchanges, 0 left.
  Gombe: never reached; `run_gombe.py --legacy-kerb-order` gives 0 metric differences and the 78 checks
  identical (results/gombe_regression_g13).
- **New rule G13 (services clear of other buildings).** Services are straight tap-to-facade segments and
  nothing checked them against the buildings they pass. `validate.check_services_buildings`; reported, not
  repaired. Stored design of P1-P4: 4,345 of 12,273 services cross another building (35.4 %; dev/audit_g13.py);
  Bukavu 48 %. P1 was already submitted when this was found: to be added in its revision.
- **Pressure zones.** One zone needs z_max - z_min <= p_max - p_min (60 m). Northern Ibanda spans ~105 m,
  so rule H (max pressure 80 m) fails by necessity; the engine designs one zone.
