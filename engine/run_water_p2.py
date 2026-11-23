"""
Water analyses of paper P2 on a stored design (default: the Gombe design, results/gombe/model.pkl).

    python engine/run_water_p2.py [--results DIR]

Writes water_p2.json and water_p2.pkl next to model.pkl.
"""
import copy
import json
import sys
import time
from pathlib import Path

import networkx as nx
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from engine import water, water_analysis as wa  # noqa: E402
from engine import store  # noqa: E402
from engine.config import RESULTS, load_rules  # noqa: E402
from engine.pipeline import model_context  # noqa: E402


def _arg_results():
    return Path(sys.argv[sys.argv.index("--results") + 1]) if "--results" in sys.argv else RESULTS


def main(results=None, log=None):
    RESULTS = Path(results) if results else _arg_results()
    t0 = time.time()
    log = log or (lambda s: print(f"[{time.strftime('%H:%M:%S')}] {s}", flush=True))
    m = store.load(RESULTS / "model.pkl")
    G, S, b, rows, packs = m["G"], m["S"], m["b"], m["rows"], m["packs"]
    # rule keys added after model.pkl was written (e.g. fire_residual_mwc) come from the live pack;
    # keys present in the design keep their design values
    for k, v in load_rules()[0]["water"]["capacity"].items():
        packs["water"]["capacity"].setdefault(k, v)
    bi = b.set_index("bldg_id")
    tv = G.graph["tap_vertex"]
    served = [r["bldg_id"] for r in rows]
    cands, elev = model_context(m)

    def prepare(N):
        q, taps = {}, {}
        for bid in served:
            t = tv[bid]
            q[t] = q.get(t, 0.0) + bi.loc[bid, "water_peak_ls"] / 1000
            taps[t] = taps.get(t, 0) + 1
        for v in N:
            N.nodes[v]["q_req"] = q.get(v, 0.0)
        return taps

    out = {}
    # -- looped design (from run_gombe) and two tree designs ---------------------------------
    Nw, wi = m["Nw"], m["wi"]
    taps = prepare(Nw)
    sources = wi["sources"]
    terminals = sorted({tv[x] for x in served}, key=str)
    p_src = packs["water"]["capacity"]["source_pressure_mwc"]
    # (1) the Steiner tree without loop closures, repaired by adding inflows (rule G4) up to K = 10:
    #     record every K, and the supply path of its worst tap (pressure bound of the tree topology)
    packs_st = copy.deepcopy(packs)
    packs_st["water"]["graph"]["bridgeless_trunk"] = False
    log("Steiner tree, K = 3 ... 10")
    rep = []
    for K in range(len(sources), 11):
        srcK = sources if K == len(sources) else water.choose_inflows(G, cands, K, elev, terminals, p_src)[0]
        NK, iK = water.design_water(G, S, served, bi, packs_st, srcK, elev, log)
        rep.append(dict(K=K, feasible=iK["feasible"], violations=iK["violations"], length_m=iK["length_m"],
                        min_pressure=iK["min_pressure_mwc"], taps_below_min=iK["taps_below_min"]))
        if K == 7:
            out["steiner_tree_worst_path"] = wa.worst_supply_path(NK, srcK, set(terminals))
        log(f"  K={K}: {rep[-1]}")
    out["steiner_tree"] = rep
    # (2) the shortest-path forest from the inflows of the looped design: the tree with the shortest
    #     supply paths, repaired by adding inflows if needed
    packs_tree = copy.deepcopy(packs_st)
    packs_tree["water"]["graph"]["tree"] = "spt"
    log("shortest-path forest")
    Nt, ti = water.design_water(G, S, served, bi, packs_tree, sources, elev, log)
    if not ti["feasible"]:
        Nt, ti = water.design_water_auto(G, S, served, bi, packs_tree, cands, elev, log)
    prepare(Nt)
    sources_of = {"looped": sources, "tree": ti["sources"]}
    designs = {"looped": (Nw, wi), "tree": (Nt, ti)}

    for name, (N, info) in designs.items():
        log(f"== {name}")
        src = sources_of[name]
        hyd = wa.Hydraulics(N, src, packs)
        r = dict(K=len(src), feasible=info["feasible"], violations=info["violations"],
                 max_velocity_ms=info["max_velocity_ms"], max_gradient_m_per_km=info["max_gradient_m_per_km"])
        base = hyd.pda()
        r["base_served_fraction"] = base["served_total"] / base["required_total"]
        r["length_m"] = sum(d["length"] for _, _, d in N.edges(data=True))
        r["mu"] = N.number_of_edges() - N.number_of_nodes() + nx.number_connected_components(N)
        r["min_pressure"] = info["min_pressure_mwc"]
        r["boq"] = wa.boq(N, packs, info["components"])
        # valve segments
        ve = wa.valve_edges(N, info["components"])
        seg_of, segs = wa.segments(N, ve, taps)
        sl = np.array([s["length"] for s in segs.values()])
        sc = np.array([s["customers"] for s in segs.values()])
        r["segments"] = dict(n=len(segs), n_valves=len(ve), length_p50=float(np.median(sl)), length_max=float(sl.max()),
                             customers_p50=float(np.median(sc)), customers_p95=float(np.percentile(sc, 95)),
                             customers_max=int(sc.max()))
        # N-1 over the segments that hold a main of at least 150 mm (and no inflow)
        big = {sid for sid, s in segs.items()
               if any(N.edges[a, b_]["diameter_mm"] >= 150 for a in s["nodes"] for b_ in N.neighbors(a)
                      if seg_of.get(b_) == sid)
               and not any(x in s["nodes"] for x in src)}
        log(f"  N-1 over {len(big)} segments")
        n1 = wa.n_minus_one(hyd, segs, taps, lambda sid, s: sid in big, log=log)
        r["n_minus_one"] = n1
        L = np.array([x["length"] for x in n1])
        cut = np.array([x["customers_cut_off"] for x in n1])
        sf = np.array([x["served_fraction"] for x in n1])
        r["n1_summary"] = dict(n=len(n1), cut_off_any=int((cut > 0).sum()), cut_off_max=int(cut.max()) if len(cut) else 0,
                               expected_cut_off=float((cut * L).sum() / L.sum()) if len(L) else 0.0,
                               unserved_mean_pct=float(100 * (1 - sf).mean()) if len(sf) else 0.0,
                               unserved_expected_pct=float(100 * ((1 - sf) * L).sum() / L.sum()) if len(L) else 0.0,
                               unserved_max_pct=float(100 * (1 - sf).max()) if len(sf) else 0.0)
        log(f"  N-1: {r['n1_summary']}")
        # customers cut off read from the segment graph alone, checked against the hydraulic N-1
        Sg = wa.segment_graph(N, seg_of, ve, src)
        pred = wa.predicted_cut_off(Sg, segs)
        mism = [x["segment"] for x in n1 if pred[x["segment"]] != x["customers_cut_off"]]
        pc = np.array(list(pred.values()))
        r["segment_graph"] = dict(n=Sg.number_of_nodes() - 1, edges=Sg.number_of_edges(),
                                  articulation=len(set(nx.articulation_points(Sg)) - {"w"}),
                                  bridges=sum(1 for _ in nx.bridges(Sg)),
                                  cut_off_any=int((pc > 0).sum()), cut_off_max=int(pc.max()),
                                  checked=len(n1), mismatches=len(mism))
        log(f"  segment graph: {r['segment_graph']}")
        r["seg_of"] = {str(k): v for k, v in seg_of.items()}
        out[name] = r
        designs[name] = (N, info, hyd, seg_of, segs)

    # -- analyses of the looped design ----------------------------------------------------
    N, info, hyd, seg_of, segs = designs["looped"]
    cap = packs["water"]["capacity"]
    log("EPANET cross-check")
    try:
        out["epanet"] = wa.epanet_check(N, sources, packs, RESULTS / "epanet")
        log(f"  max head difference {out['epanet']['max_head_diff_m']:.4f} m, "
            f"max flow difference {out['epanet']['max_flow_diff_ls']:.4f} L/s")
    except Exception as exc:  # report, do not hide
        out["epanet"] = dict(error=str(exc))
        log(f"  {exc}" if str(exc).startswith("EPANET cross-check skipped") else f"  EPANET check failed: {exc}")

    log("DMAs")
    dm = wa.dmas(N, sources, taps, packs["water"]["optimization"]["dma_target_connections"], hyd, log=log)
    out["dma"] = {k: v for k, v in dm.items() if k in ("sizes", "extra_inlets", "served_fraction", "min_pressure")}
    out["dma"].update(n=len(set(dm["dma"].values())), inlets=len(dm["inlets"]), closed=len(dm["closed"]),
                      quotient_edges=dm["quotient"].number_of_edges(), roots=len(dm["roots"]))

    log("fire flow at every hydrant")
    hn = wa.hydrant_nodes(N, info["components"])
    fire = wa.fire_all(hyd, hn, cap["fire_flow_ls"] / 1000, log=log)
    res = np.array([f["residual"] for f in fire])
    oth = np.array([f["min_other"] for f in fire])
    log("fire capacity at every hydrant (bisection)")
    qcap = []
    for k, h in enumerate(hn):
        qcap.append(wa.fire_capacity(hyd, h, cap["fire_residual_mwc"]) * 1000)
        if (k + 1) % 200 == 0:
            log(f"    capacity: {k + 1}/{len(hn)}")
    for f, qc in zip(fire, qcap):
        f["capacity_ls"] = qc
    qc = np.array(qcap)
    out["fire"] = dict(n=len(fire), residual_min=float(res.min()), residual_p05=float(np.percentile(res, 5)),
                       residual_median=float(np.median(res)),
                       below_residual=int((res < cap["fire_residual_mwc"]).sum()),
                       below_min=int((res < cap["min_pressure_mwc"]).sum()),
                       others_below_min=int((oth < cap["min_pressure_mwc"]).sum()),
                       others_min=float(oth.min()),
                       capacity_min_ls=float(qc.min()), capacity_p05_ls=float(np.percentile(qc, 5)),
                       capacity_median_ls=float(np.median(qc)),
                       capacity_below_fire_flow=int((qc < cap["fire_flow_ls"] - 0.5).sum()), all=fire)
    log(f"  fire: { {k: v for k, v in out['fire'].items() if k != 'all'} }")

    log("sensitivity")
    sens = []
    for kind, values in (("head", [-10, -5, 0, 5, 10]), ("C", [100, 110, 120, 130, 140]),
                         ("demand", [0.8, 1.0, 1.2, 1.4, 1.6])):
        for x in values:
            kw = dict(head_add=x) if kind == "head" else dict(C=x) if kind == "C" else dict(scale=x)
            p, flows, keep, d = hyd.dda(**kw)
            tp = [p[v] for v in keep if hyd.q_req.get(v, 0) > 0]
            pr = hyd.pda(**kw)
            sens.append(dict(kind=kind, value=x, min_pressure=float(min(tp)),
                             below_min_pct=float(100 * np.mean(np.array(tp) < cap["min_pressure_mwc"])),
                             served_pct=float(100 * pr["served_total"] / pr["required_total"])))
    out["sensitivity"] = sens

    out["runtime_s"] = time.time() - t0
    (RESULTS / "water_p2.json").write_text(json.dumps({k: v for k, v in out.items()}, indent=1, default=str))
    store.save(dict(out=out, Nt=Nt, ti=ti, dma=dm, hydrants=hn, sources=sources), RESULTS / "water_p2.pkl")
    log(f"done in {out['runtime_s']:.0f} s")


if __name__ == "__main__":
    main()
