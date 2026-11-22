"""
Electricity analyses of paper P3 on a stored design (default: the Gombe design, results/gombe/model.pkl).

    python engine/run_electricity_p3.py [--results DIR]

Writes elec_p3.json and elec_p3.pkl next to model.pkl.
"""
import copy
import json
import pickle
import sys
import time
from collections import Counter
from pathlib import Path

import networkx as nx
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from engine import electricity_analysis as ea  # noqa: E402
from engine.config import RESULTS, load_rules  # noqa: E402
from engine.pipeline import model_context  # noqa: E402

EMERGENCY = (100, 120, 150)  # transformer loading allowed while feeding a neighbour's feeder (%)


def _arg_results():
    return Path(sys.argv[sys.argv.index("--results") + 1]) if "--results" in sys.argv else RESULTS


def main(results=None, log=None):
    RESULTS = Path(results) if results else _arg_results()
    t0 = time.time()
    log = log or (lambda s: print(f"[{time.strftime('%H:%M:%S')}] {s}", flush=True))
    m = pickle.load(open(RESULTS / "model.pkl", "rb"))
    G, b, rows, packs, F, er, infeed = m["G"], m["b"], m["rows"], m["packs"], m["Fe"], m["er"], m["infeed"]
    # rule keys added after model.pkl was written come from the live pack; design values are kept
    live = load_rules()[0]["electricity"]
    for sec in ("capacity", "mv_cable"):
        for k, v in live.get(sec, {}).items():
            packs["electricity"].setdefault(sec, {}).setdefault(k, v)
    bi = b.set_index("bldg_id")
    tr = er["transformers"]
    demand, svc_len, mv = ea.loads(G, rows, bi, packs)
    lv = ea.LV(F, er, demand, svc_len, packs)
    out = {}

    log("LV: diversity, cables, losses")
    out["inversions"] = ea.diversity_inversions(lv)
    F2, res = ea.resize_monotone(lv, packs, bi)
    out["monotone_resize"] = {k: (v if not isinstance(v, dict) else {kk: vv for kk, vv in v.items() if kk != "by_size"})
                              for k, v in res.items()}
    out["cables"] = ea.cable_cost(F, packs)
    out["cables"]["by_size"] = {str(k): v for k, v in out["cables"]["by_size"].items()}
    out["lv_losses_kw"] = ea.lv_losses(lv)
    drops = np.array([F.nodes[t]["drop_pct"] for t in demand if t in F and "drop_pct" in F.nodes[t]])
    out["lv"] = dict(customers=len(drops), drop_median=float(np.median(drops)), drop_p95=float(np.percentile(drops, 95)),
                     drop_max=float(drops.max()), above4=int((drops > 4).sum()),
                     transformers=len(lv.sites), loading=[t["loading"] for t in tr.values() if "mv_customer" not in t],
                     diversified_kva=sum(t["kva"] for t in tr.values() if "mv_customer" not in t),
                     mv_customers=len([t for t in tr.values() if "mv_customer" in t]))
    log(f"  {out['inversions']}, monotone re-size {out['monotone_resize']}")

    log("transformer outages: whole cell through one tie")
    cell = ea.tie_restoration(lv, packs)
    out["cell_restoration"] = dict(n=len(cell), ok=sum(x["ok"] for x in cell),
                                   customers=sum(x["customers"] for x in cell),
                                   customers_ok=sum(x["customers"] for x in cell if x["ok"]),
                                   no_tie=sum(x["neighbour"] is None for x in cell),
                                   fail_loading=sum(1 for x in cell if not x["ok"] and x["loading"] and x["loading"] > 1),
                                   fail_ampacity=sum(1 for x in cell if not x["ok"] and x["ampacity_ok"] is False),
                                   fail_drop=sum(1 for x in cell if not x["ok"] and x["max_drop"] and x["max_drop"] > 10))
    log(f"  {out['cell_restoration']}")

    log("transformer outages: feeder by feeder")
    costs = {float(k): v for k, v in packs["electricity"]["materials"]["costs_per_m"]["XLPE_Cu"].items()}
    fr_all = {}
    out["feeder_restoration"] = []
    for up in (False, True):
        for pct in EMERGENCY:
            p = copy.deepcopy(packs)
            p["electricity"]["capacity"]["emergency_transformer_loading_percent"] = pct
            req = {}
            fr = ea.feeder_restoration(lv, p, upsize=up, req=req)
            why = Counter()
            for x in fr:
                why.update(x["why"])
            extra = sum((costs[A] - costs[F.edges[tuple(k)]["size_mm2"]]) * F.edges[tuple(k)]["length"] for k, A in req.items())
            row = dict(upsize=up, emergency_pct=pct, restored=sum(x["restored"] for x in fr),
                       customers=sum(x["customers"] for x in fr), cells_none=sum(x["restored"] == 0 for x in fr),
                       cells_full=sum(x["restored"] == x["customers"] for x in fr), rejected=dict(why),
                       cables_upsized=len(req), upsized_km=sum(F.edges[tuple(k)]["length"] for k in req) / 1000,
                       extra_cost=extra)
            out["feeder_restoration"].append(row)
            fr_all[(up, pct)] = fr
            log(f"  {row}")
    fk = np.array([k for x in fr_all[(False, 100)] for k in x["feeder_kva"]])
    spare = np.array([t["rating_kva"] - t["kva"] for t in tr.values() if "mv_customer" not in t])
    out["feeders"] = dict(n=len(fk), kva_median=float(np.median(fk)), kva_p90=float(np.percentile(fk, 90)),
                          spare_median=float(np.median(spare)), spare_p90=float(np.percentile(spare, 90)))

    log("MV open rings")
    ra = ea.ring_analysis(er["rings"], tr, packs, G)
    out["rings"] = ra
    two = ea.two_edge_component(G, infeed)
    subs = list(dict.fromkeys(s for r in er["rings"] for s in r["tour"][1:-1]))  # each substation once
    out["ring_visits"] = dict(visits=sum(len(r["tour"]) - 2 for r in er["rings"]), substations=len(subs),
                              rings_repeating=sum(len(r["tour"]) - 2 > len(set(r["tour"][1:-1])) for r in er["rings"]),
                              in_two_groups=sum(len(r["transformers"]) for r in er["rings"])
                              - len({s for r in er["rings"] for s in r["transformers"]}))
    out["bridge_bound"] = dict(substations=len(subs), behind_bridge=sum(s not in two for s in subs),
                               kva_behind=sum(tr[s]["kva"] for s in subs if s not in two))
    log(f"  N-1 ok {sum(r['n1_ok'] for r in ra)}/{len(ra)}, bridge bound {out['bridge_bound']}")
    # the same bound with a second infeed at the best of the boundary candidates of P1
    cands, _ = model_context(m)
    out["second_infeed"] = ea.second_infeed(G, infeed, cands, subs, tr)
    out["second_infeed"]["distance_m"] = float(nx.shortest_path_length(G, infeed, out["second_infeed"]["site"],
                                                                        weight="length"))
    log(f"  second infeed: {out['second_infeed']}")

    log("trench cuts: design, disjoint spokes, trail rings")
    variants = {"design": er["rings"],
                "spokes": [ea.reroute(r, G, tr, infeed=infeed, inner=False) for r in er["rings"]],
                "trail": [ea.reroute(r, G, tr, infeed=infeed) for r in er["rings"]]}
    out["trench"] = {}
    for name, rings in variants.items():
        tc = [ea.trench_cuts(r, tr, G, infeed) for r in rings]
        rs = ea.ring_analysis(rings, tr, packs, G)
        out["trench"][name] = dict(
            cable_km=sum(r["length_m"] for r in rings) / 1000, street_km=sum(x["street_m"] for x in tc) / 1000,
            shared_km=sum(x["shared_m"] for x in tc) / 1000, exposed_km=sum(x["exposed_m"] for x in tc) / 1000,
            exp_lost_kva=float(np.mean([x["exp_lost_kva"] for x in tc])),
            exp_lost_subs=float(np.mean([x["exp_lost_subs"] for x in tc])),
            max_lost_subs=max(x["max_lost_subs"] for x in tc),
            rings_safe=sum(x["exposed_m"] == 0 for x in tc), n1_ok=sum(r["n1_ok"] for r in rs),
            I_full_max=max(r["I_full"] for r in rs), drop_full_max=max(r["drop_full"] for r in rs),
            loss_kw=sum(r["loss_balance_kw"] for r in rs))
        log(f"  {name}: {out['trench'][name]}")

    log("load growth")
    gh = ea.growth_headroom(lv, packs, ra)
    out["growth"] = {k: dict(min=float(min(v)), p05=float(np.percentile(v, 5)), median=float(np.median(v)), all=list(v))
                     for k, v in gh.items()}
    out["runtime_s"] = time.time() - t0
    (RESULTS / "elec_p3.json").write_text(json.dumps(out, indent=1, default=str))
    with open(RESULTS / "elec_p3.pkl", "wb") as fh:
        pickle.dump(dict(out=out, F2=F2, variants=variants, cell=cell, fr=fr_all, lv_sites=lv.sites), fh)
    log(f"done in {out['runtime_s']:.0f} s")


if __name__ == "__main__":
    main()
