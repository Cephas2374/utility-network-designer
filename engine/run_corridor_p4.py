"""
Corridor analyses of paper P4 on the Gombe design (results/gombe/model.pkl): the vertical design at
crossings for the kerb order of P1 and for the order compatible with the depth order, and the losses of
the three utilities when a street is dug up.

    python engine/run_corridor_p4.py [--results DIR]

Writes corridor_p4.json and corridor_p4.pkl next to model.pkl.
"""
import copy
import json
import sys
import time
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from engine import electricity_analysis as ea, streetcuts, vertical  # noqa: E402
from engine import store  # noqa: E402
from engine.config import RESULTS, load_rules  # noqa: E402

LAYOUTS = {"p1": None, "compatible": ["gas", "water", "electricity"]}
GRADIENTS = (0.05, 0.1, 0.2)


def _arg_results():
    return Path(sys.argv[sys.argv.index("--results") + 1]) if "--results" in sys.argv else RESULTS


def main(results=None, log=None):
    RESULTS = Path(results) if results else _arg_results()
    t0 = time.time()
    log = log or (lambda s: print(f"[{time.strftime('%H:%M:%S')}] {s}", flush=True))
    m = store.load(RESULTS / "model.pkl")
    live = load_rules()
    m["cr"] = dict(m["cr"], vertical=live[1]["vertical"])  # rule section added for P4
    for sec in ("capacity", "mv_cable"):
        for k, v in live[0]["electricity"].get(sec, {}).items():
            m["packs"]["electricity"].setdefault(sec, {}).setdefault(k, v)
    out, keep = {}, {}

    for name, ko in LAYOUTS.items():
        if ko is not None and list(ko) == list(m["cr"]["corridor_layout"]["kerb_order"]):
            log(f"layout {name}: the design already has this kerb order")  # designs made after P4
            continue
        log(f"layout {name}: kerb order {ko or m['cr']['corridor_layout']['kerb_order']}")
        cor = vertical.corridor(m, kerb_order=ko, log=lambda s: None)
        rp = cor.report()
        lay = dict(kerb_order=list(cor.kerb),
                   road={u: rp["road_service"].get(u, 0) + rp["road_junction"].get(u, 0) for u in cor.U},
                   inter=rp["inter"], verge_over=rp["lane_outside_verge_sides"],
                   gas_violations=rp["gas_building_clearance_violations"],
                   dual_km=sum(rp["dual_mains_length_m"].values()) / 1000,
                   service_km={u: sum(cor.service_length(u, r) for r in m["rows"]) / 1000 for u in cor.U})
        for g in GRADIENTS:
            cr = copy.deepcopy(m["cr"])
            cr["vertical"]["max_gradient"] = g
            res = vertical.analyse(cor, m, m["packs"], cr)
            lay[f"grad_{g}"] = dict(totals=res["totals"], events=res["events"], legs=res["legs"])
            if g == 0.1:
                keep[name] = dict(prof={k: v for k, v in res["prof"].items() if k[0] == "main"},
                                  req={k: v for k, v in res["req"].items() if k[0] == "main"},
                                  sigma=cor.sigma, dual=cor.dual, lanes=cor.lane, kerb=list(cor.kerb))
            log(f"  gradient {g}: " + ", ".join(
                f"{u} mains {t['mains_lowered_m'] / 1000:.1f} km {t['mains_extra_m3']:.0f} m3"
                for u, t in res["totals"].items()))
        out[name] = lay

    log("street cuts")
    G = m["G"]
    gate = m["gate"]  # the city gate of the design
    demand_lv, _, mv_taps = ea.loads(G, m["rows"], m["b"].set_index("bldg_id"), m["packs"])
    sc = streetcuts.street_cuts(m, demand_lv, mv_taps, gate)
    taps = sc["taps"]
    S = m["S"]
    eids = set(sc["water"]) | set(sc["electricity"]) | set(sc["gas"])
    rows = []
    for e in eids:
        L = S.edges[e][3]["length"]
        sets = {u: sc[u].get(e, set()) for u in ("water", "electricity", "gas")}
        n = {u: sum(taps[t] for t in s) for u, s in sets.items()}
        both = lambda a, b: sum(taps[t] for t in sets[a] & sets[b])
        anyset = sets["water"] | sets["electricity"] | sets["gas"]
        three = sets["water"] & sets["electricity"] & sets["gas"]
        two = {t for t in anyset if sum(t in s for s in sets.values()) >= 2}
        rows.append(dict(eid=e, length=L, water=n["water"], electricity=n["electricity"], gas=n["gas"],
                         we=both("water", "electricity"), wg=both("water", "gas"), eg=both("electricity", "gas"),
                         two=sum(taps[t] for t in two), three=sum(taps[t] for t in three),
                         any=sum(taps[t] for t in anyset)))
    Ltot = sum(r["length"] for r in rows)
    ex = lambda k: sum(r[k] * r["length"] for r in rows) / Ltot
    out["street_cuts"] = dict(edges=len(rows), street_km=Ltot / 1000, connections=sum(taps.values()),
                              expected={k: ex(k) for k in ("water", "electricity", "gas", "we", "wg", "eg",
                                                           "two", "three", "any")},
                              max={k: max(r[k] for r in rows) for k in ("water", "electricity", "gas", "two", "three", "any")},
                              edges_two=sum(r["two"] > 0 for r in rows), km_two=sum(r["length"] for r in rows if r["two"] > 0) / 1000,
                              edges_three=sum(r["three"] > 0 for r in rows),
                              km_three=sum(r["length"] for r in rows if r["three"] > 0) / 1000,
                              gate=str(gate))
    log(f"  {out['street_cuts']}")
    out["gas"] = {k: v for k, v in json.loads((RESULTS / "metrics.json").read_text())["gas"].items()
                  if not isinstance(v, (list, dict))}
    out["runtime_s"] = time.time() - t0
    (RESULTS / "corridor_p4.json").write_text(json.dumps(out, indent=1, default=str))
    store.save(dict(out=out, keep=keep, cuts=rows, gate=gate), RESULTS / "corridor_p4.pkl")
    log(f"done in {out['runtime_s']:.0f} s")


if __name__ == "__main__":
    main()
