"""Rule G13 (services clear of other buildings, added 2026-11-23) applied to stored designs, without
redesigning them: the Gombe design of papers P1-P4 (results/gombe) and the Bukavu test of P5.

    python dev/audit_g13.py

Writes results/g13_audit.json.
"""
import json
import pickle
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from engine.validate import service_building_crossings  # noqa: E402

out = {}
for name in ["gombe", "bukavu_ibanda_north"]:
    f = ROOT / "results" / name / "model.pkl"
    if not f.exists():
        continue
    m = pickle.load(open(f, "rb"))
    bad, n = service_building_crossings(m["rows"], m["b"])
    out[name] = dict(services=n, crossing_other_building=bad, pct=round(100 * bad / n, 1))
    print(name, out[name])
(ROOT / "results" / "g13_audit.json").write_text(json.dumps(out, indent=1))
