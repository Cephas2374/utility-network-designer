"""EPANET 2.2 for the papers, on the authors' machine: runs dev/epanet_run.py in the separate environment
dev/.venv_epanet (WNTR is not installed in the Python of the engine).  Not part of the QGIS plugin;
engine.water_analysis imports it only when WNTR is missing and this file exists.
"""
import subprocess  # nosec B404 - fixed local interpreter and script
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EPANET_PY = ROOT / "dev" / ".venv_epanet" / "Scripts" / "python.exe"
SCRIPT = ROOT / "dev" / "epanet_run.py"


def run(inp, prefix, workdir):
    if not (EPANET_PY.exists() and SCRIPT.exists()):
        raise RuntimeError("EPANET cross-check skipped: install the Python package 'wntr' to enable it")
    pr = subprocess.run([str(EPANET_PY), str(SCRIPT), str(inp), prefix],  # nosec B603 - list args, no shell
                        capture_output=True, text=True, cwd=str(workdir), shell=False)
    if pr.returncode != 0:
        raise RuntimeError(pr.stderr[-2000:])
