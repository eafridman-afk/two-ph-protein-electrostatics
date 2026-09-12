import csv
import json
import math
from pathlib import Path

OCC = Path("/app/occupancy.csv")
EN = Path("/app/energy.csv")
SUM = Path("/app/summary.json")
RANK = Path("/app/rank.csv")
DATA = Path("/app/data")
ABS = 1e-3
REL = 1e-4
CENTER = {"HIS": "ND1", "LYS": "NZ", "ARG": "CZ", "GLU": "CD", "ASP": "CG", "SER": "OG", "ALA": "CB"}
CATIONIC = {"HIS", "LYS", "ARG"}
ANIONIC = {"GLU", "ASP"}
COULOMB = 332.0636


def close(a, b):
    return abs(a - b) <= max(ABS, REL * max(abs(a), abs(b), 1.0))


def test_artifacts_exist():
    assert OCC.is_file()
    assert EN.is_file()
    assert SUM.is_file()
    assert RANK.is_file()


def test_his_occupancy_moves():
    rows = list(csv.DictReader(OCC.open()))
    his = [r for r in rows if r["structure"] == "HHHHHH" and r["resname"] == "HIS"]
    endo = [float(r["charge"]) for r in his if r["compartment"] == "endosome"]
    cyto = [float(r["charge"]) for r in his if r["compartment"] == "cytosol"]
    assert endo and cyto
    assert min(endo) > 0.70
    assert max(cyto) < 0.08
    for r in his:
        pH = float(r["pH"])
        expect = 1.0 / (1.0 + 10 ** (pH - 6.0))
        assert close(float(r["charge"]), expect)


def test_glu_is_anionic_and_kappa_matches_job():
    summary = json.loads(SUM.read_text())
    job = json.loads((DATA / "job.json").read_text())
    assert close(summary["kappa"], 0.329 * math.sqrt(float(job["ionic_strength_M"])))
    rows = list(csv.DictReader(OCC.open()))
    glu = [r for r in rows if r["structure"] == "EEEEEE" and r["resname"] == "GLU"]
    assert glu and all(float(r["charge"]) < 0 for r in glu)


def test_rank_and_delta_u_require_two_pH_propagation():
    summary = json.loads(SUM.read_text())
    d_his = summary["peptides"]["HHHHHH"]["delta_U_endosome_minus_cytosol"]
    d_arg = summary["peptides"]["KSRRRAR"]["delta_U_endosome_minus_cytosol"]
    assert d_his < d_arg
    ranks = list(csv.DictReader(RANK.open()))
    assert ranks[0]["peptide"] == "HHHHHH"
    assert int(ranks[0]["rank"]) == 1
    energy = {(r["peptide"], r["compartment"]): float(r["U"]) for r in csv.DictReader(EN.open())}
    assert close(energy[("HHHHHH", "endosome")], summary["peptides"]["HHHHHH"]["U_endosome"])
    # Frozen single-pH charges cannot produce a large His-driven ΔU.
    assert d_his < -1.0
