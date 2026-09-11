import csv
import json
from pathlib import Path

OCC = Path("/app/occupancy.csv")
SUM = Path("/app/summary.json")
RANK = Path("/app/rank.csv")
PIPE = Path("/app/pipeline.py")

PKA = {"H": 6.00, "K": 10.40, "R": 12.50}
PH = {"endosome": 5.5, "cytosol": 7.4}


def frac(pH: float, pKa: float) -> float:
    return 1.0 / (1.0 + 10 ** (pH - pKa))


def test_files_exist():
    assert OCC.is_file()
    assert SUM.is_file()
    assert RANK.is_file()


def test_his_occupancy_moves_between_compartments():
    rows = list(csv.DictReader(OCC.open()))
    his = [r for r in rows if r["peptide"] == "HHHHHH" and r["aa"] == "H"]
    endo = [float(r["frac_protonated"]) for r in his if r["compartment"] == "endosome"]
    cyto = [float(r["frac_protonated"]) for r in his if r["compartment"] == "cytosol"]
    assert endo and cyto
    assert min(endo) > 0.70
    assert max(cyto) < 0.08
    assert abs(endo[0] - frac(5.5, 6.0)) < 1e-6
    assert abs(cyto[0] - frac(7.4, 6.0)) < 1e-6


def test_polyarg_delta_is_small():
    summary = json.loads(SUM.read_text())
    d_h = summary["peptides"]["HHHHHH"]["delta_charge_endosome_minus_cytosol"]
    d_k = summary["peptides"]["KSRRRAR"]["delta_charge_endosome_minus_cytosol"]
    assert d_h > 3.0
    assert d_k < 0.05


def test_rank_order():
    rows = list(csv.DictReader(RANK.open()))
    assert rows[0]["peptide"] == "HHHHHH"
    assert rows[1]["peptide"] == "KSRRRAR"


def test_occupancy_uses_compartment_ph_not_a_single_default():
    rows = list(csv.DictReader(OCC.open()))
    his = [r for r in rows if r["peptide"] == "HHHHHH" and r["aa"] == "H"]
    by_comp = {}
    for r in his:
        by_comp.setdefault(r["compartment"], set()).add(round(float(r["frac_protonated"]), 6))
    assert by_comp["endosome"] != by_comp["cytosol"]
