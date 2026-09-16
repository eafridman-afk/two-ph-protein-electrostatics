import csv
import json
import math
import os
from pathlib import Path

from recompute import expected

APP = Path(os.environ.get("APP_ROOT", "/app"))
OCC = APP / "occupancy.csv"
EN = APP / "energy.csv"
SUM = APP / "summary.json"
RANK = APP / "rank.csv"
ABS = 1e-3
REL = 1e-4


def close(a, b):
    return abs(a - b) <= max(ABS, REL * max(abs(a), abs(b), 1.0))


def test_artifacts_exist():
    assert OCC.is_file()
    assert EN.is_file()
    assert SUM.is_file()
    assert RANK.is_file()


def test_headers():
    assert list(csv.DictReader(OCC.open()).fieldnames) == [
        "structure",
        "chain",
        "resseq",
        "icode",
        "resname",
        "atom",
        "compartment",
        "pH",
        "pKa",
        "frac",
        "charge",
        "missing_center",
    ]
    assert list(csv.DictReader(EN.open()).fieldnames) == ["structure", "compartment", "pH", "U"]
    assert list(csv.DictReader(RANK.open()).fieldnames) == ["rank", "structure", "delta_U"]


def test_prep_drops_water_and_counts_1lyz_residues():
    rows = list(csv.DictReader(OCC.open()))
    assert all(r["resname"] not in {"HOH", "WAT"} for r in rows)
    lyz = [r for r in rows if r["structure"] == "1LYZ" and r["compartment"] == "endosome"]
    assert len(lyz) == 129
    his = [r for r in lyz if r["resname"] == "HIS"]
    assert len(his) == 1
    assert his[0]["resseq"] == "15"
    assert his[0]["atom"] == "ND1"
    assert his[0]["missing_center"] == "false"


def test_missing_center_is_zero_for_unmapped_residues():
    rows = list(csv.DictReader(OCC.open()))
    cys = [r for r in rows if r["structure"] == "1LYZ" and r["resname"] == "CYS"]
    tyr = [r for r in rows if r["structure"] == "1LYZ" and r["resname"] == "TYR"]
    assert cys and tyr
    assert all(r["missing_center"] == "true" for r in cys + tyr)
    assert all(close(float(r["charge"]), 0.0) for r in cys + tyr)


def test_his_occupancy_moves():
    rows = list(csv.DictReader(OCC.open()))
    his = [r for r in rows if r["structure"] == "1LYZ" and r["resname"] == "HIS"]
    endo = [float(r["charge"]) for r in his if r["compartment"] == "endosome"]
    cyto = [float(r["charge"]) for r in his if r["compartment"] == "cytosol"]
    assert endo and cyto
    assert min(endo) > 0.70
    assert max(cyto) < 0.08
    for r in his:
        pH = float(r["pH"])
        expect = 1.0 / (1.0 + 10 ** (pH - 6.0))
        assert close(float(r["charge"]), expect)


def test_independent_recompute_matches_artifacts():
    job, occ, energies, summary, ranked = expected()
    got_sum = json.loads(SUM.read_text())
    assert close(got_sum["kappa"], summary["kappa"])
    assert close(got_sum["kappa"], 0.329 * math.sqrt(float(job["ionic_strength_M"])))

    got_occ = list(csv.DictReader(OCC.open()))
    assert len(got_occ) == len(occ)
    by_key = {
        (r["structure"], r["chain"], int(r["resseq"]), r["icode"], r["resname"], r["compartment"]): r
        for r in got_occ
    }
    for row in occ:
        key = (row["structure"], row["chain"], row["resseq"], row["icode"], row["resname"], row["compartment"])
        got = by_key[key]
        assert got["missing_center"] == ("true" if row["missing_center"] else "false")
        assert close(float(got["charge"]), row["charge"])
        assert close(float(got["frac"]), row["frac"])

    got_en = {(r["structure"], r["compartment"]): float(r["U"]) for r in csv.DictReader(EN.open())}
    for row in energies:
        assert close(got_en[(row["structure"], row["compartment"])], row["U"])
        assert close(got_sum["structures"][row["structure"]][f"U_{row['compartment']}"], row["U"])

    ranks = list(csv.DictReader(RANK.open()))
    assert [r["structure"] for r in ranks] == ranked
    assert int(ranks[0]["rank"]) == 1
    d0 = summary["structures"][ranked[0]]["delta_U_endosome_minus_cytosol"]
    d1 = summary["structures"][ranked[1]]["delta_U_endosome_minus_cytosol"]
    assert d0 < d1
    # Frozen single-pH His charge cannot produce this endosome shift on 1LYZ vs EEEEEE.
    assert summary["structures"]["1LYZ"]["delta_U_endosome_minus_cytosol"] < -5.0
