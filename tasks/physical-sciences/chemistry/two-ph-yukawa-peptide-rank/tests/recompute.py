"""Independent SPEC recompute. Verifier-only. Do not import solution/."""
from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path

DATA = Path(os.environ.get("APP_ROOT", "/app")) / "data"
CENTER = {"HIS": "ND1", "LYS": "NZ", "ARG": "CZ", "GLU": "CD"}
CATIONIC = {"HIS", "LYS", "ARG"}
ANIONIC = {"GLU"}
PREFACTOR = 332.0636
KAPPA_COEF = 0.329


def load_job() -> dict:
    return json.loads((DATA / "job.json").read_text())


def load_pka() -> dict[str, float]:
    table: dict[str, float] = {}
    with (DATA / "pka_table.csv").open() as fh:
        for row in csv.DictReader(fh):
            name = (row.get("resn") or "").strip()
            raw = (row.get("pKa") or "").strip()
            if name and raw:
                table[name] = float(raw)
    return table


def keep_altloc(alt: str, keep: list[str]) -> bool:
    allowed = set(keep)
    if "" in allowed:
        allowed.add(" ")
    return alt in allowed or (alt == " " and "" in keep)


def load_residues(path: Path, job: dict, chain: str) -> dict[int, dict]:
    prep = job["prep"]
    drop = set(prep["drop_resnames"])
    wanted = set(prep["records"])
    keep = list(prep["altloc_keep"])
    grouped: dict[int, dict] = {}
    for line in path.read_text().splitlines():
        tag = line[:6].strip()
        if tag not in wanted or len(line) < 54:
            continue
        if not keep_altloc(line[16], keep):
            continue
        if chain and line[21] != chain:
            continue
        resn = line[17:20].strip()
        if resn in drop or resn not in CENTER:
            continue
        resseq = int(line[22:26])
        atom = line[12:16].strip()
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        grouped.setdefault(resseq, {"resname": resn, "atoms": {}})
        grouped[resseq]["atoms"][atom] = xyz
    return grouped


def frac_charge(resn: str, pH: float, pka: dict[str, float]) -> tuple[float, float, str]:
    if resn in CATIONIC and resn in pka:
        frac = 1.0 / (1.0 + 10 ** (pH - pka[resn]))
        return frac, frac, f"{pka[resn]:.2f}"
    if resn in ANIONIC and resn in pka:
        frac = 1.0 / (1.0 + 10 ** (pka[resn] - pH))
        return frac, -frac, f"{pka[resn]:.2f}"
    return 0.0, 0.0, ""


def site_list(residues: dict[int, dict], pH: float, pka: dict[str, float]):
    rows = []
    for resseq in sorted(residues):
        rec = residues[resseq]
        resn = rec["resname"]
        atom = CENTER[resn]
        missing = atom not in rec["atoms"]
        if missing:
            frac, q, pka_s = 0.0, 0.0, f"{pka[resn]:.2f}" if resn in pka else ""
            xyz = None
        else:
            frac, q, pka_s = frac_charge(resn, pH, pka)
            xyz = rec["atoms"][atom]
        rows.append(
            {
                "resseq": resseq,
                "resname": resn,
                "atom": atom,
                "missing": missing,
                "frac": frac,
                "charge": q,
                "pKa": pka_s,
                "xyz": xyz,
            }
        )
    return rows


def yukawa(left, right, kappa: float) -> float:
    acc = 0.0
    for a in left:
        if a["missing"] or a["charge"] == 0.0 or a["xyz"] is None:
            continue
        for b in right:
            if b["missing"] or b["charge"] == 0.0 or b["xyz"] is None:
                continue
            r = math.dist(a["xyz"], b["xyz"])
            if r < 0.1:
                continue
            acc += a["charge"] * b["charge"] * math.exp(-kappa * r) / r
    return PREFACTOR * acc


def expected():
    job = load_job()
    pka = load_pka()
    kappa = KAPPA_COEF * math.sqrt(float(job["ionic_strength_M"]))
    chains = job.get("chains", {})
    loaded = {
        name: load_residues(DATA / fname, job, chains.get(name, "A"))
        for name, fname in job["pdb_files"].items()
    }
    partner = job["partner"]
    candidates = list(job.get("candidates") or job["peptides"])
    occupancy = []
    charges = {}
    for name, residues in loaded.items():
        for comp, spec in job["compartments"].items():
            pH = float(spec["pH"])
            sites = site_list(residues, pH, pka)
            charges[(name, comp)] = sites
            for s in sites:
                occupancy.append(
                    {
                        "structure": name,
                        "resseq": s["resseq"],
                        "resname": s["resname"],
                        "atom": s["atom"],
                        "compartment": comp,
                        "pH": pH,
                        "pKa": s["pKa"],
                        "frac": s["frac"],
                        "charge": s["charge"],
                        "missing_center": s["missing"],
                    }
                )
    energies = []
    summary = {"kappa": kappa, "peptides": {}}
    for pep in candidates:
        rec = {}
        for comp, spec in job["compartments"].items():
            pH = float(spec["pH"])
            U = yukawa(charges[(pep, comp)], charges[(partner, comp)], kappa)
            energies.append({"peptide": pep, "compartment": comp, "pH": pH, "U": U})
            rec[f"U_{comp}"] = U
        rec["delta_U_endosome_minus_cytosol"] = rec["U_endosome"] - rec["U_cytosol"]
        summary["peptides"][pep] = rec
    ranked = sorted(
        candidates,
        key=lambda pep: (summary["peptides"][pep]["delta_U_endosome_minus_cytosol"], candidates.index(pep)),
    )
    return job, occupancy, energies, summary, ranked
