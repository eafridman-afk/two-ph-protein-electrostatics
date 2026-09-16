#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path

APP = Path(os.environ.get("APP_ROOT", "/app"))
DATA = APP / "data"
JOB = json.loads((DATA / "job.json").read_text())
CENTER = {"HIS": "ND1", "LYS": "NZ", "ARG": "CZ", "GLU": "CD"}
CATIONIC = {"HIS", "LYS", "ARG"}
ANIONIC = {"GLU"}
COULOMB = 332.0636
TITRATABLE = set(CENTER)


def load_pka() -> dict[str, float]:
    table: dict[str, float] = {}
    with (DATA / "pka_table.csv").open() as fh:
        for row in csv.DictReader(fh):
            key = (row.get("resn") or row.get("res") or row.get("aa") or "").strip()
            raw = (row.get("pKa") or "").strip()
            if key and raw:
                table[key] = float(raw)
    return table


def _altloc_ok(alt: str, keep: list[str]) -> bool:
    allowed = {a if a != "" else " " for a in keep}
    allowed.add(" ")
    return alt in allowed


def parse_pdb(path: Path, chain: str) -> dict[int, dict]:
    prep = JOB["prep"]
    drop = set(prep["drop_resnames"])
    records = set(prep["records"])
    keep_alt = list(prep["altloc_keep"])
    residues: dict[int, dict] = {}
    for line in path.read_text().splitlines():
        rec = line[:6].strip()
        if rec not in records:
            continue
        if len(line) < 54:
            continue
        alt = line[16] if len(line) > 16 else " "
        if not _altloc_ok(alt, keep_alt):
            continue
        this_chain = line[21] if len(line) > 21 else ""
        if chain and this_chain != chain:
            continue
        name = line[12:16].strip()
        resname = line[17:20].strip()
        if resname in drop:
            continue
        resseq = int(line[22:26])
        recd = residues.setdefault(resseq, {"resname": resname, "atoms": {}})
        recd["atoms"][name] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return residues


def occupancy_charge(resname: str, pH: float, pka: dict[str, float]) -> tuple[float, float, str]:
    if resname in CATIONIC and resname in pka:
        frac = 1.0 / (1.0 + 10 ** (pH - pka[resname]))
        return frac, frac, f"{pka[resname]:.2f}"
    if resname in ANIONIC and resname in pka:
        frac = 1.0 / (1.0 + 10 ** (pka[resname] - pH))
        return frac, -frac, f"{pka[resname]:.2f}"
    return 0.0, 0.0, ""


def sites(residues: dict[int, dict]):
    out = []
    for resseq in sorted(residues):
        rec = residues[resseq]
        resname = rec["resname"]
        if resname not in TITRATABLE:
            continue
        atom = CENTER[resname]
        missing = atom not in rec["atoms"]
        if missing:
            out.append((resseq, resname, atom, None, True))
        else:
            out.append((resseq, resname, atom, rec["atoms"][atom], False))
    return out


def pair_energy(a, qa, b, qb, kappa: float) -> float:
    total = 0.0
    for site_i, qi in zip(a, qa):
        xyz_i = site_i[3]
        if site_i[4] or xyz_i is None or qi == 0.0:
            continue
        for site_j, qj in zip(b, qb):
            xyz_j = site_j[3]
            if site_j[4] or xyz_j is None or qj == 0.0:
                continue
            r = math.dist(xyz_i, xyz_j)
            if r < 0.1:
                continue
            total += qi * qj * math.exp(-kappa * r) / r
    return COULOMB * total


def _out_path(key: str) -> Path:
    raw = JOB["outputs"][key]
    path = Path(raw)
    if path.is_absolute() and os.environ.get("APP_ROOT"):
        return APP / path.name
    return path


def main() -> None:
    pka = load_pka()
    I = float(JOB["ionic_strength_M"])
    kappa = 0.329 * math.sqrt(I)
    compartments = JOB["compartments"]
    pdb_files = JOB["pdb_files"]
    chains = JOB.get("chains", {})
    candidates = list(JOB.get("candidates") or JOB["peptides"])
    mapped = {}
    for name, fname in pdb_files.items():
        mapped[name] = sites(parse_pdb(DATA / fname, chains.get(name, "A")))

    occ_rows = []
    charges: dict[tuple[str, str], list[float]] = {}
    for name in pdb_files:
        for comp, spec in compartments.items():
            pH = float(spec["pH"])
            qlist = []
            for resseq, resname, atom, xyz, missing in mapped[name]:
                if missing:
                    frac, charge, pka_s = 0.0, 0.0, ""
                    if resname in pka:
                        pka_s = f"{pka[resname]:.2f}"
                else:
                    frac, charge, pka_s = occupancy_charge(resname, pH, pka)
                qlist.append(0.0 if missing else charge)
                occ_rows.append(
                    {
                        "structure": name,
                        "resseq": resseq,
                        "resname": resname,
                        "atom": atom,
                        "compartment": comp,
                        "pH": pH,
                        "pKa": pka_s,
                        "frac": f"{frac:.10f}",
                        "charge": f"{charge:.10f}",
                        "missing_center": "true" if missing else "false",
                    }
                )
            charges[(name, comp)] = qlist

    with _out_path("occupancy").open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "structure",
                "resseq",
                "resname",
                "atom",
                "compartment",
                "pH",
                "pKa",
                "frac",
                "charge",
                "missing_center",
            ],
        )
        writer.writeheader()
        writer.writerows(occ_rows)

    energy_rows = []
    summary = {"kappa": kappa, "peptides": {}}
    partner_name = JOB["partner"]
    partner = mapped[partner_name]
    for pep in candidates:
        rec = {}
        for comp, spec in compartments.items():
            pH = float(spec["pH"])
            U = pair_energy(mapped[pep], charges[(pep, comp)], partner, charges[(partner_name, comp)], kappa)
            energy_rows.append({"peptide": pep, "compartment": comp, "pH": pH, "U": f"{U:.10f}"})
            rec[f"U_{comp}"] = U
        rec["delta_U_endosome_minus_cytosol"] = rec["U_endosome"] - rec["U_cytosol"]
        summary["peptides"][pep] = rec

    with _out_path("energy").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["peptide", "compartment", "pH", "U"])
        writer.writeheader()
        writer.writerows(energy_rows)

    _out_path("summary").write_text(json.dumps(summary, indent=2) + "\n")

    ranked = sorted(
        candidates,
        key=lambda pep: (summary["peptides"][pep]["delta_U_endosome_minus_cytosol"], candidates.index(pep)),
    )
    with _out_path("rank").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["rank", "peptide", "delta_U"])
        writer.writeheader()
        for i, pep in enumerate(ranked, start=1):
            writer.writerow(
                {
                    "rank": i,
                    "peptide": pep,
                    "delta_U": f"{summary['peptides'][pep]['delta_U_endosome_minus_cytosol']:.10f}",
                }
            )


if __name__ == "__main__":
    main()
