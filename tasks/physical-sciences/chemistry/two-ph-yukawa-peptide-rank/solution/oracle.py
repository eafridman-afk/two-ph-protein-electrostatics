#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from pathlib import Path

APP = Path("/app")
DATA = APP / "data"
JOB = json.loads((DATA / "job.json").read_text())
CENTER = {"HIS": "ND1", "LYS": "NZ", "ARG": "CZ", "GLU": "CD", "ASP": "CG", "SER": "OG", "ALA": "CB"}
CATIONIC = {"HIS", "LYS", "ARG"}
ANIONIC = {"GLU", "ASP"}
COULOMB = 332.0636


def load_pka() -> dict[str, float]:
    table = {}
    with (DATA / "pka_table.csv").open() as fh:
        for row in csv.DictReader(fh):
            key = (row.get("resn") or row.get("res") or row.get("aa") or "").strip()
            raw = (row.get("pKa") or "").strip()
            if key and raw:
                table[key] = float(raw)
    return table


def parse_pdb(path: Path) -> dict[int, dict]:
    residues: dict[int, dict] = {}
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        name = line[12:16].strip()
        resname = line[17:20].strip()
        resseq = int(line[22:26])
        x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        rec = residues.setdefault(resseq, {"resname": resname, "atoms": {}})
        rec["atoms"][name] = (x, y, z)
    return residues


def occupancy_charge(resname: str, pH: float, pka: dict[str, float]) -> tuple[float, float, str]:
    if resname in CATIONIC and resname in pka:
        frac = 1.0 / (1.0 + 10 ** (pH - pka[resname]))
        return frac, frac, f"{pka[resname]:.2f}"
    if resname in ANIONIC and resname in pka:
        frac = 1.0 / (1.0 + 10 ** (pka[resname] - pH))
        return frac, -frac, f"{pka[resname]:.2f}"
    return 0.0, 0.0, ""


def centers(residues: dict[int, dict]) -> list[tuple[str, int, str, float, float, float, bool]]:
    out = []
    for resseq, rec in sorted(residues.items()):
        resname = rec["resname"]
        atom = CENTER.get(resname)
        missing = atom is None or atom not in rec["atoms"]
        if missing:
            out.append((resname, resseq, atom or "", 0.0, 0.0, 0.0, True))
        else:
            x, y, z = rec["atoms"][atom]
            out.append((resname, resseq, atom, x, y, z, False))
    return out


def pair_energy(a, qa, b, qb, kappa: float) -> float:
    total = 0.0
    for (res_i, i, atom_i, xi, yi, zi, miss_i), qi in zip(a, qa):
        if miss_i or qi == 0.0:
            continue
        for (res_j, j, atom_j, xj, yj, zj, miss_j), qj in zip(b, qb):
            if miss_j or qj == 0.0:
                continue
            r = math.sqrt((xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2)
            if r < 0.1:
                continue
            total += qi * qj * math.exp(-kappa * r) / r
    return COULOMB * total


def main() -> None:
    pka = load_pka()
    I = float(JOB["ionic_strength_M"])
    kappa = 0.329 * math.sqrt(I)
    compartments = JOB["compartments"]
    structures = {
        "HHHHHH": parse_pdb(DATA / "HHHHHH.pdb"),
        "KSRRRAR": parse_pdb(DATA / "KSRRRAR.pdb"),
        "EEEEEE": parse_pdb(DATA / "EEEEEE.pdb"),
    }
    mapped = {name: centers(res) for name, res in structures.items()}

    occ_rows = []
    charges: dict[tuple[str, str], list[float]] = {}
    for name, sites in mapped.items():
        for comp, spec in compartments.items():
            pH = float(spec["pH"])
            qlist = []
            for resname, resseq, atom, x, y, z, missing in sites:
                if missing:
                    frac, charge, pka_s = 0.0, 0.0, pka.get(resname, "")
                    pka_s = f"{pka_s:.2f}" if isinstance(pka_s, float) else ""
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

    with Path(JOB["outputs"]["occupancy"]).open("w", newline="") as fh:
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
    partner = mapped["EEEEEE"]
    for pep in JOB["peptides"]:
        rec = {}
        for comp, spec in compartments.items():
            pH = float(spec["pH"])
            U = pair_energy(mapped[pep], charges[(pep, comp)], partner, charges[("EEEEEE", comp)], kappa)
            energy_rows.append({"peptide": pep, "compartment": comp, "pH": pH, "U": f"{U:.10f}"})
            rec[f"U_{comp}"] = U
        rec["delta_U_endosome_minus_cytosol"] = rec["U_endosome"] - rec["U_cytosol"]
        summary["peptides"][pep] = rec

    with Path(JOB["outputs"]["energy"]).open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["peptide", "compartment", "pH", "U"])
        writer.writeheader()
        writer.writerows(energy_rows)

    Path(JOB["outputs"]["summary"]).write_text(json.dumps(summary, indent=2))

    ranked = sorted(
        JOB["peptides"],
        key=lambda pep: (
            summary["peptides"][pep]["delta_U_endosome_minus_cytosol"],
            JOB["peptides"].index(pep),
        ),
    )
    with Path(JOB["outputs"]["rank"]).open("w", newline="") as fh:
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
