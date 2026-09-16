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
CENTER = {"HIS": "ND1", "LYS": "NZ", "ARG": "CZ", "GLU": "CD", "ASP": "CG", "SER": "OG", "ALA": "CB"}
CATIONIC = {"HIS", "LYS", "ARG"}
ANIONIC = {"GLU", "ASP"}
COULOMB = 332.0636


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


def parse_pdb(path: Path) -> dict[tuple[str, int, str, str], dict]:
    prep = JOB["prep"]
    drop = set(prep["drop_resnames"])
    records = set(prep["records"])
    keep_alt = list(prep["altloc_keep"])
    residues: dict[tuple[str, int, str, str], dict] = {}
    for line in path.read_text().splitlines():
        rec = line[:6].strip()
        if rec not in records:
            continue
        if len(line) < 54:
            continue
        alt = line[16] if len(line) > 16 else " "
        if not _altloc_ok(alt, keep_alt):
            continue
        name = line[12:16].strip()
        resname = line[17:20].strip()
        if resname in drop:
            continue
        chain = line[21] if len(line) > 21 else ""
        resseq = int(line[22:26])
        icode = line[26] if len(line) > 26 else " "
        if icode == " ":
            icode = ""
        x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        key = (chain, resseq, icode, resname)
        recd = residues.setdefault(key, {"resname": resname, "atoms": {}})
        recd["atoms"][name] = (x, y, z)
    return residues


def occupancy_charge(resname: str, pH: float, pka: dict[str, float]) -> tuple[float, float, str]:
    if resname in CATIONIC and resname in pka:
        frac = 1.0 / (1.0 + 10 ** (pH - pka[resname]))
        return frac, frac, f"{pka[resname]:.2f}"
    if resname in ANIONIC and resname in pka:
        frac = 1.0 / (1.0 + 10 ** (pka[resname] - pH))
        return frac, -frac, f"{pka[resname]:.2f}"
    return 0.0, 0.0, ""


def centers(residues: dict[tuple[str, int, str, str], dict]):
    out = []
    for key in sorted(residues, key=lambda k: (k[0], k[1], k[2], k[3])):
        chain, resseq, icode, resname = key
        rec = residues[key]
        atom = CENTER.get(resname, "")
        missing = atom == "" or atom not in rec["atoms"]
        if missing:
            out.append((chain, resseq, icode, resname, atom, 0.0, 0.0, 0.0, True))
        else:
            x, y, z = rec["atoms"][atom]
            out.append((chain, resseq, icode, resname, atom, x, y, z, False))
    return out


def pair_energy(a, qa, b, qb, kappa: float) -> float:
    total = 0.0
    for site_i, qi in zip(a, qa):
        miss_i = site_i[8]
        if miss_i or qi == 0.0:
            continue
        xi, yi, zi = site_i[5], site_i[6], site_i[7]
        for site_j, qj in zip(b, qb):
            miss_j = site_j[8]
            if miss_j or qj == 0.0:
                continue
            xj, yj, zj = site_j[5], site_j[6], site_j[7]
            r = math.sqrt((xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2)
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
    structures = {name: parse_pdb(DATA / fname) for name, fname in pdb_files.items()}
    mapped = {name: centers(res) for name, res in structures.items()}

    occ_rows = []
    charges: dict[tuple[str, str], list[float]] = {}
    struct_order = list(pdb_files)
    for name in struct_order:
        sites = mapped[name]
        for comp, spec in compartments.items():
            pH = float(spec["pH"])
            qlist = []
            for chain, resseq, icode, resname, atom, x, y, z, missing in sites:
                if missing:
                    frac, charge = 0.0, 0.0
                    pka_s = pka.get(resname, "")
                    pka_s = f"{pka_s:.2f}" if isinstance(pka_s, float) else ""
                else:
                    frac, charge, pka_s = occupancy_charge(resname, pH, pka)
                qlist.append(0.0 if missing else charge)
                occ_rows.append(
                    {
                        "structure": name,
                        "chain": chain,
                        "resseq": resseq,
                        "icode": icode,
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
            ],
        )
        writer.writeheader()
        writer.writerows(occ_rows)

    energy_rows = []
    summary = {"kappa": kappa, "structures": {}}
    partner_name = JOB["partner"]
    partner = mapped[partner_name]
    for pep in JOB["peptides"]:
        rec = {}
        for comp, spec in compartments.items():
            pH = float(spec["pH"])
            U = pair_energy(mapped[pep], charges[(pep, comp)], partner, charges[(partner_name, comp)], kappa)
            energy_rows.append({"structure": pep, "compartment": comp, "pH": pH, "U": f"{U:.10f}"})
            rec[f"U_{comp}"] = U
        rec["delta_U_endosome_minus_cytosol"] = rec["U_endosome"] - rec["U_cytosol"]
        summary["structures"][pep] = rec

    with _out_path("energy").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["structure", "compartment", "pH", "U"])
        writer.writeheader()
        writer.writerows(energy_rows)

    _out_path("summary").write_text(json.dumps(summary, indent=2) + "\n")

    ranked = sorted(
        JOB["peptides"],
        key=lambda pep: (
            summary["structures"][pep]["delta_U_endosome_minus_cytosol"],
            JOB["peptides"].index(pep),
        ),
    )
    with _out_path("rank").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["rank", "structure", "delta_U"])
        writer.writeheader()
        for i, pep in enumerate(ranked, start=1):
            writer.writerow(
                {
                    "rank": i,
                    "structure": pep,
                    "delta_U": f"{summary['structures'][pep]['delta_U_endosome_minus_cytosol']:.10f}",
                }
            )


if __name__ == "__main__":
    main()
