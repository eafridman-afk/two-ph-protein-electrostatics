#!/usr/bin/env python3
"""Prep + rank public peptides. Occupancy uses a single working pH."""
from __future__ import annotations

import csv
import json
from pathlib import Path

JOB = Path("/app/data/job.json")
PKA = Path("/app/data/pka_table.csv")
WORKING_PH = 7.4


def load_pka(path: Path) -> dict[str, float | None]:
    table: dict[str, float | None] = {}
    with path.open() as fh:
        for row in csv.DictReader(fh):
            raw = row["pKa"]
            table[row["aa"]] = None if raw == "null" else float(raw)
    return table


def frac_protonated(pH: float, pKa: float) -> float:
    return 1.0 / (1.0 + 10 ** (pH - pKa))


def main() -> None:
    job = json.loads(JOB.read_text())
    pka = load_pka(PKA)
    peptides = job["peptides"]
    compartments = job["compartments"]
    rows = []
    net = {}
    for pep in peptides:
        net[pep] = {}
        for name, spec in compartments.items():
            _ = spec["pH"]
            q = 0.0
            for i, aa in enumerate(pep, start=1):
                pk = pka.get(aa)
                if pk is None:
                    frac = 0.0
                    charge = 0.0
                else:
                    frac = frac_protonated(WORKING_PH, pk)
                    charge = frac
                q += charge
                rows.append(
                    {
                        "peptide": pep,
                        "res_index": i,
                        "aa": aa,
                        "compartment": name,
                        "pH": spec["pH"],
                        "pKa": "" if pk is None else f"{pk:.2f}",
                        "frac_protonated": f"{frac:.8f}",
                        "charge": f"{charge:.8f}",
                    }
                )
            net[pep][name] = q
    occupancy_path = Path(job["outputs"]["occupancy"])
    occupancy_path.parent.mkdir(parents=True, exist_ok=True)
    with occupancy_path.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=[
                "peptide",
                "res_index",
                "aa",
                "compartment",
                "pH",
                "pKa",
                "frac_protonated",
                "charge",
            ],
        )
        w.writeheader()
        w.writerows(rows)
    summary = {"peptides": {}}
    rank_rows = []
    for pep in peptides:
        q_end = net[pep]["endosome"]
        q_cyt = net[pep]["cytosol"]
        delta = q_end - q_cyt
        summary["peptides"][pep] = {
            "net_charge": {"endosome": q_end, "cytosol": q_cyt},
            "delta_charge_endosome_minus_cytosol": delta,
        }
        rank_rows.append({"peptide": pep, "delta_charge": delta})
    rank_rows.sort(key=lambda r: (-r["delta_charge"], peptides.index(r["peptide"])))
    Path(job["outputs"]["summary"]).write_text(json.dumps(summary, indent=2) + "\n")
    with Path(job["outputs"]["rank"]).open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["rank", "peptide", "delta_charge"])
        w.writeheader()
        for i, row in enumerate(rank_rows, start=1):
            w.writerow(
                {
                    "rank": i,
                    "peptide": row["peptide"],
                    "delta_charge": f"{row['delta_charge']:.8f}",
                }
            )


if __name__ == "__main__":
    main()
