# Two-compartment peptide ranking

A prep pipeline in `/app` is supposed to protonate two public peptides at two compartment pH values and rank them by how much net side-chain charge moves between those compartments.

The files already on disk look complete. They are not trustworthy.
Treat pH as a compartment address. Do not assume a single physiological default.

## Inputs

- `/app/data/job.json` — peptides, compartments, output paths
- `/app/data/pka_table.csv` — locked model pKa values; use only this table
- `/app/pipeline.py` — existing prep + rank script

Peptides (public): `HHHHHH`, `KSRRRAR`.
Compartments: endosome and cytosol, with pH given in the job file.
Charge scope: side chains only. No N- or C-terminal charges.
Ionic strength is recorded for provenance. It must not enter occupancy.

## Locked occupancy

For cationic side chains (His, Lys, Arg):

```
frac_protonated = 1 / (1 + 10^(pH - pKa))
charge          = +1 * frac_protonated
```

Ser and Ala: frac = 0, charge = 0.

## Required outputs

1. `/app/occupancy.csv`
   Header: `peptide,res_index,aa,compartment,pH,pKa,frac_protonated,charge`
2. `/app/summary.json` with `net_charge` and `delta_charge_endosome_minus_cytosol`
3. `/app/rank.csv` ranked by `delta_charge` descending

You may edit `/app/pipeline.py` or replace the computation.
Do not call PROPKA, MD, or a web pKa service.
Do not add biological or therapeutic claims.

You have 1800 seconds.
