# Two-compartment peptide ranking

A prep pipeline in `/app` protonates two public peptides at two compartment pH values and ranks them by net side-chain charge change.

Use `/app/data/job.json`, `/app/data/pka_table.csv`, and `/app/pipeline.py`. Peptides: `HHHHHH`, `KSRRRAR`. Side chains only. Ionic strength is provenance; it must not enter occupancy.

Cationic: `frac_protonated = 1/(1+10^(pH-pKa))`, `charge = +frac`. Ser and Ala: 0.

Write `/app/occupancy.csv` with header `peptide,res_index,aa,compartment,pH,pKa,frac_protonated,charge`.
Write `/app/summary.json` with `net_charge` and `delta_charge_endosome_minus_cytosol`.
Write `/app/rank.csv` ordered by `delta_charge` descending.

You may edit `/app/pipeline.py`. No PROPKA, MD, or web pKa. 1800 seconds.
