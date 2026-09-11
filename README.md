# two-ph-protein-electrostatics

Harbor draft for [Terminal-Bench Science](https://www.terminal-bench-science.ai/) proposal **#979**.

Frontier agents treat proteins as pH-invariant. This task ranks public peptides by how much net side-chain charge moves between endosome and cytosol. The seed pipeline looks complete; occupancy must actually follow the compartment pH in `job.json`.

Apache-2.0. Public peptides only (`HHHHHH`, `KSRRRAR`). Not MoleculoSphere. Not a therapeutic program.

## Task

`tasks/physical-sciences/chemistry/two-ph-rank-public-peptides/`

- `instruction.md` — does not name a leaked default
- `environment/pipeline.py` — starter that writes occupancy + rank
- `environment/data/job.json` + `pka_table.csv`
- `solution/solve.sh` + `tests/`

## Status

Proposal #979 is on the TB-Science board. Human reviewer: `@MrtinoRG`. Do not open a PR on `harbor-framework/terminal-bench-science` until they say build.
