# two-ph-protein-electrostatics

Harbor draft for [Terminal-Bench Science](https://www.terminal-bench-science.ai/).

Frontier agents treat proteins as pH-invariant. This packet embeds histidine occupancy in a two-compartment Yukawa rank versus a public acidic partner.

Apache-2.0. Public coordinates only. Not MoleculoSphere. Not a therapeutic program.

## Review this (Harbor PR packet)

`tasks/physical-sciences/chemistry/two-ph-yukawa-peptide-rank/`

Deposited public PDB `1LYZ` as the input. Parse ATOM/HETATM → drop water/ligands/altLoc per SPEC → occupancy at endosome 5.5 and cytosol 7.4 → intermolecular Yukawa vs `EEEEEE` → ΔU rank. `SPEC.md` states the ranking rule only.

Approved from [Discussion #1815](https://github.com/harbor-framework/terminal-bench-science/discussions/1815) (`/approve`). Do not file another Airtable. Do not PR the old hexapeptide-only files.

## Earlier draft (do not PR)

`tasks/physical-sciences/chemistry/two-ph-rank-public-peptides/` — occupancy-only draft.
